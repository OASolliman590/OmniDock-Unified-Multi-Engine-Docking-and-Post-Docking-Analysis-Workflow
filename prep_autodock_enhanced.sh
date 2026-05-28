#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Enhanced AutoDock Preparation Script
# =============================================================================
# This script prepares ligands and receptors for AutoDock Vina with advanced
# features including comprehensive error handling,
# and flexible input/output configurations.
# =============================================================================

# ── Configuration and Setup ────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${1:-$SCRIPT_DIR/autodock_config.json}"
LOG_FILE="${SCRIPT_DIR}/autodock_prep_$(date +%Y%m%d_%H%M%S).log"

# Default configuration
DEFAULT_CONFIG='{
  "input": {
    "ligands": {
      "path": "./ligands_raw",
      "formats": ["sdf", "mol", "mol2", "pdb"],
      "in_same_folder": false
    },
    "receptors": {
      "path": "./receptors_raw", 
      "formats": ["pdb"],
      "in_same_folder": false
    }
  },
  "output": {
    "ligands": "./ligands_prep",
    "receptors": "./receptors_prep",
    "logs": "./logs"
  },
  "preparation": {
    "force_field": "AMBER",
    "ph": 7.4,
    "allow_bad_res": true,
    "default_altloc": "A",
    "receptor_use_pdb2pqr": false,
    "ligand_preparation_backend": "engine_aware_full",
    "ligand_preparation_profile": "engine_aware_full",
    "selected_engines": [],
    "autodocktools_prepare_ligand4": "",
    "autodocktools_prepare_receptor4": "",
    "autodocktools_python": ""
  },
  "quality_control": {
    "validate_outputs": true,
    "check_file_sizes": true,
    "min_file_size_kb": 1
  }
}'

# ── Logging Functions ──────────────────────────────────────────────────────
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" | tee -a "$LOG_FILE" >&2
}

log_success() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] SUCCESS: $1" | tee -a "$LOG_FILE"
}

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $1" | tee -a "$LOG_FILE"
}

# ── Reporting Helpers ──────────────────────────────────────────────────────
LIGAND_FAILURE_REPORT=""
LIGAND_SKIPPED_REPORT=""
RECEPTOR_FAILURE_REPORT=""
RECEPTOR_SKIPPED_REPORT=""

csv_escape() {
    local value="${1:-}"
    value="${value//$'\n'/ }"
    value="${value//$'\r'/ }"
    value="${value//\"/\"\"}"
    printf '%s' "$value"
}

init_report_file() {
    local report_file="$1"
    mkdir -p "$(dirname "$report_file")"
    printf '"stage","input_file","status","reason","details"\n' > "$report_file"
}

append_report_row() {
    local report_file="$1"
    local stage="$2"
    local input_file="$3"
    local status="$4"
    local reason="$5"
    local details="${6:-}"
    if [[ -z "$report_file" ]]; then
        return 0
    fi
    printf '"%s","%s","%s","%s","%s"\n' \
        "$(csv_escape "$stage")" \
        "$(csv_escape "$input_file")" \
        "$(csv_escape "$status")" \
        "$(csv_escape "$reason")" \
        "$(csv_escape "$details")" >> "$report_file"
}

pdb_content_metrics() {
    local pdb_file="$1"
    awk '
    BEGIN {
        atom=0; hetatm=0; non_water_hetatm=0;
    }
    {
        rec=substr($0,1,6); gsub(/ /,"",rec);
        chain=substr($0,22,1); gsub(/ /,"",chain);
        resseq=substr($0,23,4); gsub(/ /,"",resseq);
        resname=substr($0,18,3); gsub(/ /,"",resname);
        if (chain == "") chain="_";
        if (resseq == "") resseq="_";
        if (resname == "") resname="_";

        if (rec == "ATOM") {
            atom++;
            aa_res[chain ":" resseq]=1;
        } else if (rec == "HETATM") {
            hetatm++;
            het_res[chain ":" resseq ":" resname]=1;
            if (resname != "HOH" && resname != "WAT") {
                non_water_hetatm++;
            }
        }
    }
    END {
        aa_count=0; het_count=0;
        for (k in aa_res) aa_count++;
        for (k in het_res) het_count++;
        printf "%d|%d|%d|%d|%d\n", atom, hetatm, non_water_hetatm, aa_count, het_count;
    }' "$pdb_file"
}

classify_pdb_role() {
    local pdb_file="$1"
    local metrics
    metrics=$(pdb_content_metrics "$pdb_file")

    local atom_count=0
    local hetatm_count=0
    local non_water_hetatm_count=0
    local aa_residue_count=0
    local het_residue_count=0
    IFS='|' read -r atom_count hetatm_count non_water_hetatm_count aa_residue_count het_residue_count <<< "$metrics"

    if [[ $atom_count -ge 200 ]] || [[ $aa_residue_count -ge 30 ]]; then
        printf "receptor|%s" "$metrics"
        return 0
    fi

    if [[ $atom_count -eq 0 ]] && [[ $non_water_hetatm_count -gt 0 ]] && [[ $het_residue_count -le 8 ]]; then
        printf "ligand|%s" "$metrics"
        return 0
    fi

    if [[ $atom_count -le 40 ]] && [[ $aa_residue_count -le 3 ]] && [[ $non_water_hetatm_count -gt 0 ]] && [[ $het_residue_count -le 8 ]]; then
        printf "ligand|%s" "$metrics"
        return 0
    fi

    printf "unknown|%s" "$metrics"
}

# ── Progress Tracking ──────────────────────────────────────────────────────
show_progress() {
    local current=$1
    local total=$2
    local task=$3
    local percent=$((current * 100 / total))
    printf "\r[%3d%%] %s (%d/%d)" "$percent" "$task" "$current" "$total"
}

# ── Configuration Management ───────────────────────────────────────────────
load_config() {
    if [[ ! -f "$CONFIG_FILE" ]]; then
        log_info "Creating default configuration file: $CONFIG_FILE"
        echo "$DEFAULT_CONFIG" > "$CONFIG_FILE"
        log_info "Please edit $CONFIG_FILE and run again"
        exit 0
    fi
    
    # Load configuration using jq (install if not available)
    if ! command -v jq &> /dev/null; then
        log_error "jq is required for JSON configuration. Install with: brew install jq (macOS) or apt-get install jq (Ubuntu)"
        exit 1
    fi
    
    # Validate JSON
    if ! jq empty "$CONFIG_FILE" 2>/dev/null; then
        log_error "Invalid JSON configuration file: $CONFIG_FILE"
        exit 1
    fi
    
    log_info "Loaded configuration from: $CONFIG_FILE"
}

normalize_ligand_profile() {
    local raw="${1:-}"
    local token
    token=$(printf '%s' "$raw" | tr '[:upper:]' '[:lower:]')
    case "$token" in
        ""|auto|default) echo "engine_aware_full" ;;
        obabel|openbabel) echo "openbabel_only" ;;
        meeko) echo "meeko_only" ;;
        autodocktools|adt) echo "autodocktools_only" ;;
        openbabel_only|meeko_only|autodocktools_only|openbabel_meeko|openbabel_meeko_autodock|openbabel_autodocktools|engine_aware_full)
            echo "$token"
            ;;
        *)
            echo "engine_aware_full"
            ;;
    esac
}

profile_requires_autodocktools_ligand() {
    local profile="$1"
    local selected_engines_csv="${2:-}"
    case "$profile" in
        autodocktools_only|openbabel_autodocktools|openbabel_meeko_autodock)
            return 0
            ;;
        engine_aware_full)
            if [[ "$selected_engines_csv" == *"autodock4"* ]]; then
                return 0
            fi
            ;;
    esac
    return 1
}

derive_receptor_script_from_ligand_script() {
    local ligand_script="$1"
    if [[ -z "$ligand_script" ]]; then
        return 1
    fi
    local ligand_dir
    ligand_dir=$(dirname "$ligand_script")
    local receptor_candidate="$ligand_dir/prepare_receptor4.py"
    if [[ -f "$receptor_candidate" ]]; then
        printf '%s' "$receptor_candidate"
        return 0
    fi
    return 1
}

resolve_autodocktools_prepare_ligand4() {
    local configured_path
    configured_path=$(jq -r '.preparation.autodocktools_prepare_ligand4 // ""' "$CONFIG_FILE")
    if [[ -n "$configured_path" && -f "$configured_path" ]]; then
        printf '%s' "$configured_path"
        return 0
    fi

    for env_key in AUTODOCKTOOLS_PREPARE_LIGAND4 ADT_PREPARE_LIGAND4; do
        local env_value="${!env_key:-}"
        if [[ -n "$env_value" && -f "$env_value" ]]; then
            printf '%s' "$env_value"
            return 0
        fi
    done

    for command_name in prepare_ligand4.py prepare_ligand4; do
        if command -v "$command_name" >/dev/null 2>&1; then
            command -v "$command_name"
            return 0
        fi
    done

    local mgltools_root="${MGLTOOLS_PATH:-}"
    if [[ -n "$mgltools_root" ]]; then
        local candidate="$mgltools_root/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
        if [[ -f "$candidate" ]]; then
            printf '%s' "$candidate"
            return 0
        fi
    fi

    local candidates=(
        "$SCRIPT_DIR/.workflow/tools/autodocktools-prepare-py3k/AutoDockTools/Utilities24/prepare_ligand4.py"
        "$SCRIPT_DIR/tools/autodocktools-prepare-py3k/AutoDockTools/Utilities24/prepare_ligand4.py"
        "$HOME/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
        "/opt/mgltools/1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
        "/usr/local/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
    )
    local candidate=""
    for candidate in "${candidates[@]}"; do
        if [[ -f "$candidate" ]]; then
            printf '%s' "$candidate"
            return 0
        fi
    done

    return 1
}

resolve_autodocktools_prepare_receptor4() {
    local ligand_script_hint="${1:-}"
    local configured_path
    configured_path=$(jq -r '.preparation.autodocktools_prepare_receptor4 // ""' "$CONFIG_FILE")
    if [[ -n "$configured_path" && -f "$configured_path" ]]; then
        printf '%s' "$configured_path"
        return 0
    fi

    for env_key in AUTODOCKTOOLS_PREPARE_RECEPTOR4 ADT_PREPARE_RECEPTOR4; do
        local env_value="${!env_key:-}"
        if [[ -n "$env_value" && -f "$env_value" ]]; then
            printf '%s' "$env_value"
            return 0
        fi
    done

    if derive_receptor_script_from_ligand_script "$ligand_script_hint" >/dev/null 2>&1; then
        derive_receptor_script_from_ligand_script "$ligand_script_hint"
        return 0
    fi

    local configured_ligand_path
    configured_ligand_path=$(jq -r '.preparation.autodocktools_prepare_ligand4 // ""' "$CONFIG_FILE")
    if derive_receptor_script_from_ligand_script "$configured_ligand_path" >/dev/null 2>&1; then
        derive_receptor_script_from_ligand_script "$configured_ligand_path"
        return 0
    fi

    for env_key in AUTODOCKTOOLS_PREPARE_LIGAND4 ADT_PREPARE_LIGAND4; do
        local ligand_env="${!env_key:-}"
        if derive_receptor_script_from_ligand_script "$ligand_env" >/dev/null 2>&1; then
            derive_receptor_script_from_ligand_script "$ligand_env"
            return 0
        fi
    done

    for command_name in prepare_receptor4.py prepare_receptor4; do
        if command -v "$command_name" >/dev/null 2>&1; then
            command -v "$command_name"
            return 0
        fi
    done

    local mgltools_root="${MGLTOOLS_PATH:-}"
    if [[ -n "$mgltools_root" ]]; then
        local candidate="$mgltools_root/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
        if [[ -f "$candidate" ]]; then
            printf '%s' "$candidate"
            return 0
        fi
    fi

    local candidates=(
        "$SCRIPT_DIR/.workflow/tools/autodocktools-prepare-py3k/AutoDockTools/Utilities24/prepare_receptor4.py"
        "$SCRIPT_DIR/tools/autodocktools-prepare-py3k/AutoDockTools/Utilities24/prepare_receptor4.py"
        "$HOME/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
        "/opt/mgltools/1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
        "/usr/local/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
    )
    local candidate=""
    for candidate in "${candidates[@]}"; do
        if [[ -f "$candidate" ]]; then
            printf '%s' "$candidate"
            return 0
        fi
    done

    return 1
}

resolve_autodocktools_python() {
    local configured_python
    configured_python=$(jq -r '.preparation.autodocktools_python // ""' "$CONFIG_FILE")
    if [[ -n "$configured_python" ]]; then
        if [[ -x "$configured_python" ]]; then
            printf '%s' "$configured_python"
            return 0
        fi
        if command -v "$configured_python" >/dev/null 2>&1; then
            command -v "$configured_python"
            return 0
        fi
    fi

    for env_key in AUTODOCKTOOLS_PYTHON ADT_PYTHON; do
        local env_value="${!env_key:-}"
        if [[ -n "$env_value" ]]; then
            if [[ -x "$env_value" ]]; then
                printf '%s' "$env_value"
                return 0
            fi
            if command -v "$env_value" >/dev/null 2>&1; then
                command -v "$env_value"
                return 0
            fi
        fi
    done

    if command -v python3 >/dev/null 2>&1; then
        command -v python3
        return 0
    fi
    if command -v python >/dev/null 2>&1; then
        command -v python
        return 0
    fi
    return 1
}

autodocktools_pythonpath_root_from_script() {
    local script_path="$1"
    (
        local script_dir
        script_dir=$(cd "$(dirname "$script_path")" && pwd)
        cd "$script_dir/../.." && pwd
    )
}

prepare_receptor_with_autodocktools() {
    local input_file="$1"
    local output_file="$2"
    local adt_script="${AUTODOCKTOOLS_PREPARE_RECEPTOR4:-}"
    local adt_python="${AUTODOCKTOOLS_PYTHON:-python3}"

    if [[ -z "$adt_script" || ! -f "$adt_script" ]]; then
        log_error "AutoDockTools receptor script is not configured or missing: $adt_script"
        return 1
    fi

    local adt_root
    adt_root=$(autodocktools_pythonpath_root_from_script "$adt_script")
    local adt_output=""

    if ! adt_output=$(
        PYTHONPATH="$adt_root${PYTHONPATH:+:$PYTHONPATH}" \
        "$adt_python" "$adt_script" \
            -r "$input_file" \
            -o "$output_file" \
            -A checkhydrogens \
            -U nphs_lps_waters_deleteAltB 2>&1
    ); then
        log_error "AutoDockTools receptor preparation failed for $(basename "$input_file")"
        if [[ -n "$adt_output" ]]; then
            local first_error_line=""
            first_error_line=$(printf '%s\n' "$adt_output" | tail -n 1)
            log_error "prepare_receptor4.py output: $first_error_line"
        fi
        return 1
    fi

    if [[ -n "$adt_output" ]]; then
        local last_line=""
        last_line=$(printf '%s\n' "$adt_output" | tail -n 1)
        if [[ -n "$last_line" ]]; then
            log_info "prepare_receptor4.py: $last_line"
        fi
    fi
    return 0
}

# ── Dependency Checking ────────────────────────────────────────────────────
check_dependencies() {
    local missing_deps=()
    local ligand_profile_raw
    ligand_profile_raw=$(jq -r '.preparation.ligand_preparation_profile // .preparation.ligand_preparation_backend // "engine_aware_full"' "$CONFIG_FILE")
    local ligand_profile
    ligand_profile=$(normalize_ligand_profile "$ligand_profile_raw")
    local selected_engines_csv
    selected_engines_csv=$(jq -r '.preparation.selected_engines // [] | join(",")' "$CONFIG_FILE")
    local adt_ligand_script=""
    local adt_receptor_script=""
    local adt_python=""
    
    # Check strictly required tools
    local required_tools=("obabel" "jq")
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_deps+=("$tool")
        fi
    done

    # Check optional preparation tools and log their status.
    if ! command -v pdb2pqr30 &> /dev/null; then
        log_info "pdb2pqr30 not found. Receptor preparation will fall back to Open Babel."
    fi
    if ! mk_prepare_ligand.py --help >/dev/null 2>&1; then
        log_info "mk_prepare_ligand.py is not usable. Ligand preparation will rely on profile-specific fallbacks."
    fi
    if ! mk_prepare_receptor.py --help >/dev/null 2>&1; then
        log_info "mk_prepare_receptor.py is not usable. Receptor preparation will fall back to Open Babel."
    fi

    if profile_requires_autodocktools_ligand "$ligand_profile" "$selected_engines_csv"; then
        if adt_ligand_script=$(resolve_autodocktools_prepare_ligand4); then
            log_info "AutoDockTools ligand script resolved: $adt_ligand_script"
        else
            missing_deps+=("prepare_ligand4.py")
        fi
        if adt_python=$(resolve_autodocktools_python); then
            log_info "AutoDockTools python resolved: $adt_python"
        else
            missing_deps+=("python3")
        fi
        if [[ "$ligand_profile" == "autodocktools_only" ]]; then
            if adt_receptor_script=$(resolve_autodocktools_prepare_receptor4 "$adt_ligand_script"); then
                log_info "AutoDockTools receptor script resolved: $adt_receptor_script"
            else
                missing_deps+=("prepare_receptor4.py")
            fi
        fi
    fi
    
    if [[ ${#missing_deps[@]} -gt 0 ]]; then
        log_error "Missing required dependencies: ${missing_deps[*]}"
        log_error "Install AutoDockTools and OpenBabel to continue"
        exit 1
    fi

    if [[ -n "$adt_ligand_script" ]]; then
        export AUTODOCKTOOLS_PREPARE_LIGAND4="$adt_ligand_script"
    fi
    if [[ -n "$adt_receptor_script" ]]; then
        export AUTODOCKTOOLS_PREPARE_RECEPTOR4="$adt_receptor_script"
    fi
    if [[ -n "$adt_python" ]]; then
        export AUTODOCKTOOLS_PYTHON="$adt_python"
    fi
    if [[ "$ligand_profile" == "autodocktools_only" ]]; then
        log_info "AutoDockTools-only profile active: receptor prep uses prepare_receptor4.py."
    fi
    
    log_success "All required dependencies found"
}

has_working_meeko_ligand() {
    mk_prepare_ligand.py --help >/dev/null 2>&1
}

has_working_meeko_receptor() {
    mk_prepare_receptor.py --help >/dev/null 2>&1
}

prepare_ligand_to_pdbqt() {
    local input_file="$1"
    local output_file="$2"
    local profile_raw
    profile_raw=$(jq -r '.preparation.ligand_preparation_profile // .preparation.ligand_preparation_backend // "engine_aware_full"' "$CONFIG_FILE")
    local profile
    profile=$(normalize_ligand_profile "$profile_raw")
    local selected_engines
    selected_engines=$(jq -r '.preparation.selected_engines // [] | join(",")' "$CONFIG_FILE")
    local ph_value
    ph_value=$(jq -r '.preparation.ph // 7.4' "$CONFIG_FILE")

    if command -v python3 >/dev/null 2>&1; then
        if PYTHONPATH="$SCRIPT_DIR${PYTHONPATH:+:$PYTHONPATH}" \
            python3 -m docking.preparation.ligand_preparation \
                --input "$input_file" \
                --output "$output_file" \
                --backend-profile "$profile" \
                --ph "$ph_value" \
                --engines "$selected_engines" >/dev/null 2>&1; then
            return 0
        fi
        if [[ "$profile" != "engine_aware_full" ]]; then
            log_error "Central ligand preparation failed in forced profile mode ($profile) for $(basename "$input_file")."
            return 1
        fi
        log_info "Central ligand normalization failed for $(basename "$input_file"). Falling back to direct CLI conversion."
    fi

    if [[ "$profile" != "engine_aware_full" ]]; then
        log_error "python3-based central ligand preparation is required for profile '$profile'."
        return 1
    fi

    local fallback_input="$input_file"
    local tmp_workdir=""
    if command -v obabel >/dev/null 2>&1; then
        tmp_workdir=$(mktemp -d "${TMPDIR:-/tmp}/pdbwiz_ligprep_fallback_XXXXXX")
        local source_sdf="$tmp_workdir/source.sdf"
        local with_3d_sdf="$tmp_workdir/with_3d.sdf"
        local with_h_sdf="$tmp_workdir/with_h.sdf"
        local minimized_sdf="$tmp_workdir/minimized.sdf"
        if obabel "$input_file" -O "$source_sdf" >/dev/null 2>&1 \
            && obabel "$source_sdf" -O "$with_3d_sdf" --gen3d >/dev/null 2>&1 \
            && obabel "$with_3d_sdf" -O "$with_h_sdf" -h >/dev/null 2>&1; then
            fallback_input="$with_h_sdf"
            if obabel "$with_h_sdf" -O "$minimized_sdf" --minimize --steps 250 --ff MMFF94 >/dev/null 2>&1 \
                || obabel "$with_h_sdf" -O "$minimized_sdf" --minimize --steps 250 --ff UFF >/dev/null 2>&1; then
                fallback_input="$minimized_sdf"
            fi
            log_info "Fallback ligand normalization (Open Babel 3D) succeeded for $(basename "$input_file")."
        else
            log_info "Fallback ligand normalization failed for $(basename "$input_file"). Using original input for direct conversion."
        fi
    fi

    local rc=1
    if has_working_meeko_ligand; then
        if mk_prepare_ligand.py -i "$fallback_input" -o "$output_file" 2>/dev/null; then
            rc=0
        fi
        if [[ $rc -ne 0 ]]; then
            log_info "Meeko ligand preparation failed for $(basename "$input_file"). Falling back to Open Babel."
        fi
    fi

    if [[ $rc -ne 0 ]]; then
        if obabel "$fallback_input" -O "$output_file" -h >/dev/null 2>&1; then
            rc=0
        fi
    fi

    if [[ -n "$tmp_workdir" && -d "$tmp_workdir" ]]; then
        rm -rf "$tmp_workdir"
    fi
    return $rc
}

# ── File Format Conversion ─────────────────────────────────────────────────
convert_ligand_format() {
    local input_file="$1"
    local output_format="$2"
    local output_file="$3"
    
    case "$output_format" in
        "sdf")
            obabel "$input_file" -O "$output_file" 2>/dev/null || return 1
            ;;
        "mol2")
            obabel "$input_file" -O "$output_file" 2>/dev/null || return 1
            ;;
        "pdb")
            obabel "$input_file" -O "$output_file" 2>/dev/null || return 1
            ;;
        *)
            log_error "Unsupported output format: $output_format"
            return 1
            ;;
    esac
}

# ── Quality Control ────────────────────────────────────────────────────────
validate_file() {
    local file="$1"
    local min_size_kb="${2:-1}"
    
    if [[ ! -f "$file" ]]; then
        log_error "File not found: $file"
        return 1
    fi
    
    local file_size_kb=$(($(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null) / 1024))
    if [[ $file_size_kb -lt $min_size_kb ]]; then
        log_error "File too small: $file (${file_size_kb}KB < ${min_size_kb}KB)"
        return 1
    fi
    
    # Check if PDBQT file has proper format
    if [[ "$file" == *.pdbqt ]]; then
        if ! grep -q "ATOM\|HETATM" "$file"; then
            log_error "Invalid PDBQT file: $file (no ATOM/HETATM records)"
            return 1
        fi
    fi
    
    return 0
}

validate_ligand_pdbqt() {
    local file="$1"

    if ! validate_file "$file"; then
        return 1
    fi

    local root_count
    local torsdof_count
    root_count=$(grep -c '^ROOT$' "$file" 2>/dev/null || echo 0)
    torsdof_count=$(grep -c '^TORSDOF' "$file" 2>/dev/null || echo 0)

    if [[ "$root_count" -ne 1 ]]; then
        log_error "Invalid ligand PDBQT: $file (expected 1 ROOT block, found $root_count)"
        return 1
    fi
    if [[ "$torsdof_count" -ne 1 ]]; then
        log_error "Invalid ligand PDBQT: $file (expected 1 TORSDOF record, found $torsdof_count)"
        return 1
    fi

    return 0
}

# ── Ligand Preparation ─────────────────────────────────────────────────────
prepare_ligands() {
    local input_dir="$1"
    local output_dir="$2"
    local formats=("${@:3}")
    
    log_info "Starting ligand preparation..."
    mkdir -p "$output_dir"
    
    local candidate_files=()
    local skipped_receptor_like=0
    local processed_files=0
    
    # Collect ligand candidates while excluding receptor-like PDB artifacts that
    # may be present in the raw ligand directory.
    for format in "${formats[@]}"; do
        while IFS= read -r -d '' mol_file; do
            local format_lower
            format_lower=$(printf '%s' "$format" | tr '[:upper:]' '[:lower:]')
            if [[ "$format_lower" == "pdb" ]]; then
                local classification
                classification=$(classify_pdb_role "$mol_file")
                local role="${classification%%|*}"
                local metrics="${classification#*|}"
                if [[ "$role" == "receptor" ]]; then
                    skipped_receptor_like=$((skipped_receptor_like + 1))
                    log_info "Skipping receptor-like PDB in ligand input: $(basename "$mol_file")"
                    append_report_row "$LIGAND_SKIPPED_REPORT" "ligand" "$mol_file" "skipped" "receptor_like_pdb" "metrics=$metrics"
                    continue
                fi
                if [[ "$role" == "unknown" ]]; then
                    log_info "Ligand PDB role uncertain for $(basename "$mol_file"); proceeding with preparation (metrics=$metrics)."
                fi
            fi
            candidate_files+=("$mol_file")
        done < <(find "$input_dir" -iname "*.${format}" -type f -print0 2>/dev/null)
    done

    local total_files=${#candidate_files[@]}
    
    if [[ $total_files -eq 0 ]]; then
        if [[ $skipped_receptor_like -gt 0 ]]; then
            log_info "No valid ligand files found after filtering. Excluded $skipped_receptor_like receptor-like file(s) from: $input_dir"
        else
            log_info "No ligand files found in: $input_dir"
        fi
        return 0
    fi
    
    log_info "Found $total_files ligand files to process"
    if [[ $skipped_receptor_like -gt 0 ]]; then
        log_info "Excluded $skipped_receptor_like receptor-like file(s) from ligand preparation"
    fi
    
    # Process filtered ligand candidates
    for mol_file in "${candidate_files[@]}"; do
            processed_files=$((processed_files + 1))
            show_progress "$processed_files" "$total_files" "Processing ligands"
            
            local base_name
            base_name=$(basename "$mol_file")
            base_name="${base_name%.*}"
            local output_file="$output_dir/${base_name}.pdbqt"
            local started_epoch
            started_epoch=$(date +%s)
            log_info "[$processed_files/$total_files] Preparing ligand: $(basename "$mol_file")"
            
            # Skip if already exists and valid
            if [[ -f "$output_file" ]] && validate_file "$output_file"; then
                log_info "Skipping existing valid ligand output: $(basename "$output_file")"
                append_report_row "$LIGAND_SKIPPED_REPORT" "ligand" "$mol_file" "skipped" "already_prepared" "$output_file"
                continue
            fi
            
            if prepare_ligand_to_pdbqt "$mol_file" "$output_file"; then
                if validate_ligand_pdbqt "$output_file"; then
                    local elapsed
                    elapsed=$(( $(date +%s) - started_epoch ))
                    log_success "Prepared ligand: $base_name (${elapsed}s)"
                else
                    log_error "Invalid output file: $output_file"
                    append_report_row "$LIGAND_FAILURE_REPORT" "ligand" "$mol_file" "failed" "invalid_pdbqt_output" "$output_file"
                    rm -f "$output_file"
                fi
            else
                log_error "Failed to prepare ligand: $mol_file"
                append_report_row "$LIGAND_FAILURE_REPORT" "ligand" "$mol_file" "failed" "backend_conversion_failed" ""
            fi

    done
    
    echo # New line after progress
    local success_count=$(find "$output_dir" -name "*.pdbqt" -type f | wc -l)
    log_success "Ligand preparation completed: $success_count/$total_files files"
}

# ── Receptor Preparation ───────────────────────────────────────────────────
prepare_receptors() {
    local input_dir="$1"
    local output_dir="$2"
    local formats=("${@:3}")
    
    log_info "Starting receptor preparation..."
    mkdir -p "$output_dir"
    
    local candidate_files=()
    local processed_files=0
    local skipped_ligand_like=0
    local ligand_profile_raw
    ligand_profile_raw=$(jq -r '.preparation.ligand_preparation_profile // .preparation.ligand_preparation_backend // "engine_aware_full"' "$CONFIG_FILE")
    local ligand_profile
    ligand_profile=$(normalize_ligand_profile "$ligand_profile_raw")
    local receptor_use_pdb2pqr
    receptor_use_pdb2pqr=$(jq -r '.preparation.receptor_use_pdb2pqr // false' "$CONFIG_FILE")

    # Collect receptor candidates while excluding ligand-like PDB artifacts that
    # may be present in the raw protein directory.
    for format in "${formats[@]}"; do
        while IFS= read -r -d '' pdb_file; do
            local file_name_lower
            file_name_lower=$(basename "$pdb_file" | tr '[:upper:]' '[:lower:]')
            if [[ "$file_name_lower" == *"_ligand_"* ]] || [[ "$file_name_lower" == ligand_* ]] || [[ "$file_name_lower" == *"_ligand."* ]]; then
                skipped_ligand_like=$((skipped_ligand_like + 1))
                log_info "Skipping ligand-like file in receptor input: $(basename "$pdb_file")"
                append_report_row "$RECEPTOR_SKIPPED_REPORT" "receptor" "$pdb_file" "skipped" "ligand_like_filename" ""
                continue
            fi
            local classification
            classification=$(classify_pdb_role "$pdb_file")
            local role="${classification%%|*}"
            local metrics="${classification#*|}"
            if [[ "$role" == "ligand" ]]; then
                skipped_ligand_like=$((skipped_ligand_like + 1))
                log_info "Skipping ligand-like PDB by content in receptor input: $(basename "$pdb_file")"
                append_report_row "$RECEPTOR_SKIPPED_REPORT" "receptor" "$pdb_file" "skipped" "ligand_like_content" "metrics=$metrics"
                continue
            fi
            if [[ "$role" == "unknown" ]]; then
                log_info "Receptor PDB role uncertain for $(basename "$pdb_file"); proceeding with preparation (metrics=$metrics)."
            fi
            candidate_files+=("$pdb_file")
        done < <(find "$input_dir" -name "*.${format}" -type f -print0 2>/dev/null)
    done

    local total_files=${#candidate_files[@]}
    
    if [[ $total_files -eq 0 ]]; then
        if [[ $skipped_ligand_like -gt 0 ]]; then
            log_info "No valid receptor files found after filtering. Excluded $skipped_ligand_like ligand-like file(s) from: $input_dir"
        else
            log_info "No receptor files found in: $input_dir"
        fi
        return 0
    fi
    
    log_info "Found $total_files receptor files to process"
    if [[ $skipped_ligand_like -gt 0 ]]; then
        log_info "Excluded $skipped_ligand_like ligand-like file(s) from receptor preparation"
    fi
    
    # Process filtered receptor candidates
    for pdb_file in "${candidate_files[@]}"; do
            processed_files=$((processed_files + 1))
            show_progress "$processed_files" "$total_files" "Processing receptors"
            
            local base_name
            base_name=$(basename "$pdb_file")
            base_name="${base_name%.*}"
            local pqr_file="$output_dir/${base_name}.pqr"
            local clean_pdb="$output_dir/${base_name}_clean.pdb"
            local pdbqt_file="$output_dir/${base_name}.pdbqt"
            local started_epoch
            started_epoch=$(date +%s)
            log_info "[$processed_files/$total_files] Preparing receptor: $(basename "$pdb_file")"
            
            # Skip if already exists and valid
            if [[ -f "$pdbqt_file" ]] && validate_file "$pdbqt_file"; then
                log_info "Skipping existing valid receptor output: $(basename "$pdbqt_file")"
                append_report_row "$RECEPTOR_SKIPPED_REPORT" "receptor" "$pdb_file" "skipped" "already_prepared" "$pdbqt_file"
                continue
            fi
            
            local prep_input="$pdb_file"
            local used_clean_intermediate="false"
            local prepared_ok="false"

            if [[ "$ligand_profile" == "autodocktools_only" ]]; then
                if prepare_receptor_with_autodocktools "$pdb_file" "$pdbqt_file"; then
                    prepared_ok="true"
                fi
            else
                # Optional PDB2PQR path (disabled by default for better geometry stability,
                # especially in protein+nucleic-acid assemblies).
                if [[ "$receptor_use_pdb2pqr" == "true" ]]; then
                    if command -v pdb2pqr30 >/dev/null 2>&1; then
                        if pdb2pqr30 --ff "$(jq -r '.preparation.force_field' "$CONFIG_FILE")" \
                                    --with-ph "$(jq -r '.preparation.ph' "$CONFIG_FILE")" \
                                    "$pdb_file" "$pqr_file" >/dev/null 2>&1; then
                            if obabel "$pqr_file" -O "$clean_pdb" >/dev/null 2>&1; then
                                prep_input="$clean_pdb"
                                used_clean_intermediate="true"
                            else
                                log_info "Open Babel clean step failed for $(basename "$pdb_file"). Falling back to the original PDB."
                            fi
                        else
                            log_info "PDB2PQR failed for $(basename "$pdb_file"). Falling back to the original PDB."
                        fi
                    fi
                fi

                # Step 2: PDB → PDBQT via Meeko when usable, else Open Babel.
                if has_working_meeko_receptor; then
                    local meeko_args=("--read_pdb" "$prep_input" "-p" "$pdbqt_file")
                    
                    if [[ "$(jq -r '.preparation.allow_bad_res' "$CONFIG_FILE")" == "true" ]]; then
                        meeko_args+=("--allow_bad_res")
                    fi
                    
                    meeko_args+=("--default_altloc" "$(jq -r '.preparation.default_altloc' "$CONFIG_FILE")")
                    
                    if mk_prepare_receptor.py "${meeko_args[@]}" 2>/dev/null; then
                        prepared_ok="true"
                    else
                        log_info "Meeko receptor preparation failed for $(basename "$pdb_file"). Falling back to Open Babel."
                    fi
                fi

                if [[ "$prepared_ok" != "true" ]]; then
                    if obabel "$prep_input" -O "$pdbqt_file" -xr >/dev/null 2>&1; then
                        prepared_ok="true"
                    fi
                fi
            fi

            if [[ "$prepared_ok" == "true" ]]; then
                if validate_file "$pdbqt_file"; then
                    local elapsed
                    elapsed=$(( $(date +%s) - started_epoch ))
                    log_success "Prepared receptor: $base_name (${elapsed}s)"
                    
                else
                    log_error "Invalid output file: $pdbqt_file"
                    append_report_row "$RECEPTOR_FAILURE_REPORT" "receptor" "$pdb_file" "failed" "invalid_pdbqt_output" "$pdbqt_file"
                    rm -f "$pdbqt_file"
                fi
            else
                log_error "Failed to prepare receptor: $pdb_file"
                append_report_row "$RECEPTOR_FAILURE_REPORT" "receptor" "$pdb_file" "failed" "backend_conversion_failed" ""
            fi
            
            # Clean up intermediate files
            if [[ "$used_clean_intermediate" == "true" ]]; then
                rm -f "$pqr_file" "$clean_pdb"
            else
                rm -f "$pqr_file"
            fi
    done
    
    echo # New line after progress
    local success_count=$(find "$output_dir" -name "*.pdbqt" -type f | wc -l)
    log_success "Receptor preparation completed: $success_count/$total_files files"
}

# ── Main Execution ─────────────────────────────────────────────────────────
main() {
    log_info "Starting Enhanced AutoDock Preparation Script"
    log_info "Configuration file: $CONFIG_FILE"
    log_info "Log file: $LOG_FILE"
    
    # Load configuration
    load_config
    
    # Check dependencies
    check_dependencies
    
    # Get configuration values
    local ligands_input=$(jq -r '.input.ligands.path' "$CONFIG_FILE")
    local receptors_input=$(jq -r '.input.receptors.path' "$CONFIG_FILE")
    local ligands_output=$(jq -r '.output.ligands' "$CONFIG_FILE")
    local receptors_output=$(jq -r '.output.receptors' "$CONFIG_FILE")
    local ligand_profile_raw
    ligand_profile_raw=$(jq -r '.preparation.ligand_preparation_profile // .preparation.ligand_preparation_backend // "engine_aware_full"' "$CONFIG_FILE")
    local ligand_profile
    ligand_profile=$(normalize_ligand_profile "$ligand_profile_raw")
    local selected_engines
    selected_engines=$(jq -r '.preparation.selected_engines // [] | join(",")' "$CONFIG_FILE")
    local adt_script=$(jq -r '.preparation.autodocktools_prepare_ligand4 // ""' "$CONFIG_FILE")
    local adt_receptor_script=$(jq -r '.preparation.autodocktools_prepare_receptor4 // ""' "$CONFIG_FILE")
    local adt_python=$(jq -r '.preparation.autodocktools_python // ""' "$CONFIG_FILE")

    export PDBWIZARD_LIGAND_PREP_BACKEND="$ligand_profile"
    export PDBWIZARD_LIGAND_PREP_PROFILE="$ligand_profile"
    export PDBWIZARD_LIGAND_PREP_PH="$(jq -r '.preparation.ph // 7.4' "$CONFIG_FILE")"
    export PDBWIZARD_SELECTED_ENGINES="$selected_engines"
    if [[ -n "$adt_script" && -f "$adt_script" ]]; then
        export AUTODOCKTOOLS_PREPARE_LIGAND4="$adt_script"
    fi
    if [[ -n "$adt_receptor_script" && -f "$adt_receptor_script" ]]; then
        export AUTODOCKTOOLS_PREPARE_RECEPTOR4="$adt_receptor_script"
    fi
    if [[ -n "$adt_python" ]]; then
        export AUTODOCKTOOLS_PYTHON="$adt_python"
    fi
    
    # Handle same folder scenario
    local in_same_folder=$(jq -r '.input.ligands.in_same_folder' "$CONFIG_FILE")
    if [[ "$in_same_folder" == "true" ]]; then
        ligands_input="$receptors_input"
        log_info "Using same folder for ligands and receptors: $ligands_input"
    fi
    
    # Create output directories
    mkdir -p "$ligands_output" "$receptors_output" "$(dirname "$LOG_FILE")"
    local ligand_step_report_dir="$ligands_output/preparation_steps"
    mkdir -p "$ligand_step_report_dir"
    export PDBWIZARD_LIGAND_PREP_REPORT_DIR="$ligand_step_report_dir"

    LIGAND_FAILURE_REPORT="$ligands_output/ligand_preparation_failures.csv"
    LIGAND_SKIPPED_REPORT="$ligands_output/ligand_preparation_skipped.csv"
    RECEPTOR_FAILURE_REPORT="$receptors_output/receptor_preparation_failures.csv"
    RECEPTOR_SKIPPED_REPORT="$receptors_output/receptor_preparation_skipped.csv"
    init_report_file "$LIGAND_FAILURE_REPORT"
    init_report_file "$LIGAND_SKIPPED_REPORT"
    init_report_file "$RECEPTOR_FAILURE_REPORT"
    init_report_file "$RECEPTOR_SKIPPED_REPORT"
    
    # Prepare ligands
    local ligand_formats=($(jq -r '.input.ligands.formats[]' "$CONFIG_FILE"))
    prepare_ligands "$ligands_input" "$ligands_output" "${ligand_formats[@]}"
    
    # Prepare receptors
    local receptor_formats=($(jq -r '.input.receptors.formats[]' "$CONFIG_FILE"))
    prepare_receptors "$receptors_input" "$receptors_output" "${receptor_formats[@]}"
    
    # Final summary
    local total_ligands=$(find "$ligands_output" -name "*.pdbqt" -type f 2>/dev/null | wc -l)
    local total_receptors=$(find "$receptors_output" -name "*.pdbqt" -type f 2>/dev/null | wc -l)
    
    log_success "Preparation completed successfully!"
    log_info "Ligands prepared: $total_ligands"
    log_info "Receptors prepared: $total_receptors"
    log_info "Log file saved: $LOG_FILE"
    log_info "Ligand failure report: $LIGAND_FAILURE_REPORT"
    log_info "Ligand skipped report: $LIGAND_SKIPPED_REPORT"
    log_info "Ligand step reports: $ligand_step_report_dir"
    log_info "Receptor failure report: $RECEPTOR_FAILURE_REPORT"
    log_info "Receptor skipped report: $RECEPTOR_SKIPPED_REPORT"
    
    # Generate summary report
    local summary_file="$receptors_output/preparation_summary.txt"
    cat > "$summary_file" << EOF
AutoDock Preparation Summary
============================
Date: $(date)
Configuration: $CONFIG_FILE
Log File: $LOG_FILE

Results:
--------
Ligands prepared: $total_ligands
Receptors prepared: $total_receptors

Output Directories:
------------------
Ligands: $ligands_output
Receptors: $receptors_output

Reports:
--------
Ligand failures: $LIGAND_FAILURE_REPORT
Ligand skipped: $LIGAND_SKIPPED_REPORT
Ligand step reports: $ligand_step_report_dir
Receptor failures: $RECEPTOR_FAILURE_REPORT
Receptor skipped: $RECEPTOR_SKIPPED_REPORT

Ligand preparation profile: $ligand_profile
Selected engines: ${selected_engines:-none}
Force Field: $(jq -r '.preparation.force_field' "$CONFIG_FILE")
pH: $(jq -r '.preparation.ph' "$CONFIG_FILE")
EOF
    
    log_info "Summary report saved: $summary_file"
}

# ── Error Handling ─────────────────────────────────────────────────────────
trap 'log_error "Script interrupted at line $LINENO"' INT TERM
trap 'log_error "Script failed at line $LINENO"' ERR

# ── Execute Main Function ──────────────────────────────────────────────────
main "$@"
