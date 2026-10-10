from __future__ import annotations

import csv
import json
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from .models import ProjectManifest


# Single source of the pairlist header (pairlist_builder, pair_curation and project_builder import it).
PAIRLIST_COLUMNS = [
    "receptor",
    "site_id",
    "ligand",
    "center_x",
    "center_y",
    "center_z",
    "size_x",
    "size_y",
    "size_z",
    "protein_display_name",
    "ligand_display_name",
    "pdb_id",
    "pair_source",
    "is_cocrystal_benchmark",
    "cocrystal_ligand_name",
    "cocrystal_ligand_display_name",
    "selection_mode",
    # Spec 036 R5a box provenance. Written by every pairlist writer; readers tolerate old pairlists
    # that lack these columns (missing values load as empty or None).
    "box_method",
    "ligand_rg_angstrom",
    "edge_angstrom",
    "box_containment_status",
    "box_warnings",
    # Spec 037 R2b: where the ligand Rg came from (ligand_preparation or computed_at_box_step).
    "rg_source",
    # Spec 037 R1c pair marks (written by the pairlist builder; marks, not blocks).
    "receptor_policy_level",
    "ligand_policy_level",
    "receptor_ph",
    "ligand_ph",
    "pair_protonation_status",
]

LAYOUT_CANONICAL = "canonical"
LAYOUT_DOCKING_LEGACY = "docking_legacy"
LAYOUT_PROFILES = {LAYOUT_CANONICAL, LAYOUT_DOCKING_LEGACY}

DOCKING_LEGACY_DIRS = {
    "raw_ligands": "1-Raw_Ligand",
    # Compatibility alias: normalized SDF copies now live beside the source
    # ligands in the single raw-ligands staging folder.
    "raw_ligands_sdf": "1-Raw_Ligand",
    "raw_proteins": "2-Raw_Protien",
    "preparation": "3-Preparation",
    "prepared_proteins": "3-Preparation/2-Prepared_Protiens",
    "prepared_ligands": "3-Preparation/4-Prepared_Ligand",
    "preparation_logs": "3-Preparation/preparation_logs",
    "docking_root": "4-Docking",
    "receptors": "4-Docking/receptors",
    "ligands": "4-Docking/ligands",
    "metadata": "4-Docking/metadata",
    "post_docking_root": "5-Analysis",
    "analysis": "5-Analysis",
    "raw_data": "5-Analysis/raw_data",
    "reports": "7-Reports",
    "visualizations": "5-Analysis/visualizations",
    "interactions": "5-Analysis/interactions",
    "rmsd_analysis": "5-Analysis/rmsd_analysis",
    "complexes": "5-Analysis/complexes",
    "best_poses": "5-Analysis/best_poses",
}

GNINA_HPC_COMPAT_DIRS = [
    "results",
    "scripts",
    "scripts_hpc",
]

NUMBERED_OUTPUT_TOPOLOGY = {
    "meta": ".meta",
    "input_root": "0-Input",
    "input_receptors": "0-Input/receptors",
    "input_ligands": "0-Input/ligands",
    "input_references": "0-Input/references",
    "preparation_root": "1-Preparation",
    "prep_receptors_prepared": "1-Preparation/receptors/prepared",
    "prep_receptors_qc": "1-Preparation/receptors/qc",
    "prep_ligands_prepared": "1-Preparation/ligands/prepared",
    "prep_ligands_flagged": "1-Preparation/ligands/flagged",
    "gridboxes_root": "2-GridBoxes",
    "grid_definitions": "2-GridBoxes/definitions",
    "grid_validation": "2-GridBoxes/validation",
    "docking_root_numbered": "3-Docking",
    "docking_gnina_configs": "3-Docking/gnina/configs",
    "docking_gnina_poses": "3-Docking/gnina/poses",
    "docking_gnina_logs": "3-Docking/gnina/logs",
    "docking_vina_configs": "3-Docking/vina/configs",
    "docking_vina_poses": "3-Docking/vina/poses",
    "docking_vina_logs": "3-Docking/vina/logs",
    "docking_smina_configs": "3-Docking/smina/configs",
    "docking_smina_poses": "3-Docking/smina/poses",
    "docking_smina_logs": "3-Docking/smina/logs",
    "docking_autodock4_configs": "3-Docking/autodock4/configs",
    "docking_autodock4_poses": "3-Docking/autodock4/poses",
    "docking_autodock4_logs": "3-Docking/autodock4/logs",
    "post_docking_root_numbered": "4-Working",
    "post_metadata": "4-Working/metadata",
    "post_scores_raw": "4-Working/scores/raw",
    "post_scores_unified": "4-Working/scores/unified",
    "post_scores_consensus": "4-Working/scores/consensus",
    "post_filtering": "4-Working/filtering",
    "post_filter_reports": "4-Working/filtering/filter_reports",
    "post_rescoring_candidates": "4-Working/rescoring",
    "post_rescoring_results": "4-Working/rescoring/rescored",
    "post_rerun_manifests": "4-Working/rerun_manifests",
    "post_storage": "4-Working/storage",
    "analysis_root_numbered": "5-Analysis",
    "analysis_sessions_root": "5-Analysis/sessions",
    "analysis_comparative": "5-Analysis/comparative",
    "analysis_raw_data": "5-Analysis/raw_data",
    "analysis_rmsd": "5-Analysis/rmsd",
    "analysis_polypharmacology": "5-Analysis/polypharmacology",
    "analysis_structure_quality": "5-Analysis/structure_quality",
    "analysis_top_pose": "5-Analysis/top_pose_ligand_performance",
    "visualizations_root_numbered": "6-Visualizations",
    "visualizations_2d": "6-Visualizations/2d",
    "visualizations_2d_heatmaps": "6-Visualizations/2d/heatmaps",
    "visualizations_2d_score_plots": "6-Visualizations/2d/score_plots",
    "visualizations_2d_interaction_maps": "6-Visualizations/2d/interaction_maps",
    "visualizations_3d": "6-Visualizations/3d",
    "visualizations_3d_py3dmol": "6-Visualizations/3d/py3dmol",
    "visualizations_3d_pymol": "6-Visualizations/3d/pymol",
    "reports_root_numbered": "7-Reports",
    "reports_per_protein": "7-Reports/per_protein",
    "reports_figures": "7-Reports/figures",
    "reports_tables": "7-Reports/tables",
    "sessions_root_numbered": "sessions",
}


def _recorded_layout_profile(project_root: Path) -> Optional[str]:
    """
    Return the layout profile recorded inside an existing project manifest.

    Reads the manifest files directly (not through manifest_path) so that this
    helper can be used by detect_layout_profile without recursion. A manifest
    that sits where its own recorded profile would place it (4-Docking/ for
    docking_legacy, the project root for canonical) takes precedence. Otherwise
    the first recorded profile found (checking 4-Docking/ before the root) is
    used, so a stray manifest left by an earlier run cannot override the one
    written at the recorded location.
    """
    root = Path(project_root).expanduser().resolve()
    locations = {
        LAYOUT_DOCKING_LEGACY: root / DOCKING_LEGACY_DIRS["docking_root"] / "project_manifest.json",
        LAYOUT_CANONICAL: root / "project_manifest.json",
    }
    recorded: Dict[Path, str] = {}
    for location in locations.values():
        if not location.is_file():
            continue
        try:
            with open(location, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
        except (OSError, ValueError):
            continue
        value = payload.get("layout_profile") if isinstance(payload, dict) else None
        if value in LAYOUT_PROFILES:
            recorded[location] = value
    for profile, location in locations.items():
        if recorded.get(location) == profile:
            return profile
    return next(iter(recorded.values()), None)


def detect_layout_profile(project_root: Path, requested: Optional[str] = None) -> str:
    if requested:
        if requested not in LAYOUT_PROFILES:
            raise ValueError(f"Unsupported layout profile: {requested}")
        return requested

    root = Path(project_root).expanduser().resolve()
    recorded = _recorded_layout_profile(root)
    if recorded:
        return recorded
    if (root / DOCKING_LEGACY_DIRS["docking_root"]).exists():
        return LAYOUT_DOCKING_LEGACY
    return LAYOUT_CANONICAL


def resolve_layout_profile(project_root: Path) -> str:
    """
    Layout profile that workflow steps must use for a project.

    The profile recorded in the project manifest wins. Projects without a
    recorded manifest keep the legacy default that workflow steps have always
    used for them.
    """
    return _recorded_layout_profile(Path(project_root)) or LAYOUT_DOCKING_LEGACY


def docking_root(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    if profile == LAYOUT_DOCKING_LEGACY:
        return root / DOCKING_LEGACY_DIRS["docking_root"]
    return root


def post_docking_root(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    if profile == LAYOUT_DOCKING_LEGACY:
        return root / DOCKING_LEGACY_DIRS["post_docking_root"]
    return root / "5-Analysis"


def _canonical_layout(project_root: Path) -> Dict[str, Path]:
    return {
        "project_root": project_root,
        "raw_ligands": project_root / "ligands_raw",
        "raw_ligands_sdf": project_root / "ligands_raw",
        "raw_proteins": project_root / "proteins_raw",
        "preparation": project_root / "preparation",
        "prepared_proteins": project_root / "prepared_proteins",
        "prepared_ligands": project_root / "prepared_ligands",
        "preparation_logs": project_root / "preparation_logs",
        "docking_root": project_root,
        "receptors": project_root / "receptors",
        "ligands": project_root / "ligands",
        "metadata": project_root / "metadata",
        "post_docking_root": project_root / "5-Analysis",
        "analysis": project_root / "5-Analysis",
        "raw_data": project_root / "5-Analysis" / "raw_data",
        "reports": project_root / "7-Reports",
        "visualizations": project_root / "5-Analysis" / "visualizations",
        "interactions": project_root / "5-Analysis" / "interactions",
        "rmsd_analysis": project_root / "5-Analysis" / "rmsd_analysis",
        "complexes": project_root / "5-Analysis" / "complexes",
        "best_poses": project_root / "5-Analysis" / "best_poses",
    }


def _legacy_layout(project_root: Path) -> Dict[str, Path]:
    layout = {"project_root": project_root}
    for key, rel_path in DOCKING_LEGACY_DIRS.items():
        layout[key] = project_root / rel_path
    return layout


def project_layout_paths(project_root: Path, layout_profile: Optional[str] = None) -> Dict[str, Path]:
    """Layout dict of a project under its own root, without creating any directory.

    Use this to resolve where a stage's files live (Spec 036 R6). Use ``ensure_project_layout``
    when the directories must exist.
    """
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    layout = _legacy_layout(root) if profile == LAYOUT_DOCKING_LEGACY else _canonical_layout(root)
    layout["layout_profile"] = profile
    return layout


def ensure_project_layout(project_root: Path, layout_profile: str = LAYOUT_CANONICAL) -> Dict[str, Path]:
    project_root = Path(project_root).expanduser().resolve()
    layout = project_layout_paths(project_root, layout_profile)
    profile = layout["layout_profile"]
    for key, path in layout.items():
        if key in {"project_root", "layout_profile"}:
            continue
        path.mkdir(parents=True, exist_ok=True)
    ensure_post_docking_compat_shim(project_root)
    return layout


def ensure_engine_layout(project_root: Path, engine: str, layout_profile: Optional[str] = None) -> Dict[str, Path]:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    dock_root = docking_root(root, profile)

    if profile == LAYOUT_DOCKING_LEGACY:
        if engine == "gnina":
            layout = {
                "root": dock_root,
                "poses": dock_root / "gnina_out",
                "logs": dock_root / "logs",
                "scores": dock_root / "results",
            }
            for path in layout.values():
                path.mkdir(parents=True, exist_ok=True)
            for dirname in GNINA_HPC_COMPAT_DIRS:
                (dock_root / dirname).mkdir(parents=True, exist_ok=True)
            return layout

        engine_root = dock_root / f"{engine}_out"
        layout = {
            "root": engine_root,
            "poses": engine_root / "poses",
            "logs": engine_root / "logs",
            "scores": engine_root / "scores",
        }
        for path in layout.values():
            path.mkdir(parents=True, exist_ok=True)
        return layout

    engine_root = dock_root / "engines" / engine
    layout = {
        "root": engine_root,
        "poses": engine_root / "poses",
        "logs": engine_root / "logs",
        "scores": engine_root / "scores",
    }
    for path in layout.values():
        path.mkdir(parents=True, exist_ok=True)
    return layout


def ensure_gnina_hpc_compat_layout(project_root: Path, layout_profile: Optional[str] = None) -> Dict[str, str]:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    gnina_layout = ensure_engine_layout(root, "gnina", layout_profile=profile)
    dock_root = docking_root(root, profile)

    created_dirs: Dict[str, str] = {
        "gnina_out": str(gnina_layout["poses"]),
        "logs": str(gnina_layout["logs"]),
    }
    for dirname in GNINA_HPC_COMPAT_DIRS:
        path = dock_root / dirname
        path.mkdir(parents=True, exist_ok=True)
        created_dirs[dirname] = str(path)
    return created_dirs


def ensure_numbered_output_layout(project_root: Path, layout_profile: Optional[str] = None) -> Dict[str, Path]:
    """
    Ensure canonical numbered DockForge output folders exist.

    Numbered DockForge folders are always created at project root. This keeps
    human-facing analysis outputs in `5-Analysis/` and intermediate scratch data
    in `4-Working/`, regardless of whether the project started from a legacy dock
    layout.
    """
    root = Path(project_root).expanduser().resolve()
    detect_layout_profile(root, layout_profile)
    numbered_base = root

    layout: Dict[str, Path] = {}
    for key, relative in NUMBERED_OUTPUT_TOPOLOGY.items():
        path = numbered_base / relative
        path.mkdir(parents=True, exist_ok=True)
        layout[key] = path
    layout["numbered_base"] = numbered_base
    ensure_post_docking_compat_shim(root)
    return layout


def ensure_post_docking_compat_shim(project_root: Path) -> Dict[str, object]:
    root = Path(project_root).expanduser().resolve()
    canonical = root / "5-Analysis"
    legacy = root / "5-Post_Docking_Analysis"
    canonical.mkdir(parents=True, exist_ok=True)

    status = "created"
    if legacy.exists() or legacy.is_symlink():
        if legacy.is_symlink():
            status = "existing_symlink"
        else:
            status = "legacy_dir_exists"
            warnings.warn(
                f"Legacy post-docking directory exists at {legacy}. New runs write to {canonical}.",
                RuntimeWarning,
                stacklevel=2,
            )
    else:
        try:
            legacy.symlink_to(canonical, target_is_directory=True)
            warnings.warn(
                f"Created compatibility symlink {legacy} -> {canonical}. Use {canonical} for new analysis outputs.",
                RuntimeWarning,
                stacklevel=2,
            )
        except (OSError, NotImplementedError) as exc:
            try:
                legacy.mkdir(parents=True, exist_ok=True)
                status = "fallback_dir_created"
            except Exception:
                status = "symlink_skipped_unsupported"
            warnings.warn(
                f"Could not create symlink {legacy} -> {canonical} ({exc}). Using directory fallback.",
                RuntimeWarning,
                stacklevel=2,
            )
    return {
        "status": status,
        "canonical_root": canonical,
        "legacy_root": legacy,
    }


def pairlist_path(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    return docking_root(project_root, layout_profile) / "pairlist.csv"


def pair_intent_path(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    return docking_root(project_root, layout_profile) / "metadata" / "pair_intent.csv"


def pair_curation_state_path(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    return docking_root(project_root, layout_profile) / "metadata" / "pair_curation_state.json"


def manifest_path(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    return docking_root(project_root, layout_profile) / "project_manifest.json"


def deployment_root(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    return docking_root(project_root, layout_profile) / "deployments"


def shared_receptors_dir(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    return ensure_project_layout(project_root, detect_layout_profile(project_root, layout_profile))["receptors"]


def shared_ligands_dir(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    return ensure_project_layout(project_root, detect_layout_profile(project_root, layout_profile))["ligands"]


def load_pairlist(project_root: Path, layout_profile: Optional[str] = None) -> pd.DataFrame:
    return pd.read_csv(pairlist_path(project_root, layout_profile))


def load_manifest(project_root: Path, layout_profile: Optional[str] = None) -> Dict[str, object]:
    root = Path(project_root).expanduser().resolve()
    primary = manifest_path(root, layout_profile)
    candidates = [primary]
    root_candidate = root / "project_manifest.json"
    legacy_candidate = root / DOCKING_LEGACY_DIRS["docking_root"] / "project_manifest.json"

    for candidate in (root_candidate, legacy_candidate):
        if candidate not in candidates:
            candidates.append(candidate)

    for candidate in candidates:
        if not candidate.exists():
            continue
        with open(candidate, "r", encoding="utf-8") as handle:
            return json.load(handle)

    raise FileNotFoundError(f"No project manifest found. Checked: {', '.join(str(path) for path in candidates)}")


def save_manifest(project_root: Path, payload: Dict[str, object], layout_profile: Optional[str] = None) -> Path:
    path = manifest_path(project_root, layout_profile)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    return path


# ---------------------------------------------------------------- Spec 036 R1a: one protonation policy per project

PROTONATION_POLICY_KEY = "protonation_policy"
PROTONATION_LIGAND_POLICIES = ("ph_model", "explicit_state", "as_input")
PROTONATION_CONFLICT_REASON = "protonation_policy_conflict"


class ProtonationPolicyConflict(ValueError):
    """An explicit value disagrees with the project's stored protonation policy (Spec 036 R1a)."""

    def __init__(self, field: str, stored: object, requested: object) -> None:
        self.field = str(field)
        self.stored = stored
        self.requested = requested
        self.reason = PROTONATION_CONFLICT_REASON
        super().__init__(
            f"{PROTONATION_CONFLICT_REASON}:{self.field}: project policy has {stored!r}, "
            f"this command requested {requested!r}. Use the stored value, or change the project policy "
            "with `workflow protonation-policy --replace`."
        )


def _manifest_file_for_update(project_root: Path) -> Optional[Path]:
    """The manifest file that load_manifest would read (so a write goes back to the same file)."""
    root = Path(project_root).expanduser().resolve()
    for candidate in (
        manifest_path(root),
        root / "project_manifest.json",
        root / DOCKING_LEGACY_DIRS["docking_root"] / "project_manifest.json",
    ):
        if candidate.is_file():
            return candidate
    return None


def load_protonation_policy(project_root: Path) -> Optional[Dict[str, object]]:
    """The project's stored ``protonation_policy`` block, or None when the project has none."""
    root = Path(project_root).expanduser().resolve()
    if _manifest_file_for_update(root) is None:
        return None
    block = load_manifest(root).get(PROTONATION_POLICY_KEY)
    return dict(block) if isinstance(block, dict) and block else None


def _policy_path_text(project_root: Path, path: Path) -> str:
    resolved = Path(path).expanduser().resolve()
    try:
        return resolved.relative_to(Path(project_root).expanduser().resolve()).as_posix()
    except ValueError:
        return str(resolved)


def _file_sha256(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def protonation_state_map_path(project_root: Path, block: Optional[Dict[str, object]]) -> Optional[Path]:
    """Resolve the stored ligand state map and verify its SHA-256 (fails closed when it changed)."""
    if not block or not block.get("ligand_state_map"):
        return None
    record = block["ligand_state_map"]
    raw = Path(str(record.get("path") or "")).expanduser()
    candidate = raw if raw.is_absolute() else (Path(project_root).expanduser().resolve() / raw)
    if not candidate.is_file():
        raise ValueError(f"protonation_state_map_missing:{candidate}")
    if _file_sha256(candidate) != str(record.get("sha256") or ""):
        raise ValueError(f"protonation_state_map_hash_mismatch:{candidate}")
    return candidate.resolve()


def build_protonation_policy_block(
    project_root: Path,
    *,
    receptor_ph: float,
    receptor_force_field: str,
    ligand_policy: str,
    ligand_ph: Optional[float] = None,
    ligand_state_map: Optional[Path] = None,
    source: str = "user_entered",
    set_at: Optional[str] = None,
) -> Dict[str, object]:
    """Validate explicit protonation values and return the manifest block (no file is written)."""
    from datetime import datetime, timezone

    from .models import validate_preparation_ph
    from .preparation.receptor_preparation import PDB2PQR_FORCE_FIELDS

    ok, normalized_receptor_ph, error = validate_preparation_ph(receptor_ph)
    if not ok:
        raise ValueError(f"invalid_receptor_ph:{error}")
    force_field = str(receptor_force_field or "").strip().upper()
    if force_field not in PDB2PQR_FORCE_FIELDS:
        raise ValueError(f"invalid_receptor_force_field:{receptor_force_field!r}")
    policy = str(ligand_policy or "").strip().lower()
    if policy not in PROTONATION_LIGAND_POLICIES:
        raise ValueError(f"invalid_ligand_policy:{ligand_policy!r}")
    normalized_ligand_ph: Optional[float] = None
    if ligand_ph is not None:
        ok, normalized_ligand_ph, error = validate_preparation_ph(ligand_ph)
        if not ok:
            raise ValueError(f"invalid_ligand_ph:{error}")
    if policy == "ph_model" and normalized_ligand_ph is None:
        raise ValueError("ligand_ph_required_for_ph_model")
    state_map_record = None
    if policy == "explicit_state":
        if ligand_state_map is None:
            raise ValueError("ligand_state_map_required_for_explicit_state")
        state_map = Path(ligand_state_map).expanduser().resolve()
        if not state_map.is_file():
            raise ValueError(f"ligand_state_map_missing:{state_map}")
        state_map_record = {"path": _policy_path_text(project_root, state_map), "sha256": _file_sha256(state_map)}
    elif ligand_state_map is not None:
        raise ValueError(f"ligand_state_map_only_for_explicit_state:{policy}")
    return {
        "receptor_ph": float(normalized_receptor_ph),
        "receptor_force_field": force_field,
        "ligand_policy": policy,
        "ligand_ph": float(normalized_ligand_ph) if normalized_ligand_ph is not None else None,
        "ligand_state_map": state_map_record,
        "source": str(source or "user_entered"),
        "set_at": set_at or datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
    }


def check_protonation_policy_conflicts(
    stored: Optional[Dict[str, object]],
    requested: Dict[str, object],
) -> List[Tuple[str, object, object]]:
    """Explicit values that disagree with the stored policy. ``None`` requested values are not checked."""
    if not stored:
        return []
    conflicts: List[Tuple[str, object, object]] = []
    for field, value in requested.items():
        if value is None:
            continue
        if field == "ligand_state_map":
            record = stored.get("ligand_state_map") or {}
            stored_value: object = record.get("sha256") if isinstance(record, dict) else None
            requested_value: object = _file_sha256(Path(str(value)).expanduser()) if Path(str(value)).expanduser().is_file() else value
            same = stored_value == requested_value
        elif field in {"receptor_ph", "ligand_ph"}:
            stored_value = stored.get(field)
            requested_value = float(value)
            same = stored_value is not None and abs(float(stored_value) - requested_value) <= 1e-9
        elif field == "receptor_force_field":
            stored_value = stored.get(field)
            requested_value = str(value).strip().upper()
            same = str(stored_value or "").upper() == requested_value
        else:
            stored_value = stored.get(field)
            requested_value = str(value).strip().lower()
            same = str(stored_value or "").lower() == requested_value
        if not same:
            conflicts.append((field, stored_value, requested_value))
    return conflicts


def require_matching_protonation_policy(stored: Optional[Dict[str, object]], requested: Dict[str, object]) -> None:
    conflicts = check_protonation_policy_conflicts(stored, requested)
    if conflicts:
        field, stored_value, requested_value = conflicts[0]
        raise ProtonationPolicyConflict(field, stored_value, requested_value)


def resolve_protonation_inputs(
    project_root: Path,
    *,
    receptor_ph: Optional[float] = None,
    receptor_force_field: Optional[str] = None,
    ligand_policy: Optional[str] = None,
    ligand_ph: Optional[float] = None,
    ligand_state_map: Optional[str] = None,
) -> Dict[str, object]:
    """Explicit values, falling back to the project policy (R1a). Conflicts raise ``ProtonationPolicyConflict``."""
    stored = load_protonation_policy(project_root)
    require_matching_protonation_policy(
        stored,
        {
            "receptor_ph": receptor_ph,
            "receptor_force_field": receptor_force_field,
            "ligand_policy": ligand_policy,
            "ligand_ph": ligand_ph,
            "ligand_state_map": ligand_state_map,
        },
    )
    resolved: Dict[str, object] = {
        "project_policy": stored is not None,
        "receptor_ph": receptor_ph,
        "receptor_ph_source": "user_entered" if receptor_ph is not None else None,
        "receptor_force_field": receptor_force_field,
        "force_field_source": "user_entered" if receptor_force_field else None,
        "ligand_policy": ligand_policy,
        "ligand_ph": ligand_ph,
        "ligand_ph_source": "user_entered" if ligand_ph is not None else None,
        "ligand_state_map": ligand_state_map,
    }
    if stored:
        source = str(stored.get("source") or "user_entered")
        if resolved["receptor_ph"] is None and stored.get("receptor_ph") is not None:
            resolved["receptor_ph"] = float(stored["receptor_ph"])
            resolved["receptor_ph_source"] = source
        if resolved["receptor_force_field"] is None and stored.get("receptor_force_field"):
            resolved["receptor_force_field"] = str(stored["receptor_force_field"])
            resolved["force_field_source"] = source
        if resolved["ligand_policy"] is None and stored.get("ligand_policy"):
            resolved["ligand_policy"] = str(stored["ligand_policy"])
        if resolved["ligand_ph"] is None and stored.get("ligand_ph") is not None:
            resolved["ligand_ph"] = float(stored["ligand_ph"])
            resolved["ligand_ph_source"] = source
        if resolved["ligand_state_map"] is None:
            state_map = protonation_state_map_path(project_root, stored)
            if state_map is not None:
                resolved["ligand_state_map"] = str(state_map)
    return resolved


def set_protonation_policy(
    project_root: Path,
    *,
    receptor_ph: float,
    receptor_force_field: str,
    ligand_policy: str,
    ligand_ph: Optional[float] = None,
    ligand_state_map: Optional[Path] = None,
    source: str = "user_entered",
    replace: bool = False,
) -> Dict[str, object]:
    """Store the project protonation policy. A different stored policy needs ``replace=True``."""
    root = Path(project_root).expanduser().resolve()
    block = build_protonation_policy_block(
        root,
        receptor_ph=receptor_ph,
        receptor_force_field=receptor_force_field,
        ligand_policy=ligand_policy,
        ligand_ph=ligand_ph,
        ligand_state_map=ligand_state_map,
        source=source,
    )
    stored = load_protonation_policy(root)
    if stored and not replace:
        comparable = {k: block[k] for k in ("receptor_ph", "receptor_force_field", "ligand_policy", "ligand_ph")}
        conflicts = check_protonation_policy_conflicts(stored, comparable)
        if block.get("ligand_state_map"):
            conflicts += check_protonation_policy_conflicts(stored, {"ligand_state_map": str(ligand_state_map)})
        if conflicts:
            field, stored_value, requested_value = conflicts[0]
            raise ProtonationPolicyConflict(field, stored_value, requested_value)
        return stored
    file_path = _manifest_file_for_update(root)
    if file_path is None:
        raise FileNotFoundError(f"No project manifest found under {root}; run workflow init first.")
    payload = json.loads(file_path.read_text(encoding="utf-8"))
    previous = payload.get(PROTONATION_POLICY_KEY) or {}
    # Spec 037 R1a: replacing the project default keeps the per-entity overrides.
    for section in (PROTONATION_RECEPTOR_OVERRIDES, PROTONATION_LIGAND_OVERRIDES):
        if isinstance(previous.get(section), dict) and previous.get(section):
            block[section] = dict(previous[section])
    payload[PROTONATION_POLICY_KEY] = block
    file_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return block


# ---------------------------------------------------------------- Spec 037 R1: layered protonation policy
#
# Three levels. The project default is the flat Spec 036 block (``receptor_ph``, ``receptor_force_field``,
# ``ligand_policy``, ``ligand_ph``, ``ligand_state_map``). ``receptors.<key>`` and ``ligands.<key>`` are
# per-entity overrides that sit above it. ``<key>`` is the file stem of the receptor or ligand.

PROTONATION_LEVEL_PROJECT = "project"
PROTONATION_LEVEL_RECEPTOR = "receptor"
PROTONATION_LEVEL_LIGAND = "ligand"
PROTONATION_LEVEL_NONE = "none"
PROTONATION_RECEPTOR_OVERRIDES = "receptors"
PROTONATION_LIGAND_OVERRIDES = "ligands"
POLICY_HASH_VERSION = "protonation_effective_v1"
STALE_PREPARATION_REASON = "stale_preparation_policy_changed"
UNRECORDED_PREPARATION_REASON = "preparation_policy_unrecorded"
PAIR_STATUS_COHERENT = "coherent"
PAIR_STATUS_PH_MISMATCH = "pair_ph_mismatch"
PAIR_STATUS_EXPLICIT_STATE = "explicit_state_no_ph"
PAIR_STATUS_AS_INPUT = "as_input_no_ph"
PAIR_STATUS_RECEPTOR_PH_MISSING = "receptor_ph_missing"
PAIR_STATUS_NO_POLICY = "no_stored_policy"


class ProtonationPolicyError(ValueError):
    """A layered protonation value cannot be resolved, or a prepared entity is stale (Spec 037 R1)."""

    def __init__(self, reason: str, details: str = "") -> None:
        self.reason = str(reason)
        self.details = str(details)
        super().__init__(f"{self.reason}: {self.details}" if self.details else self.reason)


def entity_policy_key(name: object) -> str:
    """Override key of a receptor or ligand: its file stem (1IEP.pdbqt, 1IEP.pdb and 1IEP all key as 1IEP)."""
    return Path(str(name or "").strip()).stem


def _normalize_net_charge(value: object) -> Optional[object]:
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.lower() == "none":
        return None
    try:
        return int(text)
    except ValueError:
        return text


def protonation_policy_hash(effective: Dict[str, object]) -> str:
    """Content hash of the effective values that determine the prepared state (Spec 037 R1d)."""
    import hashlib

    payload = json.dumps({"version": POLICY_HASH_VERSION, "effective": effective}, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _state_map_rows(path: Path) -> Dict[str, Dict[str, object]]:
    """Raw state-map rows keyed by entity key. No chemistry here: the ligand module validates the rows."""
    rows: Dict[str, Dict[str, object]] = {}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            fields = {str(key or "").strip().lower(): str(value or "").strip() for key, value in raw.items()}
            ligand = fields.get("ligand", "")
            if ligand:
                rows.setdefault(
                    entity_policy_key(ligand),
                    {"smiles": fields.get("smiles", ""), "net_charge": _normalize_net_charge(fields.get("net_charge", ""))},
                )
    return rows


def resolve_entity_protonation(project_root: Path, kind: str, name: object) -> Dict[str, object]:
    """Effective protonation policy of one receptor or ligand and the level it came from (Spec 037 R1a).

    Precedence: ``receptors.<key>`` or ``ligands.<key>`` over the project default. A field an override does
    not set is inherited from the project default and marked ``project`` in ``field_levels``. Without any
    stored policy the level is ``none`` and no hash is recorded. A missing explicit state raises
    ``ProtonationPolicyError`` (fails closed).
    """
    if kind not in {"receptor", "ligand"}:
        raise ValueError(f"unknown protonation entity kind: {kind!r}")
    root = Path(project_root).expanduser().resolve()
    key = entity_policy_key(name)
    result: Dict[str, object] = {
        "entity_kind": kind,
        "entity": key,
        "policy_level": PROTONATION_LEVEL_NONE,
        "field_levels": {},
        "effective": None,
        "policy_hash": None,
        "policy_hash_version": POLICY_HASH_VERSION,
    }
    block = load_protonation_policy(root)
    if not block:
        return result
    if kind == "receptor":
        override = dict((block.get(PROTONATION_RECEPTOR_OVERRIDES) or {}).get(key) or {})
        ph = block.get("receptor_ph")
        force_field = block.get("receptor_force_field")
        levels = {"receptor_ph": PROTONATION_LEVEL_PROJECT, "receptor_force_field": PROTONATION_LEVEL_PROJECT}
        if override.get("receptor_ph") is not None:
            ph = override["receptor_ph"]
            levels["receptor_ph"] = PROTONATION_LEVEL_RECEPTOR
        if override.get("receptor_force_field"):
            force_field = override["receptor_force_field"]
            levels["receptor_force_field"] = PROTONATION_LEVEL_RECEPTOR
        effective: Dict[str, object] = {
            "receptor_ph": float(ph) if ph is not None else None,
            "receptor_force_field": str(force_field).strip().upper() if force_field else None,
        }
        result["policy_level"] = (
            PROTONATION_LEVEL_RECEPTOR if PROTONATION_LEVEL_RECEPTOR in levels.values() else PROTONATION_LEVEL_PROJECT
        )
    else:
        override = dict((block.get(PROTONATION_LIGAND_OVERRIDES) or {}).get(key) or {})
        policy = str(block.get("ligand_policy") or "")
        ph = block.get("ligand_ph")
        levels = {"ligand_policy": PROTONATION_LEVEL_PROJECT, "ligand_ph": PROTONATION_LEVEL_PROJECT, "explicit_state": PROTONATION_LEVEL_PROJECT}
        if override.get("ligand_policy"):
            policy = str(override["ligand_policy"])
            levels["ligand_policy"] = PROTONATION_LEVEL_LIGAND
        if override.get("ligand_ph") is not None:
            ph = override["ligand_ph"]
            levels["ligand_ph"] = PROTONATION_LEVEL_LIGAND
        explicit: Optional[Dict[str, object]] = None
        if policy == "explicit_state":
            has_inline = bool(str(override.get("smiles") or "").strip()) or _normalize_net_charge(override.get("net_charge")) is not None
            if has_inline:
                explicit = {
                    "smiles": str(override.get("smiles") or "").strip(),
                    "net_charge": _normalize_net_charge(override.get("net_charge")),
                }
                levels["explicit_state"] = PROTONATION_LEVEL_LIGAND
            else:
                if block.get("ligand_policy") != "explicit_state" or not block.get("ligand_state_map"):
                    raise ProtonationPolicyError(
                        "explicit_state_missing_entry",
                        f"{key}: no explicit state (ligands.{key} SMILES or net_charge, or a project explicit_state map)",
                    )
                state_map = protonation_state_map_path(root, block)
                row = _state_map_rows(state_map).get(key) if state_map is not None else None
                if row is None:
                    raise ProtonationPolicyError("explicit_state_missing_entry", f"{key} has no row in the project state map")
                if not row["smiles"] and row["net_charge"] is None:
                    raise ProtonationPolicyError("explicit_state_missing_entry", f"{key}: the state map row has neither smiles nor net_charge")
                explicit = {"smiles": str(row["smiles"]), "net_charge": row["net_charge"]}
        effective = {"ligand_policy": policy, "ligand_ph": None, "explicit_state": None}
        if policy == "ph_model":
            if ph is None:
                raise ProtonationPolicyError("ligand_ph_required_for_ph_model", key)
            effective["ligand_ph"] = float(ph)
        if policy == "explicit_state":
            effective["explicit_state"] = explicit
        result["ligand_policy"] = policy
        result["ligand_ph"] = effective["ligand_ph"]
        result["explicit_state"] = explicit
        result["policy_level"] = (
            PROTONATION_LEVEL_LIGAND if PROTONATION_LEVEL_LIGAND in levels.values() else PROTONATION_LEVEL_PROJECT
        )
    result["field_levels"] = levels
    if kind == "receptor":
        result["receptor_ph"] = effective["receptor_ph"]
        result["receptor_force_field"] = effective["receptor_force_field"]
    result["effective"] = effective
    result["policy_hash"] = protonation_policy_hash(effective)
    return result


def check_preparation_policy_current(provenance: Optional[Dict[str, object]], resolution: Dict[str, object]) -> None:
    """Refuse a prepared entity whose provenance does not record the current policy hash (Spec 037 R1d).

    ``provenance`` None means no record was found and nothing is checked here. A record without a
    ``policy_hash`` key was written before Spec 037; it is refused when a project policy is stored.
    """
    if provenance is None:
        return
    entity = f"{resolution.get('entity_kind')}:{resolution.get('entity')}"
    current = resolution.get("policy_hash")
    if "policy_hash" not in provenance:
        if current is None:
            return
        raise ProtonationPolicyError(
            UNRECORDED_PREPARATION_REASON,
            f"{entity}: the preparation record has no policy hash; re-prepare under the stored policy",
        )
    recorded = provenance.get("policy_hash")
    if recorded != current:
        raise ProtonationPolicyError(
            STALE_PREPARATION_REASON,
            f"{entity}: prepared under {recorded or 'no policy'}, now resolves to {current or 'no policy'}"
            f" ({resolution.get('policy_level')}); re-prepare this entity",
        )


def pair_protonation_columns(receptor: Dict[str, object], ligand: Dict[str, object]) -> Dict[str, object]:
    """Spec 037 R1c pair marks from the resolved receptor and ligand policies. Marks, not blocks."""
    ligand_policy = ligand.get("ligand_policy")
    receptor_ph = receptor.get("receptor_ph")
    ligand_ph = ligand.get("ligand_ph")
    if ligand.get("policy_level") == PROTONATION_LEVEL_NONE:
        status = PAIR_STATUS_NO_POLICY
    elif ligand_policy == "explicit_state":
        status = PAIR_STATUS_EXPLICIT_STATE
    elif ligand_policy == "as_input":
        status = PAIR_STATUS_AS_INPUT
    elif receptor_ph is None:
        status = PAIR_STATUS_RECEPTOR_PH_MISSING
    elif abs(float(receptor_ph) - float(ligand_ph)) <= 1e-9:
        status = PAIR_STATUS_COHERENT
    else:
        status = PAIR_STATUS_PH_MISMATCH
    return {
        "receptor_policy_level": receptor.get("policy_level"),
        "ligand_policy_level": ligand.get("policy_level"),
        "receptor_ph": float(receptor_ph) if receptor_ph is not None else None,
        "ligand_ph": float(ligand_ph) if ligand_policy == "ph_model" and ligand_ph is not None else None,
        "pair_protonation_status": status,
    }


def set_protonation_override(
    project_root: Path,
    *,
    kind: str,
    name: object,
    replace: bool = False,
    receptor_ph: Optional[float] = None,
    receptor_force_field: Optional[str] = None,
    ligand_policy: Optional[str] = None,
    ligand_ph: Optional[float] = None,
    smiles: Optional[str] = None,
    net_charge: Optional[int] = None,
    ligand_state_map: Optional[Path] = None,
    source: str = "user_entered",
) -> Dict[str, object]:
    """Store a receptor or ligand override under ``receptors.<key>`` or ``ligands.<key>`` (Spec 037 R1a).

    The project default must be stored first. A value that differs from the stored override is refused
    unless ``replace`` is true, as for the project default (Spec 036 R1a). A ligand state map CSV gives
    the row for this ligand, which is stored inline together with the map's SHA-256.
    """
    from datetime import datetime, timezone

    from .models import validate_preparation_ph
    from .preparation.receptor_preparation import PDB2PQR_FORCE_FIELDS

    root = Path(project_root).expanduser().resolve()
    block = load_protonation_policy(root)
    if not block:
        raise ProtonationPolicyError(
            "protonation_default_missing",
            "store the project default first (workflow protonation-policy without --receptor or --ligand)",
        )
    key = entity_policy_key(name)
    if not key:
        raise ValueError("protonation_override_name_required")
    provided: Dict[str, object] = {}
    if kind == "receptor":
        if receptor_ph is None and receptor_force_field is None:
            raise ValueError("receptor_override_needs_ph_or_force_field")
        if receptor_ph is not None:
            ok, normalized, error = validate_preparation_ph(receptor_ph)
            if not ok:
                raise ValueError(f"invalid_receptor_ph:{error}")
            provided["receptor_ph"] = float(normalized)
        if receptor_force_field is not None:
            force_field = str(receptor_force_field).strip().upper()
            if force_field not in PDB2PQR_FORCE_FIELDS:
                raise ValueError(f"invalid_receptor_force_field:{receptor_force_field!r}")
            provided["receptor_force_field"] = force_field
        section = PROTONATION_RECEPTOR_OVERRIDES
    elif kind == "ligand":
        if ligand_policy is not None:
            policy = str(ligand_policy).strip().lower()
            if policy not in PROTONATION_LIGAND_POLICIES:
                raise ValueError(f"invalid_ligand_policy:{ligand_policy!r}")
            provided["ligand_policy"] = policy
        if ligand_ph is not None:
            ok, normalized, error = validate_preparation_ph(ligand_ph)
            if not ok:
                raise ValueError(f"invalid_ligand_ph:{error}")
            provided["ligand_ph"] = float(normalized)
        if smiles is not None and str(smiles).strip():
            from docking.preparation.ligand_preparation import canonical_smiles_from_text

            text = str(smiles).strip()
            try:
                canonical_smiles_from_text(text)
            except Exception as exc:
                raise ValueError(f"invalid_ligand_smiles:{exc}") from exc
            provided["smiles"] = text
        if net_charge is not None:
            provided["net_charge"] = _normalize_net_charge(net_charge)
        if provided.get("smiles") and provided.get("net_charge") is not None:
            from docking.preparation.ligand_preparation import ProtonationStateError, _inline_state_entry

            try:
                _inline_state_entry(Path(key), {"smiles": provided["smiles"], "net_charge": provided["net_charge"]})
            except ProtonationStateError as exc:
                raise ValueError(f"invalid_ligand_explicit_state:{exc.reason}:{exc.details}") from exc
        if ligand_state_map is not None:
            from docking.preparation.ligand_preparation import load_protonation_state_map

            state_map = Path(ligand_state_map).expanduser().resolve()
            entries = load_protonation_state_map(state_map)
            match = [entry for ligand_text, entry in entries.items() if entity_policy_key(ligand_text) == key]
            if not match:
                raise ProtonationPolicyError("explicit_state_missing_entry", f"{key} has no row in {state_map.name}")
            provided["smiles"] = str(match[0].get("smiles") or "")
            provided["net_charge"] = _normalize_net_charge(match[0].get("net_charge"))
            provided["state_map"] = {"path": _policy_path_text(root, state_map), "sha256": _file_sha256(state_map)}
        if not provided:
            raise ValueError("ligand_override_needs_a_value")
        if provided.get("ligand_policy") == "explicit_state" and not (provided.get("smiles") or provided.get("net_charge") is not None):
            raise ProtonationPolicyError(
                "ligand_explicit_state_missing", f"{key}: explicit_state needs --ligand-smiles, --ligand-net-charge or --ligand-state-map"
            )
        effective_policy = str(provided.get("ligand_policy") or block.get("ligand_policy") or "")
        state_fields = {"smiles", "net_charge", "state_map"} & set(provided)
        if state_fields and effective_policy != "explicit_state":
            # A state given here would be silently ignored under ph_model or as_input: refuse instead.
            raise ProtonationPolicyError(
                "explicit_state_fields_need_explicit_policy",
                f"{key}: {', '.join(sorted(state_fields))} given but the effective ligand policy is {effective_policy!r}; "
                "add --ligand-policy explicit_state",
            )
        if effective_policy == "ph_model" and provided.get("ligand_ph") is None and block.get("ligand_ph") is None:
            raise ProtonationPolicyError("ligand_ph_required_for_ph_model", key)
        section = PROTONATION_LIGAND_OVERRIDES
    else:
        raise ValueError(f"unknown protonation entity kind: {kind!r}")

    overrides = dict(block.get(section) or {})
    stored = dict(overrides.get(key) or {})
    if stored and not replace:
        for field, value in provided.items():
            if field == "state_map" or stored.get(field) is None:
                continue
            old, new = stored[field], value
            same = abs(float(old) - float(new)) <= 1e-9 if isinstance(new, float) and isinstance(old, (int, float)) else old == new
            if not same:
                raise ProtonationPolicyConflict(f"{section}.{key}.{field}", old, new)
    record = {**stored, **provided, "source": str(source or "user_entered"), "set_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")}
    overrides[key] = record
    file_path = _manifest_file_for_update(root)
    if file_path is None:
        raise FileNotFoundError(f"No project manifest found under {root}; run workflow init first.")
    payload = json.loads(file_path.read_text(encoding="utf-8"))
    saved = dict(payload.get(PROTONATION_POLICY_KEY) or {})
    saved[section] = overrides
    payload[PROTONATION_POLICY_KEY] = saved
    file_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return record


def ensure_pairlist_stub(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    path = pairlist_path(project_root, layout_profile)
    if path.exists():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(columns=PAIRLIST_COLUMNS).to_csv(path, index=False)
    return path


def bootstrap_project_layout(
    project_root: Path,
    engines: List[str],
    project_name: str = "",
    favorite_engine: str = "",
    asset_mode: str = "symlink",
    source_paths: Optional[Dict[str, str]] = None,
    notes: Optional[List[str]] = None,
    layout_profile: str = LAYOUT_CANONICAL,
    pairlist_mode: str = "",
    has_cocrystal_benchmark_rows: bool = False,
) -> Dict[str, object]:
    project_root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(project_root, layout_profile)
    layout = ensure_project_layout(project_root, profile)
    pairlist_file = ensure_pairlist_stub(project_root, profile)
    existing = load_manifest(project_root, profile) if manifest_path(project_root, profile).exists() else {}

    for engine in engines:
        ensure_engine_layout(project_root, engine, layout_profile=profile)

    compatibility_profiles: List[str] = list(existing.get("compatibility_profiles", []))
    if "gnina" in engines:
        ensure_gnina_hpc_compat_layout(project_root, layout_profile=profile)
        if "gnina_hpc" not in compatibility_profiles:
            compatibility_profiles.append("gnina_hpc")

    merged_engines = sorted(set(existing.get("engines", [])) | set(engines))
    merged_source_paths = dict(existing.get("source_paths", {}))
    merged_source_paths.update(source_paths or {})
    merged_notes = list(existing.get("notes", []))
    for note in notes or []:
        if note not in merged_notes:
            merged_notes.append(note)
    if "gnina" in merged_engines:
        compat_note = "compatibility_profile=gnina_hpc"
        if compat_note not in merged_notes:
            merged_notes.append(compat_note)

    manifest = ProjectManifest(
        project_name=project_name or existing.get("project_name") or project_root.name,
        project_root=str(project_root),
        created_at=existing.get("created_at") or datetime.now(timezone.utc).isoformat(),
        asset_mode=existing.get("asset_mode") or asset_mode,
        engines=merged_engines,
        pairlist_file=str(pairlist_file),
        receptors_dir=str(layout["receptors"]),
        ligands_dir=str(layout["ligands"]),
        metadata_dir=str(layout["metadata"]),
        layout_profile=profile,
        docking_root=str(layout["docking_root"]),
        post_docking_root=str(layout["post_docking_root"]),
        raw_proteins_dir=str(layout["raw_proteins"]),
        raw_ligands_dir=str(layout["raw_ligands"]),
        raw_ligands_sdf_dir=str(layout["raw_ligands_sdf"]),
        prepared_proteins_dir=str(layout["prepared_proteins"]),
        prepared_ligands_dir=str(layout["prepared_ligands"]),
        enabled_panels=list(existing.get("enabled_panels", []) or merged_engines),
        pairlist_mode=pairlist_mode or str(existing.get("pairlist_mode", "")),
        has_cocrystal_benchmark_rows=bool(
            existing.get("has_cocrystal_benchmark_rows", False) or has_cocrystal_benchmark_rows
        ),
        source_paths=merged_source_paths,
        favorite_engine=favorite_engine or existing.get("favorite_engine", ""),
        engine_settings=existing.get("engine_settings", {engine: {} for engine in merged_engines}),
        pair_curation_state_file=str(existing.get("pair_curation_state_file", "") or pair_curation_state_path(project_root, profile)),
        latest_pair_round=existing.get("latest_pair_round", ""),
        deployment_root=str(existing.get("deployment_root", "") or deployment_root(project_root, profile)),
        latest_rerun_manifest_file=existing.get("latest_rerun_manifest_file", ""),
        top_pose_selection_policy=str(existing.get("top_pose_selection_policy", "best_affinity") or "best_affinity"),
        top_pose_global_aggregation=str(existing.get("top_pose_global_aggregation", "best_target") or "best_target"),
        hpc_profile=dict(existing.get("hpc_profile", {})),
        notes=merged_notes,
        protonation_policy=dict(existing.get(PROTONATION_POLICY_KEY) or {}),
        compatibility_profiles=compatibility_profiles,
    )
    manifest.engine_settings = {
        **{engine: {} for engine in merged_engines},
        **dict(existing.get("engine_settings", {})),
        **dict(manifest.engine_settings),
    }
    # Keys this function does not model (for example protonation_policy written by a newer command,
    # or a key from another writer) are kept. Modelled fields take the values computed above.
    payload = {**existing, **manifest.to_dict()}
    save_manifest(project_root, payload, layout_profile=profile)
    return {
        "project_root": str(project_root),
        "docking_root": str(layout["docking_root"]),
        "post_docking_root": str(layout["post_docking_root"]),
        "engines": merged_engines,
        "manifest_path": str(manifest_path(project_root, profile)),
        "pairlist_file": str(pairlist_file),
        "layout_type": profile,
        "compatibility_profiles": compatibility_profiles,
        "enabled_panels": manifest.enabled_panels,
    }
