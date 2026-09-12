"""
Enhanced pose extractor for post-docking analysis pipeline.

This module handles the extraction of best poses from docking results,
organizing them into a centralized folder structure, and saving the
best pose for each complex into a separate folder.
"""
import pandas as pd
from pathlib import Path
import shutil
from typing import Dict, Optional, List
import csv
import re
import math

from post_docking_analysis.complex_validation import validate_complex_pdb_structure
from post_docking_analysis.score_semantics import sort_value
from post_docking_analysis.pose_geometry import selected_record


_BEST_POSE_CRITERIA_ALIASES = {
    "affinity": "vina_affinity",
    "vina": "vina_affinity",
    "vina_affinity": "vina_affinity",
    "cnn": "cnn_affinity",
    "cnn_affinity": "cnn_affinity",
}


def _infer_element_symbol(atom_name: str, fallback: str = "C") -> str:
    letters = "".join(ch for ch in str(atom_name or "").strip() if ch.isalpha())
    if not letters:
        return fallback
    if len(letters) >= 2 and letters[1].islower():
        return letters[:2]
    return letters[:2].capitalize()


def _sanitize_resname(value: str, fallback: str = "LIG") -> str:
    token = "".join(ch for ch in str(value or "").strip().upper() if ch.isalnum())
    if not token:
        token = fallback
    return token[:3].ljust(3, "X")


def _infer_ligand_resname_from_tag(tag: str, fallback: str = "LIG") -> str:
    token = str(tag or "").strip()
    # Common pattern: receptor_site_1_ligand[_...]
    match = re.search(r"_site_\d+_([^_]+)", token)
    if match:
        return _sanitize_resname(match.group(1), fallback=fallback)
    parts = [p for p in token.split("_") if p]
    return _sanitize_resname(parts[-1] if parts else fallback, fallback=fallback)


def _infer_receptor_file_from_tag(tag: str, receptors_dir: Path) -> Path:
    """
    Resolve receptor file from a canonical tag like receptor_site_1_ligand.
    """
    token = str(tag or "").strip()
    candidates: List[Path] = []

    site_match = re.match(r"^(?P<receptor>.+?)_(site_\d+)_(?P<ligand>.+)$", token)
    if site_match:
        receptor_name = str(site_match.group("receptor")).strip()
        if receptor_name:
            candidates.append(receptors_dir / f"{receptor_name}_prep.pdbqt")
            candidates.append(receptors_dir / f"{receptor_name}.pdbqt")
            candidates.append(receptors_dir / f"{receptor_name}.pdb")

    if "_prep_" in token:
        prefix = token.split("_prep_")[0]
        candidates.append(receptors_dir / f"{prefix}_prep.pdbqt")
        candidates.append(receptors_dir / f"{prefix}.pdbqt")
        candidates.append(receptors_dir / f"{prefix}.pdb")

    parts = [part for part in token.split("_") if part]
    if len(parts) >= 2:
        prefix2 = f"{parts[0]}_{parts[1]}"
        candidates.append(receptors_dir / f"{prefix2}_prep.pdbqt")
        candidates.append(receptors_dir / f"{prefix2}.pdbqt")
    elif parts:
        candidates.append(receptors_dir / f"{parts[0]}_prep.pdbqt")
        candidates.append(receptors_dir / f"{parts[0]}.pdbqt")

    for path in candidates:
        if path.exists():
            return path

    return candidates[0] if candidates else receptors_dir / f"{token}_prep.pdbqt"


def _renumber_atom_serials(lines: List[str]) -> List[str]:
    renumbered: List[str] = []
    serial = 1
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            line = line.rstrip("\n").ljust(80)
            line = f"{line[:6]}{serial:5d}{line[11:]}"
            serial += 1
        renumbered.append(line.rstrip("\n"))
    return renumbered


def _parse_pose_row(row: Dict[str, object]) -> Dict[str, object]:
    parsed = dict(row)
    parsed["tag"] = str(row.get("tag", ""))
    parsed["mode"] = int(float(row.get("mode", 0)))
    parsed["vina_affinity"] = float(row.get("vina_affinity", "nan"))
    try:
        parsed["cnn_affinity"] = float(row.get("cnn_affinity", "nan"))
    except Exception:
        parsed["cnn_affinity"] = float("nan")
    try:
        parsed["cnn_score"] = float(row.get("cnn_score", "nan"))
    except Exception:
        parsed["cnn_score"] = float("nan")
    return parsed


def _normalize_best_pose_criterion(value: object) -> str:
    token = str(value or "").strip().lower()
    return _BEST_POSE_CRITERIA_ALIASES.get(token, "vina_affinity")


def _numeric_or_inf(value: object) -> float:
    try:
        numeric = float(value)
    except Exception:
        return float("inf")
    return numeric if math.isfinite(numeric) else float("inf")


def _is_better_pose(
    candidate: Dict[str, object],
    incumbent: Dict[str, object],
    *,
    criterion: str = "vina_affinity",
) -> bool:
    """
    Deterministic tie-breaker:
    1) lower primary criterion wins
    2) lower vina_affinity wins
    3) lower pose mode wins
    """
    cand_primary = sort_value(candidate.get(criterion), criterion)
    inc_primary = sort_value(incumbent.get(criterion), criterion)
    if cand_primary < inc_primary:
        return True
    if cand_primary > inc_primary:
        return False

    cand_aff = _numeric_or_inf(candidate.get("vina_affinity"))
    inc_aff = _numeric_or_inf(incumbent.get("vina_affinity"))
    if cand_aff < inc_aff:
        return True
    if cand_aff > inc_aff:
        return False

    return int(candidate.get("mode", 0)) < int(incumbent.get("mode", 0))


def _extract_sdf_record_text(sdf_file: Path, pose_number: int) -> str:
    """
    Return one SDF mol record (1-based pose_number) from a possibly multi-pose SDF.
    """
    try:
        return selected_record(sdf_file, pose_number).rstrip() + "\n$$$$\n"
    except (ValueError, OSError):
        return ""


def _pdbqt_record_to_pdb_line(line: str) -> Optional[str]:
    """
    Convert one ATOM/HETATM PDBQT record to a simple PDB line.

    Prefer fixed-column parsing when the line is well-formed, but fall back to
    token parsing when spacing is irregular.
    """
    raw = str(line).rstrip("\n")
    if not raw.startswith(("ATOM", "HETATM")):
        return None

    record = raw[0:6].strip() or "ATOM"
    atom_num = raw[6:11].strip()
    atom_name = raw[12:16].strip()
    res_name = raw[17:20].strip() or "UNK"
    chain_id = raw[21:22].strip() or "A"
    res_num = raw[22:26].strip() or "1"
    occupancy = raw[54:60].strip() if len(raw) > 54 else "1.00"
    temp_factor = raw[60:66].strip() if len(raw) > 60 else "20.00"
    element = raw[76:78].strip() if len(raw) > 76 else ""

    try:
        x = float(raw[30:38].strip())
        y = float(raw[38:46].strip())
        z = float(raw[46:54].strip())
    except Exception:
        parts = raw.split()
        if len(parts) < 8:
            return None
        record = parts[0]
        atom_num = parts[1]
        atom_name = parts[2]
        res_name = parts[3] if len(parts) > 3 else "UNK"
        cursor = 4
        chain_id = "A"
        res_num = "1"

        # Identify where xyz coordinates begin.
        coord_start = None
        for idx in range(4, len(parts) - 2):
            try:
                float(parts[idx])
                float(parts[idx + 1])
                float(parts[idx + 2])
                coord_start = idx
                break
            except Exception:
                continue
        if coord_start is None:
            return None

        identity_tokens = parts[4:coord_start]
        if len(identity_tokens) >= 2:
            if len(identity_tokens[0]) == 1 and identity_tokens[0].isalpha():
                chain_id = identity_tokens[0]
                res_num = identity_tokens[1]
            else:
                res_num = identity_tokens[0]
        elif len(identity_tokens) == 1:
            token = identity_tokens[0]
            if len(token) == 1 and token.isalpha():
                chain_id = token
                res_num = "1"
            else:
                res_num = token
        cursor = coord_start
        if len(parts) < cursor + 3:
            return None
        try:
            x = float(parts[cursor])
            y = float(parts[cursor + 1])
            z = float(parts[cursor + 2])
        except Exception:
            return None
        occupancy = parts[cursor + 3] if len(parts) > cursor + 3 else "1.00"
        temp_factor = parts[cursor + 4] if len(parts) > cursor + 4 else "20.00"
        element = parts[-1] if parts[-1].isalpha() and len(parts[-1]) <= 2 else ""

    element = element or _infer_element_symbol(atom_name)
    return (
        f"{record:6s}{atom_num:>5s} {atom_name:<4s}{res_name:>3s} {chain_id:1s}{res_num:>4s}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{float(occupancy):6.2f}{float(temp_factor):6.2f}          {element:>2s}"
    )

def extract_best_poses_from_gnina(
    input_dir: Path,
    output_dir: Path,
    config: dict = None,
    gnina_dir: Optional[Path] = None,
    receptors_dir: Optional[Path] = None,
    scores_csv: Optional[Path] = None,
    best_pose_criterion: Optional[str] = None,
) -> int:
    """
    Extract best poses as PDB files using GNINA outputs in input_dir.
    
    Parameters
    ----------
    input_dir : Path
        Input directory containing GNINA outputs
    output_dir : Path
        Output directory for extracted poses
    config : dict, optional
        Configuration dictionary
        
    Returns
    -------
    int
        Number of PDBs written
    """
    # Use configuration or defaults
    if config is None:
        config = {}

    extract_all = config.get("pose_extraction", {}).get("extract_all_poses", False)
    criterion_token = (
        best_pose_criterion
        if best_pose_criterion is not None
        else config.get("pose_extraction", {}).get("best_pose_criteria", "affinity")
    )
    criterion = _normalize_best_pose_criterion(criterion_token)

    # Prefer caller-provided paths, then fall back to legacy discovery.
    possible_gnina_dirs = []
    if gnina_dir is not None:
        possible_gnina_dirs.append(Path(gnina_dir))
    possible_gnina_dirs.extend([
        input_dir / "gnina_out",
        input_dir / "gnina_out_cox2",
        input_dir / "gnina_out_inha"
    ])
    
    gnina_dir = None
    for possible_dir in possible_gnina_dirs:
        if possible_dir.exists():
            gnina_dir = possible_dir
            break
    
    if gnina_dir is None:
        print(f"❌ GNINA output directory not found in {input_dir}")
        return 0
    
    receptors_dir = Path(receptors_dir) if receptors_dir is not None else input_dir / "receptors"
    scores_csv = Path(scores_csv) if scores_csv is not None else gnina_dir / "all_scores.csv"

    if not scores_csv.exists():
        print(f"❌ Scores CSV not found: {scores_csv}")
        return 0

    # Create output directories
    best_poses_dir = output_dir / "best_poses"
    all_poses_dir = output_dir / "all_poses" if extract_all else None
    
    best_poses_dir.mkdir(parents=True, exist_ok=True)
    if all_poses_dir:
        all_poses_dir.mkdir(parents=True, exist_ok=True)

    # Read CSV and pick best mode per tag according to the configured criterion.
    rows: List[Dict[str, object]] = []
    with scores_csv.open() as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                parsed = _parse_pose_row(r)
            except Exception:
                continue
            rows.append(parsed)

    if not rows:
        print(f"⚠️  No rows in {scores_csv}")
        return 0

    # Group by tag
    if extract_all:
        # For extracting all poses, we'll process all rows
        poses_to_extract = rows
    else:
        # For best poses only, we'll pick the best per tag
        best_by_tag: Dict[str, Dict[str, object]] = {}
        rows.sort(key=lambda rec: (str(rec.get("tag", "")), int(rec.get("mode", 0))))
        for r in rows:
            tag = str(r.get("tag", ""))
            if tag not in best_by_tag or _is_better_pose(r, best_by_tag[tag], criterion=criterion):
                best_by_tag[tag] = r
        poses_to_extract = list(best_by_tag.values())

    written = 0
    manifest_rows: List[Dict[str, object]] = []
    for r in poses_to_extract:
        tag = r['tag']
        pose_number = int(r.get("mode", 0) or 0)
        ligand_resname = _infer_ligand_resname_from_tag(tag)
        sdf_file = gnina_dir / f"{tag}_top.sdf"
        if not sdf_file.exists():
            fallback_sdf = gnina_dir / f"{tag}.sdf"
            if fallback_sdf.exists():
                sdf_file = fallback_sdf

        # Determine output directory based on extraction type
        if extract_all:
            out_dir = all_poses_dir
        else:
            # Create a separate folder for each complex
            complex_dir = best_poses_dir / tag
            complex_dir.mkdir(exist_ok=True)
            out_dir = complex_dir

        out_pdb = out_dir / f"{tag}_pose{pose_number}.pdb"

        receptor_file = _infer_receptor_file_from_tag(tag, receptors_dir)

        manifest_row: Dict[str, object] = {
            "tag": tag,
            "selected_pose": pose_number,
            "selection_criterion": criterion,
            "selected_score": r.get(criterion),
            "vina_affinity": r.get("vina_affinity"),
            "cnn_affinity": r.get("cnn_affinity"),
            "cnn_score": r.get("cnn_score"),
            "input_scores_csv": str(scores_csv),
            "input_sdf_file": str(sdf_file),
            "input_receptor_file": str(receptor_file),
            "output_pdb": str(out_pdb),
            "status": "",
            "validation_is_valid": False,
            "validation_errors": "",
            "validation_warnings": "",
            "receptor_atom_count": 0,
            "ligand_atom_count": 0,
        }

        if pose_number <= 0:
            manifest_row["status"] = "invalid_pose_number"
            manifest_rows.append(manifest_row)
            print(f"⚠️  Invalid pose index for tag {tag}: mode={r.get('mode')}")
            continue

        if not sdf_file.exists():
            manifest_row["status"] = "missing_sdf"
            manifest_rows.append(manifest_row)
            print(f"⚠️  SDF not found for tag {tag}: {sdf_file}")
            continue

        if not receptor_file.exists():
            manifest_row["status"] = "missing_receptor"
            manifest_rows.append(manifest_row)
            print(f"⚠️  Receptor file not found: {receptor_file}")
            continue

        pose_record_text = _extract_sdf_record_text(sdf_file, pose_number)
        if not pose_record_text:
            manifest_row["status"] = "missing_sdf_pose_record"
            manifest_rows.append(manifest_row)
            print(f"⚠️  Pose {pose_number} not found in SDF for tag {tag}: {sdf_file}")
            continue
        out_pdb.with_suffix(".ligand.sdf").write_text(pose_record_text, encoding="utf-8")

        # Try to get docking center coordinates from log file
        log_file = gnina_dir / f"{tag}.log"
        manifest_row["input_log_file"] = str(log_file)
        if log_file.exists():
            try:
                with open(log_file, 'r') as f:
                    log_content = f.read()
                    # Extract center coordinates from command line
                    import re
                    center_match = re.search(r'--center_x\s+([\d.-]+)\s+--center_y\s+([\d.-]+)\s+--center_z\s+([\d.-]+)', log_content)
                    if center_match:
                        docking_center = (
                            float(center_match.group(1)),
                            float(center_match.group(2)),
                            float(center_match.group(3)),
                        )
                        manifest_row["docking_center_xyz"] = ",".join(f"{value:.3f}" for value in docking_center)
                        print(f"📍 Found docking center: {docking_center}")
            except Exception as e:
                print(f"⚠️  Could not extract docking center: {e}")

        wrote_output = False
        fallback_reason = ""

        # Combine receptor and ligand to create complex using OpenBabel
        try:
            from openbabel import pybel

            # Read receptor PDBQT file
            receptor_lines = []
            receptor_mol = next(pybel.readfile("pdbqt", str(receptor_file)))
            receptor_pdb = receptor_mol.write("pdb")
            for line in receptor_pdb.split('\n'):
                if line.startswith(('ATOM', 'HETATM')):
                    # Preserve receptor chain IDs and cofactors.
                    line = line.ljust(80)
                    new_line = line
                    receptor_lines.append(new_line)

            # Read ligand pose record from selected SDF conformer block
            out_pdb.with_suffix(".receptor.pdb").write_text("\n".join(receptor_lines + ["END"]) + "\n", encoding="utf-8")
            ligand_lines = []
            ligand_mol = pybel.readstring("sdf", pose_record_text)
            ligand_pdb = ligand_mol.write("pdb")
            for line in ligand_pdb.split('\n'):
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # Fix the line format and assign chain B
                    line = line.ljust(80)
                    new_line = f"HETATM{line[6:21]}B{line[22:]}"
                    new_line = new_line[:17] + ligand_resname + new_line[20:]
                    ligand_lines.append(new_line)

            # Preserve the selected ligand graph for downstream interactions.
            out_pdb.with_suffix(".ligand.sdf").write_text(pose_record_text, encoding="utf-8")
            # Combine receptor and ligand
            all_lines = _renumber_atom_serials(receptor_lines + ligand_lines) + ["END"]
            combined_content = '\n'.join(all_lines)

            # Write combined complex
            with open(out_pdb, 'w') as f:
                f.write(combined_content)

            wrote_output = True
            print(f"✅ Extracted complex {out_pdb.name} (receptor + ligand)")

        except ImportError:
            fallback_reason = "openbabel_unavailable"
            print("⚠️  OpenBabel not available, using fallback method")
            # Fallback: use simple SDF to PDB conversion
            try:
                ligand_pdb_content = _convert_sdf_to_pdb_simple(
                    sdf_file,
                    resname=ligand_resname,
                    chain_id="B",
                    pose_number=pose_number,
                )
                if ligand_pdb_content:
                    # Read receptor PDBQT file manually
                    with open(receptor_file, 'r') as f:
                        receptor_content = f.read()
                    
                    # Convert PDBQT to PDB format (remove Q and T columns)
                    receptor_pdb_lines = []
                    for line in receptor_content.split('\n'):
                        if line.startswith(('ATOM', 'HETATM')):
                            pdb_line = _pdbqt_record_to_pdb_line(line)
                            if pdb_line:
                                receptor_pdb_lines.append(pdb_line)
                        elif line.startswith(('REMARK', 'HEADER', 'TITLE', 'COMPND', 'SOURCE', 'AUTHOR', 'REVDAT', 'JRNL', 'SEQRES', 'HET', 'FORMUL', 'HELIX', 'SHEET', 'SSBOND', 'LINK', 'CISPEP', 'SITE', 'CRYST1', 'ORIGX1', 'ORIGX2', 'ORIGX3', 'SCALE1', 'SCALE2', 'SCALE3', 'MTRIX1', 'MTRIX2', 'MTRIX3', 'TVECT', 'MODEL', 'ENDMDL')):
                            receptor_pdb_lines.append(line)
                    
                    # Combine receptor and ligand
                    out_pdb.with_suffix(".receptor.pdb").write_text("\n".join(receptor_pdb_lines + ["END"]) + "\n", encoding="utf-8")
                    combined_content = []
                    combined_content.extend(receptor_pdb_lines)
                    combined_content.append("")  # Empty line separator
                    combined_content.extend(ligand_pdb_content.split('\n'))
                    combined_content = _renumber_atom_serials(combined_content)
                    
                    # Write combined complex
                    with open(out_pdb, 'w') as f:
                        f.write('\n'.join(combined_content))

                    wrote_output = True
                    print(f"✅ Extracted complex {out_pdb.name} (receptor + ligand, fallback method)")
                else:
                    print(f"⚠️  Failed to convert ligand from {sdf_file}")
            except Exception as e:
                print(f"⚠️  Error creating complex for {tag}: {e}")
                # Final fallback: just copy the SDF file
                try:
                    shutil.copy2(sdf_file, out_pdb)
                    wrote_output = True
                    fallback_reason = fallback_reason or "fallback_copy_sdf"
                    print(f"✅ Copied SDF as PDB: {out_pdb.name}")
                except Exception as e2:
                    print(f"❌ Failed to copy {sdf_file}: {e2}")
        except Exception as e:
            fallback_reason = fallback_reason or "openbabel_or_merge_error"
            print(f"⚠️  Error creating complex for {tag}: {e}")
            # Fallback: just copy the SDF file
            try:
                shutil.copy2(sdf_file, out_pdb)
                wrote_output = True
                fallback_reason = "fallback_copy_sdf"
                print(f"✅ Copied SDF as PDB: {out_pdb.name}")
            except Exception as e2:
                print(f"❌ Failed to copy {sdf_file}: {e2}")

        if wrote_output:
            validation = validate_complex_pdb_structure(out_pdb, receptor_reference=receptor_file)
            validation_errors = [str(item) for item in (validation.get("errors") or [])]
            validation_warnings = [str(item) for item in (validation.get("warnings") or [])]
            manifest_row["validation_is_valid"] = bool(validation.get("is_valid", False))
            manifest_row["validation_errors"] = "; ".join(validation_errors)
            manifest_row["validation_warnings"] = "; ".join(validation_warnings)
            manifest_row["receptor_atom_count"] = int(validation.get("receptor_atom_count", 0) or 0)
            manifest_row["ligand_atom_count"] = int(validation.get("ligand_atom_count", 0) or 0)
            if validation.get("is_valid", False):
                manifest_row["status"] = "extracted"
            else:
                manifest_row["status"] = "invalid_complex_output"
            if fallback_reason:
                manifest_row["fallback_reason"] = fallback_reason
            written += 1
        else:
            if not manifest_row["status"]:
                manifest_row["status"] = "extraction_failed"
            if fallback_reason:
                manifest_row["fallback_reason"] = fallback_reason

        manifest_rows.append(manifest_row)

    if manifest_rows:
        manifest_file = output_dir / "pose_extraction_manifest.csv"
        output_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(manifest_rows).to_csv(manifest_file, index=False)
        print(f"🧾 Pose extraction manifest saved: {manifest_file}")

    print(f"✅ Extracted {written} poses to: {best_poses_dir}")
    if all_poses_dir:
        print(f"   (All poses saved to: {all_poses_dir})")
    return written

def organize_poses_by_affinity(best_poses_dir: Path, threshold: float = -8.0):
    """
    Organize extracted poses into folders based on binding affinity.
    
    Parameters
    ----------
    best_poses_dir : Path
        Directory containing extracted best poses
    threshold : float
        Affinity threshold for strong binders (kcal/mol)
    """
    strong_binders_dir = best_poses_dir / "strong_binders"
    moderate_binders_dir = best_poses_dir / "moderate_binders"
    weak_binders_dir = best_poses_dir / "weak_binders"
    
    strong_binders_dir.mkdir(exist_ok=True)
    moderate_binders_dir.mkdir(exist_ok=True)
    weak_binders_dir.mkdir(exist_ok=True)
    
    # Read the scores CSV to get affinity values
    # Try multiple possible locations for the CSV
    possible_csv_paths = [
        best_poses_dir.parent.parent / "gnina_out" / "all_scores.csv",
        best_poses_dir.parent / "reports" / "full_data.csv",
        best_poses_dir.parent.parent / "test_docking_data" / "gnina_out" / "all_scores.csv"
    ]
    
    scores_csv = None
    for path in possible_csv_paths:
        if path.exists():
            scores_csv = path
            break
    
    if not scores_csv:
        print("⚠️  Scores CSV not found for organizing poses")
        return
    
    df = pd.read_csv(scores_csv)
    
    # Filter out failed docking attempts (positive values)
    original_count = len(df)
    df = df[df['vina_affinity'] < 0]
    failed_count = original_count - len(df)
    if failed_count > 0:
        print(f"🚫 Filtered out {failed_count} failed docking attempts for pose organization")
    
    # Move pose files based on affinity
    for pdb_file in best_poses_dir.rglob("*.pdb"):
        if pdb_file.is_file() and pdb_file.parent != strong_binders_dir and \
           pdb_file.parent != moderate_binders_dir and pdb_file.parent != weak_binders_dir:
            
            # Extract tag from filename
            filename = pdb_file.stem
            # Remove _poseXX part (e.g., _pose1, _pose2, etc.)
            if filename.endswith('_pose1'):
                tag = filename[:-6]  # Remove '_pose1'
            else:
                # Fallback: remove last two parts if they look like pose numbers
                parts = filename.split("_")
                if len(parts) >= 2 and parts[-1].startswith('pose'):
                    tag = "_".join(parts[:-1])
                else:
                    tag = filename
            
            # Find affinity for this tag - try both 'tag' and 'complex_name' columns
            affinity_row = None
            if 'tag' in df.columns:
                affinity_row = df[df['tag'] == tag]
            elif 'complex_name' in df.columns:
                affinity_row = df[df['complex_name'] == tag]
            
            if affinity_row is not None and not affinity_row.empty:
                affinity = affinity_row.iloc[0]['vina_affinity']
                
                # Move to appropriate directory
                if affinity <= threshold:
                    target_dir = strong_binders_dir
                elif affinity <= -6.0:
                    target_dir = moderate_binders_dir
                else:
                    target_dir = weak_binders_dir
                    
                target_file = target_dir / pdb_file.name
                shutil.move(str(pdb_file), str(target_file))
                print(f"📁 Organized {pdb_file.name} to {target_dir.name}")

def create_pose_summary_report(best_poses_dir: Path, output_dir: Path):
    """
    Create a summary report of extracted poses.
    
    Parameters
    ----------
    best_poses_dir : Path
        Directory containing extracted best poses
    output_dir : Path
        Output directory for reports
    """
    # Read the scores CSV - try multiple possible locations
    possible_csv_paths = [
        best_poses_dir.parent.parent / "gnina_out" / "all_scores.csv",
        best_poses_dir.parent / "reports" / "full_data.csv",
        best_poses_dir.parent.parent / "test_docking_data" / "gnina_out" / "all_scores.csv"
    ]
    
    scores_csv = None
    for path in possible_csv_paths:
        if path.exists():
            scores_csv = path
            break
    
    if not scores_csv:
        print("⚠️  Scores CSV not found for creating summary report")
        return
    
    df = pd.read_csv(scores_csv)
    
    # Create summary report
    summary_data = []
    for pdb_file in best_poses_dir.rglob("*.pdb"):
        if pdb_file.is_file():
            # Extract tag from filename
            filename = pdb_file.stem
            # Remove _poseXX part (e.g., _pose1, _pose2, etc.)
            if filename.endswith('_pose1'):
                tag = filename[:-6]  # Remove '_pose1'
            else:
                # Fallback: remove last two parts if they look like pose numbers
                parts = filename.split("_")
                if len(parts) >= 2 and parts[-1].startswith('pose'):
                    tag = "_".join(parts[:-1])
                else:
                    tag = filename
            
            # Find data for this tag - try both 'tag' and 'complex_name' columns
            tag_data = None
            if 'tag' in df.columns:
                tag_data = df[df['tag'] == tag]
            elif 'complex_name' in df.columns:
                tag_data = df[df['complex_name'] == tag]
            
            if tag_data is not None and not tag_data.empty:
                row = tag_data.iloc[0]
                summary_entry = {
                    'complex': tag,
                    'vina_affinity': row['vina_affinity'],
                    'pdb_file': str(pdb_file.relative_to(best_poses_dir))
                }
                
                # Add optional columns if they exist
                if 'cnn_affinity' in row:
                    summary_entry['cnn_affinity'] = row['cnn_affinity']
                if 'cnn_score' in row:
                    summary_entry['cnn_score'] = row['cnn_score']
                if 'mode' in row:
                    summary_entry['mode'] = row['mode']
                
                summary_data.append(summary_entry)
    
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df = summary_df.sort_values('vina_affinity')
        
        # Save to CSV
        summary_file = output_dir / "pose_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        print(f"✅ Pose summary report saved to: {summary_file}")
        
        # Save to Excel if pandas Excel support is available
        try:
            excel_file = output_dir / "pose_summary.xlsx"
            summary_df.to_excel(excel_file, index=False)
            print(f"✅ Pose summary Excel saved to: {excel_file}")
        except ImportError:
            print("⚠️  Excel support not available, skipping Excel report")
    else:
        print("⚠️  No pose data found for summary report")

def _convert_sdf_to_pdb_simple(
    sdf_file: Path,
    resname: str = "LIG",
    chain_id: str = "A",
    pose_number: int = 1,
) -> str:
    """
    Simple SDF to PDB converter that extracts coordinates and creates basic PDB format.
    
    Parameters
    ----------
    sdf_file : Path
        Path to SDF file
        
    Returns
    -------
    str
        PDB content as string, or empty string if conversion fails
    """
    try:
        pose_record = _extract_sdf_record_text(sdf_file, pose_number)
        if not pose_record:
            return ""
        lines = pose_record.splitlines(keepends=True)
        
        # Find the counts line (line 4 in SDF format)
        if len(lines) < 4:
            return ""
        
        counts_line = lines[3].strip()
        if not counts_line or len(counts_line.split()) < 3:
            return ""
        
        # Parse atom count
        try:
            atom_count = int(counts_line.split()[0])
        except (ValueError, IndexError):
            return ""
        
        # Extract atom coordinates (lines 5 to 5+atom_count-1)
        pdb_lines = []
        atom_num = 1
        safe_resname = _sanitize_resname(resname, fallback="LIG")
        safe_chain = str(chain_id or "A").strip()[:1] or "A"
        
        for i in range(4, 4 + atom_count):
            if i >= len(lines):
                break
                
            line = lines[i].strip()
            if not line:
                continue
                
            parts = line.split()
            if len(parts) >= 4:
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    z = float(parts[2])
                    element = parts[3] if len(parts) > 3 else "C"
                    
                    # Create PDB ATOM line with proper formatting
                    pdb_line = (
                        f"HETATM{atom_num:5d}  {element:2s}  {safe_resname:3s} {safe_chain}{atom_num:4d}    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           {element:2s}"
                    )
                    pdb_lines.append(pdb_line)
                    atom_num += 1
                    
                except (ValueError, IndexError):
                    continue
        
        if pdb_lines:
            return '\n'.join(pdb_lines)
        else:
            return ""
            
    except Exception as e:
        print(f"⚠️  Simple SDF conversion failed: {e}")
        return ""
