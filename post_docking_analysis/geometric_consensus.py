"""Cross-engine agreement of mapped poses in the same receptor frame."""
from __future__ import annotations
import json
from itertools import combinations
from pathlib import Path
from typing import List, Optional
import numpy as np
import pandas as pd
from .pose_geometry import load_pose_molecule, fixed_frame_rmsd, selected_record

_GEOMETRIC_COLUMNS = ["tag", "geometric_agreement", "geometric_engine_count", "geometric_expected_pairs", "geometric_valid_pairs", "geometric_pass_pairs", "geometric_pairwise_rmsd_min", "geometric_pairwise_rmsd_mean", "geometric_pairwise_rmsd_max", "geometric_cutoff_angstrom", "geometric_reason", "geometric_pairwise_rmsd_json"]

def _normalize_element(value: str) -> str:
    token = str(value).strip()
    autodock = {"A":"C", "NA":"N", "OA":"O", "SA":"S", "HD":"H", "HS":"H"}
    return autodock.get(token, token.title()).upper()

def _extract_sdf_pose(path: Path, pose_index: int):
    try:
        mol = load_pose_molecule(path, pose_index)
        return mol.GetConformer().GetPositions(), [a.GetSymbol().upper() for a in mol.GetAtoms()], ""
    except Exception as exc:
        return np.empty((0, 3)), [], str(exc)

def _extract_pdb_like_pose(path: Path, pose_index: int):
    """Coordinate reader only; chemical mapping cannot be inferred from order."""
    try:
        coords, elements = [], []
        for line in selected_record(path, pose_index).splitlines():
            if not line.startswith(("ATOM", "HETATM")):
                continue
            token = line.split()[-1] if path.suffix.lower() == ".pdbqt" else line[76:78].strip()
            element = _normalize_element(token)
            if not element or element == "H":
                continue
            coords.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
            elements.append(element)
        if not coords:
            raise ValueError("no_heavy_atoms")
        return np.asarray(coords), elements, ""
    except Exception as exc:
        return np.empty((0,3)), [], str(exc)

def _extract_pose_signature(path: Path, pose_index: int):
    return _extract_sdf_pose(path, pose_index) if path.suffix.lower() in {".sdf", ".mol"} else _extract_pdb_like_pose(path, pose_index)

def _kabsch_rmsd(coords_a, coords_b):
    """Legacy name retained; docking coordinates are deliberately never fitted."""
    if coords_a.shape != coords_b.shape or not coords_a.size:
        return float("nan")
    return float(np.sqrt(np.mean(np.sum((coords_a-coords_b)**2,axis=1))))

def _pairwise_rmsd(engine_pose_data, rmsd_cutoff, pose_formats=None):
    """Coordinate arrays alone cannot establish chemical identity."""
    rows = [{"engine_a":a,"engine_b":b,"rmsd":None,"reason":"missing_authoritative_atom_mapping","soft_alignment_used":False,"soft_alignment_rmsd":None} for a,b in combinations(sorted(engine_pose_data),2)]
    return rows, len(rows), 0, 0, "pairwise_mapping_incomplete"

def compute_geometric_consensus(best_by_engine: pd.DataFrame, *, rmsd_cutoff: float = 2.0, expected_engines: Optional[List[str]] = None) -> pd.DataFrame:
    required = {"engine", "tag", "pose_file", "pose"}
    if best_by_engine is None or best_by_engine.empty or not required.issubset(best_by_engine.columns):
        return pd.DataFrame(columns=_GEOMETRIC_COLUMNS)
    if not np.isfinite(rmsd_cutoff) or rmsd_cutoff <= 0:
        raise ValueError("RMSD cutoff must be finite and positive")
    frame = best_by_engine.copy()
    frame["engine"] = frame["engine"].astype(str).str.strip().str.lower()
    engines = sorted(set(expected_engines or frame["engine"].unique()))
    cache, rows = {}, []
    for tag, group in frame.groupby("tag", dropna=False):
        poses, errors, frames = {}, {}, {}
        for engine in engines:
            records = group[group["engine"] == engine]
            if len(records) != 1:
                errors[engine] = "missing_engine_pose" if records.empty else "ambiguous_engine_pose"
                continue
            record = records.iloc[0]
            frame_id = record.get("receptor_frame_id")
            if pd.isna(frame_id) or not str(frame_id).strip():
                errors[engine] = "missing_receptor_frame_id"
                continue
            frames[engine] = str(frame_id)
            try:
                key = (str(record["pose_file"]), int(record["pose"]))
                if key not in cache:
                    cache[key] = load_pose_molecule(Path(key[0]),key[1])
                poses[engine] = cache[key]
            except Exception as exc:
                errors[engine] = str(exc)
        pair_rows = []
        for a,b in combinations(engines,2):
            value, reason = None, ""
            if a in errors or b in errors:
                reason = ";".join(f"{e}:{errors[e]}" for e in (a,b) if e in errors)
            elif frames[a] != frames[b]:
                reason = "receptor_frame_mismatch"
            else:
                try:
                    value = fixed_frame_rmsd(poses[a],poses[b])
                except Exception as exc:
                    reason = str(exc)
            pair_rows.append({"engine_a":a,"engine_b":b,"rmsd":value,"reason":reason,"soft_alignment_used":False,"soft_alignment_rmsd":None})
        finite = [p["rmsd"] for p in pair_rows if p["rmsd"] is not None and np.isfinite(p["rmsd"])]
        passed = sum(v <= rmsd_cutoff for v in finite)
        agree = bool(pair_rows) and len(finite) == len(pair_rows) and passed == len(pair_rows)
        reason = "all_pairs_within_cutoff" if agree else "insufficient_engine_poses" if not pair_rows else "pairwise_mapping_incomplete" if len(finite) < len(pair_rows) else "rmsd_above_cutoff"
        rows.append({"tag":str(tag),"geometric_agreement":agree,"geometric_engine_count":len(poses),"geometric_expected_pairs":len(pair_rows),"geometric_valid_pairs":len(finite),"geometric_pass_pairs":passed,"geometric_pairwise_rmsd_min":min(finite) if finite else np.nan,"geometric_pairwise_rmsd_mean":float(np.mean(finite)) if finite else np.nan,"geometric_pairwise_rmsd_max":max(finite) if finite else np.nan,"geometric_cutoff_angstrom":float(rmsd_cutoff),"geometric_reason":reason,"geometric_pairwise_rmsd_json":json.dumps(pair_rows,sort_keys=True)})
    return pd.DataFrame(rows,columns=_GEOMETRIC_COLUMNS)
