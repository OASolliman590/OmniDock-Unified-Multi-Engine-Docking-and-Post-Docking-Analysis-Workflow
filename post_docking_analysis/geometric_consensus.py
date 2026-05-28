"""
DockBox-style geometric consensus helpers.

This module computes cross-engine pose agreement by comparing best poses for
the same tag using pairwise heavy-atom RMSD after Kabsch alignment.
"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


_GEOMETRIC_COLUMNS = [
    "tag",
    "geometric_agreement",
    "geometric_engine_count",
    "geometric_expected_pairs",
    "geometric_valid_pairs",
    "geometric_pass_pairs",
    "geometric_pairwise_rmsd_min",
    "geometric_pairwise_rmsd_mean",
    "geometric_pairwise_rmsd_max",
    "geometric_cutoff_angstrom",
    "geometric_reason",
    "geometric_pairwise_rmsd_json",
]


def _normalize_element(value: str) -> str:
    letters = "".join(ch for ch in str(value or "").strip() if ch.isalpha())
    if not letters:
        return ""
    if len(letters) >= 2 and letters[1].islower():
        return letters[:2].upper()
    return letters[:1].upper()


def _parse_pdb_atom_line(line: str) -> Optional[Tuple[float, float, float, str]]:
    if not str(line).startswith(("ATOM", "HETATM")):
        return None
    try:
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        element = _normalize_element(line[76:78].strip() or line[12:16].strip())
        if not element:
            return None
        return (x, y, z, element)
    except Exception:
        parts = str(line).split()
        if len(parts) < 8:
            return None
        coord_idx = None
        for idx in range(0, len(parts) - 2):
            try:
                float(parts[idx])
                float(parts[idx + 1])
                float(parts[idx + 2])
                coord_idx = idx
                break
            except Exception:
                continue
        if coord_idx is None:
            return None
        try:
            x = float(parts[coord_idx])
            y = float(parts[coord_idx + 1])
            z = float(parts[coord_idx + 2])
        except Exception:
            return None
        element = _normalize_element(parts[-1])
        if not element:
            atom_name = parts[2] if len(parts) > 2 else ""
            element = _normalize_element(atom_name)
        if not element:
            return None
        return (x, y, z, element)


def _extract_pdb_like_pose(path: Path, pose_index: int) -> Tuple[np.ndarray, List[str], str]:
    try:
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception as exc:
        return np.empty((0, 3), dtype=float), [], f"read_error:{exc}"

    models: List[List[str]] = []
    current: List[str] = []
    saw_model = False
    base_atoms: List[str] = []

    for line in lines:
        if line.startswith("MODEL"):
            saw_model = True
            current = []
            continue
        if line.startswith("ENDMDL"):
            models.append(current)
            current = []
            continue
        if line.startswith(("ATOM", "HETATM")):
            if saw_model:
                current.append(line)
            else:
                base_atoms.append(line)

    if saw_model and current:
        models.append(current)
    if not saw_model:
        models = [base_atoms]

    if not models:
        return np.empty((0, 3), dtype=float), [], "no_models"

    idx = max(int(pose_index or 1) - 1, 0)
    if idx >= len(models):
        idx = 0
    atom_lines = models[idx]
    if not atom_lines:
        return np.empty((0, 3), dtype=float), [], "empty_model"

    coords: List[List[float]] = []
    elements: List[str] = []
    for line in atom_lines:
        parsed = _parse_pdb_atom_line(line)
        if parsed is None:
            continue
        x, y, z, element = parsed
        if element.startswith("H"):
            continue
        coords.append([x, y, z])
        elements.append(element)

    if not coords:
        return np.empty((0, 3), dtype=float), [], "no_heavy_atoms"
    return np.array(coords, dtype=float), elements, ""


def _extract_sdf_pose(path: Path, pose_index: int) -> Tuple[np.ndarray, List[str], str]:
    # Semi-public parser utility: reused by redocking_validation.py for
    # consistent multi-pose SDF block handling.
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except Exception as exc:
        return np.empty((0, 3), dtype=float), [], f"read_error:{exc}"

    blocks = [block for block in text.split("$$$$") if block.strip()]
    if not blocks:
        return np.empty((0, 3), dtype=float), [], "no_molecule_blocks"

    idx = max(int(pose_index or 1) - 1, 0)
    if idx >= len(blocks):
        idx = 0
    lines = blocks[idx].splitlines()
    if len(lines) < 4:
        return np.empty((0, 3), dtype=float), [], "invalid_sdf_header"

    counts = lines[3]
    atom_count = 0
    try:
        atom_count = int(counts[0:3].strip())
    except Exception:
        try:
            atom_count = int(counts.split()[0])
        except Exception:
            return np.empty((0, 3), dtype=float), [], "invalid_sdf_counts"

    start = 4
    end = min(start + atom_count, len(lines))
    coords: List[List[float]] = []
    elements: List[str] = []
    for line in lines[start:end]:
        try:
            x = float(line[0:10].strip())
            y = float(line[10:20].strip())
            z = float(line[20:30].strip())
            element = _normalize_element(line[31:34].strip())
        except Exception:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                x = float(parts[0])
                y = float(parts[1])
                z = float(parts[2])
            except Exception:
                continue
            element = _normalize_element(parts[3])
        if not element or element.startswith("H"):
            continue
        coords.append([x, y, z])
        elements.append(element)

    if not coords:
        return np.empty((0, 3), dtype=float), [], "no_heavy_atoms"
    return np.array(coords, dtype=float), elements, ""


def _extract_pose_signature(path: Path, pose_index: int) -> Tuple[np.ndarray, List[str], str]:
    suffix = path.suffix.lower()
    if suffix == ".sdf":
        return _extract_sdf_pose(path, pose_index)
    if suffix in {".pdb", ".pdbqt"}:
        return _extract_pdb_like_pose(path, pose_index)
    # Conservative fallback to pdb-like parser for unknown text formats.
    return _extract_pdb_like_pose(path, pose_index)


def _kabsch_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    if coords_a.shape != coords_b.shape or coords_a.shape[0] <= 0:
        return float("nan")
    a = coords_a - np.mean(coords_a, axis=0)
    b = coords_b - np.mean(coords_b, axis=0)
    covariance = np.dot(a.T, b)
    u, _, vt = np.linalg.svd(covariance)
    rot = np.dot(vt.T, u.T)
    if np.linalg.det(rot) < 0:
        vt[-1, :] *= -1
        rot = np.dot(vt.T, u.T)
    aligned = np.dot(a, rot)
    diff = aligned - b
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def _centroid_sorted(coords: np.ndarray) -> np.ndarray:
    if coords.shape[0] <= 1:
        return coords
    centroid = np.mean(coords, axis=0)
    distances = np.linalg.norm(coords - centroid, axis=1)
    order = np.argsort(distances)
    return coords[order]


def _pairwise_rmsd(
    engine_pose_data: Dict[str, Tuple[np.ndarray, List[str]]],
    rmsd_cutoff: float,
    pose_formats: Optional[Dict[str, str]] = None,
) -> Tuple[List[Dict[str, object]], int, int, int, str]:
    pair_rows: List[Dict[str, object]] = []
    engines = sorted(engine_pose_data.keys())
    expected_pairs = len(engines) * (len(engines) - 1) // 2
    valid_pairs = 0
    pass_pairs = 0
    normalized_formats = {
        str(engine).strip().lower(): str(fmt or "").strip().lower()
        for engine, fmt in dict(pose_formats or {}).items()
        if str(engine).strip()
    }

    for engine_a, engine_b in combinations(engines, 2):
        coords_a, elements_a = engine_pose_data[engine_a]
        coords_b, elements_b = engine_pose_data[engine_b]
        rmsd = np.nan
        reason = ""
        soft_alignment_used = False
        soft_alignment_rmsd: Optional[float] = None
        count_a = int(coords_a.shape[0])
        count_b = int(coords_b.shape[0])
        format_a = normalized_formats.get(str(engine_a).strip().lower(), "")
        format_b = normalized_formats.get(str(engine_b).strip().lower(), "")
        formats_differ = bool(format_a and format_b and format_a != format_b)
        mismatch = ""
        if count_a != count_b:
            mismatch = "atom_count_mismatch"
        elif elements_a != elements_b:
            mismatch = "element_sequence_mismatch"

        if not mismatch:
            rmsd = _kabsch_rmsd(coords_a, coords_b)
            if np.isfinite(rmsd):
                valid_pairs += 1
                if rmsd <= float(rmsd_cutoff):
                    pass_pairs += 1
            else:
                reason = "rmsd_not_finite"
        else:
            reason = mismatch
            count_delta = abs(count_a - count_b)
            if mismatch == "atom_count_mismatch" and formats_differ and count_delta <= 2:
                n_common = min(count_a, count_b)
                if n_common > 0:
                    soft_a = _centroid_sorted(coords_a)[:n_common]
                    soft_b = _centroid_sorted(coords_b)[:n_common]
                    soft_rmsd = _kabsch_rmsd(soft_a, soft_b)
                    soft_alignment_used = True
                    if np.isfinite(soft_rmsd):
                        soft_alignment_rmsd = float(soft_rmsd)
                        rmsd = float(soft_rmsd)
                        reason = "soft_aligned_count_mismatch"
                        valid_pairs += 1
                        if float(soft_rmsd) <= float(rmsd_cutoff):
                            pass_pairs += 1
                    else:
                        reason = "soft_alignment_rmsd_not_finite"
        pair_rows.append(
            {
                "engine_a": engine_a,
                "engine_b": engine_b,
                "rmsd": float(rmsd) if np.isfinite(rmsd) else None,
                "reason": reason,
                "soft_alignment_used": bool(soft_alignment_used),
                "soft_alignment_rmsd": soft_alignment_rmsd,
            }
        )

    if len(engines) < 2:
        status_reason = "insufficient_engine_poses"
    elif expected_pairs == 0:
        status_reason = "insufficient_pairs"
    elif valid_pairs < expected_pairs:
        status_reason = "pairwise_mapping_incomplete"
    elif pass_pairs == expected_pairs:
        status_reason = "all_pairs_within_cutoff"
    else:
        status_reason = "rmsd_above_cutoff"

    return pair_rows, expected_pairs, valid_pairs, pass_pairs, status_reason


def compute_geometric_consensus(
    best_by_engine: pd.DataFrame,
    *,
    rmsd_cutoff: float = 2.0,
) -> pd.DataFrame:
    """
    Compute DockBox-style geometric consensus per tag.

    Parameters
    ----------
    best_by_engine : pd.DataFrame
        One best pose per (engine, tag). Must contain columns:
        engine, tag, pose_file, pose.
    rmsd_cutoff : float
        RMSD cutoff in Angstrom for pairwise geometric agreement.
    """
    required = {"engine", "tag", "pose_file", "pose"}
    if best_by_engine is None or best_by_engine.empty or not required.issubset(best_by_engine.columns):
        return pd.DataFrame(columns=_GEOMETRIC_COLUMNS)

    frame = best_by_engine.copy()
    frame["engine"] = frame["engine"].astype(str).str.strip().str.lower()
    frame["tag"] = frame["tag"].astype(str)
    frame["pose"] = pd.to_numeric(frame["pose"], errors="coerce").fillna(1).astype(int)
    frame["pose_file"] = frame["pose_file"].astype(str)

    signature_cache: Dict[Tuple[str, int], Tuple[np.ndarray, List[str], str]] = {}
    rows: List[Dict[str, object]] = []

    for tag, group in frame.groupby("tag", dropna=False):
        engine_pose_data: Dict[str, Tuple[np.ndarray, List[str]]] = {}
        pose_formats: Dict[str, str] = {}
        parse_errors: Dict[str, str] = {}

        for _, record in group.iterrows():
            engine = str(record.get("engine") or "").strip().lower()
            pose_file = Path(str(record.get("pose_file") or "")).expanduser()
            pose_index = int(record.get("pose") or 1)
            if engine and pose_file.suffix:
                pose_formats[engine] = str(pose_file.suffix).strip().lower()
            key = (str(pose_file), pose_index)
            if key not in signature_cache:
                if not pose_file.exists() or not pose_file.is_file():
                    signature_cache[key] = (np.empty((0, 3), dtype=float), [], "missing_pose_file")
                else:
                    signature_cache[key] = _extract_pose_signature(pose_file, pose_index)
            coords, elements, error = signature_cache[key]
            if coords.size == 0 or not elements:
                parse_errors[engine] = error or "pose_parse_failed"
                continue
            engine_pose_data[engine] = (coords, elements)

        pair_rows, expected_pairs, valid_pairs, pass_pairs, status_reason = _pairwise_rmsd(
            engine_pose_data,
            rmsd_cutoff=float(rmsd_cutoff),
            pose_formats=pose_formats,
        )
        finite_rmsd = [float(row["rmsd"]) for row in pair_rows if row.get("rmsd") is not None]
        hard_failure_reasons = {
            "missing_pose_file",
            "atom_count_mismatch",
            "element_sequence_mismatch",
            "rmsd_not_finite",
            "soft_alignment_rmsd_not_finite",
        }
        all_pairs_within_cutoff = bool(pair_rows) and all(
            (
                (
                    row.get("rmsd") is not None
                    and float(row.get("rmsd")) <= float(rmsd_cutoff)
                )
                or (
                    bool(row.get("soft_alignment_used"))
                    and row.get("soft_alignment_rmsd") is not None
                    and float(row.get("soft_alignment_rmsd")) <= float(rmsd_cutoff)
                )
            )
            for row in pair_rows
        )
        no_hard_failures = all(
            not (
                str(row.get("reason") or "") in hard_failure_reasons
                and not bool(row.get("soft_alignment_used"))
            )
            for row in pair_rows
        )
        geometric_agreement = bool(expected_pairs > 0 and all_pairs_within_cutoff and no_hard_failures)

        reason = status_reason
        if parse_errors and reason in {"insufficient_engine_poses", "pairwise_mapping_incomplete"}:
            reason = f"{reason};parse_errors={','.join(sorted(parse_errors.keys()))}"

        rows.append(
            {
                "tag": str(tag),
                "geometric_agreement": geometric_agreement,
                "geometric_engine_count": int(len(engine_pose_data)),
                "geometric_expected_pairs": int(expected_pairs),
                "geometric_valid_pairs": int(valid_pairs),
                "geometric_pass_pairs": int(pass_pairs),
                "geometric_pairwise_rmsd_min": float(min(finite_rmsd)) if finite_rmsd else np.nan,
                "geometric_pairwise_rmsd_mean": float(sum(finite_rmsd) / len(finite_rmsd)) if finite_rmsd else np.nan,
                "geometric_pairwise_rmsd_max": float(max(finite_rmsd)) if finite_rmsd else np.nan,
                "geometric_cutoff_angstrom": float(rmsd_cutoff),
                "geometric_reason": reason,
                "geometric_pairwise_rmsd_json": json.dumps(pair_rows, sort_keys=True),
            }
        )

    return pd.DataFrame(rows, columns=_GEOMETRIC_COLUMNS)
