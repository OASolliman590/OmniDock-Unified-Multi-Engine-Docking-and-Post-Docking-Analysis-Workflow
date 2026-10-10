"""
DockBox-style geometric consensus helpers.

Pose-file parsers and the Spec 031 graph-based geometric QC (consensus v2).

Spec 036 R3d: the legacy element-order pairwise RMSD (`_pairwise_rmsd`, centroid soft
alignment, `compute_geometric_consensus`) is retired. Geometry is computed only by
`compute_geometric_consensus_v2`, which uses the graph/symmetry mapping in atom_mapping.py.
"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from post_docking_analysis.atom_mapping import (
    MAPPING_METHOD,
    attach_coordinates_to_topology,
    compare_graph_poses,
    load_sdf_graph_pose,
    not_comparable_result,
)


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
    """Extract one 1-based MODEL ordinal from a PDB/PDBQT pose file.

    Files without MODEL records are treated as a single-pose file. A final
    MODEL that reaches EOF without ENDMDL is accepted, as produced by some
    docking tools. Other malformed MODEL boundaries and out-of-range pose
    requests are reported explicitly; the parser never substitutes model 1.
    """
    try:
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception as exc:
        return np.empty((0, 3), dtype=float), [], f"read_error:{exc}"

    try:
        numeric_pose = float(pose_index)
        selected_pose = int(numeric_pose)
    except (TypeError, ValueError, OverflowError):
        return np.empty((0, 3), dtype=float), [], "invalid_pose_index"
    if not np.isfinite(numeric_pose) or numeric_pose != selected_pose or selected_pose < 1:
        return np.empty((0, 3), dtype=float), [], "invalid_pose_index"

    models: List[List[str]] = []
    current: Optional[List[str]] = None
    saw_model = False
    base_atoms: List[str] = []

    for line in lines:
        if line.startswith("MODEL"):
            if current is not None:
                return np.empty((0, 3), dtype=float), [], "malformed_model:nested_model"
            if base_atoms:
                return np.empty((0, 3), dtype=float), [], "malformed_model:atom_outside_model"
            saw_model = True
            current = []
            continue
        if line.startswith("ENDMDL"):
            if current is None:
                return np.empty((0, 3), dtype=float), [], "malformed_model:unexpected_endmdl"
            models.append(current)
            current = None
            continue
        if line.startswith(("ATOM", "HETATM")):
            if current is not None:
                current.append(line)
            elif saw_model:
                return np.empty((0, 3), dtype=float), [], "malformed_model:atom_outside_model"
            else:
                base_atoms.append(line)

    # Accept a final unterminated model at EOF. This is the only missing
    # ENDMDL case that can be interpreted unambiguously.
    if current is not None:
        models.append(current)
    if not saw_model:
        models = [base_atoms]

    if not models:
        return np.empty((0, 3), dtype=float), [], "no_models"

    idx = selected_pose - 1
    if idx >= len(models):
        return (
            np.empty((0, 3), dtype=float),
            [],
            f"pose_index_out_of_range:{selected_pose}>{len(models)}",
        )
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


def _strict_graph_pose(row: pd.Series):
    pose_path = Path(str(row.get("pose_file") or "")).expanduser()
    if not pose_path.is_file():
        return None, "missing_pose_file"
    try:
        pose_index = int(row.get("pose", 1) or 1)
    except (TypeError, ValueError, OverflowError):
        return None, "invalid_pose_index"
    try:
        if pose_path.suffix.lower() in {".sdf", ".mol"}:
            return load_sdf_graph_pose(pose_path, pose_index=pose_index), ""
        topology_value = row.get("topology_file", row.get("ligand_topology_file", ""))
        topology_path = Path(str(topology_value or "")).expanduser()
        if not topology_path.is_file():
            return None, "missing_explicit_topology"
        coords, elements, error = _extract_pose_signature(pose_path, pose_index)
        if error or coords.size == 0:
            return None, f"pose_parse_failed:{error or 'no_heavy_atoms'}"
        return attach_coordinates_to_topology(coords, elements, topology_path), ""
    except (OSError, RuntimeError, ValueError) as exc:
        return None, f"topology_mapping_failed:{exc}"


def compute_geometric_consensus_v2(
    selected_by_engine: pd.DataFrame,
    *,
    rmsd_cutoff: float = 2.0,
) -> pd.DataFrame:
    """Graph/symmetry-aware QC; never uses legacy centroid alignment."""
    required = {"engine", "protein", "tag", "pose_file", "pose"}
    if selected_by_engine is None or selected_by_engine.empty or not required.issubset(selected_by_engine.columns):
        return pd.DataFrame()
    records: List[Dict[str, object]] = []
    for (protein, tag), group in selected_by_engine.groupby(["protein", "tag"], dropna=False, sort=True):
        poses: Dict[str, object] = {}
        load_errors: Dict[str, str] = {}
        for _, row in group.sort_values("engine", kind="mergesort").iterrows():
            engine = str(row.get("engine", "")).strip().lower()
            graph_pose, error = _strict_graph_pose(row)
            if graph_pose is None:
                load_errors[engine] = error
            else:
                poses[engine] = graph_pose
        pairs: List[Dict[str, object]] = []
        for engine_a, engine_b in combinations(sorted(set(group["engine"].astype(str).str.lower())), 2):
            if engine_a not in poses or engine_b not in poses:
                result = not_comparable_result(
                    load_errors.get(engine_a) or load_errors.get(engine_b) or "missing_graph_pose"
                )
            else:
                result = compare_graph_poses(poses[engine_a], poses[engine_b])
            pair = result.to_dict()
            pair.update({"engine_a": engine_a, "engine_b": engine_b, "cutoff_angstrom": float(rmsd_cutoff)})
            pair["within_cutoff"] = bool(result.comparable and float(result.rmsd_angstrom) <= rmsd_cutoff)
            pairs.append(pair)
        comparable_pairs = [pair for pair in pairs if pair["status"] == "comparable"]
        expected_pairs = len(pairs)
        fully_comparable = expected_pairs > 0 and len(comparable_pairs) == expected_pairs
        agreement = fully_comparable and all(bool(pair["within_cutoff"]) for pair in comparable_pairs)
        numeric = [float(pair["rmsd_angstrom"]) for pair in comparable_pairs]
        records.append({
            "protein": str(protein),
            "tag": str(tag),
            "geometry_status": "comparable" if fully_comparable else "not_comparable",
            "geometric_agreement": bool(agreement),
            "geometric_engine_count": int(group["engine"].nunique()),
            "geometric_expected_pairs": expected_pairs,
            "geometric_valid_pairs": len(comparable_pairs),
            "geometric_pass_pairs": sum(bool(pair["within_cutoff"]) for pair in comparable_pairs),
            "geometric_pairwise_rmsd_min": min(numeric) if numeric else np.nan,
            "geometric_pairwise_rmsd_mean": float(np.mean(numeric)) if numeric else np.nan,
            "geometric_pairwise_rmsd_max": max(numeric) if numeric else np.nan,
            "geometric_cutoff_angstrom": float(rmsd_cutoff),
            "geometric_reason": (
                "all_pairs_within_cutoff" if agreement
                else "rmsd_above_cutoff" if fully_comparable
                else "pairwise_mapping_incomplete"
            ),
            "geometric_pairwise_rmsd_json": json.dumps(pairs, sort_keys=True),
            "geometric_mapping_method": MAPPING_METHOD,
        })
    return pd.DataFrame(records)
