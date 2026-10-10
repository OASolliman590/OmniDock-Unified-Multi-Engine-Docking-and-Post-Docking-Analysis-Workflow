"""Replicate handling for post-docking analysis (Spec 036 R5b, replicate_contract.md).

- Pose files of one pair are ``<pair tag>__rep<NN><ext>``. ``split_pose_stem`` returns the pair tag and
  the replicate id. A stem without the suffix is a legacy single pose (replicate id None).
- ``load_replicate_seeds`` reads the seed of each replicate from the engine run manifest.
- ``pose_reproducibility_table`` compares the top pose of every replicate of one (engine, protein, tag)
  pairwise with the Spec 031/034 in-place heavy-atom comparison (``compare_graph_poses_in_place``).
  Lineage poses come from the Meeko REMARK reader (PDBQT) or the SDF graph loader (GNINA).

The metric is a reproducibility descriptor of the docking runs. It is not a binding-mode validation.
"""
from __future__ import annotations

import json
import math
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .replicate_names import split_replicate_stem

from .atom_mapping import (
    LineageError,
    compare_graph_poses_in_place,
    load_pdbqt_lineage_pose,
    load_sdf_graph_pose,
)
from .pose_selection import select_best_pose_rows

POSE_REPRODUCIBILITY_METHOD = "spec031_in_place_heavy_atom_rmsd_v1"
POSE_REPRODUCIBILITY_THRESHOLD_ANGSTROM = 2.0
REPRODUCIBILITY_COLUMNS = [
    "engine",
    "protein",
    "tag",
    "pose_reproducibility_status",
    "pose_reproducibility_replicates",
    "pose_reproducibility_pairs",
    "pose_reproducibility_pairs_comparable",
    "pose_reproducibility_max_rmsd_angstrom",
    "pose_reproducibility_median_rmsd_angstrom",
    "pose_reproducibility_fraction_within_2A",
    "pose_reproducibility_method",
    "pose_reproducibility_reasons",
]


def split_pose_stem(stem: str) -> Tuple[str, Optional[int]]:
    """(pair tag, replicate id) of a pose or log file stem. Legacy stems give (stem, None)."""
    return split_replicate_stem(stem)


def _as_replicate_id(value: object) -> Optional[int]:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return int(numeric)


def load_replicate_seeds(run_manifest: Path) -> Dict[int, Optional[int]]:
    """Map replicate id -> seed from ``engines/<engine>/run_manifest.json``. Empty when unavailable."""
    path = Path(run_manifest)
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    seeds: Dict[int, Optional[int]] = {}
    for job in payload.get("jobs", []) or []:
        replicate_id = _as_replicate_id(job.get("replicate_id"))
        if replicate_id is None:
            continue
        seed = job.get("seed")
        seeds.setdefault(replicate_id, _as_replicate_id(seed) if seed is not None else None)
    effective_seeds = (payload.get("effective") or {}).get("seeds") or []
    for index, seed in enumerate(effective_seeds, start=1):
        seeds.setdefault(index, _as_replicate_id(seed) if seed is not None else None)
    return seeds


def _load_pose(pose_file: str, pose_index: int):
    """GraphPose for one pose, or (None, reason). PDBQT uses the Meeko lineage reader."""
    path = Path(str(pose_file or "")).expanduser()
    if not path.is_file():
        return None, "missing_pose_file"
    try:
        if path.suffix.lower() in {".sdf", ".mol"}:
            return load_sdf_graph_pose(path, pose_index=pose_index), ""
        if path.suffix.lower() == ".pdbqt":
            return load_pdbqt_lineage_pose(path, model_index=pose_index), ""
    except LineageError as exc:
        return None, f"lineage:{exc.reason}"
    except (OSError, RuntimeError, ValueError) as exc:
        return None, f"pose_load_failed:{exc}"
    return None, "unsupported_pose_format"


def _top_pose_per_replicate(group: pd.DataFrame) -> List[pd.Series]:
    """One selected row per replicate of one (engine, protein, tag), by the shared v2 selector.

    Rows without a replicate id (legacy single poses) form one replicate.
    """
    if group.empty:
        return []
    keys = [_as_replicate_id(value) for value in group.get("replicate_id", pd.Series([None] * len(group), index=group.index))]
    buckets: Dict[Optional[int], List[object]] = {}
    for index, key in zip(group.index, keys):
        buckets.setdefault(key, []).append(index)
    rows: List[pd.Series] = []
    for key in sorted(buckets, key=lambda value: (value is None, value or 0)):
        chosen = select_best_pose_rows(group.loc[buckets[key]])
        if not chosen.empty:
            rows.append(chosen.iloc[0])
    return rows


def _summarize(
    status: str,
    replicates: int,
    rmsd_values: List[float],
    pairs: int,
    reasons: List[str],
) -> Dict[str, object]:
    comparable = len(rmsd_values)
    if comparable:
        values = np.asarray(rmsd_values, dtype=float)
        max_value = float(np.max(values))
        median_value = float(np.median(values))
        fraction = float(np.mean(values <= POSE_REPRODUCIBILITY_THRESHOLD_ANGSTROM))
    else:
        max_value = median_value = fraction = float("nan")
    return {
        "pose_reproducibility_status": status,
        "pose_reproducibility_replicates": replicates,
        "pose_reproducibility_pairs": pairs,
        "pose_reproducibility_pairs_comparable": comparable,
        "pose_reproducibility_max_rmsd_angstrom": max_value,
        "pose_reproducibility_median_rmsd_angstrom": median_value,
        "pose_reproducibility_fraction_within_2A": fraction,
        "pose_reproducibility_method": POSE_REPRODUCIBILITY_METHOD,
        "pose_reproducibility_reasons": ";".join(sorted({reason for reason in reasons if reason})),
    }


def pose_reproducibility_table(normalized: pd.DataFrame) -> pd.DataFrame:
    """One row per (engine, protein, tag) with the replicate pose-reproducibility metric.

    ``normalized`` holds the pooled rows of every replicate (``replicate_id`` and ``pose_file`` columns).
    A group with one replicate is ``single_replicate`` and carries no RMSD. Pairwise values use the
    replicate's own top pose, selected by the shared v2 selector.
    """
    if normalized is None or normalized.empty:
        return pd.DataFrame(columns=REPRODUCIBILITY_COLUMNS)
    frame = normalized.copy()
    for column in ("engine", "protein", "tag"):
        if column not in frame.columns:
            frame[column] = ""
        frame[column] = frame[column].astype(str)
    if "replicate_id" not in frame.columns:
        frame["replicate_id"] = np.nan

    pose_cache: Dict[Tuple[str, int], object] = {}
    out: List[Dict[str, object]] = []
    for (engine, protein, tag), group in frame.groupby(["engine", "protein", "tag"], sort=True, dropna=False):
        tops = _top_pose_per_replicate(group)
        replicates = len(tops)
        if replicates == 0:
            record = _summarize("no_poses", 0, [], 0, [])
        elif replicates == 1:
            record = _summarize("single_replicate", 1, [], 0, [])
        else:
            rmsd_values: List[float] = []
            reasons: List[str] = []
            pairs = 0
            loaded: List[Tuple[object, object]] = []
            for top in tops:
                pose_value = pd.to_numeric(top.get("pose"), errors="coerce")
                key = (str(top.get("pose_file", "")), int(pose_value) if pd.notna(pose_value) else 1)
                if key not in pose_cache:
                    pose_cache[key] = _load_pose(key[0], key[1])
                loaded.append(pose_cache[key])
            for (pose_a, reason_a), (pose_b, reason_b) in combinations(loaded, 2):
                pairs += 1
                if pose_a is None or pose_b is None:
                    reasons.append(str(reason_a or reason_b))
                    continue
                result = compare_graph_poses_in_place(pose_a, pose_b)
                if result.comparable:
                    rmsd_values.append(float(result.rmsd_angstrom))
                else:
                    reasons.append(result.reason)
            if rmsd_values and len(rmsd_values) == pairs:
                status = "completed"
            elif rmsd_values:
                status = "partial_not_comparable"
            else:
                status = "not_comparable"
            record = _summarize(status, replicates, rmsd_values, pairs, reasons)
        record.update({"engine": engine, "protein": protein, "tag": tag})
        out.append(record)
    return pd.DataFrame(out, columns=REPRODUCIBILITY_COLUMNS)
