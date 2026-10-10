"""
Single pose selector shared by every best-pose writer (Spec 036 R4, D4).

The rules are the approved ``consensus_rank_geometry_qc_v2`` pose-selection rules
(``specs/031-appraisal-remediation-p0/decisions/consensus-v2.md``):

- GNINA: highest ``cnn_score`` selects the representative pose. Ties use higher
  ``cnn_affinity``, then lower Vina affinity (``affinity_kcal_mol``), then lower pose
  index, then the stable tag. A GNINA row without a finite ``cnn_score`` is not
  replaced by a Vina-affinity pick; it is kept with ``missing_required_cnn_score``.
- Vina / Smina / AutoDock4: lowest ``affinity_kcal_mol`` selects the pose. Ties use lower
  pose index, then the stable tag. Missing affinity is kept with
  ``missing_required_affinity``.
- Replicates (Spec 036 R5b): rows of every replicate of one (engine, protein, tag) are pooled.
  After the v2 keys, a final tie-break on ``replicate_id`` (lower first) makes the choice
  independent of row order. It only matters when two replicates tie exactly on every v2 key.
  The chosen row keeps its ``replicate_id``, ``seed`` and ``pose_file``.

Writers must not re-implement these rules (a ``min(cnn_affinity)`` or ``min(affinity)``
pick selects the worst GNINA pose). Use ``select_best_pose_rows``.
"""
from __future__ import annotations

from typing import List

import numpy as np
import pandas as pd

POSE_SELECTION_RULE_ID = "consensus_rank_geometry_qc_v2/pose_selection"
GROUP_COLUMNS = ("engine", "protein", "tag")
_POSE_SELECTION_NUMERIC = ("affinity_kcal_mol", "cnn_score", "cnn_affinity", "pose")


def select_best_pose_rows(scores: pd.DataFrame) -> pd.DataFrame:
    """Select one pose per (engine, protein, tag) using the approved v2 metrics.

    Returns one row per group, sorted by group key. Each row carries the source columns
    (numeric columns coerced) plus ``v2_pose_selection_status``, ``v2_ranking_metric``,
    ``v2_ranking_metric_name`` and ``v2_ranking_direction``.
    """
    if scores is None or scores.empty:
        return pd.DataFrame(columns=list(scores.columns) if scores is not None else [])
    frame = scores.copy()
    for column in ("engine", "protein", "tag", "ligand", "site_id"):
        if column not in frame.columns:
            frame[column] = ""
        frame[column] = frame[column].astype(str)
    if "replicate_id" not in frame.columns:
        frame["replicate_id"] = np.nan
    frame["_rep_order"] = pd.to_numeric(frame["replicate_id"], errors="coerce").fillna(np.inf)
    frame["engine"] = frame["engine"].str.strip().str.lower()
    for column in _POSE_SELECTION_NUMERIC:
        frame[column] = pd.to_numeric(frame.get(column), errors="coerce")

    selected: List[pd.Series] = []
    for (_, _, _), group in frame.groupby(list(GROUP_COLUMNS), dropna=False, sort=True):
        engine = str(group["engine"].iloc[0])
        working = group.copy()
        if engine == "gnina":
            working = working[np.isfinite(working["cnn_score"])].copy()
            if working.empty:
                # Preserve one auditable incomplete row instead of substituting affinity.
                row = group.sort_values(["tag"], kind="mergesort").iloc[0].copy()
                row["v2_pose_selection_status"] = "missing_required_cnn_score"
                row["v2_ranking_metric"] = np.nan
                row["v2_ranking_metric_name"] = "cnn_affinity"
                row["v2_ranking_direction"] = "higher_is_better"
                selected.append(row)
                continue
            working["_pose_order"] = working["pose"].fillna(np.inf)
            working.sort_values(
                ["cnn_score", "cnn_affinity", "affinity_kcal_mol", "_pose_order", "tag", "_rep_order"],
                ascending=[False, False, True, True, True, True],
                kind="mergesort",
                inplace=True,
            )
            row = working.iloc[0].copy()
            row["v2_ranking_metric"] = row.get("cnn_affinity")
            row["v2_ranking_metric_name"] = "cnn_affinity"
            row["v2_ranking_direction"] = "higher_is_better"
        else:
            working = working[np.isfinite(working["affinity_kcal_mol"])].copy()
            if working.empty:
                row = group.sort_values(["tag"], kind="mergesort").iloc[0].copy()
                row["v2_pose_selection_status"] = "missing_required_affinity"
                row["v2_ranking_metric"] = np.nan
                row["v2_ranking_metric_name"] = "affinity_kcal_mol"
                row["v2_ranking_direction"] = "lower_is_better"
                selected.append(row)
                continue
            working["_pose_order"] = working["pose"].fillna(np.inf)
            working.sort_values(
                ["affinity_kcal_mol", "_pose_order", "tag", "_rep_order"],
                ascending=[True, True, True, True],
                kind="mergesort",
                inplace=True,
            )
            row = working.iloc[0].copy()
            row["v2_ranking_metric"] = row.get("affinity_kcal_mol")
            row["v2_ranking_metric_name"] = "affinity_kcal_mol"
            row["v2_ranking_direction"] = "lower_is_better"
        row["v2_pose_selection_status"] = "completed"
        selected.append(row)
    return (
        pd.DataFrame(selected)
        .drop(columns=["_pose_order", "_rep_order"], errors="ignore")
        .reset_index(drop=True)
    )

