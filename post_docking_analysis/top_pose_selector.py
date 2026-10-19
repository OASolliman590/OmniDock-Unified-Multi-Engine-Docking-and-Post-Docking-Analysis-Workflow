"""
Top-pose ligand performance atlas utilities.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from post_docking_analysis.consensus import (
    CONSENSUS_V2_METHOD,
    SINGLE_ENGINE_NATIVE_METHOD,
    build_consensus_rankings,
)
from post_docking_analysis.pose_selection import select_best_pose_rows


SUPPORTED_SELECTION_POLICIES = {
    "best_affinity",
    "best_consensus",
    "hybrid",
}
# Spec 037 R4b: the default ranks by the v2 consensus rank (single_engine_native_v1 for one engine).
DEFAULT_SELECTION_POLICY = "best_consensus"
# Columns whose values are raw affinities. They are descriptive only and never compared across engines.
DESCRIPTIVE_NOT_COMPARABLE_COLUMNS = (
    "ligand_performance_summary.descriptive_best_affinity_kcal_mol",
    "top_pose_per_ligand_per_protein.best_affinity_kcal_mol",
)
DESCRIPTIVE_SCOPE = "descriptive_not_comparable_across_engines"

SUPPORTED_GLOBAL_AGGREGATIONS = {
    "best_target",
}

# Spec 036 R3: only v2 and the single-engine native rule remain. Both are higher-is-better.
_CONSENSUS_SCORE_HIGHER_IS_BETTER_MODES = frozenset(
    {
        "consensus_rank_geometry_qc_v2",
        "single_engine_native_v1",
    }
)


def normalize_selection_policy(value: Optional[str]) -> str:
    token = str(value or "").strip().lower()
    return token if token in SUPPORTED_SELECTION_POLICIES else DEFAULT_SELECTION_POLICY


def normalize_global_aggregation(value: Optional[str]) -> str:
    token = str(value or "").strip().lower()
    return token if token in SUPPORTED_GLOBAL_AGGREGATIONS else "best_target"


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _ensure_columns(frame: pd.DataFrame, columns: Iterable[str], default: object = None) -> pd.DataFrame:
    working = frame.copy()
    for column in columns:
        if column not in working.columns:
            working[column] = default
    return working


def _policy_sort_fields(selection_policy: str) -> Tuple[str, str, str]:
    """Primary and secondary sort keys. The secondary key is always the v2 rank (Spec 037 R4b), never a raw affinity.

    ``best_affinity`` sorts by the raw affinity and is therefore only valid within one engine (checked in
    ``build_top_pose_atlas``).
    """
    normalized = normalize_selection_policy(selection_policy)
    if normalized == "best_affinity":
        return (
            "affinity_kcal_mol",
            "consensus_rank_global",
            "affinity_kcal_mol -> consensus_rank_global -> tag",
        )
    return (
        "consensus_score",
        "consensus_rank_global",
        "consensus_score -> consensus_rank_global -> tag",
    )


def _ranking_mode_for_engines(engines: List[str]) -> str:
    """Spec 037 R3b/R4b: single_engine_native_v1 for one engine, consensus_rank_geometry_qc_v2 for two or more."""
    return SINGLE_ENGINE_NATIVE_METHOD if len(engines) == 1 else CONSENSUS_V2_METHOD


def _identity_rows(best_frame: pd.DataFrame, consensus_work: pd.DataFrame) -> pd.DataFrame:
    """One pose per ranked tag: the pose of its consensus winner engine (Spec 037 R3b, R4b).

    A tag without a complete v2 case (or without consensus data) keeps every engine row. No winner is invented,
    and the row carries its consensus status (see ``_best_by_tag`` in multi_engine_pipeline_impl).
    """
    winners: Dict[str, str] = {}
    if not consensus_work.empty and {"consensus_status", "winner_engine"}.issubset(consensus_work.columns):
        status = consensus_work["consensus_status"].astype(str)
        winner = consensus_work["winner_engine"].fillna("").astype(str).str.lower()
        ranked = consensus_work[status.isin(["completed", "single_engine"]) & (winner.str.len() > 0)]
        winners = dict(zip(ranked["tag"].astype(str), winner.loc[ranked.index]))
    if not winners:
        return best_frame.copy()
    winner_for_row = best_frame["tag"].map(winners)
    keep = winner_for_row.isna() | (best_frame["engine"] == winner_for_row)
    return best_frame[keep].copy()


def _consensus_score_ascending(mode: str) -> bool:
    """False for the v2 and single-engine modes (higher consensus_score is better). A blank mode means the v2 default."""
    token = str(mode or "").strip().lower() or CONSENSUS_V2_METHOD
    if token in _CONSENSUS_SCORE_HIGHER_IS_BETTER_MODES:
        return False
    return True


def _build_empty_contract() -> Dict[str, object]:
    return {
        "top_pose_per_ligand_per_protein": pd.DataFrame(),
        "top_pose_per_ligand_global": pd.DataFrame(),
        "ligand_performance_summary": pd.DataFrame(),
        "top_pose_confidence_metrics": pd.DataFrame(),
        "manifest": {
            "selection_policy": DEFAULT_SELECTION_POLICY,
            "global_aggregation": "best_target",
            "tie_break_rule": "consensus_score -> consensus_rank_global -> tag",
            "comparison_scope": DESCRIPTIVE_SCOPE,
            "descriptive_not_comparable_across_engines": list(DESCRIPTIVE_NOT_COMPARABLE_COLUMNS),
            "generated_at": _utc_now_iso(),
            "rows": {
                "top_pose_per_ligand_per_protein": 0,
                "top_pose_per_ligand_global": 0,
                "ligand_performance_summary": 0,
                "top_pose_confidence_metrics": 0,
            },
        },
    }


def build_top_pose_atlas(
    best_by_engine: pd.DataFrame,
    consensus_df: Optional[pd.DataFrame] = None,
    *,
    selection_policy: str = DEFAULT_SELECTION_POLICY,
    consensus_mode: str = "",
    global_aggregation: str = "best_target",
    run_id: str = "",
    generated_at: Optional[str] = None,
    context: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    """
    Build deterministic top-pose ligand performance tables.
    """
    if best_by_engine is None or best_by_engine.empty:
        payload = _build_empty_contract()
        payload["manifest"]["selection_policy"] = normalize_selection_policy(selection_policy)
        payload["manifest"]["global_aggregation"] = normalize_global_aggregation(global_aggregation)
        payload["manifest"]["run_id"] = str(run_id or "")
        payload["manifest"]["generated_at"] = str(generated_at or _utc_now_iso())
        payload["manifest"]["context"] = context or {}
        return payload

    policy = normalize_selection_policy(selection_policy)
    aggregation = normalize_global_aggregation(global_aggregation)
    primary_field, secondary_field, tie_break_rule = _policy_sort_fields(policy)
    score_ascending = _consensus_score_ascending(consensus_mode)
    primary_sort_ascending = score_ascending if primary_field == "consensus_score" else True

    best_frame = best_by_engine.copy()
    best_frame = _ensure_columns(
        best_frame,
        [
            "tag",
            "protein",
            "ligand",
            "site_id",
            "engine",
            "pose",
            "pose_file",
            "log_file",
            "affinity_kcal_mol",
            "normalized_affinity_score",
        ],
    )
    best_frame["tag"] = best_frame["tag"].astype(str)
    best_frame["protein"] = best_frame["protein"].astype(str)
    best_frame["ligand"] = best_frame["ligand"].astype(str)
    best_frame["site_id"] = best_frame["site_id"].astype(str)
    best_frame["engine"] = best_frame["engine"].astype(str).str.lower()
    best_frame["pose"] = pd.to_numeric(best_frame["pose"], errors="coerce").fillna(0).astype(int)
    best_frame["affinity_kcal_mol"] = pd.to_numeric(best_frame["affinity_kcal_mol"], errors="coerce")
    best_frame["normalized_affinity_score"] = pd.to_numeric(
        best_frame["normalized_affinity_score"],
        errors="coerce",
    )
    # Spec 036 R4: one pose per (engine, protein, tag) by the shared v2 selector, so that a GNINA
    # tag is represented by its highest-cnn_score pose and not by its lowest-affinity pose.
    best_frame = select_best_pose_rows(best_frame)
    best_frame["pose"] = pd.to_numeric(best_frame["pose"], errors="coerce").fillna(0).astype(int)

    # Spec 037 R4b: a raw affinity is compared only within one engine.
    engines = sorted(set(best_frame["engine"].dropna().astype(str)))
    if policy == "best_affinity" and len(engines) > 1:
        raise ValueError(
            "selection_policy 'best_affinity' compares raw affinities across engines, which Spec 037 R4b "
            f"removed (engines found: {', '.join(engines)}). Use 'best_consensus' (the default, consensus rank "
            "consensus_rank_geometry_qc_v2) or 'hybrid'. 'best_affinity' is available for a single engine only."
        )

    # Spec 037 R4b: without a consensus table the rank is computed from the same pose rows (never from affinities).
    if consensus_df is None or consensus_df.empty:
        consensus_df = build_consensus_rankings(
            best_by_engine=best_frame,
            consensus_mode=_ranking_mode_for_engines(engines),
        )

    engine_support = (
        best_frame.groupby("tag", dropna=False)["engine"]
        .nunique()
        .rename("engine_support_count")
        .reset_index()
    )

    consensus_work = pd.DataFrame(columns=["tag"])
    bio_columns: List[str] = []
    if consensus_df is not None and not consensus_df.empty:
        consensus_work = consensus_df.copy()
        consensus_work["tag"] = consensus_work["tag"].astype(str)
        bio_columns = [column for column in consensus_work.columns if str(column).startswith("bio_")]
        keep_columns = [
            "tag",
            "consensus_score",
            "consensus_rank_global",
            "consensus_status",
            "agreement_count",
            "agreement_fraction",
            "winner_engine",
            "best_affinity_kcal_mol",
            "mean_normalized_affinity_score",
            "docking_quality_class",
            "qc_status",
            "admet_status",
        ] + bio_columns
        keep_columns = [column for column in keep_columns if column in consensus_work.columns]
        consensus_work = consensus_work[keep_columns].drop_duplicates("tag", keep="first")

    # Spec 037 R4b: the pose of each ranked tag is its consensus winner's (no raw affinity picks the engine).
    tag_identity = _identity_rows(best_frame, consensus_work)
    atlas = tag_identity.merge(consensus_work, on="tag", how="left")
    atlas = atlas.merge(engine_support, on="tag", how="left")
    atlas["engine_support_count"] = pd.to_numeric(atlas["engine_support_count"], errors="coerce").fillna(
        pd.to_numeric(atlas.get("agreement_count"), errors="coerce")
    )
    atlas["engine_support_count"] = pd.to_numeric(atlas["engine_support_count"], errors="coerce").fillna(0).astype(int)
    atlas["single_engine_mode"] = atlas["engine_support_count"].le(1)
    atlas["consensus_score"] = pd.to_numeric(atlas.get("consensus_score"), errors="coerce")
    atlas["consensus_rank_global"] = pd.to_numeric(atlas.get("consensus_rank_global"), errors="coerce")
    atlas["agreement_fraction"] = pd.to_numeric(atlas.get("agreement_fraction"), errors="coerce")
    atlas["best_affinity_kcal_mol"] = pd.to_numeric(atlas.get("best_affinity_kcal_mol"), errors="coerce")
    atlas["mean_normalized_affinity_score"] = pd.to_numeric(
        atlas.get("mean_normalized_affinity_score"),
        errors="coerce",
    )
    atlas["affinity_kcal_mol"] = pd.to_numeric(atlas.get("affinity_kcal_mol"), errors="coerce")
    atlas = _ensure_columns(
        atlas,
        ["docking_quality_class", "qc_status", "admet_status", "winner_engine"],
        default="",
    )
    atlas["winner_engine"] = atlas["winner_engine"].astype(str).str.lower()
    atlas["winner_engine"] = atlas["winner_engine"].where(atlas["winner_engine"].str.len() > 0, atlas["engine"])

    missing_primary = np.inf if primary_sort_ascending else -np.inf
    atlas["_sort_primary"] = pd.to_numeric(atlas.get(primary_field), errors="coerce").fillna(missing_primary)
    atlas["_sort_secondary"] = pd.to_numeric(atlas.get(secondary_field), errors="coerce").fillna(np.inf)
    atlas["_sort_tag"] = atlas["tag"].astype(str)
    atlas["_sort_engine"] = atlas["engine"].astype(str)

    grouped = atlas.sort_values(
        ["ligand", "protein", "_sort_primary", "_sort_secondary", "_sort_tag", "_sort_engine"],
        ascending=[True, True, primary_sort_ascending, True, True, True],
    )
    per_protein = grouped.groupby(["ligand", "protein"], dropna=False, sort=False).head(1).copy()
    per_protein.sort_values(["ligand", "protein"], inplace=True)
    per_protein["selection_policy"] = policy
    per_protein["global_aggregation"] = aggregation
    per_protein["tie_break_rule"] = tie_break_rule

    global_sorted = per_protein.sort_values(
        ["ligand", "_sort_primary", "_sort_secondary", "_sort_tag", "_sort_engine"],
        ascending=[True, primary_sort_ascending, True, True, True],
    )
    global_table = global_sorted.groupby("ligand", dropna=False, sort=False).head(1).copy()
    protein_counts = per_protein.groupby("ligand", dropna=False)["protein"].nunique()
    global_table["protein_count_evaluated"] = global_table["ligand"].map(protein_counts).fillna(0).astype(int)
    global_table["selected_protein"] = global_table["protein"]
    global_table["selection_policy"] = policy
    global_table["global_aggregation"] = aggregation
    global_table["tie_break_rule"] = tie_break_rule
    global_table.sort_values(
        ["_sort_primary", "_sort_secondary", "ligand", "_sort_engine"],
        ascending=[primary_sort_ascending, True, True, True],
        inplace=True,
    )

    summary = (
        per_protein.groupby("ligand", dropna=False)
        .agg(
            proteins_evaluated=("protein", "nunique"),
            # Spec 037 R4b: min over the proteins' winners, which can be different engines. Descriptive only.
            descriptive_best_affinity_kcal_mol=("affinity_kcal_mol", "min"),
            best_consensus_score=("consensus_score", "min"),
            strong_hits=("docking_quality_class", lambda s: int((s.astype(str) == "Strong").sum())),
            moderate_hits=("docking_quality_class", lambda s: int((s.astype(str) == "Moderate").sum())),
            weak_hits=("docking_quality_class", lambda s: int((s.astype(str) == "Weak").sum())),
            qc_fail_count=("qc_status", lambda s: int(s.astype(str).str.startswith("fail").sum())),
            qc_warn_count=("qc_status", lambda s: int(s.astype(str).str.startswith("warn").sum())),
            admet_flagged_count=("admet_status", lambda s: int((s.astype(str) == "flagged").sum())),
        )
        .reset_index()
    )
    best_consensus_series = (
        per_protein.groupby("ligand", dropna=False)["consensus_score"].min()
        if score_ascending
        else per_protein.groupby("ligand", dropna=False)["consensus_score"].max()
    )
    summary["best_consensus_score"] = summary["ligand"].map(best_consensus_series)
    best_global_subset = global_table[["ligand", "protein", "engine", "tag"]].rename(
        columns={
            "protein": "global_best_protein",
            "engine": "global_best_engine",
            "tag": "global_best_tag",
        }
    )
    summary = summary.merge(best_global_subset, on="ligand", how="left")
    summary["selection_policy"] = policy
    summary["global_aggregation"] = aggregation
    # Spec 037 R4b: ordered by the consensus rank (higher score is better under v2), never by a raw affinity.
    summary.sort_values(["best_consensus_score", "ligand"], ascending=[score_ascending, True], inplace=True, na_position="last")

    confidence_columns = [
        "ligand",
        "protein",
        "tag",
        "engine",
        "engine_support_count",
        "single_engine_mode",
        "agreement_fraction",
        "consensus_score",
        "affinity_kcal_mol",
        "docking_quality_class",
        "qc_status",
        "admet_status",
        "selection_policy",
        "tie_break_rule",
    ] + bio_columns
    confidence_columns = [column for column in confidence_columns if column in per_protein.columns]
    confidence = per_protein[confidence_columns].copy()

    for frame in (per_protein, global_table, summary, confidence):
        frame.drop(columns=["_sort_primary", "_sort_secondary", "_sort_tag", "_sort_engine"], inplace=True, errors="ignore")

    manifest = {
        "selection_policy": policy,
        "global_aggregation": aggregation,
        "tie_break_rule": tie_break_rule,
        "comparison_scope": DESCRIPTIVE_SCOPE,
        "descriptive_not_comparable_across_engines": list(DESCRIPTIVE_NOT_COMPARABLE_COLUMNS),
        "engines": engines,
        "generated_at": str(generated_at or _utc_now_iso()),
        "run_id": str(run_id or ""),
        "rows": {
            "top_pose_per_ligand_per_protein": int(len(per_protein)),
            "top_pose_per_ligand_global": int(len(global_table)),
            "ligand_performance_summary": int(len(summary)),
            "top_pose_confidence_metrics": int(len(confidence)),
        },
        "source_rows": {
            "best_by_engine": int(len(best_by_engine)),
            "consensus": int(len(consensus_df)) if consensus_df is not None else 0,
        },
        "columns": {
            "top_pose_per_ligand_per_protein": list(per_protein.columns),
            "top_pose_per_ligand_global": list(global_table.columns),
            "ligand_performance_summary": list(summary.columns),
            "top_pose_confidence_metrics": list(confidence.columns),
        },
        "context": context or {},
    }

    return {
        "top_pose_per_ligand_per_protein": per_protein,
        "top_pose_per_ligand_global": global_table,
        "ligand_performance_summary": summary,
        "top_pose_confidence_metrics": confidence,
        "manifest": manifest,
    }


def write_top_pose_atlas(atlas_payload: Dict[str, object], output_dir: Path) -> Dict[str, str]:
    output_root = Path(output_dir).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    per_protein = atlas_payload.get("top_pose_per_ligand_per_protein", pd.DataFrame())
    global_table = atlas_payload.get("top_pose_per_ligand_global", pd.DataFrame())
    summary = atlas_payload.get("ligand_performance_summary", pd.DataFrame())
    confidence = atlas_payload.get("top_pose_confidence_metrics", pd.DataFrame())
    manifest = dict(atlas_payload.get("manifest", {}) or {})

    per_protein_csv = output_root / "top_pose_per_ligand_per_protein.csv"
    global_csv = output_root / "top_pose_per_ligand_global.csv"
    summary_csv = output_root / "ligand_performance_summary.csv"
    confidence_csv = output_root / "top_pose_confidence_metrics.csv"
    manifest_json = output_root / "top_pose_selection_manifest.json"

    per_protein.to_csv(per_protein_csv, index=False)
    global_table.to_csv(global_csv, index=False)
    summary.to_csv(summary_csv, index=False)
    confidence.to_csv(confidence_csv, index=False)

    per_protein.to_json(output_root / "top_pose_per_ligand_per_protein.json", orient="records", indent=2)
    global_table.to_json(output_root / "top_pose_per_ligand_global.json", orient="records", indent=2)
    summary.to_json(output_root / "ligand_performance_summary.json", orient="records", indent=2)

    manifest["files"] = {
        "top_pose_per_ligand_per_protein_csv": str(per_protein_csv),
        "top_pose_per_ligand_global_csv": str(global_csv),
        "ligand_performance_summary_csv": str(summary_csv),
        "top_pose_confidence_metrics_csv": str(confidence_csv),
    }
    manifest_json.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    return {
        "top_pose_per_ligand_per_protein_file": str(per_protein_csv),
        "top_pose_per_ligand_global_file": str(global_csv),
        "ligand_performance_summary_file": str(summary_csv),
        "top_pose_confidence_metrics_file": str(confidence_csv),
        "top_pose_selection_manifest_file": str(manifest_json),
    }
