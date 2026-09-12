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


SUPPORTED_SELECTION_POLICIES = {
    "best_affinity",
    "best_consensus",
    "hybrid",
}

SUPPORTED_GLOBAL_AGGREGATIONS = {
    "best_target",
}

_CONSENSUS_SCORE_ASCENDING_MODES = frozenset(
    {
        "",
        "strict_consensus",
        "weighted_hybrid",
    }
)


def normalize_selection_policy(value: Optional[str]) -> str:
    token = str(value or "").strip().lower()
    return token if token in SUPPORTED_SELECTION_POLICIES else "best_affinity"


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
    normalized = normalize_selection_policy(selection_policy)
    if normalized == "best_affinity":
        return (
            "affinity_kcal_mol",
            "consensus_score",
            "affinity_kcal_mol -> consensus_score -> tag",
        )
    return (
        "consensus_score",
        "affinity_kcal_mol",
        "consensus_score -> affinity_kcal_mol -> tag",
    )


def _consensus_score_ascending(mode: str) -> bool:
    # Every producer defines consensus_score as higher-is-better.
    return False


def _build_empty_contract() -> Dict[str, object]:
    return {
        "top_pose_per_ligand_per_protein": pd.DataFrame(),
        "top_pose_per_ligand_global": pd.DataFrame(),
        "ligand_performance_summary": pd.DataFrame(),
        "top_pose_confidence_metrics": pd.DataFrame(),
        "manifest": {
            "selection_policy": "best_affinity",
            "global_aggregation": "best_target",
            "tie_break_rule": "affinity_kcal_mol -> consensus_score -> tag",
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
    selection_policy: str = "best_affinity",
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

    # Compare only within engine/target/site score distributions. Raw energies
    # from different engines have no common calibration.
    if best_frame["normalized_affinity_score"].isna().all():
        from post_docking_analysis.consensus import normalize_engine_scores
        best_frame = normalize_engine_scores(best_frame, group_keys=["engine", "protein", "site_id"])
    multi_engine = best_frame["engine"].nunique() > 1
    effective_primary = "normalized_affinity_score" if multi_engine and policy == "best_affinity" else primary_field
    if effective_primary == "normalized_affinity_score":
        primary_sort_ascending = False
    tie_break_rule = f"{effective_primary} -> tag (raw scores never break cross-engine ties)"
    identity_sorted = best_frame.copy()
    identity_sorted["_affinity_sort"] = -identity_sorted["normalized_affinity_score"].fillna(-np.inf)
    if consensus_df is not None and not consensus_df.empty and "winner_engine" in consensus_df:
        winners = consensus_df.drop_duplicates("tag").set_index("tag")["winner_engine"]
        identity_sorted["_selected_winner"] = identity_sorted["engine"].eq(identity_sorted["tag"].map(winners))
        identity_sorted.loc[identity_sorted["_selected_winner"], "_affinity_sort"] = -np.inf
    identity_sorted.sort_values(
        ["tag", "_affinity_sort", "engine", "pose"],
        ascending=[True, True, True, True],
        inplace=True,
    )
    tag_identity = identity_sorted.groupby("tag", dropna=False, sort=False).head(1).copy()
    tag_identity.drop(columns=["_affinity_sort"], inplace=True, errors="ignore")

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

    atlas = tag_identity.merge(consensus_work, on="tag", how="left")
    atlas = atlas.merge(engine_support, on="tag", how="left")
    atlas["engine_support_count"] = pd.to_numeric(atlas["engine_support_count"], errors="coerce").fillna(
        pd.to_numeric(atlas.get("agreement_count"), errors="coerce")
    )
    atlas["engine_support_count"] = pd.to_numeric(atlas["engine_support_count"], errors="coerce").fillna(0).astype(int)
    atlas["single_engine_mode"] = atlas["engine_support_count"].le(1)
    atlas["consensus_score"] = pd.to_numeric(atlas.get("consensus_score"), errors="coerce")
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

    atlas["_sort_primary"] = pd.to_numeric(atlas.get(effective_primary), errors="coerce").fillna(np.inf if primary_sort_ascending else -np.inf)
    atlas["_sort_secondary"] = 0.0
    atlas["interpretation"] = "exploratory_relative_prioritization_not_potency_or_selectivity"
    atlas["_sort_tag"] = atlas["tag"].astype(str)

    grouped = atlas.sort_values(
        ["ligand", "protein", "_sort_primary", "_sort_secondary", "_sort_tag"],
        ascending=[True, True, primary_sort_ascending, True, True],
    )
    per_protein = grouped.groupby(["ligand", "protein"], dropna=False, sort=False).head(1).copy()
    per_protein.sort_values(["ligand", "protein"], inplace=True)
    per_protein["selection_policy"] = policy
    per_protein["global_aggregation"] = aggregation
    per_protein["tie_break_rule"] = tie_break_rule

    # Cross-target best energy is not an affinity/selectivity measurement. Select
    # an explicitly exploratory representative by relative within-target rank.
    per_protein["_global_relative"] = per_protein["consensus_score"].fillna(per_protein["normalized_affinity_score"])
    global_sorted = per_protein.sort_values(["ligand", "_global_relative", "_sort_tag"], ascending=[True, False, True])
    global_table = global_sorted.groupby("ligand", dropna=False, sort=False).head(1).copy()
    protein_counts = per_protein.groupby("ligand", dropna=False)["protein"].nunique()
    global_table["protein_count_evaluated"] = global_table["ligand"].map(protein_counts).fillna(0).astype(int)
    global_table["selected_protein"] = global_table["protein"]
    global_table["selection_policy"] = policy
    global_table["global_aggregation"] = aggregation
    global_table["tie_break_rule"] = tie_break_rule
    global_table.sort_values(
        ["_sort_primary", "_sort_secondary", "ligand"],
        ascending=[primary_sort_ascending, True, True],
        inplace=True,
    )

    summary = (
        per_protein.groupby("ligand", dropna=False)
        .agg(
            proteins_evaluated=("protein", "nunique"),
            best_affinity_kcal_mol=("affinity_kcal_mol", "min"),
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
    summary["best_affinity_kcal_mol"] = np.nan
    summary["interpretation"] = "exploratory_relative_prioritization_not_potency_or_selectivity"
    summary.sort_values(["ligand"], inplace=True)

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
        frame.drop(columns=["_sort_primary", "_sort_secondary", "_sort_tag", "_global_relative", "_selected_winner"], inplace=True, errors="ignore")

    manifest = {
        "selection_policy": policy,
        "effective_primary": effective_primary,
        "interpretation": "exploratory_relative_prioritization_not_potency_or_selectivity",
        "global_aggregation": aggregation,
        "tie_break_rule": tie_break_rule,
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
