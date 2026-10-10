"""
Consensus ranking utilities for multi-engine post-docking analysis.
"""
from __future__ import annotations

import logging
import json
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from post_docking_analysis.geometric_consensus import compute_geometric_consensus_v2
from post_docking_analysis.pose_selection import select_best_pose_rows

logger = logging.getLogger(__name__)


CONSENSUS_V2_METHOD = "consensus_rank_geometry_qc_v2"
SINGLE_ENGINE_NATIVE_METHOD = "single_engine_native_v1"
DEFAULT_CONSENSUS_MODE = CONSENSUS_V2_METHOD

# Spec 036 R3a: these legacy modes were removed. They used a method forbidden by Spec 031
# (file-order / element-order matching, centroid soft alignment, mapping-free RMSD) or
# duplicated consensus v2 without a distinct purpose. Requests for them raise.
REMOVED_CONSENSUS_MODES = frozenset(
    {
        "dockbox_geometric",
        "weighted_hybrid",
        "strict_consensus",
        "favorite_guardrails",
    }
)
SUPPORTED_CONSENSUS_MODES = {
    CONSENSUS_V2_METHOD,
    SINGLE_ENGINE_NATIVE_METHOD,
}

SUPPORTED_RESCORING_SCOPES = {
    "top_n_per_protein",
    "top_n_global",
}

SUPPORTED_NORMALIZATION_METHODS = {
    "per_engine_rank",
    "per_engine_minmax",
    "per_engine_zscore",
}

SUPPORTED_HIT_CLASS_POLICIES = {
    "target_percentile",
    "reference_anchor",
}

def normalize_consensus_mode(value: Optional[str]) -> str:
    """Return a supported consensus mode; blank means the v2 default.

    Spec 036 R3: removed and unknown names raise. There is no silent fallback to a legacy mode.
    """
    token = str(value or "").strip().lower()
    if not token:
        return DEFAULT_CONSENSUS_MODE
    if token in SUPPORTED_CONSENSUS_MODES:
        return token
    if token in REMOVED_CONSENSUS_MODES:
        raise ValueError(
            f"consensus mode '{token}' was removed by Spec 036 R3 (it relied on a method forbidden by "
            f"Spec 031 or duplicated consensus v2). Use '{CONSENSUS_V2_METHOD}', which is the default."
        )
    raise ValueError(
        f"unknown consensus mode '{token}'. Use '{CONSENSUS_V2_METHOD}' (default). "
        f"'{SINGLE_ENGINE_NATIVE_METHOD}' is selected automatically for single-engine projects."
    )


def select_representative_poses_v2(scores: pd.DataFrame) -> pd.DataFrame:
    """Select one pose per engine/target/tag using the approved v2 metrics.

    Delegates to the shared selector (Spec 036 R4). Behaviour is unchanged.
    """
    return select_best_pose_rows(scores)


def _build_consensus_v2(scores: pd.DataFrame, requested_engines: Optional[List[str]] = None) -> pd.DataFrame:
    selected = select_representative_poses_v2(scores)
    if selected.empty:
        return selected
    geometry = compute_geometric_consensus_v2(selected, rmsd_cutoff=2.0)
    if not geometry.empty:
        selected = selected.drop(
            columns=[column for column in geometry.columns if column in selected.columns and column not in {"protein", "tag"}],
            errors="ignore",
        ).merge(geometry, on=["protein", "tag"], how="left")
    outputs: List[Dict[str, object]] = []
    requested = sorted({str(engine).strip().lower() for engine in (requested_engines or []) if str(engine).strip()})
    for protein, target in selected.groupby("protein", dropna=False, sort=True):
        eligible = sorted(
            str(engine)
            for engine, group in target.groupby("engine", dropna=False)
            if pd.to_numeric(group["v2_ranking_metric"], errors="coerce").notna().any()
        )
        requested_for_target = requested or sorted(target["engine"].astype(str).unique().tolist())
        absent_for_target = sorted(set(requested_for_target) - set(target["engine"].astype(str)))
        tag_engine = target.pivot_table(
            index="tag", columns="engine", values="v2_ranking_metric", aggfunc="first"
        )
        complete_tags = [
            str(tag) for tag in tag_engine.index
            if eligible and all(engine in tag_engine.columns and pd.notna(tag_engine.at[tag, engine]) for engine in eligible)
        ]
        percentiles: Dict[str, pd.Series] = {}
        for engine in eligible:
            values = pd.to_numeric(tag_engine.loc[complete_tags, engine], errors="coerce")
            direction = str(
                target.loc[target["engine"] == engine, "v2_ranking_direction"].iloc[0]
            )
            ranks = values.rank(method="average", ascending=(direction == "lower_is_better"))
            count = len(values)
            percentiles[engine] = (
                pd.Series(1.0, index=values.index, dtype=float)
                if count == 1
                else (float(count) - ranks) / float(count - 1)
            )
        for tag, tag_rows in target.groupby("tag", dropna=False, sort=True):
            first = tag_rows.sort_values(["engine", "tag"], kind="mergesort").iloc[0]
            complete = str(tag) in complete_tags
            components = {
                engine: float(percentiles[engine].loc[str(tag)])
                for engine in eligible if complete
            }
            geometry_status = str(first.get("geometry_status", first.get("geometric_reason", "not_comparable")))
            geometry_agreement = bool(first.get("geometric_agreement", False))
            finite_affinities = pd.to_numeric(tag_rows.get("affinity_kcal_mol"), errors="coerce")
            best_affinity = float(finite_affinities.min()) if finite_affinities.notna().any() else np.nan
            winner_engine = max(sorted(components), key=lambda key: components[key]) if components else ""
            outputs.append({
                "tag": str(tag),
                "protein": str(protein),
                "ligand": str(first.get("ligand", "")),
                "site_id": str(first.get("site_id", "")),
                "consensus_score": float(np.mean(list(components.values()))) if components else np.nan,
                "mean_rank_pct": float(np.mean(list(components.values()))) if components else np.nan,
                "mean_normalized_affinity_score": float(np.mean(list(components.values()))) if components else np.nan,
                "best_affinity_kcal_mol": best_affinity,
                "agreement_count": len(components),
                "agreement_fraction": 1.0 if components else np.nan,
                "single_engine": bool(complete and len(eligible) == 1),
                "winner_engine": winner_engine,
                "score_name_primary": "equal_weight_engine_rank_percentiles",
                "score_primary": float(np.mean(list(components.values()))) if components else np.nan,
                "score_name_secondary": "geometry_qc",
                "score_secondary": np.nan,
                "consensus_status": (
                    "single_engine" if complete and len(eligible) == 1
                    else "completed" if complete
                    else "consensus_incomplete"
                ),
                "consensus_mode": CONSENSUS_V2_METHOD,
                "consensus_direction": "higher_is_better",
                "eligible_engines": ",".join(eligible),
                "requested_engines": ",".join(requested_for_target),
                "target_absent_engines": ",".join(absent_for_target),
                "included_engines": ",".join(sorted(components)),
                "missing_engines": ",".join(sorted(set(eligible) - set(components))),
                "engine_count": len(eligible),
                "complete_case_ligand_count": len(complete_tags),
                "engine_percentiles_json": json.dumps(components, sort_keys=True),
                "engine_weights_json": json.dumps(
                    {engine: 1.0 / len(eligible) for engine in eligible} if eligible else {}, sort_keys=True
                ),
                "geometry_status": geometry_status,
                "geometric_agreement": geometry_agreement,
                "geometry_tie_key": 1 if geometry_status == "comparable" and geometry_agreement else 0,
                "missing_data_policy": "complete_case_no_imputation_no_renormalization",
            })
    result = pd.DataFrame(outputs)
    result["_score"] = pd.to_numeric(result["consensus_score"], errors="coerce").fillna(-np.inf)

    def _rank_target(group: pd.DataFrame) -> pd.DataFrame:
        ranked = group.copy()
        # Geometry only orders an exact numeric tie when every member of that tie is comparable.
        ranked["_geometry_order"] = 0
        for _, tie in ranked.groupby("_score", dropna=False):
            if len(tie) > 1 and (tie["geometry_status"] == "comparable").all():
                ranked.loc[tie.index, "_geometry_order"] = tie["geometry_tie_key"]
        return ranked.sort_values(
            ["_score", "_geometry_order", "tag"], ascending=[False, False, True], kind="mergesort"
        )

    result = pd.concat(
        [_rank_target(group) for _, group in result.groupby("protein", dropna=False, sort=True)],
        ignore_index=True,
    )
    result["consensus_rank_within_protein"] = result.groupby("protein", dropna=False).cumcount() + 1
    result.sort_values(["_score", "protein", "consensus_rank_within_protein"], ascending=[False, True, True], inplace=True)
    result["consensus_rank_global"] = np.arange(1, len(result) + 1, dtype=int)
    return result.drop(columns=["_score", "_geometry_order"], errors="ignore").reset_index(drop=True)


def normalize_rescoring_scope(value: Optional[str]) -> str:
    token = str(value or "").strip().lower()
    return token if token in SUPPORTED_RESCORING_SCOPES else "top_n_per_protein"


def normalize_top_n(value: Optional[int], default: int = 3) -> int:
    try:
        parsed = int(value if value is not None else default)
    except (TypeError, ValueError):
        parsed = default
    return max(parsed, 1)


def normalize_normalization_method(value: Optional[str]) -> str:
    token = str(value or "").strip().lower()
    return token if token in SUPPORTED_NORMALIZATION_METHODS else "per_engine_rank"


def normalize_hit_class_policy(value: Optional[str]) -> str:
    token = str(value or "").strip().lower()
    return token if token in SUPPORTED_HIT_CLASS_POLICIES else "target_percentile"


def _safe_rank_pct(series: pd.Series) -> pd.Series:
    ranked = pd.to_numeric(series, errors="coerce")
    if ranked.isna().all():
        return pd.Series(np.nan, index=series.index, dtype=float)
    return ranked.rank(pct=True, method="average")


def _safe_rank_score(series: pd.Series, *, lower_is_better: bool = True) -> pd.Series:
    """
    Return a deterministic [0,1] rank score where 1.0 is best and 0.0 is worst.
    """
    values = pd.to_numeric(series, errors="coerce")
    finite_mask = np.isfinite(values.to_numpy(dtype=float))
    if finite_mask.sum() <= 0:
        return pd.Series(np.nan, index=series.index, dtype=float)
    if finite_mask.sum() == 1:
        out = pd.Series(np.nan, index=series.index, dtype=float)
        out.loc[finite_mask] = 1.0
        return out
    valid = values[finite_mask]
    # pandas rank: ascending=False gives the lowest value the largest rank number.
    ascending = not lower_is_better
    ranks = valid.rank(method="average", ascending=ascending)
    scaled = (ranks - 1.0) / float(len(valid) - 1)
    out = pd.Series(np.nan, index=series.index, dtype=float)
    out.loc[valid.index] = scaled
    return out


def _safe_minmax_scale(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    finite_mask = np.isfinite(values.to_numpy(dtype=float))
    if finite_mask.sum() <= 1:
        return pd.Series(np.zeros(len(series), dtype=float), index=series.index)
    valid = values[finite_mask]
    span = float(valid.max() - valid.min())
    if span <= 0:
        return pd.Series(np.zeros(len(series), dtype=float), index=series.index)
    # Docking convention: lower affinity is better, so invert min-max.
    scaled = (float(valid.max()) - values) / span
    return scaled.fillna(0.0)


def _safe_zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    finite_mask = np.isfinite(values.to_numpy(dtype=float))
    if finite_mask.sum() <= 1:
        return pd.Series(np.zeros(len(series), dtype=float), index=series.index)
    valid = values[finite_mask]
    mean_value = float(valid.mean())
    std_value = float(valid.std(ddof=0))
    if std_value <= 0:
        return pd.Series(np.zeros(len(series), dtype=float), index=series.index)
    zscores = (values - mean_value) / std_value
    return zscores.fillna(0.0)


def normalize_engine_scores(
    frame: pd.DataFrame,
    *,
    method: str = "per_engine_rank",
    score_column: str = "affinity_kcal_mol",
    normalized_column: str = "normalized_affinity_score",
    group_keys: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Normalize scores inside each engine/protein group before cross-engine consensus.
    All supported methods return higher normalized values for stronger binders.
    """
    normalized_method = normalize_normalization_method(method)
    keys = list(group_keys or ["engine", "protein"])
    working = frame.copy()
    working[score_column] = pd.to_numeric(working.get(score_column), errors="coerce")

    def _normalize_group(series: pd.Series) -> pd.Series:
        if normalized_method == "per_engine_minmax":
            return _safe_minmax_scale(series)
        if normalized_method == "per_engine_zscore":
            zscores = _safe_zscore(series)
            # Docking z-score: lower is better; negate before [0,1] scaling.
            return _safe_minmax_scale(-zscores)
        # Default: rank-based normalization with 1.0 = strongest.
        return _safe_rank_score(series, lower_is_better=True)

    working[normalized_column] = (
        working.groupby(keys, dropna=False)[score_column].transform(_normalize_group)
    )
    working[normalized_column] = pd.to_numeric(working[normalized_column], errors="coerce")
    working["normalization_method"] = normalized_method
    return working


def _build_single_engine_native(best_by_engine: pd.DataFrame, requested_engines: Optional[List[str]] = None) -> pd.DataFrame:
    """Explicit single-engine ranking (Spec 036 R3c): that engine's native metric with v2 pose selection.

    Not multi-engine consensus. Each ligand's composite is its within-target percentile of the
    native metric (GNINA: cnn_affinity, higher is better; Vina/Smina/AD4: affinity, lower is better).
    """
    if best_by_engine is None or best_by_engine.empty or "engine" not in best_by_engine.columns:
        return _build_consensus_v2(best_by_engine, requested_engines=requested_engines)
    engines = sorted({str(value).strip().lower() for value in best_by_engine["engine"].dropna().astype(str)})
    if len(engines) > 1:
        raise ValueError(
            f"{SINGLE_ENGINE_NATIVE_METHOD} requires exactly one engine; found {engines}. "
            f"Use '{CONSENSUS_V2_METHOD}' for multi-engine projects."
        )
    result = _build_consensus_v2(best_by_engine, requested_engines=requested_engines or engines)
    if result.empty:
        return result
    result["consensus_mode"] = SINGLE_ENGINE_NATIVE_METHOD
    # No engine agreement exists for one engine; never report 1.0 agreement.
    result["agreement_fraction"] = np.nan
    # Name and value of the native metric that ranks the ligands (GNINA: cnn_affinity; others: affinity).
    selected = select_representative_poses_v2(best_by_engine)
    native = selected.drop_duplicates("tag").set_index("tag")
    result["score_name_primary"] = result["tag"].map(native["v2_ranking_metric_name"]).map(
        lambda name: "cnn_affinity" if name == "cnn_affinity" else "vina_affinity"
    )
    result["score_primary"] = result["tag"].map(native["v2_ranking_metric"]).astype(float)
    return result


def build_consensus_rankings(
    best_by_engine: pd.DataFrame,
    consensus_mode: str = DEFAULT_CONSENSUS_MODE,
    favorite_engine: Optional[str] = None,
    normalization_method: str = "per_engine_rank",
    score_column: str = "affinity_kcal_mol",
    normalized_column: str = "normalized_affinity_score",
    requested_engines: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Build the consensus ranking table from the per-(engine, tag) pose rows.

    Modes (Spec 036 R3):
    - ``consensus_rank_geometry_qc_v2`` (default): approved Spec 031 policy.
    - ``single_engine_native_v1``: one engine only, native metric with v2 pose selection.

    ``favorite_engine``, ``normalization_method``, ``score_column`` and ``normalized_column``
    are accepted for call compatibility. They do not change the v2 ranking.
    """
    mode = normalize_consensus_mode(consensus_mode)
    if mode == CONSENSUS_V2_METHOD:
        return _build_consensus_v2(best_by_engine, requested_engines=requested_engines)
    return _build_single_engine_native(best_by_engine, requested_engines=requested_engines)


def select_rescoring_candidates(
    consensus_df: pd.DataFrame,
    rescoring_scope: str = "top_n_per_protein",
    rescoring_top_n: int = 3,
) -> pd.DataFrame:
    scope = normalize_rescoring_scope(rescoring_scope)
    top_n = normalize_top_n(rescoring_top_n, default=3)

    if consensus_df is None or consensus_df.empty:
        return pd.DataFrame(columns=list(consensus_df.columns) if consensus_df is not None else [])

    frame = consensus_df.copy()
    frame = frame.sort_values(
        ["consensus_score", "best_affinity_kcal_mol", "tag"],
        ascending=[False, True, True],
    )
    if scope == "top_n_global":
        selected = frame.head(top_n).copy()
    else:
        selected = frame.groupby("protein", dropna=False).head(top_n).copy()

    selected["rescoring_scope"] = scope
    selected["rescoring_top_n"] = int(top_n)
    return selected


def build_consensus_explainability(
    consensus_df: pd.DataFrame,
    rescoring_df: pd.DataFrame,
    consensus_mode: str,
    rescoring_scope: str,
    rescoring_top_n: int,
    favorite_engine: Optional[str] = None,
    normalization_method: str = "per_engine_rank",
) -> Dict[str, object]:
    mode = normalize_consensus_mode(consensus_mode)
    scope = normalize_rescoring_scope(rescoring_scope)
    top_n = normalize_top_n(rescoring_top_n, default=3)

    if consensus_df is None or consensus_df.empty:
        return {
            "consensus_mode": mode,
            "normalization_method": normalize_normalization_method(normalization_method),
            "rescoring_scope": scope,
            "rescoring_top_n": top_n,
            "favorite_engine": str(favorite_engine or ""),
            "candidate_count": 0,
            "rescoring_candidate_count": 0,
            "proteins": {},
            "top_global": [],
        }

    proteins: Dict[str, Dict[str, object]] = {}
    for protein, group in consensus_df.groupby("protein", dropna=False):
        protein_key = str(protein)
        proteins[protein_key] = {
            "candidate_count": int(len(group)),
            "best_consensus_score": float(group["consensus_score"].max()),
            "best_affinity_kcal_mol": float(group["best_affinity_kcal_mol"].min()),
            "mean_normalized_affinity_score": float(
                pd.to_numeric(group.get("mean_normalized_affinity_score"), errors="coerce").dropna().mean()
            )
            if "mean_normalized_affinity_score" in group.columns
            else None,
            "top_tags": group.sort_values(
                ["consensus_score", "best_affinity_kcal_mol", "tag"],
                ascending=[False, True, True],
            )["tag"].astype(str).head(5).tolist(),
        }

    top_global_rows: List[Dict[str, object]] = []
    for _, row in (
        consensus_df.sort_values(
            ["consensus_score", "best_affinity_kcal_mol", "tag"],
            ascending=[False, True, True],
        )
        .head(20)
        .iterrows()
    ):
        top_global_rows.append(
            {
                "tag": str(row.get("tag", "")),
                "protein": str(row.get("protein", "")),
                "ligand": str(row.get("ligand", "")),
                "consensus_score": float(row.get("consensus_score", np.nan)),
                "agreement_count": int(row.get("agreement_count", 0)),
                "winner_engine": str(row.get("winner_engine", "")),
            }
        )

    return {
        "consensus_mode": mode,
        "normalization_method": normalize_normalization_method(normalization_method),
        "rescoring_scope": scope,
        "rescoring_top_n": top_n,
        "favorite_engine": str(favorite_engine or ""),
        "candidate_count": int(len(consensus_df)),
        "rescoring_candidate_count": int(len(rescoring_df) if rescoring_df is not None else 0),
        "hit_class_distribution": (
            {
                str(label): int(count)
                for label, count in consensus_df["docking_quality_class"].value_counts(dropna=False).items()
            }
            if "docking_quality_class" in consensus_df.columns
            else {}
        ),
        "qc_status_distribution": (
            {
                str(label): int(count)
                for label, count in consensus_df["qc_status"].value_counts(dropna=False).items()
            }
            if "qc_status" in consensus_df.columns
            else {}
        ),
        "proteins": proteins,
        "top_global": top_global_rows,
    }


def classify_hits_target_aware(
    consensus_df: pd.DataFrame,
    *,
    policy: str = "target_percentile",
    strong_percentile: float = 0.10,
    moderate_percentile: float = 0.35,
    reference_baselines: Optional[pd.DataFrame] = None,
    reference_delta_kcal_mol: float = 1.0,
    require_reference_pass: bool = True,
) -> pd.DataFrame:
    """
    Assign target-aware Strong/Moderate/Weak classes while keeping QC and ADMET separated.
    """
    normalized_policy = normalize_hit_class_policy(policy)
    if consensus_df is None or consensus_df.empty:
        return consensus_df.copy()

    strong_pct = float(strong_percentile)
    moderate_pct = float(moderate_percentile)
    if strong_pct <= 0 or strong_pct >= 1:
        strong_pct = 0.10
    if moderate_pct <= strong_pct or moderate_pct >= 1:
        moderate_pct = max(0.35, strong_pct + 0.10)
        moderate_pct = min(moderate_pct, 0.90)

    frame = consensus_df.copy()
    frame["consensus_score"] = pd.to_numeric(frame.get("consensus_score"), errors="coerce")
    frame["best_affinity_kcal_mol"] = pd.to_numeric(frame.get("best_affinity_kcal_mol"), errors="coerce")
    frame["agreement_count"] = pd.to_numeric(frame.get("agreement_count"), errors="coerce").fillna(0).astype(int)
    single_engine_series = frame.get("single_engine", False)
    if isinstance(single_engine_series, pd.Series):
        single_engine_series = single_engine_series.fillna(False).astype(bool)
    else:
        single_engine_series = pd.Series(bool(single_engine_series), index=frame.index, dtype=bool)
    frame["single_engine"] = single_engine_series

    frame["target_percentile"] = (
        frame.groupby("protein", dropna=False)["consensus_score"].rank(
            pct=True,
            method="average",
            ascending=False,
        )
    )
    frame["target_ligand_count"] = frame.groupby("protein", dropna=False)["protein"].transform("size").astype(int)

    frame["docking_quality_class"] = "Weak"
    frame.loc[frame["target_percentile"] <= moderate_pct, "docking_quality_class"] = "Moderate"
    frame.loc[frame["target_percentile"] <= strong_pct, "docking_quality_class"] = "Strong"
    frame["classification_basis"] = normalized_policy

    frame["reference_affinity"] = np.nan
    frame["reference_tag"] = ""
    frame["delta_vs_reference"] = np.nan
    frame["reference_anchor_available"] = False
    frame["reference_anchor_reason"] = "not_requested"

    if normalized_policy == "reference_anchor":
        baseline_df = pd.DataFrame(columns=["protein", "reference_affinity", "reference_tag", "redocking_classification"])
        if reference_baselines is not None and not reference_baselines.empty:
            baseline_df = reference_baselines.copy()
        for required in ("protein", "reference_affinity", "reference_tag", "redocking_classification"):
            if required not in baseline_df.columns:
                if required == "reference_affinity":
                    baseline_df[required] = np.nan
                else:
                    baseline_df[required] = ""
        baseline_df = baseline_df[["protein", "reference_affinity", "reference_tag", "redocking_classification"]].copy()
        baseline_df["reference_affinity"] = pd.to_numeric(baseline_df["reference_affinity"], errors="coerce")
        baseline_df["redocking_classification"] = baseline_df["redocking_classification"].fillna("").astype(str).str.lower()
        baseline_df = baseline_df.drop_duplicates("protein")

        frame = frame.merge(
            baseline_df,
            on="protein",
            how="left",
            suffixes=("", "_baseline"),
        )
        if "reference_affinity_baseline" in frame.columns:
            frame["reference_affinity"] = pd.to_numeric(frame.get("reference_affinity_baseline"), errors="coerce")
        frame["reference_tag"] = (
            frame["reference_tag_baseline"].fillna("").astype(str)
            if "reference_tag_baseline" in frame.columns
            else frame.get("reference_tag", "").astype(str)
        )
        if "redocking_classification_baseline" in frame.columns:
            redocking_class = frame["redocking_classification_baseline"].fillna("").astype(str).str.lower()
        elif "redocking_classification" in frame.columns:
            redocking_class = frame["redocking_classification"].fillna("").astype(str).str.lower()
        else:
            redocking_class = pd.Series("", index=frame.index, dtype=str)
        frame.drop(
            columns=[
                "reference_affinity_baseline",
                "reference_tag_baseline",
                "redocking_classification_baseline",
                "redocking_classification",
            ],
            inplace=True,
            errors="ignore",
        )
        frame["delta_vs_reference"] = frame["best_affinity_kcal_mol"] - frame["reference_affinity"]

        valid_reference = frame["reference_affinity"].notna()
        if require_reference_pass:
            valid_reference &= redocking_class.eq("pass")
        frame["reference_anchor_available"] = valid_reference
        frame.loc[valid_reference, "reference_anchor_reason"] = "validated_reference"
        frame.loc[~valid_reference & frame["reference_affinity"].isna(), "reference_anchor_reason"] = "missing_reference_affinity"
        frame.loc[~valid_reference & frame["reference_affinity"].notna(), "reference_anchor_reason"] = "reference_validation_not_passed"

        delta_cutoff = float(reference_delta_kcal_mol)
        if delta_cutoff < 0:
            delta_cutoff = 1.0
        anchor_strong = valid_reference & (frame["delta_vs_reference"] <= 0.0)
        anchor_moderate = valid_reference & (frame["delta_vs_reference"] <= delta_cutoff)
        frame.loc[valid_reference, "docking_quality_class"] = "Weak"
        frame.loc[anchor_moderate, "docking_quality_class"] = "Moderate"
        frame.loc[anchor_strong, "docking_quality_class"] = "Strong"
        frame.loc[valid_reference, "classification_basis"] = "reference_anchor"
        frame.loc[~valid_reference, "classification_basis"] = "target_percentile_fallback"

    insufficient_n_mask = frame["target_ligand_count"] < 3
    insufficient_n_mask &= ~frame["reference_anchor_available"].astype(bool)
    frame.loc[insufficient_n_mask, "docking_quality_class"] = "Unclassified"
    frame.loc[insufficient_n_mask, "classification_basis"] = "insufficient_n"

    frame["qc_status"] = "pass"
    warn_low_engine = (~frame["single_engine"]) & (frame["agreement_count"] < 2)
    frame.loc[warn_low_engine, "qc_status"] = "warn_low_engine_agreement"
    frame.loc[frame["best_affinity_kcal_mol"] > 0, "qc_status"] = "fail_positive_affinity"

    # Keep ADMET interpretation independent from docking quality.
    frame["admet_status"] = "unknown"
    if "admet_flag_count" in frame.columns:
        flags = pd.to_numeric(frame["admet_flag_count"], errors="coerce").fillna(0).astype(int)
        frame.loc[flags == 0, "admet_status"] = "pass"
        frame.loc[flags > 0, "admet_status"] = "flagged"
    if "admet_status" in consensus_df.columns:
        frame["admet_status"] = frame["admet_status"].where(
            frame["admet_status"] != "unknown",
            consensus_df["admet_status"].astype(str),
        )

    # Atomic QC downgrade: evaluate once, then apply one final class.
    fail_mask = frame["qc_status"].astype(str).str.startswith("fail")
    warn_mask = frame["qc_status"].astype(str).str.startswith("warn")
    weak_mask = frame["docking_quality_class"] == "Weak"
    moderate_mask = frame["docking_quality_class"] == "Moderate"
    strong_mask = frame["docking_quality_class"] == "Strong"
    frame.loc[fail_mask, "docking_quality_class"] = "Weak"
    frame.loc[warn_mask & strong_mask, "docking_quality_class"] = "Moderate"
    frame.loc[warn_mask & moderate_mask, "docking_quality_class"] = "Weak"
    frame.loc[warn_mask & weak_mask, "docking_quality_class"] = "Weak"

    frame["hit_class_policy"] = normalized_policy
    frame["hit_class_strong_percentile"] = strong_pct
    frame["hit_class_moderate_percentile"] = moderate_pct
    return frame
