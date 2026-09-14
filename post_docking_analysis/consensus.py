"""
Consensus ranking utilities for multi-engine post-docking analysis.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from post_docking_analysis.geometric_consensus import compute_geometric_consensus
from post_docking_analysis.score_semantics import score_spec

logger = logging.getLogger(__name__)


SUPPORTED_CONSENSUS_MODES = {
    "dockbox_geometric",
    "weighted_hybrid",
    "strict_consensus",
    "favorite_guardrails",
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

GEOMETRIC_COLUMNS = [
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


def normalize_consensus_mode(value: Optional[str]) -> str:
    token = str(value or "").strip().lower()
    return token if token in SUPPORTED_CONSENSUS_MODES else "dockbox_geometric"


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
    valid = values[finite_mask]
    out = pd.Series(np.nan, index=series.index, dtype=float)
    if len(valid) <= 1:
        out.loc[valid.index] = 0.5
        return out
    span = float(valid.max() - valid.min())
    if span <= 0:
        out.loc[valid.index] = 0.5
        return out
    # Docking convention: lower affinity is better, so invert min-max.
    scaled = (float(valid.max()) - values) / span
    return scaled.where(finite_mask)


def _safe_zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    finite_mask = np.isfinite(values.to_numpy(dtype=float))
    valid = values[finite_mask]
    out = pd.Series(np.nan, index=series.index, dtype=float)
    if len(valid) <= 1:
        out.loc[valid.index] = 0.0
        return out
    mean_value = float(valid.mean())
    std_value = float(valid.std(ddof=0))
    if std_value <= 0:
        out.loc[valid.index] = 0.0
        return out
    zscores = (values - mean_value) / std_value
    return zscores.where(finite_mask)


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
        # All normalizers operate on a lower-is-better internal score.
        series = series if score_spec(score_column).lower_is_better else -series
        if normalized_method == "per_engine_minmax":
            return _safe_minmax_scale(series)
        if normalized_method == "per_engine_zscore":
            zscores = _safe_zscore(series)
            # Docking z-score remains lower-is-better before inverse min-max.
            return _safe_minmax_scale(zscores)
        # Default: rank-based normalization with 1.0 = strongest.
        return _safe_rank_score(series, lower_is_better=True)

    working[normalized_column] = (
        working.groupby(keys, dropna=False)[score_column].transform(_normalize_group)
    )
    working[normalized_column] = pd.to_numeric(working[normalized_column], errors="coerce")
    working["normalization_method"] = normalized_method
    return working


def _validate_consensus_input(
    frame: pd.DataFrame,
    *,
    score_column: str,
    expected_engines: Optional[List[str]],
) -> tuple[pd.DataFrame, List[str]]:
    """Validate the one-row-per-engine/tag consensus boundary.

    Consensus receives one already-selected pose per engine and biological tag.
    Duplicate rows are ambiguous because they can represent seeds, retries,
    modes, or incompatible score records; silently selecting or pivoting one
    would create an undeclared scientific aggregation. The caller may still
    provide an expected engine that has no row, which remains a missing engine
    in the denominator.
    """
    required = {"engine", "tag", "protein", "ligand", "site_id", score_column}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"missing_consensus_columns: {','.join(missing)}")

    normalized = frame.copy()
    for column in ("engine", "tag", "protein", "ligand", "site_id"):
        values = []
        for value in normalized[column].tolist():
            if pd.isna(value):
                raise ValueError(f"missing_consensus_identity: {column}")
            token = str(value).strip()
            if not token or token.lower() in {"nan", "none", "null"}:
                raise ValueError(f"missing_consensus_identity: {column}")
            values.append(token.lower() if column == "engine" else token)
        normalized[column] = values

    if "scoring_function" not in normalized.columns:
        normalized["scoring_function"] = normalized["engine"]
    else:
        scoring_values = []
        for engine, value in zip(normalized["engine"], normalized["scoring_function"]):
            if pd.isna(value) or not str(value).strip():
                scoring_values.append(engine)
            else:
                scoring_values.append(str(value).strip().lower())
        normalized["scoring_function"] = scoring_values

    duplicate_mask = normalized.duplicated(subset=["engine", "tag"], keep=False)
    if duplicate_mask.any():
        duplicate_keys = sorted(
            {
                (str(row.engine), str(row.tag))
                for row in normalized.loc[duplicate_mask, ["engine", "tag"]].itertuples(index=False)
            }
        )
        raise ValueError(f"duplicate_engine_tag_rows: {duplicate_keys}")

    identity_counts = normalized.groupby("tag", dropna=False)[
        ["protein", "ligand", "site_id"]
    ].nunique(dropna=False)
    if (identity_counts > 1).any().any():
        raise ValueError("inconsistent_tag_identity")

    scoring_counts = normalized.groupby(
        ["engine", "protein", "site_id"], dropna=False
    )["scoring_function"].nunique(dropna=False)
    if (scoring_counts > 1).any():
        raise ValueError("mixed_scoring_functions")

    if expected_engines is None:
        engines = sorted(normalized["engine"].unique().tolist())
    else:
        raw_engines = [expected_engines] if isinstance(expected_engines, str) else list(expected_engines)
        # Existing pipeline callers use [] as the legacy sentinel meaning
        # "infer the observed engine scope". Keep that compatibility while
        # validating every declared, non-empty scope below.
        if not raw_engines:
            engines = sorted(normalized["engine"].unique().tolist())
            return normalized, engines
        engines = []
        for value in raw_engines:
            if pd.isna(value):
                raise ValueError("invalid_expected_engine_scope")
            token = str(value).strip().lower()
            if not token or token in {"nan", "none", "null"}:
                raise ValueError("invalid_expected_engine_scope")
            engines.append(token)
        if not engines:
            raise ValueError("empty_expected_engine_scope")
        if len(set(engines)) != len(engines):
            raise ValueError("duplicate_expected_engines")
        extras = sorted(set(normalized["engine"]) - set(engines))
        if extras:
            raise ValueError(f"engine_outside_expected_scope: {extras}")
    if not engines:
        raise ValueError("empty_consensus_engine_scope")
    return normalized, engines


def build_consensus_rankings(
    best_by_engine: pd.DataFrame,
    consensus_mode: str = "dockbox_geometric",
    favorite_engine: Optional[str] = None,
    normalization_method: str = "per_engine_rank",
    score_column: str = "affinity_kcal_mol",
    normalized_column: str = "normalized_affinity_score",
    expected_engines: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Build DockBox-style consensus ranking table from one best pose per (engine, tag).

    Note:
    `strict_consensus` requires multi-engine agreement. In single-engine mode the
    strict agreement filter cannot be applied; behavior matches weighted/hybrid
    scoring and a warning is emitted.
    """
    mode = normalize_consensus_mode(consensus_mode)
    favorite = str(favorite_engine or "").strip().lower()

    if best_by_engine is None or best_by_engine.empty:
        return pd.DataFrame(
            columns=[
                "tag",
                "protein",
                "ligand",
                "site_id",
                "score_name_primary",
                "score_primary",
                "score_name_secondary",
                "score_secondary",
                "agreement_count",
                "agreement_fraction",
                "winner_engine",
                "best_affinity_kcal_mol",
                "mean_affinity_kcal_mol",
                "affinity_spread_kcal_mol",
                "mean_rank_pct",
                "best_affinity_rank_pct",
                "favorite_rank_pct",
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
                "consensus_score",
                "consensus_mode",
                "consensus_rank_within_protein",
                "consensus_rank_global",
                "engine_votes",
            ]
        )

    frame, engines = _validate_consensus_input(
        best_by_engine,
        score_column=score_column,
        expected_engines=expected_engines,
    )
    frame[score_column] = pd.to_numeric(frame[score_column], errors="coerce")
    frame = frame[np.isfinite(frame[score_column])].copy()
    if frame.empty:
        return pd.DataFrame()

    frame = normalize_engine_scores(
        frame,
        method=normalization_method,
        score_column=score_column,
        normalized_column=normalized_column,
        group_keys=["engine", "protein", "site_id", "scoring_function"],
    )
    total_engines = max(len(engines), 1)

    frame["engine_rank_pct"] = (
        frame.groupby(["engine", "protein", "site_id", "scoring_function"], dropna=False)[normalized_column]
        .transform(lambda s: _safe_rank_score(s, lower_is_better=False))
    )

    affinity_matrix = frame.pivot(index="tag", columns="engine", values=score_column)
    agreement_count = affinity_matrix.notna().sum(axis=1).astype(int)
    single_engine = total_engines <= 1
    if single_engine:
        agreement_fraction = pd.Series(np.nan, index=agreement_count.index, dtype=float)
    else:
        agreement_fraction = agreement_count / float(total_engines)
    normalized_matrix = frame.pivot(index="tag", columns="engine", values=normalized_column)
    winner_engine = normalized_matrix.idxmax(axis=1, skipna=True).fillna("").astype(str)
    spread = normalized_matrix.max(axis=1, skipna=True) - normalized_matrix.min(axis=1, skipna=True)

    by_tag = (
        frame.sort_values(["tag", "affinity_kcal_mol"])
        .groupby("tag", as_index=False)
        .first()[["tag", "protein", "ligand", "site_id"]]
    )
    ranking_columns = ["score_name_primary", "score_primary", "score_name_secondary", "score_secondary"]
    available_ranking_columns = [column for column in ranking_columns if column in frame.columns]
    if available_ranking_columns:
        primary_rows = (
            frame.sort_values(["tag", score_column, "engine"], ascending=[True, True, True])
            .groupby("tag", as_index=False)
            .first()[["tag"] + available_ranking_columns]
        )
        by_tag = by_tag.merge(primary_rows, on="tag", how="left")
    for column in ranking_columns:
        if column not in by_tag.columns:
            by_tag[column] = np.nan if "score" in column and not column.startswith("score_name") else ""
    by_tag = by_tag.merge(
        frame.groupby("tag", as_index=False).agg(
            best_affinity_kcal_mol=("affinity_kcal_mol", "min"),
            mean_affinity_kcal_mol=("affinity_kcal_mol", "mean"),
            mean_normalized_affinity_score=(normalized_column, "mean"),
            mean_rank_pct=("engine_rank_pct", "mean"),
        ),
        on="tag",
        how="left",
    )
    by_tag["agreement_count"] = by_tag["tag"].map(agreement_count).fillna(0).astype(int)
    by_tag["agreement_fraction"] = by_tag["tag"].map(agreement_fraction)
    by_tag["single_engine"] = bool(single_engine)
    by_tag["winner_engine"] = by_tag["tag"].map(winner_engine).fillna("").astype(str)
    by_tag["affinity_spread_kcal_mol"] = by_tag["tag"].map(spread).fillna(0.0).astype(float)
    by_tag["engine_votes"] = by_tag["tag"].map(
        affinity_matrix.apply(lambda row: ",".join([col for col in affinity_matrix.columns if pd.notna(row[col])]), axis=1)
    ).fillna("")

    by_tag["best_affinity_rank_pct"] = (
        by_tag.groupby("protein", dropna=False)["best_affinity_kcal_mol"]
        .transform(lambda s: _safe_rank_score(s, lower_is_better=True))
    )
    # Dimensionless disagreement: larger spread must never improve a score.
    by_tag["spread_norm"] = by_tag["tag"].map(spread).fillna(0.0)
    selected = frame.sort_values(["tag", normalized_column, "engine"], ascending=[True, False, True]).drop_duplicates("tag")
    selected = selected.set_index("tag")
    by_tag["best_affinity_kcal_mol"] = by_tag["tag"].map(selected["affinity_kcal_mol"])
    by_tag["winner_scoring_function"] = by_tag["tag"].map(selected["scoring_function"])
    for column in available_ranking_columns:
        by_tag[column] = by_tag["tag"].map(selected[column])
    by_tag["best_affinity_rank_pct"] = np.nan
    # Raw means and ranges across different scoring functions have no shared scale.
    by_tag.loc[by_tag["agreement_count"] > 1, ["mean_affinity_kcal_mol", "affinity_spread_kcal_mol"]] = np.nan

    by_tag["favorite_rank_pct"] = np.nan
    if favorite:
        favorite_rows = frame[frame["engine"] == favorite].copy()
        if not favorite_rows.empty:
            favorite_rows["favorite_rank_pct"] = (
                favorite_rows.groupby("protein", dropna=False)[normalized_column]
                .transform(lambda s: _safe_rank_score(s, lower_is_better=False))
            )
            by_tag = by_tag.merge(
                favorite_rows[["tag", "favorite_rank_pct"]].drop_duplicates("tag"),
                on="tag",
                how="left",
                suffixes=("", "_fav"),
            )
            by_tag["favorite_rank_pct"] = pd.to_numeric(
                by_tag.get("favorite_rank_pct_fav", by_tag["favorite_rank_pct"]),
                errors="coerce",
            )
            by_tag.drop(columns=["favorite_rank_pct_fav"], inplace=True, errors="ignore")

    if mode == "strict_consensus" and not single_engine:
        favorable = normalized_matrix.ge(0.5).sum(axis=1)
        by_tag = by_tag[(by_tag["agreement_count"] == total_engines) & by_tag["tag"].map(favorable).ge(2)].copy()
    elif mode == "strict_consensus" and single_engine:
        logger.warning(
            "strict_consensus requested in single-engine mode; agreement filtering is skipped."
        )

    if by_tag.empty:
        return by_tag

    for column in GEOMETRIC_COLUMNS:
        if column not in by_tag.columns:
            by_tag[column] = np.nan
    by_tag["geometric_agreement"] = False
    by_tag["geometric_reason"] = "not_computed"
    by_tag["geometric_pairwise_rmsd_json"] = "[]"

    if mode == "dockbox_geometric":
        geometric_df = compute_geometric_consensus(
            frame.copy(),
            rmsd_cutoff=2.0,
            expected_engines=engines,
        )
        if geometric_df is not None and not geometric_df.empty:
            by_tag = by_tag.merge(geometric_df, on="tag", how="left", suffixes=("", "_geom"))
            for column in GEOMETRIC_COLUMNS:
                geom_col = f"{column}_geom"
                if geom_col in by_tag.columns:
                    by_tag[column] = by_tag[geom_col]
                    by_tag.drop(columns=[geom_col], inplace=True, errors="ignore")
            by_tag["geometric_agreement"] = by_tag["geometric_agreement"].fillna(False).astype(bool)
            by_tag["geometric_reason"] = by_tag["geometric_reason"].fillna("not_available")
            by_tag["geometric_pairwise_rmsd_json"] = by_tag["geometric_pairwise_rmsd_json"].fillna("[]")

    if mode == "dockbox_geometric":
        agreement_component = by_tag["geometric_agreement"].astype(float)
        by_tag["consensus_score"] = (
            agreement_component
            + by_tag["mean_rank_pct"].fillna(0.0) * 1e-3
            + (1.0 - by_tag["spread_norm"].fillna(0.0)) * 1e-4
        )
    elif mode == "favorite_guardrails":
        base_favorite = pd.to_numeric(by_tag["favorite_rank_pct"], errors="coerce")
        base_favorite = base_favorite.fillna(by_tag["mean_rank_pct"])
        guardrail_penalty = (
            (
                (~by_tag["single_engine"].astype(bool))
                & (by_tag["agreement_count"] < 2)
            ).astype(float) * 0.25
            + by_tag["spread_norm"].fillna(0.0) * 0.10
            + (by_tag["winner_engine"] != favorite).astype(float) * (0.10 if favorite else 0.0)
        )
        by_tag["consensus_score"] = base_favorite - guardrail_penalty
    elif mode == "strict_consensus":
        by_tag["consensus_score"] = (
            by_tag["mean_rank_pct"].fillna(0.0) * 0.90
            + (1.0 - by_tag["spread_norm"].fillna(0.0)) * 0.10
        )
    else:
        agreement_component = pd.to_numeric(by_tag["agreement_fraction"], errors="coerce").fillna(1.0)
        by_tag["consensus_score"] = (
            by_tag["mean_rank_pct"].fillna(0.0) * 0.85
            + agreement_component * 0.10
            + (1.0 - by_tag["spread_norm"].fillna(0.0)) * 0.05
        )

    by_tag["consensus_mode"] = mode
    by_tag["normalization_method"] = normalize_normalization_method(normalization_method)
    by_tag.sort_values(
        ["protein", "consensus_score", "tag"],
        ascending=[True, False, True],
        inplace=True,
    )
    by_tag["consensus_rank_within_protein"] = by_tag.groupby("protein", dropna=False).cumcount() + 1
    global_order = by_tag.sort_values(["consensus_score", "tag"], ascending=[False, True]).index
    by_tag["consensus_rank_global"] = pd.Series(np.arange(1, len(by_tag) + 1), index=global_order)

    ordered_columns = [
        "tag",
        "protein",
        "ligand",
        "site_id",
        "score_name_primary",
        "score_primary",
        "score_name_secondary",
        "score_secondary",
        "agreement_count",
        "agreement_fraction",
        "single_engine",
        "winner_engine",
        "winner_scoring_function",
        "best_affinity_kcal_mol",
        "mean_affinity_kcal_mol",
        "mean_normalized_affinity_score",
        "affinity_spread_kcal_mol",
        "mean_rank_pct",
        "best_affinity_rank_pct",
        "favorite_rank_pct",
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
        "consensus_score",
        "consensus_mode",
        "normalization_method",
        "consensus_rank_within_protein",
        "consensus_rank_global",
        "engine_votes",
    ]
    return by_tag[ordered_columns].copy()


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
        ["consensus_score", "tag"],
        ascending=[False, True],
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
                ["consensus_score", "tag"],
                ascending=[False, True],
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
        for required in ("protein", "engine", "scoring_function", "reference_affinity", "reference_tag", "redocking_classification"):
            if required not in baseline_df.columns:
                if required == "reference_affinity":
                    baseline_df[required] = np.nan
                else:
                    baseline_df[required] = ""
        baseline_df = baseline_df[["protein", "engine", "scoring_function", "reference_affinity", "reference_tag", "redocking_classification"]].copy()
        baseline_df["reference_affinity"] = pd.to_numeric(baseline_df["reference_affinity"], errors="coerce")
        baseline_df["redocking_classification"] = baseline_df["redocking_classification"].fillna("").astype(str).str.lower()
        keys = ["protein", "engine", "scoring_function"]
        if baseline_df.duplicated(keys).any():
            raise ValueError("Reference baselines must be unique per target, engine and scoring function")
        for column in ("winner_engine", "winner_scoring_function"):
            if column not in frame:
                frame[column] = ""

        frame = frame.merge(
            baseline_df,
            left_on=["protein", "winner_engine", "winner_scoring_function"],
            right_on=keys,
            how="left",
            suffixes=("", "_baseline"),
            validate="many_to_one",
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

    frame["qc_status"] = "not_evaluated_physical_validity"
    frame["physical_validity_status"] = "not_evaluated"
    frame["interpretation"] = "relative_docking_prioritization_not_binding_evidence"
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
