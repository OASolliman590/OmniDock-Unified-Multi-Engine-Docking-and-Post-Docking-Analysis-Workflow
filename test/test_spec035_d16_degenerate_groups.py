"""Spec 035 D16 follow-up: a min-max or z-score group with one finite value, or zero span, gives NaN.

Not 0.0 (the branch's old value) and not 0.5 (main's midpoint). A missing score also stays NaN, and
groups with spread are unchanged.

per_engine_rank is deliberately not covered here. Its single-value and tied-group outputs (1.0 and 0.5)
are not named by the Lead's decision, so they are left as they are and listed as an open question in
specs/035-main-reconciliation/port_log.md (section 5).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from post_docking_analysis.consensus import _safe_minmax_scale, _safe_zscore, normalize_engine_scores

DEGENERATE_METHODS = ("per_engine_minmax", "per_engine_zscore")


def _frame(rows):
    return pd.DataFrame(
        [{"engine": "vina", "protein": protein, "tag": tag, "affinity_kcal_mol": value} for protein, tag, value in rows]
    )


def _normalized(rows, method):
    out = normalize_engine_scores(_frame(rows), method=method, score_column="affinity_kcal_mol",
                                  normalized_column="normalized_affinity_score", group_keys=["engine", "protein"])
    return out.set_index("tag")["normalized_affinity_score"]


@pytest.mark.parametrize("method", DEGENERATE_METHODS)
def test_group_with_one_finite_value_is_nan(method):
    scores = _normalized([("P1", "alone", -7.0), ("P2", "b", -9.0), ("P2", "c", -6.0), ("P2", "d", -8.0)], method)
    assert np.isnan(scores["alone"])
    # The control group, which has spread, still gets finite scores in [0, 1].
    assert scores[["b", "c", "d"]].notna().all()
    assert ((scores[["b", "c", "d"]] >= 0.0) & (scores[["b", "c", "d"]] <= 1.0)).all()


@pytest.mark.parametrize("method", DEGENERATE_METHODS)
def test_group_with_zero_span_is_nan_not_zero_or_midpoint(method):
    scores = _normalized([("P1", "a", -8.0), ("P1", "b", -8.0), ("P1", "c", -8.0),
                          ("P2", "d", -9.0), ("P2", "e", -5.0)], method)
    for tag in ("a", "b", "c"):
        assert np.isnan(scores[tag]), f"{method}: zero-span row {tag} must be NaN, got {scores[tag]}"
    assert scores[["d", "e"]].notna().all()


@pytest.mark.parametrize("method", DEGENERATE_METHODS)
def test_one_finite_value_beside_a_missing_value_is_nan(method):
    scores = _normalized([("P1", "a", -8.0), ("P1", "b", np.nan)], method)
    assert np.isnan(scores["a"])
    assert np.isnan(scores["b"])


@pytest.mark.parametrize("method", DEGENERATE_METHODS)
def test_nonfinite_values_do_not_count_toward_group_spread(method):
    # Only two finite values remain, so the group has spread and the finite rows are scored; the
    # infinite and missing rows are NaN.
    scores = _normalized([("P1", "a", -8.0), ("P1", "b", -6.0), ("P1", "c", np.inf), ("P1", "d", np.nan)], method)
    assert scores[["a", "b"]].notna().all()
    assert np.isnan(scores["c"]) and np.isnan(scores["d"])


def test_minmax_helper_is_nan_for_single_value_and_zero_span():
    assert _safe_minmax_scale(pd.Series([-7.0])).isna().all()
    assert _safe_minmax_scale(pd.Series([-7.0, -7.0, -7.0])).isna().all()
    assert _safe_minmax_scale(pd.Series([-8.0, -4.0])).notna().all()


def test_zscore_helper_is_nan_for_single_value_and_zero_spread():
    assert _safe_zscore(pd.Series([-7.0])).isna().all()
    assert _safe_zscore(pd.Series([-7.0, -7.0, -7.0])).isna().all()
    assert _safe_zscore(pd.Series([-8.0, -4.0, -6.0])).notna().all()
