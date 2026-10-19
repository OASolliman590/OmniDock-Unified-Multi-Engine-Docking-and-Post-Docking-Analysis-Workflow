import numpy as np
import pandas as pd

from post_docking_analysis.consensus import normalize_engine_scores


def _scores(method):
    frame = pd.DataFrame({"engine": ["vina"] * 3, "protein": ["P"] * 3, "affinity_kcal_mol": [-9.0, -7.0, -4.0]})
    return normalize_engine_scores(frame, method=method)["normalized_affinity_score"].to_numpy()


def test_all_normalisers_give_the_strongest_binder_the_highest_score():
    for method in ("per_engine_rank", "per_engine_minmax", "per_engine_zscore"):
        scores = _scores(method)
        assert np.argmax(scores) == 0, method
        assert np.argmin(scores) == 2, method


def test_zscore_matches_minmax_order():
    np.testing.assert_allclose(_scores("per_engine_zscore"), _scores("per_engine_minmax"))
