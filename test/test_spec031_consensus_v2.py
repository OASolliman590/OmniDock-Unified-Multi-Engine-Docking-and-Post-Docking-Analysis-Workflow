import json

import numpy as np
import pandas as pd

from post_docking_analysis.consensus import build_consensus_rankings, select_representative_poses_v2
from post_docking_analysis.top_pose_selector import build_top_pose_atlas


def _scores() -> pd.DataFrame:
    rows = []
    for tag, ligand, vina, smina, ad4, cnn_affinity, cnn_score in [
        ("p__winner", "winner", -9.0, -8.0, -7.0, 7.5, 0.80),
        ("p__other", "other", -7.0, -6.0, -5.0, 5.0, 0.70),
    ]:
        for engine, affinity in (("vina", vina), ("smina", smina), ("ad4", ad4)):
            rows.append({"engine": engine, "protein": "p", "tag": tag, "ligand": ligand,
                         "site_id": "site1", "pose": 1, "affinity_kcal_mol": affinity})
        rows.append({"engine": "gnina", "protein": "p", "tag": tag, "ligand": ligand,
                     "site_id": "site1", "pose": 1, "affinity_kcal_mol": vina,
                     "cnn_affinity": cnn_affinity, "cnn_score": cnn_score})
    return pd.DataFrame(rows)


def test_v2_respects_all_engine_score_directions() -> None:
    result = build_consensus_rankings(_scores(), consensus_mode="consensus_rank_geometry_qc_v2")
    assert result.iloc[0]["ligand"] == "winner"
    assert result.iloc[0]["consensus_score"] == 1.0
    assert result.iloc[1]["consensus_score"] == 0.0
    assert json.loads(result.iloc[0]["engine_weights_json"]) == {
        "ad4": 0.25, "gnina": 0.25, "smina": 0.25, "vina": 0.25
    }


def test_gnina_pose_selection_uses_cnn_score_then_ranks_cnn_affinity() -> None:
    rows = pd.DataFrame([
        {"engine": "gnina", "protein": "p", "tag": "x", "ligand": "x", "pose": 1,
         "affinity_kcal_mol": -10.0, "cnn_score": 0.2, "cnn_affinity": 9.0},
        {"engine": "gnina", "protein": "p", "tag": "x", "ligand": "x", "pose": 2,
         "affinity_kcal_mol": -7.0, "cnn_score": 0.9, "cnn_affinity": 6.0},
    ])
    selected = select_representative_poses_v2(rows)
    assert int(selected.iloc[0]["pose"]) == 2
    assert float(selected.iloc[0]["v2_ranking_metric"]) == 6.0


def test_missing_eligible_engine_score_is_incomplete_without_imputation() -> None:
    scores = _scores()
    scores = scores[~((scores["tag"] == "p__other") & (scores["engine"] == "smina"))]
    result = build_consensus_rankings(scores, consensus_mode="consensus_rank_geometry_qc_v2")
    incomplete = result[result["tag"] == "p__other"].iloc[0]
    assert incomplete["consensus_status"] == "consensus_incomplete"
    assert np.isnan(incomplete["consensus_score"])
    assert incomplete["included_engines"] == ""


def test_input_permutation_does_not_change_v2_output() -> None:
    scores = _scores()
    first = build_consensus_rankings(scores, consensus_mode="consensus_rank_geometry_qc_v2")
    second = build_consensus_rankings(
        scores.sample(frac=1.0, random_state=31), consensus_mode="consensus_rank_geometry_qc_v2"
    )
    columns = ["tag", "consensus_score", "consensus_rank_within_protein", "engine_percentiles_json"]
    pd.testing.assert_frame_equal(first[columns].reset_index(drop=True), second[columns].reset_index(drop=True))


def test_geometry_breaks_only_fully_comparable_exact_tie() -> None:
    scores = pd.DataFrame([
        {"engine": "vina", "protein": "p", "tag": "a", "ligand": "a", "pose": 1,
         "affinity_kcal_mol": -7.0, "geometry_status": "comparable", "geometric_agreement": False},
        {"engine": "vina", "protein": "p", "tag": "b", "ligand": "b", "pose": 1,
         "affinity_kcal_mol": -7.0, "geometry_status": "comparable", "geometric_agreement": True},
    ])
    result = build_consensus_rankings(scores, consensus_mode="consensus_rank_geometry_qc_v2")
    assert list(result["tag"]) == ["b", "a"]


def test_single_engine_is_disclosed() -> None:
    scores = _scores()
    scores = scores[scores["engine"] == "vina"]
    result = build_consensus_rankings(scores, consensus_mode="consensus_rank_geometry_qc_v2")
    assert set(result["consensus_status"]) == {"single_engine"}


def test_top_pose_selector_treats_v2_higher_score_as_better() -> None:
    best = pd.DataFrame([
        {"engine": "vina", "protein": "p", "ligand": "lig", "site_id": "a", "tag": "low",
         "pose": 1, "affinity_kcal_mol": -9.0},
        {"engine": "vina", "protein": "p", "ligand": "lig", "site_id": "b", "tag": "high",
         "pose": 1, "affinity_kcal_mol": -7.0},
    ])
    consensus = pd.DataFrame([
        {"tag": "low", "consensus_score": 0.1, "agreement_count": 2},
        {"tag": "high", "consensus_score": 0.9, "agreement_count": 2},
    ])
    atlas = build_top_pose_atlas(
        best, consensus, selection_policy="best_consensus",
        consensus_mode="consensus_rank_geometry_qc_v2",
    )
    assert atlas["top_pose_per_ligand_per_protein"].iloc[0]["tag"] == "high"
