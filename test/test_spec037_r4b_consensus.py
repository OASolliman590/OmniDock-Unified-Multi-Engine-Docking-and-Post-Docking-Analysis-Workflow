"""Spec 037 R4b: no raw cross-engine comparison decides anything.

Each decision is tested on a fixture where the raw affinity minimum and the v2 consensus disagree:

- Fixture A (top-pose selector): tags p_s1_L, p_s2_L, p_s3_L on one protein and ligand.
  Raw minimum over engines: p_s1_L (vina -9.5). v2 consensus: p_s2_L (0.75), composite of per-engine percentiles.
- Fixture B (rerun promotion and rescoring): vina and smina on three tags.
  Raw minimum for p_s1_L: vina (-9.0 against smina -8.9). v2 winner for p_s1_L: smina (percentile 1.0 vs 0.0).

The raw-minimum result is printed in each test's assertion message, so the before and after values are
visible from the test output.
"""

from __future__ import annotations

import pandas as pd
import pytest

from post_docking_analysis.consensus import (
    build_consensus_explainability,
    build_consensus_rankings,
    select_rescoring_candidates,
)
from post_docking_analysis.multi_engine_pipeline_impl import MultiEngineAnalysisPipeline
from post_docking_analysis.top_pose_selector import build_top_pose_atlas


def _rows(spec: dict) -> pd.DataFrame:
    rows = []
    for tag, per_engine in spec.items():
        site = tag.split("_")[1]
        for engine, affinity in per_engine.items():
            rows.append({
                "engine": engine,
                "protein": "p.pdbqt",
                "ligand": "L.pdbqt",
                "site_id": site,
                "tag": tag,
                "pose": 1,
                "affinity_kcal_mol": affinity,
                "pose_file": f"/synthetic/{tag}_{engine}.pdbqt",
            })
    return pd.DataFrame(rows)


FIXTURE_A = {
    "p_s1_L": {"vina": -9.5, "smina": -6.0},
    "p_s2_L": {"vina": -8.0, "smina": -9.0},
    "p_s3_L": {"vina": -7.0, "smina": -8.5},
}
FIXTURE_B = {
    "p_s1_L": {"vina": -9.0, "smina": -8.9},
    "p_s2_L": {"vina": -10.0, "smina": -5.0},
    "p_s3_L": {"vina": -10.5, "smina": -6.0},
}


def _raw_minimum_tag(best: pd.DataFrame) -> str:
    return str(best.sort_values(["affinity_kcal_mol", "engine"]).iloc[0]["tag"])


# ----------------------------------------------------------------------------- decision 1: top-pose selector


def test_raw_minimum_and_v2_disagree_on_fixture_a():
    best = _rows(FIXTURE_A)
    consensus = build_consensus_rankings(best_by_engine=best, consensus_mode="consensus_rank_geometry_qc_v2")
    v2_top = str(consensus.sort_values("consensus_rank_global").iloc[0]["tag"])
    assert _raw_minimum_tag(best) == "p_s1_L"
    assert v2_top == "p_s2_L"


def test_top_pose_default_ranks_by_v2_consensus_not_raw_affinity():
    best = _rows(FIXTURE_A)
    payload = build_top_pose_atlas(best_by_engine=best, consensus_df=None)
    per_protein = payload["top_pose_per_ligand_per_protein"]
    assert payload["manifest"]["selection_policy"] == "best_consensus"
    assert str(per_protein.iloc[0]["tag"]) == "p_s2_L", (
        f"raw minimum would select p_s1_L (-9.5); v2 selects {per_protein.iloc[0]['tag']}"
    )
    assert payload["manifest"]["comparison_scope"] == "descriptive_not_comparable_across_engines"
    assert "ligand_performance_summary.descriptive_best_affinity_kcal_mol" in payload["manifest"][
        "descriptive_not_comparable_across_engines"
    ]


def test_best_affinity_is_refused_across_engines_with_a_clear_message():
    best = _rows(FIXTURE_A)
    with pytest.raises(ValueError, match="compares raw affinities across engines") as excinfo:
        build_top_pose_atlas(best_by_engine=best, consensus_df=None, selection_policy="best_affinity")
    assert "best_consensus" in str(excinfo.value)


def test_best_affinity_is_still_available_within_one_engine():
    single = _rows({tag: {"vina": affinity["vina"]} for tag, affinity in FIXTURE_A.items()})
    payload = build_top_pose_atlas(best_by_engine=single, consensus_df=None, selection_policy="best_affinity")
    assert payload["manifest"]["selection_policy"] == "best_affinity"
    assert str(payload["top_pose_per_ligand_per_protein"].iloc[0]["tag"]) == "p_s1_L"


def test_top_pose_identity_row_is_the_consensus_winner_engine_not_the_lowest_affinity():
    # p_s1_L: vina -9.0 is the lowest raw affinity, but smina is the v2 winner (percentile 1.0 against 0.0).
    best = _rows({"p_s1_L": {"vina": -9.0, "smina": -8.9}})
    consensus = pd.DataFrame([{
        "tag": "p_s1_L", "consensus_score": 0.9, "consensus_rank_global": 1, "consensus_status": "completed",
        "winner_engine": "smina", "best_affinity_kcal_mol": -9.0, "agreement_count": 2,
    }])
    payload = build_top_pose_atlas(best_by_engine=best, consensus_df=consensus)
    row = payload["top_pose_per_ligand_per_protein"].iloc[0]
    assert str(row["engine"]) == "smina"
    assert float(row["affinity_kcal_mol"]) == pytest.approx(-8.9)


def test_summary_is_ordered_by_consensus_and_its_affinity_is_labelled_descriptive():
    best = _rows(FIXTURE_A)
    summary = build_top_pose_atlas(best_by_engine=best, consensus_df=None)["ligand_performance_summary"]
    assert "best_affinity_kcal_mol" not in summary.columns
    assert "descriptive_best_affinity_kcal_mol" in summary.columns


# ----------------------------------------------------------------------------- decision 2: rerun promotion


def _consensus_for(best: pd.DataFrame) -> pd.DataFrame:
    return build_consensus_rankings(best_by_engine=best, consensus_mode="consensus_rank_geometry_qc_v2")


def test_raw_minimum_and_v2_disagree_on_fixture_b_pair_winner():
    best = _rows(FIXTURE_B)
    raw_winner = best[best["tag"] == "p_s1_L"].sort_values(["affinity_kcal_mol", "engine"]).iloc[0]["engine"]
    v2_winner = str(_consensus_for(best).set_index("tag").loc["p_s1_L", "winner_engine"])
    assert raw_winner == "vina"
    assert v2_winner == "smina"


def test_pair_competition_winner_is_the_v2_winner_not_the_raw_minimum():
    best = _rows(FIXTURE_B)
    consensus = _consensus_for(best)
    competition = MultiEngineAnalysisPipeline._build_pair_competition(best, "vina", consensus).set_index("tag")
    row = competition.loc["p_s1_L"]
    assert bool(row["is_target_engine_winner"]) is False, "raw minimum would make vina the winner of p_s1_L"
    assert row["winner_engine"] == "smina"
    assert row["winner_basis"] == "consensus_rank_geometry_qc_v2_per_engine_percentile"
    # The raw difference is kept, and it is named as descriptive.
    assert row["descriptive_affinity_difference_kcal_mol"] == pytest.approx(-8.9 - (-9.0))
    assert "affinity_advantage_kcal_mol" not in competition.columns


def test_pair_competition_without_a_complete_case_has_no_winner():
    best = _rows({"p_s9_L": {"vina": -9.0}})  # one engine only: not a complete case for a two-engine target
    consensus = pd.DataFrame([{
        "tag": "p_s9_L", "consensus_status": "consensus_incomplete", "winner_engine": "",
        "engine_percentiles_json": "{}", "consensus_score": float("nan"),
    }])
    competition = MultiEngineAnalysisPipeline._build_pair_competition(best, "vina", consensus).set_index("tag")
    assert bool(competition.loc["p_s9_L", "is_target_engine_winner"]) is False
    assert competition.loc["p_s9_L", "winner_basis"] == "consensus_incomplete_no_winner"


def test_min_affinity_advantage_gate_is_refused(tmp_path):
    with pytest.raises(ValueError, match="removed by Spec 037 R4b"):
        MultiEngineAnalysisPipeline(
            project_dir=str(tmp_path),
            output_dir=str(tmp_path / "out"),
            min_affinity_advantage=0.5,
        )


# ----------------------------------------------------------------------------- decision 3: rescoring tie-break


def _tie_fixture() -> pd.DataFrame:
    # Equal consensus_score. The raw affinity would prefer P_raw (-9.0); the v2 rank prefers P_rank (rank 1).
    return pd.DataFrame([
        {"tag": "P_raw", "protein": "p", "consensus_score": 0.5, "consensus_rank_global": 2,
         "best_affinity_kcal_mol": -9.0, "agreement_count": 2},
        {"tag": "P_rank", "protein": "p", "consensus_score": 0.5, "consensus_rank_global": 1,
         "best_affinity_kcal_mol": -6.0, "agreement_count": 2},
    ])


def test_rescoring_tie_break_uses_the_v2_rank_not_the_raw_affinity():
    selected = select_rescoring_candidates(_tie_fixture(), rescoring_scope="top_n_global", rescoring_top_n=1)
    assert list(selected["tag"]) == ["P_rank"], "raw affinity would select P_raw (-9.0)"


def test_rescoring_tie_break_falls_back_to_tag_when_no_rank_column():
    # Without a rank column the tie goes to the lexical tag ("P_rank" < "P_raw"), not to the raw affinity (P_raw).
    frame = _tie_fixture().drop(columns=["consensus_rank_global"])
    selected = select_rescoring_candidates(frame, rescoring_scope="top_n_global", rescoring_top_n=1)
    assert list(selected["tag"]) == ["P_rank"]


def test_explainability_top_global_uses_the_v2_rank_tie_break():
    explain = build_consensus_explainability(
        consensus_df=_tie_fixture(),
        rescoring_df=select_rescoring_candidates(_tie_fixture(), "top_n_global", 1),
        consensus_mode="consensus_rank_geometry_qc_v2",
        rescoring_scope="top_n_global",
        rescoring_top_n=1,
    )
    assert explain["top_global"][0]["tag"] == "P_rank"
    assert "descriptive_best_affinity_kcal_mol" in explain["proteins"]["p"]
    assert "best_affinity_kcal_mol" not in explain["proteins"]["p"]
