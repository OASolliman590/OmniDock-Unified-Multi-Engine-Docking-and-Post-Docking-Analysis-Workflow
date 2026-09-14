"""Consensus input contract tests that do not require chemistry backends."""

from pathlib import Path

import pandas as pd
import pytest

from post_docking_analysis.consensus import build_consensus_rankings


def _row(
    *,
    engine: str = "vina",
    tag: str = "P_site_L",
    protein: str = "P",
    ligand: str = "L",
    site_id: str = "site",
    score: float = -8.0,
    scoring_function: str | None = None,
) -> dict[str, object]:
    return {
        "engine": engine,
        "tag": tag,
        "protein": protein,
        "ligand": ligand,
        "site_id": site_id,
        "pose": 1,
        "pose_file": "missing.sdf",
        "affinity_kcal_mol": score,
        "scoring_function": scoring_function or engine,
        "receptor_frame_id": "frame",
    }


def test_duplicate_engine_tag_rows_are_rejected_before_aggregation() -> None:
    frame = pd.DataFrame([_row(score=-7.0), _row(score=-11.0)])

    with pytest.raises(ValueError, match="duplicate_engine_tag"):
        build_consensus_rankings(frame, expected_engines=["vina"])


def test_tag_identity_must_match_across_engines() -> None:
    frame = pd.DataFrame(
        [
            _row(engine="vina"),
            _row(engine="gnina", ligand="different_ligand"),
        ]
    )

    with pytest.raises(ValueError, match="inconsistent_tag_identity"):
        build_consensus_rankings(frame, expected_engines=["vina", "gnina"])


def test_mixed_scoring_functions_in_one_engine_target_site_are_rejected() -> None:
    frame = pd.DataFrame(
        [
            _row(tag="P_site_L1", ligand="L1", scoring_function="vina"),
            _row(tag="P_site_L2", ligand="L2", scoring_function="custom-vina"),
        ]
    )

    with pytest.raises(ValueError, match="mixed_scoring_functions"):
        build_consensus_rankings(frame, expected_engines=["vina"])


def test_expected_engine_scope_is_normalized_and_extras_are_rejected() -> None:
    frame = pd.DataFrame([_row(engine="VINA")])

    result = build_consensus_rankings(frame, expected_engines=[" vina ", "GNINA"])
    assert result.iloc[0].agreement_count == 1
    assert result.iloc[0].agreement_fraction == pytest.approx(0.5)

    with pytest.raises(ValueError, match="engine_outside_expected_scope"):
        build_consensus_rankings(frame, expected_engines=["gnina"])


def test_duplicate_expected_engines_are_rejected() -> None:
    with pytest.raises(ValueError, match="duplicate_expected_engines"):
        build_consensus_rankings(pd.DataFrame([_row()]), expected_engines=["vina", "VINA"])


@pytest.mark.parametrize("invalid", [None, float("nan"), ""])
def test_invalid_expected_engine_scope_is_rejected(invalid: object) -> None:
    with pytest.raises(ValueError, match="invalid_expected_engine_scope|empty_expected_engine_scope"):
        build_consensus_rankings(pd.DataFrame([_row()]), expected_engines=["vina", invalid])


def test_unique_rows_preserve_score_direction_and_missing_engine_coverage() -> None:
    frame = pd.DataFrame(
        [
            _row(engine="vina", tag="P_site_best", ligand="best", score=-10.0),
            _row(engine="vina", tag="P_site_worst", ligand="worst", score=-6.0),
        ]
    )

    result = build_consensus_rankings(
        frame,
        expected_engines=["vina", "gnina"],
        consensus_mode="weighted_hybrid",
    )
    assert result.iloc[0].tag == "P_site_best"
    assert result["agreement_fraction"].eq(0.5).all()
    assert result["agreement_count"].eq(1).all()
    assert (result["agreement_fraction"] <= 1.0).all()
