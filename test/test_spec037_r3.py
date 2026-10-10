"""Spec 037 R3: basic-mode exhaustiveness, summaries on the consensus rank, and the GNINA summary label.

- R3a: basic mode refuses an explicit --exhaustiveness (as it refuses --num-modes). The effective value and its
  source are printed and recorded in run_manifest.json. Covered for the docking CLI and the workflow `dock run` path.
- R3b: best-pose-per-tag, top overall, best per protein and per ligand, and the protein/ligand summaries follow the
  consensus rank (consensus_rank_geometry_qc_v2 for two or more engines, single_engine_native_v1 for one). A raw
  affinity is never compared across engines to order a summary.
- R3c: the single-engine GNINA summary names cnn_score as the pose selection key.

Expected values are closed-form (percentiles of three tags are 0, 0.5 and 1), not read back from the code under test.
No docking is run: dry runs and synthetic score tables only.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import docking.cli as docking_cli
from docking.parameter_schema import resolve_parameter_schema, validate_parameter_schema
from post_docking_analysis.consensus import CONSENSUS_V2_METHOD, SINGLE_ENGINE_NATIVE_METHOD
from post_docking_analysis.multi_engine_pipeline_impl import (
    RANKING_UNAVAILABLE_CONSENSUS_INCOMPLETE,
    UNIFIED_COMPAT_COLUMNS,
    MultiEngineAnalysisPipeline,
    ranking_method_for_scores,
)
from post_docking_analysis.simplified_pipeline_impl import SimplifiedPostDockingPipeline
from workflow import execution
from workflow.cli import main as workflow_main

from test_spec036_docking_params import _dock_args, _dock_project, _run_manifest


# ----------------------------------------------------------------------------- fixtures


def _schema(mode: str, *, preset: str = "exhaustive", exhaustiveness: int = 32, explicit: bool = False):
    return resolve_parameter_schema(
        mode=mode,
        preset=preset,
        exhaustiveness=exhaustiveness,
        num_modes=20,
        seed=1,
        box_scale=1.0,
        box_padding=0.0,
        runtime_by_engine={"vina": {}},
        exhaustiveness_explicit=explicit,
    )


def _two_engine_disagreement() -> pd.DataFrame:
    """Three tags on one target, two engines. The raw minimum across engines and the v2 rank disagree.

    Raw minimum: ligand A, vina, -10.0 kcal/mol (the lowest affinity in the table).
    v2 per-engine percentiles over the three complete tags (vina: lower is better; gnina: cnn_affinity, higher is better):
      vina   A=-10.0 -> 1.0, B=-9.0 -> 0.5, C=-8.0 -> 0.0
      gnina  A=3.0   -> 0.0, B=6.0   -> 1.0, C=5.0  -> 0.5
    Composite (mean of percentiles): A=0.50, B=0.75, C=0.25, so the v2 order is B, A, C.
    Winner engines: A vina, B gnina, C gnina. GNINA pose choice is the highest cnn_score (one pose per tag here).
    """
    rows = [
        {"engine": "vina", "protein": "P1", "ligand": "A", "site_id": "site_1", "tag": "P1_site_1_A", "pose": 1,
         "affinity_kcal_mol": -10.0},
        {"engine": "vina", "protein": "P1", "ligand": "B", "site_id": "site_1", "tag": "P1_site_1_B", "pose": 1,
         "affinity_kcal_mol": -9.0},
        {"engine": "vina", "protein": "P1", "ligand": "C", "site_id": "site_1", "tag": "P1_site_1_C", "pose": 1,
         "affinity_kcal_mol": -8.0},
        {"engine": "gnina", "protein": "P1", "ligand": "A", "site_id": "site_1", "tag": "P1_site_1_A", "pose": 1,
         "affinity_kcal_mol": -6.0, "cnn_score": 0.5, "cnn_affinity": 3.0},
        {"engine": "gnina", "protein": "P1", "ligand": "B", "site_id": "site_1", "tag": "P1_site_1_B", "pose": 1,
         "affinity_kcal_mol": -7.0, "cnn_score": 0.8, "cnn_affinity": 6.0},
        {"engine": "gnina", "protein": "P1", "ligand": "C", "site_id": "site_1", "tag": "P1_site_1_C", "pose": 1,
         "affinity_kcal_mol": -6.5, "cnn_score": 0.6, "cnn_affinity": 5.0},
    ]
    return pd.DataFrame(rows)


def _with_unified_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """The writers need the full unified score schema. Missing columns are blank; ranking does not read them."""
    frame = frame.copy()
    for column in UNIFIED_COMPAT_COLUMNS:
        if column not in frame.columns:
            frame[column] = None
    if "pose_file" in frame.columns:
        frame["pose_file"] = frame["pose_file"].fillna("")
    return frame


def _single_gnina_selection_fixture() -> pd.DataFrame:
    """One engine. X's pose 1 has the highest cnn_score (picked) but not the highest cnn_affinity.

    Ranking metric cnn_affinity (higher is better) for the picked poses: Y 9.0, Z 6.0, X 4.0, so the order is Y, Z, X.
    """
    rows = [
        {"engine": "gnina", "protein": "P1", "ligand": "X", "site_id": "site_1", "tag": "P1_site_1_X", "pose": 1,
         "affinity_kcal_mol": -7.0, "cnn_score": 0.9, "cnn_affinity": 4.0, "pose_file": ""},
        {"engine": "gnina", "protein": "P1", "ligand": "X", "site_id": "site_1", "tag": "P1_site_1_X", "pose": 2,
         "affinity_kcal_mol": -9.0, "cnn_score": 0.3, "cnn_affinity": 8.0, "pose_file": ""},
        {"engine": "gnina", "protein": "P1", "ligand": "Y", "site_id": "site_1", "tag": "P1_site_1_Y", "pose": 1,
         "affinity_kcal_mol": -8.0, "cnn_score": 0.5, "cnn_affinity": 9.0, "pose_file": ""},
        {"engine": "gnina", "protein": "P1", "ligand": "Z", "site_id": "site_1", "tag": "P1_site_1_Z", "pose": 1,
         "affinity_kcal_mol": -6.0, "cnn_score": 0.7, "cnn_affinity": 6.0, "pose_file": ""},
    ]
    return pd.DataFrame(rows)


def _single_vina_fixture() -> pd.DataFrame:
    rows = [
        {"engine": "vina", "protein": "P1", "ligand": name, "site_id": "site_1", "tag": f"P1_site_1_{name}", "pose": 1,
         "affinity_kcal_mol": affinity}
        for name, affinity in (("A", -10.0), ("B", -8.0), ("C", -6.0))
    ]
    return pd.DataFrame(rows)


def _pipeline(project_dir: Path, output_dir: Path) -> MultiEngineAnalysisPipeline:
    return MultiEngineAnalysisPipeline(
        project_dir=str(project_dir),
        output_dir=str(output_dir),
        analysis_mode="comparative_all_engines",
        analysis_scope="comparison_only",
        normalization_method="per_engine_rank",
        favorite_engine="gnina",
        hit_class_policy="target_percentile",
        hit_class_strong_percentile=0.10,
        hit_class_moderate_percentile=0.40,
        top_pose_selection_policy="best_affinity",
        top_pose_global_aggregation="best_target",
    )


def _stub() -> MultiEngineAnalysisPipeline:
    """A pipeline object for the classmethods and static helpers that need no project state."""
    return object.__new__(MultiEngineAnalysisPipeline)


def _workflow_dock_argv(root: Path, *extra: str) -> list[str]:
    return [
        "dock", "run",
        "--project-dir", str(root),
        "--engines", "vina",
        "--vina-binary", "vina",
        "--seed", "100",
        "--dry-run",
        "--no-ligand-qc-gate",
        "--no-receptor-qc-gate",
        *extra,
    ]


# ----------------------------------------------------------------------------- R3a: schema


def test_basic_mode_refuses_an_explicit_exhaustiveness_with_a_clear_message() -> None:
    schema = _schema("basic", explicit=True, exhaustiveness=16)
    errors, _warnings = validate_parameter_schema(schema, ["vina"])
    assert any(
        "`--exhaustiveness` is not used in basic mode" in error
        and "'exhaustive'" in error
        and "--parameter-mode advanced" in error
        for error in errors
    )


def test_basic_mode_default_exhaustiveness_comes_from_the_preset_and_is_accepted() -> None:
    schema = _schema("basic", preset="screening_fast", exhaustiveness=32, explicit=False)
    assert schema.common.exhaustiveness == 8
    assert schema.common.exhaustiveness_source == "basic_preset:screening_fast"
    assert validate_parameter_schema(schema, ["vina"])[0] == []


def test_advanced_mode_honours_an_explicit_exhaustiveness_and_records_the_source() -> None:
    schema = _schema("advanced", exhaustiveness=8, explicit=True)
    assert schema.common.exhaustiveness == 8
    assert schema.common.exhaustiveness_source == "user_explicit"
    assert schema.common.exhaustiveness_overridden is True
    assert validate_parameter_schema(schema, ["vina"])[0] == []


def test_advanced_mode_default_exhaustiveness_is_recorded_as_the_advanced_default() -> None:
    schema = _schema("advanced", exhaustiveness=32, explicit=False)
    assert schema.common.exhaustiveness == 32
    assert schema.common.exhaustiveness_source == "advanced_default"
    assert schema.common.exhaustiveness_overridden is False


# ----------------------------------------------------------------------------- R3a: docking CLI


def test_dock_cli_refuses_an_explicit_exhaustiveness_in_basic_mode(tmp_path, monkeypatch, capsys) -> None:
    root = _dock_project(tmp_path, monkeypatch)
    assert docking_cli.dock_main(_dock_args(root, "--seed", "100", "--exhaustiveness", "16")) == 1
    out = capsys.readouterr().out
    assert "`--exhaustiveness` is not used in basic mode" in out
    assert "--parameter-mode advanced" in out


def test_dock_cli_records_the_basic_preset_source_of_exhaustiveness(tmp_path, monkeypatch, capsys) -> None:
    root = _dock_project(tmp_path, monkeypatch)
    assert docking_cli.dock_main(_dock_args(root, "--seed", "100")) == 0
    assert "exhaustiveness=32 (source: basic_preset:exhaustive)" in capsys.readouterr().out
    manifest = _run_manifest(root)
    assert manifest["effective"]["exhaustiveness"] == 32
    assert manifest["effective"]["exhaustiveness_source"] == "basic_preset:exhaustive"
    assert manifest["jobs"][0]["exhaustiveness"] == 32
    project = json.loads((root / "project_manifest.json").read_text(encoding="utf-8"))
    assert project["parameter_schema"]["common"]["exhaustiveness_source"] == "basic_preset:exhaustive"


def test_dock_cli_advanced_explicit_exhaustiveness_is_honoured_and_recorded(tmp_path, monkeypatch, capsys) -> None:
    root = _dock_project(tmp_path, monkeypatch)
    rc = docking_cli.dock_main(_dock_args(root, "--seed", "100", "--parameter-mode", "advanced", "--exhaustiveness", "8"))
    assert rc == 0
    assert "exhaustiveness=8 (source: user_explicit)" in capsys.readouterr().out
    manifest = _run_manifest(root)
    assert manifest["effective"]["exhaustiveness"] == 8
    assert manifest["effective"]["exhaustiveness_source"] == "user_explicit"
    assert manifest["jobs"][0]["exhaustiveness"] == 8
    project = json.loads((root / "project_manifest.json").read_text(encoding="utf-8"))
    assert project["parameter_schema"]["common"]["exhaustiveness_source"] == "user_explicit"


# ----------------------------------------------------------------------------- R3a: workflow dock run


def test_workflow_dock_run_refuses_an_explicit_exhaustiveness_in_basic_mode(tmp_path, monkeypatch, capsys) -> None:
    root = _dock_project(tmp_path, monkeypatch)
    assert workflow_main(_workflow_dock_argv(root, "--exhaustiveness", "16")) == 1
    assert "`--exhaustiveness` is not used in basic mode" in capsys.readouterr().out


def test_workflow_dock_run_without_the_flag_uses_the_preset_and_records_it(tmp_path, monkeypatch, capsys) -> None:
    root = _dock_project(tmp_path, monkeypatch)
    from workflow.cli import _build_dock_argv, build_parser

    argv = _build_dock_argv(build_parser().parse_args(_workflow_dock_argv(root)))
    assert "--exhaustiveness" not in argv
    assert workflow_main(_workflow_dock_argv(root)) == 0
    assert _run_manifest(root)["effective"]["exhaustiveness_source"] == "basic_preset:exhaustive"


def test_workflow_dock_run_advanced_exhaustiveness_is_forwarded_and_recorded(tmp_path, monkeypatch, capsys) -> None:
    root = _dock_project(tmp_path, monkeypatch)
    rc = workflow_main(_workflow_dock_argv(root, "--parameter-mode", "advanced", "--exhaustiveness", "8"))
    assert rc == 0
    assert "exhaustiveness=8 (source: user_explicit)" in capsys.readouterr().out
    manifest = _run_manifest(root)
    assert manifest["effective"]["exhaustiveness"] == 8
    assert manifest["effective"]["exhaustiveness_source"] == "user_explicit"


# ----------------------------------------------------------------------------- R3b: ranking


def test_ranking_method_follows_the_engine_count() -> None:
    assert ranking_method_for_scores(_two_engine_disagreement()) == CONSENSUS_V2_METHOD
    assert ranking_method_for_scores(_single_vina_fixture()) == SINGLE_ENGINE_NATIVE_METHOD


def test_two_engine_v2_composites_are_the_expected_percentiles() -> None:
    best = MultiEngineAnalysisPipeline._best_by_tag(_two_engine_disagreement()).set_index("ligand")
    assert best.loc["A", "consensus_score"] == pytest.approx(0.50)
    assert best.loc["B", "consensus_score"] == pytest.approx(0.75)
    assert best.loc["C", "consensus_score"] == pytest.approx(0.25)
    assert best.loc["A", "winner_engine"] == "vina"
    assert best.loc["B", "winner_engine"] == "gnina"
    assert best.loc["C", "winner_engine"] == "gnina"


def test_best_by_tag_orders_by_the_v2_rank_not_the_raw_minimum() -> None:
    fixture = _two_engine_disagreement()
    # The raw minimum across engines is ligand A (vina -10.0). The v2 rank does not put A first.
    assert fixture.sort_values("affinity_kcal_mol").iloc[0]["ligand"] == "A"
    best = MultiEngineAnalysisPipeline._best_by_tag(fixture)
    assert list(best["ligand"]) == ["B", "A", "C"]
    assert list(best["engine"]) == ["gnina", "vina", "gnina"]
    assert list(best["consensus_rank_global"]) == [1, 2, 3]
    assert set(best["ranking_method"]) == {CONSENSUS_V2_METHOD}


def test_downstream_summaries_follow_the_v2_rank() -> None:
    results = MultiEngineAnalysisPipeline._build_downstream_results(_stub(), _two_engine_disagreement())

    top = results["top_overall"]
    assert list(top["ligand"]) == ["B", "A", "C"]
    assert top.iloc[0]["engine"] == "gnina"
    assert top.iloc[0]["consensus_rank_global"] == 1
    assert set(top["ranking_method"]) == {CONSENSUS_V2_METHOD}

    per_protein = results["best_per_protein"]
    assert len(per_protein) == 1
    assert per_protein.iloc[0]["ligand"] == "B"
    assert per_protein.iloc[0]["vina_affinity"] == pytest.approx(-7.0)

    assert list(results["best_per_ligand"]["ligand"]) == ["B", "A", "C"]

    protein_summary = results["protein_summary"].iloc[0]
    # best_affinity is the affinity of the top-ranked row (-7.0), not the raw minimum across engines (-10.0).
    assert protein_summary["best_affinity"] == pytest.approx(-7.0)
    assert protein_summary["top_consensus_rank_global"] == 1
    assert protein_summary["ranking_method"] == CONSENSUS_V2_METHOD

    assert list(results["ligand_summary"]["ligand"]) == ["B", "A", "C"]


def test_a_tag_without_a_complete_v2_case_is_marked_and_never_ranked_by_raw_affinity() -> None:
    fixture = pd.concat(
        [
            _two_engine_disagreement(),
            # Vina only, with the lowest raw affinity of all. It has no gnina row, so its v2 case is incomplete.
            pd.DataFrame([{"engine": "vina", "protein": "P1", "ligand": "D", "site_id": "site_1",
                           "tag": "P1_site_1_D", "pose": 1, "affinity_kcal_mol": -11.0}]),
        ],
        ignore_index=True,
    )
    best = MultiEngineAnalysisPipeline._best_by_tag(fixture)
    row = best[best["ligand"] == "D"].iloc[0]
    assert row["ranking_method"] == RANKING_UNAVAILABLE_CONSENSUS_INCOMPLETE
    assert pd.isna(row["consensus_rank_global"])
    assert pd.isna(row["consensus_score"])
    # Ranked rows come first; the unranked row is last.
    assert list(best["ligand"]) == ["B", "A", "C", "D"]

    results = MultiEngineAnalysisPipeline._build_downstream_results(_stub(), fixture)
    assert "D" not in set(results["top_overall"]["ligand"])
    assert "D" not in set(results["best_per_ligand"]["ligand"])
    assert results["best_per_protein"].iloc[0]["ligand"] == "B"


def test_single_engine_vina_ranks_by_single_engine_native_v1() -> None:
    best = MultiEngineAnalysisPipeline._best_by_tag(_single_vina_fixture())
    assert list(best["ligand"]) == ["A", "B", "C"]
    assert set(best["ranking_method"]) == {SINGLE_ENGINE_NATIVE_METHOD}
    assert set(best["consensus_status"]) == {"single_engine"}


def test_single_engine_gnina_picks_by_cnn_score_and_ranks_by_cnn_affinity() -> None:
    best = MultiEngineAnalysisPipeline._best_by_tag(_single_gnina_selection_fixture())
    assert list(best["ligand"]) == ["Y", "Z", "X"]
    assert set(best["ranking_method"]) == {SINGLE_ENGINE_NATIVE_METHOD}
    x_row = best[best["ligand"] == "X"].iloc[0]
    assert int(x_row["pose"]) == 1  # highest cnn_score, not the highest cnn_affinity (pose 2)


def test_simplified_bridge_picks_the_consensus_row_per_protein_ligand() -> None:
    # Same protein and ligand, two engines: the raw minimum is vina -10.0, but the consensus rank is gnina.
    frame = pd.DataFrame(
        [
            {"protein": "P1", "ligand": "A", "tag": "P1_site_1_A", "engine": "vina", "vina_affinity": -10.0,
             "consensus_rank_global": 3.0},
            {"protein": "P1", "ligand": "A", "tag": "P1_site_1_A", "engine": "gnina", "vina_affinity": -6.0,
             "consensus_rank_global": 1.0},
        ]
    )
    picked = SimplifiedPostDockingPipeline._best_pose_per_protein_ligand(frame)
    assert list(picked["engine"]) == ["gnina"]


# ----------------------------------------------------------------------------- R3b/R3c: single-engine writer


def test_single_engine_gnina_summary_names_cnn_score_as_the_selection_key(tmp_path) -> None:
    root = tmp_path / "project"
    execution.run_workflow_init(str(root), layout_profile="canonical")
    pipeline = _pipeline(root, tmp_path / "analysis")
    target = tmp_path / "analysis" / "single_gnina"
    pipeline._write_single_engine_reports(_with_unified_columns(_single_gnina_selection_fixture()), "gnina", target)

    summary = (target / "summary.txt").read_text(encoding="utf-8")
    assert "Primary ranking score: cnn_score (pose selection key" in summary
    assert "Ligand ranking metric: cnn_affinity" in summary
    assert "Primary ranking score: cnn_affinity" not in summary
    assert "Ranking method (Spec 037 R3b): single_engine_native_v1" in summary
    assert "Top ranked complex (single_engine_native_v1): P1_site_1_Y" in summary

    best = pd.read_csv(target / "best_poses.csv")
    assert list(best["tag"]) == ["P1_site_1_Y", "P1_site_1_Z", "P1_site_1_X"]
    assert set(best["ranking_method"]) == {SINGLE_ENGINE_NATIVE_METHOD}
    assert list(best["consensus_rank_global"]) == [1, 2, 3]

    protein = pd.read_csv(target / "protein_summary.csv")
    # best_affinity is the affinity of the top-ranked complex (Y, -8.0), not the minimum over the protein.
    assert protein.iloc[0]["best_affinity"] == pytest.approx(-8.0)
    assert protein.iloc[0]["ranking_method"] == SINGLE_ENGINE_NATIVE_METHOD


def test_single_engine_vina_summary_names_the_affinity_selection_key(tmp_path) -> None:
    root = tmp_path / "project"
    execution.run_workflow_init(str(root), layout_profile="canonical")
    pipeline = _pipeline(root, tmp_path / "analysis")
    target = tmp_path / "analysis" / "single_vina"
    pipeline._write_single_engine_reports(_with_unified_columns(_single_vina_fixture()), "vina", target)
    summary = (target / "summary.txt").read_text(encoding="utf-8")
    assert "Primary ranking score: affinity_kcal_mol (pose selection key; lowest is kept)" in summary
    assert "cnn_" not in summary
