"""Spec 038 R3 (T004): one post-docking path, checked against the current code.

Status (2026-10-19): R3 is STOPPED for Lead decision. The checks below pass where the two paths
already agree, and each known divergence is a strict xfail that names its STOP point in
specs/038-entry-point-consolidation/postdocking_inventory.md (section 4). When a divergence is
resolved, its xfail turns into a failure and must be flipped to a plain test.
"""
from __future__ import annotations

import inspect
import io
import json
from pathlib import Path

import pandas as pd
import pytest

import post_docking_analysis.simplified_pipeline_impl as spi
from docking.project_layout import bootstrap_project_layout, ensure_engine_layout, pairlist_path, post_docking_root
from post_docking_analysis.multi_engine_pipeline import MultiEngineAnalysisPipeline
from post_docking_analysis.unified_pipeline import UnifiedPostDockingPipeline
from test_dockforge_smoke import (
    _seed_engine_outputs,
    _stub_dag_optional_interactions,
    _write_test_pairlist,
)

ENGINES3 = ("gnina", "vina", "smina")
CONSENSUS_FILES = (
    "consensus_ranked_hits.csv",
    "engine_summary.csv",
    "best_engine_per_complex.csv",
    "combined_engine_scores.csv",
)
# Columns that describe how the engine scope was resolved, not the score or pose itself.
SCOPE_COLUMNS = {"engines_in_scope", "scope_source", "scoped_engine_count", "detected_engine_count", "excluded_engines"}
STOP_SUMMARY = "STOP D-A1 (Spec 038 R3): DAG reports/summary.txt is a placeholder; legacy writes the full summary."
STOP_ATLAS = "STOP D-A2 (Spec 038 R3): DAG top-pose atlas drops normalized_affinity_score and the status columns."
STOP_FAVORITE = "STOP D-A3 (Spec 038 R3): the DAG path has no favorite_engine_continue stage (legacy writes favorite_engine/)."
STOP_ROUTE = "STOP D-A (Spec 038 R3): run() is not yet routed through the DAG; that changes legacy outputs."
STOP_FILTER = "STOP D-B (Spec 038 R3): the DAG ignores the problematic-ligand filter; legacy applies it when enabled."
STOP_WRITER = "Single-writer violation (Spec 038 R3): DAG nodes share reports/ paths with the comparative writer."


def _seed_three_engines(root: Path) -> None:
    bootstrap_project_layout(root, list(ENGINES3), layout_profile="docking_legacy")
    _write_test_pairlist(root)
    _seed_engine_outputs(root, ENGINES3)


def _seed_two_tags_vina(root: Path) -> None:
    """One vina project with a usable ligand (L1) and a positive-affinity ligand (L4)."""
    bootstrap_project_layout(root, ["vina"], layout_profile="docking_legacy")
    rows = []
    for ligand in ("L1", "L4"):
        rows.append({
            "receptor": "R1", "site_id": "site_1", "ligand": ligand,
            "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
            "size_x": 20.0, "size_y": 20.0, "size_z": 20.0,
        })
    pd.DataFrame(rows).to_csv(pairlist_path(root), index=False)
    layout = ensure_engine_layout(root, "vina")
    for ligand, affinity in (("L1", -7.2), ("L4", 0.4)):
        tag = f"R1_site_1_{ligand}"
        (layout["poses"] / f"{tag}.pdbqt").write_text(
            f"MODEL 1\nREMARK VINA RESULT: {affinity:.2f} 0.00 0.00\nENDMDL\n", encoding="utf-8"
        )
        (layout["logs"] / f"{tag}.log").write_text(
            "mode | affinity | rmsd_lb | rmsd_ub\n"
            f"    1      {affinity:.2f}      0.00      0.00\n",
            encoding="utf-8",
        )


def _dag(root: Path, session: Path, **kwargs) -> UnifiedPostDockingPipeline:
    pipeline = UnifiedPostDockingPipeline(
        project_dir=str(root), output_dir=str(session), rmsd_workers=1, **kwargs
    )
    _stub_dag_optional_interactions(pipeline)
    return pipeline


def _legacy(root: Path, session: Path, engines: list, **kwargs) -> MultiEngineAnalysisPipeline:
    return MultiEngineAnalysisPipeline(
        project_dir=str(root), output_dir=str(session), rmsd_workers=1, engines_in_scope=engines, **kwargs
    )


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    """Each run gets its own project copy, because the DAG writes project-level metadata."""
    base = tmp_path_factory.mktemp("spec038_postdoc")
    results = {}

    def project(name: str, seeder) -> tuple[Path, Path]:
        root = base / name / "project"
        seeder(root)
        return root, post_docking_root(root) / "sessions" / name

    def record(name: str, fn):
        results[name]["error"] = ""
        try:
            fn()
        except Exception as exc:  # the legacy smina crash is a result to test, not a fixture failure
            results[name]["error"] = f"{type(exc).__name__}: {exc}"

    comparative_kwargs = dict(analysis_mode="comparative_all_engines", analysis_scope="full")

    root, session = project("dag_comparative", _seed_three_engines)
    results["dag_comparative"] = {"root": root, "session": session}
    record("dag_comparative", lambda: _dag(root, session, engines="gnina,vina,smina", **comparative_kwargs).run())

    root, session = project("legacy_comparative", _seed_three_engines)
    results["legacy_comparative"] = {"root": root, "session": session}
    record("legacy_comparative", lambda: _legacy(root, session, list(ENGINES3), **comparative_kwargs).run())

    root, session = project("dag_single_vina", _seed_three_engines)
    results["dag_single_vina"] = {"root": root, "session": session}
    record("dag_single_vina", lambda: _dag(root, session, analysis_mode="single_engine", engine="vina", engines="vina").run())

    root, session = project("legacy_single_vina", _seed_three_engines)
    results["legacy_single_vina"] = {"root": root, "session": session}
    record("legacy_single_vina", lambda: _legacy(root, session, ["vina"], analysis_mode="single_engine", engine="vina").run())

    root, session = project("legacy_single_smina", _seed_three_engines)
    results["legacy_single_smina"] = {"root": root, "session": session}
    record("legacy_single_smina", lambda: _legacy(root, session, ["smina"], analysis_mode="single_engine", engine="smina").run())

    root, session = project("dag_favorite", _seed_three_engines)
    results["dag_favorite"] = {"root": root, "session": session}
    record("dag_favorite", lambda: _dag(root, session, analysis_mode="favorite_engine_continue", favorite_engine="vina",
                                        engines="gnina,vina,smina", analysis_scope="full").run())

    root, session = project("legacy_favorite", _seed_three_engines)
    results["legacy_favorite"] = {"root": root, "session": session}
    record("legacy_favorite", lambda: _legacy(root, session, list(ENGINES3), analysis_mode="favorite_engine_continue",
                                              favorite_engine="vina", analysis_scope="full").run())

    root, session = project("legacy_filter_on", _seed_two_tags_vina)
    results["legacy_filter_on"] = {"root": root, "session": session}
    record("legacy_filter_on", lambda: _legacy(root, session, ["vina"], analysis_mode="single_engine", engine="vina",
                                               exclude_problematic_ligands=True).run())

    root, session = project("legacy_filter_off", _seed_two_tags_vina)
    results["legacy_filter_off"] = {"root": root, "session": session}
    record("legacy_filter_off", lambda: _legacy(root, session, ["vina"], analysis_mode="single_engine", engine="vina").run())

    root, session = project("dag_filter_on", _seed_two_tags_vina)
    results["dag_filter_on"] = {"root": root, "session": session}
    record("dag_filter_on", lambda: _dag(root, session, analysis_mode="single_engine", engine="vina", engines="vina",
                                         exclude_problematic_ligands=True).run())
    return results


def _normalized_text(path: Path, root: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace").replace(str(root), "<P>")


def _frame(path: Path, root: Path) -> pd.DataFrame:
    frame = pd.read_csv(io.StringIO(_normalized_text(path, root)), dtype=str, keep_default_na=False)
    return frame.drop(columns=[c for c in frame.columns if c in SCOPE_COLUMNS])


# ----------------------------------------------------------------- passing: what already agrees

def test_compat_entry_points_import():
    assert issubclass(UnifiedPostDockingPipeline, MultiEngineAnalysisPipeline)
    assert hasattr(spi, "SimplifiedPostDockingPipeline")
    from post_docking_analysis import simplified_cli, cli

    assert callable(simplified_cli.main) and callable(cli.main)


def test_problematic_filter_default_is_off_on_every_entry_point():
    # UnifiedPostDockingPipeline forwards its keyword arguments to MultiEngineAnalysisPipeline.
    for cls in (MultiEngineAnalysisPipeline, spi.SimplifiedPostDockingPipeline):
        default = inspect.signature(cls.__init__).parameters["exclude_problematic_ligands"].default
        assert default is False, cls.__name__


def test_consensus_outputs_identical_between_run_and_dag(runs):
    for name in CONSENSUS_FILES:
        dag_file = runs["dag_comparative"]["session"] / "reports" / name
        legacy_file = runs["legacy_comparative"]["session"] / "reports" / name
        assert dag_file.exists() and legacy_file.exists(), name
        dag_frame = _frame(dag_file, runs["dag_comparative"]["root"])
        legacy_frame = _frame(legacy_file, runs["legacy_comparative"]["root"])
        pd.testing.assert_frame_equal(dag_frame, legacy_frame, check_like=False, obj=name)


def test_single_engine_values_identical_between_run_and_dag(runs):
    for rel in ("vina/best_poses.csv", "vina/scores.csv", "best_pose_table/best_pose_per_tag_by_engine.csv"):
        dag_file = runs["dag_single_vina"]["session"] / rel
        legacy_file = runs["legacy_single_vina"]["session"] / rel
        assert dag_file.exists() and legacy_file.exists(), rel
        pd.testing.assert_frame_equal(
            _frame(dag_file, runs["dag_single_vina"]["root"]),
            _frame(legacy_file, runs["legacy_single_vina"]["root"]),
            obj=rel,
        )


def test_legacy_favorite_engine_stage_writes_its_directory(runs):
    assert (runs["legacy_favorite"]["session"] / "favorite_engine" / "summary.txt").exists()


def test_legacy_filter_applied_when_enabled_and_not_when_disabled(runs):
    excluded = runs["legacy_filter_on"]["session"] / "excluded_problematic_ligands.csv"
    assert excluded.exists()
    assert "R1_site_1_L4" in excluded.read_text(encoding="utf-8")
    assert not (runs["legacy_filter_off"]["session"] / "excluded_problematic_ligands.csv").exists()


# ----------------------------------------------------------------- STOP points (strict xfail)

def test_legacy_single_engine_smina_completes(runs):
    assert runs["legacy_single_smina"]["error"] == ""
    assert (runs["legacy_single_smina"]["session"] / "smina" / "best_poses.csv").exists()


def test_dag_reports_summary_is_the_full_summary(runs):
    summary = (runs["dag_comparative"]["session"] / "reports" / "summary.txt").read_text(encoding="utf-8")
    assert "Multi-engine docking summary" in summary


@pytest.mark.xfail(strict=True, reason=STOP_ATLAS)
def test_top_pose_normalized_score_matches_legacy(runs):
    name = "top_pose_ligand_performance/top_pose_per_ligand_global.csv"
    dag = pd.read_csv(runs["dag_comparative"]["session"] / name)
    legacy = pd.read_csv(runs["legacy_comparative"]["session"] / name)
    pd.testing.assert_series_equal(dag["normalized_affinity_score"], legacy["normalized_affinity_score"], check_names=False)


@pytest.mark.xfail(strict=True, reason=STOP_ATLAS)
def test_top_pose_status_columns_present_on_dag(runs):
    name = "top_pose_ligand_performance/top_pose_per_ligand_global.csv"
    dag_cols = set(pd.read_csv(runs["dag_comparative"]["session"] / name, nrows=0).columns)
    legacy_cols = set(pd.read_csv(runs["legacy_comparative"]["session"] / name, nrows=0).columns)
    assert {"redocking_validation_status", "pose_reproducibility_status", "docked_microspecies_smiles"} <= dag_cols
    assert legacy_cols - SCOPE_COLUMNS <= dag_cols


@pytest.mark.xfail(strict=True, reason=STOP_FAVORITE)
def test_favorite_engine_stage_present_on_dag(runs):
    assert (runs["dag_favorite"]["session"] / "favorite_engine" / "summary.txt").exists()


@pytest.mark.xfail(strict=True, reason=STOP_FILTER)
def test_problematic_filter_recorded_on_dag(runs):
    analysis_parameters = runs["dag_filter_on"]["root"] / "4-Working" / "metadata" / "analysis_parameters.json"
    assert analysis_parameters.exists()
    payload = json.loads(analysis_parameters.read_text(encoding="utf-8"))
    assert payload.get("exclude_problematic_ligands") is True


@pytest.mark.xfail(strict=True, reason=STOP_ROUTE)
def test_run_routes_through_the_dag():
    assert "request_artifact_scope" in inspect.getsource(MultiEngineAnalysisPipeline.run)


def test_report_nodes_do_not_call_the_comparative_writer():
    from post_docking_analysis import multi_engine_pipeline_impl as mep

    for method in ("_dag_compute_polypharmacology_node", "_dag_compute_biology_correlation_node"):
        assert "_write_comparative_reports(" not in inspect.getsource(getattr(mep.MultiEngineAnalysisPipeline, method))


def test_dag_run_retires_no_shared_report_files(runs):
    root = runs["dag_comparative"]["root"]
    retired = [
        p.relative_to(root).as_posix()
        for p in root.rglob("*")
        if p.is_file() and ".artifact_history" in p.parts
        and (p.name == "summary.txt" or "polypharmacology" in p.parts or "biology_correlation" in p.parts)
    ]
    assert retired == []


@pytest.mark.xfail(strict=True, reason=STOP_ROUTE)
def test_run_and_dag_produce_the_same_session_outputs(runs):
    """Every text file the legacy run writes must match the DAG output (run_tracking excluded)."""
    legacy_root, legacy_session = runs["legacy_comparative"]["root"], runs["legacy_comparative"]["session"]
    dag_root, dag_session = runs["dag_comparative"]["root"], runs["dag_comparative"]["session"]
    mismatched = []
    for path in sorted(legacy_session.rglob("*")):
        rel = path.relative_to(legacy_session)
        if not path.is_file() or rel.parts[0] in {"run_tracking", "favorite_engine"} or path.suffix not in {".csv", ".json", ".txt"}:
            continue
        other = dag_session / rel
        if not other.exists() or _normalized_text(path, legacy_root) != _normalized_text(other, dag_root):
            mismatched.append(rel.as_posix())
    assert mismatched == []


# ----------------------------------------------------------------- D-C: GNINA best pose via the shared selector

def _simplified_with_scores(scores: pd.DataFrame):
    pipeline = object.__new__(spi.SimplifiedPostDockingPipeline)
    pipeline.scores_df = scores
    return pipeline


def test_simplified_gnina_best_pose_is_highest_cnn_score_not_lowest_vina():
    scores = pd.DataFrame([
        {"engine": "gnina", "tag": "R1_site_1_L2", "protein": "R1", "pose": 1, "vina_affinity": -7.5, "cnn_score": 0.55, "cnn_affinity": -7.1, "affinity_kcal_mol": -7.5},
        {"engine": "gnina", "tag": "R1_site_1_L2", "protein": "R1", "pose": 2, "vina_affinity": -7.0, "cnn_score": 0.80, "cnn_affinity": -6.6, "affinity_kcal_mol": -7.0},
        {"engine": "gnina", "tag": "R1_site_1_L1", "protein": "R1", "pose": 1, "vina_affinity": -8.6, "cnn_score": 0.91, "cnn_affinity": -8.2, "affinity_kcal_mol": -8.6},
    ])
    best = _simplified_with_scores(scores)._best_poses_from_scores()
    chosen = {row["tag"]: int(row["pose"]) for _, row in best.iterrows()}
    assert chosen == {"R1_site_1_L2": 2, "R1_site_1_L1": 1}
    # The output keeps the input rows and columns unchanged.
    assert list(best.columns) == list(scores.columns)


def test_simplified_vina_best_pose_is_lowest_affinity():
    scores = pd.DataFrame([
        {"engine": "vina", "tag": "R1_site_1_L2", "protein": "R1", "pose": 1, "vina_affinity": -7.5, "affinity_kcal_mol": -7.5},
        {"engine": "vina", "tag": "R1_site_1_L2", "protein": "R1", "pose": 2, "vina_affinity": -7.0, "affinity_kcal_mol": -7.0},
    ])
    best = _simplified_with_scores(scores)._best_poses_from_scores()
    assert best["pose"].tolist() == [1]
