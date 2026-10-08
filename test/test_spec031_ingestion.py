from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
import pytest

from post_docking_analysis.engine_detector import detect_engines
from post_docking_analysis.engine_hpc_adapter import detect_engine_layout
from post_docking_analysis.multi_engine_pipeline import MultiEngineAnalysisPipeline
from post_docking_analysis.unified_pipeline import UnifiedPostDockingPipeline


def _tag() -> str:
    return "REC.pdbqt_site_1_LIG.pdbqt"


def _write_project(project: Path, engine: str) -> None:
    dock_root = project / "4-Docking"
    dock_root.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                "receptor": "REC.pdbqt",
                "site_id": "site_1",
                "ligand": "LIG.pdbqt",
                "center_x": 1.0,
                "center_y": 2.0,
                "center_z": 3.0,
                "size_x": 20.0,
                "size_y": 20.0,
                "size_z": 20.0,
            }
        ]
    ).to_csv(dock_root / "pairlist.csv", index=False)
    (dock_root / "project_manifest.json").write_text(
        json.dumps(
            {
                "project_name": "spec031-ingestion",
                "project_root": str(project),
                "layout_profile": "docking_legacy",
                "engines": [engine],
                "pairlist_file": str(dock_root / "pairlist.csv"),
            }
        ),
        encoding="utf-8",
    )


def _vina_payload(affinity: float) -> str:
    return "\n".join(
        [
            "MODEL 1",
            f"REMARK VINA RESULT: {affinity:.2f} 0.000 0.000",
            "ATOM      1  C1  LIG A   1       0.000   0.000   0.000  0.00  0.00     0.000 C",
            "ENDMDL",
            "",
        ]
    )


def _gnina_log_payload(affinity: float) -> str:
    return "\n".join(
        [
            "mode | affinity | intramol | CNN pose score | CNN affinity",
            "-----+------------+------------+------------+----------",
            f"   1      {affinity:.2f}       -0.50       0.8000      -7.50",
            "",
        ]
    )


def _pipeline(project: Path, engine: str, *, force: bool = False) -> MultiEngineAnalysisPipeline:
    return MultiEngineAnalysisPipeline(
        project_dir=str(project),
        output_dir=str(project / "analysis"),
        analysis_mode="single_engine",
        engine=engine,
        engines_in_scope=[engine],
        force_score_rebuild=force,
    )


def _make_newer(path: Path, reference: Path) -> None:
    future_ns = max(path.stat().st_mtime_ns, reference.stat().st_mtime_ns) + 2_000_000_000
    os.utime(path, ns=(future_ns, future_ns))


def test_read_only_resolver_selects_populated_artifacts_and_reports_ambiguity(tmp_path: Path) -> None:
    project = tmp_path / "project"
    nested = project / "4-Docking" / "vina_out" / "poses"
    flat = project / "4-Docking" / "vina_out"
    nested.mkdir(parents=True)
    (nested / "nested.pdbqt").write_text(_vina_payload(-7.0), encoding="utf-8")
    (flat / "flat.pdbqt").write_text(_vina_payload(-8.0), encoding="utf-8")
    before = sorted(str(path.relative_to(project)) for path in project.rglob("*"))

    detected = detect_engine_layout(project, "vina")

    after = sorted(str(path.relative_to(project)) for path in project.rglob("*"))
    assert before == after
    assert detected["pose_folder"] == nested.resolve()
    assert detected["ambiguous"] is True
    pose_resolution = detected["artifact_resolution"]["poses"]
    assert pose_resolution["status"] == "ambiguous"
    assert pose_resolution["populated_candidates"] == [str(nested.resolve()), str(flat.resolve())]
    assert "Multiple populated poses directories" in detected["warnings"][0]


def test_unified_ingestion_reads_flat_legacy_vina_output(tmp_path: Path) -> None:
    project = tmp_path / "project"
    _write_project(project, "vina")
    flat_pose = project / "4-Docking" / "vina_out" / f"{_tag()}.pdbqt"
    flat_pose.parent.mkdir(parents=True)
    flat_pose.write_text(_vina_payload(-7.25), encoding="utf-8")

    scores = _pipeline(project, "vina")._load_or_build_scores()

    assert len(scores) == 1
    assert scores.iloc[0]["engine"] == "vina"
    assert scores.iloc[0]["tag"] == _tag()
    assert scores.iloc[0]["affinity_kcal_mol"] == pytest.approx(-7.25)
    assert Path(scores.iloc[0]["pose_file"]) == flat_pose.resolve()


def test_unified_ingestion_resolves_numbered_pose_with_legacy_log(tmp_path: Path) -> None:
    project = tmp_path / "project"
    _write_project(project, "vina")
    numbered_pose = project / "3-Docking" / "vina" / "poses" / f"{_tag()}.pdbqt"
    legacy_log = project / "4-Docking" / "logs" / f"{_tag()}.log"
    numbered_pose.parent.mkdir(parents=True)
    legacy_log.parent.mkdir(parents=True)
    numbered_pose.write_text(_vina_payload(-7.75), encoding="utf-8")
    legacy_log.write_text("Vina execution log\n", encoding="utf-8")

    detected = detect_engine_layout(project, "vina")
    scores = _pipeline(project, "vina")._load_or_build_scores()

    assert detected["pose_folder"] == numbered_pose.parent.resolve()
    assert detected["log_folder"] == legacy_log.parent.resolve()
    assert scores.iloc[0]["affinity_kcal_mol"] == pytest.approx(-7.75)


def test_ingestion_fails_explicitly_when_raw_outputs_parse_to_zero_rows(tmp_path: Path) -> None:
    project = tmp_path / "project"
    _write_project(project, "vina")
    pose = project / "4-Docking" / "vina_out" / f"{_tag()}.pdbqt"
    pose.parent.mkdir(parents=True)
    pose.write_text("MODEL 1\nREMARK no score here\nENDMDL\n", encoding="utf-8")

    with pytest.raises(ValueError, match="produced zero valid rows") as exc_info:
        _pipeline(project, "vina")._load_or_build_scores()

    message = str(exc_info.value)
    assert "engine=vina" in message
    assert "parser=parse_vina_pdbqt" in message
    assert str(pose.resolve()) in message


def test_shared_log_from_another_engine_is_not_treated_as_vina_pose_output(tmp_path: Path) -> None:
    project = tmp_path / "project"
    _write_project(project, "vina")
    shared_log = project / "4-Docking" / "logs" / "gnina_only.log"
    shared_log.parent.mkdir(parents=True)
    shared_log.write_text(_gnina_log_payload(-7.0), encoding="utf-8")

    scores = _pipeline(project, "vina")._load_or_build_scores()

    assert scores.empty


def test_normalized_vina_cache_rebuilds_when_pose_changes(tmp_path: Path) -> None:
    project = tmp_path / "project"
    _write_project(project, "vina")
    pose = project / "4-Docking" / "vina_out" / "poses" / f"{_tag()}.pdbqt"
    pose.parent.mkdir(parents=True)
    pose.write_text(_vina_payload(-7.0), encoding="utf-8")
    pipeline = _pipeline(project, "vina")
    first = pipeline._load_or_build_scores()
    normalized = project / "4-Docking" / "vina_out" / "scores" / "normalized_scores.csv"
    assert first.iloc[0]["affinity_kcal_mol"] == pytest.approx(-7.0)
    assert normalized.exists()

    pose.write_text(_vina_payload(-8.5), encoding="utf-8")
    _make_newer(pose, normalized)
    second = pipeline._load_or_build_scores()

    assert second.iloc[0]["affinity_kcal_mol"] == pytest.approx(-8.5)
    assert pd.read_csv(normalized).iloc[0]["affinity_kcal_mol"] == pytest.approx(-8.5)

    pose.write_text(_vina_payload(-9.5), encoding="utf-8")
    older_ns = max(1, normalized.stat().st_mtime_ns - 2_000_000_000)
    os.utime(pose, ns=(older_ns, older_ns))
    cached = pipeline._load_or_build_scores()
    forced = _pipeline(project, "vina", force=True)._load_or_build_scores()
    assert cached.iloc[0]["affinity_kcal_mol"] == pytest.approx(-8.5)
    assert forced.iloc[0]["affinity_kcal_mol"] == pytest.approx(-9.5)


def test_gnina_all_scores_is_refreshed_with_changed_log(tmp_path: Path) -> None:
    project = tmp_path / "project"
    _write_project(project, "gnina")
    pose = project / "4-Docking" / "gnina_out" / f"{_tag()}.sdf"
    log = project / "4-Docking" / "logs" / f"{_tag()}.log"
    pose.parent.mkdir(parents=True)
    log.parent.mkdir(parents=True)
    pose.write_text("$$$$\n", encoding="utf-8")
    log.write_text(_gnina_log_payload(-7.0), encoding="utf-8")
    pipeline = _pipeline(project, "gnina")
    first = pipeline._load_or_build_scores()
    normalized = project / "4-Docking" / "results" / "normalized_scores.csv"
    all_scores = project / "4-Docking" / "results" / "all_scores.csv"
    assert first.iloc[0]["affinity_kcal_mol"] == pytest.approx(-7.0)
    assert all_scores.exists()

    log.write_text(_gnina_log_payload(-9.0), encoding="utf-8")
    _make_newer(log, normalized)
    second = pipeline._load_or_build_scores()

    assert second.iloc[0]["affinity_kcal_mol"] == pytest.approx(-9.0)
    assert pd.read_csv(all_scores).iloc[0]["vina_affinity"] == pytest.approx(-9.0)


def test_dag_raw_scores_declares_resolved_pose_log_and_pair_metadata_inputs(tmp_path: Path) -> None:
    project = tmp_path / "project"
    _write_project(project, "vina")
    pose = project / "4-Docking" / "vina_out" / f"{_tag()}.pdbqt"
    log = project / "4-Docking" / "logs" / f"{_tag()}.log"
    intent = project / "4-Docking" / "metadata" / "pair_intent.csv"
    pose.parent.mkdir(parents=True)
    log.parent.mkdir(parents=True)
    intent.parent.mkdir(parents=True)
    pose.write_text(_vina_payload(-7.0), encoding="utf-8")
    log.write_text("Vina execution log\n", encoding="utf-8")
    pd.read_csv(project / "4-Docking" / "pairlist.csv").to_csv(intent, index=False)

    graph = _pipeline(project, "vina").build_artifact_graph()
    inputs = set(graph.nodes["raw_scores"].inputs)
    initial_cache_key = graph._compute_cache_key(graph.nodes["raw_scores"])

    assert str(pose.resolve()) in inputs
    assert str(log.resolve()) in inputs
    assert str(intent.resolve()) in inputs
    assert str((project / "4-Docking" / "pairlist.csv").resolve()) in inputs
    pose.write_text(_vina_payload(-8.0), encoding="utf-8")
    _make_newer(pose, intent)
    assert graph._compute_cache_key(graph.nodes["raw_scores"]) != initial_cache_key


def test_engine_detection_cache_invalidates_on_raw_output_change_and_force_propagates(tmp_path: Path) -> None:
    project = tmp_path / "project"
    _write_project(project, "vina")
    pose = project / "4-Docking" / "vina_out" / f"{_tag()}.pdbqt"
    pose.parent.mkdir(parents=True)
    pose.write_text(_vina_payload(-7.0), encoding="utf-8")

    first = detect_engines(project, redetect=True)
    cached = detect_engines(project)
    assert first["cache_hit"] is False
    assert cached["cache_hit"] is True
    assert cached["source_signature"] == first["source_signature"]

    pose.write_text(_vina_payload(-8.0), encoding="utf-8")
    _make_newer(pose, project / "4-Docking" / "project_manifest.json")
    invalidated = detect_engines(project)
    assert invalidated["cache_hit"] is False
    assert invalidated["source_signature"] != first["source_signature"]

    forced = UnifiedPostDockingPipeline(
        project_dir=str(project),
        output_dir=str(project / "forced-analysis"),
        analysis_mode="single_engine",
        engine="vina",
        dag_force=True,
    )
    assert forced.dag_force is True
    assert forced.force_score_rebuild is True
    assert forced.engine_detection_report["cache_hit"] is False
