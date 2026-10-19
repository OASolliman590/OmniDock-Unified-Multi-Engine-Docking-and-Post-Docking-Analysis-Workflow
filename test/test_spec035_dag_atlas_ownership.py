"""Spec 035 follow-up 1 (O6/S9): the analysis DAG's top_pose_atlas node is the only writer of the atlas.

The comparative, polypharmacology and biology report nodes read the atlas and never rewrite it.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest

import post_docking_analysis.multi_engine_pipeline_impl as mep
from docking.project_layout import bootstrap_project_layout, post_docking_root
from post_docking_analysis.unified_pipeline import UnifiedPostDockingPipeline
from test_dockforge_smoke import _seed_engine_outputs, _stub_dag_optional_interactions, _write_test_pairlist

ATLAS_NODE = "top_pose_atlas"
REPORT_NODES = ("comparative", "polypharmacology", "biology_correlation")


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _atlas_files(pipeline: UnifiedPostDockingPipeline) -> list[Path]:
    paths = pipeline._dag_artifact_paths()
    return [paths[key] for key in ("best_poses", "ligand_performance_summary", "best_poses_per_protein",
                                   "top_pose_confidence", "top_pose_manifest")]


def _dag_pipeline(root: Path, name: str) -> UnifiedPostDockingPipeline:
    bootstrap_project_layout(root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
    _write_test_pairlist(root)
    _seed_engine_outputs(root, ("gnina", "vina", "smina"))
    pipeline = UnifiedPostDockingPipeline(
        project_dir=str(root),
        output_dir=str(post_docking_root(root) / "sessions" / name),
        analysis_mode="comparative_all_engines",
        analysis_scope="full",
        engines="gnina,vina,smina",
        rmsd_workers=4,
    )
    _stub_dag_optional_interactions(pipeline)
    return pipeline


def test_atlas_outputs_have_one_registered_owner(tmp_path):
    pipeline = _dag_pipeline(tmp_path, "owner_registry")
    graph = pipeline.build_artifact_graph()
    for path in _atlas_files(pipeline):
        assert graph.output_to_node[str(path)] == ATLAS_NODE
    for name in REPORT_NODES:
        declared_outputs = set(graph.nodes[name].outputs)
        assert not declared_outputs & {str(path) for path in _atlas_files(pipeline)}
        # Readers depend on every atlas file, so they run after the owner and see its final bytes.
        assert {str(path) for path in _atlas_files(pipeline)} <= set(graph.nodes[name].inputs)


def test_report_nodes_run_after_the_atlas_node(tmp_path):
    pipeline = _dag_pipeline(tmp_path, "owner_order")
    graph = pipeline.build_artifact_graph()
    tiers = graph._topo_sort(str(pipeline._dag_artifact_paths()["reports"]))
    tier_of = {node.name: index for index, tier in enumerate(tiers) for node in tier}
    for name in (ATLAS_NODE, *REPORT_NODES):
        assert name in tier_of
    assert all(tier_of[ATLAS_NODE] < tier_of[name] for name in REPORT_NODES)


def test_only_the_atlas_node_calls_the_atlas_writer(tmp_path, monkeypatch):
    pipeline = _dag_pipeline(tmp_path, "owner_writer")
    real_writer = mep.write_top_pose_atlas
    callers: list[str] = []

    def recording_writer(payload, output_dir):
        # The immediate caller is the function that performs the write.
        callers.append(sys._getframe(1).f_code.co_name)
        return real_writer(payload, output_dir)

    monkeypatch.setattr(mep, "write_top_pose_atlas", recording_writer)
    assert pipeline.run() is True
    # One write for the session root and one for the canonical root, both from the atlas node.
    assert callers == [f"_dag_compute_{ATLAS_NODE}_node", f"_dag_compute_{ATLAS_NODE}_node"]


def test_report_nodes_read_the_atlas_without_rewriting_it(tmp_path, monkeypatch):
    pipeline = _dag_pipeline(tmp_path, "owner_readers")
    assert pipeline.run() is True
    before = {path: _digest(path) for path in _atlas_files(pipeline)}
    mtimes = {path: path.stat().st_mtime_ns for path in _atlas_files(pipeline)}

    def refuse(*_args, **_kwargs):
        raise AssertionError("a report node must not call write_top_pose_atlas")

    monkeypatch.setattr(mep, "write_top_pose_atlas", refuse)
    scores = pipeline._load_or_build_scores()
    pipeline._write_comparative_reports(scores, analysis_scope=pipeline.analysis_scope, preserve_dag_atlas=True)
    assert {path: _digest(path) for path in _atlas_files(pipeline)} == before
    assert {path: path.stat().st_mtime_ns for path in _atlas_files(pipeline)} == mtimes


def test_preserve_mode_refuses_a_missing_atlas(tmp_path):
    pipeline = _dag_pipeline(tmp_path, "owner_missing")
    scores = pipeline._load_or_build_scores()
    # Nothing has produced the atlas yet, so a reader must fail loudly rather than pass.
    with pytest.raises(FileNotFoundError, match="DAG atlas node"):
        pipeline._write_comparative_reports(scores, analysis_scope=pipeline.analysis_scope, preserve_dag_atlas=True)
