#!/usr/bin/env python3
"""DockForge smoke checks for workflow scaffolding and preparation contracts."""

from __future__ import annotations

import argparse
import builtins
import contextlib
import io
import inspect
import hashlib
import json
import logging
import os
import shutil
import sys
import tempfile
import time
import types
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from docking.models import (
    PairlistRow,
    resolve_effective_ligand_preparation_profile,
    validate_ligand_preparation_profile,
)
from docking.execution_environment import (
    ExecutionEnvironmentConfig,
    apply_environment_runtime,
    validate_environment_selection,
)
from docking.preparation.ligand_preparation import validate_ligand_preparation_output_contract
from docking.parameter_schema import (
    apply_schema_to_runtime,
    resolve_parameter_schema,
    validate_parameter_schema,
)
from docking.engine_registry import build_runner
from docking.preflight import run_docking_preflight
from docking.project_layout import (
    ensure_engine_layout,
    ensure_numbered_output_layout,
    ensure_post_docking_compat_shim,
    ensure_project_layout,
    load_manifest,
    load_pairlist,
    pairlist_path,
    post_docking_root,
    bootstrap_project_layout,
    shared_ligands_dir,
    shared_receptors_dir,
)
from post_docking_analysis.multi_engine_pipeline import MultiEngineAnalysisPipeline
from post_docking_analysis.storage_sqlite import read_dataset_registry
from post_docking_analysis.top_pose_selector import build_top_pose_atlas
from post_docking_analysis.consensus import (
    build_consensus_rankings,
    classify_hits_target_aware,
    normalize_engine_scores,
)
from post_docking_analysis.geometric_consensus import compute_geometric_consensus
from post_docking_analysis.engine_detector import detect_engines
from post_docking_analysis.correlation_analyzer import (
    MIN_CORRELATION_N,
    analyze_vina_cnn_correlation,
    compute_cross_engine_rank_correlations,
)
from post_docking_analysis.docking_parser import parse_autodock4_dlg
from post_docking_analysis.engine_hpc_adapter import detect_engine_layout
from post_docking_analysis.gnina_hpc_adapter import detect_gnina_layout
from post_docking_analysis.vina_hpc_adapter import detect_vina_layout
from post_docking_analysis.smina_hpc_adapter import detect_smina_layout
from post_docking_analysis.autodock4_hpc_adapter import detect_autodock4_layout
from post_docking_analysis.generate_scores_csv import generate_all_scores_csv
from post_docking_analysis.pose_extractor import _pdbqt_record_to_pdb_line, extract_best_poses_from_gnina
from post_docking_analysis.report_generator import generate_analysis_root_index, generate_dashboard_index
from post_docking_analysis.unified_pipeline import UnifiedPostDockingPipeline
from post_docking_analysis.visualization_suite import (
    _build_engine_agreement_summary,
    _build_plot_context,
    _ligand_label,
    _pose_diversity_diagnostic,
    _validation_joined_table,
    generate_visualization_suite,
)
from post_docking_analysis.cli import _prompt_cli_engine_fallback
from post_docking_analysis.artifact_graph import ArtifactGraph, ArtifactNode
from post_docking_analysis.simplified_cli import (
    _is_dockforge_manifest_project,
    _resolve_project_favorite_engine,
)
from post_docking_analysis.simplified_pipeline_impl import (
    SimplifiedPostDockingPipeline,
    _normalize_rmsd_scopes as normalize_simplified_rmsd_scopes,
)
from post_docking_analysis.redocking_validation import (
    _parse_pdb_like_heavy_atoms,
    run_redocking_validation,
)
from docking.cli import deploy_main
from docking.deployment import auto_execute_rerun_manifest
from docking.hpc_profiles import load_hpc_profile
from test_pipeline import run_smoke_test
from workflow.execution import (
    BackgroundTaskManager,
    run_analysis_comparative,
    run_analysis_favorite,
    run_analysis_target,
    run_autodock_prepare,
)
from workflow.ids import build_run_id, format_run_timestamp
from workflow.interactive import _PromptNavigation, _confirm, _maybe_navigation, _optional_text, _select, _text
from workflow.interactive import (
    _analysis_engine_scope_picker,
    _analysis_progress_file,
    _refresh_latest_analysis_shortcuts,
    _run_hpc_stage,
    _write_analysis_progress,
)
import workflow.interactive as interactive_module
from workflow.models import WorkflowStepResult
from workflow.state import (
    ensure_state,
    get_feature_flags,
    list_background_tasks,
    list_checkpoint_metadata,
    load_state,
    record_checkpoint_metadata,
    summarize_state,
    update_feature_flags,
    write_meta_config_freeze,
    write_meta_env_lock,
    write_meta_run_manifest,
)


def _assert(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


class _FakePrompt:
    def __init__(self, value):
        self._value = value

    def ask(self):
        return self._value


class _FakeQuestionary:
    class Choice:
        def __init__(self, title, value, checked=False):
            self.title = title
            self.value = value
            self.checked = checked

    def __init__(self, *, select_values=None, text_values=None, confirm_values=None, checkbox_values=None):
        self._select_values = list(select_values or [])
        self._text_values = list(text_values or [])
        self._confirm_values = list(confirm_values or [])
        self._checkbox_values = list(checkbox_values or [])
        self.last_select_choices = []

    def select(self, _message, choices, default=None):
        self.last_select_choices = list(choices)
        if self._select_values:
            return _FakePrompt(self._select_values.pop(0))
        return _FakePrompt(default)

    def text(self, _message, default=""):
        if self._text_values:
            return _FakePrompt(self._text_values.pop(0))
        return _FakePrompt(default)

    def confirm(self, _message, default=True):
        if self._confirm_values:
            return _FakePrompt(self._confirm_values.pop(0))
        return _FakePrompt(default)

    def checkbox(self, _message, choices):
        if self._checkbox_values:
            return _FakePrompt(self._checkbox_values.pop(0))
        return _FakePrompt([getattr(choice, "value", "") for choice in choices if getattr(choice, "checked", False)])


def _smoke_deterministic_run_id_utility() -> None:
    fixed = datetime(2026, 3, 31, 12, 0, 0, tzinfo=timezone.utc)
    _assert(
        format_run_timestamp(fixed) == "20260331_120000",
        "format_run_timestamp should produce deterministic UTC run tokens",
    )
    token = build_run_id("Task Runner", value=fixed, suffix="abc-123")
    _assert(
        token == "task_runner_20260331_120000_abc_123",
        f"Unexpected deterministic run id format: {token}",
    )


def _smoke_feature_flags_and_checkpoint_scaffold() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_state_flags_"))
    try:
        ensure_state(temp_root)
        flags = get_feature_flags(temp_root)
        _assert("enable_timeline" in flags, "feature flags should include enable_timeline")
        _assert("enable_checkpoint_revise" in flags, "feature flags should include enable_checkpoint_revise")

        update_feature_flags(temp_root, enable_sqlite_dual_write=True, enable_timeline=False)
        updated = get_feature_flags(temp_root)
        _assert(updated["enable_sqlite_dual_write"] is True, "feature flag update should persist")
        _assert(updated["enable_timeline"] is False, "feature flag update should accept disabled toggles")

        checkpoint = record_checkpoint_metadata(
            temp_root,
            source_project_dir="/tmp/source_project",
            target_project_dir="/tmp/target_project",
            layout_profile="canonical",
            marker_file="/tmp/target_project/CHECKPOINT_REVISE_WORKSPACE.md",
            note="smoke checkpoint",
            metadata={"origin": "smoke"},
            checkpoint_id="checkpoint_20260331_120000_smoke",
        )
        records = list_checkpoint_metadata(temp_root)
        _assert(records, "checkpoint metadata records should not be empty")
        _assert(
            str(checkpoint.get("checkpoint_id", "")) == "checkpoint_20260331_120000_smoke",
            "checkpoint metadata should keep provided checkpoint id",
        )
        _assert(
            str(records[0].get("target_project_dir", "")) == "/tmp/target_project",
            "checkpoint metadata should persist target project path",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_background_task_non_blocking() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_bg_task_"))
    try:
        ensure_state(temp_root)
        manager = BackgroundTaskManager(max_workers=1)

        def _slow_job(delay_seconds: float) -> WorkflowStepResult:
            time.sleep(delay_seconds)
            return WorkflowStepResult(target="dock.run", status="completed", root=str(temp_root))

        task_id = manager.submit(
            root=temp_root,
            label="background-smoke",
            target="dock.run",
            fn=_slow_job,
            delay_seconds=0.4,
        )
        _assert(task_id.startswith("task_"), "background task id should use deterministic task prefix")

        # While the task is running, other state operations should remain responsive.
        update_feature_flags(temp_root, enable_post_docking_first_class=False)
        summary_lines = summarize_state(temp_root)
        _assert(
            any("Active background tasks" in line for line in summary_lines),
            "state summary should report active background task counts during execution",
        )

        timeout = time.time() + 5.0
        final_status = ""
        final_result_status = ""
        while time.time() < timeout:
            rows = list_background_tasks(temp_root, include_finished=True)
            if rows:
                row = next((entry for entry in rows if str(entry.get("task_id")) == task_id), rows[0])
                final_status = str(row.get("status", ""))
                final_result_status = str(row.get("result_status", ""))
                if final_status in {"completed", "failed"}:
                    break
            time.sleep(0.05)

        _assert(final_status == "completed", f"background task should finish as completed, got: {final_status}")
        _assert(
            final_result_status == "completed",
            f"background task result status should be completed, got: {final_result_status}",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_artifact_graph_executor() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_artifact_graph_"))
    try:
        metadata_root = temp_root / "4-Working" / "metadata"
        metadata_root.mkdir(parents=True, exist_ok=True)
        source_a = temp_root / "input_a.txt"
        source_b = temp_root / "input_b.txt"
        source_a.write_text("A\n", encoding="utf-8")
        source_b.write_text("B\n", encoding="utf-8")
        parsed = temp_root / "parsed.txt"
        branch_ok = temp_root / "branch_ok.txt"
        branch_optional = temp_root / "branch_optional.txt"
        report_out = temp_root / "report.txt"

        def _write_parsed():
            parsed.write_text(source_a.read_text() + source_b.read_text(), encoding="utf-8")

        def _write_branch_ok():
            time.sleep(0.2)
            branch_ok.write_text(parsed.read_text(encoding="utf-8").upper(), encoding="utf-8")

        def _write_branch_optional():
            time.sleep(0.2)
            raise RuntimeError("synthetic optional failure")

        def _write_report():
            lines = ["ok"]
            if branch_optional.exists():
                lines.append("optional")
            report_out.write_text("\n".join(lines), encoding="utf-8")

        graph = ArtifactGraph(
            cache_file=metadata_root / "dag_cache.json",
            report_file=metadata_root / "dag_execution_report.json",
            max_workers=4,
        )
        graph.register(
            ArtifactNode(
                name="parse",
                inputs=[str(source_a), str(source_b)],
                outputs=[str(parsed)],
                compute=_write_parsed,
            )
        )
        graph.register(
            ArtifactNode(
                name="branch_ok",
                inputs=[str(parsed)],
                outputs=[str(branch_ok)],
                compute=_write_branch_ok,
            )
        )
        graph.register(
            ArtifactNode(
                name="branch_optional",
                inputs=[str(parsed)],
                outputs=[str(branch_optional)],
                compute=_write_branch_optional,
                optional=True,
            )
        )
        graph.register(
            ArtifactNode(
                name="report",
                inputs=[str(branch_ok), str(branch_optional)],
                optional_inputs=[str(branch_optional)],
                outputs=[str(report_out)],
                compute=_write_report,
            )
        )

        report = graph.request(str(report_out))
        _assert(report_out.exists(), "artifact graph should materialize requested report artifact")
        _assert(report["nodes"]["parse"]["status"] == "completed", "upstream node should complete")
        _assert(report["nodes"]["branch_ok"]["status"] == "completed", "successful parallel branch should complete")
        _assert(
            report["nodes"]["branch_optional"]["status"] == "skipped_optional",
            "optional failing branch should be marked skipped_optional",
        )
        _assert(report["nodes"]["report"]["status"] == "completed", "downstream node should continue after optional failure")
        _assert(
            isinstance(report.get("status_counts"), dict) and report["status_counts"].get("completed", 0) >= 1,
            "artifact graph report should include node status counts",
        )
        _assert(
            float(report.get("parallel_time_saved_s", 0.0) or 0.0) > 0.0,
            "parallel artifact graph tiers should report positive time saved for concurrent branches",
        )
        _assert((metadata_root / "dag_cache.json").exists(), "artifact graph should persist dag cache metadata")
        _assert((metadata_root / "dag_execution_report.json").exists(), "artifact graph should persist execution report")
        persisted_report = json.loads((metadata_root / "dag_execution_report.json").read_text(encoding="utf-8"))
        _assert(
            str(Path(str(persisted_report.get("report_file", ""))).resolve()) == str((metadata_root / "dag_execution_report.json").resolve()),
            "persisted DAG report should record its own file path",
        )
        _assert(
            str(Path(str(persisted_report.get("cache_file", ""))).resolve()) == str((metadata_root / "dag_cache.json").resolve()),
            "persisted DAG report should record cache file path",
        )
        _assert(
            isinstance(persisted_report["nodes"]["branch_ok"].get("duration_s"), (int, float)),
            "persisted DAG report should include per-node duration",
        )
        _assert(
            isinstance(persisted_report["nodes"]["branch_ok"].get("cache_hit"), bool),
            "persisted DAG report should include per-node cache_hit flag",
        )
        print("✅ Artifact graph executor checks passed")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_artifact_graph_cache_hit() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_artifact_graph_cache_"))
    try:
        metadata_root = temp_root / "4-Working" / "metadata"
        metadata_root.mkdir(parents=True, exist_ok=True)
        source = temp_root / "input.txt"
        middle = temp_root / "middle.txt"
        final = temp_root / "final.txt"
        source.write_text("dockforge\n", encoding="utf-8")
        call_counts = {"prepare": 0, "finalize": 0}

        def _prepare():
            call_counts["prepare"] += 1
            middle.write_text(source.read_text(encoding="utf-8").strip().upper(), encoding="utf-8")

        def _finalize():
            call_counts["finalize"] += 1
            final.write_text(f"FINAL:{middle.read_text(encoding='utf-8')}", encoding="utf-8")

        graph = ArtifactGraph(
            cache_file=metadata_root / "dag_cache.json",
            report_file=metadata_root / "dag_execution_report.json",
            max_workers=2,
        )
        graph.register(
            ArtifactNode(
                name="prepare",
                inputs=[str(source)],
                outputs=[str(middle)],
                compute=_prepare,
            )
        )
        graph.register(
            ArtifactNode(
                name="finalize",
                inputs=[str(middle)],
                outputs=[str(final)],
                compute=_finalize,
            )
        )

        first_report = graph.request(str(final))
        second_report = graph.request(str(final))
        _assert(first_report["nodes"]["prepare"]["status"] == "completed", "first run should execute prepare node")
        _assert(second_report["nodes"]["prepare"]["status"] == "cache_hit", "unchanged rerun should hit prepare cache")
        _assert(second_report["nodes"]["finalize"]["status"] == "cache_hit", "unchanged rerun should hit finalize cache")
        _assert(call_counts["prepare"] == 1, "prepare compute should run once before cache hits")
        _assert(call_counts["finalize"] == 1, "finalize compute should run once before cache hits")
        cache_payload = json.loads((metadata_root / "dag_cache.json").read_text(encoding="utf-8"))
        _assert(
            "prepare" in (cache_payload.get("nodes") or {}) and "cache_key" in cache_payload["nodes"]["prepare"],
            "dag cache should store cache keys for completed nodes",
        )
        _assert(
            "finalize" in (cache_payload.get("nodes") or {}) and "cache_key" in cache_payload["nodes"]["finalize"],
            "dag cache should store cache keys for downstream completed nodes",
        )
        print("✅ Artifact graph cache checks passed")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_artifact_graph_registration_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_artifact_registry_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina"))
        pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(post_docking_root(temp_root) / "sessions" / "dag_registration"),
            analysis_mode="comparative_all_engines",
            engines="gnina,vina",
        )
        graph = pipeline.build_artifact_graph()
        expected_nodes = {
            "raw_scores",
            "normalized_scores",
            "validation_gate",
            "consensus_ranked",
            "classified_hits",
            "top_pose_atlas",
            "complexes",
            "prolif",
            "pandamap",
            "poseview",
            "pymol",
            "interactions",
            "polypharmacology",
            "comparative",
            "biology_correlation",
            "reports",
        }
        _assert(expected_nodes.issubset(set(graph.nodes)), "artifact graph should register the full node scaffold")
        _assert(graph.nodes["prolif"].optional is True, "prolif node should be optional in DAG scaffold")
        _assert(graph.nodes["poseview"].optional is True, "poseview node should be optional in DAG scaffold")
        _assert(pipeline.resolve_dag_scope_artifact("comparison_only") == "classified_hits", "comparison_only should map to classified_hits")
        _assert(pipeline.resolve_dag_scope_artifact("top_pose_only") == "best_poses", "top_pose_only should map to best_poses")
        print("✅ Artifact graph registration contract checks passed")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_artifact_graph_scope_request_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_artifact_scope_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina"))
        pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(post_docking_root(temp_root) / "sessions" / "dag_scope"),
            analysis_mode="comparative_all_engines",
            engines="gnina,vina",
        )
        report = pipeline.request_artifact_scope("comparison_only", force=True)
        node_names = set((report.get("nodes") or {}).keys())
        _assert(
            {"raw_scores", "normalized_scores", "validation_gate", "consensus_ranked", "classified_hits"}.issubset(node_names),
            "comparison_only scope should execute through classified_hits dependencies",
        )
        _assert("complexes" not in node_names, "comparison_only scope should not execute complexes node")
        _assert("reports" not in node_names, "comparison_only scope should not execute reports node")
        dag_paths = pipeline._dag_artifact_paths()
        _assert(dag_paths["classified_hits"].exists(), "comparison_only scope should materialize classified_hits artifact")
        _assert(not dag_paths["complexes"].exists(), "comparison_only scope should not materialize complexes output")
        _assert(not dag_paths["reports"].exists(), "comparison_only scope should not materialize reports output")

        top_pose_output = post_docking_root(temp_root) / "sessions" / "dag_scope_top_pose"
        top_pose_pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(top_pose_output),
            analysis_mode="comparative_all_engines",
            engines="gnina,vina",
        )
        top_pose_report = top_pose_pipeline.request_artifact_scope("top_pose_only", force=True)
        top_pose_nodes = set((top_pose_report.get("nodes") or {}).keys())
        _assert("top_pose_atlas" in top_pose_nodes, "top_pose_only scope should include top_pose_atlas node")
        _assert("complexes" not in top_pose_nodes, "top_pose_only scope should stop before complexes node")

        interactions_output = post_docking_root(temp_root) / "sessions" / "dag_scope_interactions"
        interactions_pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(interactions_output),
            analysis_mode="comparative_all_engines",
            engines="gnina,vina",
        )
        interactions_report = interactions_pipeline.request_artifact_scope("interactions", force=True)
        interaction_nodes = set((interactions_report.get("nodes") or {}).keys())
        _assert("interactions" in interaction_nodes, "interactions scope should include interaction aggregate node")
        _assert("reports" not in interaction_nodes, "interactions scope should not execute reports node")

        structures_output = post_docking_root(temp_root) / "sessions" / "dag_scope_structures"
        structures_pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(structures_output),
            analysis_mode="comparative_all_engines",
            engines="gnina,vina",
        )
        structures_report = structures_pipeline.request_artifact_scope("structures_only", force=True)
        structure_nodes = set((structures_report.get("nodes") or {}).keys())
        _assert("complexes" in structure_nodes, "structures_only scope should include complexes node")
        _assert("reports" not in structure_nodes, "structures_only scope should not execute reports node")
        _assert(
            structures_pipeline._dag_artifact_paths()["complexes"].exists(),
            "structures_only scope should materialize complexes output",
        )

        report_only_output = post_docking_root(temp_root) / "sessions" / "dag_scope_report_only"
        report_only_pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(report_only_output),
            analysis_mode="comparative_all_engines",
            engines="gnina,vina",
        )
        report_only_report = report_only_pipeline.request_artifact_scope("report_only", force=True)
        report_only_nodes = set((report_only_report.get("nodes") or {}).keys())
        _assert("reports" in report_only_nodes, "report_only scope should include reports node")
        # complexes/prolif/pandamap are intentionally decoupled from reports (Spec-026 fix):
        # reports no longer blocks on interaction tools; they run independently via interactions scope.
        _assert("complexes" not in report_only_nodes, "report_only scope should NOT pull in complexes (decoupled)")
        _assert(
            report_only_pipeline._dag_artifact_paths()["reports"].exists(),
            "report_only scope should materialize summary report artifact",
        )

        full_output = post_docking_root(temp_root) / "sessions" / "dag_scope_full"
        full_pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(full_output),
            analysis_mode="comparative_all_engines",
            engines="gnina,vina",
        )
        full_report = full_pipeline.request_artifact_scope("full", force=True)
        full_nodes = set((full_report.get("nodes") or {}).keys())
        _assert("reports" in full_nodes, "full scope should include reports node")
        # prolif/pandamap are optional interaction nodes, not in reports dependency chain (Spec-026 fix).
        # They are only executed when scope == "interactions".
        _assert("prolif" not in full_nodes, "full scope should NOT include prolif (decoupled from reports)")
        _assert("pandamap" not in full_nodes, "full scope should NOT include pandamap (decoupled from reports)")
        _assert("comparative" in full_nodes, "full scope should include comparative node")
        print("✅ Artifact graph scope request checks passed")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_artifact_graph_pipeline_run_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_artifact_run_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina"))
        output_dir = post_docking_root(temp_root) / "sessions" / "dag_pipeline_run"
        pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(output_dir),
            analysis_mode="comparative_all_engines",
            engines="gnina,vina",
            dag_scope="comparison_only",
            dag_force=True,
        )
        _assert(pipeline.run() is True, "pipeline.run should honor dag_scope and succeed")
        _assert(
            bool(pipeline.dag_execution_report),
            "pipeline.run with dag_scope should persist last dag execution report in-memory",
        )
        _assert(
            pipeline.dag_execution_report.get("artifact_requested_key") == "classified_hits",
            "dag execution report should record requested artifact key",
        )
        _assert(
            Path(pipeline.dag_execution_report_file).exists(),
            "pipeline.run with dag_scope should write dag execution report file",
        )
        print("✅ Artifact graph pipeline-run contract checks passed")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_artifact_graph_cli_scope_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_artifact_cli_"))
    argv_before = list(sys.argv)
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina"))
        from post_docking_analysis import cli as cli_mod

        sys.argv = [
            "post_docking_analysis",
            "--project-dir",
            str(temp_root),
            "--analysis-mode",
            "comparative_all_engines",
            "--engines",
            "gnina,vina",
            "--scope",
            "comparison_only",
            "--force",
        ]
        try:
            cli_mod.main()
        except SystemExit as exc:
            _assert(int(exc.code or 0) == 0, "CLI artifact DAG scope run should exit successfully")
        report_file = ensure_numbered_output_layout(temp_root)["post_metadata"] / "dag_execution_report.json"
        _assert(report_file.exists(), "CLI --scope should emit dag execution report")
        payload = json.loads(report_file.read_text(encoding="utf-8"))
        _assert(
            payload.get("artifact_requested_key") == "classified_hits",
            "CLI --scope comparison_only should request classified_hits artifact",
        )
        print("✅ Artifact graph CLI scope contract checks passed")
    finally:
        sys.argv = argv_before
        shutil.rmtree(temp_root, ignore_errors=True)


def _assert_valid_dag_report_file(pipeline: UnifiedPostDockingPipeline) -> dict:
    report_file = Path(str(pipeline.dag_execution_report_file))
    _assert(report_file.exists(), "dag execution report file should exist")
    payload = json.loads(report_file.read_text(encoding="utf-8"))
    _assert(isinstance(payload, dict), "dag execution report should be valid JSON object")
    _assert(isinstance(payload.get("nodes"), dict), "dag execution report should include node records")
    return payload


def _stub_dag_optional_interactions(
    pipeline: UnifiedPostDockingPipeline,
    *,
    sleep_s: float = 0.0,
    fail_prolif: bool = False,
) -> None:
    def _make_stub(node_key: str, *, fail: bool = False):
        def _compute(paths):
            if sleep_s > 0:
                time.sleep(sleep_s)
            if fail:
                raise RuntimeError(f"{node_key} smoke failure")
            return pipeline._write_dag_placeholder_outputs(
                [paths[node_key]],
                node_name=node_key,
                details="smoke optional interaction stub",
            )

        return _compute

    pipeline._dag_compute_prolif_node = _make_stub("prolif", fail=fail_prolif)
    pipeline._dag_compute_pandamap_node = _make_stub("pandamap")
    pipeline._dag_compute_poseview_node = _make_stub("poseview")
    pipeline._dag_compute_pymol_node = _make_stub("pymol")


def _smoke_artifact_graph_default_runner_migration_contract() -> None:
    source = inspect.getsource(UnifiedPostDockingPipeline.run)
    _assert("request_artifact_scope" in source, "default unified runner should dispatch through DAG scope requests")
    _assert("super().run" not in source, "default unified runner should not fall back to legacy step loop")
    print("✅ Artifact graph default unified-runner migration checks passed")


def _smoke_artifact_graph_full_parallel_execution_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_artifact_full_parallel_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina"))
        output_dir = post_docking_root(temp_root) / "sessions" / "dag_default_full_parallel"
        pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(output_dir),
            analysis_mode="comparative_all_engines",
            analysis_scope="full",
            engines="gnina,vina,smina",
            rmsd_workers=4,
        )
        _stub_dag_optional_interactions(pipeline, sleep_s=0.15)
        run_succeeded = pipeline.run()
        report = pipeline.dag_execution_report
        failed_or_blocked_nodes = {
            str(node_name): {
                "status": str((node_payload or {}).get("status") or ""),
                "details": str((node_payload or {}).get("details") or ""),
            }
            for node_name, node_payload in (report.get("nodes") or {}).items()
            if str((node_payload or {}).get("status") or "") in {"failed", "blocked_by_failure"}
        }
        _assert(
            run_succeeded is True,
            "default full unified run should succeed through DAG execution; "
            f"status_counts={report.get('status_counts')!r}; "
            f"failed_or_blocked_nodes={failed_or_blocked_nodes!r}; "
            f"report_file={getattr(pipeline, 'dag_execution_report_file', '')}",
        )
        _assert(
            report.get("artifact_requested_key") == "reports",
            "default full unified run should request reports artifact",
        )
        _assert(
            float(report.get("parallel_time_saved_s", 0.0) or 0.0) > 0.0,
            "full DAG run should report positive parallel_time_saved_s",
        )
        persisted = _assert_valid_dag_report_file(pipeline)
        _assert(
            float(persisted.get("parallel_time_saved_s", 0.0) or 0.0) > 0.0,
            "persisted DAG report should record positive parallel_time_saved_s for full runs",
        )
        print("✅ Artifact graph full-run parallel execution checks passed")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_artifact_graph_unchanged_rerun_cache_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_artifact_rerun_cache_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina"))
        output_dir = post_docking_root(temp_root) / "sessions" / "dag_rerun_cache"

        first = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(output_dir),
            analysis_mode="comparative_all_engines",
            analysis_scope="report_only",
            engines="gnina,vina,smina",
            rmsd_workers=4,
        )
        _stub_dag_optional_interactions(first, sleep_s=0.01)
        _assert(first.run() is True, "initial report-only run should succeed")

        second = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(output_dir),
            analysis_mode="comparative_all_engines",
            analysis_scope="report_only",
            engines="gnina,vina,smina",
            rmsd_workers=4,
        )
        _stub_dag_optional_interactions(second, sleep_s=0.01)
        started = time.perf_counter()
        _assert(second.run() is True, "unchanged rerun should succeed")
        elapsed = time.perf_counter() - started
        report = second.dag_execution_report
        for node_name, payload in (report.get("nodes") or {}).items():
            status = str((payload or {}).get("status") or "")
            if node_name == "reports":
                _assert(status != "cache_hit", "reports node should remain non-cacheable on unchanged reruns")
            else:
                _assert(status == "cache_hit", f"unchanged rerun should hit cache for {node_name}")
        wall_clock = float(report.get("wall_clock_s", elapsed) or elapsed)
        _assert(wall_clock < 5.0, "report_only rerun with cached upstream should finish quickly")
        _assert_valid_dag_report_file(second)
        print("✅ Artifact graph unchanged-rerun cache checks passed")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_artifact_graph_parameter_invalidation_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_artifact_param_invalidation_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina"))
        output_dir = post_docking_root(temp_root) / "sessions" / "dag_param_invalidation"

        first = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(output_dir),
            analysis_mode="comparative_all_engines",
            analysis_scope="report_only",
            normalization_method="per_engine_rank",
            engines="gnina,vina,smina",
            rmsd_workers=4,
        )
        _stub_dag_optional_interactions(first)
        _assert(first.run() is True, "initial report-only run should succeed")

        second = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(output_dir),
            analysis_mode="comparative_all_engines",
            analysis_scope="report_only",
            normalization_method="per_engine_minmax",
            engines="gnina,vina,smina",
            rmsd_workers=4,
        )
        _stub_dag_optional_interactions(second)
        _assert(second.run() is True, "normalization-change rerun should succeed")
        report = second.dag_execution_report
        node_rows = report.get("nodes") or {}
        _assert(
            str((node_rows.get("raw_scores") or {}).get("status") or "") in {"completed", "completed_with_warnings"},
            "changed analysis parameters must invalidate the conservatively scoped raw-score cache",
        )
        for node_name in (
            "normalized_scores",
            "validation_gate",
            "consensus_ranked",
            "classified_hits",
            "top_pose_atlas",
            "complexes",
            "comparative",
        ):
            _assert(
                str((node_rows.get(node_name) or {}).get("status") or "") != "cache_hit",
                f"changing normalization should invalidate downstream node {node_name}",
            )
        _assert_valid_dag_report_file(second)
        print("✅ Artifact graph parameter-invalidation checks passed")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_artifact_graph_optional_failure_isolation_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_artifact_optional_failure_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina"))
        output_dir = post_docking_root(temp_root) / "sessions" / "dag_optional_failure"
        pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(output_dir),
            analysis_mode="comparative_all_engines",
            analysis_scope="interactions",
            engines="gnina,vina,smina",
            rmsd_workers=4,
        )
        _stub_dag_optional_interactions(pipeline, fail_prolif=True)
        report = pipeline.request_artifact_scope("interactions", force=True)
        node_rows = report.get("nodes") or {}
        _assert(
            str((node_rows.get("prolif") or {}).get("status") or "") == "skipped_optional",
            "failed optional prolif node should be marked skipped_optional",
        )
        for node_name in ("pandamap", "poseview", "pymol", "interactions"):
            _assert(
                str((node_rows.get(node_name) or {}).get("status") or "") in {"completed", "completed_with_warnings"},
                f"{node_name} should still complete when optional prolif fails",
            )
        _assert("reports" not in node_rows, "interactions scope should not execute reports node")
        report_file = Path(str(report.get("report_file", ""))).expanduser()
        _assert(report_file.exists(), "dag execution report file should exist for interactions scope request")
        payload = json.loads(report_file.read_text(encoding="utf-8"))
        _assert(isinstance(payload.get("nodes"), dict), "dag execution report should include node records")
        print("✅ Artifact graph optional-failure isolation checks passed")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_artifact_graph_default_comparison_only_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_artifact_comparison_only_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina"))
        output_dir = post_docking_root(temp_root) / "sessions" / "dag_default_comparison_only"
        pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(output_dir),
            analysis_mode="comparative_all_engines",
            analysis_scope="comparison_only",
            engines="gnina,vina,smina",
        )
        _assert(pipeline.run() is True, "comparison_only unified run should succeed")
        report = pipeline.dag_execution_report
        _assert(
            report.get("artifact_requested_key") == "classified_hits",
            "comparison_only unified run should request classified_hits artifact",
        )
        dag_paths = pipeline._dag_artifact_paths()
        _assert(not dag_paths["complexes"].exists(), "comparison_only should not materialize complexes")
        _assert(not dag_paths["best_poses"].exists(), "comparison_only should not materialize top-pose outputs")
        _assert_valid_dag_report_file(pipeline)
        print("✅ Artifact graph default comparison-only checks passed")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_artifact_graph_gnina_solo_schema_contract() -> None:
    multi_root = Path(tempfile.mkdtemp(prefix="dockforge_artifact_multi_schema_"))
    solo_root = Path(tempfile.mkdtemp(prefix="dockforge_artifact_gnina_schema_"))
    try:
        bootstrap_project_layout(multi_root, ["gnina", "vina"], layout_profile="docking_legacy")
        _write_test_pairlist(multi_root)
        _seed_engine_outputs(multi_root, ("gnina", "vina"))
        multi_pipeline = UnifiedPostDockingPipeline(
            project_dir=str(multi_root),
            output_dir=str(post_docking_root(multi_root) / "sessions" / "multi_schema"),
            analysis_mode="comparative_all_engines",
            analysis_scope="comparison_only",
            engines="gnina,vina",
        )
        _assert(multi_pipeline.run() is True, "multi-engine comparison run should succeed for schema comparison")
        multi_consensus = pd.read_csv(multi_pipeline._dag_artifact_paths()["consensus_ranked"])

        bootstrap_project_layout(solo_root, ["gnina"], layout_profile="docking_legacy")
        _write_test_pairlist(solo_root)
        _seed_engine_outputs(solo_root, ("gnina",))
        solo_pipeline = UnifiedPostDockingPipeline(
            project_dir=str(solo_root),
            output_dir=str(post_docking_root(solo_root) / "sessions" / "gnina_schema"),
            analysis_mode="auto",
            analysis_scope="comparison_only",
        )
        _assert(solo_pipeline.run() is True, "GNINA solo comparison run should succeed")
        solo_consensus = pd.read_csv(solo_pipeline._dag_artifact_paths()["consensus_ranked"])
        _assert("score_name_primary" in solo_consensus.columns, "GNINA solo consensus output should expose score_name_primary")
        _assert(
            str(solo_consensus.loc[0, "score_name_primary"]) == "cnn_affinity",
            "GNINA solo consensus output should prioritize cnn_affinity",
        )
        _assert(
            set(solo_consensus.columns) == set(multi_consensus.columns),
            "GNINA solo consensus output should keep the same schema as multi-engine consensus",
        )
        _assert_valid_dag_report_file(solo_pipeline)
        print("✅ Artifact graph GNINA-solo schema checks passed")
    finally:
        shutil.rmtree(multi_root, ignore_errors=True)
        shutil.rmtree(solo_root, ignore_errors=True)


def _smoke_artifact_graph_all_strategy_schema_contract() -> None:
    roots = {
        "multi": Path(tempfile.mkdtemp(prefix="dockforge_artifact_strategy_multi_")),
        "gnina": Path(tempfile.mkdtemp(prefix="dockforge_artifact_strategy_gnina_")),
        "vina": Path(tempfile.mkdtemp(prefix="dockforge_artifact_strategy_vina_")),
        "smina": Path(tempfile.mkdtemp(prefix="dockforge_artifact_strategy_smina_")),
    }
    try:
        bootstrap_project_layout(roots["multi"], ["gnina", "vina"], layout_profile="docking_legacy")
        _write_test_pairlist(roots["multi"])
        _seed_engine_outputs(roots["multi"], ("gnina", "vina"))
        multi_pipeline = UnifiedPostDockingPipeline(
            project_dir=str(roots["multi"]),
            output_dir=str(post_docking_root(roots["multi"]) / "sessions" / "strategy_multi"),
            analysis_mode="comparative_all_engines",
            analysis_scope="comparison_only",
            engines="gnina,vina",
        )
        _assert(multi_pipeline.run() is True, "multi strategy run should succeed")
        multi_df = pd.read_csv(multi_pipeline._dag_artifact_paths()["consensus_ranked"])

        expected_primary = {
            "gnina": "cnn_affinity",
            "vina": "vina_affinity",
            "smina": "vina_affinity",
        }
        strategy_frames = {"multi": multi_df}
        for engine_name in ("gnina", "vina", "smina"):
            bootstrap_project_layout(roots[engine_name], [engine_name], layout_profile="docking_legacy")
            _write_test_pairlist(roots[engine_name])
            _seed_engine_outputs(roots[engine_name], (engine_name,))
            pipeline = UnifiedPostDockingPipeline(
                project_dir=str(roots[engine_name]),
                output_dir=str(post_docking_root(roots[engine_name]) / "sessions" / f"strategy_{engine_name}"),
                analysis_mode="auto",
                analysis_scope="comparison_only",
            )
            _assert(pipeline.run() is True, f"{engine_name} solo strategy run should succeed")
            strategy_frames[engine_name] = pd.read_csv(pipeline._dag_artifact_paths()["consensus_ranked"])

        multi_columns = set(strategy_frames["multi"].columns)
        for engine_name, frame in strategy_frames.items():
            _assert(set(frame.columns) == multi_columns, f"{engine_name} strategy should keep identical consensus schema")
        _assert(str(strategy_frames["multi"].loc[0, "engine_mode"]) == "multi", "multi strategy should label engine_mode=multi")
        for engine_name in ("gnina", "vina", "smina"):
            frame = strategy_frames[engine_name]
            _assert(bool(frame.loc[0, "single_engine_mode"]) is True, f"{engine_name} strategy should mark single_engine_mode")
            _assert(
                str(frame.loc[0, "score_name_primary"]) == expected_primary[engine_name],
                f"{engine_name} strategy should expose expected primary score",
            )
        _assert(
            str(strategy_frames["vina"].loc[0, "score_name_secondary"]) == "rmsd_lb",
            "vina solo strategy should expose rmsd_lb as secondary score",
        )
        _assert(
            str(strategy_frames["smina"].loc[0, "score_name_secondary"]) == "smina_scoring_function",
            "smina solo strategy should expose scoring provenance as secondary score",
        )
        print("✅ Artifact graph all-strategy schema checks passed")
    finally:
        for root in roots.values():
            shutil.rmtree(root, ignore_errors=True)


def _smoke_artifact_graph_per_complex_warning_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_artifact_per_complex_warn_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina",))
        pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(post_docking_root(temp_root) / "sessions" / "warning_contract"),
            analysis_mode="auto",
            analysis_scope="comparison_only",
        )
        complexes_dir = pipeline._dag_artifact_paths()["complexes"]
        complexes_dir.mkdir(parents=True, exist_ok=True)
        for stem in ("complex_ok", "complex_fail"):
            (complexes_dir / f"{stem}.pdb").write_text("ATOM      1  C   UNK A   1      0.000   0.000   0.000\nEND\n", encoding="utf-8")

        result = pipeline._run_per_complex_interaction_tasks(
            complexes_dir,
            node_name="warning_smoke",
            output_root=pipeline._dag_artifact_paths()["prolif"],
            worker=lambda pdb_file: {
                "complex": pdb_file.stem,
                "success": pdb_file.stem.endswith("ok"),
                "error": "" if pdb_file.stem.endswith("ok") else "synthetic_failure",
            },
        )
        _assert(any("synthetic_failure" in warning for warning in result.get("warnings", [])), "per-complex failure should surface in warnings")
        _assert(
            (pipeline._dag_artifact_paths()["prolif"] / "warning_smoke_per_complex_results.csv").exists(),
            "per-complex interaction helper should materialize a per-complex results table",
        )
        print("✅ Artifact graph per-complex warning checks passed")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_meta_artifact_writers() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_meta_artifacts_"))
    try:
        ensure_state(temp_root)
        config_file = write_meta_config_freeze(
            temp_root,
            config={"analysis_scope": "full", "consensus_mode": "weighted_hybrid"},
            source="smoke",
            note="meta artifact smoke test",
        )
        _assert(config_file.exists(), ".meta/config.yaml should be created")

        manifest_file = write_meta_run_manifest(
            temp_root,
            input_paths=[config_file],
            run_config={"run_id": "smoke_meta_001"},
        )
        _assert(manifest_file.exists(), ".meta/run_manifest.json should be created")
        manifest_payload = json.loads(manifest_file.read_text(encoding="utf-8"))
        _assert(int(manifest_payload.get("input_count", 0)) >= 1, "run manifest should include input checksum rows")
        _assert(
            isinstance(manifest_payload.get("tool_versions", {}), dict),
            "run manifest should include tool_versions payload",
        )

        env_lock_file = write_meta_env_lock(temp_root, preferred="pip")
        _assert(env_lock_file.exists(), ".meta/env.lock should be created")
        env_text = env_lock_file.read_text(encoding="utf-8")
        _assert("# source:" in env_text, "env.lock should document snapshot source header")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_checkpoint_revise_lineage_integrity() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_checkpoint_lineage_"))
    try:
        ensure_state(temp_root)
        cp1 = record_checkpoint_metadata(
            temp_root,
            source_project_dir=str(temp_root / "source_a"),
            target_project_dir=str(temp_root / "checkpoint_a"),
            layout_profile="canonical",
            marker_file=str(temp_root / "checkpoint_a" / "CHECKPOINT_REVISE_WORKSPACE.md"),
            note="lineage smoke - first checkpoint",
            metadata={"lineage_depth": 1},
            checkpoint_id="checkpoint_20260331_120000_a",
        )
        cp2 = record_checkpoint_metadata(
            temp_root,
            source_project_dir=str(temp_root / "checkpoint_a"),
            target_project_dir=str(temp_root / "checkpoint_b"),
            layout_profile="canonical",
            marker_file=str(temp_root / "checkpoint_b" / "CHECKPOINT_REVISE_WORKSPACE.md"),
            note="lineage smoke - revise/rerun checkpoint",
            metadata={"lineage_depth": 2, "parent_checkpoint_id": cp1.get("checkpoint_id", "")},
            checkpoint_id="checkpoint_20260331_120100_b",
        )

        rows = list_checkpoint_metadata(temp_root)
        _assert(len(rows) == 2, "lineage smoke should persist both checkpoint records")
        _assert(
            str(rows[0].get("checkpoint_id", "")) == str(cp2.get("checkpoint_id", "")),
            "checkpoint listing should be newest-first",
        )
        _assert(
            str(rows[1].get("checkpoint_id", "")) == str(cp1.get("checkpoint_id", "")),
            "checkpoint listing should retain older checkpoint",
        )
        latest_state = load_state(temp_root)
        latest_id = (
            latest_state.get("checkpoint_metadata", {}).get("latest_checkpoint_id", "")
            if isinstance(latest_state.get("checkpoint_metadata", {}), dict)
            else ""
        )
        _assert(
            str(latest_id) == "checkpoint_20260331_120100_b",
            "checkpoint store should track latest lineage checkpoint id",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_manifest_completeness_and_deterministic_paths() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_manifest_det_"))
    try:
        ensure_state(temp_root)
        inputs_dir = temp_root / "inputs"
        inputs_dir.mkdir(parents=True, exist_ok=True)
        file_a = inputs_dir / "a.txt"
        file_b = inputs_dir / "b.txt"
        file_a.write_text("alpha\n", encoding="utf-8")
        file_b.write_text("beta\n", encoding="utf-8")
        missing_file = inputs_dir / "missing.txt"

        manifest_1 = write_meta_run_manifest(
            temp_root,
            input_paths=[file_b, file_a, missing_file],
            run_config={"run_id": "determinism_001", "scope": "smoke"},
        )
        payload_1 = json.loads(manifest_1.read_text(encoding="utf-8"))
        _assert(int(payload_1.get("input_count", 0)) == 2, "manifest should include exactly two existing inputs")
        _assert(int(payload_1.get("missing_input_count", 0)) == 1, "manifest should track missing inputs")
        input_paths = list(payload_1.get("input_paths", []))
        _assert(input_paths == sorted(input_paths), "manifest input paths should be sorted deterministically")
        _assert(
            all(not str(path).startswith(str(temp_root)) for path in input_paths),
            "manifest should store project-relative paths when possible",
        )
        digest_1 = str(payload_1.get("input_checksum_digest", ""))
        _assert(bool(digest_1), "manifest should include checksum digest")

        manifest_2 = write_meta_run_manifest(
            temp_root,
            input_paths=[file_a, file_b, missing_file],
            run_config={"run_id": "determinism_001", "scope": "smoke"},
        )
        payload_2 = json.loads(manifest_2.read_text(encoding="utf-8"))
        digest_2 = str(payload_2.get("input_checksum_digest", ""))
        _assert(digest_1 == digest_2, "checksum digest should be deterministic regardless of input ordering")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _mark_score_import(layout) -> None:
    """These score-only fixtures are deliberate imports, not engine executions."""
    for path in layout["poses"].iterdir():
        if path.is_file() and path.read_text(encoding="utf-8") in {"$$$$\n", "REMARK\n"}:
            path.unlink()
    scores = layout["scores"] / "normalized_scores.csv"
    scores.with_suffix(".import.json").write_text(json.dumps({
        "kind": "imported_scores", "sha256": hashlib.sha256(scores.read_bytes()).hexdigest(),
    }), encoding="utf-8")


def _build_minimal_multi_engine_project(root: Path) -> Tuple[Path, Path]:
    bootstrap_project_layout(
        root,
        engines=["gnina", "vina"],
        project_name="dockforge-sqlite-smoke",
        favorite_engine="gnina",
        layout_profile="canonical",
    )
    pair_df = pd.DataFrame(
        [
            {
                "receptor": "R1.pdbqt",
                "site_id": "site_1",
                "ligand": "L1.pdbqt",
                "center_x": 0.0,
                "center_y": 0.0,
                "center_z": 0.0,
                "size_x": 20.0,
                "size_y": 20.0,
                "size_z": 20.0,
                "protein_display_name": "R1",
                "ligand_display_name": "L1",
            }
        ]
    )
    pair_df.to_csv(pairlist_path(root), index=False)

    normalized_rows = []
    for engine, offset in [("gnina", 0.0), ("vina", 0.2)]:
        normalized_rows.append(
            {
                "engine": engine,
                "tag": "R1.pdbqt_site_1_L1.pdbqt",
                "protein": "R1.pdbqt",
                "ligand": "L1.pdbqt",
                "site_id": "site_1",
                "pose": 1,
                "affinity_kcal_mol": -8.0 + offset,
                "score_name_primary": "vina_affinity",
                "score_primary": -8.0 + offset,
                "score_name_secondary": "",
                "score_secondary": None,
                "rmsd_lb": None,
                "rmsd_ub": None,
                "pose_file": f"/tmp/{engine}_R1_L1.pdbqt",
                "log_file": f"/tmp/{engine}_R1_L1.log",
            }
        )
    normalized_df = pd.DataFrame(normalized_rows)
    for engine in ["gnina", "vina"]:
        layout = ensure_engine_layout(root, engine, layout_profile="canonical")
        for _, row in normalized_df[normalized_df["engine"] == engine].iterrows():
            tag = str(row["tag"])
            if engine == "gnina":
                (layout["poses"] / f"{tag}.sdf").write_text("$$$$\n", encoding="utf-8")
            else:
                (layout["poses"] / f"{tag}.pdbqt").write_text("REMARK\n", encoding="utf-8")
            (layout["logs"] / f"{tag}.log").write_text("REMARK\n", encoding="utf-8")
        normalized_df[normalized_df["engine"] == engine].to_csv(
            layout["scores"] / "normalized_scores.csv",
            index=False,
        )
        _mark_score_import(layout)
    output_dir = root / "analysis" / "sessions" / "sqlite_smoke"
    return root, output_dir


def _smoke_sqlite_parity_and_skip_rationale() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_sqlite_parity_"))
    try:
        project_dir, output_dir = _build_minimal_multi_engine_project(temp_root)

        # Enabled path: expect SQLite DB + OK parity rows.
        update_feature_flags(project_dir, enable_sqlite_dual_write=True)
        pipeline = MultiEngineAnalysisPipeline(
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
        _assert(pipeline.run() is True, "SQLite-enabled comparative run should complete")

        numbered_layout = ensure_numbered_output_layout(project_dir, "canonical")
        db_file = numbered_layout["post_storage"] / "results.db"
        parity_file = numbered_layout["post_storage"] / "csv_sqlite_parity.csv"
        start_here = numbered_layout["reports_root_numbered"] / "START_HERE.md"
        dashboard_index = numbered_layout["reports_root_numbered"] / "dashboard_export_index.json"
        consolidated_json = numbered_layout["reports_root_numbered"] / "consolidated_run_summary.json"

        _assert(db_file.exists(), "SQLite DB should exist when dual-write is enabled")
        _assert(parity_file.exists(), "CSV-vs-SQLite parity report should be generated")
        _assert(start_here.exists(), "Final START_HERE report index should exist")
        _assert(dashboard_index.exists(), "Dashboard JSON export index should exist")
        _assert(consolidated_json.exists(), "Consolidated summary JSON should exist")

        registry = read_dataset_registry(db_file, pipeline.run_id)
        _assert(not registry.empty, "SQLite dataset registry should contain run datasets")
        parity_df = pd.read_csv(parity_file)
        _assert(
            set(parity_df["status"].astype(str)).issubset({"ok"}),
            "Parity report should be fully OK when SQLite dual-write is enabled",
        )

        # Disabled path: expect explicit skip rationale in parity CSV.
        update_feature_flags(project_dir, enable_sqlite_dual_write=False)
        output_dir_skip = project_dir / "analysis" / "sessions" / "sqlite_smoke_disabled"
        pipeline_skip = MultiEngineAnalysisPipeline(
            project_dir=str(project_dir),
            output_dir=str(output_dir_skip),
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
        _assert(pipeline_skip.run() is True, "SQLite-disabled comparative run should still complete")
        parity_skip = pd.read_csv(numbered_layout["post_storage"] / "csv_sqlite_parity.csv")
        _assert(
            "skipped" in set(parity_skip["status"].astype(str)),
            "When SQLite is disabled, parity file should explicitly report skipped status",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_extension_stubs() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_extensions_"))
    try:
        bootstrap_project_layout(
            temp_root,
            engines=["vina"],
            project_name="dockforge-extension-smoke",
            layout_profile="canonical",
        )
        runner = build_runner("vina", temp_root, runtime={"ensemble_receptors": {"enabled": True}})
        plan = runner.build_ensemble_plan([])
        _assert(plan.get("enabled") is True, "ensemble plan should reflect enabled scaffold mode")
        _assert(plan.get("supported") is False, "vina scaffold should report unsupported ensemble execution currently")

        rerun_payload = auto_execute_rerun_manifest(
            temp_root,
            rerun_manifest_file=temp_root / "rerun_manifest.csv",
            execution_mode="plan_only",
            submit=False,
        )
        _assert(
            str(rerun_payload.get("status", "")) == "not_implemented",
            "rerun auto-execution stub should return not_implemented status",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_simplified_cli_unified_wrapper_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_simplified_cli_wrapper_"))
    try:
        bootstrap_project_layout(
            temp_root,
            engines=["gnina", "vina"],
            project_name="dockforge-simplified-cli-wrapper",
            favorite_engine="vina",
            layout_profile="canonical",
        )
        _assert(
            _is_dockforge_manifest_project(temp_root),
            "simplified_cli should detect manifest-backed DockForge projects",
        )
        _assert(
            _resolve_project_favorite_engine(temp_root) == "vina",
            "simplified_cli unified wrapper should resolve favorite engine from manifest",
        )

        non_project = temp_root / "plain_folder"
        non_project.mkdir(parents=True, exist_ok=True)
        _assert(
            not _is_dockforge_manifest_project(non_project),
            "simplified_cli should not treat plain folders as manifest projects",
        )
        _assert(
            _resolve_project_favorite_engine(non_project) is None,
            "simplified_cli should return no favorite engine for non-project folders",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_run_analysis_target_unified_delegation() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_target_delegate_"))
    try:
        bootstrap_project_layout(
            temp_root,
            engines=["gnina"],
            project_name="dockforge-target-delegate",
            favorite_engine="gnina",
            layout_profile="canonical",
        )
        pair_df = pd.DataFrame(
            [
                {
                    "receptor": "R1.pdbqt",
                    "site_id": "site_1",
                    "ligand": "L1.pdbqt",
                    "center_x": 0.0,
                    "center_y": 0.0,
                    "center_z": 0.0,
                    "size_x": 20.0,
                    "size_y": 20.0,
                    "size_z": 20.0,
                }
            ]
        )
        pair_df.to_csv(pairlist_path(temp_root), index=False)

        layout = ensure_engine_layout(temp_root, "gnina", layout_profile="canonical")
        (layout["poses"] / "R1.pdbqt_site_1_L1.pdbqt.sdf").write_text("$$$$\n", encoding="utf-8")
        (layout["logs"] / "R1.pdbqt_site_1_L1.pdbqt.log").write_text("REMARK\n", encoding="utf-8")
        pd.DataFrame(
            [
                {
                    "engine": "gnina",
                    "tag": "R1.pdbqt_site_1_L1.pdbqt",
                    "protein": "R1.pdbqt",
                    "ligand": "L1.pdbqt",
                    "site_id": "site_1",
                    "pose": 1,
                    "affinity_kcal_mol": -8.1,
                    "score_name_primary": "vina_affinity",
                    "score_primary": -8.1,
                    "score_name_secondary": "",
                    "score_secondary": None,
                    "rmsd_lb": None,
                    "rmsd_ub": None,
                    "pose_file": "/tmp/R1_site_1_L1_pose1.pdbqt",
                    "log_file": "/tmp/R1_site_1_L1.log",
                }
            ]
        ).to_csv(layout["scores"] / "normalized_scores.csv", index=False)
        _mark_score_import(layout)

        output_dir = temp_root / "analysis" / "sessions" / "delegation_smoke"
        result = run_analysis_target(
            target="analyze.stage.reports",
            project_dir=str(temp_root),
            engine="gnina",
            output_dir=str(output_dir),
            speed_profile="fast",
        )
        _assert(
            str(result.status) == "completed",
            "run_analysis_target should complete through unified delegation for canonical projects",
        )
        _assert(
            any("delegated to unified" in str(note).lower() for note in result.notes),
            "delegated stage run should emit delegation note in WorkflowStepResult",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_config_override_overlay_contract() -> None:
    from post_docking_analysis.config_manager import load_config_overrides

    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_config_overlay_contract_"))
    try:
        config_file = temp_root / "analysis_config.yaml"
        config_file.write_text(
            "\n".join(
                [
                    "quality_control:",
                    "  hit_classification:",
                    "    strong_percentile: \"0.20\"",
                    "clean_interactions:",
                    "  max_targets: 5",
                    "  run_plip: false",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        overrides = load_config_overrides(str(config_file))
        _assert(isinstance(overrides, dict), "load_config_overrides should return a dictionary")
        _assert("analysis" not in overrides, "override overlay should not include untouched default sections")
        _assert(
            float(
                overrides.get("quality_control", {})
                .get("hit_classification", {})
                .get("strong_percentile", 0.0)
            )
            == 0.20,
            "override overlay should keep normalized explicit quality-control values only",
        )
        _assert(
            int(overrides.get("clean_interactions", {}).get("max_targets", 0)) == 5,
            "override overlay should preserve explicit clean_interactions numeric override",
        )
        _assert(
            bool(overrides.get("clean_interactions", {}).get("run_plip", True)) is False,
            "override overlay should preserve explicit clean_interactions boolean override",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_run_analysis_target_delegation_config_forwarding_contract() -> None:
    import workflow.execution as workflow_execution

    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_delegate_config_contract_"))
    captured: dict[str, object] = {}
    original_runner = workflow_execution.run_analysis_favorite
    try:
        bootstrap_project_layout(
            temp_root,
            engines=["gnina"],
            project_name="dockforge-delegate-config",
            favorite_engine="gnina",
            layout_profile="canonical",
        )
        config_file = temp_root / "analysis_config.yaml"
        config_file.write_text(
            "\n".join(
                [
                    "quality_control:",
                    "  hit_classification:",
                    "    policy: target_aware",
                    "",
                ]
            ),
            encoding="utf-8",
        )

        def _fake_run_analysis_favorite(project_dir, favorite_engine, output_dir=None, **kwargs):
            captured["project_dir"] = project_dir
            captured["favorite_engine"] = favorite_engine
            captured["config_file"] = kwargs.get("config_file")
            return WorkflowStepResult(
                target="analyze.favorite_engine",
                status="completed",
                root=str(temp_root),
                inputs={"analysis_config_file": str(config_file.resolve())},
                outputs={"output_dir": str(output_dir or "")},
            )

        workflow_execution.run_analysis_favorite = _fake_run_analysis_favorite
        result = workflow_execution.run_analysis_target(
            target="analyze.stage.reports",
            project_dir=str(temp_root),
            engine="gnina",
            output_dir=str(temp_root / "analysis" / "sessions" / "delegate_config"),
            config_file=str(config_file),
        )
        _assert(
            str(result.status) == "completed",
            "stage target should still complete when delegated favorite runner succeeds",
        )
        _assert(
            str(captured.get("config_file", "")) == str(config_file),
            "delegated stage target should forward config_file to run_analysis_favorite",
        )
        _assert(
            str(result.inputs.get("analysis_config_file", "")) == str(config_file.resolve()),
            "delegated stage result should surface resolved analysis_config_file from favorite runner",
        )
    finally:
        workflow_execution.run_analysis_favorite = original_runner
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_simplified_cli_wrapper_config_forwarding_contract() -> None:
    import post_docking_analysis.simplified_cli as simplified_cli_module
    import workflow.execution as workflow_execution

    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_simplified_cli_config_forward_"))
    captured: dict[str, object] = {}
    original_runner = workflow_execution.run_analysis_favorite
    original_layout_detector = simplified_cli_module._detect_layout
    original_argv = list(sys.argv)
    try:
        bootstrap_project_layout(
            temp_root,
            engines=["gnina"],
            project_name="dockforge-simplified-cli-config-forward",
            favorite_engine="gnina",
            layout_profile="canonical",
        )
        mock_inputs = temp_root / "mock_inputs"
        (mock_inputs / "sdf").mkdir(parents=True, exist_ok=True)
        (mock_inputs / "logs").mkdir(parents=True, exist_ok=True)
        (mock_inputs / "receptors").mkdir(parents=True, exist_ok=True)
        config_file = temp_root / "analysis_config.yaml"
        config_file.write_text("clean_interactions:\n  max_targets: 2\n", encoding="utf-8")

        def _fake_layout(_project_dir):
            return {
                "layout": "mock",
                "sdf_folder": mock_inputs / "sdf",
                "log_folder": mock_inputs / "logs",
                "receptors_folder": mock_inputs / "receptors",
                "pairlist_file": None,
            }

        def _fake_run_analysis_favorite(project_dir, favorite_engine, output_dir=None, **kwargs):
            captured["project_dir"] = project_dir
            captured["favorite_engine"] = favorite_engine
            captured["config_file"] = kwargs.get("config_file")
            return WorkflowStepResult(
                target="analyze.favorite_engine",
                status="completed",
                root=str(temp_root),
                outputs={"output_dir": str(output_dir or "")},
            )

        simplified_cli_module._detect_layout = _fake_layout
        workflow_execution.run_analysis_favorite = _fake_run_analysis_favorite
        sys.argv = [
            "simplified_cli.py",
            "--project-dir",
            str(temp_root),
            "--output",
            str(temp_root / "analysis" / "sessions" / "simplified_cli_config_forward"),
            "--config-file",
            str(config_file),
        ]
        exit_code = simplified_cli_module.main()
        _assert(exit_code == 0, "simplified_cli wrapper should succeed when delegated runner succeeds")
        _assert(
            str(captured.get("project_dir", "")) == str(temp_root),
            "simplified_cli wrapper should forward project_dir into unified runner call",
        )
        _assert(
            str(captured.get("config_file", "")) == str(config_file),
            "simplified_cli wrapper should forward --config-file into unified runner call",
        )
    finally:
        sys.argv = original_argv
        workflow_execution.run_analysis_favorite = original_runner
        simplified_cli_module._detect_layout = original_layout_detector
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_user_facing_analysis_unified_path_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_unified_path_contract_"))
    try:
        bootstrap_project_layout(
            temp_root,
            engines=["gnina", "vina"],
            project_name="dockforge-unified-path-contract",
            favorite_engine="gnina",
            layout_profile="canonical",
        )
        import post_docking_analysis.unified_pipeline as unified_module

        captured_calls = []
        original_unified = unified_module.UnifiedPostDockingPipeline

        class _FakeUnifiedPipeline:
            def __init__(self, *args, **kwargs):
                self.analysis_mode = str(kwargs.get("analysis_mode", ""))
                self.favorite_engine = str(kwargs.get("favorite_engine", "") or "")
                self.top_pose_selection_policy = str(
                    kwargs.get("top_pose_selection_policy", "best_affinity") or "best_affinity"
                )
                self.top_pose_global_aggregation = str(
                    kwargs.get("top_pose_global_aggregation", "best_target") or "best_target"
                )
                self.top_pose_outputs = {}
                self.rerun_manifest_file = ""
                captured_calls.append(
                    {
                        "analysis_mode": self.analysis_mode,
                        "favorite_engine": self.favorite_engine,
                    }
                )

            def run(self) -> bool:
                return True

        unified_module.UnifiedPostDockingPipeline = _FakeUnifiedPipeline
        try:
            comparative_result = run_analysis_comparative(
                project_dir=str(temp_root),
                output_dir=str(temp_root / "analysis" / "sessions" / "unified_comparative"),
                analysis_scope="report_only",
            )
            favorite_result = run_analysis_favorite(
                project_dir=str(temp_root),
                favorite_engine="gnina",
                output_dir=str(temp_root / "analysis" / "sessions" / "unified_favorite"),
                analysis_scope="report_only",
            )
            stage_result = run_analysis_target(
                target="analyze.stage.reports",
                project_dir=str(temp_root),
                engine="gnina",
                output_dir=str(temp_root / "analysis" / "sessions" / "unified_stage"),
                speed_profile="fast",
            )
        finally:
            unified_module.UnifiedPostDockingPipeline = original_unified

        _assert(
            str(comparative_result.status) == "completed",
            "Comparative user-facing command should complete through unified pipeline path",
        )
        _assert(
            str(favorite_result.status) == "completed",
            "Favorite-engine user-facing command should complete through unified pipeline path",
        )
        _assert(
            str(stage_result.status) == "completed",
            "Stage-target user-facing command should complete through unified pipeline path",
        )
        _assert(
            len(captured_calls) == 3,
            f"Expected 3 unified pipeline calls, observed {len(captured_calls)}",
        )
        analysis_modes = [str(entry.get("analysis_mode", "")) for entry in captured_calls]
        _assert(
            "comparative_all_engines" in analysis_modes,
            "run_analysis_comparative should instantiate unified pipeline in comparative mode",
        )
        _assert(
            analysis_modes.count("favorite_engine_continue") == 2,
            "run_analysis_favorite and canonical stage target should both instantiate unified favorite mode",
        )
        _assert(
            any("delegated to unified" in str(note).lower() for note in stage_result.notes),
            "Canonical stage run should report unified delegation in notes",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_clean_interaction_cli_dispatch_contract() -> None:
    import workflow.cli as workflow_cli
    import workflow.execution as workflow_execution

    captured: dict[str, object] = {}
    original_runner = workflow_execution.run_analysis_target

    def _fake_run_analysis_target(target: str, **kwargs):
        captured["target"] = target
        captured["kwargs"] = kwargs
        return WorkflowStepResult(
            target=target,
            status="completed",
            root=str(Path.cwd()),
            outputs={"output_dir": str(Path.cwd() / "tmp_clean_interactions")},
        )

    workflow_execution.run_analysis_target = _fake_run_analysis_target
    try:
        exit_code = workflow_cli.main(
            [
                "analyze",
                "interactions",
                "clean",
                "--dataset-root",
                "/tmp/sertaline_dataset",
                "--execution-mode",
                "local_dry_run",
                "--max-targets",
                "2",
                "--max-ligands",
                "7",
                "--max-poses",
                "4",
                "--target-filter",
                "COX2",
                "--ligand-filter",
                "Sertraline",
                "--skip-chem-fixes",
                "--skip-layered-plip",
                "--no-prolif",
                "--dry-run",
            ]
        )
    finally:
        workflow_execution.run_analysis_target = original_runner

    _assert(exit_code == 0, "Clean interaction CLI command should return success when runner succeeds")
    _assert(
        str(captured.get("target", "")) == "analyze.interactions.clean",
        "Clean interaction command should dispatch to analyze.interactions.clean target",
    )
    kwargs = captured.get("kwargs", {})
    _assert(isinstance(kwargs, dict), "Captured dispatch kwargs should be a dictionary")
    kwargs = kwargs if isinstance(kwargs, dict) else {}
    _assert(
        str(kwargs.get("dataset_root", "")) == "/tmp/sertaline_dataset",
        "Dispatch should forward --dataset-root",
    )
    _assert(
        str(kwargs.get("pipeline_execution_mode", "")) == "local_dry_run",
        "Dispatch should forward execution mode override",
    )
    _assert(
        int(kwargs.get("max_targets", 0)) == 2 and int(kwargs.get("max_ligands", 0)) == 7 and int(kwargs.get("max_poses", 0)) == 4,
        "Dispatch should forward target/ligand/pose caps",
    )
    _assert(
        bool(kwargs.get("run_chemistry_fixes", True)) is False and bool(kwargs.get("run_layered_plip", True)) is False,
        "Dispatch should map skip flags to execution booleans",
    )
    _assert(
        bool(kwargs.get("run_plip", False)) is True and bool(kwargs.get("run_prolif", True)) is False,
        "Dispatch should preserve PLIP while disabling ProLIF when requested",
    )
    _assert(bool(kwargs.get("dry_run", False)) is True, "Dispatch should forward dry-run flag")


def _smoke_clean_interaction_cli_shortcut_contract() -> None:
    import workflow.cli as workflow_cli
    import workflow.execution as workflow_execution

    captured: dict[str, object] = {}
    original_runner = workflow_execution.run_analysis_target

    def _fake_run_analysis_target(target: str, **kwargs):
        captured["target"] = target
        captured["kwargs"] = kwargs
        return WorkflowStepResult(
            target=target,
            status="completed",
            root=str(Path.cwd()),
            outputs={"output_dir": str(Path.cwd() / "tmp_clean_interactions")},
        )

    workflow_execution.run_analysis_target = _fake_run_analysis_target
    try:
        exit_code = workflow_cli.main(
            [
                "analyze",
                "clean",
                "--dry-run",
            ]
        )
    finally:
        workflow_execution.run_analysis_target = original_runner

    _assert(exit_code == 0, "Clean interaction shortcut should return success when runner succeeds")
    _assert(
        str(captured.get("target", "")) == "analyze.interactions.clean",
        "analyze clean shortcut should dispatch to analyze.interactions.clean target",
    )
    kwargs = captured.get("kwargs", {})
    _assert(isinstance(kwargs, dict), "Captured dispatch kwargs should be a dictionary")
    kwargs = kwargs if isinstance(kwargs, dict) else {}
    _assert(
        int(kwargs.get("max_targets", 0)) == 3 and int(kwargs.get("max_poses", 0)) == 3,
        "Shortcut should preserve safe clean-pipeline defaults (max_targets=3, max_poses=3)",
    )
    _assert(bool(kwargs.get("run_chemistry_fixes", False)) is True, "Shortcut should keep chemistry-fix stage enabled by default")
    _assert(bool(kwargs.get("run_plip", False)) is True and bool(kwargs.get("run_prolif", False)) is True, "Shortcut should keep PLIP and ProLIF enabled by default")
    _assert(bool(kwargs.get("dry_run", False)) is True, "Shortcut should forward dry-run flag")


def _smoke_non_gnina_rmsd_enforcement_contract() -> None:
    import docking.project_layout as project_layout
    import post_docking_analysis.protein_naming as protein_naming_module
    import post_docking_analysis.simplified_input_handler as simplified_input_module
    import post_docking_analysis.simplified_pipeline as simplified_pipeline_module
    import post_docking_analysis.unified_pipeline as unified_module
    import workflow.execution as workflow_execution

    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_non_gnina_rmsd_contract_"))
    captured_kwargs: dict[str, object] = {}

    original_ensure_engine_layout = project_layout.ensure_engine_layout
    original_pairlist_path = project_layout.pairlist_path
    original_shared_receptors_dir = project_layout.shared_receptors_dir
    original_unified_pipeline = unified_module.UnifiedPostDockingPipeline
    original_compat_cols = list(unified_module.UNIFIED_COMPAT_COLUMNS)
    original_simplified_pipeline = simplified_pipeline_module.SimplifiedPostDockingPipeline
    original_build_protein_map = protein_naming_module.build_protein_name_mapping
    original_find_receptors = simplified_input_module.find_receptor_files

    class _FakeUnifiedPipeline:
        def __init__(self, *args, **kwargs):
            self.pairlist_df = pd.DataFrame([{"receptor": "R1.pdbqt", "ligand": "L1.pdbqt", "site_id": "site_1"}])

        def _load_or_build_scores(self) -> pd.DataFrame:
            return pd.DataFrame(
                [
                    {
                        "engine": "vina",
                        "tag": "R1_site_1_L1",
                        "protein": "R1.pdbqt",
                        "ligand": "L1.pdbqt",
                        "site_id": "site_1",
                        "pose": 1,
                        "affinity_kcal_mol": -7.4,
                    }
                ]
            )

        def _build_downstream_results(self, engine_scores: pd.DataFrame) -> dict[str, pd.DataFrame]:
            best = engine_scores.copy()
            return {"best_poses": best, "full_data": engine_scores.copy()}

        def _build_simplified_bridge_scores(self, engine_scores: pd.DataFrame) -> pd.DataFrame:
            return engine_scores.copy()

        def _build_simplified_bridge_complexes(self, best_poses: pd.DataFrame) -> list[dict[str, object]]:
            return [{"tag": str(best_poses.iloc[0]["tag"])}] if not best_poses.empty else []

        def _extract_best_pose_complexes(self, manifest_df: pd.DataFrame, complexes_dir: Path) -> int:
            complexes_dir.mkdir(parents=True, exist_ok=True)
            return int(len(manifest_df))

        def _mirror_complexes_to_best_poses(self, complexes_dir: Path, best_poses_dir: Path) -> int:
            best_poses_dir.mkdir(parents=True, exist_ok=True)
            return 1

    class _FakeSimplifiedPipeline:
        def __init__(self, *args, **kwargs):
            captured_kwargs.update(kwargs)
            self.run_rmsd = bool(kwargs.get("run_rmsd", False))
            self.protein_name_override_file = temp_root / "protein_name_overrides.csv"
            self.pairlist_df = pd.DataFrame()
            self.receptor_files = []
            self.protein_name_map = {}
            self.ligand_name_map = {}
            self.scores_df = pd.DataFrame()
            self.complexes = []

        def _prompt_for_protein_name_overrides(self) -> None:
            return None

        def _persist_protein_name_mapping(self) -> None:
            return None

        def _build_ligand_name_map(self) -> dict[str, str]:
            return {}

        def _prompt_for_ligand_name_overrides(self) -> None:
            return None

        def _persist_ligand_name_mapping(self) -> None:
            return None

    def _fake_ensure_engine_layout(project_dir: str, engine: str, *args, **kwargs) -> dict[str, Path]:
        poses_dir = temp_root / "poses"
        logs_dir = temp_root / "logs"
        poses_dir.mkdir(parents=True, exist_ok=True)
        logs_dir.mkdir(parents=True, exist_ok=True)
        return {"poses": poses_dir, "logs": logs_dir}

    def _fake_pairlist_path(project_dir: str, *args, **kwargs) -> Path:
        pairlist_file = temp_root / "pairlist.csv"
        if not pairlist_file.exists():
            pd.DataFrame(
                [{"receptor": "R1.pdbqt", "site_id": "site_1", "ligand": "L1.pdbqt"}]
            ).to_csv(pairlist_file, index=False)
        return pairlist_file

    def _fake_shared_receptors_dir(project_dir: str, *args, **kwargs) -> Path:
        receptors_dir = temp_root / "receptors"
        receptors_dir.mkdir(parents=True, exist_ok=True)
        return receptors_dir

    project_layout.ensure_engine_layout = _fake_ensure_engine_layout
    project_layout.pairlist_path = _fake_pairlist_path
    project_layout.shared_receptors_dir = _fake_shared_receptors_dir
    unified_module.UnifiedPostDockingPipeline = _FakeUnifiedPipeline
    unified_module.UNIFIED_COMPAT_COLUMNS = ["engine", "tag", "protein", "ligand", "site_id", "pose", "affinity_kcal_mol"]
    simplified_pipeline_module.SimplifiedPostDockingPipeline = _FakeSimplifiedPipeline
    protein_naming_module.build_protein_name_mapping = lambda *_args, **_kwargs: {}
    simplified_input_module.find_receptor_files = lambda *_args, **_kwargs: []

    try:
        output_dir = temp_root / "analysis" / "sessions" / "non_gnina_contract"
        result = workflow_execution._prepare_non_gnina_context(
            project_dir=str(temp_root),
            engine="vina",
            output_dir=output_dir,
            ligplus_root=None,
            enable_poseview=False,
            prompt_protein_names=False,
            prompt_ligand_names=False,
            complex_query=None,
            rmsd_scopes="per_complex",
            analysis_config={},
            resume_rmsd=True,
            force_global_rmsd=False,
            global_rmsd_defer_threshold=150,
            shared_mapping_dir=None,
            exclude_problematic_ligands=False,
            positive_affinity_threshold=0.0,
            minimum_pose_count=1,
            rmsd_workers=0,
            speed_profile="standard",
        )
        _assert(isinstance(result, tuple) and len(result) == 5, "non-GNINA bridge context should return expected tuple contract")
        _assert(bool(captured_kwargs.get("run_rmsd", False)) is True, "non-GNINA bridge must force run_rmsd=True")
        _assert(bool(captured_kwargs.get("run_visualizations", False)) is True, "non-GNINA bridge should keep visualization stages enabled")
    finally:
        project_layout.ensure_engine_layout = original_ensure_engine_layout
        project_layout.pairlist_path = original_pairlist_path
        project_layout.shared_receptors_dir = original_shared_receptors_dir
        unified_module.UnifiedPostDockingPipeline = original_unified_pipeline
        unified_module.UNIFIED_COMPAT_COLUMNS = original_compat_cols
        simplified_pipeline_module.SimplifiedPostDockingPipeline = original_simplified_pipeline
        protein_naming_module.build_protein_name_mapping = original_build_protein_map
        simplified_input_module.find_receptor_files = original_find_receptors
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_legacy_stage_routes_blocked_contract() -> None:
    import workflow.execution as workflow_execution

    for target in ["analyze.stage.structure_quality", "analyze.visuals.pymol"]:
        try:
            workflow_execution.run_analysis_target(target=target, project_dir="/tmp/legacy_route_block_test")
            raise RuntimeError(f"Expected run_analysis_target to block retired target: {target}")
        except ValueError as exc:
            message = str(exc)
            _assert("retired" in message, f"{target} block message should include retirement notice")
            _assert(
                "PostDockingAnalysisPipeline" in message,
                f"{target} block message should explain legacy PostDockingAnalysisPipeline dependency",
            )

    try:
        workflow_execution._run_gnina_analysis_target("analyze.stage.structure_quality")
        raise RuntimeError("_run_gnina_analysis_target should block structure_quality before legacy route execution")
    except ValueError as exc:
        _assert("PostDockingAnalysisPipeline" in str(exc), "GNINA helper block message should reference legacy pipeline dependency")

    try:
        workflow_execution._run_non_gnina_analysis_target(
            target="analyze.visuals.pymol",
            project_dir="/tmp/legacy_route_block_test",
            engine="vina",
        )
        raise RuntimeError("_run_non_gnina_analysis_target should block pymol before legacy route execution")
    except ValueError as exc:
        _assert("PostDockingAnalysisPipeline" in str(exc), "non-GNINA helper block message should reference legacy pipeline dependency")


def _smoke_non_gnina_bridge_no_legacy_stage_execution_contract() -> None:
    source = inspect.getsource(MultiEngineAnalysisPipeline._run_non_gnina_downstream_bridge)
    _assert(
        "PostDockingAnalysisPipeline(" not in source,
        "non-GNINA downstream bridge should not instantiate legacy PostDockingAnalysisPipeline",
    )
    _assert(
        "skipped_retired_route" in source,
        "non-GNINA downstream bridge should mark retired structure-quality/PyMOL stages explicitly",
    )


def _smoke_private_analysis_helpers_reject_interaction_alias_bypass_contract() -> None:
    import workflow.execution as workflow_execution

    alias_targets = [
        "analyze.interactions.prolif",
        "analyze.interactions.ligplot",
        "analyze.interactions.pandamap",
    ]
    for target in alias_targets:
        try:
            workflow_execution._run_gnina_analysis_target(target=target, project_dir="/tmp/private_helper_alias_block")
            raise RuntimeError(f"_run_gnina_analysis_target should reject clean-contract alias bypass for {target}")
        except ValueError as exc:
            _assert(
                "normalized to analyze.interactions.clean" in str(exc),
                "GNINA helper alias bypass error should direct callers to clean-contract dispatch",
            )

        try:
            workflow_execution._run_non_gnina_analysis_target(
                target=target,
                project_dir="/tmp/private_helper_alias_block",
                engine="vina",
            )
            raise RuntimeError(f"_run_non_gnina_analysis_target should reject clean-contract alias bypass for {target}")
        except ValueError as exc:
            _assert(
                "normalized to analyze.interactions.clean" in str(exc),
                "non-GNINA helper alias bypass error should direct callers to clean-contract dispatch",
            )


def _smoke_interactive_clean_alias_collapse_contract() -> None:
    normalized, collapsed = interactive_module._collapse_interaction_alias_targets(
        [
            "analyze.interactions.prolif",
            "analyze.interactions.ligplot",
            "analyze.stage.reports",
            "analyze.interactions.clean",
        ]
    )
    _assert(
        normalized == ["analyze.interactions.clean", "analyze.stage.reports"],
        f"interactive alias collapse should dedupe to clean+reports, got {normalized}",
    )
    _assert(
        collapsed == ["analyze.interactions.prolif", "analyze.interactions.ligplot"],
        f"interactive alias collapse should track collapsed targets in order, got {collapsed}",
    )

    normalized_no_alias, collapsed_no_alias = interactive_module._collapse_interaction_alias_targets(
        ["analyze.stage.rmsd", "analyze.interactions.poseview"]
    )
    _assert(
        normalized_no_alias == ["analyze.stage.rmsd", "analyze.interactions.poseview"],
        "non-collapsed targets should remain unchanged",
    )
    _assert(collapsed_no_alias == [], "non-collapsed targets should report empty collapsed list")


def _smoke_interactive_run_tracking_summary_surface_contract() -> None:
    import post_docking_analysis.engine_detector as engine_detector_module

    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_interactive_tracking_surface_"))
    original_run_analysis_target = interactive_module.run_analysis_target
    original_detect_engines = engine_detector_module.detect_engines

    captured_manifest_path = {"path": None}

    def _fake_run_analysis_target(target: str, **kwargs):
        output_dir = Path(str(kwargs.get("output_dir", ""))).expanduser().resolve()
        run_tracking_dir = output_dir / "run_tracking"
        run_tracking_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = run_tracking_dir / "run_manifest.json"
        manifest_payload = {
            "pipeline": "interactive_smoke",
            "input": {"project_dir": str(temp_root)},
            "output_dir": str(output_dir),
            "run_started_at": datetime.now(timezone.utc).isoformat(),
            "run_completed_at": datetime.now(timezone.utc).isoformat(),
            "run_status": "completed",
            "error_message": "",
            "stage_contract": {
                "requested_run_rmsd": True,
                "requested_run_visualizations": True,
                "effective_run_rmsd": True,
                "effective_run_visualizations": True,
                "requested_run_prolif": True,
                "requested_run_ligplot": True,
                "effective_run_prolif": True,
                "effective_run_ligplot": True,
            },
            "step_status_file": str(run_tracking_dir / "step_status.csv"),
            "outputs_index_file": str(run_tracking_dir / "outputs_index.csv"),
            "steps_total": 1,
            "step_status_counts": {"completed": 1},
            "outputs_count": 0,
        }
        manifest_path.write_text(json.dumps(manifest_payload, indent=2), encoding="utf-8")
        captured_manifest_path["path"] = str(manifest_path.resolve())
        return WorkflowStepResult(
            target=target,
            status="completed",
            root=str(temp_root),
            outputs={"output_dir": str(output_dir)},
        )

    fake_q = _FakeQuestionary(
        checkbox_values=[["analyze.stage.reports"]],
        select_values=[
            "full",
            "per_engine_rank",
            "auto",
            "target_percentile",
            "best_affinity",
            "best_target",
            "standard",
            "run",
        ],
        confirm_values=[False, False, False, True],
    )

    try:
        bootstrap_project_layout(
            temp_root,
            engines=["gnina"],
            project_name="dockforge-interactive-tracking-surface",
            favorite_engine="gnina",
            layout_profile="canonical",
        )
        interactive_module.run_analysis_target = _fake_run_analysis_target
        engine_detector_module.detect_engines = lambda *_args, **_kwargs: {
            "valid_engines": ["gnina"],
            "engines": {"gnina": {"coverage_pct": 100.0}},
        }

        stdout_buffer = io.StringIO()
        with contextlib.redirect_stdout(stdout_buffer):
            interactive_module._run_analysis_stage(fake_q, temp_root)
        output = stdout_buffer.getvalue()
    finally:
        interactive_module.run_analysis_target = original_run_analysis_target
        engine_detector_module.detect_engines = original_detect_engines
        shutil.rmtree(temp_root, ignore_errors=True)

    expected_manifest_path = str(captured_manifest_path.get("path") or "")
    _assert(expected_manifest_path != "", "interactive smoke should capture a run manifest path")
    _assert("Run tracking manifest:" in output, "interactive completion should print run tracking manifest path")
    _assert(expected_manifest_path in output, "interactive completion should print absolute run manifest path")
    _assert("Stage contract:" in output, "interactive completion should print a stage contract summary line")
    _assert(
        "run_rmsd requested=True enforced=True" in output,
        "interactive completion should include one-line stage contract diff entries",
    )


def _smoke_interactive_favorite_flow_clean_followup_contract() -> None:
    import post_docking_analysis.engine_detector as engine_detector_module

    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_interactive_favorite_clean_followup_"))
    original_run_analysis_favorite = interactive_module.run_analysis_favorite
    original_run_analysis_target = interactive_module.run_analysis_target
    original_detect_engines = engine_detector_module.detect_engines

    captured_followup_calls: list[dict[str, object]] = []

    def _fake_run_analysis_favorite(project_dir, favorite_engine, output_dir=None, **_kwargs):
        return WorkflowStepResult(
            target="analyze.favorite_engine",
            status="completed",
            root=str(temp_root),
            outputs={
                "output_dir": str(output_dir or ""),
                "top_pose_per_ligand_global_file": str(temp_root / "favorite_top_pose_global.csv"),
            },
        )

    def _fake_run_analysis_target(target: str, **kwargs):
        captured_followup_calls.append({"target": target, **dict(kwargs)})
        output_dir = str(kwargs.get("output_dir", ""))
        return WorkflowStepResult(
            target=target,
            status="completed",
            root=str(temp_root),
            outputs={"output_dir": output_dir},
        )

    fake_q = _FakeQuestionary(
        checkbox_values=[["analyze.favorite_engine"]],
        select_values=[
            "full",
            "per_engine_rank",
            "auto",
            "target_percentile",
            "best_affinity",
            "best_target",
            "standard",
            "run",
            "run_clean",
        ],
        confirm_values=[False, False, False, True],
    )

    try:
        bootstrap_project_layout(
            temp_root,
            engines=["gnina"],
            project_name="dockforge-interactive-favorite-clean-followup",
            favorite_engine="gnina",
            layout_profile="canonical",
        )
        interactive_module.run_analysis_favorite = _fake_run_analysis_favorite
        interactive_module.run_analysis_target = _fake_run_analysis_target
        engine_detector_module.detect_engines = lambda *_args, **_kwargs: {
            "valid_engines": ["gnina"],
            "engines": {"gnina": {"coverage_pct": 100.0}},
        }

        stdout_buffer = io.StringIO()
        with contextlib.redirect_stdout(stdout_buffer):
            interactive_module._run_analysis_stage(fake_q, temp_root)
        output = stdout_buffer.getvalue()
    finally:
        interactive_module.run_analysis_favorite = original_run_analysis_favorite
        interactive_module.run_analysis_target = original_run_analysis_target
        engine_detector_module.detect_engines = original_detect_engines
        shutil.rmtree(temp_root, ignore_errors=True)

    _assert(
        len(captured_followup_calls) == 1,
        f"favorite-engine continuation should trigger one clean-followup call when requested, got {len(captured_followup_calls)}",
    )
    followup_call = captured_followup_calls[0]
    _assert(
        str(followup_call.get("target", "")) == "analyze.interactions.clean",
        "favorite-flow clean follow-up should dispatch through analyze.interactions.clean",
    )
    _assert(
        str(followup_call.get("project_dir", "")) == str(temp_root),
        "favorite-flow clean follow-up should forward project_dir",
    )
    clean_output = str(followup_call.get("output_dir", ""))
    _assert(
        Path(clean_output).parts[-2:] == ("favorite_engine", "clean_interactions"),
        f"favorite-flow clean follow-up should use favorite_engine/clean_interactions output suffix, got: {clean_output}",
    )
    _assert(
        "Favorite clean-interactions follow-up: completed" in output,
        "interactive completion should report favorite-flow clean follow-up status",
    )


def _smoke_noninteractive_interaction_alias_to_clean_contract() -> None:
    import workflow.execution as workflow_execution

    captured_calls: list[dict[str, object]] = []
    original_clean_runner = workflow_execution._run_clean_interaction_pipeline_target

    def _fake_clean_runner(**kwargs):
        captured_calls.append(dict(kwargs))
        return WorkflowStepResult(
            target="analyze.interactions.clean",
            status="completed",
            root=str(Path.cwd()),
            outputs={"output_dir": str(Path.cwd() / "tmp_clean_interactions")},
        )

    workflow_execution._run_clean_interaction_pipeline_target = _fake_clean_runner
    try:
        for alias_target in [
            "analyze.interactions.prolif",
            "analyze.interactions.ligplot",
            "analyze.interactions.pandamap",
        ]:
            result = workflow_execution.run_analysis_target(
                alias_target,
                project_dir="/tmp/noninteractive_alias_project",
            )
            _assert(
                str(result.target) == "analyze.interactions.clean",
                f"{alias_target} should dispatch via clean interaction target",
            )
            _assert(
                str(result.status) == "completed",
                f"{alias_target} alias dispatch should preserve clean runner status",
            )
    finally:
        workflow_execution._run_clean_interaction_pipeline_target = original_clean_runner

    _assert(len(captured_calls) == 3, f"expected 3 clean-runner alias calls, got {len(captured_calls)}")
    for index, payload in enumerate(captured_calls, start=1):
        _assert(
            str(payload.get("project_dir", "")) == "/tmp/noninteractive_alias_project",
            f"alias call {index} should forward project_dir",
        )
        _assert(
            str(payload.get("dataset_root", "")) == "/tmp/noninteractive_alias_project",
            f"alias call {index} should default dataset_root to project_dir when omitted",
        )
        _assert(
            bool(payload.get("run_plip", False)) is True and bool(payload.get("run_prolif", False)) is True,
            f"alias call {index} should preserve mandatory PLIP+ProLIF defaults",
        )


def _smoke_noninteractive_interaction_alias_cli_contract() -> None:
    import workflow.cli as workflow_cli
    import workflow.execution as workflow_execution

    captured_calls: list[dict[str, object]] = []
    original_clean_runner = workflow_execution._run_clean_interaction_pipeline_target

    def _fake_clean_runner(**kwargs):
        captured_calls.append(dict(kwargs))
        return WorkflowStepResult(
            target="analyze.interactions.clean",
            status="completed",
            root=str(Path.cwd()),
            outputs={"output_dir": str(Path.cwd() / "tmp_clean_interactions")},
        )

    workflow_execution._run_clean_interaction_pipeline_target = _fake_clean_runner
    try:
        for alias_command in ["prolif", "ligplot", "pandamap"]:
            exit_code = workflow_cli.main(
                [
                    "analyze",
                    "interactions",
                    alias_command,
                    "--project-dir",
                    "/tmp/noninteractive_alias_cli_project",
                ]
            )
            _assert(
                exit_code == 0,
                f"CLI analyze interactions {alias_command} should succeed when clean runner succeeds",
            )
    finally:
        workflow_execution._run_clean_interaction_pipeline_target = original_clean_runner

    _assert(len(captured_calls) == 3, f"expected 3 clean-runner alias calls via CLI, got {len(captured_calls)}")
    for index, payload in enumerate(captured_calls, start=1):
        _assert(
            str(payload.get("project_dir", "")) == "/tmp/noninteractive_alias_cli_project",
            f"CLI alias call {index} should forward project_dir",
        )
        _assert(
            str(payload.get("dataset_root", "")) == "/tmp/noninteractive_alias_cli_project",
            f"CLI alias call {index} should default dataset_root to project_dir when omitted",
        )
        _assert(
            bool(payload.get("run_plip", False)) is True and bool(payload.get("run_prolif", False)) is True,
            f"CLI alias call {index} should preserve mandatory PLIP+ProLIF defaults",
        )


def _smoke_prompt_back_navigation_contract() -> None:
    # Prompt token guardrail.
    try:
        _maybe_navigation("back")
        raise RuntimeError("Expected _maybe_navigation('back') to raise _PromptNavigation")
    except _PromptNavigation as exc:
        _assert(exc.action == "back", "back token should map to PromptNavigation(back)")

    try:
        _maybe_navigation("exit")
        raise RuntimeError("Expected _maybe_navigation('exit') to raise _PromptNavigation")
    except _PromptNavigation as exc:
        _assert(exc.action == "exit", "exit token should map to PromptNavigation(exit)")

    # _select should inject Back choice and propagate back selection.
    fake_select = _FakeQuestionary(select_values=["opt_a"])
    selected = _select(fake_select, "Choose option", [("opt_a", "Option A")])
    _assert(selected == "opt_a", "_select should return selected option value")
    _assert(
        any(getattr(choice, "value", "") == "back" for choice in fake_select.last_select_choices),
        "_select should always append Back navigation choice by default",
    )

    fake_back = _FakeQuestionary(select_values=["back"])
    try:
        _select(fake_back, "Choose option", [("opt_a", "Option A")])
        raise RuntimeError("Expected selecting Back to raise _PromptNavigation")
    except _PromptNavigation as exc:
        _assert(exc.action == "back", "Back selection should raise PromptNavigation(back)")

    # Text/optional/confirm prompt helpers should also support universal Back behavior.
    for helper_name, helper_fn, fake in [
        ("_text", _text, _FakeQuestionary(text_values=["back"])),
        ("_optional_text", _optional_text, _FakeQuestionary(text_values=["back"])),
        ("_confirm", _confirm, _FakeQuestionary(confirm_values=[None])),
    ]:
        try:
            if helper_name == "_text":
                helper_fn(fake, "Required text", required=True)
            elif helper_name == "_optional_text":
                helper_fn(fake, "Optional text")
            else:
                helper_fn(fake, "Confirm?")
            raise RuntimeError(f"Expected {helper_name} to raise _PromptNavigation for back navigation")
        except _PromptNavigation as exc:
            _assert(exc.action == "back", f"{helper_name} should raise PromptNavigation(back)")


def _smoke_hpc_interactive_qc_bypass_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_hpc_interactive_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")

        fake = _FakeQuestionary(
            select_values=["submit", "screen"],
            text_values=["alex-bibalex", "/cluster/users/test/project", "user@host", "hpc_screen_003"],
            confirm_values=[False, False, False, True],
            checkbox_values=[["gnina", "vina", "smina"]],
        )

        captured = {"submit_argv": None, "ligand_audits": 0, "receptor_audits": 0}
        original_load_hpc_profile = interactive_module.load_hpc_profile
        original_resolve_remote_project_dir = interactive_module.resolve_remote_project_dir
        original_resolve_ssh_target = interactive_module.resolve_ssh_target
        original_run_docking_submit = interactive_module.run_docking_submit
        original_audit_project_ligands = interactive_module.audit_project_ligands
        original_audit_project_receptors = interactive_module.audit_project_receptors
        try:
            interactive_module.load_hpc_profile = lambda *_args, **_kwargs: {}
            interactive_module.resolve_remote_project_dir = lambda *_args, **_kwargs: "/cluster/users/test/project"
            interactive_module.resolve_ssh_target = lambda *_args, **_kwargs: "user@host"

            def _fake_submit(argv):
                captured["submit_argv"] = list(argv)
                return types.SimpleNamespace(status="completed")

            def _fake_ligand_audit(*_args, **_kwargs):
                captured["ligand_audits"] += 1
                return {"issue_count": 99}

            def _fake_receptor_audit(*_args, **_kwargs):
                captured["receptor_audits"] += 1
                return {"issue_count": 99}

            interactive_module.run_docking_submit = _fake_submit
            interactive_module.audit_project_ligands = _fake_ligand_audit
            interactive_module.audit_project_receptors = _fake_receptor_audit

            _run_hpc_stage(fake, temp_root)
        finally:
            interactive_module.load_hpc_profile = original_load_hpc_profile
            interactive_module.resolve_remote_project_dir = original_resolve_remote_project_dir
            interactive_module.resolve_ssh_target = original_resolve_ssh_target
            interactive_module.run_docking_submit = original_run_docking_submit
            interactive_module.audit_project_ligands = original_audit_project_ligands
            interactive_module.audit_project_receptors = original_audit_project_receptors

        _assert(captured["ligand_audits"] == 0, "interactive HPC submit should skip ligand preflight when ligand QC gate is disabled")
        _assert(captured["receptor_audits"] == 0, "interactive HPC submit should skip receptor preflight when receptor QC gate is disabled")
        _assert(isinstance(captured["submit_argv"], list), "interactive HPC submit should invoke run_docking_submit")
        _assert("--no-ligand-qc-gate" in captured["submit_argv"], "interactive HPC submit should pass --no-ligand-qc-gate through to the submit command")
        _assert("--no-receptor-qc-gate" in captured["submit_argv"], "interactive HPC submit should pass --no-receptor-qc-gate through to the submit command")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_hpc_interactive_admet_bypass_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_hpc_interactive_admet_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")

        fake = _FakeQuestionary(
            select_values=["submit", "screen"],
            text_values=["alex-bibalex", "/cluster/users/test/project", "user@host", "hpc_screen_003"],
            confirm_values=[False, True, False, False, True],
            checkbox_values=[["gnina", "vina", "smina"]],
        )

        captured = {"submit_argv": None, "ligand_audit_kwargs": None, "receptor_audits": 0}
        original_load_hpc_profile = interactive_module.load_hpc_profile
        original_resolve_remote_project_dir = interactive_module.resolve_remote_project_dir
        original_resolve_ssh_target = interactive_module.resolve_ssh_target
        original_run_docking_submit = interactive_module.run_docking_submit
        original_audit_project_ligands = interactive_module.audit_project_ligands
        original_audit_project_receptors = interactive_module.audit_project_receptors
        try:
            interactive_module.load_hpc_profile = lambda *_args, **_kwargs: {}
            interactive_module.resolve_remote_project_dir = lambda *_args, **_kwargs: "/cluster/users/test/project"
            interactive_module.resolve_ssh_target = lambda *_args, **_kwargs: "user@host"

            def _fake_submit(argv):
                captured["submit_argv"] = list(argv)
                return types.SimpleNamespace(status="completed")

            def _fake_ligand_audit(*_args, **_kwargs):
                captured["ligand_audit_kwargs"] = dict(_kwargs)
                return {"issue_count": 0}

            def _fake_receptor_audit(*_args, **_kwargs):
                captured["receptor_audits"] += 1
                return {"issue_count": 99}

            interactive_module.run_docking_submit = _fake_submit
            interactive_module.audit_project_ligands = _fake_ligand_audit
            interactive_module.audit_project_receptors = _fake_receptor_audit

            _run_hpc_stage(fake, temp_root)
        finally:
            interactive_module.load_hpc_profile = original_load_hpc_profile
            interactive_module.resolve_remote_project_dir = original_resolve_remote_project_dir
            interactive_module.resolve_ssh_target = original_resolve_ssh_target
            interactive_module.run_docking_submit = original_run_docking_submit
            interactive_module.audit_project_ligands = original_audit_project_ligands
            interactive_module.audit_project_receptors = original_audit_project_receptors

        _assert(isinstance(captured["ligand_audit_kwargs"], dict), "interactive HPC submit should still run ligand preflight when ligand QC gate is enabled")
        _assert(
            captured["ligand_audit_kwargs"].get("enable_admet_filters") is False,
            "interactive HPC submit should disable ADMET filtering inside ligand QC when requested",
        )
        _assert(captured["receptor_audits"] == 0, "interactive HPC submit should skip receptor preflight when receptor QC gate is disabled")
        _assert(isinstance(captured["submit_argv"], list), "interactive HPC submit should invoke run_docking_submit when ADMET filters are bypassed")
        _assert("--no-ligand-admet-filters" in captured["submit_argv"], "interactive HPC submit should pass --no-ligand-admet-filters through to the submit command")
        _assert("--no-ligand-qc-gate" not in captured["submit_argv"], "interactive HPC submit should keep ligand QC enabled while only bypassing ADMET filtering")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_hpc_profile_yaml_loader_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_hpc_yaml_profile_"))
    try:
        profile_dir = temp_root / ".workflow" / "hpc_profiles"
        profile_dir.mkdir(parents=True, exist_ok=True)
        profile_file = profile_dir / "nmrbox-local.yaml"
        profile_file.write_text(
            "\n".join(
                [
                    "name: nmrbox-local",
                    "scheduler: condor",
                    "remote:",
                    "  ssh_target: user@sulfur.nmrbox.org",
                    "  project_root_base: /home/user/docking",
                    "engines:",
                    "  gnina:",
                    "    runtime:",
                    "      binary: /home/user/gnina",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        by_name = load_hpc_profile(temp_root, profile_name="nmrbox-local")
        by_file = load_hpc_profile(temp_root, profile_file=str(profile_file))
        _assert(str(by_name.get("scheduler", "")) == "condor", "YAML profile-by-name should preserve scheduler")
        _assert(
            str(by_name.get("remote", {}).get("ssh_target", "")) == "user@sulfur.nmrbox.org",
            "YAML profile-by-name should preserve remote SSH target",
        )
        _assert(str(by_file.get("source", "")) == "profile-file", "YAML profile-file load should mark source=profile-file")
        _assert(
            str(by_file.get("engines", {}).get("gnina", {}).get("runtime", {}).get("binary", "")) == "/home/user/gnina",
            "YAML profile-file load should preserve nested engine runtime fields",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_interactive_ssh_systems_yaml_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_ssh_systems_yaml_"))
    try:
        systems_file = temp_root / ".workflow" / "ssh_systems.yaml"
        systems_file.parent.mkdir(parents=True, exist_ok=True)
        systems_file.write_text(
            "\n".join(
                [
                    "systems:",
                    "  nmrbox:",
                    "    kind: nmrbox",
                    "    hpc_profile: nmrbox-condor",
                    "    ssh_target: user@sulfur.nmrbox.org",
                    "    remote_project_dir: /home/user/docking",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        systems, source = interactive_module._load_optional_ssh_systems(temp_root)
        _assert(source.endswith("ssh_systems.yaml"), "interactive SSH systems loader should report source file path")
        _assert(len(systems) == 1, f"interactive SSH systems loader should return one entry, got {len(systems)}")
        entry = systems[0]
        _assert(str(entry.get("name", "")) == "nmrbox", "interactive SSH systems loader should preserve system name")
        _assert(str(entry.get("kind", "")) == "nmrbox", "interactive SSH systems loader should preserve system kind")
        _assert(
            str(entry.get("hpc_profile", "")) == "nmrbox-condor",
            "interactive SSH systems loader should preserve mapped hpc_profile",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_nmrbox_condor_slot_parser_contract() -> None:
    raw = "\n".join(
        [
            "slot1@hostA Unclaimed Idle 16 64000 83886080 1 8.9 1",
            "slot1@hostB Claimed Busy 16 64000 83886080 1 8.9 1",
            "slot1@hostC Unclaimed Idle 8 32768 41943040 0 0.0 0",
        ]
    )
    rows = interactive_module._parse_condor_slot_snapshot(raw)
    _assert(len(rows) == 3, f"expected 3 parsed condor rows, got {len(rows)}")
    _assert(rows[0]["state"] == "unclaimed" and rows[0]["activity"] == "idle", "row 1 state/activity parse mismatch")
    _assert(int(rows[0]["cpus"]) == 16 and int(rows[2]["cpus"]) == 8, "CPU parse mismatch for condor rows")
    _assert(
        abs(interactive_module._recommended_gpu_capability([7.5, 8.9, 8.6]) - 8.9) < 1e-6,
        "GPU capability policy should prefer 8.9 when available",
    )


def _smoke_nmrbox_condor_fit_probe_contract() -> None:
    import workflow.interactive as workflow_interactive

    original_run = workflow_interactive.subprocess.run

    class _FakeCompleted:
        def __init__(self, returncode: int, stdout: str, stderr: str = ""):
            self.returncode = returncode
            self.stdout = stdout
            self.stderr = stderr

    def _fake_run(*_args, **_kwargs):
        return _FakeCompleted(
            0,
            "\n".join(
                [
                    "slot1@gpuA Unclaimed Idle 16 65536 104857600 1 8.9 1",
                    "slot1@gpuB Unclaimed Idle 12 49152 83886080 1 7.5 1",
                    "slot1@cpuA Unclaimed Idle 8 32768 41943040 0 0.0 0",
                ]
            ),
        )

    workflow_interactive.subprocess.run = _fake_run
    try:
        ok, payload, message = workflow_interactive._probe_nmrbox_condor_fit("user@sulfur.nmrbox.org")
    finally:
        workflow_interactive.subprocess.run = original_run

    _assert(ok is True, f"expected successful condor fit probe, got ok={ok} message={message}")
    _assert(int(payload.get("idle_slots_total", 0)) == 3, "idle slot count should include all unclaimed+idle rows")
    _assert(int(payload.get("idle_gpu_slots_total", 0)) == 2, "idle GPU slot count should include GPU idle rows")
    _assert(int(payload.get("recommended_cpus", 0)) >= 8, "recommended CPUs should be derived from idle-slot median")
    _assert(str(payload.get("recommended_memory", "")).endswith("GB"), "recommended memory should be formatted in GB")
    _assert(str(payload.get("recommended_disk", "")).endswith("GB"), "recommended disk should be formatted in GB")
    _assert(
        str(payload.get("recommended_requirements", "")).startswith("(GPUs >= 1)"),
        "recommended requirements should include GPU predicate when idle GPU slots exist",
    )


def _smoke_analysis_session_progress_file_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_session_progress_"))
    try:
        analysis_root = ensure_project_layout(temp_root, "docking_legacy")["analysis"]
        sessions_root = analysis_root / "sessions"
        session_root = sessions_root / "interactive_contract"
        stage_root = session_root / "stage_targets"
        stage_root.mkdir(parents=True, exist_ok=True)

        # Refresh should not prune active sessions.
        _refresh_latest_analysis_shortcuts(temp_root, session_root)
        _assert(session_root.exists(), "active analysis session root should not be pruned")
        _assert(stage_root.exists(), "active stage_targets root should not be pruned")

        # Progress writer should create parent if needed.
        progress_file = _analysis_progress_file(session_root)
        if progress_file.exists():
            progress_file.unlink()
        _write_analysis_progress(
            progress_file,
            session_root,
            selected_targets=["analyze.comparative", "analyze.stage.rmsd"],
            completed_targets=[],
            running_target="analyze.comparative",
            status="in_progress",
            note="session progress smoke contract",
        )
        _assert(progress_file.exists(), "analysis progress file should be created successfully")
        payload = json.loads(progress_file.read_text(encoding="utf-8"))
        _assert(int(payload.get("total_targets", 0)) == 2, "progress payload should include selected target count")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_profile_validator() -> None:
    invalid = validate_ligand_preparation_profile("meeko_only", ["autodock4"])
    _assert(not invalid.is_valid, "Expected meeko_only + autodock4 to be rejected")
    _assert(
        any("AutoDock4 selected" in error for error in invalid.errors),
        "Expected AutoDock4 compatibility guidance in validation errors",
    )

    valid = validate_ligand_preparation_profile("engine_aware_full", ["gnina", "vina", "smina", "autodock4"])
    _assert(valid.is_valid, "Expected engine_aware_full to be valid with all engines selected")
    _assert(
        valid.effective_profile == "openbabel_meeko_autodock",
        "Expected engine_aware_full to resolve to openbabel_meeko_autodock when AutoDock4 is selected",
    )


def _smoke_top_pose_selector_determinism() -> None:
    best_by_engine = pd.DataFrame(
        [
            {
                "engine": "gnina",
                "tag": "P1_site_1_L1_A",
                "protein": "P1.pdbqt",
                "ligand": "L1.pdbqt",
                "site_id": "site_1",
                "pose": 1,
                "affinity_kcal_mol": -8.0,
                "pose_file": "/tmp/P1_site_1_L1_A.pdbqt",
            },
            {
                "engine": "vina",
                "tag": "P1_site_1_L1_B",
                "protein": "P1.pdbqt",
                "ligand": "L1.pdbqt",
                "site_id": "site_1",
                "pose": 1,
                "affinity_kcal_mol": -8.0,
                "pose_file": "/tmp/P1_site_1_L1_B.pdbqt",
            },
            {
                "engine": "smina",
                "tag": "P1_site_1_L1_C",
                "protein": "P1.pdbqt",
                "ligand": "L1.pdbqt",
                "site_id": "site_1",
                "pose": 1,
                "affinity_kcal_mol": -9.5,
                "pose_file": "/tmp/P1_site_1_L1_C.pdbqt",
            },
        ]
    )
    consensus_df = pd.DataFrame(
        [
            {"tag": "P1_site_1_L1_A", "consensus_score": 0.20, "agreement_count": 2},
            {"tag": "P1_site_1_L1_B", "consensus_score": 0.20, "agreement_count": 2},
            {"tag": "P1_site_1_L1_C", "consensus_score": 0.80, "agreement_count": 2},
        ]
    )

    hybrid_payload = build_top_pose_atlas(
        best_by_engine=best_by_engine,
        consensus_df=consensus_df,
        selection_policy="hybrid",
    )
    hybrid_row = hybrid_payload["top_pose_per_ligand_per_protein"].iloc[0]
    _assert(
        str(hybrid_row.get("tag")) == "P1_site_1_L1_C",
        "Hybrid policy should select the highest consensus score",
    )

    affinity_payload = build_top_pose_atlas(
        best_by_engine=best_by_engine,
        consensus_df=consensus_df,
        selection_policy="best_affinity",
    )
    affinity_row = affinity_payload["top_pose_per_ligand_per_protein"].iloc[0]
    _assert(
        str(affinity_row.get("tag")) == "P1_site_1_L1_A",
        "uncalibrated single-candidate engine scores must tie deterministically without raw-energy comparison",
    )


def _smoke_top_pose_consensus_direction_contracts() -> None:
    best_by_engine = pd.DataFrame(
        [
            {
                "engine": "gnina",
                "tag": "P1_site_1_L1_A",
                "protein": "P1.pdbqt",
                "ligand": "L1.pdbqt",
                "site_id": "site_1",
                "pose": 1,
                "affinity_kcal_mol": -7.0,
            },
            {
                "engine": "vina",
                "tag": "P1_site_1_L1_B",
                "protein": "P1.pdbqt",
                "ligand": "L1.pdbqt",
                "site_id": "site_1",
                "pose": 1,
                "affinity_kcal_mol": -9.0,
            },
        ]
    )
    consensus_df = pd.DataFrame(
        [
            {"tag": "P1_site_1_L1_A", "consensus_score": 0.90, "agreement_count": 2},
            {"tag": "P1_site_1_L1_B", "consensus_score": 0.10, "agreement_count": 2},
        ]
    )

    geometric_payload = build_top_pose_atlas(
        best_by_engine=best_by_engine,
        consensus_df=consensus_df,
        selection_policy="best_consensus",
        consensus_mode="dockbox_geometric",
    )
    geometric_row = geometric_payload["top_pose_per_ligand_per_protein"].iloc[0]
    _assert(
        str(geometric_row.get("tag")) == "P1_site_1_L1_A",
        "dockbox_geometric best_consensus should treat higher consensus_score as better",
    )

    hybrid_payload = build_top_pose_atlas(
        best_by_engine=best_by_engine,
        consensus_df=consensus_df,
        selection_policy="best_consensus",
        consensus_mode="weighted_hybrid",
    )
    hybrid_row = hybrid_payload["top_pose_per_ligand_per_protein"].iloc[0]
    _assert(
        str(hybrid_row.get("tag")) == "P1_site_1_L1_A",
        "weighted_hybrid best_consensus should treat higher consensus_score as better",
    )


def _smoke_consensus_direction_and_single_engine_qc() -> None:
    frame = pd.DataFrame(
        [
            {"engine": "gnina", "protein": "P1", "tag": "P1_site_1_A", "affinity_kcal_mol": -10.0},
            {"engine": "gnina", "protein": "P1", "tag": "P1_site_1_B", "affinity_kcal_mol": -8.0},
            {"engine": "gnina", "protein": "P1", "tag": "P1_site_1_C", "affinity_kcal_mol": -6.0},
        ]
    )
    normalized = normalize_engine_scores(frame, method="per_engine_minmax")
    best_norm = float(normalized.loc[normalized["tag"] == "P1_site_1_A", "normalized_affinity_score"].iloc[0])
    worst_norm = float(normalized.loc[normalized["tag"] == "P1_site_1_C", "normalized_affinity_score"].iloc[0])
    _assert(best_norm > worst_norm, "min-max normalization should score strongest binder higher than weakest binder")

    consensus_input = frame.copy()
    consensus_input["ligand"] = ["A", "B", "C"]
    consensus_input["site_id"] = ["site_1", "site_1", "site_1"]
    ranked = build_consensus_rankings(
        consensus_input,
        consensus_mode="weighted_hybrid",
        favorite_engine="gnina",
        normalization_method="per_engine_minmax",
    )
    _assert(str(ranked.iloc[0]["tag"]) == "P1_site_1_A", "consensus ranking should keep strongest binder first")
    _assert(bool(ranked["single_engine"].all()), "single_engine flag should be true for one-engine runs")
    _assert(ranked["agreement_fraction"].isna().all(), "agreement_fraction should be NaN for single-engine runs")

    classes = classify_hits_target_aware(
        ranked,
        strong_percentile=0.4,
        moderate_percentile=0.8,
    )
    _assert(
        "warn_low_engine_agreement" not in set(classes["qc_status"].astype(str)),
        "single-engine runs should not emit warn_low_engine_agreement",
    )
    singleton = classes.head(1).copy()
    singleton["protein"] = "P_single"
    singleton_classes = classify_hits_target_aware(
        singleton,
        strong_percentile=0.4,
        moderate_percentile=0.8,
    )
    _assert(
        str(singleton_classes.iloc[0]["classification_basis"]) == "insufficient_n",
        "N<3 target groups should be marked as insufficient_n",
    )
    _assert(
        str(singleton_classes.iloc[0]["docking_quality_class"]) == "Unclassified",
        "N<3 target groups should suppress Strong/Moderate/Weak class assignment",
    )

    anchored = classify_hits_target_aware(
        ranked,
        policy="reference_anchor",
        reference_baselines=pd.DataFrame(
            [
                {
                    "protein": "P1",
                    "reference_affinity": -9.0,
                    "engine": "gnina",
                    "scoring_function": "gnina",
                    "reference_tag": "P1_ref",
                    "redocking_classification": "pass",
                }
            ]
        ),
        reference_delta_kcal_mol=1.0,
    )
    anchored_by_tag = anchored.set_index("tag")
    _assert(
        str(anchored_by_tag.loc["P1_site_1_A", "docking_quality_class"]) == "Strong",
        "reference-anchor policy should mark better-than-reference affinity as Strong",
    )
    _assert(
        str(anchored_by_tag.loc["P1_site_1_B", "docking_quality_class"]) == "Moderate",
        "reference-anchor policy should mark near-reference affinity as Moderate",
    )
    _assert(
        str(anchored_by_tag.loc["P1_site_1_C", "docking_quality_class"]) == "Weak",
        "reference-anchor policy should mark worse-than-reference affinity as Weak",
    )
    _assert(
        str(anchored_by_tag.loc["P1_site_1_A", "classification_basis"]) == "reference_anchor",
        "reference-anchor policy should emit reference_anchor classification basis when validation passes",
    )


def _smoke_dag_reference_baselines_wiring_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_ref_baseline_wiring_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        output_dir = post_docking_root(temp_root) / "sessions" / "ref_baseline_wiring"
        pipeline = MultiEngineAnalysisPipeline(
            project_dir=str(temp_root),
            output_dir=str(output_dir),
            analysis_mode="comparative_all_engines",
            hit_class_policy="reference_anchor",
        )
        paths = pipeline._dag_artifact_paths()
        Path(paths["consensus_ranked"]).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            [
                {
                    "protein": "P1",
                    "tag": "P1_site_1_A",
                    "ligand": "L1",
                    "site_id": "site_1",
                    "consensus_score": 0.05,
                    "best_affinity_kcal_mol": -10.0,
                    "winner_engine": "gnina",
                    "winner_scoring_function": "gnina",
                    "agreement_count": 2,
                    "single_engine": False,
                }
            ]
        ).to_csv(paths["consensus_ranked"], index=False)
        pd.DataFrame(
            [
                {
                    "protein": "P1",
                    "reference_affinity": -9.0,
                    "engine": "gnina",
                    "scoring_function": "gnina",
                    "reference_tag": "P1_ref",
                    "redocking_classification": "pass",
                }
            ]
        ).to_csv(paths["reference_baselines"], index=False)

        paths["validation_gate"].parent.mkdir(parents=True, exist_ok=True)
        paths["validation_gate"].write_text(json.dumps({"allow_reference_anchor": True}), encoding="utf-8")
        node_result = pipeline._dag_compute_classified_hits_node(paths)
        _assert("classified_rows=1" in str(node_result.get("details", "")), "classified-hits node should process test row")
        classified = pd.read_csv(paths["classified_hits"])
        _assert(not classified.empty, "classified hits output should be present")
        _assert(
            str(classified.loc[0, "classification_basis"]) == "reference_anchor",
            "DAG classified-hits node should apply reference baselines when present",
        )
        _assert(
            bool(classified.loc[0, "reference_anchor_available"]) is True,
            "reference anchor should be marked available when validated baseline exists",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_geometric_consensus_soft_alignment_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_geo_soft_align_"))
    try:
        sdf_file = temp_root / "gnina_pose.sdf"
        pdbqt_file = temp_root / "vina_pose.pdbqt"
        coords_44 = _compact_test_coords(44, offset=0.0)
        coords_43 = [(x + 0.004, y - 0.003, z + 0.002) for x, y, z in coords_44[:43]]
        _write_test_sdf(sdf_file, [coords_44])
        _write_test_pdbqt_pose(pdbqt_file, coords_43)

        geometric = compute_geometric_consensus(
            pd.DataFrame(
                [
                    {"engine": "gnina", "tag": "P1_site_1_L1", "pose_file": str(sdf_file), "pose": 1},
                    {"engine": "vina", "tag": "P1_site_1_L1", "pose_file": str(pdbqt_file), "pose": 1},
                ]
            ),
            rmsd_cutoff=2.0,
        )
        _assert(len(geometric) == 1, "geometric consensus should produce one row for one tag")
        row = geometric.iloc[0]
        _assert(bool(row.get("geometric_agreement", False)) is False, "atom-count mismatch must not produce geometric agreement")
        pair_rows = json.loads(str(row.get("geometric_pairwise_rmsd_json", "[]")))
        _assert(bool(pair_rows), "pairwise geometric JSON should include one engine-pair row")
        _assert(bool(pair_rows[0].get("soft_alignment_used", False)) is False, "atom truncation must never repair a mismatch")
        _assert(
            pair_rows[0].get("rmsd") is None,
            "unsupported atom mapping must not have a fabricated RMSD",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_geometric_consensus_hard_mismatch_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_geo_hard_mismatch_"))
    try:
        sdf_file = temp_root / "gnina_pose.sdf"
        pdbqt_file = temp_root / "vina_pose.pdbqt"
        coords_44 = _compact_test_coords(44, offset=0.0)
        coords_30 = _compact_test_coords(30, offset=0.01)
        _write_test_sdf(sdf_file, [coords_44])
        _write_test_pdbqt_pose(pdbqt_file, coords_30)

        geometric = compute_geometric_consensus(
            pd.DataFrame(
                [
                    {"engine": "gnina", "tag": "P1_site_1_L1", "pose_file": str(sdf_file), "pose": 1},
                    {"engine": "vina", "tag": "P1_site_1_L1", "pose_file": str(pdbqt_file), "pose": 1},
                ]
            ),
            rmsd_cutoff=2.0,
        )
        _assert(len(geometric) == 1, "geometric consensus should produce one row for one tag")
        row = geometric.iloc[0]
        _assert(bool(row.get("geometric_agreement", True)) is False, "count mismatch >2 should remain a hard failure")
        pair_rows = json.loads(str(row.get("geometric_pairwise_rmsd_json", "[]")))
        _assert(bool(pair_rows), "pairwise geometric JSON should include one engine-pair row")
        _assert(
            bool(pair_rows[0].get("soft_alignment_used", False)) is False,
            "soft alignment should not run when count delta exceeds safeguard threshold",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_redocking_multi_pose_sdf_parser_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_redock_parser_"))
    try:
        sdf_file = temp_root / "gnina_multi_pose.sdf"
        coords_pose_1 = _compact_test_coords(44, offset=0.0)
        coords_pose_2 = _compact_test_coords(44, offset=0.02)
        # Deliberately malformed first counts line (88) mirrors the historical
        # overcount risk when a non-block parser spans multiple poses.
        _write_test_sdf(
            sdf_file,
            [coords_pose_1, coords_pose_2],
            declared_atom_counts=[88, 44],
        )
        coords, elements, error = _parse_pdb_like_heavy_atoms(sdf_file, pose_index=1)
        _assert(bool(error), "malformed counts must fail, never borrow atoms from another record")
        coords, elements, error = _parse_pdb_like_heavy_atoms(sdf_file, pose_index=2)
        _assert(error == "", f"the valid second record remains independently readable: {error}")
        _assert(int(coords.shape[0]) == 44, "only the selected record's atoms should be parsed")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_redocking_validation_multi_pose_end_to_end_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_redock_e2e_"))
    try:
        docked_sdf = temp_root / "gnina_best_pose.sdf"
        reference_pose = temp_root / "reference_pose.sdf"
        coords_pose_1 = _compact_test_coords(44, offset=0.0)
        coords_pose_2 = _compact_test_coords(44, offset=0.20)
        _write_test_sdf(
            docked_sdf,
            [coords_pose_1, coords_pose_2],
            declared_atom_counts=[44, 44],
        )
        _write_test_sdf(reference_pose, [coords_pose_1])

        best_by_engine = pd.DataFrame(
            [
                {
                    "protein": "P1",
                    "ligand": "REF1",
                    "tag": "P1_site_1_REF1",
                    "engine": "gnina",
                    "pose_file": str(docked_sdf),
                    "affinity_kcal_mol": -8.2,
                    "is_cocrystal_benchmark": True,
                    "pair_source": "reference",
                    "reference_pose_file": str(reference_pose),
                    "reference_source": "experimental:1ABC",
                    "reference_frame_id": "native-1ABC",
                    "receptor_frame_id": "native-1ABC",
                    "reference_pdb_id": "1ABC",
                }
            ]
        )
        outputs = run_redocking_validation(
            project_dir=temp_root,
            best_by_engine=best_by_engine,
            output_dir=temp_root / "reports",
        )
        validation_df = outputs.get("validation_df")
        _assert(isinstance(validation_df, pd.DataFrame), "redocking validation should return validation_df")
        _assert(not validation_df.empty, "redocking validation should include a row for reference candidate input")
        classification = str(validation_df.iloc[0]["redocking_classification"])
        _assert(
            classification in {"pass", "warn", "fail"},
            f"multi-pose SDF redocking should be evaluable (not not_evaluable), got {classification}",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _build_favorite_engine_guard_pipeline(
    root: Path,
    *,
    favorite_engine: str,
    engines_in_scope: List[str],
) -> MultiEngineAnalysisPipeline:
    bootstrap_project_layout(
        root,
        engines=["vina"],
        project_name="dockforge-favorite-engine-guard",
        favorite_engine="vina",
        layout_profile="canonical",
    )
    _write_test_pairlist(root)
    pipeline = MultiEngineAnalysisPipeline(
        project_dir=str(root),
        output_dir=str(root / "analysis" / "sessions" / "favorite_guard"),
        analysis_mode="favorite_engine_continue",
        analysis_scope="comparison_only",
        favorite_engine=favorite_engine,
        engines_in_scope=engines_in_scope,
        normalization_method="per_engine_rank",
    )
    seed_scores = pd.DataFrame(
        [
            {
                "engine": "vina",
                "tag": "R1_site_1_L1",
                "protein": "R1",
                "ligand": "L1",
                "affinity_kcal_mol": -7.2,
            }
        ]
    )
    pipeline._load_or_build_scores = lambda: seed_scores.copy()
    pipeline._apply_problematic_ligand_filters = lambda frame: frame
    pipeline._write_comparative_reports = lambda *_args, **_kwargs: None
    pipeline._write_run_tracking_artifacts = lambda **_kwargs: None
    return pipeline


def _smoke_favorite_engine_scope_guard_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_favorite_guard_scope_"))
    try:
        pipeline = _build_favorite_engine_guard_pipeline(
            temp_root,
            favorite_engine="gnina",
            engines_in_scope=["vina"],
        )
        try:
            pipeline.run()
            raise RuntimeError("favorite_engine scope guard should raise for out-of-scope engine")
        except ValueError as exc:
            message = str(exc).lower()
            _assert("gnina" in message, "favorite-engine scope guard message should include requested engine")
            _assert("vina" in message, "favorite-engine scope guard message should list valid scoped engines")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_favorite_engine_empty_scope_guard_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_favorite_guard_empty_"))
    try:
        pipeline = _build_favorite_engine_guard_pipeline(
            temp_root,
            favorite_engine="gnina",
            engines_in_scope=[],
        )
        pipeline._write_single_engine_reports = (
            lambda _scores, _engine, output_dir: Path(output_dir).mkdir(parents=True, exist_ok=True)
        )
        pipeline._run_downstream_bridge = lambda *_args, **_kwargs: None
        _assert(
            pipeline.run() is True,
            "empty engines_in_scope should not trigger favorite-engine membership guard",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_pairlist_directional_matching_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_pairlist_match_"))
    try:
        logs_dir = temp_root / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        pairlist_file = temp_root / "pairlist.csv"
        out_file = temp_root / "all_scores.csv"

        pd.DataFrame(
            [
                {"receptor": "3KK6", "site_id": "site_1", "ligand": "CEL"},
            ]
        ).to_csv(pairlist_file, index=False)

        log_payload = "\n".join(
            [
                "mode | affinity | intramol | CNN pose score | CNN affinity",
                "-----+------------+------------+------------+----------",
                "   1      -8.00       -0.50       0.8000      -7.50",
                "",
            ]
        )
        (logs_dir / "3KK6_site_1_CEL.log").write_text(log_payload, encoding="utf-8")
        (logs_dir / "3KK6_site_1_CEL_A_701.log").write_text(log_payload, encoding="utf-8")
        # Legacy score tables require a corresponding pose; a log alone is not
        # evidence of a completed docking result.
        for stem in ("3KK6_site_1_CEL", "3KK6_site_1_CEL_A_701"):
            _write_test_sdf(logs_dir / f"{stem}.sdf", [[(0.0, 0.0, 0.0)]])

        success = generate_all_scores_csv(
            gnina_out_dir=logs_dir,
            output_file=out_file,
            pairlist_file=pairlist_file,
            log_dir=logs_dir,
        )
        _assert(success is True, "generate_all_scores_csv should complete with strict pairlist matching")
        scores = pd.read_csv(out_file)
        _assert(
            "3KK6_site_1_CEL_A_701" not in set(scores["tag"].astype(str)),
            "unmatched log files must be reported but excluded from scored candidates",
        )

        unmatched_file = temp_root / "unmatched_log_files.txt"
        _assert(unmatched_file.exists(), "unmatched_log_files.txt should be generated for unmatched logs")
        unmatched = unmatched_file.read_text(encoding="utf-8")
        _assert(
            "3KK6_site_1_CEL_A_701.log" in unmatched,
            "unmatched report should include strict-mismatch logfile",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_pose_parser_chain_and_dashboard_contract() -> None:
    raw = "HETATM    1  C1  LIG A   5      1.000   2.000   3.000  1.00  0.00           C"
    parsed = _pdbqt_record_to_pdb_line(raw)
    _assert(parsed is not None, "pdbqt fallback parser should return a valid pdb line")
    parsed_tokens = str(parsed).split()
    _assert(len(parsed_tokens) >= 6, "parsed pdb line should retain chain and residue fields")
    _assert(parsed_tokens[4] == "A", "chain ID should be preserved as A")
    _assert(parsed_tokens[5] == "5", "residue number should remain 5")

    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_dashboard_contract_"))
    try:
        outputs = generate_dashboard_index(
            output_dir=temp_root,
            run_context={"analysis_scope": "full"},
            artifact_paths={"example": temp_root / "dummy.csv"},
        )
        contract_file = Path(str(outputs["dashboard_contract_file"]))
        content = contract_file.read_text(encoding="utf-8")
        _assert("\\n" not in content, "dashboard contract should contain actual newlines, not escaped literals")
        _assert("\n" in content, "dashboard contract should contain real newline separators")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_multiligand_residue_identity_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_multiligand_residue_"))
    try:
        input_dir = temp_root / "input"
        gnina_dir = input_dir / "gnina_out"
        receptors_dir = input_dir / "receptors"
        output_dir = temp_root / "analysis"
        gnina_dir.mkdir(parents=True, exist_ok=True)
        receptors_dir.mkdir(parents=True, exist_ok=True)
        output_dir.mkdir(parents=True, exist_ok=True)

        tags = [
            ("2XCT_cleaned_site_1_CPF", "CPF", -8.2),
            ("2XCT_cleaned_site_1_EVP", "EVP", -7.9),
        ]
        scores = pd.DataFrame(
            [
                {"tag": tag, "mode": 1, "vina_affinity": affinity, "cnn_score": 0.80}
                for tag, _resname, affinity in tags
            ]
        )
        scores.to_csv(gnina_dir / "all_scores.csv", index=False)

        sdf_template = (
            "Ligand\n"
            "  DockForge\n"
            "\n"
            "  1  0  0  0  0  0            999 V2000\n"
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "M  END\n"
            "$$$$\n"
        )
        for tag, _resname, _affinity in tags:
            (gnina_dir / f"{tag}_top.sdf").write_text(sdf_template, encoding="utf-8")

        receptor_pdbqt = (
            "ATOM      1  N   ALA A   1      11.104  13.207   2.100  1.00  0.00    -0.300 N\n"
            "ATOM      2  CA  ALA A   1      12.000  12.200   2.600  1.00  0.00     0.100 C\n"
            "END\n"
        )
        (receptors_dir / "2XCT_cleaned_prep.pdbqt").write_text(receptor_pdbqt, encoding="utf-8")

        extracted_count = extract_best_poses_from_gnina(
            input_dir=input_dir,
            output_dir=output_dir,
            gnina_dir=gnina_dir,
            receptors_dir=receptors_dir,
        )
        _assert(extracted_count == 2, f"Expected two extracted complexes, observed {extracted_count}")

        observed_ligand_resnames = set()
        for tag, expected_resname, _affinity in tags:
            complex_pdb = output_dir / "best_poses" / tag / f"{tag}_pose1.pdb"
            _assert(complex_pdb.exists(), f"Expected extracted complex file: {complex_pdb}")
            hetatm_lines = [
                line for line in complex_pdb.read_text(encoding="utf-8").splitlines()
                if line.startswith("HETATM")
            ]
            _assert(hetatm_lines, f"Expected ligand HETATM lines in {complex_pdb.name}")
            ligand_resnames = {line[17:20].strip() for line in hetatm_lines}
            observed_ligand_resnames.update(ligand_resnames)
            _assert(
                expected_resname in ligand_resnames,
                f"Ligand residue name should be preserved as {expected_resname} in {complex_pdb.name}",
            )
            _assert(
                "UNK" not in ligand_resnames,
                f"Ligand residue name should not regress to UNK in {complex_pdb.name}",
            )

        _assert(
            observed_ligand_resnames == {"CPF", "EVP"},
            "Multi-ligand extraction should preserve distinct residue identities end-to-end",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_vina_cnn_correlation_stats_contract() -> None:
    # n < 5 -> insufficient_n contract
    small_df = pd.DataFrame(
        {
            "vina_affinity": [-8.0, -7.8, -7.2, -6.9],
            "cnn_affinity": [-7.5, -7.2, -6.7, -6.5],
            "cnn_score": [0.91, 0.89, 0.81, 0.78],
        }
    )
    small_res = analyze_vina_cnn_correlation(small_df)
    small_pair = small_res["pearson_correlations"]["vina_cnn_affinity"]
    _assert(bool(small_pair.get("insufficient_n", False)), "n<5 should be flagged as insufficient_n")
    _assert(pd.isna(small_pair.get("correlation")), "n<5 should suppress correlation value")

    # n >= 5 -> p/q values present
    full_df = pd.DataFrame(
        {
            "vina_affinity": [-9.2, -8.8, -8.1, -7.9, -7.3, -6.8, -6.2],
            "cnn_affinity": [-8.6, -8.2, -7.9, -7.4, -7.0, -6.5, -6.0],
            "cnn_score": [0.96, 0.93, 0.90, 0.87, 0.83, 0.80, 0.76],
        }
    )
    full_res = analyze_vina_cnn_correlation(full_df)
    for block in ("pearson_correlations", "spearman_correlations"):
        for key in ("vina_cnn_affinity", "vina_cnn_score", "cnn_affinity_score"):
            payload = full_res[block][key]
            _assert("p_value" in payload and "q_value" in payload, "correlation payload must include p and q values")
            _assert(pd.notna(payload["p_value"]), "p_value should be numeric for n>=5")
            _assert(pd.notna(payload["q_value"]), "q_value should be numeric for n>=5")


def _smoke_cross_engine_correlation_contract_sparse_and_mixed() -> None:
    rows = []
    engines = [("gnina", 0.00), ("vina", 0.20), ("smina", 0.35)]

    # Sparse protein: only 4 paired tags (below MIN_CORRELATION_N), should be suppressed.
    for idx in range(1, 5):
        for engine, offset in engines:
            rows.append(
                {
                    "engine": engine,
                    "protein": "P_sparse",
                    "tag": f"P_sparse_L{idx}",
                    "affinity_kcal_mol": -8.5 + (0.15 * idx) + offset,
                }
            )

    # Mixed protein: 6 paired tags (>= minimum, < LOW_POWER_N), should compute with low_power flag.
    for idx in range(1, 7):
        for engine, offset in engines:
            rows.append(
                {
                    "engine": engine,
                    "protein": "P_mixed",
                    "tag": f"P_mixed_L{idx}",
                    "affinity_kcal_mol": -9.2 + (0.12 * idx) + offset,
                }
            )

    best_by_engine = pd.DataFrame(rows)
    per_protein_df, global_df = compute_cross_engine_rank_correlations(best_by_engine)

    _assert(not per_protein_df.empty, "per-protein correlation table should not be empty")
    _assert(not global_df.empty, "global correlation table should not be empty")
    _assert(
        set(str(value) for value in per_protein_df["protein"].dropna().unique()) == {"P_sparse", "P_mixed"},
        "per-protein output should include both sparse and mixed proteins",
    )
    _assert(
        set(str(value) for value in per_protein_df["normalization_method"].dropna().unique()) == {"rank_normalized"},
        "per-protein correlation rows should declare rank_normalized method",
    )
    _assert(
        set(str(value) for value in global_df["normalization_method"].dropna().unique()) == {"rank_normalized"},
        "global correlation rows should declare rank_normalized method",
    )

    sparse_rows = per_protein_df[per_protein_df["protein"] == "P_sparse"]
    _assert(not sparse_rows.empty, "sparse protein rows should be present")
    _assert(
        sparse_rows["insufficient_n"].astype(bool).all(),
        "sparse protein rows should be marked insufficient_n when paired tags are below minimum",
    )
    _assert(
        pd.to_numeric(sparse_rows["spearman_q_value"], errors="coerce").isna().all(),
        "sparse protein rows should not produce q-values",
    )
    _assert(
        pd.to_numeric(sparse_rows["pearson_q_value"], errors="coerce").isna().all(),
        "sparse protein rows should not produce q-values",
    )

    mixed_rows = per_protein_df[per_protein_df["protein"] == "P_mixed"]
    _assert(not mixed_rows.empty, "mixed protein rows should be present")
    _assert(
        (~mixed_rows["insufficient_n"].astype(bool)).all(),
        "mixed protein rows should compute correlations when paired tags meet minimum N",
    )
    _assert(
        (pd.to_numeric(mixed_rows["paired_tags"], errors="coerce") >= MIN_CORRELATION_N).all(),
        "mixed protein paired tag counts should satisfy MIN_CORRELATION_N",
    )
    _assert(
        mixed_rows["low_power"].astype(bool).all(),
        "mixed protein rows should be flagged low_power when N is between minimum and LOW_POWER_N",
    )

    _assert(
        (~global_df["insufficient_n"].astype(bool)).all(),
        "global correlation rows should compute when aggregated paired tags satisfy minimum N",
    )
    _assert(
        pd.to_numeric(global_df["spearman_q_value"], errors="coerce").notna().all(),
        "global rows should include BH-corrected spearman q-values",
    )
    _assert(
        pd.to_numeric(global_df["pearson_q_value"], errors="coerce").notna().all(),
        "global rows should include BH-corrected pearson q-values",
    )


def _smoke_consensus_mode_isolation_contract() -> None:
    import post_docking_analysis.consensus as consensus_module

    frame = pd.DataFrame(
        [
            {
                "engine": "gnina",
                "protein": "P1",
                "ligand": "L1",
                "site_id": "site_1",
                "tag": "P1_site_1_A",
                "pose": 1,
                "affinity_kcal_mol": -10.0,
                "pose_file": "/tmp/P1_A_gnina.pdbqt",
            },
            {
                "engine": "vina",
                "protein": "P1",
                "ligand": "L1",
                "site_id": "site_1",
                "tag": "P1_site_1_A",
                "pose": 1,
                "affinity_kcal_mol": -9.6,
                "pose_file": "/tmp/P1_A_vina.pdbqt",
            },
            {
                "engine": "gnina",
                "protein": "P1",
                "ligand": "L1",
                "site_id": "site_1",
                "tag": "P1_site_1_B",
                "pose": 1,
                "affinity_kcal_mol": -7.2,
                "pose_file": "/tmp/P1_B_gnina.pdbqt",
            },
            {
                "engine": "vina",
                "protein": "P1",
                "ligand": "L1",
                "site_id": "site_1",
                "tag": "P1_site_1_B",
                "pose": 1,
                "affinity_kcal_mol": -7.0,
                "pose_file": "/tmp/P1_B_vina.pdbqt",
            },
        ]
    )

    original_geometric = consensus_module.compute_geometric_consensus
    call_counter = {"count": 0}

    def _fake_geometric(_frame: pd.DataFrame, rmsd_cutoff: float = 2.0, **_kwargs) -> pd.DataFrame:
        call_counter["count"] += 1
        return pd.DataFrame(
            [
                {
                    "tag": "P1_site_1_A",
                    "geometric_agreement": False,
                    "geometric_engine_count": 2,
                    "geometric_expected_pairs": 1,
                    "geometric_valid_pairs": 1,
                    "geometric_pass_pairs": 0,
                    "geometric_pairwise_rmsd_min": 3.5,
                    "geometric_pairwise_rmsd_mean": 3.5,
                    "geometric_pairwise_rmsd_max": 3.5,
                    "geometric_cutoff_angstrom": rmsd_cutoff,
                    "geometric_reason": "outside_cutoff",
                    "geometric_pairwise_rmsd_json": "[]",
                },
                {
                    "tag": "P1_site_1_B",
                    "geometric_agreement": True,
                    "geometric_engine_count": 2,
                    "geometric_expected_pairs": 1,
                    "geometric_valid_pairs": 1,
                    "geometric_pass_pairs": 1,
                    "geometric_pairwise_rmsd_min": 1.1,
                    "geometric_pairwise_rmsd_mean": 1.1,
                    "geometric_pairwise_rmsd_max": 1.1,
                    "geometric_cutoff_angstrom": rmsd_cutoff,
                    "geometric_reason": "within_cutoff",
                    "geometric_pairwise_rmsd_json": "[]",
                },
            ]
        )

    consensus_module.compute_geometric_consensus = _fake_geometric
    try:
        geometric_ranked = build_consensus_rankings(
            frame,
            consensus_mode="dockbox_geometric",
            normalization_method="per_engine_rank",
        )
        weighted_ranked = build_consensus_rankings(
            frame,
            consensus_mode="weighted_hybrid",
            normalization_method="per_engine_rank",
        )
    finally:
        consensus_module.compute_geometric_consensus = original_geometric

    _assert(
        call_counter["count"] == 1,
        "Geometric consensus computation should be called only for dockbox_geometric mode",
    )
    _assert(
        str(geometric_ranked.iloc[0]["tag"]) == "P1_site_1_B",
        "dockbox_geometric mode should prioritize geometric agreement over raw affinity rank",
    )
    _assert(
        str(weighted_ranked.iloc[0]["tag"]) == "P1_site_1_A",
        "weighted_hybrid mode should remain rank-based and not inherit geometric overrides",
    )
    _assert(
        set(str(value) for value in weighted_ranked["geometric_reason"].astype(str).unique()) == {"not_computed"},
        "Non-geometric modes should keep geometric fields in not_computed state",
    )


def _smoke_prepare_guard() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_prepare_guard_"))
    original_cwd = Path.cwd()
    try:
        ligands_in = temp_root / "ligands_raw"
        ligands_out = temp_root / "ligands_prep"
        receptors_out = temp_root / "receptors_prep"
        ligands_in.mkdir(parents=True, exist_ok=True)
        ligands_out.mkdir(parents=True, exist_ok=True)
        receptors_out.mkdir(parents=True, exist_ok=True)
        os.chdir(temp_root)

        result = run_autodock_prepare(
            target="pdb.prepare_ligand",
            receptors_input=None,
            ligands_input=str(ligands_in),
            receptors_output=str(receptors_out),
            ligands_output=str(ligands_out),
            ligand_preparation_backend="meeko_only",
            selected_engines=["autodock4"],
        )
        _assert(result.status == "blocked", f"Expected blocked status for invalid profile combo, got: {result.status}")
        _assert(
            any("AutoDock4 selected" in note for note in result.notes),
            "Expected blocked result notes to include AutoDock4 compatibility guidance",
        )
    finally:
        os.chdir(original_cwd)
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_prepare_ph_guard() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_prepare_ph_guard_"))
    original_cwd = Path.cwd()
    try:
        ligands_in = temp_root / "ligands_raw"
        ligands_out = temp_root / "ligands_prep"
        receptors_out = temp_root / "receptors_prep"
        ligands_in.mkdir(parents=True, exist_ok=True)
        ligands_out.mkdir(parents=True, exist_ok=True)
        receptors_out.mkdir(parents=True, exist_ok=True)
        os.chdir(temp_root)

        result = run_autodock_prepare(
            target="pdb.prepare_ligand",
            receptors_input=None,
            ligands_input=str(ligands_in),
            receptors_output=str(receptors_out),
            ligands_output=str(ligands_out),
            ph=22.0,
            ligand_preparation_backend="openbabel_only",
            selected_engines=["vina"],
        )
        _assert(result.status == "blocked", f"Expected blocked status for invalid pH, got: {result.status}")
        ph_validation = result.outputs.get("ph_validation", {})
        _assert(ph_validation.get("is_valid") is False, "Invalid pH should report failed pH validation")
        _assert(
            "between 0.0 and 14.0" in str(ph_validation.get("error", "")),
            "Invalid pH should include range guidance",
        )
    finally:
        os.chdir(original_cwd)
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_ligand_output_contract_validator() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_ligand_contract_"))
    try:
        bad_output = temp_root / "bad_output.pdbqt"
        bad_output.write_text("REMARK malformed\n", encoding="utf-8")
        summary = {
            "input_file": str(temp_root / "ligA.sdf"),
            "output_file": str(bad_output),
            "requested_profile": "openbabel_meeko",
            "effective_profile": "openbabel_meeko",
            "selected_engines": ["vina"],
            "preparation_method": "autodocktools_direct",
            "protonation_ph": 7.4,
            "normalization": {"normalization_backend": "openbabel"},
        }
        contract = validate_ligand_preparation_output_contract(summary)
        _assert(contract.get("is_valid") is False, "Malformed/mismatched ligand contract should be invalid")
        error_blob = " | ".join([str(item) for item in contract.get("errors", [])])
        _assert(
            "preparation_method" in error_blob,
            "Contract validator should catch method/profile mismatches",
        )
        _assert(
            "no ATOM/HETATM records" in error_blob or "ROOT block" in error_blob,
            "Contract validator should catch malformed output PDBQT contents",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_preparation_preflight_artifacts() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_prepare_preflight_"))
    original_cwd = Path.cwd()
    try:
        import autodock_preparation as adprep

        ligands_in = temp_root / "ligands_raw"
        ligands_out = temp_root / "ligands_prep"
        receptors_out = temp_root / "receptors_prep"
        ligands_in.mkdir(parents=True, exist_ok=True)
        ligands_out.mkdir(parents=True, exist_ok=True)
        receptors_out.mkdir(parents=True, exist_ok=True)
        (ligands_in / "ligA.sdf").write_text("ligand", encoding="utf-8")
        os.chdir(temp_root)

        original_check_dependencies = adprep.AutoDockPreparationPipeline.check_dependencies
        original_run_enhanced_preparation = adprep.AutoDockPreparationPipeline.run_enhanced_preparation

        def _fake_check_dependencies(self):
            return True, []

        def _fake_run_enhanced_preparation(self, config_path=None):
            lig_out = Path(self.config.ligands_output).expanduser().resolve()
            lig_out.mkdir(parents=True, exist_ok=True)
            prepared = lig_out / "ligA.pdbqt"
            prepared.write_text(
                "\n".join(
                    [
                        "ROOT",
                        "HETATM    1  C   UNK A   1       0.000   0.000   0.000  1.00  0.00           C",
                        "ENDROOT",
                        "TORSDOF 0",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            step_dir = lig_out / "preparation_steps"
            step_dir.mkdir(parents=True, exist_ok=True)
            step_payload = {
                "input_file": str(ligands_in / "ligA.sdf"),
                "output_file": str(prepared),
                "preparation_method": "openbabel_pdbqt",
                "requested_profile": "openbabel_only",
                "effective_profile": "openbabel_only",
                "selected_engines": ["vina"],
                "protonation_ph": 7.4,
                "normalization": {"normalization_backend": "openbabel"},
            }
            (step_dir / "ligA.json").write_text(json.dumps(step_payload, indent=2), encoding="utf-8")
            return True

        adprep.AutoDockPreparationPipeline.check_dependencies = _fake_check_dependencies
        adprep.AutoDockPreparationPipeline.run_enhanced_preparation = _fake_run_enhanced_preparation
        try:
            result = run_autodock_prepare(
                target="pdb.prepare_ligand",
                receptors_input=None,
                ligands_input=str(ligands_in),
                receptors_output=str(receptors_out),
                ligands_output=str(ligands_out),
                ligand_preparation_backend="openbabel_only",
                selected_engines=["vina"],
            )
        finally:
            adprep.AutoDockPreparationPipeline.check_dependencies = original_check_dependencies
            adprep.AutoDockPreparationPipeline.run_enhanced_preparation = original_run_enhanced_preparation

        _assert(result.status in {"completed", "partial"}, f"Expected completed/partial status, got: {result.status}")
        for key in [
            "ligand_output_contract_validation_csv",
            "ligand_output_contract_validation_json",
            "preparation_preflight_summary_csv",
            "preparation_preflight_summary_json",
        ]:
            artifact_path = Path(str(result.outputs.get(key, "")))
            _assert(artifact_path.exists(), f"Expected preparation artifact file to exist: {key}")
    finally:
        os.chdir(original_cwd)
        shutil.rmtree(temp_root, ignore_errors=True)


def _autodocktools_script_available() -> bool:
    env_candidates = [
        Path(str(value).strip()).expanduser()
        for value in (
            str(os.environ.get("AUTODOCKTOOLS_PREPARE_LIGAND4", "") or ""),
            str(os.environ.get("ADT_PREPARE_LIGAND4", "") or ""),
        )
        if str(value).strip()
    ]
    for candidate in env_candidates:
        if candidate.exists() and candidate.is_file():
            return True

    if shutil.which("prepare_ligand4.py"):
        return True
    if shutil.which("prepare_ligand4"):
        return True

    bundled = REPO_ROOT / ".workflow" / "tools" / "autodocktools-prepare-py3k" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"
    return bundled.exists() and bundled.is_file()


def _smoke_preparation_mode_matrix() -> None:
    required_tools_missing = [tool for tool in ("obabel", "jq") if not shutil.which(tool)]
    adt_available = _autodocktools_script_available()
    prep_logger = logging.getLogger("autodock_preparation")
    previous_prep_level = prep_logger.level
    prep_logger.setLevel(logging.ERROR)
    cases = [
        ("engine_aware_full", ["gnina", "vina", "smina"]),
        ("engine_aware_full", ["gnina", "vina", "smina", "autodock4"]),
        ("openbabel_only", ["gnina"]),
        ("meeko_only", ["gnina"]),
        ("autodocktools_only", ["autodock4"]),
        ("openbabel_meeko", ["gnina", "vina", "smina"]),
        ("openbabel_meeko_autodock", ["gnina", "vina", "smina", "autodock4"]),
        ("openbabel_autodocktools", ["autodock4"]),
        ("meeko_only", ["autodock4"]),  # expected compatibility block
        ("openbabel_meeko", ["autodock4"]),  # expected compatibility block
    ]

    try:
        for profile, engines in cases:
            temp_root = Path(tempfile.mkdtemp(prefix=f"dockforge_profile_{profile}_"))
            original_cwd = Path.cwd()
            try:
                ligands_in = temp_root / "ligands_raw"
                ligands_out = temp_root / "ligands_prep"
                receptors_out = temp_root / "receptors_prep"
                ligands_in.mkdir(parents=True, exist_ok=True)
                ligands_out.mkdir(parents=True, exist_ok=True)
                receptors_out.mkdir(parents=True, exist_ok=True)
                os.chdir(temp_root)

                compatibility = validate_ligand_preparation_profile(profile, engines)
                result = run_autodock_prepare(
                    target="pdb.prepare_ligand",
                    receptors_input=None,
                    ligands_input=str(ligands_in),
                    receptors_output=str(receptors_out),
                    ligands_output=str(ligands_out),
                    ligand_preparation_backend=profile,
                    selected_engines=engines,
                )

                if not compatibility.is_valid:
                    _assert(result.status == "blocked", f"{profile}/{engines}: expected blocked for incompatibility")
                    _assert(
                        any("AutoDock4 selected" in note for note in result.notes),
                        f"{profile}/{engines}: expected AutoDock4 compatibility guidance",
                    )
                    continue

                if required_tools_missing:
                    _assert(result.status == "blocked", f"{profile}/{engines}: expected blocked when required tools are missing")
                    missing = result.outputs.get("missing_dependencies", [])
                    _assert(
                        any(tool in ",".join(missing) for tool in required_tools_missing),
                        f"{profile}/{engines}: expected missing dependency report for required tools",
                    )
                    continue

                effective_profile = resolve_effective_ligand_preparation_profile(profile, engines)
                adt_required = effective_profile in {
                    "autodocktools_only",
                    "openbabel_autodocktools",
                    "openbabel_meeko_autodock",
                }
                if adt_required and not adt_available:
                    _assert(result.status == "blocked", f"{profile}/{engines}: expected blocked when ADT script is unavailable")
                    missing = result.outputs.get("missing_dependencies", [])
                    _assert(
                        any("prepare_ligand4.py" in item for item in missing),
                        f"{profile}/{engines}: expected missing prepare_ligand4.py report",
                    )
                    continue

                _assert(
                    result.status in {"blocked", "failed", "partial", "completed"},
                    f"{profile}/{engines}: unexpected status {result.status}",
                )
                _assert(
                    result.outputs.get("ligand_input_count", 0) == 0,
                    f"{profile}/{engines}: expected zero ligand inputs in matrix smoke setup",
                )
            finally:
                os.chdir(original_cwd)
                shutil.rmtree(temp_root, ignore_errors=True)
    finally:
        prep_logger.setLevel(previous_prep_level)


def _smoke_execution_environment_and_schema() -> None:
    engines = ["gnina", "vina", "smina", "autodock4"]
    runtime_by_engine = {
        "gnina": {"binary": "gnina", "cnn_scoring": "rescore", "use_gpu": None},
        "vina": {"binary": "vina"},
        "smina": {"binary": "smina"},
        "autodock4": {"binary": "autodock4", "autogrid_binary": "autogrid4"},
    }
    config = ExecutionEnvironmentConfig(environment="local_gpu")
    env_errors, _ = validate_environment_selection(engines, config, runtime_by_engine)
    _assert(not env_errors, f"Unexpected environment selection errors: {env_errors}")
    runtime_env = apply_environment_runtime(engines, runtime_by_engine, config)
    _assert(runtime_env["gnina"].get("use_gpu") is True, "local_gpu should force GNINA GPU usage by default")

    schema = resolve_parameter_schema(
        mode="basic",
        preset="screening_fast",
        exhaustiveness=16,
        num_modes=20,
        seed=42,
        box_scale=1.0,
        box_padding=0.0,
        runtime_by_engine=runtime_env,
    )
    schema_errors, _ = validate_parameter_schema(schema, engines)
    _assert(not schema_errors, f"Unexpected schema validation errors: {schema_errors}")
    _assert(schema.common.exhaustiveness == 8, "screening_fast preset should set exhaustiveness=8")
    _assert(schema.common.num_modes == 10, "screening_fast preset should set num_modes=10")
    runtime_final = apply_schema_to_runtime(schema, runtime_env)
    _assert(runtime_final["gnina"].get("seed") == 42, "shared seed should propagate to runtime")


def _smoke_docking_preflight() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_preflight_"))
    try:
        pair_rows = [
            PairlistRow(
                receptor="missing_receptor.pdbqt",
                site_id="site_1",
                ligand="missing_ligand.pdbqt",
                center_x=0.0,
                center_y=0.0,
                center_z=0.0,
                size_x=20.0,
                size_y=20.0,
                size_z=20.0,
            )
        ]
        runtime_by_engine = {
            "gnina": {"binary": "gnina", "use_gpu": False},
            "vina": {"binary": "vina"},
        }
        preflight = run_docking_preflight(
            project_dir=temp_root,
            engines=["gnina", "vina"],
            pairlist_rows=pair_rows,
            runtime_by_engine=runtime_by_engine,
            dry_run=True,
        )
        _assert(not preflight["ok"], "preflight should fail when required receptor/ligand files are missing")
        joined_errors = " | ".join(preflight.get("errors", []))
        _assert("Missing receptor assets" in joined_errors, "missing receptor assets should be reported")
        _assert("Missing ligand assets" in joined_errors, "missing ligand assets should be reported")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_docking_preflight_receptor_qc_gate() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_preflight_receptor_qc_"))
    try:
        bootstrap_project_layout(
            temp_root,
            engines=["vina"],
            project_name="dockforge-preflight-receptor-qc",
            layout_profile="canonical",
        )
        receptors_dir = shared_receptors_dir(temp_root, "canonical")
        ligands_dir = shared_ligands_dir(temp_root, "canonical")
        receptors_dir.mkdir(parents=True, exist_ok=True)
        ligands_dir.mkdir(parents=True, exist_ok=True)

        receptor_name = "R_small.pdbqt"
        ligand_name = "L_ok.pdbqt"
        (receptors_dir / receptor_name).write_text(
            "\n".join(
                [
                    "ATOM      1  C1  UNL A   1       0.000   0.000   0.000  1.00  0.00           C",
                    "ATOM      2  O1  UNL A   1       1.200   0.000   0.000  1.00  0.00           O",
                    "END",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        (ligands_dir / ligand_name).write_text("REMARK ligand placeholder\nEND\n", encoding="utf-8")

        rows = [
            PairlistRow(
                receptor=receptor_name,
                site_id="site_1",
                ligand=ligand_name,
                center_x=0.0,
                center_y=0.0,
                center_z=0.0,
                size_x=20.0,
                size_y=20.0,
                size_z=20.0,
            )
        ]
        report_path = temp_root / "metadata" / "receptor_validation_report.json"
        preflight = run_docking_preflight(
            project_dir=temp_root,
            engines=["vina"],
            pairlist_rows=rows,
            runtime_by_engine={"vina": {"binary": "vina"}},
            dry_run=True,
            enable_ligand_qc=False,
            enable_receptor_qc=True,
            receptor_qc_report_path=str(report_path),
            receptor_min_atom_count=10,
            receptor_min_heavy_atom_count=5,
            receptor_min_chain_count=1,
            receptor_max_coordinate_span=500.0,
        )
        _assert(not preflight["ok"], "preflight should fail when receptor QC thresholds are violated")
        joined_errors = " | ".join(preflight.get("errors", []))
        _assert("Receptor QC gate failed" in joined_errors, "receptor QC failure should be reported")
        details = preflight.get("details", {})
        receptor_qc = details.get("receptor_qc", {}) if isinstance(details, dict) else {}
        _assert(
            int(receptor_qc.get("issue_count", 0) or 0) > 0,
            "receptor QC report should contain at least one issue",
        )
        _assert(report_path.exists(), "receptor QC report file should be written")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_hpc_deploy_receptor_qc_gate() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_deploy_receptor_qc_"))
    try:
        bootstrap_project_layout(
            temp_root,
            engines=["vina"],
            project_name="dockforge-deploy-receptor-qc",
            layout_profile="canonical",
        )
        receptors_dir = shared_receptors_dir(temp_root, "canonical")
        ligands_dir = shared_ligands_dir(temp_root, "canonical")
        receptors_dir.mkdir(parents=True, exist_ok=True)
        ligands_dir.mkdir(parents=True, exist_ok=True)

        receptor_name = "R_small.pdbqt"
        ligand_name = "L_ok.pdbqt"
        (receptors_dir / receptor_name).write_text(
            "\n".join(
                [
                    "ATOM      1  C1  UNL A   1       0.000   0.000   0.000  1.00  0.00           C",
                    "ATOM      2  O1  UNL A   1       1.200   0.000   0.000  1.00  0.00           O",
                    "END",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        (ligands_dir / ligand_name).write_text(
            "\n".join(
                [
                    "ROOT",
                    "ATOM      1  C   LIG A   1       0.200   0.100   0.000  1.00  0.00           C",
                    "ENDROOT",
                    "TORSDOF 0",
                    "",
                ]
            ),
            encoding="utf-8",
        )

        pair_df = pd.DataFrame(
            [
                {
                    "receptor": receptor_name,
                    "site_id": "site_1",
                    "ligand": ligand_name,
                    "center_x": 0.0,
                    "center_y": 0.0,
                    "center_z": 0.0,
                    "size_x": 20.0,
                    "size_y": 20.0,
                    "size_z": 20.0,
                }
            ]
        )
        pair_df.to_csv(pairlist_path(temp_root, "canonical"), index=False)

        status = deploy_main(
            [
                "--project-dir",
                str(temp_root),
                "--engines",
                "vina",
                "--mode",
                "screen",
                "--round",
                "qc_gate_smoke",
            ]
        )
        _assert(status == 1, f"deploy_main should fail receptor QC gate, got status={status}")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_post_docking_allscore_contracts() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_allscore_"))
    try:
        bootstrap_project_layout(
            temp_root,
            engines=["gnina", "vina"],
            project_name="dockforge-smoke",
            favorite_engine="gnina",
            layout_profile="canonical",
        )
        pairlist_file = pairlist_path(temp_root)
        pair_df = pd.DataFrame(
            [
                {
                    "receptor": "R1.pdbqt",
                    "site_id": "site_1",
                    "ligand": "L1.pdbqt",
                    "center_x": 0.0,
                    "center_y": 0.0,
                    "center_z": 0.0,
                    "size_x": 20.0,
                    "size_y": 20.0,
                    "size_z": 20.0,
                    "protein_display_name": "R1",
                    "ligand_display_name": "L1",
                },
                {
                    "receptor": "R2.pdbqt",
                    "site_id": "site_1",
                    "ligand": "L2.pdbqt",
                    "center_x": 1.0,
                    "center_y": 1.0,
                    "center_z": 1.0,
                    "size_x": 20.0,
                    "size_y": 20.0,
                    "size_z": 20.0,
                    "protein_display_name": "R2",
                    "ligand_display_name": "L2",
                },
            ]
        )
        pair_df.to_csv(pairlist_file, index=False)
        _assert(len(load_pairlist(temp_root)) == 2, "pairlist should contain two rows")

        rows = []
        for engine, offset in [("gnina", 0.0), ("vina", 0.2)]:
            for tag, protein, ligand, affinity in [
                ("R1.pdbqt_site_1_L1.pdbqt", "R1.pdbqt", "L1.pdbqt", -8.0 + offset),
                ("R2.pdbqt_site_1_L2.pdbqt", "R2.pdbqt", "L2.pdbqt", -7.2 + offset),
            ]:
                rows.append(
                    {
                        "engine": engine,
                        "tag": tag,
                        "protein": protein,
                        "ligand": ligand,
                        "site_id": "site_1",
                        "pose": 1,
                        "affinity_kcal_mol": affinity,
                        "score_name_primary": "vina_affinity",
                        "score_primary": affinity,
                        "score_name_secondary": "",
                        "score_secondary": None,
                        "rmsd_lb": None,
                        "rmsd_ub": None,
                        "pose_file": f"/tmp/{tag}.pdbqt",
                        "log_file": f"/tmp/{tag}.log",
                    }
                )
        normalized_df = pd.DataFrame(rows)
        for engine in ["gnina", "vina"]:
            layout = ensure_engine_layout(temp_root, engine, layout_profile="canonical")
            # Explicit score-only imports; placeholder poses are removed below.
            for _, row in normalized_df[normalized_df["engine"] == engine].iterrows():
                tag = str(row["tag"])
                if engine == "gnina":
                    (layout["poses"] / f"{tag}.sdf").write_text("$$$$\n", encoding="utf-8")
                else:
                    (layout["poses"] / f"{tag}.pdbqt").write_text("REMARK\n", encoding="utf-8")
                (layout["logs"] / f"{tag}.log").write_text("REMARK\n", encoding="utf-8")
            normalized_df[normalized_df["engine"] == engine].to_csv(
                layout["scores"] / "normalized_scores.csv",
                index=False,
            )
            _mark_score_import(layout)

        biology_file = temp_root / "biology_annotations.csv"
        pd.DataFrame(
            [
                {"protein": "R1.pdbqt", "ligand": "L1.pdbqt", "expression_score": 2.1, "essentiality_score": 0.8},
                {"protein": "R2.pdbqt", "ligand": "L2.pdbqt", "expression_score": 1.5, "essentiality_score": 0.4},
            ]
        ).to_csv(biology_file, index=False)

        output_dir = temp_root / "analysis" / "sessions" / "smoke_allscore"
        pipeline = MultiEngineAnalysisPipeline(
            project_dir=str(temp_root),
            output_dir=str(output_dir),
            analysis_mode="comparative_all_engines",
            analysis_scope="comparison_only",
            normalization_method="per_engine_rank",
            favorite_engine="gnina",
            biology_file=str(biology_file),
            biology_mapping_mode="protein_ligand",
            hit_class_policy="reference_anchor",
            hit_class_strong_percentile=0.10,
            hit_class_moderate_percentile=0.40,
            top_pose_selection_policy="best_affinity",
            top_pose_global_aggregation="best_target",
        )
        _assert(pipeline.run() is True, "comparative allscore pipeline should complete")

        reports_dir = output_dir / "reports"
        _assert((reports_dir / "consensus_ranked_hits.csv").exists(), "consensus ranked hits report missing")
        _assert((reports_dir / "redocking_validation.csv").exists(), "redocking validation report missing")
        _assert((reports_dir / "reference_baselines.csv").exists(), "reference baseline report missing")
        _assert((reports_dir / "VALIDATION_SUMMARY.md").exists(), "validation summary report missing")
        validation_gate_file = reports_dir / "validation_gate_status.json"
        _assert(validation_gate_file.exists(), "validation gate status report missing")
        gate_payload = json.loads(validation_gate_file.read_text(encoding="utf-8"))
        _assert(
            bool(gate_payload.get("allow_reference_anchor", True)) is False,
            "reference-anchor should be disallowed when no validated reference rows exist",
        )
        _assert(
            (reports_dir / "consensus_ranked_hits_with_classes.csv").exists(),
            "consensus ranked hits with classes report missing",
        )
        _assert(
            (reports_dir / "engine_rank_correlation_per_protein.csv").exists(),
            "per-protein engine rank correlation report missing",
        )
        _assert(
            (reports_dir / "engine_rank_correlation_global.csv").exists(),
            "global engine rank correlation report missing",
        )
        per_protein_corr_df = pd.read_csv(reports_dir / "engine_rank_correlation_per_protein.csv")
        global_corr_df = pd.read_csv(reports_dir / "engine_rank_correlation_global.csv")
        for required_col in [
            "spearman_p_value",
            "spearman_q_value",
            "pearson_p_value",
            "pearson_q_value",
            "insufficient_n",
            "normalization_method",
        ]:
            _assert(required_col in per_protein_corr_df.columns, f"per-protein correlation missing {required_col}")
            _assert(required_col in global_corr_df.columns, f"global correlation missing {required_col}")
        _assert(
            (reports_dir / "biology_mapping_report.json").exists(),
            "biology mapping report missing",
        )
        _assert(
            (reports_dir / "biology_correlation_global.csv").exists(),
            "biology correlation global report missing",
        )
        run_tracking_dir = output_dir / "run_tracking"
        _assert((run_tracking_dir / "run_manifest.json").exists(), "run tracking manifest missing")
        _assert((run_tracking_dir / "outputs_index.csv").exists(), "run outputs index missing")
        _assert((run_tracking_dir / "step_status.csv").exists(), "run step status missing")
        numbered_layout = ensure_numbered_output_layout(temp_root, "canonical")
        consolidated_summary = json.loads(
            (numbered_layout["reports_root_numbered"] / "consolidated_run_summary.json").read_text(encoding="utf-8")
        )
        run_context = dict(consolidated_summary.get("run_context", {}))
        _assert(
            str(run_context.get("hit_class_policy_requested", "")) == "reference_anchor",
            "consolidated run context should record requested hit-class policy",
        )
        _assert(
            str(run_context.get("hit_class_policy_effective", "")) == "target_percentile",
            "effective hit-class policy should fall back when validation gate disallows reference anchors",
        )
        class_df = pd.read_csv(reports_dir / "consensus_ranked_hits_with_classes.csv")
        _assert("docking_quality_class" in class_df.columns, "hit class column missing from comparative report")
        _assert("qc_status" in class_df.columns, "qc status column missing from comparative report")
        _assert("admet_status" in class_df.columns, "admet status column missing from comparative report")
        _assert("bio_expression_score" in class_df.columns, "biology annotation columns were not merged")
        step_df = pd.read_csv(run_tracking_dir / "step_status.csv")
        _assert("write_comparative_reports" in set(step_df["step"].astype(str)), "expected comparative step tracking entry")
        top_pose_dir = output_dir / "top_pose_ligand_performance"
        top_pose_per_protein_file = top_pose_dir / "top_pose_per_ligand_per_protein.csv"
        top_pose_global_file = top_pose_dir / "top_pose_per_ligand_global.csv"
        top_pose_manifest_file = top_pose_dir / "top_pose_selection_manifest.json"
        _assert(top_pose_per_protein_file.exists(), "top-pose per-protein atlas missing")
        _assert(top_pose_global_file.exists(), "top-pose global atlas missing")
        _assert(top_pose_manifest_file.exists(), "top-pose selection manifest missing")
        top_pose_manifest = json.loads(top_pose_manifest_file.read_text(encoding="utf-8"))
        _assert(
            str(top_pose_manifest.get("selection_policy", "")) == "best_affinity",
            "top-pose manifest should reflect selected policy",
        )
        _assert(
            str(top_pose_manifest.get("global_aggregation", "")) == "best_target",
            "top-pose manifest should reflect selected aggregation",
        )

        top_pose_per_protein = pd.read_csv(top_pose_per_protein_file)
        top_pose_global = pd.read_csv(top_pose_global_file)
        _assert(
            int(top_pose_global["ligand"].nunique()) == int(len(top_pose_global)),
            "top-pose global atlas should contain exactly one row per ligand",
        )
        pair_counts = (
            top_pose_per_protein.groupby(["ligand", "protein"], dropna=False)
            .size()
            .reset_index(name="n")
        )
        _assert(
            bool(pair_counts["n"].max() == 1),
            "top-pose per-protein atlas should contain exactly one row per ligand/protein pair",
        )

        outputs_index = pd.read_csv(run_tracking_dir / "outputs_index.csv")
        indexed_paths = set(outputs_index["relative_path"].astype(str))
        _assert(
            "top_pose_ligand_performance/top_pose_per_ligand_global.csv" in indexed_paths,
            "run outputs index should include top-pose atlas files",
        )

        top_pose_only_output = temp_root / "analysis" / "sessions" / "smoke_top_pose_only"
        top_pose_only_pipeline = MultiEngineAnalysisPipeline(
            project_dir=str(temp_root),
            output_dir=str(top_pose_only_output),
            analysis_mode="comparative_all_engines",
            analysis_scope="top_pose_only",
            normalization_method="per_engine_rank",
            favorite_engine="gnina",
            biology_file=str(biology_file),
            biology_mapping_mode="protein_ligand",
            hit_class_policy="target_percentile",
            hit_class_strong_percentile=0.10,
            hit_class_moderate_percentile=0.40,
            top_pose_selection_policy="best_affinity",
            top_pose_global_aggregation="best_target",
        )
        _assert(top_pose_only_pipeline.run() is True, "top_pose_only scope run should complete")
        top_pose_only_global = pd.read_csv(
            top_pose_only_output / "top_pose_ligand_performance" / "top_pose_per_ligand_global.csv"
        )
        pd.testing.assert_frame_equal(
            top_pose_global.sort_values(["ligand", "protein", "tag"]).reset_index(drop=True),
            top_pose_only_global.sort_values(["ligand", "protein", "tag"]).reset_index(drop=True),
            check_dtype=False,
        )

        canonical_scores = ensure_numbered_output_layout(temp_root)["post_scores_consensus"].parent
        _assert((canonical_scores / "raw" / "all_scores_raw.csv").exists(), "working raw allscore output missing")
        _assert(
            (canonical_scores / "unified" / "engine_rank_correlation_global.csv").exists(),
            "working unified correlation output missing",
        )
        canonical_top_pose = post_docking_root(temp_root, "canonical") / "top_pose_ligand_performance"
        _assert(
            (canonical_top_pose / "top_pose_per_ligand_global.csv").exists(),
            "canonical top-pose global atlas missing",
        )
        _assert(
            (canonical_scores / "consensus" / "consensus_ranked_hits_with_classes.csv").exists(),
            "canonical consensus hit classes output missing",
        )
        _assert(
            (canonical_scores / "consensus" / "validation_gate_status.json").exists(),
            "canonical validation gate mirror output missing",
        )
        explainability_file = canonical_scores / "consensus" / "consensus_explainability.json"
        _assert(explainability_file.exists(), "canonical consensus explainability output missing")
        explainability = json.loads(explainability_file.read_text(encoding="utf-8"))
        _assert(
            str(explainability.get("normalization_method", "")) == "per_engine_rank",
            "explainability should record normalization method",
        )
        _assert(
            isinstance(explainability.get("hit_class_distribution", {}), dict),
            "explainability should include hit class distribution",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_biology_unresolved_mapping_reporting() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_biology_unresolved_"))
    try:
        bootstrap_project_layout(
            temp_root,
            engines=["gnina", "vina"],
            project_name="dockforge-biology-unresolved",
            favorite_engine="gnina",
            layout_profile="canonical",
        )
        pair_df = pd.DataFrame(
            [
                {
                    "receptor": "R1.pdbqt",
                    "site_id": "site_1",
                    "ligand": "L1.pdbqt",
                    "center_x": 0.0,
                    "center_y": 0.0,
                    "center_z": 0.0,
                    "size_x": 20.0,
                    "size_y": 20.0,
                    "size_z": 20.0,
                    "protein_display_name": "R1",
                    "ligand_display_name": "L1",
                }
            ]
        )
        pair_df.to_csv(pairlist_path(temp_root), index=False)

        normalized_rows = []
        for engine, offset in [("gnina", 0.0), ("vina", 0.2)]:
            normalized_rows.append(
                {
                    "engine": engine,
                    "tag": "R1.pdbqt_site_1_L1.pdbqt",
                    "protein": "R1.pdbqt",
                    "ligand": "L1.pdbqt",
                    "site_id": "site_1",
                    "pose": 1,
                    "affinity_kcal_mol": -8.0 + offset,
                    "score_name_primary": "vina_affinity",
                    "score_primary": -8.0 + offset,
                    "score_name_secondary": "",
                    "score_secondary": None,
                    "rmsd_lb": None,
                    "rmsd_ub": None,
                    "pose_file": f"/tmp/{engine}_R1_L1.pdbqt",
                    "log_file": f"/tmp/{engine}_R1_L1.log",
                }
            )
        normalized_df = pd.DataFrame(normalized_rows)
        for engine in ["gnina", "vina"]:
            layout = ensure_engine_layout(temp_root, engine, layout_profile="canonical")
            for _, row in normalized_df[normalized_df["engine"] == engine].iterrows():
                tag = str(row["tag"])
                if engine == "gnina":
                    (layout["poses"] / f"{tag}.sdf").write_text("$$$$\n", encoding="utf-8")
                else:
                    (layout["poses"] / f"{tag}.pdbqt").write_text("REMARK\n", encoding="utf-8")
                (layout["logs"] / f"{tag}.log").write_text("REMARK\n", encoding="utf-8")
            normalized_df[normalized_df["engine"] == engine].to_csv(
                layout["scores"] / "normalized_scores.csv",
                index=False,
            )
            _mark_score_import(layout)

        biology_file = temp_root / "biology_annotations_unresolved.csv"
        pd.DataFrame(
            [
                {"protein": "R1.pdbqt", "ligand": "L1.pdbqt", "expression_score": 2.0},
                {"protein": "R9.pdbqt", "ligand": "L9.pdbqt", "expression_score": 9.0},
            ]
        ).to_csv(biology_file, index=False)

        output_dir = temp_root / "analysis" / "sessions" / "smoke_biology_unresolved"
        pipeline = MultiEngineAnalysisPipeline(
            project_dir=str(temp_root),
            output_dir=str(output_dir),
            analysis_mode="comparative_all_engines",
            analysis_scope="comparison_only",
            normalization_method="per_engine_rank",
            favorite_engine="gnina",
            biology_file=str(biology_file),
            biology_mapping_mode="protein_ligand",
            hit_class_policy="target_percentile",
            hit_class_strong_percentile=0.10,
            hit_class_moderate_percentile=0.40,
            top_pose_selection_policy="best_affinity",
            top_pose_global_aggregation="best_target",
        )
        _assert(pipeline.run() is True, "comparative run with unresolved biology keys should complete")

        mapping_report_file = output_dir / "reports" / "biology_mapping_report.json"
        _assert(mapping_report_file.exists(), "biology mapping report should be written for unresolved-key scenario")
        mapping_report = json.loads(mapping_report_file.read_text(encoding="utf-8"))
        _assert(
            int(mapping_report.get("unmatched_biology_keys", 0)) >= 1,
            "biology mapping report should count unmatched biology keys",
        )
        _assert(
            isinstance(mapping_report.get("unmatched_biology_examples", []), list)
            and len(mapping_report.get("unmatched_biology_examples", [])) >= 1,
            "biology mapping report should include unresolved biology key examples",
        )

        session_poly_root = output_dir / "reports" / "polypharmacology"
        _assert(
            (session_poly_root / "biology_unmatched_biology_keys.csv").exists(),
            "session polypharmacology should include unresolved biology key export",
        )
        unresolved_rows = pd.read_csv(session_poly_root / "biology_unmatched_biology_keys.csv")
        _assert(len(unresolved_rows) >= 1, "unmatched biology key export should contain at least one row")

        canonical_poly_root = post_docking_root(temp_root, "canonical") / "polypharmacology"
        _assert(
            (canonical_poly_root / "biology_unmatched_biology_keys.csv").exists(),
            "canonical polypharmacology should mirror unresolved biology key export",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_plugin_system_retired_contract() -> None:
    try:
        import post_docking_analysis.plugin_manager  # noqa: F401
        raise RuntimeError("plugin_manager module should be removed under FR-011 cleanup")
    except ModuleNotFoundError:
        pass


def _smoke_complex_export_index_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_complex_index_"))
    try:
        bridge_dir = temp_root / "deep_analysis"
        complexes_dir = bridge_dir / "complexes"
        best_poses_dir = bridge_dir / "best_poses_pdb"
        complexes_dir.mkdir(parents=True, exist_ok=True)
        best_poses_dir.mkdir(parents=True, exist_ok=True)
        (complexes_dir / "pair_A.pdb").write_text(
            "\n".join(
                [
                    "ATOM      1  N   ALA A   1      11.104   8.177   2.255  1.00 20.00           N",
                    "HETATM    2  C   UNK X   1      12.000   9.000   3.000  1.00 20.00           C",
                    "END",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        (complexes_dir / "pair_B.pdb").write_text("END\n", encoding="utf-8")
        (best_poses_dir / "pair_A.pdb").write_text("END\n", encoding="utf-8")

        source_table = pd.DataFrame(
            [
                {
                    "tag": "pair_A",
                    "protein": "R1.pdbqt",
                    "ligand": "L1.pdbqt",
                    "site_id": "site_1",
                    "pose": 1,
                    "affinity_kcal_mol": -8.1,
                }
            ]
        )
        index_file = MultiEngineAnalysisPipeline._write_complex_export_index(
            output_dir=bridge_dir,
            complexes_dir=complexes_dir,
            best_poses_dir=best_poses_dir,
            favorite_engine="vina",
            source_table=source_table,
        )
        _assert(index_file.exists(), "complex export index file should exist")
        index_df = pd.read_csv(index_file)
        _assert(len(index_df) == 2, "complex export index should contain one row per complex")
        _assert(str(index_df.loc[0, "tag"]) == "pair_A", "complex export index tag mismatch")
        valid_map = {str(row["tag"]): bool(row["is_valid_complex"]) for _, row in index_df.iterrows()}
        _assert(valid_map.get("pair_A") is True, "pair_A should be marked as a valid complex export")
        _assert(valid_map.get("pair_B") is False, "pair_B should be marked as an invalid complex export")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_unified_output_topology_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_output_topology_"))
    try:
        ensure_project_layout(temp_root, "docking_legacy")
        numbered = ensure_numbered_output_layout(temp_root)
        analysis_root = post_docking_root(temp_root)
        session_root = analysis_root / "sessions" / "interactive_20260401_120000"
        reports_root = temp_root / "7-Reports"
        for child in (
            session_root / "complexes",
            session_root / "scores",
            session_root / "best_poses",
            session_root / "reports",
            reports_root,
        ):
            child.mkdir(parents=True, exist_ok=True)

        (session_root / "complexes" / "pair_A.pdb").write_text("END\n", encoding="utf-8")
        (session_root / "reports" / "consensus_ranked_hits.csv").write_text("tag\npair_A\n", encoding="utf-8")

        shim = ensure_post_docking_compat_shim(temp_root)
        _assert(shim["legacy_root"].exists(), "legacy compatibility root should exist")
        _assert(
            shim["legacy_root"].is_symlink() or shim["status"] == "legacy_dir_exists",
            "legacy compatibility root should be a symlink unless a prior legacy dir exists",
        )
        _assert(numbered["post_docking_root_numbered"].name == "4-Working", "working root should be renamed to 4-Working")

        _refresh_latest_analysis_shortcuts(temp_root, session_root)
        index_file = Path(
            generate_analysis_root_index(
                analysis_root,
                project_root=temp_root,
                session_root=session_root,
                reports_root=reports_root,
            )
        )
        _assert(index_file.exists(), "5-Analysis/START_HERE.md should be generated")
        contents = index_file.read_text(encoding="utf-8")
        _assert("complexes" in contents, "analysis index should mention complexes path")
        _assert("top_pose_ligand_performance" in contents, "analysis index should mention top-pose atlas path")
        _assert("7-Reports" in contents, "analysis index should mention consolidated reports path")

        latest_link = analysis_root / "LATEST_SESSION"
        if latest_link.is_symlink():
            _assert(latest_link.resolve() == session_root.resolve(), "LATEST_SESSION should resolve to current session")
        else:
            pointer = latest_link.with_suffix(".txt")
            _assert(pointer.exists(), "LATEST_SESSION should provide a pointer when symlinks are unavailable")
            _assert(Path(pointer.read_text(encoding="utf-8").strip()) == session_root.resolve(), "pointer must identify the current session")

        raw_copy = analysis_root / "raw_data" / "START_HERE.md"
        _assert(raw_copy.exists(), "raw_data duplicate START_HERE.md should exist")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _write_test_pairlist(project_root: Path) -> None:
    pairlist = pd.DataFrame(
        [
            {
                "receptor": "R1",
                "site_id": "site_1",
                "ligand": "L1",
                "center_x": 0.0,
                "center_y": 0.0,
                "center_z": 0.0,
                "size_x": 20.0,
                "size_y": 20.0,
                "size_z": 20.0,
            }
        ]
    )
    pairlist.to_csv(pairlist_path(project_root), index=False)


def _compact_test_coords(
    count: int,
    *,
    offset: float = 0.0,
) -> List[Tuple[float, float, float]]:
    coords: List[Tuple[float, float, float]] = []
    for idx in range(max(int(count), 0)):
        x = (idx % 7) * 0.07 + offset
        y = ((idx // 7) % 7) * 0.07 + (offset * 0.5)
        z = (idx // 49) * 0.07 + (offset * 0.25)
        coords.append((float(x), float(y), float(z)))
    return coords


def _write_test_pdbqt_pose(path: Path, coords: List[Tuple[float, float, float]]) -> None:
    lines = ["MODEL        1"]
    for index, (x, y, z) in enumerate(coords, start=1):
        atom_name = f"C{index % 10}"
        lines.append(
            f"HETATM{index:5d} {atom_name:<4} LIG A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
    lines.extend(["ENDMDL", "END"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _build_sdf_pose_block(
    coords: List[Tuple[float, float, float]],
    *,
    title: str,
    declared_atom_count: Optional[int] = None,
) -> str:
    natoms = int(declared_atom_count) if declared_atom_count is not None else len(coords)
    lines = [
        title,
        "  DockForge",
        "",
        f"{natoms:>3d}{max(len(coords) - 1, 0):>3d}  0  0  0  0            999 V2000",
    ]
    for x, y, z in coords:
        lines.append(
            f"{x:10.4f}{y:10.4f}{z:10.4f} C   0  0  0  0  0  0  0  0  0  0  0  0"
        )
    for index in range(1, len(coords)):
        lines.append(f"{index:3d}{index + 1:3d}  1  0  0  0  0")
    lines.append("M  END")
    return "\n".join(lines) + "\n$$$$\n"


def _write_test_sdf(
    path: Path,
    blocks: List[List[Tuple[float, float, float]]],
    *,
    declared_atom_counts: Optional[List[Optional[int]]] = None,
) -> None:
    chunks: List[str] = []
    counts = declared_atom_counts or [None] * len(blocks)
    for index, coords in enumerate(blocks, start=1):
        declared = counts[index - 1] if index - 1 < len(counts) else None
        chunks.append(
            _build_sdf_pose_block(
                coords,
                title=f"Pose_{index}",
                declared_atom_count=declared,
            )
        )
    path.write_text("".join(chunks), encoding="utf-8")


def _seed_engine_outputs(project_root: Path, engines: Tuple[str, ...]) -> None:
    for engine in engines:
        layout = ensure_engine_layout(project_root, engine)
        if engine == "gnina":
            (shared_receptors_dir(project_root) / "R1_site_prep.pdbqt").write_text(
                "\n".join(
                    [
                        "ATOM      1  N   GLY A   1       0.000   0.000   0.000  0.00  0.00      A    N",
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            (layout["poses"] / "R1_site_1_L1.sdf").write_text(
                "\n".join(
                    [
                        "L1",
                        "  DockForge",
                        "",
                        "  1  0  0  0  0  0            999 V2000",
                        "    1.0000    2.0000    3.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
                        "M  END",
                        "> <minimizedAffinity>",
                        "-8.10",
                        "",
                        "> <CNNscore>",
                        "0.61",
                        "",
                        "> <CNNaffinity>",
                        "7.40",
                        "",
                        "$$$$",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            (layout["logs"] / "R1_site_1_L1.log").write_text(
                "\n".join(
                    [
                        "mode | affinity | intramol | CNNscore | cnn_affinity",
                        "    1      -8.10       -0.50       0.61      -7.40",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
        elif engine == "autodock4":
            (layout["poses"] / "R1_site_1_L1.dlg").write_text(
                "\n".join(
                    [
                        "DOCKED: MODEL 1",
                        "DOCKED: USER    Run = 1",
                        "DOCKED: USER    Estimated Free Energy of Binding    =   -7.20 kcal/mol",
                        "DOCKED: MODEL 2",
                        "DOCKED: USER    Run = 2",
                        "DOCKED: USER    Estimated Free Energy of Binding    =   -6.80 kcal/mol",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            (layout["logs"] / "R1_site_1_L1.log").write_text(
                "autodock4 run log placeholder\n",
                encoding="utf-8",
            )
        else:
            (layout["poses"] / "R1_site_1_L1.pdbqt").write_text(
                "\n".join(
                    [
                        "MODEL 1",
                        "REMARK VINA RESULT: -7.20 0.00 0.00",
                        "ENDMDL",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            log_lines = [
                "mode | affinity | rmsd_lb | rmsd_ub",
                "    1      -7.20      0.00      0.00",
            ]
            if engine == "smina":
                log_lines = [
                    "smina is based off AutoDock Vina",
                    "gauss 0.035579",
                    "repulsion 0.840245",
                    "hydrophobic -0.035069",
                    "non_dir_h_bond -0.587439",
                    "num_tors_div 1.923000",
                ] + log_lines
            (layout["logs"] / "R1_site_1_L1.log").write_text("\n".join(log_lines) + "\n", encoding="utf-8")


def _write_visualization_fixture(
    root: Path,
    *,
    engines: Tuple[str, ...],
    include_reference: bool = True,
    validation_state: str = "validated",
) -> dict:
    root.mkdir(parents=True, exist_ok=True)
    proteins = ["P1", "P2"]
    ligands = ["L1", "L2"] + (["REF1"] if include_reference else [])
    base_affinity = {"L1": -8.8, "L2": -7.2, "REF1": -8.0}

    normalized_rows = []
    consensus_rows = []
    classified_rows = []
    agreement_rows = []
    top_pose_rows = []
    for protein_index, protein in enumerate(proteins):
        for ligand in ligands:
            tag = f"{protein}_site_1_{ligand}"
            best_affinity = base_affinity.get(ligand, -7.0) + (0.15 * protein_index)
            quality = "Strong" if best_affinity <= -8.3 else "Moderate" if best_affinity <= -7.2 else "Weak"
            for engine_index, engine in enumerate(engines):
                affinity = best_affinity + (0.08 * engine_index)
                normalized_rows.append(
                    {
                        "engine": engine,
                        "tag": tag,
                        "protein": protein,
                        "ligand": ligand,
                        "site_id": "site_1",
                        "pose": 1,
                        "affinity_kcal_mol": affinity,
                        "cnn_affinity": affinity + 0.5 if engine == "gnina" else float("nan"),
                        "cnn_score": 0.72 if engine == "gnina" else float("nan"),
                        "rmsd_lb": 0.6 if engine == "vina" else 0.0,
                        "rmsd_ub": 1.1 if engine == "vina" else 0.0,
                    }
                )
            consensus_rows.append(
                {
                    "tag": tag,
                    "protein": protein,
                    "ligand": ligand,
                    "consensus_score": abs(best_affinity),
                    "docking_quality_class": quality,
                }
            )
            classified_rows.append(
                {
                    "tag": tag,
                    "protein": protein,
                    "ligand": ligand,
                    "best_affinity_kcal_mol": best_affinity,
                    "affinity_kcal_mol": best_affinity,
                    "docking_quality_class": quality,
                    "ligand_type": "reference" if ligand == "REF1" else "series",
                    "is_cocrystal_benchmark": bool(ligand == "REF1"),
                }
            )
            agreement_rows.append(
                {
                    "tag": tag,
                    "engine_support_count": len(engines) if ligand == "L1" else max(1, len(engines) - 1),
                }
            )
            top_pose_rows.append(
                {
                    "tag": tag,
                    "protein": protein,
                    "ligand": ligand,
                    "best_affinity_kcal_mol": best_affinity,
                    "docking_quality_class": quality,
                }
            )

    normalized_scores = root / "normalized_scores.csv"
    consensus_ranked = root / "consensus_ranked.csv"
    classified_hits = root / "classified_hits.csv"
    engine_agreement = root / "engine_agreement.csv"
    top_pose_global = root / "top_pose_per_ligand_global.csv"
    engine_scope_config = root / "engine_scope_config.json"
    validation_gate = root / "validation_gate.json"
    protein_name_mapping = root / "protein_name_mapping.csv"
    visualizations_dir = root / "visualizations"

    pd.DataFrame(normalized_rows).to_csv(normalized_scores, index=False)
    pd.DataFrame(consensus_rows).to_csv(consensus_ranked, index=False)
    pd.DataFrame(classified_rows).to_csv(classified_hits, index=False)
    pd.DataFrame(agreement_rows).to_csv(engine_agreement, index=False)
    pd.DataFrame(top_pose_rows).to_csv(top_pose_global, index=False)
    engine_scope_config.write_text(json.dumps({"engines_in_scope": list(engines)}, indent=2), encoding="utf-8")

    gate_payload = {"status": validation_state, "validation_gate_state": validation_state}
    validation_gate.write_text(json.dumps(gate_payload, indent=2), encoding="utf-8")
    pd.DataFrame(
        [
            {"record_type": "receptor", "key": "P1", "display_name": "Cyclooxygenase-2 (3LN1)"},
            {"record_type": "receptor", "key": "P2", "display_name": "EGFR kinase (1M17)"},
        ]
    ).to_csv(protein_name_mapping, index=False)

    validation_rows = []
    reference_rows = []
    if validation_state == "validated":
        for protein in proteins:
            reference_rows.append({"protein": protein, "reference_affinity": -8.0})
            ref_ligand = "REF1" if include_reference else "L1"
            for engine in engines:
                validation_rows.append(
                    {
                        "protein": protein,
                        "ligand": ref_ligand,
                        "tag": f"{protein}_site_1_{ref_ligand}",
                        "engine": engine,
                        "redocking_rmsd_angstrom": 1.4,
                        "redocking_classification": "pass",
                    }
                )
    pd.DataFrame(
        validation_rows,
        columns=[
            "protein",
            "ligand",
            "tag",
            "engine",
            "redocking_rmsd_angstrom",
            "redocking_classification",
        ],
    ).to_csv(root / "redocking_validation.csv", index=False)
    pd.DataFrame(reference_rows, columns=["protein", "reference_affinity"]).to_csv(
        root / "reference_baselines.csv",
        index=False,
    )

    return {
        "classified_hits": classified_hits,
        "consensus_ranked": consensus_ranked,
        "engine_agreement": engine_agreement,
        "normalized_scores": normalized_scores,
        "engine_scope_config": engine_scope_config,
        "validation_gate": validation_gate,
        "top_pose_global": top_pose_global,
        "protein_name_mapping": protein_name_mapping,
        "visualizations_dir": visualizations_dir,
    }


def _manifest_figure_row(manifest: dict, figure_name: str) -> dict:
    for row in manifest.get("figures", []) or []:
        if isinstance(row, dict) and str(row.get("name", "")) == figure_name:
            return row
    return {}


def _manifest_skipped_row(manifest: dict, figure_name: str) -> dict:
    for row in manifest.get("skipped", []) or []:
        if isinstance(row, dict) and str(row.get("name", "")) == figure_name:
            return row
    return {}


def _smoke_visualization_suite_core_contracts() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_viz_suite_"))
    try:
        fixture = _write_visualization_fixture(
            temp_root,
            engines=("gnina", "vina", "smina"),
            include_reference=True,
            validation_state="validated",
        )
        manifest = generate_visualization_suite(
            classified_hits_file=fixture["classified_hits"],
            consensus_ranked_file=fixture["consensus_ranked"],
            engine_agreement_file=fixture["engine_agreement"],
            normalized_scores_file=fixture["normalized_scores"],
            engine_scope_config_file=fixture["engine_scope_config"],
            validation_gate_file=fixture["validation_gate"],
            top_pose_global_file=fixture["top_pose_global"],
            output_dir=fixture["visualizations_dir"],
        )
        manifest_file = fixture["visualizations_dir"] / "figures_manifest.json"
        _assert(manifest_file.exists(), "visualization suite should write figures_manifest.json")
        loaded = json.loads(manifest_file.read_text(encoding="utf-8"))
        _assert(isinstance(loaded, dict), "figures_manifest should be a JSON object")
        _assert(int((loaded.get("summary") or {}).get("generated", 0)) >= 14, "multi-engine validated run should generate >=14 figures")
        _assert(str(loaded.get("engine_mode", "")) == "multi_engine", "multi-engine fixture should report multi_engine mode")
        _assert(manifest.get("summary", {}).get("generated") == loaded.get("summary", {}).get("generated"), "returned manifest should match persisted manifest")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_visualization_suite_single_engine_conditionals() -> None:
    gnina_root = Path(tempfile.mkdtemp(prefix="dockforge_viz_gnina_"))
    try:
        fixture = _write_visualization_fixture(
            gnina_root,
            engines=("gnina",),
            include_reference=True,
            validation_state="validated",
        )
        manifest = generate_visualization_suite(
            classified_hits_file=fixture["classified_hits"],
            consensus_ranked_file=fixture["consensus_ranked"],
            engine_agreement_file=fixture["engine_agreement"],
            normalized_scores_file=fixture["normalized_scores"],
            engine_scope_config_file=fixture["engine_scope_config"],
            validation_gate_file=fixture["validation_gate"],
            top_pose_global_file=fixture["top_pose_global"],
            output_dir=fixture["visualizations_dir"],
        )
        cnn_row = _manifest_figure_row(manifest, "cnn_confidence_distribution.png")
        _assert(bool(cnn_row.get("generated", False)), "GNINA solo run should generate cnn_confidence_distribution.png")
        _assert(not _manifest_figure_row(manifest, "engine_rank_correlation.png"), "GNINA solo run should not include multi-engine rank-correlation figure")
        _assert(not _manifest_figure_row(manifest, "engine_agreement_per_complex.png"), "GNINA solo run should not include multi-engine agreement figure")
    finally:
        shutil.rmtree(gnina_root, ignore_errors=True)

    vina_root = Path(tempfile.mkdtemp(prefix="dockforge_viz_vina_"))
    try:
        fixture = _write_visualization_fixture(
            vina_root,
            engines=("vina",),
            include_reference=True,
            validation_state="validated",
        )
        manifest = generate_visualization_suite(
            classified_hits_file=fixture["classified_hits"],
            consensus_ranked_file=fixture["consensus_ranked"],
            engine_agreement_file=fixture["engine_agreement"],
            normalized_scores_file=fixture["normalized_scores"],
            engine_scope_config_file=fixture["engine_scope_config"],
            validation_gate_file=fixture["validation_gate"],
            top_pose_global_file=fixture["top_pose_global"],
            output_dir=fixture["visualizations_dir"],
        )
        pose_row = _manifest_figure_row(manifest, "pose_diversity.png")
        _assert(bool(pose_row.get("generated", False)), "Vina solo run should generate pose_diversity.png")
        skipped_cnn = _manifest_skipped_row(manifest, "cnn_confidence_distribution.png")
        _assert(
            str(skipped_cnn.get("reason", "")) == "gnina_not_in_scope",
            "Vina solo run should skip cnn confidence with gnina_not_in_scope",
        )
    finally:
        shutil.rmtree(vina_root, ignore_errors=True)


def _smoke_visualization_suite_skip_policies_and_no_raise() -> None:
    no_ref_root = Path(tempfile.mkdtemp(prefix="dockforge_viz_noref_"))
    try:
        fixture = _write_visualization_fixture(
            no_ref_root,
            engines=("gnina", "vina"),
            include_reference=False,
            validation_state="needs_review",
        )
        manifest = generate_visualization_suite(
            classified_hits_file=fixture["classified_hits"],
            consensus_ranked_file=fixture["consensus_ranked"],
            engine_agreement_file=fixture["engine_agreement"],
            normalized_scores_file=fixture["normalized_scores"],
            engine_scope_config_file=fixture["engine_scope_config"],
            validation_gate_file=fixture["validation_gate"],
            top_pose_global_file=fixture["top_pose_global"],
            output_dir=fixture["visualizations_dir"],
        )
        skipped_poly = _manifest_skipped_row(manifest, "polypharmacology_delta_heatmap.png")
        _assert(str(skipped_poly.get("reason", "")) == "no_references", "polypharmacology delta heatmap should be skipped when references are absent")
        for figure_name in ("rmsd_distribution.png", "reference_vs_docked.png", "rmsd_per_complex.png"):
            skipped_row = _manifest_skipped_row(manifest, figure_name)
            _assert(
                str(skipped_row.get("reason", "")) == "no_rmsd_data",
                f"{figure_name} should be skipped when RMSD data is unavailable",
            )

        def _boom(*_args, **_kwargs):
            raise RuntimeError("intentional smoke failure")

        injected = generate_visualization_suite(
            classified_hits_file=fixture["classified_hits"],
            consensus_ranked_file=fixture["consensus_ranked"],
            engine_agreement_file=fixture["engine_agreement"],
            normalized_scores_file=fixture["normalized_scores"],
            engine_scope_config_file=fixture["engine_scope_config"],
            validation_gate_file=fixture["validation_gate"],
            top_pose_global_file=fixture["top_pose_global"],
            output_dir=fixture["visualizations_dir"] / "fault_injected",
            figure_overrides={"score_distribution.png": _boom},
        )
        failed_row = _manifest_figure_row(injected, "score_distribution.png")
        _assert(not bool(failed_row.get("generated", False)), "figure override exception should mark figure as not generated")
        _assert(str(failed_row.get("name", "")) == "score_distribution.png", "injected failure should be recorded for score_distribution figure")
        generated_any = int((injected.get("summary") or {}).get("generated", 0))
        _assert(generated_any > 0, "suite should continue generating other figures after a single figure exception")
    finally:
        shutil.rmtree(no_ref_root, ignore_errors=True)


def _smoke_visualization_suite_scientific_guardrails() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_viz_guardrails_"))
    try:
        fixture = _write_visualization_fixture(
            temp_root,
            engines=("gnina", "vina", "smina"),
            include_reference=True,
            validation_state="validated",
        )
        manifest = generate_visualization_suite(
            classified_hits_file=fixture["classified_hits"],
            consensus_ranked_file=fixture["consensus_ranked"],
            engine_agreement_file=fixture["engine_agreement"],
            normalized_scores_file=fixture["normalized_scores"],
            engine_scope_config_file=fixture["engine_scope_config"],
            validation_gate_file=fixture["validation_gate"],
            top_pose_global_file=fixture["top_pose_global"],
            output_dir=fixture["visualizations_dir"],
            protein_name_mapping_file=fixture["protein_name_mapping"],
        )
        for figure_name in (
            "engine_rank_correlation.png",
            "engine_agreement_per_complex.png",
            "reference_vs_docked.png",
            "rmsd_per_complex.png",
        ):
            row = _manifest_figure_row(manifest, figure_name)
            _assert(bool(row.get("generated", False)), f"{figure_name} should generate for validated multi-engine fixtures")

        ctx = _build_plot_context(
            classified_hits_file=fixture["classified_hits"],
            consensus_ranked_file=fixture["consensus_ranked"],
            engine_agreement_file=fixture["engine_agreement"],
            normalized_scores_file=fixture["normalized_scores"],
            engine_scope_config_file=fixture["engine_scope_config"],
            validation_gate_file=fixture["validation_gate"],
            top_pose_global_file=fixture["top_pose_global"],
            output_dir=fixture["visualizations_dir"],
            protein_name_mapping_file=fixture["protein_name_mapping"],
        )
        summary = _build_engine_agreement_summary(ctx)
        _assert(not summary.empty, "engine agreement summary should not be empty for multi-engine fixture")
        _assert((summary["full_support_fraction"] < 1.0).any(), "agreement summary should preserve partial-support ligands and not collapse everything to 1.0")
        _assert("[Ref]" in _ligand_label("3LN1_ligand_CEL_A_682.pdbqt"), "co-crystal ligand labels should carry an explicit [Ref] marker")
        _assert("pose-collapse" in str(_pose_diversity_diagnostic(pd.Series([0.0, 0.0, 0.0]))), "all-zero pose diversity should produce a collapse diagnostic")

        broken_validation = pd.read_csv(fixture["validation_gate"].parent / "redocking_validation.csv")
        broken_validation["engine"] = "bogus_engine"
        broken_validation.to_csv(fixture["validation_gate"].parent / "redocking_validation.csv", index=False)
        fallback_ctx = _build_plot_context(
            classified_hits_file=fixture["classified_hits"],
            consensus_ranked_file=fixture["consensus_ranked"],
            engine_agreement_file=fixture["engine_agreement"],
            normalized_scores_file=fixture["normalized_scores"],
            engine_scope_config_file=fixture["engine_scope_config"],
            validation_gate_file=fixture["validation_gate"],
            top_pose_global_file=fixture["top_pose_global"],
            output_dir=fixture["visualizations_dir"] / "fallback_ctx",
            protein_name_mapping_file=fixture["protein_name_mapping"],
        )
        joined = _validation_joined_table(fallback_ctx)
        _assert(joined["affinity_kcal_mol"].isna().all(), "validation plots must not backfill raw affinity from a different engine")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_visualization_suite_start_here_embedding_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_viz_start_here_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina"))
        output_dir = post_docking_root(temp_root) / "sessions" / "viz_start_here"
        pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(output_dir),
            analysis_mode="comparative_all_engines",
            analysis_scope="report_only",
            engines="gnina,vina,smina",
        )
        _stub_dag_optional_interactions(pipeline, sleep_s=0.01)
        _assert(pipeline.run() is True, "report-only DAG run should succeed for visualization embedding contract")
        manifest_file = output_dir / "visualizations" / "figures_manifest.json"
        _assert(manifest_file.exists(), "visualization DAG node should produce figures_manifest.json")
        manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
        generated_rows = [row for row in (manifest.get("figures") or []) if isinstance(row, dict) and bool(row.get("generated", False))]
        _assert(bool(generated_rows), "visualization manifest should contain at least one generated figure")
        first_path = Path(str(generated_rows[0].get("path", ""))).expanduser()
        try:
            rel_path = str(first_path.resolve().relative_to(temp_root.resolve()))
        except Exception:
            rel_path = str(first_path)

        numbered = ensure_numbered_output_layout(temp_root)
        start_here = numbered["reports_root_numbered"] / "START_HERE.md"
        _assert(start_here.exists(), "START_HERE.md should exist after report-only DAG run")
        text = start_here.read_text(encoding="utf-8")
        _assert("## Visualization Suite" in text, "START_HERE should include Visualization Suite section")
        _assert(rel_path in text, "START_HERE should list at least one generated figure path from figures_manifest")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_engine_detection_and_auto_routing() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_engine_detect_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina"))

        detected = detect_engines(temp_root, redetect=True)
        _assert(detected["detection_method"] == "manifest", "engine detection should prefer manifest engines")
        _assert(detected["routing_decision"] == "multi_engine", "three valid engines should route to multi_engine")
        _assert(Path(str(detected["report_file"])).exists(), "engine detection report should be written")

        override_detected = detect_engines(temp_root, override_engine="gnina", redetect=True)
        _assert(bool(override_detected["detection_override"]) is True, "engine override should be recorded")
        _assert(str(override_detected["routed_to_engine"]) == "gnina", "override should route to gnina")

        gnina_only_root = Path(tempfile.mkdtemp(prefix="dockforge_engine_detect_gnina_"))
        try:
            bootstrap_project_layout(gnina_only_root, ["gnina"], layout_profile="docking_legacy")
            _write_test_pairlist(gnina_only_root)
            _seed_engine_outputs(gnina_only_root, ("gnina",))
            pipeline = UnifiedPostDockingPipeline(
                project_dir=str(gnina_only_root),
                output_dir=str(post_docking_root(gnina_only_root) / "sessions" / "auto_route_smoke"),
                analysis_mode="auto",
            )
            _assert(pipeline.analysis_mode == "single_engine", "GNINA-only project should auto-route to single_engine")
            _assert(str(pipeline.engine) == "gnina", "GNINA-only project should select gnina as engine")
            _assert(pipeline.run() is True, "GNINA-only auto-routed pipeline should run")
            best_df = pd.read_csv(post_docking_root(gnina_only_root) / "sessions" / "auto_route_smoke" / "gnina" / "best_poses.csv")
            _assert("cnn_confidence" in best_df.columns, "GNINA solo outputs should include cnn_confidence")
            _assert("ligand_efficiency" in best_df.columns, "GNINA solo outputs should include ligand_efficiency")
            _assert(
                str(best_df.loc[0, "primary_ranking_label"]) == "cnn_affinity",
                "GNINA solo outputs should rank by cnn_affinity",
            )
        finally:
            shutil.rmtree(gnina_only_root, ignore_errors=True)
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_engine_detection_autodock4_parity_contract() -> None:
    autodock_root = Path(tempfile.mkdtemp(prefix="dockforge_engine_detect_autodock4_"))
    try:
        bootstrap_project_layout(autodock_root, ["autodock4"], layout_profile="docking_legacy")
        _write_test_pairlist(autodock_root)
        _seed_engine_outputs(autodock_root, ("autodock4",))

        parsed = parse_autodock4_dlg(ensure_engine_layout(autodock_root, "autodock4")["poses"] / "R1_site_1_L1.dlg")
        _assert(not parsed.empty, "AutoDock4 parser should extract pose rows from DLG")
        _assert(
            {"pose", "autodock4_affinity"}.issubset(parsed.columns),
            "AutoDock4 parser output should include pose and autodock4_affinity columns",
        )

        detected = detect_engines(autodock_root, redetect=True)
        _assert(
            str(detected.get("routing_decision", "")) == "single_engine",
            "AutoDock4-only project should route to single_engine",
        )
        _assert(
            str(detected.get("routed_to_engine", "")) == "autodock4",
            "AutoDock4-only project should route to autodock4",
        )
        _assert(
            "autodock4" in list(detected.get("valid_engines", [])),
            "AutoDock4 should be included in valid_engines",
        )
        payload = dict((detected.get("engines") or {}).get("autodock4") or {})
        _assert(
            str(payload.get("pose_format", "")) == "dlg",
            "AutoDock4 detection payload should declare DLG pose format",
        )
        _assert(
            int(payload.get("valid_score_rows", 0)) > 0,
            "AutoDock4 detection payload should report positive valid_score_rows",
        )
        _assert(
            list(payload.get("score_columns", [])) == ["autodock4_affinity"],
            "AutoDock4 detection payload should expose autodock4_affinity score column",
        )

        mixed_root = Path(tempfile.mkdtemp(prefix="dockforge_engine_detect_mixed4_"))
        try:
            bootstrap_project_layout(mixed_root, ["gnina", "vina", "smina", "autodock4"], layout_profile="docking_legacy")
            _write_test_pairlist(mixed_root)
            _seed_engine_outputs(mixed_root, ("gnina", "vina", "smina", "autodock4"))
            mixed = detect_engines(mixed_root, redetect=True)
            _assert(
                str(mixed.get("routing_decision", "")) == "multi_engine",
                "four-engine project should route to multi_engine",
            )
            _assert(
                {"gnina", "vina", "smina", "autodock4"}.issubset(set(mixed.get("valid_engines", []))),
                "engine detection should include all four engines as valid when outputs exist",
            )
        finally:
            shutil.rmtree(mixed_root, ignore_errors=True)
    finally:
        shutil.rmtree(autodock_root, ignore_errors=True)


def _smoke_engine_hpc_adapter_parity_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_engine_hpc_adapter_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina", "autodock4"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina", "autodock4"))

        detectors = {
            "gnina": (detect_gnina_layout, ".sdf"),
            "vina": (detect_vina_layout, ".pdbqt"),
            "smina": (detect_smina_layout, ".pdbqt"),
            "autodock4": (detect_autodock4_layout, ".dlg"),
        }

        for engine, (detector, suffix) in detectors.items():
            layout = detector(temp_root)
            generic = detect_engine_layout(temp_root, engine)

            pose_folder = layout.get("pose_folder") or layout.get("sdf_folder")
            _assert(isinstance(pose_folder, Path) and pose_folder.is_dir(), f"{engine} adapter should resolve pose folder")
            _assert(any(pose_folder.glob(f"*{suffix}")), f"{engine} adapter should resolve pose folder with {suffix} files")

            log_folder = layout.get("log_folder")
            _assert(isinstance(log_folder, Path) and log_folder.is_dir(), f"{engine} adapter should resolve log folder")

            receptors_folder = layout.get("receptors_folder")
            _assert(isinstance(receptors_folder, Path) and receptors_folder.is_dir(), f"{engine} adapter should resolve receptors folder")

            pairlist_file = layout.get("pairlist_file")
            _assert(isinstance(pairlist_file, Path) and pairlist_file.is_file(), f"{engine} adapter should resolve pairlist file")

            _assert(
                str(layout.get("layout", "")) in {"hpc", "local"},
                f"{engine} adapter should report local/hpc layout label",
            )
            _assert(
                str(generic.get("engine", "")) == engine,
                f"generic engine adapter should label the payload with engine={engine}",
            )
            _assert(
                generic.get("pose_folder") == pose_folder and generic.get("log_folder") == log_folder,
                f"generic and wrapper adapter paths should agree for {engine}",
            )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_run_tracking_schema_consistency_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_run_tracking_schema_"))
    expected_step_columns = [
        "index",
        "step",
        "required",
        "status",
        "started_at",
        "ended_at",
        "details",
        "error",
    ]
    expected_output_columns = [
        "relative_path",
        "category",
        "extension",
        "size_bytes",
        "modified_utc",
    ]
    expected_stage_contract_keys = {
        "requested_run_rmsd",
        "requested_run_visualizations",
        "effective_run_rmsd",
        "effective_run_visualizations",
        "requested_run_prolif",
        "requested_run_ligplot",
        "effective_run_prolif",
        "effective_run_ligplot",
    }
    expected_manifest_core_keys = {
        "pipeline",
        "input",
        "output_dir",
        "run_started_at",
        "run_completed_at",
        "run_status",
        "error_message",
        "stage_contract",
        "step_status_file",
        "outputs_index_file",
        "steps_total",
        "step_status_counts",
        "outputs_count",
        "optional_features",
    }

    def _capture_schema(run_tracking_dir: Path) -> dict:
        manifest = json.loads((run_tracking_dir / "run_manifest.json").read_text(encoding="utf-8"))
        step_df = pd.read_csv(run_tracking_dir / "step_status.csv")
        outputs_df = pd.read_csv(run_tracking_dir / "outputs_index.csv")
        return {
            "manifest_keys": set(manifest.keys()),
            "step_columns": list(step_df.columns.astype(str)),
            "outputs_columns": list(outputs_df.columns.astype(str)),
            "stage_contract_keys": set((manifest.get("stage_contract") or {}).keys()),
        }

    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina", "autodock4"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina", "autodock4"))

        baseline_schema = None
        for engine in ("gnina", "vina", "smina", "autodock4"):
            output_dir = post_docking_root(temp_root) / "sessions" / f"schema_{engine}"
            pipeline = UnifiedPostDockingPipeline(
                project_dir=str(temp_root),
                output_dir=str(output_dir),
                analysis_mode="single_engine",
                engine=engine,
                analysis_scope="comparison_only",
            )
            _stub_dag_optional_interactions(pipeline)
            _assert(pipeline.run() is True, f"single-engine run should succeed for {engine}")

            schema = _capture_schema(output_dir / "run_tracking")
            _assert(schema["step_columns"] == expected_step_columns, f"{engine} step_status schema mismatch")
            _assert(schema["outputs_columns"] == expected_output_columns, f"{engine} outputs_index schema mismatch")
            _assert(
                schema["stage_contract_keys"] == expected_stage_contract_keys,
                f"{engine} stage_contract schema mismatch",
            )
            _assert(
                expected_manifest_core_keys.issubset(schema["manifest_keys"]),
                f"{engine} run manifest should include canonical run-tracking keys",
            )
            if baseline_schema is None:
                baseline_schema = schema
            else:
                _assert(
                    schema["manifest_keys"] == baseline_schema["manifest_keys"],
                    f"{engine} run manifest keys should match baseline schema",
                )

        simplified_root = temp_root / "simplified_schema_fixture"
        sdf_dir = simplified_root / "sdf"
        log_dir = simplified_root / "logs"
        receptors_dir = simplified_root / "receptors"
        output_dir = simplified_root / "output"
        sdf_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)
        receptors_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "reports").mkdir(parents=True, exist_ok=True)
        (output_dir / "reports" / "dummy.csv").write_text("ok\n", encoding="utf-8")

        simplified = SimplifiedPostDockingPipeline(
            sdf_folder=str(sdf_dir),
            log_folder=str(log_dir),
            receptors_folder=str(receptors_dir),
            output_dir=str(output_dir),
            run_rmsd=True,
            run_visualizations=True,
        )
        simplified._write_run_tracking_manifest(
            run_started_at=datetime.now(timezone.utc).isoformat(),
            run_status="completed",
            error_message="",
            step_states=[
                {
                    "index": 1,
                    "step": "Schema smoke step",
                    "required": True,
                    "status": "completed",
                    "started_at": datetime.now(timezone.utc).isoformat(),
                    "ended_at": datetime.now(timezone.utc).isoformat(),
                    "details": "schema_smoke",
                    "error": "",
                }
            ],
        )
        simplified_schema = _capture_schema(output_dir / "run_tracking")
        _assert(
            simplified_schema["step_columns"] == expected_step_columns,
            "simplified step_status schema should match canonical run_tracking schema",
        )
        _assert(
            simplified_schema["outputs_columns"] == expected_output_columns,
            "simplified outputs_index schema should match canonical run_tracking schema",
        )
        _assert(
            simplified_schema["stage_contract_keys"] == expected_stage_contract_keys,
            "simplified stage_contract schema should match canonical run_tracking schema",
        )
        _assert(
            expected_manifest_core_keys.issubset(simplified_schema["manifest_keys"]),
            "simplified run manifest should include canonical run-tracking keys",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_optional_feature_taxonomy_contract() -> None:
    import post_docking_analysis.simplified_pipeline_impl as simplified_impl_module

    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_optional_feature_taxonomy_"))
    original_py3dmol_available = bool(simplified_impl_module.PY3DMOL_AVAILABLE)
    try:
        # --- Simplified manifest taxonomy branch checks ---
        simplified_root = temp_root / "simplified_taxonomy"
        sdf_dir = simplified_root / "sdf"
        log_dir = simplified_root / "logs"
        receptors_dir = simplified_root / "receptors"
        output_dir = simplified_root / "output"
        (output_dir / "complexes").mkdir(parents=True, exist_ok=True)
        (output_dir / "complexes" / "dummy_complex.pdb").write_text("ATOM      1  C   UNK A   1\nEND\n", encoding="utf-8")
        sdf_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)
        receptors_dir.mkdir(parents=True, exist_ok=True)

        simplified_impl_module.PY3DMOL_AVAILABLE = False
        simplified = SimplifiedPostDockingPipeline(
            sdf_folder=str(sdf_dir),
            log_folder=str(log_dir),
            receptors_folder=str(receptors_dir),
            output_dir=str(output_dir),
            enable_poseview=False,
            run_rmsd=True,
            run_visualizations=True,
        )
        simplified._write_run_tracking_manifest(
            run_started_at=datetime.now(timezone.utc).isoformat(),
            run_status="completed",
            error_message="",
            step_states=[
                {
                    "index": 1,
                    "step": "Generate PoseView diagrams",
                    "required": False,
                    "status": "completed",
                    "started_at": datetime.now(timezone.utc).isoformat(),
                    "ended_at": datetime.now(timezone.utc).isoformat(),
                    "details": "disabled branch smoke",
                    "error": "",
                },
                {
                    "index": 2,
                    "step": "Generate PandaMap analysis",
                    "required": False,
                    "status": "skipped_or_unavailable",
                    "started_at": datetime.now(timezone.utc).isoformat(),
                    "ended_at": datetime.now(timezone.utc).isoformat(),
                    "details": "failed branch smoke",
                    "error": "",
                },
                {
                    "index": 3,
                    "step": "Generate py3Dmol visualizations",
                    "required": False,
                    "status": "skipped_or_unavailable",
                    "started_at": datetime.now(timezone.utc).isoformat(),
                    "ended_at": datetime.now(timezone.utc).isoformat(),
                    "details": "dependency branch smoke",
                    "error": "",
                },
            ],
        )
        simplified_manifest = json.loads((output_dir / "run_tracking" / "run_manifest.json").read_text(encoding="utf-8"))
        simplified_optional = dict(simplified_manifest.get("optional_features") or {})
        _assert(
            str(((simplified_optional.get("poseview") or {}).get("status") or "")) == "skipped_disabled",
            "simplified optional feature taxonomy should mark disabled PoseView as skipped_disabled",
        )
        _assert(
            str(((simplified_optional.get("py3dmol") or {}).get("status") or "")) == "skipped_missing_dependency",
            "simplified optional feature taxonomy should mark missing py3Dmol as skipped_missing_dependency",
        )
        _assert(
            str(((simplified_optional.get("pandamap") or {}).get("status") or "")) == "failed_error",
            "simplified optional feature taxonomy should mark unsuccessful PandaMap stage as failed_error",
        )

        # --- Multi-engine manifest taxonomy branch checks ---
        bootstrap_project_layout(temp_root, ["gnina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina",))
        multi_output = post_docking_root(temp_root) / "sessions" / "optional_taxonomy_multi"
        multi = MultiEngineAnalysisPipeline(
            project_dir=str(temp_root),
            output_dir=str(multi_output),
            analysis_mode="comparative_all_engines",
            analysis_scope="report_only",
        )
        multi.dag_execution_report = {
            "nodes": {
                "visualizations": {"status": "completed", "details": "figures_generated=2"},
                "prolif": {"status": "completed_with_warnings", "details": "ProLIF not installed"},
                "poseview": {"status": "completed_with_warnings", "details": "PoseView stage disabled"},
                "pandamap": {"status": "failed", "details": "network timeout"},
                "pymol": {"status": "failed", "details": "PyMOL crash"},
            }
        }
        multi._write_run_tracking_artifacts(
            run_started_at=datetime.now(timezone.utc).isoformat(),
            run_status="completed",
            error_message="",
            step_states=[
                {
                    "step": "load_or_build_scores",
                    "status": "completed",
                    "started_at": datetime.now(timezone.utc).isoformat(),
                    "ended_at": datetime.now(timezone.utc).isoformat(),
                    "details": "taxonomy smoke",
                }
            ],
        )
        multi_manifest = json.loads((multi_output / "run_tracking" / "run_manifest.json").read_text(encoding="utf-8"))
        multi_optional = dict(multi_manifest.get("optional_features") or {})
        _assert(
            str(((multi_optional.get("visualizations") or {}).get("status") or "")) == "completed",
            "multi-engine optional feature taxonomy should keep successful visualizations as completed",
        )
        _assert(
            str(((multi_optional.get("prolif") or {}).get("status") or "")) == "skipped_missing_dependency",
            "multi-engine optional feature taxonomy should mark unavailable ProLIF as skipped_missing_dependency",
        )
        _assert(
            str(((multi_optional.get("poseview") or {}).get("status") or "")) == "skipped_disabled",
            "multi-engine optional feature taxonomy should mark disabled PoseView as skipped_disabled",
        )
        _assert(
            str(((multi_optional.get("pandamap") or {}).get("status") or "")) == "failed_error",
            "multi-engine optional feature taxonomy should mark failed PandaMap as failed_error",
        )
    finally:
        simplified_impl_module.PY3DMOL_AVAILABLE = original_py3dmol_available
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_simplified_pipeline_interaction_contract_stages() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_simplified_interaction_contract_"))
    try:
        sdf_dir = temp_root / "sdf"
        log_dir = temp_root / "logs"
        receptors_dir = temp_root / "receptors"
        output_dir = temp_root / "output"
        for directory in (sdf_dir, log_dir, receptors_dir, output_dir):
            directory.mkdir(parents=True, exist_ok=True)

        pipeline = SimplifiedPostDockingPipeline(
            sdf_folder=str(sdf_dir),
            log_folder=str(log_dir),
            receptors_folder=str(receptors_dir),
            output_dir=str(output_dir),
            run_rmsd=True,
            run_visualizations=True,
        )

        for method_name in (
            "_find_input_files",
            "_generate_scores_csv",
            "_match_poses_to_receptors",
            "_receptor_chain_integrity_qc_and_repair",
            "_create_complexes",
            "_analyze_binding_affinity",
            "_analyze_polypharmacology",
            "_analyze_rmsd",
            "_extract_poses",
            "_generate_reports",
            "_generate_visualizations",
            "_generate_plip_interaction_outputs",
            "_generate_prolif_interaction_maps",
            "_generate_ligplot_diagrams",
            "_normalize_interaction_tables",
            "_build_source_of_truth_interaction_analytics",
            "_generate_poseview_diagrams",
            "_generate_pandamap_analysis",
            "_generate_py3dmol_visualizations",
            "_organize_output_artifacts",
        ):
            setattr(pipeline, method_name, lambda _name=method_name: True)

        _assert(pipeline.run() is True, "simplified pipeline should complete with stubbed interaction contract stages")

        step_df = pd.read_csv(output_dir / "run_tracking" / "step_status.csv")
        step_names = set(step_df.get("step", pd.Series(dtype="object")).astype(str).tolist())
        for required_step in (
            "Receptor chain-integrity QC + repair",
            "Generate PLIP interaction outputs",
            "Normalize layered PLIP/ProLIF tables",
            "Build source-of-truth interaction analytics",
        ):
            _assert(required_step in step_names, f"missing simplified contract stage in step log: {required_step}")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_pandamap_publication_artifact_contract() -> None:
    import post_docking_analysis.publication_pandamap as publication_module

    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_pandamap_artifact_contract_"))
    original_runner = publication_module.PandaMapRunner
    try:
        complexes_dir = temp_root / "complexes"
        output_dir = temp_root / "pandamap"
        complexes_dir.mkdir(parents=True, exist_ok=True)
        (complexes_dir / "COX2_ligand_SRL_A_1.pdb").write_text(
            "HETATM    1  C1  SRL A 401      11.111  11.111  11.111  1.00  0.00           C\nEND\n",
            encoding="utf-8",
        )

        class _DummyPandaMapRunner:
            def __init__(self, conda_env: str = "pandamap", logger=None):
                self.conda_env = conda_env
                self.logger = logger

            def detect_cli_mode(self) -> str:
                return "single"

            def get_supported_options(self) -> set[str]:
                return {
                    "--output",
                    "--dpi",
                    "--title",
                    "--ligand",
                    "--report",
                    "--report-file",
                    "--3d",
                    "--3d-output",
                    "--deltaG",
                }

            def describe_capabilities(self) -> Dict[str, object]:
                return {
                    "conda_env": self.conda_env,
                    "cli_mode": "single",
                    "supported_options": sorted(self.get_supported_options()),
                    "option_matrix": {
                        "text_report": True,
                        "interactive_3d": True,
                        "delta_g_estimation": True,
                    },
                }

            def run(self, args, timeout=300, cwd=None):
                if "--output" in args:
                    output_path = Path(args[args.index("--output") + 1])
                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    output_path.write_text("ok\n", encoding="utf-8")
                if "--3d-output" in args:
                    output_path = Path(args[args.index("--3d-output") + 1])
                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    output_path.write_text("<html></html>\n", encoding="utf-8")
                if "--report-file" in args:
                    report_path = Path(args[args.index("--report-file") + 1])
                    report_path.parent.mkdir(parents=True, exist_ok=True)
                    report_path.write_text("RESIDUE INTERACTIONS\n", encoding="utf-8")

                class _Result:
                    returncode = 0
                    stdout = ""
                    stderr = ""

                return _Result()

        publication_module.PandaMapRunner = _DummyPandaMapRunner

        summary = publication_module.run_publication_pandamap_analysis(
            complexes_dir=complexes_dir,
            output_dir=output_dir,
            ligand_name="UNK",
            conda_env="pandamap",
            config={"overwrite": True},
            max_complexes=1,
        )
        _assert(bool(summary), "PandaMap publication contract smoke should return summary payload")
        capability_manifest = Path(str(summary.get("capability_manifest", "")))
        _assert(capability_manifest.exists(), "PandaMap capability manifest should be written")

        results_csv = output_dir / "pandamap_analysis_results.csv"
        _assert(results_csv.exists(), "PandaMap publication stage should emit analysis results CSV")
        results_df = pd.read_csv(results_csv)
        for column in ["ligand_resolution_source", "figure_title", "report_file"]:
            _assert(column in results_df.columns, f"PandaMap results CSV should include {column}")
        _assert(
            bool(results_df["report_file"].astype(str).str.len().gt(0).any()),
            "PandaMap results CSV should include at least one non-empty report_file path",
        )
    finally:
        publication_module.PandaMapRunner = original_runner
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_pandamap_quality_gate_contract() -> None:
    import post_docking_analysis.simplified_pipeline_impl as simplified_impl_module

    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_pandamap_quality_gate_"))
    original_runner = simplified_impl_module.run_publication_pandamap_analysis
    try:
        sdf_dir = temp_root / "sdf"
        log_dir = temp_root / "logs"
        receptors_dir = temp_root / "receptors"
        output_dir = temp_root / "output"
        for directory in [sdf_dir, log_dir, receptors_dir, output_dir / "complexes"]:
            directory.mkdir(parents=True, exist_ok=True)
        (output_dir / "complexes" / "dummy_complex.pdb").write_text(
            "HETATM    1  C1  UNK A   1      11.111  11.111  11.111  1.00  0.00           C\nEND\n",
            encoding="utf-8",
        )

        def _fake_failed_quality_run(**kwargs):
            root = Path(str(kwargs.get("output_dir", "")))
            maps_2d = root / "maps_2d"
            maps_2d.mkdir(parents=True, exist_ok=True)
            (maps_2d / "dummy_complex.png").write_text("fake\n", encoding="utf-8")
            return {
                "generated_2d_maps": 1,
                "generated_3d_visualizations": 0,
                "analysis_results": [
                    {
                        "complex": "dummy_complex",
                        "ligand_resolution_source": "generic_fallback",
                        "report_file": "",
                        "2d_maps": ["png"],
                    }
                ],
                "pandamap_runtime_capabilities": {
                    "option_matrix": {"text_report": True}
                },
            }

        simplified_impl_module.run_publication_pandamap_analysis = _fake_failed_quality_run
        simplified = SimplifiedPostDockingPipeline(
            sdf_folder=str(sdf_dir),
            log_folder=str(log_dir),
            receptors_folder=str(receptors_dir),
            output_dir=str(output_dir),
            run_rmsd=True,
            run_visualizations=True,
        )
        failed = simplified._generate_pandamap_analysis()
        _assert(failed is False, "PandaMap stage should fail when quality gate checks fail")
        failed_summary = simplified.results.get("pandamap_summary", {})
        qc_failed = dict(failed_summary.get("quality_control") or {})
        _assert(
            str(qc_failed.get("status", "")) == "failed",
            "PandaMap quality gate summary should mark failed status",
        )
        failure_reasons = [str(item) for item in (qc_failed.get("failure_reasons") or [])]
        _assert(
            any("non_generic_ligand_ratio" in reason for reason in failure_reasons),
            "PandaMap quality gate should report ligand provenance failure",
        )
        _assert(
            any("text reports" in reason for reason in failure_reasons),
            "PandaMap quality gate should report missing report-file outputs when reports are supported",
        )

        def _fake_pass_quality_run(**kwargs):
            root = Path(str(kwargs.get("output_dir", "")))
            maps_2d = root / "maps_2d"
            maps_2d.mkdir(parents=True, exist_ok=True)
            (maps_2d / "dummy_complex.png").write_text("fake\n", encoding="utf-8")
            report_file = root / "maps_2d" / "dummy_complex.txt"
            report_file.write_text("ok\n", encoding="utf-8")
            return {
                "generated_2d_maps": 1,
                "generated_3d_visualizations": 0,
                "analysis_results": [
                    {
                        "complex": "dummy_complex",
                        "ligand_resolution_source": "detected_from_pdb_hetatm",
                        "report_file": str(report_file),
                        "2d_maps": ["png"],
                    }
                ],
                "pandamap_runtime_capabilities": {
                    "option_matrix": {"text_report": True}
                },
            }

        simplified_impl_module.run_publication_pandamap_analysis = _fake_pass_quality_run
        passed = simplified._generate_pandamap_analysis()
        _assert(passed is True, "PandaMap stage should pass when quality gate checks pass")
        passed_summary = simplified.results.get("pandamap_summary", {})
        qc_passed = dict(passed_summary.get("quality_control") or {})
        _assert(
            str(qc_passed.get("status", "")) == "passed",
            "PandaMap quality gate summary should mark passed status on healthy outputs",
        )
    finally:
        simplified_impl_module.run_publication_pandamap_analysis = original_runner
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_engine_detection_cache_override_contract() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_engine_detect_cache_override_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina"))

        detect_engines(temp_root, redetect=True)
        cached_override = detect_engines(temp_root, override_engine="vina", redetect=False)
        _assert(bool(cached_override.get("cache_hit", False)) is True, "second detection should use cache")
        _assert(
            bool(cached_override.get("detection_override", False)) is True,
            "cached detection should honor override_engine",
        )
        _assert(
            str(cached_override.get("routed_to_engine", "")) == "vina",
            "cached detection should route to override engine",
        )
        _assert(
            str(cached_override.get("routing_decision", "")) == "single_engine",
            "cached override should force single_engine routing",
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_rmsd_scope_normalizer_contract() -> None:
    scopes = normalize_simplified_rmsd_scopes("per_protein,global")
    _assert(
        scopes == ("per_protein", "global"),
        f"expected per_protein/global scopes, got: {scopes}",
    )
    all_scopes = normalize_simplified_rmsd_scopes("all")
    _assert(
        all_scopes == ("per_complex", "per_protein", "global"),
        f"expected full RMSD scope expansion, got: {all_scopes}",
    )


def _smoke_strict_consensus_single_engine_warning_contract() -> None:
    frame = pd.DataFrame(
        [
            {"engine": "gnina", "protein": "P1", "ligand": "L1", "site_id": "site_1", "tag": "P1_site_1_A", "affinity_kcal_mol": -8.0},
            {"engine": "gnina", "protein": "P1", "ligand": "L2", "site_id": "site_1", "tag": "P1_site_1_B", "affinity_kcal_mol": -7.0},
        ]
    )
    captured: list[str] = []
    consensus_logger = logging.getLogger("post_docking_analysis.consensus")
    original_level = consensus_logger.level

    class _Capture(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            captured.append(str(record.getMessage()))

    handler = _Capture(level=logging.WARNING)
    consensus_logger.addHandler(handler)
    consensus_logger.setLevel(logging.WARNING)
    try:
        ranked = build_consensus_rankings(
            frame,
            consensus_mode="strict_consensus",
            favorite_engine="gnina",
            normalization_method="per_engine_rank",
        )
    finally:
        consensus_logger.removeHandler(handler)
        consensus_logger.setLevel(original_level)

    _assert(len(ranked) == len(frame), "strict_consensus should retain rows in single-engine mode")
    _assert(
        any("strict_consensus requested in single-engine mode" in message for message in captured),
        "single-engine strict_consensus should emit a warning",
    )


def _smoke_engine_scope_filter_foundation() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_engine_scope_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina"))
        manifest = load_manifest(temp_root)
        manifest["engine_presets"] = {
            "cnn_only": ["gnina"],
            "vina_pair": ["gnina", "vina"],
        }
        from docking.project_layout import save_manifest
        save_manifest(temp_root, manifest)

        pair_pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(post_docking_root(temp_root) / "sessions" / "scope_pair"),
            analysis_mode="auto",
            engines="gnina,vina",
        )
        _assert(pair_pipeline.analysis_mode == "comparative_all_engines", "two-engine scope should stay multi-engine")
        _assert(pair_pipeline.engines_in_scope == ["gnina", "vina"], "scoped engines should preserve requested subset")
        _assert(
            any(item.get("engine") == "smina" for item in pair_pipeline.excluded_engines),
            "excluded engines should be tracked in scoped session metadata",
        )
        scores = pair_pipeline._load_or_build_scores()
        _assert(set(scores["engine"].astype(str)) == {"gnina", "vina"}, "scoped score load should exclude smina")

        solo_pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(post_docking_root(temp_root) / "sessions" / "scope_solo"),
            analysis_mode="auto",
            engines="gnina",
        )
        _assert(solo_pipeline.analysis_mode == "single_engine", "single scoped engine should route to solo mode")
        _assert(str(solo_pipeline.engine) == "gnina", "single scoped engine should become active engine")

        preset_pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(post_docking_root(temp_root) / "sessions" / "scope_preset"),
            analysis_mode="auto",
            engine_preset="cnn_only",
        )
        _assert(preset_pipeline.analysis_mode == "single_engine", "preset with one engine should route to solo mode")
        _assert(preset_pipeline.engines_in_scope == ["gnina"], "preset should resolve engines_in_scope")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_engine_scope_filter_acceptance() -> None:
    temp_root = Path(tempfile.mkdtemp(prefix="dockforge_engine_scope_accept_"))
    try:
        bootstrap_project_layout(temp_root, ["gnina", "vina", "smina"], layout_profile="docking_legacy")
        _write_test_pairlist(temp_root)
        _seed_engine_outputs(temp_root, ("gnina", "vina", "smina"))
        manifest = load_manifest(temp_root)
        manifest["engine_presets"] = {"cnn_only": ["gnina"]}
        from docking.project_layout import save_manifest
        save_manifest(temp_root, manifest)

        smina_pose = ensure_engine_layout(temp_root, "smina")["poses"] / "R1_site_1_L1.pdbqt"
        before_stat = smina_pose.stat().st_mtime_ns

        pair_output = post_docking_root(temp_root) / "sessions" / "scope_pair_accept"
        pair_pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(pair_output),
            analysis_mode="comparative_all_engines",
            engines="gnina,vina",
        )
        _assert(pair_pipeline.run() is True, "two-engine scoped comparative run should succeed")
        combined = pd.read_csv(pair_output / "reports" / "combined_engine_scores.csv")
        _assert(set(combined["engine"].astype(str)) == {"gnina", "vina"}, "excluded engine should not appear in scoped score outputs")
        _assert(set(combined["engines_in_scope"].astype(str)) == {"gnina,vina"}, "scoped outputs should record engine scope")
        top_pose = pd.read_csv(pair_output / "top_pose_ligand_performance" / "top_pose_per_ligand_global.csv")
        _assert(int(top_pose["engine_support_count"].iloc[0]) == 2, "engine_support_count should equal scoped engine count")
        run_manifest = json.loads((pair_output / "run_tracking" / "run_manifest.json").read_text(encoding="utf-8"))
        _assert(run_manifest.get("engines_in_scope") == ["gnina", "vina"], "run manifest should capture scoped engines")
        detection_report = json.loads((temp_root / "4-Working" / "metadata" / "engine_detection_report.json").read_text(encoding="utf-8"))
        _assert(
            bool(((detection_report.get("engines") or {}).get("smina") or {}).get("excluded_this_session")) is True,
            "persisted detection report should mark excluded engines as excluded_this_session",
        )
        after_stat = smina_pose.stat().st_mtime_ns
        _assert(before_stat == after_stat, "excluded engine files should remain untouched after scoped runs")

        solo_output = post_docking_root(temp_root) / "sessions" / "scope_solo_accept"
        solo_pipeline = UnifiedPostDockingPipeline(
            project_dir=str(temp_root),
            output_dir=str(solo_output),
            analysis_mode="auto",
            engine_preset="cnn_only",
        )
        _assert(solo_pipeline.run() is True, "preset-scoped solo run should succeed")
        start_here = (temp_root / "7-Reports" / "START_HERE.md").read_text(encoding="utf-8")
        _assert("Scope Notice" in start_here, "START_HERE should note cross-session scope discrepancy")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _smoke_engine_scope_picker_zero_guard() -> None:
    import workflow.interactive as interactive_mod

    detection = {
        "valid_engines": ["gnina", "vina"],
        "engines": {
            "gnina": {"coverage_pct": 100.0},
            "vina": {"coverage_pct": 80.0},
        },
    }
    manifest = {}
    responses = iter([[], ["gnina"]])
    original_checkbox = interactive_mod._checkbox
    try:
        interactive_mod._checkbox = lambda *args, **kwargs: next(responses)
        selected, preset = _analysis_engine_scope_picker(
            _FakeQuestionary(),
            manifest,
            detection,
            preselected=["gnina", "vina"],
            preset_name="",
        )
        _assert(selected == ["gnina"], "picker should re-prompt until at least one engine remains selected")
        _assert(preset == "", "manual picker path should not return a preset")
    finally:
        interactive_mod._checkbox = original_checkbox


def _smoke_no_valid_engine_abort_and_solo_semantics() -> None:
    empty_root = Path(tempfile.mkdtemp(prefix="dockforge_no_valid_engines_"))
    try:
        bootstrap_project_layout(empty_root, ["gnina"], layout_profile="docking_legacy")
        _write_test_pairlist(empty_root)
        try:
            UnifiedPostDockingPipeline(
                project_dir=str(empty_root),
                output_dir=str(post_docking_root(empty_root) / "sessions" / "empty"),
                analysis_mode="auto",
            )
        except ValueError as exc:
            _assert("No valid docking-engine outputs were detected" in str(exc), "no-valid-engines should raise actionable error")
        else:
            raise RuntimeError("no-valid-engine project should abort")
    finally:
        shutil.rmtree(empty_root, ignore_errors=True)

    smina_root = Path(tempfile.mkdtemp(prefix="dockforge_smina_solo_"))
    try:
        bootstrap_project_layout(smina_root, ["smina"], layout_profile="docking_legacy")
        _write_test_pairlist(smina_root)
        _seed_engine_outputs(smina_root, ("smina",))
        output_dir = post_docking_root(smina_root) / "sessions" / "smina_solo"
        pipeline = UnifiedPostDockingPipeline(
            project_dir=str(smina_root),
            output_dir=str(output_dir),
            analysis_mode="auto",
        )
        _assert(pipeline.analysis_mode == "single_engine", "Smina-only project should auto-route to single_engine")
        _assert(pipeline.run() is True, "Smina solo pipeline should run")
        _assert(
            not (output_dir / "reports").exists(),
            "Solo mode should not emit comparative multi-engine report bundles",
        )
        best_df = pd.read_csv(output_dir / "smina" / "best_poses.csv")
        _assert(
            "pose_diversity_status" in best_df.columns and str(best_df.loc[0, "pose_diversity_status"]) == "not_available",
            "Smina solo should mark RMSD-dependent diversity as not_available",
        )
        manifest = load_manifest(smina_root)
        _assert(
            isinstance(manifest.get("engine_settings", {}).get("smina", {}).get("smina_scoring_function", {}), dict),
            "Smina solo should persist scoring-function provenance into manifest",
        )
    finally:
        shutil.rmtree(smina_root, ignore_errors=True)

    vina_root = Path(tempfile.mkdtemp(prefix="dockforge_vina_solo_"))
    try:
        bootstrap_project_layout(vina_root, ["vina"], layout_profile="docking_legacy")
        _write_test_pairlist(vina_root)
        _seed_engine_outputs(vina_root, ("vina",))
        output_dir = post_docking_root(vina_root) / "sessions" / "vina_solo"
        pipeline = UnifiedPostDockingPipeline(
            project_dir=str(vina_root),
            output_dir=str(output_dir),
            analysis_mode="auto",
        )
        _assert(pipeline.analysis_mode == "single_engine", "Vina-only project should auto-route to single_engine")
        _assert(pipeline.run() is True, "Vina solo pipeline should run")
        best_df = pd.read_csv(output_dir / "vina" / "best_poses.csv")
        _assert("pose_diversity_metric" in best_df.columns, "Vina solo should expose pose_diversity_metric")
        summary_text = (output_dir / "vina" / "summary.txt").read_text(encoding="utf-8")
        _assert("docking collapse" not in summary_text.lower(), "a single best pose's zero self-RMSD is not evidence of docking collapse")
    finally:
        shutil.rmtree(vina_root, ignore_errors=True)


def _smoke_cli_interactive_fallback_and_gnina_complex_export() -> None:
    fallback_root = Path(tempfile.mkdtemp(prefix="dockforge_cli_fallback_"))
    try:
        bootstrap_project_layout(fallback_root, ["gnina", "vina"], layout_profile="docking_legacy")
        responses = iter(["2"])
        original_input = builtins.input
        original_isatty = sys.stdin.isatty
        try:
            builtins.input = lambda _prompt="": next(responses)
            sys.stdin.isatty = lambda: True
            selected = _prompt_cli_engine_fallback(fallback_root)
            _assert(selected == "vina", "CLI interactive fallback should return the selected declared engine")
        finally:
            builtins.input = original_input
            sys.stdin.isatty = original_isatty
    finally:
        shutil.rmtree(fallback_root, ignore_errors=True)

    gnina_root = Path(tempfile.mkdtemp(prefix="dockforge_gnina_complex_export_"))
    try:
        bootstrap_project_layout(gnina_root, ["gnina"], layout_profile="docking_legacy")
        _write_test_pairlist(gnina_root)
        _seed_engine_outputs(gnina_root, ("gnina",))
        output_dir = post_docking_root(gnina_root) / "sessions" / "gnina_export"
        pipeline = UnifiedPostDockingPipeline(
            project_dir=str(gnina_root),
            output_dir=str(output_dir),
            analysis_mode="auto",
        )
        _assert(pipeline.run() is True, "GNINA solo pipeline should run for complex export confirmation")
        summary_text = (output_dir / "gnina" / "summary.txt").read_text(encoding="utf-8")
        _assert("Complex generation parser: sdf_multiconformer" in summary_text, "GNINA solo summary should document SDF parser path")
        _assert("Complex exports written: 1" in summary_text, "GNINA solo summary should report exported complexes")
        complex_pdb = output_dir / "gnina" / "complex_exports" / "best_poses" / "R1_site_1_L1" / "R1_site_1_L1_pose1.pdb"
        _assert(complex_pdb.exists(), "GNINA solo should export a best-pose complex PDB from SDF input")
    finally:
        shutil.rmtree(gnina_root, ignore_errors=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run DockForge smoke checks.")
    parser.add_argument(
        "--skip-all-engines",
        action="store_true",
        help="Skip the full all-engines dry-run smoke and execute only fast contract checks.",
    )
    parser.add_argument(
        "--skip-prep-matrix",
        action="store_true",
        help="Skip the ligand preparation profile matrix smoke checks.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    print("🧪 Running DockForge smoke checks")
    print("=" * 60)
    _smoke_deterministic_run_id_utility()
    print("✅ Deterministic run-id utility checks passed")
    _smoke_feature_flags_and_checkpoint_scaffold()
    print("✅ Feature flags and checkpoint metadata scaffold checks passed")
    _smoke_meta_artifact_writers()
    print("✅ Meta reproducibility artifact writer checks passed")
    _smoke_checkpoint_revise_lineage_integrity()
    print("✅ Checkpoint->revise lineage integrity checks passed")
    _smoke_manifest_completeness_and_deterministic_paths()
    print("✅ Manifest completeness and deterministic path checks passed")
    _smoke_prompt_back_navigation_contract()
    print("✅ Universal back-navigation prompt contract checks passed")
    _smoke_hpc_interactive_qc_bypass_contract()
    print("✅ Interactive HPC QC bypass checks passed")
    _smoke_hpc_interactive_admet_bypass_contract()
    print("✅ Interactive HPC ligand ADMET bypass checks passed")
    _smoke_hpc_profile_yaml_loader_contract()
    print("✅ HPC YAML profile loader checks passed")
    _smoke_interactive_ssh_systems_yaml_contract()
    print("✅ Interactive SSH systems YAML checks passed")
    _smoke_nmrbox_condor_slot_parser_contract()
    print("✅ NMRbox Condor slot parser checks passed")
    _smoke_nmrbox_condor_fit_probe_contract()
    print("✅ NMRbox Condor fit probe checks passed")
    _smoke_analysis_session_progress_file_contract()
    print("✅ Analysis session progress-file contract checks passed")
    _smoke_background_task_non_blocking()
    print("✅ Background task non-blocking checks passed")
    _smoke_artifact_graph_executor()
    _smoke_artifact_graph_cache_hit()
    _smoke_artifact_graph_registration_contract()
    _smoke_artifact_graph_scope_request_contract()
    _smoke_artifact_graph_pipeline_run_contract()
    _smoke_artifact_graph_cli_scope_contract()
    _smoke_artifact_graph_default_runner_migration_contract()
    _smoke_artifact_graph_full_parallel_execution_contract()
    _smoke_artifact_graph_unchanged_rerun_cache_contract()
    _smoke_artifact_graph_parameter_invalidation_contract()
    _smoke_artifact_graph_optional_failure_isolation_contract()
    _smoke_artifact_graph_default_comparison_only_contract()
    _smoke_artifact_graph_gnina_solo_schema_contract()
    _smoke_artifact_graph_all_strategy_schema_contract()
    _smoke_artifact_graph_per_complex_warning_contract()
    _smoke_visualization_suite_core_contracts()
    print("✅ Visualization suite core manifest checks passed")
    _smoke_visualization_suite_single_engine_conditionals()
    print("✅ Visualization suite single-engine conditional checks passed")
    _smoke_visualization_suite_skip_policies_and_no_raise()
    print("✅ Visualization suite skip-policy + no-raise checks passed")
    _smoke_visualization_suite_scientific_guardrails()
    print("✅ Visualization suite scientific-guardrail checks passed")
    _smoke_visualization_suite_start_here_embedding_contract()
    print("✅ Visualization suite START_HERE embedding checks passed")
    _smoke_profile_validator()
    print("✅ Preparation profile validator checks passed")
    _smoke_top_pose_selector_determinism()
    print("✅ Top-pose selector determinism checks passed")
    _smoke_top_pose_consensus_direction_contracts()
    print("✅ Top-pose consensus direction checks passed")
    _smoke_consensus_direction_and_single_engine_qc()
    print("✅ Consensus normalization direction + single-engine QC checks passed")
    _smoke_strict_consensus_single_engine_warning_contract()
    print("✅ strict_consensus single-engine warning contract checks passed")
    _smoke_dag_reference_baselines_wiring_contract()
    print("✅ DAG reference-baseline wiring checks passed")
    _smoke_geometric_consensus_soft_alignment_contract()
    print("✅ Geometric consensus soft-alignment checks passed")
    _smoke_geometric_consensus_hard_mismatch_contract()
    print("✅ Geometric consensus hard-mismatch checks passed")
    _smoke_redocking_multi_pose_sdf_parser_contract()
    print("✅ Redocking multi-pose SDF parser checks passed")
    _smoke_redocking_validation_multi_pose_end_to_end_contract()
    print("✅ Redocking multi-pose end-to-end checks passed")
    _smoke_favorite_engine_scope_guard_contract()
    print("✅ favorite_engine scoped-membership guard checks passed")
    _smoke_favorite_engine_empty_scope_guard_contract()
    print("✅ favorite_engine empty-scope guard checks passed")
    _smoke_rmsd_scope_normalizer_contract()
    print("✅ RMSD scope normalizer checks passed")
    _smoke_consensus_mode_isolation_contract()
    print("✅ Consensus mode-isolation contract checks passed")
    _smoke_pairlist_directional_matching_contract()
    print("✅ Pairlist directional matching contract checks passed")
    _smoke_pose_parser_chain_and_dashboard_contract()
    print("✅ Pose parser chain + dashboard contract newline checks passed")
    _smoke_multiligand_residue_identity_contract()
    print("✅ Multi-ligand residue identity contract checks passed")
    _smoke_vina_cnn_correlation_stats_contract()
    print("✅ Vina/CNN correlation p/q + minimum-N checks passed")
    _smoke_cross_engine_correlation_contract_sparse_and_mixed()
    print("✅ Cross-engine per-protein/global correlation contract checks passed")

    _smoke_prepare_guard()
    print("✅ Invalid preparation profile guard check passed")

    _smoke_prepare_ph_guard()
    print("✅ Preparation pH contract guard checks passed")

    _smoke_ligand_output_contract_validator()
    print("✅ Ligand preparation output contract validator checks passed")

    _smoke_preparation_preflight_artifacts()
    print("✅ Preparation preflight artifact contract checks passed")

    _smoke_execution_environment_and_schema()
    print("✅ Execution environment and parameter schema checks passed")

    _smoke_docking_preflight()
    print("✅ Docking preflight checks passed")

    _smoke_docking_preflight_receptor_qc_gate()
    print("✅ Receptor QC preflight gate checks passed")

    _smoke_hpc_deploy_receptor_qc_gate()
    print("✅ HPC deploy receptor QC gate checks passed")

    _smoke_post_docking_allscore_contracts()
    print("✅ Post-docking allscore contract checks passed")

    _smoke_biology_unresolved_mapping_reporting()
    print("✅ Biology unresolved mapping reporting checks passed")

    _smoke_plugin_system_retired_contract()
    print("✅ Plugin system retirement checks passed")

    _smoke_complex_export_index_contract()
    print("✅ Complex export indexing contract checks passed")
    _smoke_extension_stubs()
    print("✅ Future extension stub contract checks passed")
    _smoke_simplified_cli_unified_wrapper_contract()
    print("✅ Simplified CLI unified-wrapper contract checks passed")
    _smoke_run_analysis_target_unified_delegation()
    print("✅ Canonical stage-target unified delegation checks passed")
    _smoke_config_override_overlay_contract()
    print("✅ Config override overlay extraction checks passed")
    _smoke_run_analysis_target_delegation_config_forwarding_contract()
    print("✅ Stage-target delegation config forwarding checks passed")
    _smoke_simplified_cli_wrapper_config_forwarding_contract()
    print("✅ Simplified CLI wrapper config forwarding checks passed")
    _smoke_user_facing_analysis_unified_path_contract()
    print("✅ User-facing post-docking commands unified-path checks passed")
    _smoke_clean_interaction_cli_dispatch_contract()
    print("✅ Clean interaction CLI dispatch contract checks passed")
    _smoke_clean_interaction_cli_shortcut_contract()
    print("✅ Clean interaction CLI shortcut checks passed")
    _smoke_non_gnina_rmsd_enforcement_contract()
    print("✅ Non-GNINA RMSD enforcement contract checks passed")
    _smoke_legacy_stage_routes_blocked_contract()
    print("✅ Legacy structure-quality/PyMOL route block checks passed")
    _smoke_non_gnina_bridge_no_legacy_stage_execution_contract()
    print("✅ Non-GNINA bridge legacy-stage detachment checks passed")
    _smoke_interactive_clean_alias_collapse_contract()
    print("✅ Interactive clean-alias collapse contract checks passed")
    _smoke_interactive_run_tracking_summary_surface_contract()
    print("✅ Interactive run-tracking summary surface checks passed")
    _smoke_interactive_favorite_flow_clean_followup_contract()
    print("✅ Interactive favorite-flow clean follow-up contract checks passed")
    _smoke_noninteractive_interaction_alias_to_clean_contract()
    print("✅ Non-interactive interaction alias-to-clean contract checks passed")
    _smoke_noninteractive_interaction_alias_cli_contract()
    print("✅ Non-interactive interaction CLI alias contract checks passed")
    _smoke_private_analysis_helpers_reject_interaction_alias_bypass_contract()
    print("✅ Private analysis helper alias-bypass rejection checks passed")
    _smoke_unified_output_topology_contract()
    print("✅ Unified output topology checks passed")
    _smoke_engine_detection_and_auto_routing()
    print("✅ Engine detection and auto-routing checks passed")
    _smoke_engine_detection_autodock4_parity_contract()
    print("✅ Engine detection AutoDock4 parity checks passed")
    _smoke_engine_hpc_adapter_parity_contract()
    print("✅ Engine HPC adapter parity checks passed")
    _smoke_run_tracking_schema_consistency_contract()
    print("✅ Run tracking schema consistency checks passed")
    _smoke_optional_feature_taxonomy_contract()
    print("✅ Optional feature taxonomy checks passed")
    _smoke_simplified_pipeline_interaction_contract_stages()
    print("✅ Simplified interaction-contract stage wiring checks passed")
    _smoke_pandamap_publication_artifact_contract()
    print("✅ PandaMap publication artifact contract checks passed")
    _smoke_pandamap_quality_gate_contract()
    print("✅ PandaMap quality gate checks passed")
    _smoke_engine_detection_cache_override_contract()
    print("✅ Engine detection cache-override checks passed")
    _smoke_engine_scope_filter_foundation()
    print("✅ Engine scope filter foundation checks passed")
    _smoke_engine_scope_picker_zero_guard()
    print("✅ Engine scope picker zero-selection guard checks passed")
    _smoke_engine_scope_filter_acceptance()
    print("✅ Engine scope filter acceptance checks passed")
    _smoke_no_valid_engine_abort_and_solo_semantics()
    print("✅ No-valid-engine abort and solo semantic checks passed")
    _smoke_cli_interactive_fallback_and_gnina_complex_export()
    print("✅ CLI fallback and GNINA SDF complex export checks passed")
    _smoke_sqlite_parity_and_skip_rationale()
    print("✅ SQLite dual-write parity and skip-rationale checks passed")

    if not args.skip_prep_matrix:
        _smoke_preparation_mode_matrix()
        print("✅ Preparation profile matrix smoke checks passed")

    if not args.skip_all_engines:
        artifacts = run_smoke_test(keep_temp=False)
        print("✅ All-engines dry-run smoke passed")
        print(f"   pairlist: {artifacts['pairlist']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
