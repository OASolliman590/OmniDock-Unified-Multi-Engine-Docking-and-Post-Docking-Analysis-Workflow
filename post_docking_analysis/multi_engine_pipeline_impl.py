"""
Engine-aware post-docking analysis for the canonical multi-engine project layout.
"""
from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import re
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from docking.project_layout import (
    ensure_engine_layout,
    ensure_numbered_output_layout,
    ensure_post_docking_compat_shim,
    load_manifest,
    load_pairlist,
    pair_intent_path,
    pairlist_path,
    post_docking_root,
    save_manifest,
    shared_receptors_dir,
)
from post_docking_analysis.complex_query import filter_frame_by_complex_query, tags_from_frame
from post_docking_analysis.consensus import (
    build_consensus_explainability,
    build_consensus_rankings,
    classify_hits_target_aware,
    normalize_hit_class_policy,
    normalize_consensus_mode,
    normalize_normalization_method,
    normalize_rescoring_scope,
    normalize_top_n,
    normalize_engine_scores,
    select_rescoring_candidates,
)
from post_docking_analysis.biology_integration import (
    attach_biology_annotations,
    compute_biology_correlations,
    load_biology_table,
)
from post_docking_analysis.correlation_analyzer import compute_cross_engine_rank_correlations
from post_docking_analysis.docking_parser import parse_autodock4_dlg, parse_vina_pdbqt
from post_docking_analysis.generate_scores_csv import generate_all_scores_csv, parse_gnina_log
from docking.runners.job_contract import parse_sdf_scores
from post_docking_analysis.ligand_naming import (
    build_ligand_name_mapping,
    ligand_mapping_to_dataframe,
    resolve_ligand_display_name,
)
from post_docking_analysis.protein_naming import (
    build_protein_name_mapping,
    extract_pdb_code,
    format_protein_label,
    mapping_to_dataframe,
)
from post_docking_analysis.report_generator import generate_hit_classification_reports
from post_docking_analysis.report_generator import (
    generate_analysis_root_index,
    generate_consolidated_run_summary,
    generate_dashboard_index,
    generate_start_here_index,
)
from post_docking_analysis.redocking_validation import run_redocking_validation
from post_docking_analysis.storage_sqlite import validate_csv_sqlite_parity, write_comparative_bundle
from post_docking_analysis.complex_validation import validate_complex_pdb_structure
from post_docking_analysis.pose_extractor import extract_best_poses_from_gnina
from post_docking_analysis.top_pose_selector import (
    build_top_pose_atlas,
    normalize_global_aggregation,
    normalize_selection_policy,
    write_top_pose_atlas,
)
from post_docking_analysis.artifact_graph import ArtifactGraph, ArtifactNode
from post_docking_analysis.visualization_suite import generate_visualization_suite
from post_docking_analysis.score_import import read_explicit_score_import
from post_docking_analysis.score_semantics import score_spec, sort_value, explicit_true
from post_docking_analysis.pose_geometry import content_hash, load_pose_molecule, selected_record


NORMALIZED_COLUMNS = [
    "engine",
    "tag",
    "protein",
    "ligand",
    "site_id",
    "pose",
    "affinity_kcal_mol",
    "score_name_primary",
    "score_primary",
    "score_name_secondary",
    "score_secondary",
    "rmsd_lb",
    "rmsd_ub",
    "pose_file",
    "log_file",
]

UNIFIED_COMPAT_COLUMNS = [
    "engine",
    "tag",
    "complex_name",
    "protein",
    "protein_label",
    "ligand",
    "site_id",
    "pose",
    "mode",
    "affinity_kcal_mol",
    "vina_affinity",
    "score_name_primary",
    "score_primary",
    "score_name_secondary",
    "score_secondary",
    "cnn_affinity",
    "cnn_score",
    "rmsd_lb",
    "rmsd_ub",
    "pose_file",
    "log_file",
]

_DEFAULT_RMSD_SCOPES = ("per_complex",)
_SUPPORTED_ANALYSIS_SCOPES = {
    "full",
    "comparison_only",
    "rescoring_only",
    "qc_only",
    "report_only",
    "top_pose_only",
}
_DAG_SCOPE_MAP = {
    "full": "reports",
    "comparison_only": "classified_hits",
    "rescoring_only": "reports",
    "qc_only": "reports",
    "structures_only": "complexes",
    "interactions": "interactions",
    "report_only": "reports",
    "top_pose_only": "best_poses",
}
_PAIR_METADATA_COLUMNS = (
    "pair_source",
    "selection_mode",
    "is_cocrystal_benchmark",
    "cocrystal_ligand_name",
    "reference_pose_file",
    "reference_ligand_file",
    "reference_source", "reference_pdb_id", "reference_frame_id", "receptor_frame_id",
    "reference_pose_index", "scoring_function",
)
NORMALIZED_COLUMNS += list(_PAIR_METADATA_COLUMNS) + ["score_provenance"]
_BEST_POSE_SELECTION_METRIC_ALIASES = {
    "auto": "auto",
    "affinity": "vina_affinity",
    "vina": "vina_affinity",
    "vina_affinity": "vina_affinity",
    "cnn": "cnn_affinity",
    "cnn_affinity": "cnn_affinity",
}
logger = logging.getLogger(__name__)
_RUN_TRACKING_STATES = {
    "not_started",
    "in_progress",
    "completed",
    "validated",
    "failed",
    "skipped",
    "needs_review",
}
_RUN_TRACKING_STEP_COLUMNS = [
    "index",
    "step",
    "required",
    "status",
    "started_at",
    "ended_at",
    "details",
    "error",
]
_RUN_TRACKING_OUTPUT_COLUMNS = [
    "relative_path",
    "category",
    "extension",
    "size_bytes",
    "modified_utc",
]
_RUN_TRACKING_STAGE_CONTRACT_KEYS = (
    "requested_run_rmsd",
    "requested_run_visualizations",
    "effective_run_rmsd",
    "effective_run_visualizations",
    "requested_run_prolif",
    "requested_run_ligplot",
    "effective_run_prolif",
    "effective_run_ligplot",
)
_RUN_TRACKING_OPTIONAL_STATUS_VALUES = {
    "completed",
    "skipped_disabled",
    "skipped_missing_dependency",
    "failed_error",
}
_OPTIONAL_DAG_FEATURES = {
    "visualizations": "visualizations",
    "prolif": "prolif",
    "pandamap": "pandamap",
    "poseview": "poseview",
    "pymol": "pymol",
}


def _normalize_rmsd_scopes(scopes: Optional[object]) -> Tuple[str, ...]:
    if scopes is None:
        return _DEFAULT_RMSD_SCOPES
    if isinstance(scopes, str):
        raw_tokens = re.split(r"[,\s]+", scopes.strip())
    else:
        raw_tokens = []
        for value in scopes:
            raw_tokens.extend(re.split(r"[,\s]+", str(value).strip()))
    normalized: List[str] = []
    for token in raw_tokens:
        value = str(token).strip().lower()
        if not value:
            continue
        if value in {"all", "*"}:
            return _DEFAULT_RMSD_SCOPES
        if value in _DEFAULT_RMSD_SCOPES and value not in normalized:
            normalized.append(value)
    return tuple(normalized) if normalized else _DEFAULT_RMSD_SCOPES


def normalize_best_pose_selection_metric(value: Optional[str]) -> str:
    token = str(value or "").strip().lower()
    return _BEST_POSE_SELECTION_METRIC_ALIASES.get(token, "auto")


def _optional_feature_record(
    status: str,
    reason: str,
    *,
    node_status: str = "",
    details: str = "",
    error: str = "",
) -> Dict[str, str]:
    normalized_status = status if status in _RUN_TRACKING_OPTIONAL_STATUS_VALUES else "failed_error"
    return {
        "status": normalized_status,
        "reason": str(reason or ""),
        "node_status": str(node_status or ""),
        "details": str(details or ""),
        "error": str(error or ""),
    }


class MultiEngineAnalysisPipeline:
    def __init__(
        self,
        project_dir: str,
        output_dir: Optional[str] = None,
        analysis_mode: str = "favorite_engine_continue",
        engine: Optional[str] = None,
        favorite_engine: Optional[str] = None,
        promote_exhaustive: bool = False,
        rerun_engine: Optional[str] = None,
        top_per_protein: int = 1,
        winner_only: bool = False,
        min_affinity_advantage: float = 0.0,
        max_rerun_pairs: int = 0,
        pair_allowlist: Optional[str] = None,
        consensus_mode: str = "dockbox_geometric",
        rescoring_scope: str = "top_n_per_protein",
        rescoring_top_n: int = 3,
        complex_query: Optional[str] = None,
        prompt_protein_names: bool = False,
        prompt_ligand_names: bool = False,
        enable_poseview: bool = False,
        rmsd_scopes: Optional[object] = None,
        resume_rmsd: bool = True,
        force_global_rmsd: bool = False,
        global_rmsd_defer_threshold: int = 150,
        shared_mapping_dir: Optional[str] = None,
        exclude_problematic_ligands: bool = False,
        positive_affinity_threshold: float = 0.0,
        minimum_pose_count: int = 1,
        rmsd_workers: int = 0,
        speed_profile: str = "standard",
        analysis_scope: str = "full",
        normalization_method: str = "per_engine_rank",
        biology_file: Optional[str] = None,
        biology_mapping_mode: str = "auto",
        hit_class_policy: str = "target_percentile",
        hit_class_strong_percentile: float = 0.10,
        hit_class_moderate_percentile: float = 0.35,
        top_pose_selection_policy: str = "best_affinity",
        top_pose_global_aggregation: str = "best_target",
        best_pose_selection_metric: str = "auto",
        engines_in_scope: Optional[List[str]] = None,
        excluded_engines: Optional[List[Dict[str, object]]] = None,
        scope_source: str = "detected_all",
        engine_preset_name: str = "",
        detected_engine_count: int = 0,
    ):
        self.project_dir = Path(project_dir).expanduser().resolve()
        self.compat_shim_info = ensure_post_docking_compat_shim(self.project_dir)
        if output_dir:
            self.output_dir = Path(output_dir).expanduser().resolve()
        else:
            base_output_root = post_docking_root(self.project_dir) / "sessions"
            self.output_dir = base_output_root / datetime.now(timezone.utc).strftime("interactive_%Y%m%d_%H%M%S")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.analysis_mode = analysis_mode
        self.engine = engine
        self.favorite_engine = favorite_engine
        self.promote_exhaustive = promote_exhaustive
        self.rerun_engine = rerun_engine
        self.top_per_protein = max(int(top_per_protein or 1), 1)
        self.winner_only = bool(winner_only)
        self.min_affinity_advantage = float(min_affinity_advantage or 0.0)
        self.max_rerun_pairs = max(int(max_rerun_pairs or 0), 0)
        self.pair_allowlist = pair_allowlist
        self.consensus_mode = normalize_consensus_mode(consensus_mode)
        self.rescoring_scope = normalize_rescoring_scope(rescoring_scope)
        self.rescoring_top_n = normalize_top_n(rescoring_top_n, default=3)
        self.complex_query = str(complex_query or "").strip()
        self.prompt_protein_names = bool(prompt_protein_names)
        self.prompt_ligand_names = bool(prompt_ligand_names)
        self.enable_poseview = bool(enable_poseview)
        self.shared_mapping_dir = (
            Path(shared_mapping_dir).expanduser().resolve() if shared_mapping_dir else None
        )
        if self.shared_mapping_dir is not None:
            self.shared_mapping_dir.mkdir(parents=True, exist_ok=True)
        self.mapping_root = self.shared_mapping_dir or self.output_dir
        self.rmsd_scopes = _normalize_rmsd_scopes(rmsd_scopes)
        self.resume_rmsd = bool(resume_rmsd)
        self.force_global_rmsd = bool(force_global_rmsd)
        self.exclude_problematic_ligands = bool(exclude_problematic_ligands)
        self.positive_affinity_threshold = float(positive_affinity_threshold)
        self.minimum_pose_count = max(int(minimum_pose_count or 1), 1)
        self.rmsd_workers = int(rmsd_workers or 0)
        speed_value = str(speed_profile or "standard").strip().lower()
        if speed_value not in {"standard", "fast"}:
            speed_value = "standard"
        self.speed_profile = speed_value
        self.fast_mode = speed_value == "fast"
        scope_value = str(analysis_scope or "full").strip().lower()
        if scope_value not in _SUPPORTED_ANALYSIS_SCOPES:
            scope_value = "full"
        self.analysis_scope = scope_value
        self.normalization_method = normalize_normalization_method(normalization_method)
        self.biology_file = str(biology_file or "").strip()
        self.biology_mapping_mode = str(biology_mapping_mode or "auto").strip().lower()
        self.hit_class_policy = normalize_hit_class_policy(hit_class_policy)
        self.hit_class_strong_percentile = float(hit_class_strong_percentile)
        self.hit_class_moderate_percentile = float(hit_class_moderate_percentile)
        self.top_pose_selection_policy = normalize_selection_policy(top_pose_selection_policy)
        self.top_pose_global_aggregation = normalize_global_aggregation(top_pose_global_aggregation)
        self.best_pose_selection_metric = normalize_best_pose_selection_metric(best_pose_selection_metric)
        self.engines_in_scope = [
            str(engine).strip().lower()
            for engine in (engines_in_scope or [])
            if str(engine).strip()
        ]
        self.excluded_engines = list(excluded_engines or [])
        self.scope_source = str(scope_source or "detected_all")
        self.engine_preset_name = str(engine_preset_name or "")
        self.detected_engine_count = max(int(detected_engine_count or 0), 0)
        try:
            self.global_rmsd_defer_threshold = max(int(global_rmsd_defer_threshold or 0), 0)
        except (TypeError, ValueError):
            self.global_rmsd_defer_threshold = 150
        self.manifest = load_manifest(self.project_dir)
        self.pairlist_df = self._load_pairlist_with_intent()
        self.selected_query_tags: Set[str] = set()
        if self.complex_query and not self.pairlist_df.empty:
            original_pair_count = len(self.pairlist_df)
            self.pairlist_df = filter_frame_by_complex_query(self.pairlist_df, self.complex_query)
            if self.pairlist_df.empty:
                raise ValueError("Complex query removed all canonical pairlist rows")
            self.selected_query_tags = tags_from_frame(self.pairlist_df)
            (self.output_dir / "complex_query.txt").write_text(self.complex_query + "\n", encoding="utf-8")
            (self.output_dir / "complex_query_summary.txt").write_text(
                f"query={self.complex_query}\noriginal_pairlist_rows={original_pair_count}\nfiltered_pairlist_rows={len(self.pairlist_df)}\n",
                encoding="utf-8",
            )
        self.rerun_manifest_file = ""
        self.run_tracking_dir = self.output_dir / "run_tracking"
        self.run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        self.top_pose_outputs: Dict[str, str] = {}
        self.sqlite_dual_write_enabled = self._resolve_sqlite_dual_write_flag()
        self.sqlite_write_summary: Dict[str, object] = {}
        if str(self.compat_shim_info.get("status", "")) == "legacy_dir_exists":
            logger.warning(
                "Legacy post-docking directory preserved at %s. New analysis outputs will be written to %s.",
                self.compat_shim_info.get("legacy_root"),
                self.compat_shim_info.get("canonical_root"),
            )

    @staticmethod
    def _utc_now_iso() -> str:
        return datetime.now(timezone.utc).isoformat()

    def _annotate_scope_columns(self, frame: pd.DataFrame) -> pd.DataFrame:
        if frame is None or frame.empty:
            return frame
        scoped = frame.copy()
        scoped["engines_in_scope"] = ",".join(self.engines_in_scope or [])
        scoped["scope_source"] = self.scope_source
        scoped["scoped_engine_count"] = int(len(self.engines_in_scope))
        scoped["detected_engine_count"] = int(self.detected_engine_count or len(self.manifest.get("engines", [])))
        if self.excluded_engines:
            scoped["excluded_engines"] = json.dumps(self.excluded_engines)
        else:
            scoped["excluded_engines"] = "[]"
        return scoped

    def _scope_discrepancy_notice(self) -> str:
        sessions_root = post_docking_root(self.project_dir) / "sessions"
        current_scope = tuple(self.engines_in_scope or [])
        if not sessions_root.exists() or not current_scope:
            return ""
        for session_dir in sorted(path for path in sessions_root.iterdir() if path.is_dir()):
            try:
                if session_dir.resolve() == self.output_dir.resolve():
                    continue
            except Exception:
                pass
            manifest_file = session_dir / "run_tracking" / "run_manifest.json"
            if not manifest_file.exists():
                continue
            try:
                payload = json.loads(manifest_file.read_text(encoding="utf-8"))
            except Exception:
                continue
            prior_scope = tuple(str(item).strip().lower() for item in (payload.get("engines_in_scope") or []) if str(item).strip())
            if prior_scope and prior_scope != current_scope:
                return (
                    f"Previous session `{session_dir.name}` used a different engine scope: "
                    f"{', '.join(prior_scope)}. Current session uses: {', '.join(current_scope)}."
                )
        return ""

    def resolve_dag_scope_artifact(self, scope: Optional[str] = None) -> str:
        normalized = str(scope or self.analysis_scope or "full").strip().lower()
        if normalized not in _DAG_SCOPE_MAP:
            raise ValueError(f"Unsupported DAG scope: {normalized}")
        return _DAG_SCOPE_MAP[normalized]

    def _dag_artifact_paths(self) -> Dict[str, Path]:
        numbered = ensure_numbered_output_layout(self.project_dir)
        session_root = self.output_dir
        reports_root = session_root / "reports"
        interactions_root = session_root / "interactions"
        visualizations_root = session_root / "visualizations"
        paths = {
            "engine_scope_config": numbered["post_metadata"] / "engine_scope_config.json",
            "analysis_parameters": numbered["post_metadata"] / "analysis_parameters.json",
            "dag_cache": numbered["post_metadata"] / "dag_cache.json",
            "dag_execution_report": numbered["post_metadata"] / "dag_execution_report.json",
            "raw_scores": numbered["post_scores_raw"] / "raw_scores.csv",
            "normalized_scores": numbered["post_scores_unified"] / "normalized_scores.csv",
            "validation_gate": numbered["post_scores_consensus"] / "validation_gate.json",
            "reference_baselines": numbered["post_scores_consensus"] / "reference_baselines.csv",
            "consensus_ranked": numbered["post_scores_consensus"] / "consensus_ranked.csv",
            "engine_agreement": numbered["post_scores_consensus"] / "engine_agreement.csv",
            "classified_hits": numbered["post_scores_consensus"] / "classified_hits.csv",
            "best_poses": session_root / "top_pose_ligand_performance" / "top_pose_per_ligand_global.csv",
            "best_poses_per_protein": session_root / "top_pose_ligand_performance" / "top_pose_per_ligand_per_protein.csv",
            "top_pose_confidence": session_root / "top_pose_ligand_performance" / "top_pose_confidence_metrics.csv",
            "top_pose_manifest": session_root / "top_pose_ligand_performance" / "top_pose_selection_manifest.json",
            "ligand_performance_summary": session_root / "top_pose_ligand_performance" / "ligand_performance_summary.csv",
            "complexes": session_root / "complexes",
            "prolif": interactions_root / "prolif",
            "pandamap": interactions_root / "pandamap",
            "poseview": interactions_root / "poseview",
            "pymol": interactions_root / "pymol",
            "interactions": interactions_root / "_interactions_complete.txt",
            "polypharmacology": reports_root / "polypharmacology",
            "comparative": reports_root / "comparative",
            "biology_correlation": reports_root / "biology_correlation",
            "visualizations": visualizations_root,
            "figures_manifest": visualizations_root / "figures_manifest.json",
            "reports": reports_root / "summary.txt",
            "protein_name_mapping": session_root / "protein_name_mapping.csv",
        }
        return paths

    def _write_engine_scope_config_artifact(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "analysis_mode": self.analysis_mode,
            "analysis_scope": self.analysis_scope,
            "engines_in_scope": list(self.engines_in_scope),
            "excluded_engines": list(self.excluded_engines),
            "scope_source": self.scope_source,
            "engine_preset_name": self.engine_preset_name,
            "favorite_engine": str(self.favorite_engine or self.engine or ""),
            "detected_engine_count": int(self.detected_engine_count or len(self.manifest.get("engines", []))),
        }
        serialized = json.dumps(payload, indent=2, sort_keys=True)
        if path.exists():
            try:
                if path.read_text(encoding="utf-8") == serialized:
                    return path
            except Exception:
                pass
        path.write_text(serialized, encoding="utf-8")
        return path

    def _write_analysis_parameters_artifact(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "analysis_scope": self.analysis_scope,
            "analysis_mode": self.analysis_mode,
            "normalization_method": self.normalization_method,
            "consensus_mode": self.consensus_mode,
            "rescoring_scope": self.rescoring_scope,
            "rescoring_top_n": int(self.rescoring_top_n),
            "hit_class_policy": self.hit_class_policy,
            "hit_class_strong_percentile": float(self.hit_class_strong_percentile),
            "hit_class_moderate_percentile": float(self.hit_class_moderate_percentile),
            "top_pose_selection_policy": self.top_pose_selection_policy,
            "top_pose_global_aggregation": self.top_pose_global_aggregation,
            "best_pose_selection_metric": self.best_pose_selection_metric,
            "favorite_engine": str(self.favorite_engine or self.engine or self.rerun_engine or ""),
            "biology_file": str(self.biology_file or ""),
            "biology_mapping_mode": self.biology_mapping_mode,
            "scientific_schema_version": 2,
            "complex_query": self.complex_query,
            "exclude_problematic_ligands": self.exclude_problematic_ligands,
            "positive_affinity_threshold": self.positive_affinity_threshold,
            "minimum_pose_count": self.minimum_pose_count,
            "enable_poseview": self.enable_poseview,
            "rmsd_scopes": list(self.rmsd_scopes),
            "speed_profile": self.speed_profile,
            "pair_allowlist": str(self.pair_allowlist or ""),
            "winner_only": self.winner_only,
            "min_affinity_advantage": self.min_affinity_advantage,
        }
        serialized = json.dumps(payload, indent=2, sort_keys=True)
        if path.exists():
            try:
                if path.read_text(encoding="utf-8") == serialized:
                    return path
            except Exception:
                pass
        path.write_text(serialized, encoding="utf-8")
        return path

    @staticmethod
    def _write_dag_placeholder_outputs(outputs: List[Path], *, node_name: str, details: str) -> Dict[str, object]:
        for output in outputs:
            if output.suffix:
                output.parent.mkdir(parents=True, exist_ok=True)
                if output.suffix.lower() == ".json":
                    output.write_text(
                        json.dumps({"node": node_name, "details": details}, indent=2),
                        encoding="utf-8",
                    )
                else:
                    output.write_text(f"{node_name}\n{details}\n", encoding="utf-8")
            else:
                output.mkdir(parents=True, exist_ok=True)
                (output / "_node_manifest.json").write_text(
                    json.dumps({"node": node_name, "details": details}, indent=2),
                    encoding="utf-8",
                )
        return {"details": details}

    def _dag_compute_raw_scores_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        scores = self._load_or_build_scores()
        output = paths["raw_scores"]
        output.parent.mkdir(parents=True, exist_ok=True)
        scores.to_csv(output, index=False)
        return {"details": f"score_rows={len(scores)}"}

    def _dag_compute_normalized_scores_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        raw_file = paths["raw_scores"]
        scores = pd.read_csv(raw_file) if raw_file.exists() else self._load_or_build_scores()
        output = paths["normalized_scores"]
        output.parent.mkdir(parents=True, exist_ok=True)
        scores.to_csv(output, index=False)
        return {"details": "current unified loader already returns normalized score rows"}

    def _dag_compute_validation_gate_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        normalized_file = paths["normalized_scores"]
        scores = pd.read_csv(normalized_file) if normalized_file.exists() else self._load_or_build_scores()
        if scores.empty:
            payload = {"validation_gate_state": "no_scores", "allow_reference_anchor": False}
            paths["validation_gate"].parent.mkdir(parents=True, exist_ok=True)
            paths["validation_gate"].write_text(json.dumps(payload, indent=2), encoding="utf-8")
            pd.DataFrame(columns=["protein", "engine", "scoring_function", "reference_affinity"]).to_csv(paths["reference_baselines"], index=False)
            return {"details": "no scores available for validation gate"}
        best_by_engine, _ = self._prepare_dag_consensus_inputs(scores, self._resolve_dag_consensus_strategy(paths))
        validation_bundle = run_redocking_validation(
            project_dir=self.project_dir,
            best_by_engine=best_by_engine,
            expected_reference_rows=self.pairlist_df,
            expected_engines=self.engines_in_scope,
            output_dir=paths["validation_gate"].parent,
        )
        payload = dict(validation_bundle.get("validation_gate") or validation_bundle.get("summary") or {})
        paths["validation_gate"].write_text(json.dumps(payload, indent=2), encoding="utf-8")
        return {"details": f"validation_gate_state={payload.get('validation_gate_state', '')}"}

    def _dag_compute_consensus_ranked_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        normalized_file = paths["normalized_scores"]
        normalized_scores = pd.read_csv(normalized_file) if normalized_file.exists() else self._load_or_build_scores()
        if normalized_scores.empty:
            pd.DataFrame().to_csv(paths["consensus_ranked"], index=False)
            pd.DataFrame().to_csv(paths["engine_agreement"], index=False)
            return {"warnings": ["no_valid_scores"], "details": "no normalized scores available"}
        strategy = self._resolve_dag_consensus_strategy(paths)
        best_by_engine, strategy_details = self._prepare_dag_consensus_inputs(normalized_scores, strategy)
        consensus_df = self._run_dag_consensus_strategy(best_by_engine, strategy=strategy)
        consensus_df = self._finalize_dag_consensus_output(consensus_df, strategy=strategy)
        consensus_df.to_csv(paths["consensus_ranked"], index=False)
        if strategy == "multi":
            engine_agreement = (
                best_by_engine.groupby("tag", dropna=False)
                .agg(
                    engine_support_count=("engine", "nunique"),
                )
                .reset_index()
            )
            engine_agreement = engine_agreement.merge(
                consensus_df[["tag", "winner_engine", "best_affinity_kcal_mol"]],
                on="tag", how="left", validate="one_to_one",
            )
        else:
            engine_agreement = pd.DataFrame(
                columns=["tag", "engine_support_count", "best_affinity_kcal_mol", "engine_mode", "engines_in_scope"]
            )
        engine_agreement.to_csv(paths["engine_agreement"], index=False)
        return {
            "details": f"strategy={strategy}; ranked_rows={len(consensus_df)}; primary={strategy_details.get('primary_score_name', '')}",
            "warnings": strategy_details.get("warnings", []),
        }

    def _load_engine_scope_config(self, paths: Dict[str, Path]) -> Dict[str, object]:
        config_file = paths.get("engine_scope_config")
        if config_file and Path(config_file).exists():
            try:
                payload = json.loads(Path(config_file).read_text(encoding="utf-8"))
                if isinstance(payload, dict):
                    return payload
            except Exception:
                pass
        return {
            "analysis_mode": self.analysis_mode,
            "engines_in_scope": list(self.engines_in_scope),
            "detected_engine_count": int(self.detected_engine_count or len(self.manifest.get("engines", []))),
        }

    def _resolve_dag_consensus_strategy(self, paths: Dict[str, Path]) -> str:
        scope_config = self._load_engine_scope_config(paths)
        engines = [str(item).strip().lower() for item in (scope_config.get("engines_in_scope") or []) if str(item).strip()]
        if len(engines) >= 2:
            return "multi"
        engine_name = engines[0] if engines else str(self._dag_preferred_engine() or "gnina").strip().lower()
        if engine_name == "gnina":
            return "single/gnina"
        if engine_name == "smina":
            return "single/smina"
        if engine_name == "autodock4":
            return "single/autodock4"
        return "single/vina"

    def _prepare_dag_consensus_inputs(
        self,
        normalized_scores: pd.DataFrame,
        strategy: str,
    ) -> Tuple[pd.DataFrame, Dict[str, object]]:
        frame = normalized_scores.copy()
        frame["engine"] = frame["engine"].astype(str).str.strip().str.lower()
        warnings: List[str] = []
        ranking_column = "affinity_kcal_mol"
        primary_name = "vina_affinity"
        secondary_name = ""
        secondary_column = ""

        if strategy != "multi":
            target_engine = strategy.split("/", 1)[-1].strip().lower()
            frame = frame[frame["engine"] == target_engine].copy()
            if frame.empty:
                return frame, {"warnings": [f"no_scores_for_{target_engine}"], "primary_score_name": primary_name}
            if target_engine == "gnina":
                cnn_affinity = pd.to_numeric(frame.get("cnn_affinity"), errors="coerce")
                if cnn_affinity.notna().any():
                    ranking_column = "cnn_affinity"
                    primary_name = "cnn_affinity"
                    secondary_name = "vina_affinity"
                    secondary_column = "affinity_kcal_mol"
                else:
                    warnings.append("gnina_cnn_affinity_unavailable_fallback_to_vina_affinity")
            elif target_engine == "smina":
                primary_name = "vina_affinity"
                secondary_name = "smina_scoring_function"
                secondary_column = "smina_scoring_function"
                smina_meta = (
                    (self.engine_detection_report.get("engines") or {}).get("smina", {})
                    if isinstance(self.engine_detection_report, dict)
                    else {}
                )
                frame["smina_scoring_function"] = json.dumps(smina_meta.get("smina_scoring_weights", {}) or {})
                if bool(smina_meta.get("inconsistent_scoring_weights", False)):
                    warnings.append("smina_inconsistent_scoring_weights")
            elif target_engine == "autodock4":
                primary_name = "autodock4_affinity"
            else:
                primary_name = "vina_affinity"
                secondary_name = "rmsd_lb"
                secondary_column = "rmsd_lb"

        best_by_engine = self._best_rows_by_group(frame, ["engine", "tag"], ranking_column).sort_values(
            ["engine", ranking_column, "tag"],
            ascending=[True, True, True],
        )
        best_by_engine["score_name_primary"] = primary_name
        best_by_engine["score_primary"] = pd.to_numeric(best_by_engine.get(ranking_column), errors="coerce")
        if secondary_name:
            best_by_engine["score_name_secondary"] = secondary_name
            if secondary_column in best_by_engine.columns:
                best_by_engine["score_secondary"] = best_by_engine.get(secondary_column)
            else:
                best_by_engine["score_secondary"] = np.nan
        else:
            best_by_engine["score_name_secondary"] = ""
            best_by_engine["score_secondary"] = np.nan
        return best_by_engine, {
            "warnings": warnings,
            "ranking_column": ranking_column,
            "primary_score_name": primary_name,
            "secondary_score_name": secondary_name,
        }

    def _run_dag_consensus_strategy(self, best_by_engine: pd.DataFrame, *, strategy: str) -> pd.DataFrame:
        if best_by_engine is None or best_by_engine.empty:
            return pd.DataFrame()
        ranking_column = "affinity_kcal_mol"
        if strategy == "single/gnina" and pd.to_numeric(best_by_engine.get("cnn_affinity"), errors="coerce").notna().any():
            ranking_column = "cnn_affinity"
        consensus_mode = self.consensus_mode if strategy == "multi" else "weighted_hybrid"
        consensus_df = build_consensus_rankings(
            best_by_engine=best_by_engine,
            consensus_mode=consensus_mode,
            favorite_engine=self.favorite_engine or self.rerun_engine or self._dag_preferred_engine(),
            normalization_method=self.normalization_method,
            score_column=ranking_column,
            expected_engines=self.engines_in_scope,
        )
        return consensus_df

    def _finalize_dag_consensus_output(self, frame: pd.DataFrame, *, strategy: str) -> pd.DataFrame:
        if frame is None:
            return pd.DataFrame()
        finalized = frame.copy()
        single_engine = strategy != "multi"
        finalized["engine_mode"] = "single" if single_engine else "multi"
        finalized["single_engine_mode"] = bool(single_engine)
        finalized["total_detected_engines"] = int(self.detected_engine_count or len(self.manifest.get("engines", [])))
        finalized["rank"] = pd.to_numeric(finalized.get("consensus_rank_global"), errors="coerce")
        finalized["normalized_score"] = pd.to_numeric(finalized.get("mean_rank_pct"), errors="coerce")
        finalized["engine_support_count"] = pd.to_numeric(finalized.get("agreement_count"), errors="coerce")
        finalized["primary_score_name"] = finalized.get("score_name_primary", "")
        finalized = self._annotate_scope_columns(finalized)
        return finalized

    def _dag_compute_classified_hits_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        consensus_file = paths["consensus_ranked"]
        if not consensus_file.exists():
            pd.DataFrame().to_csv(paths["classified_hits"], index=False)
            return {"warnings": ["no_consensus_scores"], "details": "consensus artifact missing"}
        consensus_df = pd.read_csv(consensus_file)
        if consensus_df.empty:
            consensus_df.to_csv(paths["classified_hits"], index=False)
            return {"details": "consensus artifact empty"}
        reference_baselines = pd.DataFrame()
        baselines_path = paths.get("reference_baselines")
        if baselines_path and Path(baselines_path).exists():
            try:
                reference_baselines = pd.read_csv(Path(baselines_path))
            except Exception:
                reference_baselines = pd.DataFrame()
        gate = json.loads(paths["validation_gate"].read_text(encoding="utf-8")) if paths["validation_gate"].exists() else {}
        policy = self.hit_class_policy
        if policy == "reference_anchor" and not explicit_true(gate.get("allow_reference_anchor")):
            policy = "target_percentile"
        classified = classify_hits_target_aware(
            consensus_df,
            policy=policy,
            strong_percentile=self.hit_class_strong_percentile,
            moderate_percentile=self.hit_class_moderate_percentile,
            reference_baselines=reference_baselines,
        )
        classified.to_csv(paths["classified_hits"], index=False)
        return {"details": f"classified_rows={len(classified)}"}

    def _dag_compute_top_pose_atlas_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        classified_file = paths["classified_hits"]
        normalized_file = paths["normalized_scores"]
        classified = pd.read_csv(classified_file) if classified_file.exists() else pd.DataFrame()
        normalized_scores = pd.read_csv(normalized_file) if normalized_file.exists() else self._load_or_build_scores()
        best_by_engine, _ = self._prepare_dag_consensus_inputs(normalized_scores, self._resolve_dag_consensus_strategy(paths)) if not normalized_scores.empty else (pd.DataFrame(), {})
        top_pose_payload = build_top_pose_atlas(
            best_by_engine=best_by_engine,
            consensus_df=classified,
            selection_policy=self.top_pose_selection_policy,
            consensus_mode=self.consensus_mode,
            global_aggregation=self.top_pose_global_aggregation,
            run_id=self.run_id,
            generated_at=self._utc_now_iso(),
            context={
                "analysis_scope": self.analysis_scope,
                "analysis_mode": self.analysis_mode,
                "consensus_mode": self.consensus_mode,
                "normalization_method": self.normalization_method,
            },
        )
        outputs = write_top_pose_atlas(top_pose_payload, paths["best_poses"].parent)
        details = f"top_pose_global={outputs.get('top_pose_per_ligand_global_file', '')}"
        return {"details": details}

    def _dag_best_pose_manifest_df(self, paths: Dict[str, Path]) -> pd.DataFrame:
        normalized_file = paths["normalized_scores"]
        normalized_scores = pd.read_csv(normalized_file) if normalized_file.exists() else self._load_or_build_scores()
        if normalized_scores.empty:
            return pd.DataFrame()
        strategy = self._resolve_dag_consensus_strategy(paths)
        candidates, selection = self._prepare_dag_consensus_inputs(normalized_scores, strategy=strategy)
        selection_metric = selection.get("ranking_column", "affinity_kcal_mol")
        consensus = pd.read_csv(paths["consensus_ranked"]) if paths["consensus_ranked"].is_file() else pd.DataFrame()
        if not consensus.empty and "winner_engine" in consensus:
            winners = consensus.drop_duplicates("tag").set_index("tag")["winner_engine"]
            manifest_df = candidates[candidates["engine"].eq(candidates["tag"].map(winners))].copy()
        else:
            manifest_df = candidates.sort_values(["tag", "engine"]).drop_duplicates("tag").copy()
        if manifest_df.empty:
            return manifest_df
        manifest_df["bridge_output_name"] = manifest_df["tag"].astype(str)
        manifest_df["selection_criterion"] = selection_metric
        if selection_metric == "cnn_affinity" and "cnn_affinity" in manifest_df.columns:
            manifest_df["selected_score"] = pd.to_numeric(manifest_df.get("cnn_affinity"), errors="coerce")
        else:
            manifest_df["selected_score"] = pd.to_numeric(manifest_df.get("affinity_kcal_mol"), errors="coerce")
        return manifest_df

    def _dag_preferred_engine(self) -> str:
        for candidate in (
            str(self.favorite_engine or "").strip().lower(),
            str(self.engine or "").strip().lower(),
            str(self.rerun_engine or "").strip().lower(),
        ):
            if candidate:
                return candidate
        if self.engines_in_scope:
            return str(self.engines_in_scope[0]).strip().lower()
        return "gnina"

    def _build_dag_simplified_bridge_pipeline(self, paths: Dict[str, Path]):
        from post_docking_analysis.simplified_pipeline import SimplifiedPostDockingPipeline

        engine_name = self._dag_preferred_engine()
        engine_layout = ensure_engine_layout(self.project_dir, engine_name)
        simplified_pipeline = SimplifiedPostDockingPipeline(
            sdf_folder=str(engine_layout["poses"]),
            log_folder=str(engine_layout["logs"]),
            receptors_folder=str(shared_receptors_dir(self.project_dir)),
            output_dir=str(self.output_dir),
            pairlist_file=str(pairlist_path(self.project_dir)),
            run_rmsd=False,
            run_visualizations=False,
            prompt_protein_names=False,
            prompt_ligand_names=False,
            enable_poseview=self.enable_poseview,
            complex_query=self.complex_query,
            rmsd_scopes=self.rmsd_scopes,
            resume_rmsd=self.resume_rmsd,
            force_global_rmsd=self.force_global_rmsd,
            global_rmsd_defer_threshold=self.global_rmsd_defer_threshold,
            shared_mapping_dir=str(self.mapping_root),
            exclude_problematic_ligands=self.exclude_problematic_ligands,
            positive_affinity_threshold=self.positive_affinity_threshold,
            minimum_pose_count=self.minimum_pose_count,
            rmsd_workers=self.rmsd_workers,
        )
        simplified_pipeline.pairlist_df = self.pairlist_df.copy()
        receptor_dir = shared_receptors_dir(self.project_dir)
        receptor_files = sorted(path for path in receptor_dir.glob("*") if path.is_file()) if receptor_dir.exists() else []
        protein_override_file = self.mapping_root / "protein_name_overrides.csv"
        ligand_override_file = self.mapping_root / "ligand_name_overrides.csv"
        protein_name_map = build_protein_name_mapping(
            receptor_files,
            pairlist_df=self.pairlist_df,
            overrides_file=protein_override_file if protein_override_file.exists() else None,
        )
        ligand_identifiers: List[str] = []
        if not self.pairlist_df.empty:
            for column in ("ligand", "ligand_name", "cocrystal_ligand_name"):
                if column in self.pairlist_df.columns:
                    ligand_identifiers.extend(self.pairlist_df[column].dropna().astype(str).tolist())
        ligand_name_map = build_ligand_name_mapping(
            ligand_identifiers,
            pairlist_df=self.pairlist_df,
            overrides_file=ligand_override_file if ligand_override_file.exists() else None,
        )
        simplified_pipeline.protein_name_map = protein_name_map
        simplified_pipeline.ligand_name_map = ligand_name_map
        manifest_df = self._dag_best_pose_manifest_df(paths)
        simplified_pipeline.scores_df = self._build_simplified_bridge_scores(
            pd.read_csv(paths["normalized_scores"]) if paths["normalized_scores"].exists() else self._load_or_build_scores()
        )
        simplified_pipeline.complexes = self._build_simplified_bridge_complexes(manifest_df)
        downstream_results = self._build_downstream_results(
            pd.read_csv(paths["normalized_scores"]) if paths["normalized_scores"].exists() else self._load_or_build_scores()
        )
        simplified_pipeline.results["affinity_analysis"] = {
            "best_poses": downstream_results.get("best_poses", pd.DataFrame()).copy()
        }
        return simplified_pipeline

    def _dag_compute_complexes_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        manifest_df = self._dag_best_pose_manifest_df(paths)
        if manifest_df.empty:
            paths["complexes"].mkdir(parents=True, exist_ok=True)
            pd.DataFrame(
                [{"status": "no_best_poses", "note": "No normalized scores or best poses were available for complex extraction."}]
            ).to_csv(paths["complexes"] / "pose_extraction_manifest.csv", index=False)
            return {"warnings": ["no_best_poses"], "details": "no complexes extracted"}
        extracted_count = self._extract_best_pose_complexes(manifest_df, paths["complexes"])
        best_poses_pdb_dir = self.output_dir / "best_poses_pdb"
        mirrored_count = self._mirror_complexes_to_best_poses(paths["complexes"], best_poses_pdb_dir)
        index_file = self._write_complex_export_index(
            output_dir=paths["complexes"],
            complexes_dir=paths["complexes"],
            best_poses_dir=best_poses_pdb_dir,
            favorite_engine=str(self.favorite_engine or self.engine or ""),
            source_table=manifest_df,
        )
        return {
            "details": f"complexes_extracted={extracted_count}; mirrored_best_pose_pdb={mirrored_count}",
            "complex_export_index": str(index_file),
        }

    def _build_complex_metadata(self) -> Dict[str, Dict[str, object]]:
        """Build per-complex display metadata for interaction tools.

        Builds from the top_pose_atlas / classified_hits canonical artifacts when
        available; falls back to an empty dict so callers can use file-stem labels.
        """
        metadata: Dict[str, Dict[str, object]] = {}
        try:
            paths = self._dag_artifact_paths()
            classified_file = paths.get("classified_hits")
            if classified_file and Path(classified_file).exists():
                df = pd.read_csv(classified_file)
                for _, row in df.iterrows():
                    tag = str(row.get("tag") or "").strip()
                    if not tag:
                        continue
                    protein_raw = str(row.get("protein") or "").strip()
                    ligand_raw = str(row.get("ligand") or "").strip()
                    protein_label = protein_raw
                    try:
                        from post_docking_analysis.protein_naming import resolve_protein_display_name
                        name_map = getattr(self, "protein_name_map", None) or {}
                        protein_label = resolve_protein_display_name(protein_raw, name_map) or protein_raw
                    except Exception:
                        pass
                    metadata[tag] = {
                        "protein_label": protein_label,
                        "ligand_name": ligand_raw,
                        "affinity": float(row["affinity_kcal_mol"]) if "affinity_kcal_mol" in row.index else None,
                        "hit_class": str(row.get("hit_class") or ""),
                    }
        except Exception:
            pass
        return metadata

    def _dag_compute_prolif_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        try:
            import matplotlib
            matplotlib.use("Agg", force=True)
        except Exception:
            pass
        from post_docking_analysis.prolif_interaction_maps import PROLIF_AVAILABLE, ProLifInteractionMapper

        if not PROLIF_AVAILABLE:
            paths["prolif"].mkdir(parents=True, exist_ok=True)
            return {"warnings": ["prolif_unavailable"], "details": "ProLIF not installed"}
        mapper = ProLifInteractionMapper()
        visual_root = paths["prolif"] / "best_pose" / "visuals"
        visual_root.mkdir(parents=True, exist_ok=True)
        complex_metadata = self._build_complex_metadata()

        def _worker(pdb_file: Path) -> Dict[str, object]:
            meta = complex_metadata.get(pdb_file.stem, {})
            ligand_hint, _, _ = self._derive_ligand_identity({"complex_name": pdb_file.stem})
            output_png = visual_root / f"{pdb_file.stem}__best_pose_2d.png"
            ok = mapper.create_interaction_map(
                pdb_file,
                output_png,
                ligand_resname=ligand_hint,
                dpi=300,
                figsize=(12, 8),
                write_html=True,
                display_title=f"ProLIF Best-Pose Map: {meta.get('protein_label') or pdb_file.stem}",
            )
            return {"complex": pdb_file.stem, "success": bool(ok)}

        return self._run_per_complex_interaction_tasks(
            paths["complexes"],
            node_name="prolif",
            worker=_worker,
            output_root=paths["prolif"],
        )

    def _dag_compute_pandamap_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        try:
            import matplotlib
            matplotlib.use("Agg", force=True)
        except Exception:
            pass
        from post_docking_analysis.publication_pandamap import PublicationPandaMapAnalyzer

        pandamap_root = paths["pandamap"]
        maps_2d_dir = pandamap_root / "maps_2d"
        maps_3d_dir = pandamap_root / "maps_3d"
        maps_2d_dir.mkdir(parents=True, exist_ok=True)
        maps_3d_dir.mkdir(parents=True, exist_ok=True)
        analyzer = PublicationPandaMapAnalyzer(
            conda_env="pandamap",
            config={
                "protein_name_map": getattr(self, "protein_name_map", {}) or {},
                "complex_metadata": self._build_complex_metadata(),
                "overwrite": False,
            },
            publication_mode=True,
        )

        def _worker(pdb_file: Path) -> Dict[str, object]:
            ligand_name = analyzer._resolve_ligand_name(pdb_file, "UNK")
            maps_payload = analyzer.generate_publication_2d_map(
                pdb_file,
                ligand_name=ligand_name,
                output_dir=maps_2d_dir,
                map_name=pdb_file.stem,
            )
            maps = maps_payload.get("files", {}) if isinstance(maps_payload, dict) else {}
            vis = analyzer.generate_publication_3d_visualization(
                pdb_file,
                ligand_name=ligand_name,
                output_dir=maps_3d_dir,
                vis_name=pdb_file.stem,
                interactive=True,
            )
            return {
                "complex": pdb_file.stem,
                "success": bool(maps or vis),
                "generated_2d": int(bool(maps)),
                "generated_3d": int(bool(vis)),
            }

        return self._run_per_complex_interaction_tasks(
            paths["complexes"],
            node_name="pandamap",
            worker=_worker,
            output_root=paths["pandamap"],
        )

    def _dag_compute_poseview_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        from post_docking_analysis.poseview_integration import PoseViewAnalyzer

        if not self.enable_poseview:
            paths["poseview"].mkdir(parents=True, exist_ok=True)
            return {"warnings": ["poseview_disabled"], "details": "PoseView stage disabled"}
        analyzer = PoseViewAnalyzer(
            {
                "output_formats": ["png", "svg", "pdf"],
                "max_complexes": 0,
            }
        )
        output_dir = paths["poseview"] / "2d_diagrams"
        output_dir.mkdir(parents=True, exist_ok=True)

        def _worker(pdb_file: Path) -> Dict[str, object]:
            ligand_hint, _, _ = analyzer._infer_ligand_from_complex_name(pdb_file.stem) or (None, None, None)
            result = analyzer.analyze_complex(pdb_file, output_dir, ligand_name=ligand_hint)
            return {
                "complex": pdb_file.stem,
                "success": bool(result.success),
                "error": str(result.error_message or ""),
            }

        return self._run_per_complex_interaction_tasks(
            paths["complexes"],
            node_name="poseview",
            worker=_worker,
            output_root=paths["poseview"],
        )

    def _dag_compute_pymol_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        from post_docking_analysis.pymol_visualizer import PyMOLVisualizer

        if shutil.which("pymol") is None:
            paths["pymol"].mkdir(parents=True, exist_ok=True)
            return {"warnings": ["pymol_binary_absent"], "details": "PyMOL executable not found"}
        visualizer = PyMOLVisualizer(paths["pymol"], {"visualization": {"dpi": 300}})

        def _worker(pdb_file: Path) -> Dict[str, object]:
            ligand_hint, _, _ = self._derive_ligand_identity({"complex_name": pdb_file.stem})
            session = visualizer.create_interaction_analysis(
                pdb_file,
                ligand_resname=ligand_hint,
                scene_name=pdb_file.stem,
            )
            return {"complex": pdb_file.stem, "success": bool(session and Path(session).exists())}

        return self._run_per_complex_interaction_tasks(
            paths["complexes"],
            node_name="pymol",
            worker=_worker,
            output_root=paths["pymol"],
        )

    def _interaction_worker_count(self, complex_count: int) -> int:
        configured = int(self.rmsd_workers or 0)
        if configured > 1:
            return configured
        cpu_cap = max(4, (os.cpu_count() or 4))
        return max(1, min(cpu_cap, complex_count))

    def _run_per_complex_interaction_tasks(
        self,
        complexes_dir: Path,
        *,
        node_name: str,
        worker: Callable[[Path], Dict[str, object]],
        output_root: Path,
    ) -> Dict[str, object]:
        output_root.mkdir(parents=True, exist_ok=True)
        complex_files = sorted(path for path in complexes_dir.glob("*.pdb") if not path.name.endswith(".receptor.pdb")) if complexes_dir.exists() else []
        if not complex_files:
            return {"warnings": [f"{node_name}_no_complexes"], "details": "no complexes available"}

        workers = self._interaction_worker_count(len(complex_files))
        successes = 0
        failures = 0
        warnings: List[str] = []
        per_complex_rows: List[Dict[str, object]] = []
        _per_complex_timeout = int(getattr(self, "interaction_timeout_seconds", 0) or 0) or 180
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(worker, pdb_file): pdb_file for pdb_file in complex_files}
            for future in as_completed(futures):
                pdb_file = futures[future]
                try:
                    result = future.result(timeout=_per_complex_timeout)
                except TimeoutError:
                    failures += 1
                    warnings.append(f"{node_name}:{pdb_file.stem}:timeout_after_{_per_complex_timeout}s")
                    per_complex_rows.append({"complex": pdb_file.stem, "success": False, "error": f"timeout>{_per_complex_timeout}s"})
                    continue
                except Exception as exc:
                    failures += 1
                    warnings.append(f"{node_name}:{pdb_file.stem}:{exc}")
                    per_complex_rows.append({"complex": pdb_file.stem, "success": False, "error": str(exc)})
                    continue
                per_complex_rows.append(result)
                if bool(result.get("success")):
                    successes += 1
                else:
                    failures += 1
                    warnings.append(
                        f"{node_name}:{pdb_file.stem}:{str(result.get('error') or 'failed')}"
                    )

        pd.DataFrame(per_complex_rows).to_csv(output_root / f"{node_name}_per_complex_results.csv", index=False)
        return {
            "warnings": warnings,
            "details": f"processed={len(complex_files)}; successful={successes}; failed={failures}; workers={workers}",
        }

    def _dag_compute_comparative_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        try:
            import matplotlib
            matplotlib.use("Agg", force=True)
        except Exception:
            pass
        scores = self._load_or_build_scores()
        self._write_comparative_reports(scores, analysis_scope=self.analysis_scope, preserve_dag_atlas=True)
        destination = paths["comparative"]
        destination.mkdir(parents=True, exist_ok=True)
        reports_dir = self.output_dir / "reports"
        copied = 0
        for candidate in (
            "combined_engine_scores.csv",
            "best_pose_per_tag_by_engine.csv",
            "best_engine_per_complex.csv",
            "engine_affinity_matrix.csv",
            "engine_summary.csv",
            "engine_rank_correlation_per_protein.csv",
            "engine_rank_correlation_global.csv",
            "consensus_ranked_hits.csv",
            "consensus_ranked_hits_with_classes.csv",
            "rescoring_candidates.csv",
            "consensus_explainability.json",
            "VALIDATION_SUMMARY.md",
            "validation_gate_status.json",
        ):
            source = reports_dir / candidate
            if source.exists():
                shutil.copy2(source, destination / source.name)
                copied += 1
        return {"details": f"comparative_reports_copied={copied}"}

    def _dag_compute_polypharmacology_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        source = self.output_dir / "reports" / "polypharmacology"
        if not source.exists():
            scores = self._load_or_build_scores()
            self._write_comparative_reports(scores, analysis_scope=self.analysis_scope, preserve_dag_atlas=True)
        copied = self._copy_matching_files(source, paths["polypharmacology"])
        return {"details": f"polypharmacology_files_copied={copied}"}

    def _dag_compute_biology_correlation_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        destination = paths["biology_correlation"]
        destination.mkdir(parents=True, exist_ok=True)
        reports_dir = self.output_dir / "reports"
        if not (reports_dir / "biology_mapping_report.json").exists():
            scores = self._load_or_build_scores()
            self._write_comparative_reports(scores, analysis_scope=self.analysis_scope, preserve_dag_atlas=True)
        copied = 0
        for candidate in (
            "biology_correlation_global.csv",
            "biology_correlation_per_protein.csv",
            "biology_mapping_report.json",
        ):
            source = reports_dir / candidate
            if source.exists():
                shutil.copy2(source, destination / source.name)
                copied += 1
        return {"details": f"biology_correlation_files_copied={copied}"}

    def _dag_compute_visualizations_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        output_root = paths["visualizations"]
        output_root.mkdir(parents=True, exist_ok=True)
        manifest = generate_visualization_suite(
            classified_hits_file=paths["classified_hits"],
            consensus_ranked_file=paths["consensus_ranked"],
            engine_agreement_file=paths["engine_agreement"],
            normalized_scores_file=paths["normalized_scores"],
            engine_scope_config_file=paths["engine_scope_config"],
            validation_gate_file=paths["validation_gate"],
            top_pose_global_file=paths["best_poses"],
            output_dir=output_root,
            protein_name_mapping_file=paths.get("protein_name_mapping"),
        )
        manifest_file = output_root / "figures_manifest.json"
        if manifest_file.exists() and manifest_file.resolve() != paths["figures_manifest"].resolve():
            paths["figures_manifest"].parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(manifest_file, paths["figures_manifest"])
        summary = dict(manifest.get("summary") or {})
        generated = int(summary.get("generated", 0))
        skipped = int(summary.get("skipped", 0))
        group_summary = summary.get("groups") or {}
        details = f"figures_generated={generated} skipped={skipped} groups={json.dumps(group_summary, sort_keys=True)}"
        warnings: List[str] = []
        if generated <= 0:
            warnings.append("visualization_suite_generated_zero_figures")
        return {"details": details, "warnings": warnings}

    def _dag_compute_reports_node(self, paths: Dict[str, Path]) -> Dict[str, object]:
        reports_file = paths["reports"]
        existing_summary = self.output_dir / "reports" / "summary.txt"
        if existing_summary.exists():
            reports_file.parent.mkdir(parents=True, exist_ok=True)
            if existing_summary.resolve() != reports_file.resolve():
                shutil.copy2(existing_summary, reports_file)
        elif not reports_file.exists():
            reports_file.parent.mkdir(parents=True, exist_ok=True)
            reports_file.write_text(
                "Report generation completed through comparative writer, but summary artifact was absent.\n",
                encoding="utf-8",
            )

        figures_manifest_file = paths.get("figures_manifest")
        figures_embedded = 0
        if figures_manifest_file and Path(figures_manifest_file).exists():
            try:
                manifest_payload = json.loads(Path(figures_manifest_file).read_text(encoding="utf-8"))
                figure_rows = manifest_payload.get("figures") if isinstance(manifest_payload, dict) else []
                if not isinstance(figure_rows, list):
                    figure_rows = []
                generated_rows = []
                for row in figure_rows:
                    if not isinstance(row, dict):
                        continue
                    if not bool(row.get("generated", False)):
                        continue
                    path_token = str(row.get("path", "")).strip()
                    if not path_token:
                        continue
                    generated_rows.append(
                        {
                            "name": str(row.get("name", "")).strip(),
                            "title": str(row.get("title", "")).strip(),
                            "group": str(row.get("group", "")).strip(),
                            "path": path_token,
                        }
                    )
                figures_embedded = len(generated_rows)

                numbered_layout = ensure_numbered_output_layout(self.project_dir)
                start_here_file = numbered_layout["reports_root_numbered"] / "START_HERE.md"
                if start_here_file.exists():
                    start_here_lines = start_here_file.read_text(encoding="utf-8").splitlines()
                    block_start = "<!-- BEGIN: VISUALIZATION_SUITE -->"
                    block_end = "<!-- END: VISUALIZATION_SUITE -->"
                    def _rel_for_start_here(path_value: Path) -> str:
                        try:
                            return str(path_value.resolve().relative_to(Path(self.project_dir).expanduser().resolve()))
                        except Exception:
                            return str(path_value)
                    section_lines = [
                        block_start,
                        "",
                        "## Visualization Suite",
                        f"- figures_manifest: `{_rel_for_start_here(Path(figures_manifest_file))}`",
                        f"- generated_figures: `{figures_embedded}`",
                    ]
                    if generated_rows:
                        for row in generated_rows:
                            figure_path = Path(str(row.get("path", "")).strip()).expanduser()
                            rel_path = _rel_for_start_here(figure_path)
                            label = row.get("title") or row.get("name") or "figure"
                            group = row.get("group") or "misc"
                            section_lines.append(f"- `{label}` ({group}) -> `{rel_path}`")
                    else:
                        section_lines.append("- No generated figures were recorded for this run.")
                    section_lines.extend(["", block_end])

                    joined = "\n".join(start_here_lines)
                    if block_start in joined and block_end in joined:
                        prefix, remainder = joined.split(block_start, 1)
                        _, suffix = remainder.split(block_end, 1)
                        updated = prefix.rstrip() + "\n\n" + "\n".join(section_lines) + suffix
                    else:
                        updated = joined.rstrip() + "\n\n" + "\n".join(section_lines) + "\n"
                    start_here_file.write_text(updated, encoding="utf-8")

                # Publish symlinks to 7-Reports/figures/ so figures are always one click away.
                try:
                    figures_dir = numbered_layout["reports_figures"]
                    figures_dir.mkdir(parents=True, exist_ok=True)
                    for row in generated_rows:
                        fig_path = Path(str(row.get("path", "")).strip()).expanduser().resolve()
                        if fig_path.is_file():
                            link = figures_dir / fig_path.name
                            link.unlink(missing_ok=True)
                            link.symlink_to(fig_path)
                        elif fig_path.is_dir():
                            subdir = figures_dir / fig_path.name
                            subdir.mkdir(parents=True, exist_ok=True)
                            for child in sorted(fig_path.glob("*.png")):
                                link = subdir / child.name
                                link.unlink(missing_ok=True)
                                link.symlink_to(child.resolve())
                except Exception as exc:
                    logger.warning("Could not publish figure symlinks to reports_figures: %s", exc)
            except Exception as exc:
                logger.warning("Could not embed visualization suite section into START_HERE.md: %s", exc)

        # Second pass: symlink any .png files produced by the older visualization path
        # (cross_engine_same_pair_batches/, cross_engine_pair_heatmap.png, etc.) that are
        # NOT registered in figures_manifest.json and therefore not covered by the first pass.
        try:
            vis_root = paths.get("visualizations")
            if vis_root and Path(vis_root).is_dir():
                _nl2 = ensure_numbered_output_layout(self.project_dir)
                figures_dir = _nl2["reports_figures"]
                figures_dir.mkdir(parents=True, exist_ok=True)
                already_linked = {lnk.resolve() for lnk in figures_dir.rglob("*.png")}
                for png in sorted(Path(vis_root).rglob("*.png")):
                    if png.resolve() in already_linked:
                        continue
                    # Preserve one level of sub-directory grouping.
                    rel = png.relative_to(Path(vis_root))
                    dest = figures_dir / rel
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    dest.unlink(missing_ok=True)
                    dest.symlink_to(png.resolve())
        except Exception as exc:
            logger.warning("Could not publish non-manifest figure symlinks to reports_figures: %s", exc)

        return {"details": f"summary_file={reports_file}; figures_embedded={figures_embedded}"}

    def _dag_compute_placeholder_node(self, node_name: str, outputs: List[Path], details: str) -> Dict[str, object]:
        return self._write_dag_placeholder_outputs(outputs, node_name=node_name, details=details)

    def build_artifact_graph(self, *, max_workers: Optional[int] = None) -> ArtifactGraph:
        paths = self._dag_artifact_paths()
        self._write_engine_scope_config_artifact(paths["engine_scope_config"])
        self._write_analysis_parameters_artifact(paths["analysis_parameters"])
        graph = ArtifactGraph(
            cache_file=paths["dag_cache"],
            report_file=paths["dag_execution_report"],
            max_workers=max_workers or max(1, int(self.rmsd_workers or 1)),
        )

        def _register(
            name: str,
            inputs: List[Path],
            outputs: List[Path],
            compute,
            optional: bool = False,
            cacheable: bool = True,
            optional_inputs: Optional[List[Path]] = None,
        ) -> None:
            graph.register(
                ArtifactNode(
                    name=name,
                    inputs=[str(path) for path in inputs],
                    outputs=[str(path) for path in outputs],
                    compute=compute,
                    optional=optional,
                    cacheable=cacheable,
                    optional_inputs=[str(path) for path in (optional_inputs or [])],
                )
            )

        pairlist_file = pairlist_path(self.project_dir)
        source_inputs = []
        for source_engine in (self.engines_in_scope or self.manifest.get("engines", [])):
            layout = ensure_engine_layout(self.project_dir, source_engine)
            source_inputs.extend([layout["poses"], layout["logs"], layout["scores"] / "all_scores.csv", layout["scores"] / "normalized_scores.import.json"])
            if (layout["scores"] / "normalized_scores.import.json").is_file():
                source_inputs.append(layout["scores"] / "normalized_scores.csv")
        reference_inputs = []
        for column in ("reference_pose_file", "reference_ligand_file"):
            if column in self.pairlist_df:
                for value in self.pairlist_df[column].dropna().astype(str):
                    if value.strip():
                        path = Path(value)
                        reference_inputs.append(path if path.is_absolute() else self.project_dir / path)
        receptor_inputs = [shared_receptors_dir(self.project_dir)]
        biology_inputs = [Path(self.biology_file)] if self.biology_file else []
        extra_inputs = [Path(self.pair_allowlist)] if self.pair_allowlist else []
        _register(
            "raw_scores",
            [pairlist_file, paths["engine_scope_config"], paths["analysis_parameters"], *source_inputs, *receptor_inputs, *extra_inputs],
            [paths["raw_scores"]],
            lambda: self._dag_compute_raw_scores_node(paths),
            optional_inputs=source_inputs + receptor_inputs,
        )
        _register(
            "normalized_scores",
            [paths["raw_scores"], paths["analysis_parameters"]],
            [paths["normalized_scores"]],
            lambda: self._dag_compute_normalized_scores_node(paths),
        )
        _register(
            "validation_gate",
            [paths["normalized_scores"], *reference_inputs, *receptor_inputs, paths["analysis_parameters"]],
            [paths["validation_gate"], paths["reference_baselines"]],
            lambda: self._dag_compute_validation_gate_node(paths),
            optional_inputs=reference_inputs + receptor_inputs,
        )
        _register(
            "consensus_ranked",
            [paths["normalized_scores"], paths["engine_scope_config"], paths["analysis_parameters"], *source_inputs, *receptor_inputs],
            [paths["consensus_ranked"], paths["engine_agreement"]],
            lambda: self._dag_compute_consensus_ranked_node(paths),
            optional_inputs=source_inputs + receptor_inputs,
        )
        _register(
            "classified_hits",
            [paths["consensus_ranked"], paths["validation_gate"], paths["reference_baselines"], paths["analysis_parameters"]],
            [paths["classified_hits"]],
            lambda: self._dag_compute_classified_hits_node(paths),
        )
        _register(
            "top_pose_atlas",
            [paths["classified_hits"], paths["normalized_scores"], paths["analysis_parameters"]],
            [paths["best_poses"], paths["ligand_performance_summary"], paths["best_poses_per_protein"], paths["top_pose_confidence"], paths["top_pose_manifest"]],
            lambda: self._dag_compute_top_pose_atlas_node(paths),
        )
        _register(
            "visualizations",
            [
                paths["classified_hits"],
                paths["consensus_ranked"],
                paths["engine_agreement"],
                paths["normalized_scores"],
                paths["engine_scope_config"],
                paths["validation_gate"],
                paths["best_poses"],
            ],
            [paths["figures_manifest"]],
            lambda: self._dag_compute_visualizations_node(paths),
            optional=True,
        )
        _register(
            "complexes",
            [paths["classified_hits"], paths["best_poses"], *source_inputs, *receptor_inputs],
            [paths["complexes"]],
            lambda: self._dag_compute_complexes_node(paths),
            optional_inputs=source_inputs + receptor_inputs,
        )
        for node_name, path_key in (
            ("prolif", "prolif"),
            ("pandamap", "pandamap"),
            ("poseview", "poseview"),
            ("pymol", "pymol"),
        ):
            _register(
                node_name,
                [paths["complexes"], paths["classified_hits"]],
                [paths[path_key]],
                (
                    lambda node_name=node_name: self._dag_compute_prolif_node(paths)
                    if node_name == "prolif"
                    else self._dag_compute_pandamap_node(paths)
                    if node_name == "pandamap"
                    else self._dag_compute_poseview_node(paths)
                    if node_name == "poseview"
                    else self._dag_compute_pymol_node(paths)
                ),
                optional=True,
            )
        for node_name, path_key in (
            ("polypharmacology", "polypharmacology"),
            ("comparative", "comparative"),
            ("biology_correlation", "biology_correlation"),
        ):
            _register(
                node_name,
                [paths["classified_hits"], paths["best_poses"], paths["best_poses_per_protein"], paths["ligand_performance_summary"], paths["top_pose_confidence"], paths["top_pose_manifest"], *biology_inputs, paths["analysis_parameters"]] if node_name == "comparative" else [paths["comparative"], *biology_inputs],
                [paths[path_key]],
                (
                    lambda node_name=node_name: self._dag_compute_polypharmacology_node(paths)
                    if node_name == "polypharmacology"
                    else self._dag_compute_comparative_node(paths)
                    if node_name == "comparative"
                    else self._dag_compute_biology_correlation_node(paths)
                ),
            )
        _register(
            "interactions",
            [paths["prolif"], paths["pandamap"], paths["poseview"], paths["pymol"]],
            [paths["interactions"]],
            lambda: self._dag_compute_placeholder_node(
                "interactions",
                [paths["interactions"]],
                "Registered for DAG migration; interaction branch aggregate node for scope resolution.",
            ),
            optional_inputs=[paths["prolif"], paths["pandamap"], paths["poseview"], paths["pymol"]],
        )
        _register(
            "reports",
            [
                paths["classified_hits"],
                paths["engine_agreement"],
                paths["analysis_parameters"],
                paths["polypharmacology"],
                paths["comparative"],
                paths["biology_correlation"],
                # Interaction nodes (prolif/pandamap/poseview/pymol) are intentionally
                # NOT listed here. They are optional slow nodes; the reports node reads
                # whatever output files already exist rather than blocking on them.
                paths["figures_manifest"],
            ],
            [paths["reports"]],
            lambda: self._dag_compute_reports_node(paths),
            cacheable=False,
            optional_inputs=[paths["figures_manifest"]],
        )
        return graph

    def request_artifact_scope(
        self,
        scope: Optional[str] = None,
        *,
        force: bool = False,
        max_workers: Optional[int] = None,
    ) -> Dict[str, object]:
        paths = self._dag_artifact_paths()
        artifact_key = self.resolve_dag_scope_artifact(scope)
        artifact_output_map = {
            "reports": paths["reports"],
            "classified_hits": paths["classified_hits"],
            "complexes": paths["complexes"],
            "interactions": paths["interactions"],
            "best_poses": paths["best_poses"],
        }
        target_path = artifact_output_map[artifact_key]
        graph = self.build_artifact_graph(max_workers=max_workers)
        report = graph.request(str(target_path), force=force)
        report["scope_requested"] = str(scope or self.analysis_scope or "full")
        report["artifact_requested_key"] = artifact_key
        report["artifact_requested_path"] = str(target_path)
        graph._write_report(report)
        return report

    @staticmethod
    def _count_heavy_atoms_from_sdf(sdf_file: Path) -> Optional[int]:
        try:
            lines = sdf_file.read_text(encoding="utf-8", errors="replace").splitlines()
        except Exception:
            return None
        if len(lines) < 4:
            return None
        counts_line = lines[3]
        try:
            atom_count = int(counts_line[0:3].strip())
        except Exception:
            return None
        heavy_atoms = 0
        for atom_line in lines[4 : 4 + atom_count]:
            if len(atom_line) < 34:
                continue
            element = atom_line[31:34].strip().upper()
            if element and element != "H":
                heavy_atoms += 1
        return heavy_atoms or None

    @staticmethod
    def _count_heavy_atoms_from_pdbqt(pdbqt_file: Path) -> Optional[int]:
        try:
            heavy_atoms = 0
            for raw_line in selected_record(pdbqt_file, 1).splitlines():
                if not raw_line.startswith(("ATOM", "HETATM")):
                    continue
                atom_type = raw_line.split()[-1]
                if atom_type in {"H", "HD", "HS", "G0", "G1", "G2", "G3"}:
                    continue
                if not atom_type or not atom_type[0].isalpha():
                    return None
                heavy_atoms += 1
            return heavy_atoms or None
        except Exception:
            return None

    def _infer_heavy_atom_count(self, pose_file: Path) -> Optional[int]:
        suffix = pose_file.suffix.lower()
        if suffix == ".sdf":
            return self._count_heavy_atoms_from_sdf(pose_file)
        if suffix == ".pdbqt":
            return self._count_heavy_atoms_from_pdbqt(pose_file)
        return None

    def _load_pairlist_with_intent(self) -> pd.DataFrame:
        pairlist_df = load_pairlist(self.project_dir)
        intent_file = pair_intent_path(self.project_dir, self.manifest.get("layout_profile"))
        if intent_file.exists():
            try:
                intent_df = pd.read_csv(intent_file)
                required = {"receptor", "site_id", "ligand"}
                if required.issubset(intent_df.columns):
                    return intent_df
            except Exception as exc:
                logger.warning("Could not load pair_intent.csv (%s): %s", intent_file, exc)
        return pairlist_df

    def _resolve_sqlite_dual_write_flag(self) -> bool:
        """
        Determine whether CSV+SQLite dual-write is enabled.

        Priority:
        1. Explicit environment toggle (`DOCKFORGE_SQLITE_DUAL_WRITE`).
        2. Workflow state feature flag (`enable_sqlite_dual_write`).
        3. Default `False`.
        """
        env_token = str(os.environ.get("DOCKFORGE_SQLITE_DUAL_WRITE", "") or "").strip().lower()
        if env_token in {"1", "true", "yes", "on"}:
            return True
        if env_token in {"0", "false", "no", "off"}:
            return False
        try:
            from workflow.state import get_feature_flags

            flags = get_feature_flags(self.project_dir)
            return bool(flags.get("enable_sqlite_dual_write", False))
        except Exception:
            return False

    @staticmethod
    def _classify_output_path(relative_path: Path) -> str:
        parts = [str(part).lower() for part in relative_path.parts]
        if not parts:
            return "other"
        head = parts[0]
        if head == "top_pose_ligand_performance":
            return "top_pose_atlas"
        if head == "reports":
            return "analysis_report"
        if head == "raw_data":
            return "raw_data"
        if head == "visualizations":
            if len(parts) >= 2 and parts[1] == "2d":
                return "visualization_2d"
            if len(parts) >= 2 and parts[1] == "3d":
                return "visualization_3d"
            return "visualization"
        if head == "interactions":
            return "interaction_output"
        if head == "analysis":
            return "canonical_analysis"
        if head == "favorite_engine":
            return "favorite_engine_output"
        if head == "run_tracking":
            return "run_tracking"
        if head.endswith("_engine"):
            return "engine_output"
        return "other"

    def _canonical_top_pose_root(self) -> Path:
        """
        Resolve canonical top-pose atlas location.
        """
        return post_docking_root(self.project_dir) / "top_pose_ligand_performance"

    def _canonical_polypharmacology_root(self) -> Path:
        """
        Resolve canonical docking-biology output location.
        """
        return post_docking_root(self.project_dir) / "polypharmacology"

    @staticmethod
    def _prune_empty_directories(root: Path) -> int:
        removed = 0
        directories = sorted(
            (path for path in root.rglob("*") if path.is_dir()),
            key=lambda path: len(path.parts),
            reverse=True,
        )
        for directory in directories:
            if directory == root:
                continue
            try:
                if any(directory.iterdir()):
                    continue
                directory.rmdir()
                removed += 1
            except Exception:
                continue
        return removed

    def _infer_stage_contract(self, step_rows: List[Dict[str, object]]) -> Dict[str, object]:
        step_names = {str(row.get("step", "") or "").strip().lower() for row in step_rows}
        analysis_scope = str(getattr(self, "analysis_scope", "") or "").strip().lower()

        requested_run_rmsd = "rmsd_analysis" in step_names or analysis_scope in {
            "full",
            "comparison_only",
            "rescoring_only",
            "qc_only",
            "report_only",
            "top_pose_only",
        }
        requested_run_visualizations = "visualizations" in step_names or analysis_scope in {"full", "comparison_only"}
        requested_run_prolif = "prolif" in step_names
        requested_run_ligplot = "ligplot" in step_names

        stage_contract: Dict[str, object] = {
            "requested_run_rmsd": bool(requested_run_rmsd),
            "requested_run_visualizations": bool(requested_run_visualizations),
            "effective_run_rmsd": bool(requested_run_rmsd),
            "effective_run_visualizations": bool(requested_run_visualizations),
            "requested_run_prolif": bool(requested_run_prolif),
            "requested_run_ligplot": bool(requested_run_ligplot),
            "effective_run_prolif": bool(requested_run_prolif),
            "effective_run_ligplot": bool(requested_run_ligplot),
        }
        for key in _RUN_TRACKING_STAGE_CONTRACT_KEYS:
            stage_contract.setdefault(key, False)
        return stage_contract

    @staticmethod
    def _classify_optional_node_status(
        *,
        feature: str,
        node_status: str,
        details: str,
        blocked_by: List[str],
    ) -> Dict[str, str]:
        status = str(node_status or "").strip().lower()
        detail_text = str(details or "").strip()
        detail_lower = detail_text.lower()
        blocked_text = ",".join(str(item) for item in blocked_by if str(item).strip())
        error_text = blocked_text or detail_text

        dependency_markers = (
            "not installed",
            "executable not found",
            "binary absent",
            "not found",
            "unavailable",
        )

        if status in {"completed", "validated", "cache_hit"}:
            return _optional_feature_record(
                "completed",
                f"{feature}_completed",
                node_status=node_status,
                details=details,
                error=error_text,
            )

        if status == "completed_with_warnings":
            if "disabled" in detail_lower:
                return _optional_feature_record(
                    "skipped_disabled",
                    f"{feature}_disabled",
                    node_status=node_status,
                    details=details,
                    error=error_text,
                )
            if any(marker in detail_lower for marker in dependency_markers):
                return _optional_feature_record(
                    "skipped_missing_dependency",
                    f"{feature}_missing_dependency",
                    node_status=node_status,
                    details=details,
                    error=error_text,
                )
            if "no complexes available" in detail_lower or f"{feature}_no_complexes" in detail_lower:
                return _optional_feature_record(
                    "skipped_disabled",
                    f"{feature}_no_complexes",
                    node_status=node_status,
                    details=details,
                    error=error_text,
                )

            success_match = re.search(r"successful=(\d+)", detail_lower)
            failed_match = re.search(r"failed=(\d+)", detail_lower)
            generated_match = re.search(r"figures_generated=(\d+)", detail_lower)
            if generated_match and int(generated_match.group(1)) <= 0:
                return _optional_feature_record(
                    "failed_error",
                    f"{feature}_generated_zero_outputs",
                    node_status=node_status,
                    details=details,
                    error=error_text,
                )
            if success_match and int(success_match.group(1)) <= 0:
                if failed_match and int(failed_match.group(1)) > 0:
                    return _optional_feature_record(
                        "failed_error",
                        f"{feature}_all_tasks_failed",
                        node_status=node_status,
                        details=details,
                        error=error_text,
                    )
            return _optional_feature_record(
                "completed",
                f"{feature}_completed_with_warnings",
                node_status=node_status,
                details=details,
                error=error_text,
            )

        if status in {"skipped_optional", "skipped"}:
            if "disabled" in detail_lower:
                return _optional_feature_record(
                    "skipped_disabled",
                    f"{feature}_disabled",
                    node_status=node_status,
                    details=details,
                    error=error_text,
                )
            if any(marker in detail_lower for marker in dependency_markers):
                return _optional_feature_record(
                    "skipped_missing_dependency",
                    f"{feature}_missing_dependency",
                    node_status=node_status,
                    details=details,
                    error=error_text,
                )
            if "no complexes available" in detail_lower or f"{feature}_no_complexes" in detail_lower:
                return _optional_feature_record(
                    "skipped_disabled",
                    f"{feature}_no_complexes",
                    node_status=node_status,
                    details=details,
                    error=error_text,
                )
            return _optional_feature_record(
                "failed_error",
                f"{feature}_skipped_with_error",
                node_status=node_status,
                details=details,
                error=error_text,
            )

        if status in {"failed", "blocked_by_failure"}:
            return _optional_feature_record(
                "failed_error",
                f"{feature}_failed",
                node_status=node_status,
                details=details,
                error=error_text,
            )

        return _optional_feature_record(
            "skipped_disabled",
            f"{feature}_not_requested_by_scope",
            node_status=node_status,
            details=details,
            error=error_text,
        )

    def _infer_optional_feature_statuses(
        self,
        step_rows: List[Dict[str, object]],
    ) -> Dict[str, Dict[str, str]]:
        dag_report = getattr(self, "dag_execution_report", {})
        report = dag_report if isinstance(dag_report, dict) else {}
        node_rows = report.get("nodes") if isinstance(report.get("nodes"), dict) else {}
        if not isinstance(node_rows, dict):
            node_rows = {}

        optional_payload: Dict[str, Dict[str, str]] = {}
        for feature_name, node_name in _OPTIONAL_DAG_FEATURES.items():
            payload = node_rows.get(node_name)
            if not isinstance(payload, dict):
                optional_payload[feature_name] = _optional_feature_record(
                    "skipped_disabled",
                    f"{feature_name}_not_requested_by_scope",
                )
                continue
            optional_payload[feature_name] = self._classify_optional_node_status(
                feature=feature_name,
                node_status=str(payload.get("status", "") or ""),
                details=str(payload.get("details", "") or ""),
                blocked_by=[str(item) for item in (payload.get("blocked_by") or []) if str(item).strip()],
            )
        return optional_payload

    def _write_run_tracking_artifacts(
        self,
        *,
        run_started_at: str,
        run_status: str,
        error_message: str,
        step_states: List[Dict[str, object]],
    ) -> None:
        self.run_tracking_dir.mkdir(parents=True, exist_ok=True)
        pruned_empty_dirs = 0

        output_rows: List[Dict[str, object]] = []
        for file_path in sorted(self.output_dir.rglob("*")):
            if not file_path.is_file():
                continue
            relative = file_path.relative_to(self.output_dir)
            if relative.parts and relative.parts[0] == "run_tracking":
                continue
            category = self._classify_output_path(relative)
            suffix = file_path.suffix.lower()
            output_rows.append(
                {
                    "relative_path": relative.as_posix(),
                    "category": category,
                    "extension": suffix,
                    "size_bytes": int(file_path.stat().st_size),
                    "modified_utc": datetime.fromtimestamp(file_path.stat().st_mtime, tz=timezone.utc).isoformat(),
                }
            )
        outputs_index_df = pd.DataFrame(
            output_rows,
            columns=_RUN_TRACKING_OUTPUT_COLUMNS,
        )
        outputs_index_csv = self.run_tracking_dir / "outputs_index.csv"
        outputs_index_json = self.run_tracking_dir / "outputs_index.json"
        outputs_index_df.to_csv(outputs_index_csv, index=False)
        outputs_index_df.to_json(outputs_index_json, orient="records", indent=2)

        step_rows: List[Dict[str, object]] = []
        for index, state in enumerate(step_states, start=1):
            status = str(state.get("status", "") or "")
            details = str(state.get("details", "") or "")
            error = str(state.get("error", "") or "")
            if not error and status == "failed" and details:
                error = details
            step_rows.append(
                {
                    "index": int(state.get("index", index) or index),
                    "step": str(state.get("step", "") or ""),
                    "required": bool(state.get("required", True)),
                    "status": status,
                    "started_at": str(state.get("started_at", "") or ""),
                    "ended_at": str(state.get("ended_at", "") or ""),
                    "details": details,
                    "error": error,
                }
            )

        steps_df = pd.DataFrame(
            step_rows,
            columns=_RUN_TRACKING_STEP_COLUMNS,
        )
        step_status_csv = self.run_tracking_dir / "step_status.csv"
        step_status_json = self.run_tracking_dir / "step_status.json"
        steps_df.to_csv(step_status_csv, index=False)
        steps_df.to_json(step_status_json, orient="records", indent=2)

        status_counts: Dict[str, int] = {}
        for row in step_rows:
            status = str(row.get("status", "") or "")
            if not status:
                continue
            status_counts[status] = status_counts.get(status, 0) + 1

        optional_features = self._infer_optional_feature_statuses(step_rows)

        run_manifest = {
            "pipeline": "multi_engine_unified",
            "run_id": self.run_id,
            "input": {
                "project_dir": str(self.project_dir),
                "pairlist_file": str(pairlist_path(self.project_dir)),
            },
            "project_dir": str(self.project_dir),
            "output_dir": str(self.output_dir),
            "analysis_mode": self.analysis_mode,
            "analysis_scope": self.analysis_scope,
            "consensus_mode": self.consensus_mode,
            "normalization_method": self.normalization_method,
            "rescoring_scope": self.rescoring_scope,
            "rescoring_top_n": int(self.rescoring_top_n),
            "top_pose_selection_policy": self.top_pose_selection_policy,
            "top_pose_global_aggregation": self.top_pose_global_aggregation,
            "favorite_engine": str(self.favorite_engine or self.manifest.get("favorite_engine") or self.engine or ""),
            "engines": list(self.manifest.get("engines", [])),
            "engines_in_scope": list(self.engines_in_scope),
            "excluded_engines": list(self.excluded_engines),
            "scope_source": self.scope_source,
            "engine_preset_name": self.engine_preset_name,
            "scoped_engine_count": int(len(self.engines_in_scope)),
            "detected_engine_count": int(self.detected_engine_count or len(self.manifest.get("engines", []))),
            "run_started_at": run_started_at,
            "run_completed_at": self._utc_now_iso(),
            "run_status": run_status,
            "error_message": error_message,
            "stage_contract": self._infer_stage_contract(step_rows),
            "step_status_file": str(step_status_csv),
            "outputs_index_file": str(outputs_index_csv),
            "steps_total": int(len(step_rows)),
            "step_status_counts": status_counts,
            "outputs_count": int(len(outputs_index_df)),
            "optional_features": optional_features,
            "pruned_empty_directories": int(pruned_empty_dirs),
            "sqlite_dual_write_enabled": bool(self.sqlite_dual_write_enabled),
            "sqlite_write_summary": dict(self.sqlite_write_summary or {}),
        }
        (self.run_tracking_dir / "run_manifest.json").write_text(
            json.dumps(run_manifest, indent=2),
            encoding="utf-8",
        )
        pruned_empty_dirs = self._prune_empty_directories(self.output_dir)
        if pruned_empty_dirs:
            run_manifest["pruned_empty_directories"] = int(pruned_empty_dirs)
            (self.run_tracking_dir / "run_manifest.json").write_text(
                json.dumps(run_manifest, indent=2),
                encoding="utf-8",
            )

    def run(self) -> bool:
        run_started_at = self._utc_now_iso()
        run_status = "failed"
        error_message = ""
        active_step = ""
        step_order: List[str] = []
        step_state_map: Dict[str, Dict[str, object]] = {}

        def _set_step(step: str, status: str, details: str = "") -> None:
            normalized = status if status in _RUN_TRACKING_STATES else "needs_review"
            now = self._utc_now_iso()
            if step not in step_state_map:
                step_state_map[step] = {
                    "step": step,
                    "status": "not_started",
                    "started_at": "",
                    "ended_at": "",
                    "details": "",
                }
                step_order.append(step)
            record = step_state_map[step]
            if normalized == "in_progress" and not str(record.get("started_at") or ""):
                record["started_at"] = now
            if normalized in {"completed", "validated", "failed", "skipped", "needs_review"}:
                if not str(record.get("started_at") or ""):
                    record["started_at"] = now
                record["ended_at"] = now
            record["status"] = normalized
            if details:
                record["details"] = details

        def _start_step(step: str) -> None:
            nonlocal active_step
            active_step = step
            _set_step(step, "in_progress")

        def _complete_step(step: str, details: str = "") -> None:
            nonlocal active_step
            _set_step(step, "completed", details)
            active_step = ""

        try:
            _start_step("load_or_build_scores")
            scores = self._load_or_build_scores()
            _complete_step("load_or_build_scores", f"score_rows={len(scores)}")

            if self.complex_query:
                _start_step("apply_complex_query")
                original_score_count = len(scores)
                scores = filter_frame_by_complex_query(scores, self.complex_query)
                if scores.empty:
                    raise ValueError("Complex query removed all comparative score rows")
                score_tags = tags_from_frame(scores)
                self.selected_query_tags = self.selected_query_tags.intersection(score_tags) if self.selected_query_tags else score_tags
                if self.complex_query and not self.selected_query_tags:
                    raise ValueError("Complex query matched canonical pairs, but none of those tags had usable engine scores")
                if self.selected_query_tags and not self.pairlist_df.empty:
                    pair_tags = tags_from_frame(self.pairlist_df)
                    keep_tags = pair_tags.intersection(self.selected_query_tags)
                    if not keep_tags:
                        raise ValueError("Complex query selected engine score rows, but no canonical pairlist tags overlap after filtering")
                    pair_tag_series = self.pairlist_df.apply(
                        lambda row: f"{row['receptor']}_{row['site_id']}_{row['ligand']}",
                        axis=1,
                    ).astype(str)
                    self.pairlist_df = self.pairlist_df.loc[pair_tag_series.isin(keep_tags)].copy()
                (self.output_dir / "complex_query_summary.txt").write_text(
                    (
                        f"query={self.complex_query}\n"
                        f"filtered_pairlist_rows={len(self.pairlist_df)}\n"
                        f"filtered_score_rows={len(scores)}\n"
                        f"original_score_rows={original_score_count}\n"
                    ),
                    encoding="utf-8",
                )
                _complete_step("apply_complex_query", f"filtered_score_rows={len(scores)}")
            else:
                _set_step("apply_complex_query", "skipped", "No complex query provided")

            if scores.empty:
                raise ValueError(
                    "No engine scores were available for analysis. Run docking to produce raw outputs or place normalized_scores.csv files under engines/<engine>/scores/."
                )

            _start_step("apply_problematic_ligand_filters")
            before_filter_count = len(scores)
            scores = self._apply_problematic_ligand_filters(scores)
            if scores.empty:
                raise ValueError(
                    "All engine scores were filtered out by problematic-ligand exclusion. "
                    "Disable exclusion or loosen thresholds."
                )
            _complete_step(
                "apply_problematic_ligand_filters",
                f"score_rows_before={before_filter_count};score_rows_after={len(scores)};enabled={self.exclude_problematic_ligands}",
            )

            if self.analysis_mode == "single_engine":
                _set_step("write_comparative_reports", "skipped", "analysis_mode=single_engine")
                _set_step("promote_exhaustive_rerun_manifest", "skipped", "analysis_mode=single_engine")
                _start_step("single_engine_reports")
                target_engine = self.engine or self.favorite_engine
                if not target_engine:
                    raise ValueError("--engine is required for single_engine analysis")
                self._write_single_engine_reports(scores, target_engine, self.output_dir / target_engine)
                _complete_step("single_engine_reports", f"engine={target_engine}")
                _set_step("favorite_engine_reports", "skipped", "analysis_mode=single_engine")
                _set_step("favorite_engine_downstream_bridge", "skipped", "analysis_mode=single_engine")
                run_status = "completed"
                return True

            _start_step("write_comparative_reports")
            self._write_comparative_reports(scores, analysis_scope=self.analysis_scope)
            _complete_step("write_comparative_reports", f"analysis_scope={self.analysis_scope}")

            if self.promote_exhaustive:
                _start_step("promote_exhaustive_rerun_manifest")
                self.rerun_manifest_file = self._write_exhaustive_rerun_manifest(scores)
                _complete_step(
                    "promote_exhaustive_rerun_manifest",
                    f"rerun_manifest_file={self.rerun_manifest_file}",
                )
            else:
                _set_step("promote_exhaustive_rerun_manifest", "skipped", "promote_exhaustive=False")

            if self.analysis_mode == "comparative_all_engines":
                _set_step("single_engine_reports", "skipped", "analysis_mode=comparative_all_engines")
                _set_step("favorite_engine_reports", "skipped", "analysis_mode=comparative_all_engines")
                _set_step("favorite_engine_downstream_bridge", "skipped", "analysis_mode=comparative_all_engines")
                run_status = "completed"
                return True

            _set_step("single_engine_reports", "skipped", "analysis_mode=favorite_engine_continue")
            _start_step("favorite_engine_reports")
            favorite_engine = str(
                self.favorite_engine or self.manifest.get("favorite_engine") or self.engine or ""
            ).strip().lower()
            if not favorite_engine:
                raise ValueError(
                    "favorite_engine_continue requires --favorite-engine, --engine, or project_manifest.json.favorite_engine"
                )
            if self.engines_in_scope and favorite_engine not in self.engines_in_scope:
                raise ValueError(
                    f"favorite_engine '{favorite_engine}' is not in engines_in_scope "
                    f"{self.engines_in_scope}. Use --engine to specify a valid engine."
                )
            if not self.engines_in_scope:
                detection_report = (
                    getattr(self, "engine_detection_report", {})
                    if isinstance(getattr(self, "engine_detection_report", {}), dict)
                    else {}
                )
                detection_valid_engines = [
                    str(item).strip().lower()
                    for item in (detection_report.get("valid_engines") or [])
                    if str(item).strip()
                ]
                if detection_valid_engines and favorite_engine not in detection_valid_engines:
                    raise ValueError(
                        f"favorite_engine '{favorite_engine}' is not in detected valid engines "
                        f"{detection_valid_engines}. Use --engine to specify a valid engine."
                    )
            # Avoid path doubling when output_dir was already rooted inside
            # a "favorite_engine" directory (e.g. a re-run from an existing session).
            if self.output_dir.name == "favorite_engine":
                favorite_dir = self.output_dir
            else:
                favorite_dir = self.output_dir / "favorite_engine"
            self._write_single_engine_reports(scores, favorite_engine, favorite_dir)
            _complete_step("favorite_engine_reports", f"engine={favorite_engine}")

            if self.analysis_scope in {"full", "rescoring_only"}:
                _start_step("favorite_engine_downstream_bridge")
                self._run_downstream_bridge(
                    favorite_engine,
                    scores[scores["engine"] == favorite_engine].copy(),
                    favorite_dir,
                )
                _complete_step("favorite_engine_downstream_bridge", f"analysis_scope={self.analysis_scope}")
            else:
                (favorite_dir / "downstream_bridge.txt").write_text(
                    (
                        "Downstream bridge skipped by analysis scope.\n"
                        f"analysis_scope={self.analysis_scope}\n"
                    ),
                    encoding="utf-8",
                )
                _set_step(
                    "favorite_engine_downstream_bridge",
                    "skipped",
                    f"analysis_scope={self.analysis_scope}",
                )
            run_status = "completed"
            return True
        except Exception as exc:
            error_message = str(exc)
            if active_step:
                _set_step(active_step, "failed", error_message)
            run_status = "failed"
            raise
        finally:
            _set_step("finalize_run_tracking", "in_progress")
            step_states = [step_state_map[name] for name in step_order]
            try:
                self._write_run_tracking_artifacts(
                    run_started_at=run_started_at,
                    run_status=run_status,
                    error_message=error_message,
                    step_states=step_states,
                )
                _set_step("finalize_run_tracking", "validated")
                step_states = [step_state_map[name] for name in step_order]
                self._write_run_tracking_artifacts(
                    run_started_at=run_started_at,
                    run_status=run_status,
                    error_message=error_message,
                    step_states=step_states,
                )
            except Exception as tracking_exc:
                _set_step("finalize_run_tracking", "needs_review", str(tracking_exc))
                logger.warning("Run tracking finalization failed: %s", tracking_exc)

    @staticmethod
    def _source_fingerprint(engine_layout: Dict[str, Path], pair_index: Dict[str, object]) -> str:
        digest = hashlib.sha256(json.dumps(pair_index, sort_keys=True, default=str).encode())
        source_paths = [engine_layout["poses"], engine_layout["logs"], engine_layout["scores"] / "all_scores.csv"]
        for source in source_paths:
            source = Path(source)
            digest.update(str(source).encode())
            candidates = sorted(source.rglob("*")) if source.is_dir() else [source]
            for path in candidates:
                if path.is_file():
                    digest.update(str(path).encode())
                    digest.update(content_hash(path).encode())
                elif not path.exists():
                    digest.update(b"missing")
        for code_file in sorted(Path(__file__).parent.glob("*.py")):
            digest.update(content_hash(code_file).encode())
        parser_path = Path(__file__).parents[1] / "docking" / "runners" / "job_contract.py"
        digest.update(content_hash(parser_path).encode())
        return digest.hexdigest()

    @staticmethod
    def _pose_is_usable(path: Path) -> bool:
        if not path.is_file() or path.stat().st_size == 0:
            return False
        completion = Path(str(path) + ".completion.json")
        if completion.exists():
            try:
                payload = json.loads(completion.read_text(encoding="utf-8"))
                if payload.get("status") not in {"completed", "success"}:
                    return False
                expected_hash = payload.get("pose_sha256") or payload.get("pose_hash")
                if expected_hash and expected_hash != content_hash(path):
                    return False
            except (OSError, ValueError):
                return False
        return True

    def _load_or_build_scores(self) -> pd.DataFrame:
        frames: List[pd.DataFrame] = []
        self.pairlist_df = load_pairlist(self.project_dir)
        pair_index = self._pair_index()
        candidate_engines = self.engines_in_scope or [
            str(engine).strip().lower()
            for engine in (self.manifest.get("engines", []) or [])
            if str(engine).strip()
        ]
        for engine in candidate_engines:
            engine_layout = ensure_engine_layout(self.project_dir, engine)
            normalized_path = engine_layout["scores"] / "normalized_scores.csv"
            raw_pose_files = [path for path in engine_layout["poses"].glob("*") if path.suffix.lower() in {".sdf", ".pdbqt", ".dlg"}]
            imported = read_explicit_score_import(engine_layout["scores"]) if not raw_pose_files else None
            if imported is not None:
                if not imported["engine"].eq(engine).all():
                    raise ValueError(f"Imported score engine must match folder engine {engine}")
                frame = self._merge_pair_metadata(imported, self._pair_metadata_frame())
                for column in NORMALIZED_COLUMNS:
                    if column not in frame:
                        frame[column] = None
                frame["scoring_function"] = frame["scoring_function"].fillna(engine)
                if pair_index:
                    frame = frame[frame["tag"].isin(pair_index)]
                frames.append(frame[NORMALIZED_COLUMNS].copy())
                continue
            _use_cached = False
            cache_key_file = normalized_path.with_suffix(".sources.json")
            fingerprint = self._source_fingerprint(engine_layout, pair_index)
            if normalized_path.exists() and cache_key_file.exists():
                try:
                    saved_key = json.loads(cache_key_file.read_text(encoding="utf-8"))
                    _use_cached = saved_key.get("source") == fingerprint and saved_key.get("output") == content_hash(normalized_path)
                except (ValueError, OSError):
                    pass
            if _use_cached:
                frame = pd.read_csv(normalized_path)
            else:
                frame = self._build_normalized_frame(engine, engine_layout, pair_index)
                if fingerprint != self._source_fingerprint(engine_layout, pair_index):
                    raise RuntimeError("Docking inputs changed during score import; rerun from stable completed outputs")
                if not frame.empty:
                    frame.to_csv(normalized_path, index=False)
                    cache_key_file.write_text(json.dumps({"source": fingerprint, "output": content_hash(normalized_path)}), encoding="utf-8")
            if not frame.empty:
                for column in NORMALIZED_COLUMNS:
                    if column not in frame.columns:
                        frame[column] = None
                frame = self._merge_pair_metadata(frame, self._pair_metadata_frame())
                if "scoring_function" not in frame:
                    frame["scoring_function"] = engine
                frame["scoring_function"] = frame["scoring_function"].fillna(engine)
                if "receptor_frame_id" not in frame:
                    frame["receptor_frame_id"] = None
                for protein in frame["protein"].dropna().unique():
                    mask = frame["protein"].eq(protein) & frame["receptor_frame_id"].isna()
                    receptor = shared_receptors_dir(self.project_dir) / f"{Path(str(protein)).stem}.pdbqt"
                    if receptor.is_file():
                        frame.loc[mask, "receptor_frame_id"] = "prepared-state:" + content_hash(receptor)
                frame = frame[frame["pose_file"].map(lambda value: self._pose_is_usable(Path(str(value))))]
                frames.append(frame[NORMALIZED_COLUMNS].copy())
        if not frames:
            return pd.DataFrame(columns=NORMALIZED_COLUMNS)
        return self._add_unified_compatibility_columns(pd.concat(frames, ignore_index=True))

    def _apply_problematic_ligand_filters(self, scores: pd.DataFrame) -> pd.DataFrame:
        if not self.exclude_problematic_ligands or scores.empty:
            return scores
        if "tag" not in scores.columns:
            return scores

        frame = scores.copy()
        frame["tag"] = frame["tag"].astype(str)
        frame["affinity_kcal_mol"] = pd.to_numeric(frame.get("affinity_kcal_mol"), errors="coerce")
        best_affinity = frame.groupby("tag", dropna=False)["affinity_kcal_mol"].min()
        pose_counts = frame.groupby("tag", dropna=False).size()

        score_tags = set(frame["tag"].dropna().astype(str))
        pair_tags = tags_from_frame(self.pairlist_df) if not self.pairlist_df.empty else set()
        all_tags = pair_tags.union(score_tags)

        positive_tags = {
            str(tag)
            for tag, value in best_affinity.items()
            if pd.notna(value) and float(value) > self.positive_affinity_threshold
        }
        low_pose_tags = {
            str(tag)
            for tag, count in pose_counts.items()
            if int(count) < self.minimum_pose_count
        }
        missing_pose_tags = {str(tag) for tag in all_tags.difference(score_tags)}

        excluded_tags = positive_tags.union(low_pose_tags).union(missing_pose_tags)
        if not excluded_tags:
            return frame

        keep_tags = all_tags.difference(excluded_tags)
        filtered = frame[frame["tag"].isin(keep_tags)].copy()
        if not self.pairlist_df.empty:
            pair_tag_series = self.pairlist_df.apply(
                lambda row: f"{row['receptor']}_{row['site_id']}_{row['ligand']}",
                axis=1,
            ).astype(str)
            self.pairlist_df = self.pairlist_df.loc[pair_tag_series.isin(keep_tags)].copy()

        exclusion_rows = []
        for tag in sorted(excluded_tags):
            reasons = []
            if tag in positive_tags:
                reasons.append(f"best_affinity>{self.positive_affinity_threshold}")
            if tag in low_pose_tags:
                reasons.append(f"pose_count<{self.minimum_pose_count}")
            if tag in missing_pose_tags:
                reasons.append("no_pose_rows")
            exclusion_rows.append({"tag": tag, "reason": ";".join(reasons) or "excluded"})
        exclusion_df = pd.DataFrame(exclusion_rows)
        exclusion_path = self.output_dir / "excluded_problematic_ligands.csv"
        exclusion_df.to_csv(exclusion_path, index=False)
        logger.info(
            f"🚫 Excluded {len(excluded_tags)} problematic ligands before multi-engine analysis "
            f"(report: {exclusion_path})"
        )
        return filtered

    def _add_unified_compatibility_columns(self, scores: pd.DataFrame) -> pd.DataFrame:
        frame = scores.copy()
        frame["pose"] = pd.to_numeric(frame["pose"], errors="coerce").fillna(0).astype(int)
        frame["affinity_kcal_mol"] = pd.to_numeric(frame["affinity_kcal_mol"], errors="coerce")
        frame["complex_name"] = frame["tag"]
        frame["protein_label"] = frame["protein"]
        frame["mode"] = frame["pose"]
        frame["vina_affinity"] = frame["affinity_kcal_mol"]

        frame["cnn_affinity"] = None
        primary_cnn_affinity = frame["score_name_primary"].fillna("").eq("cnn_affinity")
        secondary_cnn_affinity = frame["score_name_secondary"].fillna("").eq("cnn_affinity")
        frame.loc[primary_cnn_affinity, "cnn_affinity"] = frame.loc[primary_cnn_affinity, "score_primary"]
        frame.loc[secondary_cnn_affinity, "cnn_affinity"] = frame.loc[secondary_cnn_affinity, "score_secondary"]

        frame["cnn_score"] = None
        primary_cnn_score = frame["score_name_primary"].fillna("").eq("cnn_score")
        secondary_cnn_score = frame["score_name_secondary"].fillna("").eq("cnn_score")
        frame.loc[primary_cnn_score, "cnn_score"] = frame.loc[primary_cnn_score, "score_primary"]
        frame.loc[secondary_cnn_score, "cnn_score"] = frame.loc[secondary_cnn_score, "score_secondary"]

        for column in UNIFIED_COMPAT_COLUMNS:
            if column not in frame.columns:
                frame[column] = None
        return frame

    def _pair_index(self) -> Dict[str, Dict[str, object]]:
        index: Dict[str, Dict[str, object]] = {}
        for _, row in self.pairlist_df.iterrows():
            tag = f"{row['receptor']}_{row['site_id']}_{row['ligand']}"
            payload: Dict[str, object] = {
                "protein": row["receptor"],
                "ligand": row["ligand"],
                "site_id": row["site_id"],
            }
            for optional in _PAIR_METADATA_COLUMNS:
                if optional in row:
                    payload[optional] = row.get(optional)
            index[tag] = payload
        return index

    def _pair_metadata_frame(self) -> pd.DataFrame:
        if self.pairlist_df is None or self.pairlist_df.empty:
            return pd.DataFrame(columns=["tag", *_PAIR_METADATA_COLUMNS])
        frame = self.pairlist_df.copy()
        if "tag" not in frame.columns:
            frame["tag"] = frame.apply(
                lambda row: f"{row['receptor']}_{row['site_id']}_{row['ligand']}",
                axis=1,
            ).astype(str)
        keep = ["tag"] + [column for column in _PAIR_METADATA_COLUMNS if column in frame.columns]
        meta = frame[keep].copy().drop_duplicates("tag")
        if "is_cocrystal_benchmark" in meta.columns:
            meta["is_cocrystal_benchmark"] = meta["is_cocrystal_benchmark"].map(explicit_true)
        return meta

    @staticmethod
    def _merge_pair_metadata(frame: pd.DataFrame, pair_meta: pd.DataFrame) -> pd.DataFrame:
        if frame is None or frame.empty or pair_meta is None or pair_meta.empty:
            return frame
        merged = frame.merge(pair_meta, on="tag", how="left", suffixes=("", "_pair"))
        for column in pair_meta.columns:
            if column == "tag":
                continue
            alt = f"{column}_pair"
            if alt not in merged.columns:
                continue
            if column in merged.columns:
                merged[column] = merged[column].where(~merged[column].isna(), merged[alt])
            else:
                merged[column] = merged[alt]
            merged.drop(columns=[alt], inplace=True, errors="ignore")
        return merged

    def _build_normalized_frame(
        self,
        engine: str,
        engine_layout: Dict[str, Path],
        pair_index: Dict[str, Dict[str, object]],
    ) -> pd.DataFrame:
        if engine == "gnina":
            source_rows = []
            for pose_file in sorted(engine_layout["poses"].glob("*.sdf")):
                tag = pose_file.stem
                if (pair_index and tag not in pair_index) or not self._pose_is_usable(pose_file):
                    continue
                try:
                    parsed = parse_sdf_scores(pose_file)
                    source_rows.extend({"tag": tag, "mode": item["pose"], "vina_affinity": item["affinity"],
                                        "cnn_affinity": item.get("cnn_affinity"), "cnn_score": item.get("cnn_score")}
                                       for item in parsed)
                except (ValueError, KeyError, OSError):
                    # Legacy imports may lack SDF properties. Require a real pose
                    # and exact requested tag before accepting its engine log.
                    log_file = engine_layout["logs"] / f"{tag}.log"
                    if log_file.is_file():
                        parsed, _ = parse_gnina_log(log_file, {key:key for key in pair_index} if pair_index else None)
                        source_rows.extend(parsed)
            source = pd.DataFrame(source_rows)
            rows = []
            for _, record in source.iterrows():
                tag = str(record.get("tag", ""))
                if pair_index and tag not in pair_index:
                    continue
                affinity = pd.to_numeric(record.get("vina_affinity"), errors="coerce")
                if pd.isna(affinity):
                    continue
                pose_value = pd.to_numeric(record.get("mode", 0), errors="coerce")
                pose = int(pose_value) if pd.notna(pose_value) else 0
                pair = pair_index.get(tag, {})
                rows.append(
                    {
                        "engine": engine,
                        "tag": tag,
                        "protein": pair.get("protein", ""),
                        "ligand": pair.get("ligand", ""),
                        "site_id": pair.get("site_id", ""),
                        "pose": pose,
                        "affinity_kcal_mol": float(affinity),
                        "score_name_primary": "cnn_affinity",
                        "score_primary": float(record.get("cnn_affinity")) if pd.notna(record.get("cnn_affinity")) else None,
                        "score_name_secondary": "cnn_score",
                        "score_secondary": float(record.get("cnn_score")) if pd.notna(record.get("cnn_score")) else None,
                        "rmsd_lb": None,
                        "rmsd_ub": None,
                        "pose_file": str(engine_layout["poses"] / f"{tag}.sdf"),
                        "log_file": str(engine_layout["logs"] / f"{tag}.log"),
                        "pair_source": pair.get("pair_source"),
                        "selection_mode": pair.get("selection_mode"),
                        "is_cocrystal_benchmark": pair.get("is_cocrystal_benchmark"),
                        "cocrystal_ligand_name": pair.get("cocrystal_ligand_name"),
                        "reference_pose_file": pair.get("reference_pose_file"),
                        "reference_ligand_file": pair.get("reference_ligand_file"),
                    }
                )
            return pd.DataFrame(rows)

        rows = []
        if engine == "autodock4":
            for pose_file in sorted(engine_layout["poses"].glob("*.dlg")):
                tag = pose_file.stem
                if pair_index and tag not in pair_index:
                    continue
                parsed = parse_autodock4_dlg(pose_file)
                pair = pair_index.get(tag, {})
                for _, record in parsed.iterrows():
                    affinity = pd.to_numeric(record.get("autodock4_affinity"), errors="coerce")
                    if pd.isna(affinity):
                        continue
                    pose_value = pd.to_numeric(record.get("pose", 0), errors="coerce")
                    pose = int(pose_value) if pd.notna(pose_value) else 0
                    rows.append(
                        {
                            "engine": engine,
                            "tag": tag,
                            "protein": pair.get("protein", ""),
                            "ligand": pair.get("ligand", ""),
                            "site_id": pair.get("site_id", ""),
                            "pose": pose,
                            "affinity_kcal_mol": float(affinity),
                            "score_name_primary": "autodock4_affinity",
                            "score_primary": float(affinity),
                            "score_name_secondary": "",
                            "score_secondary": None,
                            "rmsd_lb": None,
                            "rmsd_ub": None,
                            "pose_file": str(pose_file),
                            "log_file": str(engine_layout["logs"] / f"{tag}.log"),
                            "pair_source": pair.get("pair_source"),
                            "selection_mode": pair.get("selection_mode"),
                            "is_cocrystal_benchmark": pair.get("is_cocrystal_benchmark"),
                            "cocrystal_ligand_name": pair.get("cocrystal_ligand_name"),
                            "reference_pose_file": pair.get("reference_pose_file"),
                            "reference_ligand_file": pair.get("reference_ligand_file"),
                        }
                    )
            return pd.DataFrame(rows)

        for pose_file in sorted(engine_layout["poses"].glob("*.pdbqt")):
            tag = pose_file.stem
            if pair_index and tag not in pair_index:
                continue
            parsed = parse_vina_pdbqt(pose_file)
            pair = pair_index.get(tag, {})
            for _, record in parsed.iterrows():
                affinity = pd.to_numeric(record.get("vina_affinity"), errors="coerce")
                if pd.isna(affinity):
                    continue
                pose_value = pd.to_numeric(record.get("pose", 0), errors="coerce")
                pose = int(pose_value) if pd.notna(pose_value) else 0
                rows.append(
                    {
                        "engine": engine,
                        "tag": tag,
                        "protein": pair.get("protein", ""),
                        "ligand": pair.get("ligand", ""),
                        "site_id": pair.get("site_id", ""),
                        "pose": pose,
                        "affinity_kcal_mol": float(affinity),
                        "score_name_primary": "vina_affinity",
                        "score_primary": float(affinity),
                        "score_name_secondary": "",
                        "score_secondary": None,
                        "rmsd_lb": float(record.get("rmsd_lb")) if pd.notna(record.get("rmsd_lb")) else None,
                        "rmsd_ub": float(record.get("rmsd_ub")) if pd.notna(record.get("rmsd_ub")) else None,
                        "pose_file": str(pose_file),
                        "log_file": str(engine_layout["logs"] / f"{tag}.log"),
                        "pair_source": pair.get("pair_source"),
                        "selection_mode": pair.get("selection_mode"),
                        "is_cocrystal_benchmark": pair.get("is_cocrystal_benchmark"),
                        "cocrystal_ligand_name": pair.get("cocrystal_ligand_name"),
                        "reference_pose_file": pair.get("reference_pose_file"),
                        "reference_ligand_file": pair.get("reference_ligand_file"),
                    }
                )
        return pd.DataFrame(rows)

    def _write_comparative_reports(self, scores: pd.DataFrame, analysis_scope: str = "full", *, preserve_dag_atlas: bool = False) -> None:
        scope = str(analysis_scope or "full").strip().lower()
        if scope not in _SUPPORTED_ANALYSIS_SCOPES:
            scope = "full"

        numbered_layout = ensure_numbered_output_layout(self.project_dir)
        reports_dir = self.output_dir / "reports"
        reports_dir.mkdir(parents=True, exist_ok=True)
        raw_dir = self.output_dir / "raw_data"
        raw_dir.mkdir(parents=True, exist_ok=True)
        numbered_scores_raw = numbered_layout["post_scores_raw"]
        numbered_scores_unified = numbered_layout["post_scores_unified"]
        numbered_scores_consensus = numbered_layout["post_scores_consensus"]
        numbered_reports_root = numbered_layout["reports_root_numbered"]

        normalized_scores = normalize_engine_scores(
            scores,
            method=self.normalization_method,
            score_column="affinity_kcal_mol",
            normalized_column="normalized_affinity_score",
            group_keys=["engine", "protein"],
        )
        pair_metadata = self._pair_metadata_frame()
        normalized_scores = self._merge_pair_metadata(normalized_scores, pair_metadata)
        normalized_scores = self._annotate_scope_columns(normalized_scores)
        normalized_scores.to_csv(reports_dir / "combined_engine_scores.csv", index=False)
        normalized_scores.to_csv(reports_dir / "consensus_inputs_all_engines.csv", index=False)
        normalized_scores[UNIFIED_COMPAT_COLUMNS].to_csv(raw_dir / "unified_all_scores.csv", index=False)

        best_by_engine = self._best_rows_by_group(normalized_scores, ["engine", "tag"], "affinity_kcal_mol")
        best_by_engine = self._merge_pair_metadata(best_by_engine, pair_metadata)
        best_by_engine = self._annotate_scope_columns(best_by_engine)
        best_by_engine.sort_values(["engine", "affinity_kcal_mol", "tag"], inplace=True)
        best_by_engine.to_csv(reports_dir / "best_pose_per_tag_by_engine.csv", index=False)
        best_by_engine[UNIFIED_COMPAT_COLUMNS].to_csv(raw_dir / "unified_best_poses.csv", index=False)

        engine_summary = (
            best_by_engine.groupby("engine")
            .agg(
                complexes=("tag", "nunique"),
                best_affinity=("affinity_kcal_mol", "min"),
                mean_best_affinity=("affinity_kcal_mol", "mean"),
                median_best_affinity=("affinity_kcal_mol", "median"),
            )
            .reset_index()
        )
        engine_summary = self._annotate_scope_columns(engine_summary)
        engine_summary.to_csv(reports_dir / "engine_summary.csv", index=False)

        relative_best = normalize_engine_scores(best_by_engine, group_keys=["engine", "protein", "site_id"])
        cross_engine_best = relative_best.sort_values(["tag", "normalized_affinity_score", "engine"], ascending=[True, False, True]).drop_duplicates("tag").copy()
        cross_engine_best["interpretation"] = "relative_engine_rank_not_raw_energy_winner"
        cross_engine_best = self._annotate_scope_columns(cross_engine_best)
        cross_engine_best.to_csv(reports_dir / "best_engine_per_complex.csv", index=False)

        matrix = best_by_engine.pivot_table(
            index="tag",
            columns="engine",
            values="affinity_kcal_mol",
            aggfunc="min",
        ).reset_index()
        matrix = self._annotate_scope_columns(matrix)
        matrix.to_csv(reports_dir / "engine_affinity_matrix.csv", index=False)

        validation_bundle = run_redocking_validation(
            project_dir=self.project_dir,
            best_by_engine=best_by_engine,
            expected_reference_rows=self.pairlist_df,
            expected_engines=self.engines_in_scope,
            output_dir=reports_dir,
        )
        redocking_validation_df = validation_bundle.get("validation_df", pd.DataFrame())
        reference_baselines_df = validation_bundle.get("reference_baselines_df", pd.DataFrame())
        redocking_validation_file = Path(str(validation_bundle.get("validation_file", reports_dir / "redocking_validation.csv")))
        validation_summary_file = Path(str(validation_bundle.get("summary_file", reports_dir / "VALIDATION_SUMMARY.md")))
        validation_gate_file = Path(
            str(validation_bundle.get("validation_gate_file", reports_dir / "validation_gate_status.json"))
        )
        reference_baselines_file = Path(
            str(validation_bundle.get("reference_baselines_file", reports_dir / "reference_baselines.csv"))
        )
        validation_gate = dict(validation_bundle.get("validation_gate") or validation_bundle.get("summary") or {})
        allow_reference_anchor = bool(validation_gate.get("allow_reference_anchor", False))
        effective_hit_class_policy = str(self.hit_class_policy or "target_percentile")
        hit_class_policy_note = ""
        if effective_hit_class_policy == "reference_anchor" and not allow_reference_anchor:
            effective_hit_class_policy = "target_percentile"
            hit_class_policy_note = (
                f"reference_anchor_requested_but_validation_gate_{validation_gate.get('validation_gate_state', 'needs_review')}"
            )
            logger.warning(
                "Reference-anchor classification requested but validation gate is not validated; "
                "falling back to target_percentile policy."
            )

        run_consensus = scope in {"full", "comparison_only", "rescoring_only", "qc_only", "report_only", "top_pose_only"}
        run_visualizations = scope in {"full", "comparison_only"}

        if run_consensus:
            consensus_base = build_consensus_rankings(
                best_by_engine=best_by_engine,
                consensus_mode=self.consensus_mode,
                favorite_engine=self.favorite_engine or self.rerun_engine,
                expected_engines=self.engines_in_scope,
                normalization_method=self.normalization_method,
            )
            consensus_df = classify_hits_target_aware(
                consensus_base,
                policy=effective_hit_class_policy,
                strong_percentile=self.hit_class_strong_percentile,
                moderate_percentile=self.hit_class_moderate_percentile,
                reference_baselines=reference_baselines_df,
            )
            consensus_df = self._annotate_scope_columns(consensus_df)
        else:
            consensus_df = pd.DataFrame()
        consensus_file = reports_dir / "consensus_ranked_hits.csv"
        consensus_class_file = reports_dir / "consensus_ranked_hits_with_classes.csv"

        if run_consensus:
            rescoring_df = select_rescoring_candidates(
                consensus_df=consensus_df,
                rescoring_scope=self.rescoring_scope,
                rescoring_top_n=self.rescoring_top_n,
            )
            rescoring_df = self._annotate_scope_columns(rescoring_df)
        else:
            rescoring_df = pd.DataFrame()
        rescoring_file = reports_dir / "rescoring_candidates.csv"
        rescoring_df.to_csv(rescoring_file, index=False)

        explainability = build_consensus_explainability(
            consensus_df=consensus_df,
            rescoring_df=rescoring_df,
            consensus_mode=self.consensus_mode,
            rescoring_scope=self.rescoring_scope,
            rescoring_top_n=self.rescoring_top_n,
            favorite_engine=self.favorite_engine or self.rerun_engine,
            normalization_method=self.normalization_method,
        )
        explainability_file = reports_dir / "consensus_explainability.json"
        explainability_file.write_text(json.dumps(explainability, indent=2), encoding="utf-8")

        biology_mapping_report: Dict[str, object] = {"mapping_mode": "none", "matched_rows": 0}
        biology_global_corr = pd.DataFrame(
            columns=["score_column", "biology_column", "paired_rows", "spearman_corr", "pearson_corr"]
        )
        biology_per_protein_corr = pd.DataFrame(
            columns=["protein", "score_column", "biology_column", "paired_rows", "spearman_corr", "pearson_corr"]
        )
        annotated_consensus_file = reports_dir / "consensus_ranked_hits_with_biology.csv"
        if run_consensus and self.biology_file:
            try:
                biology_df, biology_info = load_biology_table(self.biology_file)
                annotated_consensus, biology_mapping_report = attach_biology_annotations(
                    consensus_df,
                    biology_df,
                    mapping_mode=self.biology_mapping_mode,
                    biology_prefix="bio_",
                )
                biology_mapping_report["source"] = biology_info
                biology_global_corr, biology_per_protein_corr = compute_biology_correlations(
                    annotated_consensus,
                    score_columns=["consensus_score", "best_affinity_kcal_mol"],
                    biology_prefix="bio_",
                )
                consensus_df = annotated_consensus
                consensus_df.to_csv(annotated_consensus_file, index=False)
                rescoring_df = select_rescoring_candidates(
                    consensus_df=consensus_df,
                    rescoring_scope=self.rescoring_scope,
                    rescoring_top_n=self.rescoring_top_n,
                )
                rescoring_df.to_csv(rescoring_file, index=False)
            except Exception as exc:
                biology_mapping_report = {
                    "mapping_mode": "error",
                    "source": {"path": self.biology_file},
                    "error": str(exc),
                }
        else:
            consensus_df.to_csv(annotated_consensus_file, index=False)

        consensus_df.to_csv(consensus_file, index=False)
        consensus_df.to_csv(consensus_class_file, index=False)
        classification_reports = generate_hit_classification_reports(
            consensus_df=consensus_df,
            output_dir=reports_dir,
            classifier_config={
                "policy": effective_hit_class_policy,
                "strong_percentile": self.hit_class_strong_percentile,
                "moderate_percentile": self.hit_class_moderate_percentile,
                "normalization_method": self.normalization_method,
                "consensus_mode": self.consensus_mode,
            },
            output_stem="consensus_hit_classification",
        )

        biology_mapping_file = reports_dir / "biology_mapping_report.json"
        biology_mapping_file.write_text(json.dumps(biology_mapping_report, indent=2), encoding="utf-8")
        biology_global_corr_file = reports_dir / "biology_correlation_global.csv"
        biology_per_protein_corr_file = reports_dir / "biology_correlation_per_protein.csv"
        biology_global_corr.to_csv(biology_global_corr_file, index=False)
        biology_per_protein_corr.to_csv(biology_per_protein_corr_file, index=False)

        session_polypharm_outputs = self._write_biology_polypharmacology_bundle(
            self.output_dir / "reports" / "polypharmacology",
            consensus_df=consensus_df,
            biology_mapping_report=biology_mapping_report,
            biology_global_corr=biology_global_corr,
            biology_per_protein_corr=biology_per_protein_corr,
        )
        canonical_polypharm_outputs = self._write_biology_polypharmacology_bundle(
            self._canonical_polypharmacology_root(),
            consensus_df=consensus_df,
            biology_mapping_report=biology_mapping_report,
            biology_global_corr=biology_global_corr,
            biology_per_protein_corr=biology_per_protein_corr,
        )

        per_protein_corr, global_corr = compute_cross_engine_rank_correlations(best_by_engine)
        per_protein_corr = self._annotate_scope_columns(per_protein_corr)
        global_corr = self._annotate_scope_columns(global_corr)
        per_protein_corr_file = reports_dir / "engine_rank_correlation_per_protein.csv"
        global_corr_file = reports_dir / "engine_rank_correlation_global.csv"
        per_protein_corr.to_csv(per_protein_corr_file, index=False)
        global_corr.to_csv(global_corr_file, index=False)

        # Canonical AllScore outputs for project-level consumption.
        score_root = ensure_numbered_output_layout(self.project_dir)["post_docking_root_numbered"] / "scores"
        score_raw_root = score_root / "raw"
        score_unified_root = score_root / "unified"
        score_consensus_root = score_root / "consensus"
        for folder in (score_raw_root, score_unified_root, score_consensus_root):
            folder.mkdir(parents=True, exist_ok=True)
        normalized_scores.to_csv(score_raw_root / "all_scores_raw.csv", index=False)
        redocking_validation_df.to_csv(score_raw_root / "redocking_validation.csv", index=False)
        reference_baselines_df.to_csv(score_raw_root / "reference_baselines.csv", index=False)
        best_by_engine.to_csv(score_unified_root / "best_pose_per_tag_by_engine.csv", index=False)
        matrix.to_csv(score_unified_root / "engine_affinity_matrix.csv", index=False)
        engine_summary.to_csv(score_unified_root / "engine_summary.csv", index=False)
        cross_engine_best.to_csv(score_unified_root / "best_engine_per_complex.csv", index=False)
        per_protein_corr.to_csv(score_unified_root / "engine_rank_correlation_per_protein.csv", index=False)
        global_corr.to_csv(score_unified_root / "engine_rank_correlation_global.csv", index=False)
        consensus_df.to_csv(score_consensus_root / "consensus_ranked_hits.csv", index=False)
        consensus_df.to_csv(score_consensus_root / "consensus_ranked_hits_with_classes.csv", index=False)
        rescoring_df.to_csv(score_consensus_root / "rescoring_candidates.csv", index=False)
        explainability_file_target = score_consensus_root / "consensus_explainability.json"
        explainability_file_target.write_text(json.dumps(explainability, indent=2), encoding="utf-8")
        for key in (
            "classification_config_file",
            "classification_distribution_file",
            "classification_qc_breakdown_file",
        ):
            source_path = Path(str(classification_reports.get(key, ""))).expanduser()
            if source_path.exists() and source_path.is_file():
                shutil.copy2(source_path, score_consensus_root / source_path.name)
        (score_consensus_root / "biology_mapping_report.json").write_text(
            json.dumps(biology_mapping_report, indent=2),
            encoding="utf-8",
        )
        if validation_summary_file.exists():
            shutil.copy2(validation_summary_file, score_consensus_root / "VALIDATION_SUMMARY.md")
        if validation_gate_file.exists():
            shutil.copy2(validation_gate_file, score_consensus_root / "validation_gate_status.json")
        biology_global_corr.to_csv(score_consensus_root / "biology_correlation_global.csv", index=False)
        biology_per_protein_corr.to_csv(score_consensus_root / "biology_correlation_per_protein.csv", index=False)

        # Numbered project layout mirrors (4-Working/*).
        normalized_scores.to_csv(numbered_scores_raw / "all_scores_raw.csv", index=False)
        redocking_validation_df.to_csv(numbered_scores_raw / "redocking_validation.csv", index=False)
        reference_baselines_df.to_csv(numbered_scores_raw / "reference_baselines.csv", index=False)
        best_by_engine.to_csv(numbered_scores_unified / "best_pose_per_tag_by_engine.csv", index=False)
        matrix.to_csv(numbered_scores_unified / "engine_affinity_matrix.csv", index=False)
        engine_summary.to_csv(numbered_scores_unified / "engine_summary.csv", index=False)
        cross_engine_best.to_csv(numbered_scores_unified / "best_engine_per_complex.csv", index=False)
        per_protein_corr.to_csv(numbered_scores_unified / "engine_rank_correlation_per_protein.csv", index=False)
        global_corr.to_csv(numbered_scores_unified / "engine_rank_correlation_global.csv", index=False)
        consensus_df.to_csv(numbered_scores_consensus / "consensus_ranked_hits.csv", index=False)
        consensus_df.to_csv(numbered_scores_consensus / "consensus_ranked_hits_with_classes.csv", index=False)
        rescoring_df.to_csv(numbered_scores_consensus / "rescoring_candidates.csv", index=False)
        (numbered_scores_consensus / "consensus_explainability.json").write_text(
            json.dumps(explainability, indent=2),
            encoding="utf-8",
        )
        (numbered_scores_consensus / "biology_mapping_report.json").write_text(
            json.dumps(biology_mapping_report, indent=2),
            encoding="utf-8",
        )
        biology_global_corr.to_csv(numbered_scores_consensus / "biology_correlation_global.csv", index=False)
        biology_per_protein_corr.to_csv(numbered_scores_consensus / "biology_correlation_per_protein.csv", index=False)
        if validation_summary_file.exists():
            shutil.copy2(validation_summary_file, numbered_scores_consensus / "VALIDATION_SUMMARY.md")
        if validation_gate_file.exists():
            shutil.copy2(validation_gate_file, numbered_scores_consensus / "validation_gate_status.json")

        parity_file = numbered_layout["post_storage"] / "csv_sqlite_parity.csv"
        if self.sqlite_dual_write_enabled:
            sqlite_db_file = numbered_layout["post_storage"] / "results.db"
            dataset_frames = {
                "combined_engine_scores": normalized_scores,
                "redocking_validation": redocking_validation_df,
                "reference_baselines": reference_baselines_df,
                "best_pose_per_tag_by_engine": best_by_engine,
                "engine_summary": engine_summary,
                "best_engine_per_complex": cross_engine_best,
                "engine_affinity_matrix": matrix,
                "engine_rank_correlation_per_protein": per_protein_corr,
                "engine_rank_correlation_global": global_corr,
                "consensus_ranked_hits": consensus_df,
                "rescoring_candidates": rescoring_df,
                "biology_correlation_global": biology_global_corr,
                "biology_correlation_per_protein": biology_per_protein_corr,
            }
            source_csv_map = {
                "combined_engine_scores": str(reports_dir / "combined_engine_scores.csv"),
                "redocking_validation": str(redocking_validation_file),
                "reference_baselines": str(reference_baselines_file),
                "best_pose_per_tag_by_engine": str(reports_dir / "best_pose_per_tag_by_engine.csv"),
                "engine_summary": str(reports_dir / "engine_summary.csv"),
                "best_engine_per_complex": str(reports_dir / "best_engine_per_complex.csv"),
                "engine_affinity_matrix": str(reports_dir / "engine_affinity_matrix.csv"),
                "engine_rank_correlation_per_protein": str(per_protein_corr_file),
                "engine_rank_correlation_global": str(global_corr_file),
                "consensus_ranked_hits": str(consensus_file),
                "rescoring_candidates": str(rescoring_file),
                "biology_correlation_global": str(biology_global_corr_file),
                "biology_correlation_per_protein": str(biology_per_protein_corr_file),
            }
            sqlite_summary = write_comparative_bundle(
                sqlite_db_file,
                run_id=self.run_id,
                run_metadata={
                    "generated_at": self._utc_now_iso(),
                    "project_dir": str(self.project_dir),
                    "output_dir": str(self.output_dir),
                    "analysis_mode": self.analysis_mode,
                    "analysis_scope": scope,
                    "consensus_mode": self.consensus_mode,
                    "normalization_method": self.normalization_method,
                    "favorite_engine": str(self.favorite_engine or self.rerun_engine or ""),
                },
                dataset_frames=dataset_frames,
                source_csv_map=source_csv_map,
            )
            self.sqlite_write_summary = sqlite_summary.to_dict()
            parity_df = validate_csv_sqlite_parity(
                sqlite_db_file,
                run_id=self.run_id,
                dataset_csv_map={key: Path(value) for key, value in source_csv_map.items()},
            )
            parity_df.to_csv(parity_file, index=False)
        else:
            self.sqlite_write_summary = {
                "status": "disabled",
                "reason": "feature_flag_enable_sqlite_dual_write_false",
            }
            pd.DataFrame(
                [
                    {
                        "dataset_key": "all",
                        "csv_path": "",
                        "csv_exists": False,
                        "csv_rows": 0,
                        "sqlite_rows": 0,
                        "column_match": False,
                        "status": "skipped",
                        "note": "SQLite dual-write disabled by feature flag",
                    }
                ]
            ).to_csv(parity_file, index=False)

        top_pose_payload = build_top_pose_atlas(
            best_by_engine=best_by_engine,
            consensus_df=consensus_df,
            selection_policy=self.top_pose_selection_policy,
            consensus_mode=self.consensus_mode,
            global_aggregation=self.top_pose_global_aggregation,
            run_id=self.run_id,
            generated_at=self._utc_now_iso(),
            context={
                "analysis_scope": scope,
                "analysis_mode": self.analysis_mode,
                "consensus_mode": self.consensus_mode,
                "normalization_method": self.normalization_method,
                "top_pose_selection_policy": self.top_pose_selection_policy,
                "top_pose_global_aggregation": self.top_pose_global_aggregation,
                "engines_in_scope": list(self.engines_in_scope),
                "scoped_engine_count": int(len(self.engines_in_scope)),
                "detected_engine_count": int(self.detected_engine_count or len(self.manifest.get("engines", []))),
            },
        )
        top_pose_session_root = self.output_dir / "top_pose_ligand_performance"
        top_pose_canonical_root = self._canonical_top_pose_root()
        if preserve_dag_atlas:
            # The atlas node owns these files. A downstream report must consume
            # its exact representative choices and must never rewrite its inputs.
            file_keys = {
                "top_pose_per_ligand_per_protein": "top_pose_per_ligand_per_protein",
                "top_pose_per_ligand_global": "top_pose_per_ligand_global",
                "ligand_performance_summary": "ligand_performance_summary",
                "top_pose_confidence_metrics": "top_pose_confidence_metrics",
            }
            for payload_key, stem in file_keys.items():
                top_pose_payload[payload_key] = pd.read_csv(top_pose_session_root / f"{stem}.csv")
            manifest_path = top_pose_session_root / "top_pose_selection_manifest.json"
            top_pose_payload["manifest"] = json.loads(manifest_path.read_text(encoding="utf-8"))
            session_top_pose_outputs = {f"{key}_file": str(top_pose_session_root / f"{stem}.csv") for key, stem in file_keys.items()}
            session_top_pose_outputs["top_pose_selection_manifest_file"] = str(manifest_path)
        else:
            session_top_pose_outputs = write_top_pose_atlas(top_pose_payload, top_pose_session_root)
        canonical_top_pose_outputs = write_top_pose_atlas(top_pose_payload, top_pose_canonical_root)
        self.top_pose_outputs = {
            **session_top_pose_outputs,
            "top_pose_session_root": str(top_pose_session_root),
            "top_pose_canonical_root": str(top_pose_canonical_root),
            "top_pose_selection_manifest_canonical_file": str(
                canonical_top_pose_outputs.get("top_pose_selection_manifest_file", "")
            ),
            "top_pose_per_ligand_global_canonical_file": str(
                canonical_top_pose_outputs.get("top_pose_per_ligand_global_file", "")
            ),
        }

        summary_lines = [
            "Multi-engine docking summary",
            "============================",
            "",
            f"Project: {self.project_dir}",
            f"Analysis mode: {self.analysis_mode}",
            f"Analysis scope: {scope}",
            f"Consensus mode: {self.consensus_mode}",
            f"Normalization method: {self.normalization_method}",
            f"Hit class policy: {self.hit_class_policy}",
            f"Hit class policy (effective): {effective_hit_class_policy}",
            f"Hit class strong percentile: {self.hit_class_strong_percentile}",
            f"Hit class moderate percentile: {self.hit_class_moderate_percentile}",
            f"Rescoring scope: {self.rescoring_scope}",
            f"Rescoring top-N: {self.rescoring_top_n}",
            f"Top-pose policy: {self.top_pose_selection_policy}",
            f"Top-pose aggregation: {self.top_pose_global_aggregation}",
            f"Engines in scope: {', '.join(self.engines_in_scope) if self.engines_in_scope else 'all detected'}",
            f"Scope source: {self.scope_source}",
            f"Scoped engine count: {len(self.engines_in_scope) if self.engines_in_scope else len(self.manifest.get('engines', []))} of {self.detected_engine_count or len(self.manifest.get('engines', []))}",
            "",
        ]
        if self.excluded_engines:
            summary_lines.append(
                "Excluded engines: "
                + ", ".join(
                    f"{item.get('engine')} ({item.get('exclusion_reason', 'user_excluded')})"
                    for item in self.excluded_engines
                )
            )
            summary_lines.append("")
        for _, row in engine_summary.iterrows():
            summary_lines.append(
                f"- {row['engine']}: complexes={int(row['complexes'])}, best={row['best_affinity']:.3f}, mean_best={row['mean_best_affinity']:.3f}"
            )
        summary_lines.extend(
            [
                "",
                "Consensus ranking",
                "-----------------",
                f"- ranked_candidates={len(consensus_df)}",
                f"- rescoring_candidates={len(rescoring_df)}",
                f"- consensus_file={consensus_file}",
                f"- rescoring_file={rescoring_file}",
                f"- explainability_file={explainability_file}",
                f"- hit_classes_file={consensus_class_file}",
                f"- redocking_validation_file={redocking_validation_file}",
                f"- reference_baselines_file={reference_baselines_file}",
                f"- validation_summary_file={validation_summary_file}",
                f"- validation_gate_file={validation_gate_file}",
                f"- validation_gate_state={validation_gate.get('validation_gate_state', '')}",
                f"- validation_gate_allow_reference_anchor={validation_gate.get('allow_reference_anchor', False)}",
                f"- hit_class_policy_effective={effective_hit_class_policy}",
                f"- hit_class_policy_note={hit_class_policy_note}",
                f"- hit_classification_config_file={classification_reports.get('classification_config_file', '')}",
                f"- hit_classification_distribution_file={classification_reports.get('classification_distribution_file', '')}",
                f"- hit_classification_qc_breakdown_file={classification_reports.get('classification_qc_breakdown_file', '')}",
                f"- biology_mapping_file={biology_mapping_file}",
                f"- biology_correlation_global_file={biology_global_corr_file}",
                f"- biology_correlation_per_protein_file={biology_per_protein_corr_file}",
                f"- biology_unmatched_biology_keys={session_polypharm_outputs.get('unmatched_biology_file', '')}",
                f"- biology_unmatched_docking_keys={session_polypharm_outputs.get('unmatched_docking_file', '')}",
                f"- per_protein_rank_correlation_file={per_protein_corr_file}",
                f"- global_rank_correlation_file={global_corr_file}",
                f"- canonical_scores_root={score_root}",
                f"- polypharmacology_session_root={session_polypharm_outputs.get('polypharmacology_root', '')}",
                f"- polypharmacology_canonical_root={canonical_polypharm_outputs.get('polypharmacology_root', '')}",
                f"- top_pose_per_ligand_per_protein_file={session_top_pose_outputs.get('top_pose_per_ligand_per_protein_file', '')}",
                f"- top_pose_per_ligand_global_file={session_top_pose_outputs.get('top_pose_per_ligand_global_file', '')}",
                f"- top_pose_summary_file={session_top_pose_outputs.get('ligand_performance_summary_file', '')}",
                f"- top_pose_manifest_file={session_top_pose_outputs.get('top_pose_selection_manifest_file', '')}",
                f"- top_pose_canonical_root={top_pose_canonical_root}",
            ]
        )
        viz_summary = self._write_cross_engine_visualizations(best_by_engine, reports_dir) if run_visualizations else {}
        if viz_summary:
            summary_lines.extend(
                [
                    "",
                    "Cross-engine visualizations",
                    "--------------------------",
                    f"- paired_tags={viz_summary.get('paired_tags', 0)}",
                    f"- per_protein_batches={viz_summary.get('batch_count', 0)}",
                    f"- visualizations_dir={viz_summary.get('visualizations_dir', '')}",
                ]
            )

        compat_status = str(self.compat_shim_info.get("status", "") or "")
        if compat_status in {"legacy_dir_exists", "created", "existing_symlink"}:
            summary_lines.extend(
                [
                    "",
                    "Compatibility Notice",
                    "--------------------",
                    f"- compatibility_status={compat_status}",
                    f"- canonical_analysis_root={self.compat_shim_info.get('canonical_root', '')}",
                    f"- legacy_analysis_root={self.compat_shim_info.get('legacy_root', '')}",
                ]
            )
            if compat_status == "legacy_dir_exists":
                summary_lines.append("- existing legacy directory was preserved untouched; new results were written only to 5-Analysis/")
            elif compat_status == "created":
                summary_lines.append("- created compatibility symlink 5-Post_Docking_Analysis -> 5-Analysis for transition support")
            else:
                summary_lines.append("- existing compatibility symlink remains in place")

        class_dist_file = Path(str(classification_reports.get("classification_distribution_file", ""))).expanduser()
        class_qc_file = Path(str(classification_reports.get("classification_qc_breakdown_file", ""))).expanduser()
        class_dist_df = pd.read_csv(class_dist_file) if class_dist_file.exists() else pd.DataFrame()
        class_qc_df = pd.read_csv(class_qc_file) if class_qc_file.exists() else pd.DataFrame()
        scope_discrepancy_notice = self._scope_discrepancy_notice()

        run_context = {
            "run_id": self.run_id,
            "project_dir": str(self.project_dir),
            "output_dir": str(self.output_dir),
            "analysis_mode": self.analysis_mode,
            "analysis_scope": scope,
            "consensus_mode": self.consensus_mode,
            "normalization_method": self.normalization_method,
            "favorite_engine": str(self.favorite_engine or self.rerun_engine or ""),
            "engines_in_scope": list(self.engines_in_scope),
            "excluded_engines": list(self.excluded_engines),
            "scope_source": self.scope_source,
            "engine_preset_name": self.engine_preset_name,
            "scoped_engine_count": int(len(self.engines_in_scope)),
            "detected_engine_count": int(self.detected_engine_count or len(self.manifest.get("engines", []))),
            "scope_discrepancy_notice": scope_discrepancy_notice,
            "hit_class_policy_requested": str(self.hit_class_policy),
            "hit_class_policy_effective": str(effective_hit_class_policy),
            "validation_gate_state": str(validation_gate.get("validation_gate_state", "")),
            "validation_gate_allow_reference_anchor": bool(validation_gate.get("allow_reference_anchor", False)),
            "sqlite_dual_write_enabled": bool(self.sqlite_dual_write_enabled),
        }
        consolidated_outputs = generate_consolidated_run_summary(
            numbered_reports_root,
            run_context=run_context,
            engine_summary=engine_summary,
            classification_distribution=class_dist_df,
            qc_breakdown=class_qc_df,
            output_stem="consolidated_run_summary",
        )
        start_here_artifacts: Dict[str, object] = {
            "summary_txt": str(reports_dir / "summary.txt"),
            "combined_engine_scores_csv": str(reports_dir / "combined_engine_scores.csv"),
            "consensus_ranked_hits_csv": str(consensus_file),
            "consensus_ranked_hits_with_classes_csv": str(consensus_class_file),
            "rescoring_candidates_csv": str(rescoring_file),
            "consensus_explainability_json": str(explainability_file),
            "redocking_validation_csv": str(redocking_validation_file),
            "reference_baselines_csv": str(reference_baselines_file),
            "validation_summary_md": str(validation_summary_file),
            "validation_gate_status_json": str(validation_gate_file),
            "engine_rank_correlation_per_protein_csv": str(per_protein_corr_file),
            "engine_rank_correlation_global_csv": str(global_corr_file),
            "biology_mapping_report_json": str(biology_mapping_file),
            "biology_correlation_global_csv": str(biology_global_corr_file),
            "biology_correlation_per_protein_csv": str(biology_per_protein_corr_file),
            "top_pose_per_ligand_global_csv": session_top_pose_outputs.get("top_pose_per_ligand_global_file", ""),
            "top_pose_per_ligand_per_protein_csv": session_top_pose_outputs.get("top_pose_per_ligand_per_protein_file", ""),
            "sqlite_parity_csv": str(parity_file),
            "sqlite_write_summary": json.dumps(self.sqlite_write_summary),
        }
        start_here_artifacts.update(consolidated_outputs)
        dashboard_outputs = generate_dashboard_index(
            numbered_reports_root,
            run_context=run_context,
            artifact_paths=start_here_artifacts,
            output_stem="dashboard_export_index",
        )
        start_here_artifacts.update(dashboard_outputs)
        start_here_file = generate_start_here_index(
            numbered_reports_root,
            project_root=self.project_dir,
            run_context=run_context,
            artifact_paths=start_here_artifacts,
            output_name="START_HERE.md",
        )
        summary_lines.extend(
            [
                "",
                "Final Reports",
                "-------------",
                f"- numbered_reports_root={numbered_reports_root}",
                f"- start_here_file={start_here_file}",
                f"- consolidated_summary_json={consolidated_outputs.get('consolidated_summary_json_file', '')}",
                f"- consolidated_summary_csv={consolidated_outputs.get('consolidated_summary_csv_file', '')}",
                f"- dashboard_index_file={dashboard_outputs.get('dashboard_index_file', '')}",
                f"- dashboard_contract_file={dashboard_outputs.get('dashboard_contract_file', '')}",
                f"- sqlite_parity_file={parity_file}",
            ]
        )
        summary_file = reports_dir / "summary.txt"
        summary_file.write_text("\n".join(summary_lines), encoding="utf-8")
        start_here_artifacts["summary_txt"] = str(summary_file)
        generate_start_here_index(
            numbered_reports_root,
            project_root=self.project_dir,
            run_context=run_context,
            artifact_paths=start_here_artifacts,
            output_name="START_HERE.md",
        )

        self._materialize_canonical_analysis_view(
            reports_dir=reports_dir,
            raw_dir=raw_dir,
            visualizations_dir=Path(str(viz_summary.get("visualizations_dir", ""))).resolve()
            if viz_summary and viz_summary.get("visualizations_dir")
            else None,
        )
        generate_analysis_root_index(
            post_docking_root(self.project_dir),
            project_root=self.project_dir,
            session_root=self.output_dir,
            reports_root=numbered_reports_root,
        )

    @staticmethod
    def _copy_matching_files(source: Path, destination: Path, suffixes: Optional[Set[str]] = None) -> int:
        if source is None or not source.exists() or not source.is_dir():
            return 0
        destination.mkdir(parents=True, exist_ok=True)
        copied = 0
        normalized_suffixes = {value.lower() for value in (suffixes or set())}
        for file_path in sorted(source.rglob("*")):
            if not file_path.is_file():
                continue
            if normalized_suffixes and file_path.suffix.lower() not in normalized_suffixes:
                continue
            relative = file_path.relative_to(source)
            target = destination / relative
            try:
                if file_path.resolve() == target.resolve():
                    continue
            except Exception:
                if file_path == target:
                    continue
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(file_path, target)
            copied += 1
        return copied

    def _materialize_canonical_analysis_view(
        self,
        reports_dir: Path,
        raw_dir: Path,
        visualizations_dir: Optional[Path] = None,
    ) -> None:
        """
        Mirror root-level inventories and summary assets into the canonical
        `5-Analysis/` view without creating an extra nested analysis root.
        """
        canonical_root = post_docking_root(self.project_dir)
        raw_root = canonical_root / "raw_data"
        processed_root = canonical_root / "comparative"
        figures_2d_root = self.output_dir / "visualizations" / "2d"
        figures_3d_root = self.output_dir / "visualizations" / "3d"
        reports_root = canonical_root / "raw_data"

        raw_root.mkdir(parents=True, exist_ok=True)
        processed_root.mkdir(parents=True, exist_ok=True)
        figures_2d_root.mkdir(parents=True, exist_ok=True)
        figures_3d_root.mkdir(parents=True, exist_ok=True)
        reports_root.mkdir(parents=True, exist_ok=True)

        copied_raw = self._copy_matching_files(raw_dir, raw_root)
        copied_processed = self._copy_matching_files(
            reports_dir,
            processed_root,
            suffixes={".csv", ".json"},
        )
        copied_reports = self._copy_matching_files(
            reports_dir,
            reports_root,
            suffixes={".txt", ".md"},
        )

        copied_2d = 0
        copied_3d = 0
        if visualizations_dir and visualizations_dir.exists():
            for file_path in sorted(visualizations_dir.rglob("*")):
                if not file_path.is_file():
                    continue
                ext = file_path.suffix.lower()
                if ext not in {".png", ".svg", ".pdf", ".jpg", ".jpeg", ".gif", ".html", ".pse", ".pml"}:
                    continue
                relative = file_path.relative_to(visualizations_dir)
                # Route interactive/static 3D assets into figures/3d.
                is_3d = "3d" in {part.lower() for part in relative.parts} or ext in {".html", ".pse", ".pml"}
                if is_3d:
                    target = figures_3d_root / relative
                    copied_3d += 1
                else:
                    target = figures_2d_root / relative
                    copied_2d += 1
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(file_path, target)

        layout_manifest = {
            "analysis_root": str(canonical_root),
            "raw_dir": str(raw_root),
            "processed_dir": str(processed_root),
            "figures_2d_dir": str(self.output_dir / "visualizations" / "2d"),
            "figures_3d_dir": str(self.output_dir / "visualizations" / "3d"),
            "reports_dir": str(reports_root),
            "copied_counts": {
                "raw": int(copied_raw),
                "processed": int(copied_processed),
                "reports": int(copied_reports),
                "figures_2d": int(copied_2d),
                "figures_3d": int(copied_3d),
            },
        }
        (canonical_root / "layout_manifest.json").write_text(
            json.dumps(layout_manifest, indent=2),
            encoding="utf-8",
        )

    def _write_cross_engine_visualizations(
        self,
        best_by_engine: pd.DataFrame,
        reports_dir: Path,
    ) -> Dict[str, object]:
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError:
            return {}

        if best_by_engine.empty:
            return {}

        visualizations_dir = self.output_dir / "visualizations"
        visualizations_dir.mkdir(parents=True, exist_ok=True)
        batch_dir = visualizations_dir / "cross_engine_same_pair_batches"
        batch_dir.mkdir(parents=True, exist_ok=True)

        protein_name_map, ligand_name_map = self._build_comparative_name_maps()
        labelled = self._label_cross_engine_rows(best_by_engine, protein_name_map, ligand_name_map)
        labelled.to_csv(reports_dir / "cross_engine_pair_best_scores.csv", index=False)

        engine_order = [engine for engine in ("gnina", "vina", "smina") if engine in set(labelled["engine"].astype(str))]
        if not engine_order:
            engine_order = sorted(labelled["engine"].dropna().astype(str).unique().tolist())
        engine_palette = {
            "gnina": "#4c72b0",
            "vina": "#55a868",
            "smina": "#c17d11",
        }

        paired = labelled.groupby("tag", dropna=False).filter(
            lambda group: group["engine"].astype(str).nunique() >= 2
        ).copy()
        if paired.empty:
            self._save_cross_engine_placeholder(
                visualizations_dir / "cross_engine_pair_heatmap.png",
                "Cross-Engine Pair Heatmap",
                "No same-pair multi-engine rows were available for comparison.",
            )
            return {
                "paired_tags": 0,
                "batch_count": 0,
                "visualizations_dir": str(visualizations_dir),
            }

        sns.set_theme(style="whitegrid")

        wide = paired.pivot_table(
            index=["tag", "protein_label", "ligand_display_name"],
            columns="engine",
            values="affinity_kcal_mol",
            aggfunc="min",
        ).reset_index()
        present_engine_order = [engine for engine in engine_order if engine in wide.columns]
        if not present_engine_order:
            present_engine_order = [column for column in wide.columns if column not in {"tag", "protein_label", "ligand_display_name"}]
        wide["best_affinity"] = wide[present_engine_order].min(axis=1, skipna=True)
        wide["worst_affinity"] = wide[present_engine_order].max(axis=1, skipna=True)
        wide["affinity_spread"] = wide["worst_affinity"] - wide["best_affinity"]
        wide["is_reference_ligand"] = wide["ligand_display_name"].astype(str).str.contains(r"\[Ref\]", case=False, na=False)
        wide["pair_label"] = wide["protein_label"].astype(str) + " | " + wide["ligand_display_name"].astype(str)
        best_engine_indices = wide[present_engine_order].idxmin(axis=1)
        wide["best_engine"] = best_engine_indices.where(best_engine_indices.notna(), "")
        wide.sort_values(["is_reference_ligand", "best_affinity", "pair_label"], ascending=[False, True, True], inplace=True)

        heatmap_df = wide.head(min(60, len(wide))).set_index("pair_label")[present_engine_order]
        self._plot_cross_engine_pair_heatmap(
            heatmap_df,
            visualizations_dir / "cross_engine_pair_heatmap.png",
        )
        self._plot_cross_engine_distribution(
            paired,
            present_engine_order,
            engine_palette,
            visualizations_dir / "cross_engine_engine_distribution.png",
        )
        self._plot_cross_engine_disagreement(
            wide,
            engine_palette,
            visualizations_dir / "cross_engine_engine_disagreement.png",
        )

        batch_count = self._plot_cross_engine_per_protein_batches(
            wide,
            present_engine_order,
            engine_palette,
            batch_dir,
        )

        gnina_rows = labelled[labelled["engine"].astype(str) == "gnina"].copy()
        if not gnina_rows.empty and gnina_rows["cnn_affinity"].notna().any():
            gnina_rows["cnn_affinity"] = pd.to_numeric(gnina_rows["cnn_affinity"], errors="coerce")
            gnina_rows["cnn_score"] = pd.to_numeric(gnina_rows["cnn_score"], errors="coerce")
            gnina_rows["classical_strength"] = -pd.to_numeric(gnina_rows["affinity_kcal_mol"], errors="coerce")
            gnina_rows = gnina_rows[np.isfinite(gnina_rows["classical_strength"]) & np.isfinite(gnina_rows["cnn_affinity"])].copy()
            if not gnina_rows.empty:
                gnina_rows.sort_values(["classical_strength", "cnn_affinity"], ascending=[False, False], inplace=True)
                gnina_rows.to_csv(reports_dir / "gnina_cnn_pair_metrics.csv", index=False)
                self._plot_gnina_cnn_profiles(
                    gnina_rows,
                    visualizations_dir / "gnina_cnn_affinity_profile.png",
                )

        return {
            "paired_tags": int(paired["tag"].nunique()),
            "batch_count": int(batch_count),
            "visualizations_dir": str(visualizations_dir),
        }

    @staticmethod
    def _safe_subset(frame: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
        if frame is None or frame.empty:
            return pd.DataFrame(columns=columns)
        existing = [column for column in columns if column in frame.columns]
        return frame[existing].copy()

    @staticmethod
    def _write_biology_polypharmacology_bundle(
        root: Path,
        *,
        consensus_df: pd.DataFrame,
        biology_mapping_report: Dict[str, object],
        biology_global_corr: pd.DataFrame,
        biology_per_protein_corr: pd.DataFrame,
    ) -> Dict[str, str]:
        root.mkdir(parents=True, exist_ok=True)
        bio_columns = sorted([column for column in consensus_df.columns if str(column).startswith("bio_")])
        base_columns = [
            "protein",
            "ligand",
            "tag",
            "engine",
            "best_affinity_kcal_mol",
            "consensus_score",
            "consensus_rank",
            "docking_quality_class",
            "qc_status",
            "admet_status",
        ]
        pair_level = MultiEngineAnalysisPipeline._safe_subset(consensus_df, base_columns + bio_columns)
        pair_level_file = root / "polypharmacology_pair_level_with_biology.csv"
        pair_level.to_csv(pair_level_file, index=False)

        ligand_summary_cols = ["ligand", "best_affinity_kcal_mol", "consensus_score", "consensus_rank", "docking_quality_class"]
        ligand_summary = MultiEngineAnalysisPipeline._safe_subset(consensus_df, ligand_summary_cols + bio_columns)
        if not ligand_summary.empty and "ligand" in ligand_summary.columns:
            aggregations: Dict[str, str] = {}
            if "best_affinity_kcal_mol" in ligand_summary.columns:
                aggregations["best_affinity_kcal_mol"] = "min"
            if "consensus_score" in ligand_summary.columns:
                aggregations["consensus_score"] = "mean"
            if "consensus_rank" in ligand_summary.columns:
                aggregations["consensus_rank"] = "min"
            for column in bio_columns:
                if column in ligand_summary.columns:
                    aggregations[column] = "mean"
            if aggregations:
                ligand_summary = ligand_summary.groupby("ligand", dropna=False).agg(aggregations).reset_index()
            sort_columns = [column for column in ("consensus_rank", "best_affinity_kcal_mol") if column in ligand_summary.columns]
            if sort_columns:
                ligand_summary = ligand_summary.sort_values(sort_columns, ascending=[True] * len(sort_columns))
        ligand_summary_file = root / "polypharmacology_ligand_summary_with_biology.csv"
        ligand_summary.to_csv(ligand_summary_file, index=False)

        protein_summary_cols = ["protein", "best_affinity_kcal_mol", "consensus_score", "consensus_rank"] + bio_columns
        protein_summary = MultiEngineAnalysisPipeline._safe_subset(consensus_df, protein_summary_cols)
        if not protein_summary.empty and "protein" in protein_summary.columns:
            aggregations: Dict[str, str] = {}
            if "best_affinity_kcal_mol" in protein_summary.columns:
                aggregations["best_affinity_kcal_mol"] = "min"
            if "consensus_score" in protein_summary.columns:
                aggregations["consensus_score"] = "mean"
            if "consensus_rank" in protein_summary.columns:
                aggregations["consensus_rank"] = "mean"
            for column in bio_columns:
                if column in protein_summary.columns:
                    aggregations[column] = "mean"
            if aggregations:
                protein_summary = protein_summary.groupby("protein", dropna=False).agg(aggregations).reset_index()
            sort_columns = [column for column in ("consensus_score", "best_affinity_kcal_mol") if column in protein_summary.columns]
            if sort_columns:
                ascending = [False if column == "consensus_score" else True for column in sort_columns]
                protein_summary = protein_summary.sort_values(sort_columns, ascending=ascending)
        protein_summary_file = root / "polypharmacology_protein_summary_with_biology.csv"
        protein_summary.to_csv(protein_summary_file, index=False)

        biology_global_file = root / "biology_correlation_global.csv"
        biology_per_protein_file = root / "biology_correlation_per_protein.csv"
        biology_global_corr.to_csv(biology_global_file, index=False)
        biology_per_protein_corr.to_csv(biology_per_protein_file, index=False)

        mapping_file = root / "biology_mapping_report.json"
        mapping_file.write_text(json.dumps(biology_mapping_report, indent=2), encoding="utf-8")

        unresolved_biology = pd.DataFrame(biology_mapping_report.get("unmatched_biology_examples") or [])
        unresolved_docking = pd.DataFrame(biology_mapping_report.get("unmatched_docking_examples") or [])
        unresolved_biology_file = root / "biology_unmatched_biology_keys.csv"
        unresolved_docking_file = root / "biology_unmatched_docking_keys.csv"
        unresolved_biology.to_csv(unresolved_biology_file, index=False)
        unresolved_docking.to_csv(unresolved_docking_file, index=False)

        return {
            "polypharmacology_root": str(root),
            "pair_level_file": str(pair_level_file),
            "ligand_summary_file": str(ligand_summary_file),
            "protein_summary_file": str(protein_summary_file),
            "mapping_file": str(mapping_file),
            "global_correlation_file": str(biology_global_file),
            "per_protein_correlation_file": str(biology_per_protein_file),
            "unmatched_biology_file": str(unresolved_biology_file),
            "unmatched_docking_file": str(unresolved_docking_file),
        }

    def _build_comparative_name_maps(self) -> Tuple[Dict[str, Dict[str, str]], Dict[str, str]]:
        receptor_dir = shared_receptors_dir(self.project_dir)
        receptor_files = sorted(path for path in receptor_dir.glob("*") if path.is_file()) if receptor_dir.exists() else []

        protein_override_file = self.mapping_root / "protein_name_overrides.csv"
        ligand_override_file = self.mapping_root / "ligand_name_overrides.csv"

        protein_name_map = build_protein_name_mapping(
            receptor_files,
            pairlist_df=self.pairlist_df,
            overrides_file=protein_override_file if protein_override_file.exists() else None,
        )
        ligand_identifiers: List[str] = []
        if not self.pairlist_df.empty:
            for column in ("ligand", "ligand_name", "cocrystal_ligand_name"):
                if column in self.pairlist_df.columns:
                    ligand_identifiers.extend(self.pairlist_df[column].dropna().astype(str).tolist())
        ligand_name_map = build_ligand_name_mapping(
            ligand_identifiers,
            pairlist_df=self.pairlist_df,
            overrides_file=ligand_override_file if ligand_override_file.exists() else None,
        )

        if self.prompt_protein_names and sys.stdin.isatty():
            code_map = protein_name_map.get("pdb_code_to_name", {})
            receptor_map = protein_name_map.get("receptor_to_name", {})
            seen_codes = set()
            seen_receptors = set()
            prompt_targets: List[Tuple[str, str, str]] = []

            for receptor_file in receptor_files:
                receptor_name = str(receptor_file.name)
                receptor_stem = str(receptor_file.stem)
                pdb_code = extract_pdb_code(receptor_name) or extract_pdb_code(receptor_stem)
                if pdb_code:
                    if pdb_code in seen_codes:
                        continue
                    seen_codes.add(pdb_code)
                    current = str(code_map.get(pdb_code, "")).strip() or f"Protein {pdb_code}"
                    prompt_targets.append(("code", pdb_code, current))
                else:
                    if receptor_stem in seen_receptors:
                        continue
                    seen_receptors.add(receptor_stem)
                    current = str(
                        receptor_map.get(receptor_name)
                        or receptor_map.get(receptor_stem)
                        or format_protein_label(receptor_name, protein_name_map)
                    ).strip()
                    prompt_targets.append(("receptor", receptor_stem, current or receptor_stem))

            if not prompt_targets:
                for pdb_code in sorted(code_map):
                    current = str(code_map.get(pdb_code, "")).strip() or f"Protein {pdb_code}"
                    prompt_targets.append(("code", pdb_code, current))

            for target_type, target_id, current in prompt_targets:
                prompt = (
                    f"  Protein name for PDB {target_id} [{current}]: "
                    if target_type == "code"
                    else f"  Protein name for receptor {target_id} [{current}]: "
                )
                try:
                    updated = input(prompt).strip()
                except Exception:
                    updated = ""
                final_name = updated or current
                if target_type == "code":
                    code_map[target_id] = final_name
                    for receptor_key in list(receptor_map.keys()):
                        key_code = extract_pdb_code(receptor_key)
                        if key_code and key_code == target_id:
                            receptor_map[receptor_key] = final_name
                else:
                    receptor_map[target_id] = final_name
                    receptor_map[Path(target_id).stem] = final_name

        if self.prompt_ligand_names and sys.stdin.isatty():
            seen = set()
            prompt_targets: List[Tuple[str, str]] = []
            for source in sorted(ligand_name_map):
                stem = Path(str(source)).stem
                if source != stem:
                    continue
                display = str(ligand_name_map.get(source, "")).strip()
                if not display or stem in seen:
                    continue
                seen.add(stem)
                prompt_targets.append((stem, display))
            for ligand_id, current in prompt_targets:
                try:
                    updated = input(f"  Ligand name for {ligand_id} [{current}]: ").strip()
                except Exception:
                    updated = ""
                final_name = updated or current
                ligand_name_map[ligand_id] = final_name
                ligand_name_map[str(Path(ligand_id).stem)] = final_name

        protein_map_file = self.mapping_root / "protein_name_mapping.csv"
        ligand_map_file = self.mapping_root / "ligand_name_mapping.csv"
        mapping_to_dataframe(protein_name_map).to_csv(protein_map_file, index=False)
        ligand_mapping_to_dataframe(ligand_name_map).to_csv(ligand_map_file, index=False)
        if self.mapping_root != self.output_dir:
            mapping_to_dataframe(protein_name_map).to_csv(self.output_dir / "protein_name_mapping.csv", index=False)
            ligand_mapping_to_dataframe(ligand_name_map).to_csv(self.output_dir / "ligand_name_mapping.csv", index=False)

        if not protein_override_file.exists():
            rows = []
            for pdb_code, display_name in protein_name_map.get("pdb_code_to_name", {}).items():
                rows.append(
                    {
                        "pdb_code": pdb_code,
                        "receptor": "",
                        "display_name": display_name,
                        "notes": "Edit display_name if needed; keep pdb_code unchanged.",
                    }
                )
            pd.DataFrame(rows, columns=["pdb_code", "receptor", "display_name", "notes"]).to_csv(
                protein_override_file,
                index=False,
            )

        if not ligand_override_file.exists():
            rows = []
            seen_stems = set()
            for key, display_name in sorted(ligand_name_map.items()):
                stem = Path(str(key)).stem
                if stem in seen_stems:
                    continue
                seen_stems.add(stem)
                rows.append(
                    {
                        "ligand": stem,
                        "display_name": display_name,
                        "notes": "Edit display_name if needed; keep ligand unchanged.",
                    }
                )
            pd.DataFrame(rows, columns=["ligand", "display_name", "notes"]).to_csv(
                ligand_override_file,
                index=False,
            )

        return protein_name_map, ligand_name_map

    def _label_cross_engine_rows(
        self,
        rows: pd.DataFrame,
        protein_name_map: Dict[str, Dict[str, str]],
        ligand_name_map: Dict[str, str],
    ) -> pd.DataFrame:
        labelled = rows.copy()
        labelled["protein_display_name"] = labelled["protein"].astype(str).apply(
            lambda value: format_protein_label(value, protein_name_map)
        )
        labelled["protein_label"] = labelled["protein_display_name"]
        labelled["ligand_display_name"] = labelled["ligand"].astype(str).apply(
            lambda value: resolve_ligand_display_name(value, ligand_name_map)
        )
        labelled["is_reference_ligand"] = labelled["ligand"].astype(str).str.contains(r"^[0-9][A-Za-z0-9]{3}_ligand_", case=False, na=False)
        labelled.loc[labelled["is_reference_ligand"], "ligand_display_name"] = (
            labelled.loc[labelled["is_reference_ligand"], "ligand"].astype(str)
            + " [Ref]"
        )
        labelled["pair_label"] = labelled["protein_label"].astype(str) + " | " + labelled["ligand_display_name"].astype(str)
        return labelled

    @staticmethod
    def _save_cross_engine_placeholder(output_file: Path, title: str, message: str) -> None:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=11)
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.axis("off")
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)

    @staticmethod
    def _slugify_plot_name(value: str) -> str:
        slug = re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_")
        return slug[:80] or "plot"

    def _plot_cross_engine_pair_heatmap(self, matrix_df: pd.DataFrame, output_file: Path) -> None:
        import matplotlib.pyplot as plt
        import seaborn as sns

        if matrix_df.empty:
            self._save_cross_engine_placeholder(
                output_file,
                "Cross-Engine Pair Heatmap",
                "No same-pair multi-engine rows were available for heatmap rendering.",
            )
            return

        fig_height = max(7, len(matrix_df) * 0.30 + 2)
        annotate = len(matrix_df) <= 35
        fig, ax = plt.subplots(figsize=(11, fig_height))
        sns.heatmap(
            matrix_df,
            cmap="RdYlGn_r",
            annot=annotate,
            fmt=".2f",
            linewidths=0.35,
            linecolor="white",
            cbar_kws={"label": "Affinity (kcal/mol)"},
            ax=ax,
        )
        ax.set_title("Cross-Engine Affinity Heatmap (Same Protein-Ligand Pair)", fontsize=14, fontweight="bold")
        ax.set_xlabel("Docking Engine")
        ax.set_ylabel("Protein | Ligand")
        plt.setp(ax.get_yticklabels(), fontsize=8)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)

    def _plot_cross_engine_distribution(
        self,
        paired_rows: pd.DataFrame,
        engine_order: List[str],
        engine_palette: Dict[str, str],
        output_file: Path,
    ) -> None:
        import matplotlib.pyplot as plt
        import seaborn as sns

        fig, ax = plt.subplots(figsize=(10, 6))
        sns.boxplot(
            data=paired_rows,
            x="engine",
            y="affinity_kcal_mol",
            hue="engine",
            order=engine_order,
            hue_order=engine_order,
            palette=engine_palette,
            dodge=False,
            legend=False,
            ax=ax,
        )
        sns.stripplot(
            data=paired_rows,
            x="engine",
            y="affinity_kcal_mol",
            order=engine_order,
            color="black",
            alpha=0.35,
            size=3,
            jitter=0.18,
            ax=ax,
        )
        ax.set_title("Cross-Engine Affinity Distribution on Matched Pairs", fontsize=14, fontweight="bold")
        ax.set_xlabel("Docking Engine")
        ax.set_ylabel("Best-Pose Affinity (kcal/mol; positive = unfavorable)")
        ax.axhline(-7.0, color="red", linestyle="--", alpha=0.45)
        ax.axhline(0.0, color="darkorange", linestyle=":", alpha=0.65)
        stats_lines = []
        for engine in engine_order:
            sub = paired_rows[paired_rows["engine"].astype(str) == engine]
            if sub.empty:
                continue
            unfavorable = int((pd.to_numeric(sub["affinity_kcal_mol"], errors="coerce") > 0).sum())
            total = int(len(sub))
            stats_lines.append(f"{engine.upper()}: +ve {unfavorable}/{total}")
        if stats_lines:
            ax.text(
                0.98,
                0.98,
                "\n".join(stats_lines),
                ha="right",
                va="top",
                transform=ax.transAxes,
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7),
            )
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)

    def _plot_cross_engine_disagreement(
        self,
        wide: pd.DataFrame,
        engine_palette: Dict[str, str],
        output_file: Path,
    ) -> None:
        import matplotlib.pyplot as plt

        disagreement = wide.sort_values(["affinity_spread", "best_affinity"], ascending=[False, True]).head(min(30, len(wide))).copy()
        if disagreement.empty:
            self._save_cross_engine_placeholder(
                output_file,
                "Cross-Engine Disagreement",
                "No same-pair multi-engine rows were available for disagreement analysis.",
            )
            return

        colors = [engine_palette.get(str(value), "#999999") for value in disagreement["best_engine"].tolist()]
        fig_height = max(7, len(disagreement) * 0.34 + 1.5)
        fig, ax = plt.subplots(figsize=(12, fig_height))
        ax.barh(disagreement["pair_label"], disagreement["affinity_spread"], color=colors, edgecolor="black")
        ax.invert_yaxis()
        ax.set_xlabel("Worst-to-Best Affinity Spread (kcal/mol)")
        ax.set_ylabel("Protein | Ligand")
        ax.set_title("Largest Cross-Engine Disagreements (Same Pair)", fontsize=14, fontweight="bold")
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)

    def _plot_cross_engine_per_protein_batches(
        self,
        wide: pd.DataFrame,
        engine_order: List[str],
        engine_palette: Dict[str, str],
        batch_dir: Path,
    ) -> int:
        import matplotlib.pyplot as plt
        import numpy as np

        count = 0
        for protein_label, protein_df in wide.groupby("protein_label", dropna=False):
            matrix = protein_df.copy()
            matrix.sort_values(["best_affinity", "ligand_display_name"], inplace=True)
            matrix = matrix.head(min(28, len(matrix)))
            if matrix.empty:
                continue

            ligands = matrix["ligand_display_name"].astype(str).tolist()
            y_positions = np.arange(len(ligands), dtype=float)
            bar_width = 0.78 / max(len(engine_order), 1)
            fig_height = max(6, len(ligands) * 0.38 + 1.8)
            fig, ax = plt.subplots(figsize=(13, fig_height))

            for offset, engine in enumerate(engine_order):
                if engine not in matrix.columns:
                    continue
                values = pd.to_numeric(matrix[engine], errors="coerce").to_numpy(dtype=float)
                positions = y_positions - 0.39 + (offset + 0.5) * bar_width
                ax.barh(
                    positions,
                    values,
                    height=bar_width * 0.92,
                    color=engine_palette.get(engine, "#999999"),
                    edgecolor="black",
                    label=engine.upper(),
                )

            ax.set_yticks(y_positions)
            ax.set_yticklabels(ligands, fontsize=9)
            ax.axvline(-7.0, color="red", linestyle="--", alpha=0.45)
            ax.axvline(0.0, color="darkorange", linestyle=":", alpha=0.65)
            ax.set_xlabel("Best-Pose Affinity (kcal/mol)")
            ax.set_ylabel("Ligand")
            ax.set_title(f"Cross-Engine Comparison: {protein_label}", fontsize=14, fontweight="bold")
            ax.legend(loc="best")
            plt.tight_layout()
            output_file = batch_dir / f"{count + 1:02d}_{self._slugify_plot_name(str(protein_label))}_cross_engine.png"
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
            plt.close(fig)
            count += 1
        return count

    def _plot_gnina_cnn_profiles(self, gnina_rows: pd.DataFrame, output_file: Path) -> None:
        import matplotlib.pyplot as plt
        import numpy as np

        if gnina_rows.empty:
            self._save_cross_engine_placeholder(
                output_file,
                "GNINA CNN Affinity Profile",
                "No GNINA CNN affinity values were available for visualization.",
            )
            return

        ranked = gnina_rows.sort_values(["classical_strength", "cnn_affinity"], ascending=[False, False]).head(min(40, len(gnina_rows))).copy()
        ranked["plot_label"] = ranked["protein_label"].astype(str) + " | " + ranked["ligand_display_name"].astype(str)
        ranked = ranked.iloc[::-1].copy()

        fig, (ax_scatter, ax_rank) = plt.subplots(
            1,
            2,
            figsize=(18, max(7, len(ranked) * 0.22 + 3)),
            gridspec_kw={"width_ratios": [1.1, 1.2]},
        )

        size_values = pd.to_numeric(gnina_rows["cnn_score"], errors="coerce").fillna(0.0).clip(lower=0.0)
        scaled_sizes = 30 + 180 * (size_values / size_values.max()) if float(size_values.max()) > 0 else 60
        scatter = ax_scatter.scatter(
            gnina_rows["classical_strength"],
            gnina_rows["cnn_affinity"],
            c=pd.to_numeric(gnina_rows["cnn_score"], errors="coerce"),
            cmap="magma",
            s=scaled_sizes,
            alpha=0.78,
            edgecolors="white",
            linewidths=0.4,
        )
        ax_scatter.set_xlabel("Classical Affinity Strength (-vina_affinity)")
        ax_scatter.set_ylabel("GNINA CNN Affinity")
        ax_scatter.set_title("GNINA Classical vs CNN Affinity", fontsize=14, fontweight="bold")
        colorbar = fig.colorbar(scatter, ax=ax_scatter)
        colorbar.set_label("GNINA CNN Score")

        y_positions = np.arange(len(ranked), dtype=float)
        ax_rank.barh(
            y_positions,
            ranked["classical_strength"],
            color="#4c72b0",
            edgecolor="black",
            alpha=0.72,
            label="Classical strength (-vina_affinity)",
        )
        ax_rank.set_yticks(y_positions)
        ax_rank.set_yticklabels(ranked["plot_label"], fontsize=8)
        ax_rank.set_xlabel("Classical Affinity Strength")
        ax_rank.set_ylabel("Protein | Ligand")
        ax_rank.set_title("Top GNINA Pairs: Classical Strength vs CNN Affinity", fontsize=14, fontweight="bold")

        ax_rank_cnn = ax_rank.twiny()
        ax_rank_cnn.plot(
            ranked["cnn_affinity"],
            y_positions,
            color="#dd8452",
            marker="o",
            linewidth=2.0,
            markersize=4.5,
            label="CNN affinity",
        )
        ax_rank_cnn.set_xlabel("GNINA CNN Affinity")

        handles_left, labels_left = ax_rank.get_legend_handles_labels()
        handles_right, labels_right = ax_rank_cnn.get_legend_handles_labels()
        ax_rank.legend(handles_left + handles_right, labels_left + labels_right, loc="best", fontsize=8)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)

    def _write_exhaustive_rerun_manifest(self, scores: pd.DataFrame) -> str:
        self.manifest = load_manifest(self.project_dir)
        target_engine = self.rerun_engine or self.favorite_engine or self.manifest.get("favorite_engine") or ""
        if not target_engine:
            raise ValueError(
                "Comparative rerun promotion requires --rerun-engine or a favorite_engine in project_manifest.json"
            )

        best_by_engine = self._best_rows_by_group(scores, ["engine", "tag"], "affinity_kcal_mol")
        if self.min_affinity_advantage > 0 and best_by_engine["engine"].nunique() > 1:
            raise ValueError("Cross-engine raw affinity advantage is uncalibrated; use consensus rank promotion")
        engine_best = best_by_engine[best_by_engine["engine"] == target_engine].copy()
        if engine_best.empty:
            raise ValueError(f"No comparative scores were available for rerun engine '{target_engine}'")

        promoted = engine_best.merge(
            self._build_pair_competition(best_by_engine, target_engine),
            how="left",
            on="tag",
        )
        if self.winner_only:
            promoted = promoted[promoted["is_target_engine_winner"] == True].copy()
        if self.min_affinity_advantage > 0:
            promoted = promoted[
                promoted["affinity_advantage_kcal_mol"].fillna(float("-inf")) >= self.min_affinity_advantage
            ].copy()
        allowlist_tags = self._load_pair_allowlist_tags()
        if allowlist_tags:
            promoted = promoted[promoted["tag"].isin(allowlist_tags)].copy()

        # Consensus-first gating before engine-specific promotion details are finalized.
        consensus_df = build_consensus_rankings(
            best_by_engine=best_by_engine,
            consensus_mode=self.consensus_mode,
            favorite_engine=target_engine,
            normalization_method=self.normalization_method,
            expected_engines=self.engines_in_scope,
        )
        rescoring_df = select_rescoring_candidates(
            consensus_df=consensus_df,
            rescoring_scope=self.rescoring_scope,
            rescoring_top_n=self.rescoring_top_n,
        )
        consensus_tags = set(rescoring_df["tag"].astype(str).tolist()) if not rescoring_df.empty else set()
        if consensus_tags:
            promoted = promoted[promoted["tag"].astype(str).isin(consensus_tags)].copy()

        if promoted.empty:
            raise ValueError(f"No rerun candidates were selected for engine '{target_engine}'")

        promoted.sort_values(
            ["protein", "affinity_kcal_mol", "affinity_advantage_kcal_mol", "tag"],
            ascending=[True, True, False, True],
            inplace=True,
        )
        promoted = (
            promoted.groupby("protein", dropna=False)
            .head(self.top_per_protein)
            .copy()
        )
        if promoted.empty:
            raise ValueError(f"No rerun candidates were selected for engine '{target_engine}'")
        if self.max_rerun_pairs:
            promoted.sort_values(
                ["affinity_kcal_mol", "affinity_advantage_kcal_mol", "protein", "tag"],
                ascending=[True, False, True, True],
                inplace=True,
            )
            promoted = promoted.head(self.max_rerun_pairs).copy()

        pairlist = self.pairlist_df.copy()
        pairlist["tag"] = pairlist.apply(
            lambda row: f"{row['receptor']}_{row['site_id']}_{row['ligand']}",
            axis=1,
        )
        promoted["rank_within_protein"] = (
            promoted.groupby("protein", dropna=False)["affinity_kcal_mol"]
            .rank(method="first")
            .astype(int)
        )
        rerun_df = pairlist.merge(
            promoted[
                [
                    "tag",
                    "engine",
                    "protein",
                    "affinity_kcal_mol",
                    "rank_within_protein",
                    "winner_engine",
                    "winner_affinity_kcal_mol",
                    "best_other_engine",
                    "best_other_affinity_kcal_mol",
                    "affinity_advantage_kcal_mol",
                    "is_target_engine_winner",
                ]
            ],
            how="inner",
            on="tag",
        )
        if rerun_df.empty:
            raise ValueError(
                "Rerun candidate tags did not overlap the canonical pairlist during exhaustive manifest generation"
            )
        rerun_df["docking_mode"] = "exhaustive"
        rerun_df["round_id"] = str(self.manifest.get("latest_pair_round") or "round_001")
        rerun_df["promotion_rule"] = self._promotion_rule_label()
        rerun_df["consensus_mode"] = self.consensus_mode
        rerun_df["rescoring_scope"] = self.rescoring_scope
        rerun_df["rescoring_top_n"] = int(self.rescoring_top_n)
        rerun_df.rename(columns={"affinity_kcal_mol": "source_affinity_kcal_mol"}, inplace=True)
        rerun_df["promotion_reason"] = rerun_df.apply(self._promotion_reason, axis=1)
        rerun_df.sort_values(["protein", "rank_within_protein", "ligand"], inplace=True)

        reports_dir = self.output_dir / "reports"
        reports_dir.mkdir(parents=True, exist_ok=True)
        csv_path = reports_dir / "exhaustive_rerun_manifest.csv"
        json_path = reports_dir / "exhaustive_rerun_manifest.json"
        rerun_df.to_csv(csv_path, index=False)
        metadata = {
            "project_dir": str(self.project_dir),
            "analysis_mode": self.analysis_mode,
            "target_engine": target_engine,
            "top_per_protein": self.top_per_protein,
            "winner_only": self.winner_only,
            "min_affinity_advantage": self.min_affinity_advantage,
            "max_rerun_pairs": self.max_rerun_pairs,
            "pair_allowlist": self.pair_allowlist or "",
            "promotion_rule": self._promotion_rule_label(),
            "consensus_mode": self.consensus_mode,
            "rescoring_scope": self.rescoring_scope,
            "rescoring_top_n": int(self.rescoring_top_n),
            "consensus_candidate_count": int(len(consensus_df)),
            "rescoring_candidate_count": int(len(rescoring_df)),
            "candidate_count": int(len(rerun_df)),
            "round_id": str(rerun_df["round_id"].iloc[0]) if not rerun_df.empty else "",
            "csv_path": str(csv_path),
        }
        json_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        return str(csv_path)

    @staticmethod
    def _build_pair_competition(best_by_engine: pd.DataFrame, target_engine: str) -> pd.DataFrame:
        rows = []
        relative = normalize_engine_scores(best_by_engine, group_keys=["engine", "protein", "site_id"])
        for tag, group in relative.groupby("tag", dropna=False):
            ordered = group.sort_values(["normalized_affinity_score", "engine"], ascending=[False, True]).reset_index(drop=True)
            target_rows = ordered[ordered["engine"] == target_engine]
            if target_rows.empty:
                continue
            target_row = target_rows.iloc[0]
            best_row = ordered.iloc[0]
            other_rows = ordered[ordered["engine"] != target_engine]
            best_other = other_rows.iloc[0] if not other_rows.empty else None
            best_other_affinity = float(best_other["affinity_kcal_mol"]) if best_other is not None else None
            target_affinity = float(target_row["affinity_kcal_mol"])
            affinity_advantage = None
            # Raw affinity differences across engines are not meaningful.
            affinity_advantage = None
            rows.append(
                {
                    "tag": tag,
                    "winner_engine": str(best_row["engine"]),
                    "winner_affinity_kcal_mol": float(best_row["affinity_kcal_mol"]),
                    "best_other_engine": str(best_other["engine"]) if best_other is not None else "",
                    "best_other_affinity_kcal_mol": best_other_affinity,
                    "affinity_advantage_kcal_mol": affinity_advantage,
                    "is_target_engine_winner": bool(str(best_row["engine"]) == target_engine),
                    "winner_basis": "relative_within_engine_target_rank",
                }
            )
        return pd.DataFrame(rows)

    def _load_pair_allowlist_tags(self) -> Set[str]:
        if not self.pair_allowlist:
            return set()
        path = Path(self.pair_allowlist).expanduser().resolve()
        if not path.exists():
            raise ValueError(f"Pair allowlist file does not exist: {path}")
        if path.suffix.lower() in {".txt", ".list"}:
            return {line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}

        frame = pd.read_csv(path)
        if "tag" in frame.columns:
            return {str(value).strip() for value in frame["tag"].dropna().tolist() if str(value).strip()}
        required = {"receptor", "site_id", "ligand"}
        if required.issubset(frame.columns):
            tags: Set[str] = set()
            for _, row in frame.iterrows():
                receptor = str(row.get("receptor") or "").strip()
                site_id = str(row.get("site_id") or "").strip()
                ligand = str(row.get("ligand") or "").strip()
                if not receptor or not site_id or not ligand:
                    continue
                lowered = {receptor.lower(), site_id.lower(), ligand.lower()}
                if lowered.intersection({"nan", "none"}):
                    continue
                tags.add(f"{receptor}_{site_id}_{ligand}")
            return tags
        raise ValueError(
            f"Pair allowlist must contain either a 'tag' column or columns: {', '.join(sorted(required))}"
        )

    def _promotion_rule_label(self) -> str:
        labels = [
            f"top_per_protein={self.top_per_protein}",
            f"consensus_mode={self.consensus_mode}",
            f"rescoring_scope={self.rescoring_scope}",
            f"rescoring_top_n={self.rescoring_top_n}",
        ]
        if self.winner_only:
            labels.append("winner_only")
        if self.min_affinity_advantage > 0:
            labels.append(f"min_affinity_advantage={self.min_affinity_advantage:g}")
        if self.max_rerun_pairs > 0:
            labels.append(f"max_rerun_pairs={self.max_rerun_pairs}")
        if self.pair_allowlist:
            labels.append("pair_allowlist")
        return ";".join(labels)

    @staticmethod
    def _promotion_reason(row: pd.Series) -> str:
        reasons = [f"rank_within_protein={int(row['rank_within_protein'])}"]
        if bool(row.get("is_target_engine_winner", False)):
            reasons.append("target_engine_won_pair")
        if pd.notna(row.get("affinity_advantage_kcal_mol")):
            reasons.append(f"advantage={float(row['affinity_advantage_kcal_mol']):.3f}")
        return ";".join(reasons)

    def _write_single_engine_reports(self, scores: pd.DataFrame, engine: str, output_dir: Path) -> None:
        output_dir.mkdir(parents=True, exist_ok=True)
        engine_scores = scores[scores["engine"] == engine].copy()
        if engine_scores.empty:
            raise ValueError(f"No scores found for engine '{engine}'")

        ranking_column = "affinity_kcal_mol"
        ranking_label = "vina_affinity"
        if engine == "gnina" and engine_scores["cnn_affinity"].notna().any():
            engine_scores["cnn_affinity"] = pd.to_numeric(engine_scores["cnn_affinity"], errors="coerce")
            ranking_column = "cnn_affinity"
            ranking_label = "cnn_affinity"
            cnn_score = pd.to_numeric(engine_scores.get("cnn_score"), errors="coerce")
            engine_scores["cnn_confidence"] = pd.cut(
                cnn_score,
                bins=[-float("inf"), 0.3, 0.5, float("inf")],
                labels=["low", "moderate", "high"],
                right=False,
            ).astype(str)
            engine_scores["cnn_confidence"] = engine_scores["cnn_confidence"].replace("nan", "")
            engine_scores["heavy_atom_count"] = engine_scores["pose_file"].apply(
                lambda value: self._infer_heavy_atom_count(Path(str(value))) if str(value) else None
            )
            engine_scores["ligand_efficiency"] = engine_scores.apply(
                lambda row: (abs(float(row["cnn_affinity"])) / float(row["heavy_atom_count"]))
                if pd.notna(row.get("cnn_affinity")) and pd.notna(row.get("heavy_atom_count")) and float(row["heavy_atom_count"]) > 0
                else None,
                axis=1,
            )
        elif engine == "smina":
            engine_scores["rmsd_lb"] = pd.to_numeric(engine_scores.get("rmsd_lb"), errors="coerce")
            engine_scores["rmsd_ub"] = pd.to_numeric(engine_scores.get("rmsd_ub"), errors="coerce")
            invalid_rmsd = (engine_scores["rmsd_lb"] == -1.0) & (engine_scores["rmsd_ub"] == -1.0)
            engine_scores.loc[invalid_rmsd, ["rmsd_lb", "rmsd_ub"]] = None
            engine_scores["pose_diversity_status"] = "not_available"
            smina_meta = (
                (self.engine_detection_report.get("engines") or {}).get("smina", {})
                if isinstance(self.engine_detection_report, dict)
                else {}
            )
            engine_scores["smina_scoring_function"] = json.dumps(smina_meta.get("smina_scoring_weights", {}) or {})
            if bool(smina_meta.get("inconsistent_scoring_weights", False)):
                logger.warning("Smina scoring weights differ across logs for this project; provenance marked inconsistent.")
        elif engine == "vina":
            engine_scores["rmsd_lb"] = pd.to_numeric(engine_scores.get("rmsd_lb"), errors="coerce")
            engine_scores["pose_diversity_metric"] = engine_scores["rmsd_lb"]

        engine_scores["single_engine_mode"] = True
        engine_scores["primary_ranking_score"] = pd.to_numeric(engine_scores.get(ranking_column), errors="coerce")
        engine_scores["primary_ranking_label"] = ranking_label
        engine_scores = self._annotate_scope_columns(engine_scores)
        engine_scores.to_csv(output_dir / "scores.csv", index=False)
        engine_scores[UNIFIED_COMPAT_COLUMNS].to_csv(output_dir / "all_scores.csv", index=False)
        best, _selected_metric = self._select_best_pose_rows(engine_scores, requested_metric=ranking_label)
        best["pose"] = pd.to_numeric(best.get("pose"), errors="coerce")
        if _selected_metric == "cnn_affinity":
            best.sort_values(["cnn_affinity", "affinity_kcal_mol", "pose", "tag"], ascending=[False, True, True, True], inplace=True)
        else:
            best.sort_values(["affinity_kcal_mol", "pose", "tag"], inplace=True)
        best = self._annotate_scope_columns(best)
        best.to_csv(output_dir / "best_poses.csv", index=False)
        best[UNIFIED_COMPAT_COLUMNS].to_csv(output_dir / "best_poses_unified.csv", index=False)

        gnina_complex_export_count = 0
        if engine == "gnina":
            engine_layout = ensure_engine_layout(self.project_dir, "gnina")
            scores_csv = output_dir / "all_scores.csv"
            if scores_csv.exists():
                try:
                    gnina_complex_export_count = extract_best_poses_from_gnina(
                        input_dir=self.project_dir,
                        output_dir=output_dir / "complex_exports",
                        gnina_dir=engine_layout["poses"],
                        receptors_dir=shared_receptors_dir(self.project_dir),
                        scores_csv=scores_csv,
                        best_pose_criterion=ranking_label,
                    )
                except Exception as exc:
                    logger.warning("GNINA solo complex export failed: %s", exc)

        summary = (
            best.groupby("protein")
            .agg(best_affinity=("affinity_kcal_mol", "min"), ligand_count=("ligand", "nunique"))
            .reset_index()
            .sort_values(["best_affinity", "protein"])
        )
        summary = self._annotate_scope_columns(summary)
        summary.to_csv(output_dir / "protein_summary.csv", index=False)
        collapse_warning = ""
        # Vina reports distance from its best mode: mode 1 is zero by definition.
        if engine == "smina":
            smina_meta = (
                (self.engine_detection_report.get("engines") or {}).get("smina", {})
                if isinstance(self.engine_detection_report, dict)
                else {}
            )
            try:
                self.manifest.setdefault("engine_settings", {}).setdefault("smina", {})[
                    "smina_scoring_function"
                ] = smina_meta.get("smina_scoring_weights", {}) or {}
                save_manifest(self.project_dir, self.manifest)
            except Exception:
                pass
        (output_dir / "summary.txt").write_text(
            "\n".join(
                [
                    f"Engine: {engine}",
                    f"Primary ranking score: {ranking_label}",
                    f"Complexes: {best['tag'].nunique()}",
                    f"Best affinity: {best['affinity_kcal_mol'].min():.3f}",
                    f"Mean best affinity: {best['affinity_kcal_mol'].mean():.3f}",
                    *([f"Complex generation parser: sdf_multiconformer", f"Complex exports written: {gnina_complex_export_count}"] if engine == "gnina" else []),
                    *( [collapse_warning] if collapse_warning else [] ),
                ]
            ),
            encoding="utf-8",
        )
        numbered_layout = ensure_numbered_output_layout(self.project_dir)
        run_context = {
            "run_id": self.run_id,
            "project_dir": str(self.project_dir),
            "output_dir": str(output_dir),
            "analysis_mode": "single_engine",
            "analysis_scope": self.analysis_scope,
            "consensus_mode": self.consensus_mode,
            "normalization_method": self.normalization_method,
            "favorite_engine": str(engine),
            "engines_in_scope": list(self.engines_in_scope),
            "excluded_engines": list(self.excluded_engines),
            "scope_source": self.scope_source,
            "engine_preset_name": self.engine_preset_name,
            "scoped_engine_count": int(len(self.engines_in_scope)),
            "detected_engine_count": int(self.detected_engine_count or len(self.manifest.get("engines", []))),
            "scope_discrepancy_notice": self._scope_discrepancy_notice(),
        }
        start_here_artifacts = {
            "summary_txt": str(output_dir / "summary.txt"),
            "scores_csv": str(output_dir / "scores.csv"),
            "best_poses_csv": str(output_dir / "best_poses.csv"),
            "protein_summary_csv": str(output_dir / "protein_summary.csv"),
        }
        generate_start_here_index(
            numbered_layout["reports_root_numbered"],
            project_root=self.project_dir,
            run_context=run_context,
            artifact_paths=start_here_artifacts,
            output_name="START_HERE.md",
        )

    @staticmethod
    def _infer_element_from_pdbqt(atom_name: str, atom_type: str) -> str:
        """Infer a valid element symbol from a PDBQT atom name/type pair."""
        name = "".join(ch for ch in str(atom_name).strip() if ch.isalpha())
        adt = str(atom_type or "").strip()
        adt_upper = adt.upper()
        adt_lower = adt.lower()

        if adt_lower in {"cl", "br", "zn", "mg", "mn", "fe", "cu", "ni", "co", "na", "cd", "hg"}:
            return adt_lower.capitalize()
        if adt_lower == "ca":
            return "Ca" if adt == "Ca" else "C"
        if adt_upper in {"A", "C", "CA", "CG0", "CG1", "CG2"}:
            return "C"
        if adt_upper.startswith("N"):
            return "N"
        if adt_upper.startswith("O"):
            return "O"
        if adt_upper.startswith("S"):
            return "S"
        if adt_upper.startswith("P"):
            return "P"
        if adt_upper.startswith("H"):
            return "H"
        if adt_upper in {"F", "I"}:
            return adt_upper

        up = name.upper()
        if up.startswith("CL"):
            return "Cl"
        if up.startswith("BR"):
            return "Br"
        if up.startswith("ZN"):
            return "Zn"
        if up.startswith("MG"):
            return "Mg"
        if up.startswith("MN"):
            return "Mn"
        if up.startswith("FE"):
            return "Fe"
        if up.startswith("CU"):
            return "Cu"
        if up.startswith("NI"):
            return "Ni"
        if up.startswith("CO"):
            return "Co"
        if up.startswith("NA"):
            return "Na"
        if up.startswith("CD"):
            return "Cd"
        if up.startswith("HG"):
            return "Hg"
        if up:
            return up[0]
        return ""

    @classmethod
    def _pdbqt_atom_to_pdb_record(
        cls,
        line: str,
        record_name: str,
        chain_id: str,
        residue_name: Optional[str] = None,
        residue_number: Optional[str] = None,
        serial_override: Optional[int] = None,
    ) -> Optional[str]:
        """Convert one PDBQT atom record into a fixed-width PDB record."""
        if not line.startswith(("ATOM", "HETATM")):
            return None
        raw = line.rstrip("\n")
        if len(raw) < 80:
            raw = raw.ljust(80)
        raw = record_name + raw[6:]
        if serial_override is not None:
            raw = raw[:6] + f"{serial_override:5d}" + raw[11:]
        if residue_name is not None:
            raw = raw[:17] + f"{residue_name:>3}"[:3] + raw[20:]
        raw = raw[:21] + chain_id + raw[22:]
        if residue_number is not None:
            raw = raw[:22] + str(residue_number)[-4:].rjust(4) + raw[26:]
        atom_name = raw[12:16]
        atom_type = raw.split()[-1] if raw.split() else ""
        element = cls._infer_element_from_pdbqt(atom_name, atom_type)
        raw = raw[:76] + f"{element:>2}" + raw[78:]
        return raw[:80]

    @staticmethod
    def _normalize_pdb_record(
        line: str,
        record_name: str,
        chain_id: str,
        serial_override: Optional[int] = None,
    ) -> Optional[str]:
        """Normalize an ATOM/HETATM PDB line into a consistent fixed-width record."""
        if not line.startswith(("ATOM", "HETATM")):
            return None
        raw = line.rstrip("\n")
        if len(raw) < 80:
            raw = raw.ljust(80)
        raw = record_name + raw[6:]
        raw = raw[:21] + chain_id + raw[22:]
        if serial_override is not None:
            raw = raw[:6] + f"{serial_override:5d}" + raw[11:]
        return raw[:80]

    @staticmethod
    def _parse_atom_serial(pdb_line: str) -> Optional[int]:
        try:
            token = str(pdb_line)[6:11].strip()
            return int(token) if token else None
        except Exception:
            return None

    @staticmethod
    def _renumber_pdb_serials(lines: List[str]) -> List[str]:
        """Renumber all ATOM/HETATM serial numbers sequentially from 1.

        Eliminates duplicate serials that arise when receptor and ligand
        PDB lines are concatenated without coordinating their numbering.
        """
        result: List[str] = []
        serial = 0
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                serial += 1
                line = line[:6] + f"{serial:5d}" + line[11:]
            result.append(line)
        return result

    @staticmethod
    def _format_hetatm(
        serial: int,
        atom_name: str,
        res_name: str,
        chain_id: str,
        res_seq: str,
        x: float,
        y: float,
        z: float,
        element: str,
    ) -> str:
        """Format a single HETATM record in fixed-width PDB format."""
        try:
            resnum = int(str(res_seq)[-4:])
        except Exception:
            resnum = 1
        # PDB atom-name convention: 1-char elements are padded with a leading space.
        elem2 = element[:2].upper()
        if len(element) == 1:
            name_field = f" {atom_name:<3}"
        else:
            name_field = f"{atom_name:<4}"
        name_field = name_field[:4]
        line = (
            f"HETATM"
            f"{serial:5d}"
            f" "
            f"{name_field}"
            f" "  # alt loc
            f"{res_name:>3}"
            f" "  # space between resname and chain
            f"{chain_id}"
            f"{resnum:4d}"
            f"    "  # insertion code + 3 spaces (cols 27-30)
            f"{x:8.3f}"
            f"{y:8.3f}"
            f"{z:8.3f}"
            f"  1.00"
            f"  0.00"
            f"          "  # 10 spaces (cols 67-76)
            f"{elem2:>2}"
            f"  "  # charge
        )
        return line[:80]

    @staticmethod
    def _normalize_ligand_resname(raw_name: str) -> str:
        """Convert an arbitrary ligand token into a 3-character PDB residue name."""
        reserved = {
            "ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLN", "GLU", "GLX",
            "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
            "THR", "TRP", "TYR", "VAL", "SEC", "PYL",
        }
        token = "".join(ch for ch in str(raw_name).upper() if ch.isalnum())
        if not token:
            return "LIG"

        candidates: List[str] = []
        if len(token) >= 3:
            candidates.extend([token[:3], f"{token[0]}{token[-2:]}", f"{token[:2]}{token[-1]}", token[-3:]])
        elif len(token) == 2:
            candidates.extend([token + "X", f"{token[0]}X{token[1]}", f"X{token}"])
        else:
            candidates.extend([token + "XX", f"X{token}X", f"XX{token}"])

        normalized: List[str] = []
        seen = set()
        for candidate in candidates:
            code = "".join(ch for ch in candidate.upper() if ch.isalnum())[:3]
            if len(code) < 3:
                code = code.ljust(3, "X")
            if code and code not in seen:
                seen.add(code)
                normalized.append(code)

        for code in normalized:
            if code not in reserved:
                return code
        return normalized[0] if normalized else "LIG"

    def _derive_ligand_identity(self, ligand_name: str) -> Tuple[str, str, str]:
        """Derive a ligand residue identity from a prepared ligand filename."""
        source = str(ligand_name or "")
        ligand_match = re.search(r"_ligand_([A-Za-z0-9]{1,6})_([A-Za-z])_(\d+)", source)
        if ligand_match:
            return (
                self._normalize_ligand_resname(ligand_match.group(1)),
                ligand_match.group(2).upper(),
                ligand_match.group(3),
            )
        stem = Path(source).stem
        token = stem.split("_")[-1] if stem else "LIG"
        if token.lower() in {"pdbqt", "poses", "top", "out"} and "_" in stem:
            for part in reversed(stem.split("_")):
                if part.lower() not in {"pdbqt", "poses", "top", "out"}:
                    token = part
                    break
        return self._normalize_ligand_resname(token), "B", "1"

    def _read_receptor_lines(self, receptor_file: Path) -> List[str]:
        """Read receptor atoms from PDBQT or PDB into normalized PDB lines."""
        receptor_lines: List[str] = []
        if not receptor_file.exists():
            return receptor_lines
        with open(receptor_file, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if receptor_file.suffix.lower() == ".pdbqt":
                    converted = self._pdbqt_atom_to_pdb_record(line, line[:6], line[21:22] or " ")
                else:
                    converted = self._normalize_pdb_record(line, line[:6], line[21:22] or " ")
                if converted:
                    receptor_lines.append(converted)
        return receptor_lines

    def _extract_pose_lines_from_pdbqt(
        self,
        pose_file: Path,
        pose_number: int,
        ligand_name: str,
        starting_serial: int,
    ) -> List[str]:
        """Extract one ligand pose from a Vina/Smina-style PDBQT file as PDB HETATM records."""
        residue_name, chain_id, residue_number = self._derive_ligand_identity(ligand_name)
        ligand_lines: List[str] = []
        current_serial = starting_serial
        current_model = 0
        model_headers_seen = False
        capture_pose = pose_number == 1

        with open(pose_file, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith("MODEL"):
                    model_headers_seen = True
                    try:
                        current_model = int(line.split()[1])
                    except (IndexError, ValueError):
                        current_model += 1
                    capture_pose = current_model == pose_number
                    continue
                if line.startswith("ENDMDL"):
                    if capture_pose:
                        break
                    capture_pose = False
                    continue
                if model_headers_seen and not capture_pose:
                    continue
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                current_serial += 1
                converted = self._pdbqt_atom_to_pdb_record(
                    line,
                    "HETATM",
                    chain_id,
                    residue_name=residue_name,
                    residue_number=residue_number,
                    serial_override=current_serial,
                )
                if converted:
                    ligand_lines.append(converted)
        if model_headers_seen and pose_number != 1 and not ligand_lines:
            return []
        return ligand_lines

    def _extract_pose_lines_from_sdf(
        self,
        pose_file: Path,
        pose_number: int,
        ligand_name: str,
        starting_serial: int,
    ) -> List[str]:
        """Extract one ligand pose from a GNINA-style SDF file as PDB HETATM records.

        Supports V2000 and V3000 mol records. Poses are separated by `$$$$`.
        """
        residue_name, chain_id, residue_number = self._derive_ligand_identity(ligand_name)
        try:
            content = pose_file.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            return []

        # Each pose is one mol record terminated by $$$$
        records = content.split("$$$$")
        records = [r for r in records if r.strip()]
        if pose_number < 1 or pose_number > len(records):
            return []

        record_text = records[pose_number - 1]
        lines = record_text.splitlines()

        # Locate the counts line robustly.
        # SDF V2000 counts line: "aaabbblll..." — starts with two integers ≥0.
        # SDF V3000 counts line: contains "V3000" or body contains "M  V30".
        # The mol header is always: name / program / comment / counts (4 lines),
        # but GNINA SDFs have a blank name line and records after the first begin
        # with an extra blank separator line, so the offset varies.
        import re as _re
        _counts_re = _re.compile(r"^\s*(\d+)\s+(\d+)")
        counts_idx = None
        for _i, _ln in enumerate(lines):
            if "M  V30 COUNTS" in _ln or ("V3000" in _ln and _i < 10):
                counts_idx = _i
                break
            if _i > 0 and _counts_re.match(_ln):
                # Verify the preceding three lines look like a mol header
                # (name, program, comment) — all relatively short.
                counts_idx = _i
                break
        if counts_idx is None or counts_idx >= len(lines):
            return []

        current_serial = starting_serial
        ligand_lines_out: List[str] = []
        element_counts: Dict[str, int] = {}

        counts_line = lines[counts_idx]
        is_v3000 = "V3000" in counts_line or any("M  V30" in ln for ln in lines[:20])

        if is_v3000:
            # V3000: atom records look like: M  V30 <idx> <symbol> <x> <y> <z> [aamap]
            in_atom_block = False
            for line in lines:
                stripped = line.strip()
                if "M  V30 BEGIN ATOM" in stripped:
                    in_atom_block = True
                    continue
                if "M  V30 END ATOM" in stripped:
                    break
                if not in_atom_block or "M  V30" not in stripped:
                    continue
                parts = stripped.split()
                # Expected: M V30 index symbol x y z ...
                if len(parts) < 7:
                    continue
                try:
                    x = float(parts[4])
                    y = float(parts[5])
                    z = float(parts[6])
                    element = parts[3]
                except (ValueError, IndexError):
                    continue
                element_counts[element] = element_counts.get(element, 0) + 1
                atom_name = f"{element}{element_counts[element]}"
                current_serial += 1
                ligand_lines_out.append(
                    self._format_hetatm(current_serial, atom_name, residue_name,
                                        chain_id, residue_number, x, y, z, element)
                )
        else:
            # V2000: counts line cols 0-2 = num_atoms, atom block starts immediately after counts line
            try:
                num_atoms = int(counts_line[0:3].strip())
            except (ValueError, IndexError):
                return []
            atom_start = counts_idx + 1
            for i in range(num_atoms):
                atom_line = lines[atom_start + i] if (atom_start + i) < len(lines) else ""
                if len(atom_line) < 34:
                    continue
                try:
                    x = float(atom_line[0:10].strip())
                    y = float(atom_line[10:20].strip())
                    z = float(atom_line[20:30].strip())
                    element = atom_line[31:34].strip()
                except (ValueError, IndexError):
                    continue
                if not element:
                    continue
                element_counts[element] = element_counts.get(element, 0) + 1
                atom_name = f"{element}{element_counts[element]}"
                current_serial += 1
                ligand_lines_out.append(
                    self._format_hetatm(current_serial, atom_name, residue_name,
                                        chain_id, residue_number, x, y, z, element)
                )

        return ligand_lines_out

    @staticmethod
    def _display_complex_name(record: pd.Series) -> str:
        protein = Path(str(record.get("protein") or "")).stem or "protein"
        ligand = Path(str(record.get("ligand") or "")).stem or "ligand"
        site_id = str(record.get("site_id") or "site_1")
        return f"{protein}__{site_id}__{ligand}"

    @staticmethod
    def _best_rows_by_group(scores: pd.DataFrame, group_cols: List[str], value_col: str) -> pd.DataFrame:
        """Return one best row per group while skipping all-NaN affinity groups."""
        if scores.empty:
            return scores.copy()
        ranked = scores.copy()
        ranked[value_col] = pd.to_numeric(ranked[value_col], errors="coerce")
        valid = ranked[np.isfinite(ranked[value_col])]
        grouped = valid.groupby(group_cols, dropna=False)[value_col]
        idx = grouped.idxmin() if score_spec(value_col).lower_is_better else grouped.idxmax()
        idx = idx.dropna().astype(int)
        if idx.empty:
            return ranked.iloc[0:0].copy()
        best = ranked.loc[idx].copy()
        return best

    @classmethod
    def _best_by_tag(cls, scores: pd.DataFrame) -> pd.DataFrame:
        if scores.empty:
            return scores.copy()
        best = cls._best_rows_by_group(scores, ["tag"], "affinity_kcal_mol")
        best.sort_values(["affinity_kcal_mol", "tag"], inplace=True)
        return best

    @staticmethod
    def _resolve_best_pose_metric_for_frame(scores: pd.DataFrame, requested_metric: str) -> str:
        metric = normalize_best_pose_selection_metric(requested_metric)
        if metric != "auto":
            return metric
        if scores is None or scores.empty or "engine" not in scores.columns:
            return "vina_affinity"
        engine_series = scores["engine"].astype(str).str.strip().str.lower()
        engines = sorted({value for value in engine_series.tolist() if value})
        if len(engines) == 1 and engines[0] == "gnina" and "cnn_affinity" in scores.columns:
            cnn = pd.to_numeric(scores["cnn_affinity"], errors="coerce")
            if cnn.notna().any():
                return "cnn_affinity"
        return "vina_affinity"

    def _select_best_pose_rows(
        self,
        scores: pd.DataFrame,
        *,
        requested_metric: Optional[str] = None,
    ) -> Tuple[pd.DataFrame, str]:
        """
        Deterministically select one best row per tag with explicit tie-breaks.

        Tie-break priority:
        1) primary metric (vina_affinity or cnn_affinity; lower is better)
        2) vina_affinity (lower is better)
        3) pose index (lower is better)
        4) engine label (lexical) for deterministic stability
        """
        if scores is None or scores.empty:
            return pd.DataFrame(columns=scores.columns if isinstance(scores, pd.DataFrame) else []), "vina_affinity"

        metric = self._resolve_best_pose_metric_for_frame(
            scores,
            requested_metric if requested_metric is not None else self.best_pose_selection_metric,
        )

        working = scores.copy()
        working["tag"] = working["tag"].astype(str)
        working["pose"] = pd.to_numeric(working.get("pose"), errors="coerce").fillna(np.inf)
        working["affinity_kcal_mol"] = pd.to_numeric(working.get("affinity_kcal_mol"), errors="coerce")
        working["cnn_affinity"] = pd.to_numeric(working.get("cnn_affinity"), errors="coerce")
        if "engine" in working.columns:
            working["engine"] = working["engine"].astype(str)
        else:
            working["engine"] = ""

        if metric == "cnn_affinity":
            primary = pd.to_numeric(working.get("cnn_affinity"), errors="coerce")
        else:
            primary = pd.to_numeric(working.get("affinity_kcal_mol"), errors="coerce")
            metric = "vina_affinity"

        working["_primary_sort"] = primary.map(lambda value: sort_value(value, metric))
        working["_secondary_sort"] = pd.to_numeric(working.get("affinity_kcal_mol"), errors="coerce").fillna(np.inf)
        working["_pose_sort"] = pd.to_numeric(working.get("pose"), errors="coerce").fillna(np.inf)
        working["_engine_sort"] = working["engine"].astype(str)

        sorted_rows = working.sort_values(
            ["tag", "_primary_sort", "_secondary_sort", "_pose_sort", "_engine_sort"],
            ascending=[True, True, True, True, True],
        )
        best = sorted_rows.groupby("tag", dropna=False, sort=False).head(1).copy()
        best.drop(columns=["_primary_sort", "_secondary_sort", "_pose_sort", "_engine_sort"], inplace=True, errors="ignore")
        best.sort_values(["tag"], inplace=True)
        return best, metric

    def _build_downstream_results(self, engine_scores: pd.DataFrame) -> Dict[str, object]:
        """Build legacy-style result tables from the unified engine score schema."""
        full_data = engine_scores.copy()
        full_data["complex_name"] = full_data.apply(self._display_complex_name, axis=1)
        full_data["pose"] = pd.to_numeric(full_data["pose"], errors="coerce").fillna(0).astype(int)
        full_data["vina_affinity"] = pd.to_numeric(full_data["affinity_kcal_mol"], errors="coerce")
        full_data["cnn_affinity"] = pd.to_numeric(full_data["cnn_affinity"], errors="coerce")
        full_data["cnn_score"] = pd.to_numeric(full_data["cnn_score"], errors="coerce")

        best_poses = self._best_by_tag(full_data)
        summary_stats = (
            full_data.groupby(["complex_name", "protein", "ligand", "site_id"])
            .agg(
                vina_affinity_min=("vina_affinity", "min"),
                vina_affinity_max=("vina_affinity", "max"),
                vina_affinity_mean=("vina_affinity", "mean"),
                vina_affinity_std=("vina_affinity", "std"),
                pose_count=("pose", "count"),
            )
            .reset_index()
        )
        top_overall = best_poses[
            ["complex_name", "protein", "ligand", "site_id", "vina_affinity", "pose", "engine"]
        ].head(10)
        best_per_protein = self._best_rows_by_group(best_poses, ["protein"], "vina_affinity").sort_values(
            ["vina_affinity", "protein"]
        ).reset_index(drop=True)
        best_per_ligand = self._best_rows_by_group(best_poses, ["ligand"], "vina_affinity").sort_values(
            ["vina_affinity", "ligand"]
        ).reset_index(drop=True)
        protein_summary = (
            best_poses.groupby("protein")
            .agg(
                best_affinity=("vina_affinity", "min"),
                mean_affinity=("vina_affinity", "mean"),
                ligand_count=("ligand", "nunique"),
                complex_count=("complex_name", "nunique"),
            )
            .reset_index()
            .sort_values(["best_affinity", "protein"])
        )
        ligand_summary = (
            best_poses.groupby("ligand")
            .agg(
                best_affinity=("vina_affinity", "min"),
                mean_affinity=("vina_affinity", "mean"),
                protein_count=("protein", "nunique"),
                complex_count=("complex_name", "nunique"),
            )
            .reset_index()
            .sort_values(["best_affinity", "ligand"])
        )
        return {
            "full_data": full_data,
            "best_poses": best_poses,
            "summary_stats": summary_stats,
            "top_overall": top_overall,
            "best_per_protein": best_per_protein,
            "best_per_ligand": best_per_ligand,
            "protein_summary": protein_summary,
            "ligand_summary": ligand_summary,
        }

    def _extract_best_pose_complexes(self, best_poses: pd.DataFrame, poses_dir: Path) -> int:
        """Create best-pose receptor-ligand complex PDBs for Vina/Smina favorite-engine continuation."""
        poses_dir.mkdir(parents=True, exist_ok=True)
        extracted = 0
        pose_manifest_rows: List[Dict[str, object]] = []

        for _, row in best_poses.iterrows():
            pose_file = Path(str(row.get("pose_file") or ""))
            receptor_name = str(row.get("protein") or "")
            receptor_file = shared_receptors_dir(self.project_dir) / receptor_name if receptor_name else Path("")
            if not pose_file.exists() or not receptor_file.exists():
                pose_manifest_rows.append(
                    {
                        "complex_name": row.get("complex_name"),
                        "tag": row.get("tag"),
                        "pose": int(pd.to_numeric(row.get("pose"), errors="coerce") or 0),
                        "selection_criterion": str(row.get("selection_criterion") or self.best_pose_selection_metric),
                        "selected_score": pd.to_numeric(row.get("selected_score"), errors="coerce"),
                        "vina_affinity": pd.to_numeric(row.get("affinity_kcal_mol"), errors="coerce"),
                        "cnn_affinity": pd.to_numeric(row.get("cnn_affinity"), errors="coerce"),
                        "pose_file": str(pose_file),
                        "receptor_file": str(receptor_file),
                        "status": "missing_input",
                    }
                )
                continue

            receptor_lines = self._read_receptor_lines(receptor_file)
            if not receptor_lines:
                pose_manifest_rows.append(
                    {
                        "complex_name": row.get("complex_name"),
                        "tag": row.get("tag"),
                        "pose": int(pd.to_numeric(row.get("pose"), errors="coerce") or 0),
                        "selection_criterion": str(row.get("selection_criterion") or self.best_pose_selection_metric),
                        "selected_score": pd.to_numeric(row.get("selected_score"), errors="coerce"),
                        "vina_affinity": pd.to_numeric(row.get("affinity_kcal_mol"), errors="coerce"),
                        "cnn_affinity": pd.to_numeric(row.get("cnn_affinity"), errors="coerce"),
                        "pose_file": str(pose_file),
                        "receptor_file": str(receptor_file),
                        "status": "empty_receptor",
                    }
                )
                continue

            max_serial = 0
            for receptor_line in receptor_lines:
                serial = self._parse_atom_serial(receptor_line)
                if serial is not None:
                    max_serial = max(max_serial, serial)

            pose_number = int(row.get("pose", 1))
            ligand_str = str(row.get("ligand") or "")
            ext = pose_file.suffix.lower()
            if ext == ".sdf":
                ligand_lines = self._extract_pose_lines_from_sdf(
                    pose_file, pose_number, ligand_str, max_serial,
                )
            else:
                ligand_lines = self._extract_pose_lines_from_pdbqt(
                    pose_file, pose_number, ligand_str, max_serial,
                )
            if not ligand_lines:
                pose_manifest_rows.append(
                    {
                        "complex_name": row.get("complex_name"),
                        "tag": row.get("tag"),
                        "pose": int(pd.to_numeric(row.get("pose"), errors="coerce") or 0),
                        "selection_criterion": str(row.get("selection_criterion") or self.best_pose_selection_metric),
                        "selected_score": pd.to_numeric(row.get("selected_score"), errors="coerce"),
                        "vina_affinity": pd.to_numeric(row.get("affinity_kcal_mol"), errors="coerce"),
                        "cnn_affinity": pd.to_numeric(row.get("cnn_affinity"), errors="coerce"),
                        "pose_file": str(pose_file),
                        "receptor_file": str(receptor_file),
                        "status": "missing_pose_atoms",
                    }
                )
                continue

            output_name = str(row.get("bridge_output_name") or row.get("tag") or row.get("complex_name") or "").strip()
            if not output_name:
                output_name = f"complex_{extracted + 1}"
            output_file = poses_dir / f"{output_name}.pdb"
            for sidecar_suffix in (".ligand.sdf", ".geometry.json"):
                output_file.with_suffix(sidecar_suffix).unlink(missing_ok=True)
            all_lines = self._renumber_pdb_serials(receptor_lines + ligand_lines)
            output_file.write_text("\n".join(all_lines + ["END"]) + "\n", encoding="utf-8")
            output_file.with_suffix(".receptor.pdb").write_text("\n".join(receptor_lines + ["END"]) + "\n", encoding="utf-8")
            try:
                from rdkit import Chem
                ligand_molecule = load_pose_molecule(pose_file, pose_number)
                with Chem.SDWriter(str(output_file.with_suffix(".ligand.sdf"))) as writer:
                    writer.write(ligand_molecule)
                output_file.with_suffix(".geometry.json").write_text(json.dumps({"receptor_frame_id": row.get("receptor_frame_id"), "pose_file": str(pose_file), "pose": pose_number}), encoding="utf-8")
            except Exception as exc:
                logger.warning("Complex %s is display-only: authoritative chemical graph unavailable (%s)", output_name, exc)
            validation = validate_complex_pdb_structure(output_file, receptor_reference=receptor_file)
            validation_errors = [str(item) for item in (validation.get("errors") or [])]
            validation_warnings = [str(item) for item in (validation.get("warnings") or [])]
            extraction_status = "extracted" if validation.get("is_valid", False) else "invalid_complex_output"
            pose_manifest_rows.append(
                {
                    "complex_name": row.get("complex_name"),
                    "tag": row.get("tag"),
                    "pose": int(row.get("pose", 1)),
                    "selection_criterion": str(row.get("selection_criterion") or self.best_pose_selection_metric),
                    "selected_score": pd.to_numeric(row.get("selected_score"), errors="coerce"),
                    "vina_affinity": pd.to_numeric(row.get("affinity_kcal_mol"), errors="coerce"),
                    "cnn_affinity": pd.to_numeric(row.get("cnn_affinity"), errors="coerce"),
                    "bridge_output_name": output_name,
                    "pose_file": str(pose_file),
                    "receptor_file": str(receptor_file),
                    "output_pdb": str(output_file),
                    "status": extraction_status,
                    "validation_is_valid": bool(validation.get("is_valid", False)),
                    "validation_errors": "; ".join(validation_errors),
                    "validation_warnings": "; ".join(validation_warnings),
                    "receptor_atom_count": int(validation.get("receptor_atom_count", 0) or 0),
                    "ligand_atom_count": int(validation.get("ligand_atom_count", 0) or 0),
                }
            )
            if extraction_status == "extracted":
                extracted += 1

        if pose_manifest_rows:
            pd.DataFrame(pose_manifest_rows).to_csv(poses_dir / "pose_extraction_manifest.csv", index=False)
        return extracted

    def _build_simplified_bridge_scores(self, engine_scores: pd.DataFrame) -> pd.DataFrame:
        """Build a simplified-pipeline-compatible score table from unified engine scores."""
        bridge_scores = engine_scores.copy()
        bridge_scores["mode"] = pd.to_numeric(bridge_scores["pose"], errors="coerce").fillna(0).astype(int)
        bridge_scores["pose"] = bridge_scores["mode"]
        bridge_scores["vina_affinity"] = pd.to_numeric(bridge_scores["affinity_kcal_mol"], errors="coerce")
        bridge_scores["cnn_affinity"] = pd.to_numeric(bridge_scores["cnn_affinity"], errors="coerce")
        bridge_scores["cnn_score"] = pd.to_numeric(bridge_scores["cnn_score"], errors="coerce")
        bridge_scores["tag"] = bridge_scores["tag"].astype(str)
        return bridge_scores

    def _build_simplified_bridge_complexes(self, best_poses: pd.DataFrame) -> List[Dict[str, object]]:
        """Build simplified-pipeline complex manifest entries for best-pose PDB complexes."""
        complexes: List[Dict[str, object]] = []
        for _, row in best_poses.iterrows():
            tag = str(row.get("tag") or "").strip()
            if not tag:
                continue
            receptor_name = str(row.get("protein") or "").strip()
            ligand_name = str(row.get("ligand") or "").strip()
            complexes.append(
                {
                    "complex_name": tag,
                    "receptor_name": receptor_name,
                    "receptor_file": str(shared_receptors_dir(self.project_dir) / receptor_name) if receptor_name else "",
                    "ligand_name": ligand_name,
                    "pose_file": str(row.get("pose_file") or ""),
                    "site_id": str(row.get("site_id") or ""),
                    "pose": int(pd.to_numeric(row.get("pose"), errors="coerce") or 0),
                    "protein_display_name": Path(receptor_name).stem if receptor_name else "",
                    "protein_label": Path(receptor_name).stem if receptor_name else "",
                }
            )
        return complexes

    def _write_pose_complex_pdb(
        self,
        pose_file: Path,
        receptor_name: str,
        ligand_name: str,
        pose_number: int,
        output_file: Path,
    ) -> str:
        """Write one receptor-ligand complex PDB extracted from a multi-pose PDBQT file."""
        receptor_file = shared_receptors_dir(self.project_dir) / receptor_name if receptor_name else Path("")
        if not pose_file.exists() or not receptor_file.exists():
            return "missing_input"

        receptor_lines = self._read_receptor_lines(receptor_file)
        if not receptor_lines:
            return "empty_receptor"

        max_serial = 0
        for receptor_line in receptor_lines:
            serial = self._parse_atom_serial(receptor_line)
            if serial is not None:
                max_serial = max(max_serial, serial)

        ligand_lines = self._extract_pose_lines_from_pdbqt(
            pose_file,
            pose_number,
            ligand_name,
            max_serial,
        )
        if not ligand_lines:
            return "missing_pose_atoms"

        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text("\n".join(receptor_lines + ligand_lines + ["END"]) + "\n", encoding="utf-8")
        validation = validate_complex_pdb_structure(output_file, receptor_reference=receptor_file)
        if not validation.get("is_valid", False):
            return "invalid_complex_output"
        return "extracted"

    @staticmethod
    def _mirror_complexes_to_best_poses(complexes_dir: Path, best_poses_dir: Path) -> int:
        """Copy extracted best-pose complexes into the legacy best_poses_pdb location."""
        best_poses_dir.mkdir(parents=True, exist_ok=True)
        copied = 0
        for complex_file in sorted(complexes_dir.glob("*.pdb")):
            if complex_file.name.endswith(".receptor.pdb"):
                continue
            shutil.copy2(complex_file, best_poses_dir / complex_file.name)
            for suffix in (".ligand.sdf", ".receptor.pdb", ".geometry.json"):
                sidecar = complex_file.with_suffix(suffix)
                if sidecar.exists():
                    shutil.copy2(sidecar, best_poses_dir / sidecar.name)
            copied += 1
        return copied

    @staticmethod
    def _write_complex_export_index(
        *,
        output_dir: Path,
        complexes_dir: Path,
        best_poses_dir: Optional[Path] = None,
        favorite_engine: str = "",
        source_table: Optional[pd.DataFrame] = None,
    ) -> Path:
        """
        Emit an index for exported complex PDB assets used by rescoring/visualization.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        source_lookup: Dict[str, Dict[str, object]] = {}
        if source_table is not None and not source_table.empty:
            for _, row in source_table.iterrows():
                tag = str(row.get("tag") or "").strip()
                if tag:
                    source_lookup[tag] = {
                        "protein": str(row.get("protein") or ""),
                        "ligand": str(row.get("ligand") or ""),
                        "site_id": str(row.get("site_id") or ""),
                        "pose": int(pd.to_numeric(row.get("pose"), errors="coerce") or 0),
                        "affinity_kcal_mol": float(pd.to_numeric(row.get("affinity_kcal_mol"), errors="coerce"))
                        if pd.notna(pd.to_numeric(row.get("affinity_kcal_mol"), errors="coerce"))
                        else None,
                    }

        rows: List[Dict[str, object]] = []
        for pdb_file in sorted(complexes_dir.glob("*.pdb")):
            if pdb_file.name.endswith(".receptor.pdb"):
                continue
            stem = pdb_file.stem
            source = source_lookup.get(stem, {})
            validation = validate_complex_pdb_structure(pdb_file)
            validation_errors = [str(item) for item in (validation.get("errors") or [])]
            validation_warnings = [str(item) for item in (validation.get("warnings") or [])]
            rows.append(
                {
                    "tag": stem,
                    "favorite_engine": favorite_engine,
                    "complex_pdb": str(pdb_file),
                    "best_pose_pdb": str(best_poses_dir / pdb_file.name)
                    if best_poses_dir is not None and (best_poses_dir / pdb_file.name).exists()
                    else "",
                    "protein": str(source.get("protein", "")),
                    "ligand": str(source.get("ligand", "")),
                    "site_id": str(source.get("site_id", "")),
                    "pose": int(source.get("pose") or 0),
                    "affinity_kcal_mol": source.get("affinity_kcal_mol"),
                    "exists": True,
                    "is_valid_complex": bool(validation.get("is_valid", False)),
                    "validation_errors": "; ".join(validation_errors),
                    "validation_warnings": "; ".join(validation_warnings),
                    "receptor_atom_count": int(validation.get("receptor_atom_count", 0) or 0),
                    "ligand_atom_count": int(validation.get("ligand_atom_count", 0) or 0),
                }
            )
        index_file = output_dir / "complex_export_index.csv"
        pd.DataFrame(rows).to_csv(index_file, index=False)
        return index_file

    @staticmethod
    def _pairwise_values(matrix: np.ndarray) -> np.ndarray:
        """Return finite upper-triangle RMSD values."""
        if matrix.size == 0 or matrix.ndim != 2:
            return np.array([], dtype=float)
        upper = matrix[np.triu_indices(matrix.shape[0], k=1)]
        return upper[np.isfinite(upper)]

    @staticmethod
    def _summarize_rmsd_values(values: np.ndarray) -> Dict[str, Optional[float]]:
        """Summarize finite RMSD values."""
        if values.size == 0:
            return {
                "mean_pairwise_rmsd": None,
                "median_pairwise_rmsd": None,
                "std_pairwise_rmsd": None,
                "min_pairwise_rmsd": None,
                "max_pairwise_rmsd": None,
            }
        return {
            "mean_pairwise_rmsd": float(np.mean(values)),
            "median_pairwise_rmsd": float(np.median(values)),
            "std_pairwise_rmsd": float(np.std(values)),
            "min_pairwise_rmsd": float(np.min(values)),
            "max_pairwise_rmsd": float(np.max(values)),
        }

    @staticmethod
    def _write_json(output_file: Path, payload: Dict[str, object]) -> None:
        with open(output_file, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)

    @staticmethod
    def _save_rmsd_heatmap(matrix_df: pd.DataFrame, output_file: Path, title: str) -> None:
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except Exception:
            return

        fig, ax = plt.subplots(
            figsize=(max(7, len(matrix_df.columns) * 0.7), max(6, len(matrix_df.index) * 0.7))
        )
        sns.heatmap(
            matrix_df,
            cmap="viridis",
            square=True,
            cbar_kws={"label": "RMSD (Å)"},
            ax=ax,
        )
        ax.set_title(title, fontsize=11)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)

    def _run_rmsd_scope(
        self,
        scope_dir: Path,
        pose_rows: pd.DataFrame,
        scope_name: str,
        heatmap_title: str,
        label_column: str = "tag",
        inputs_filename: str = "best_pose_inputs.csv",
        matrix_filename: str = "best_pose_rmsd_matrix.csv",
        summary_basename: str = "best_pose_rmsd_summary",
        heatmap_filename: str = "rmsd_heatmap.png",
        enable_enhanced_visualizations: bool = True,
        visualizations_dir_name: Optional[str] = "visualizations",
    ) -> Dict[str, object]:
        """Run a best-pose RMSD scope using extracted complex PDB files."""
        scope_dir.mkdir(parents=True, exist_ok=True)
        export_df = pose_rows.copy()
        if "pdb_file" in export_df.columns:
            export_df["pdb_file"] = export_df["pdb_file"].astype(str)
        inputs_file = scope_dir / inputs_filename
        export_df.to_csv(inputs_file, index=False)
        visualizations_dir = scope_dir / visualizations_dir_name if visualizations_dir_name else None
        summary_csv = scope_dir / f"{summary_basename}.csv"
        summary_json = scope_dir / f"{summary_basename}.json"
        state_file = scope_dir / f"{summary_basename}_scope_state.json"
        if any(column in export_df and export_df[column].nunique() > 1 for column in ("ligand", "protein_name", "protein_label")):
            summary = {"scope": scope_name, "state": "not_evaluable", "reason": "cross_ligand_or_cross_receptor_RMSD_requires_explicit_mapping", "pose_count": len(export_df)}
            summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
            pd.DataFrame([summary]).to_csv(summary_csv, index=False)
            return summary

        def _scope_signature(frame: pd.DataFrame) -> str:
            signature_df = frame.copy()
            for column in signature_df.columns:
                signature_df[column] = signature_df[column].map(
                    lambda value: "" if pd.isna(value) else str(value)
                )
            payload = json.dumps(
                {
                    "scope_name": scope_name,
                    "force_global_rmsd": self.force_global_rmsd,
                    "global_rmsd_defer_threshold": self.global_rmsd_defer_threshold,
                    "label_column": label_column,
                    "inputs_filename": inputs_filename,
                    "matrix_filename": matrix_filename,
                    "summary_basename": summary_basename,
                    "signature_csv": signature_df.to_csv(index=False),
                    "source_content": {str(path): content_hash(path) for value in signature_df.to_numpy().ravel() for original in [Path(str(value))] if original.is_file() for path in [original, original.with_suffix(".ligand.sdf"), original.with_suffix(".geometry.json")] if path.is_file()},
                },
                sort_keys=True,
            )
            return hashlib.sha256(payload.encode("utf-8")).hexdigest()

        def _persist_scope_state(input_signature: str, state: str, summary: Dict[str, object]) -> None:
            self._write_json(
                state_file,
                {
                    "input_signature": input_signature,
                    "state": state,
                    "summary_file": str(summary_json),
                    "scope": str(summary.get("scope") or ""),
                },
            )

        input_signature = _scope_signature(export_df)
        if self.resume_rmsd and state_file.exists() and summary_json.exists():
            try:
                state_payload = json.loads(state_file.read_text(encoding="utf-8"))
                cached_summary = json.loads(summary_json.read_text(encoding="utf-8"))
                if (
                    isinstance(cached_summary, dict)
                    and str(state_payload.get("input_signature") or "") == input_signature
                ):
                    cached_summary["resumed"] = True
                    cached_summary["state"] = str(
                        cached_summary.get("state") or state_payload.get("state") or "completed"
                    )
                    return cached_summary
            except Exception:
                pass

        if export_df.empty:
            summary = {
                "scope": scope_name,
                "state": "skipped_empty",
                "pose_count": 0,
                "pairwise_count": 0,
                "matrix_file": "",
                "inputs_file": str(inputs_file),
                "visualizations_dir": str(visualizations_dir) if visualizations_dir else "",
                "input_signature": input_signature,
                "mean_pairwise_rmsd": None,
                "median_pairwise_rmsd": None,
                "std_pairwise_rmsd": None,
                "min_pairwise_rmsd": None,
                "max_pairwise_rmsd": None,
            }
            pd.DataFrame([summary]).to_csv(summary_csv, index=False)
            self._write_json(summary_json, summary)
            _persist_scope_state(input_signature, "skipped_empty", summary)
            return summary

        if (
            scope_name == "global_best_poses"
            and not self.force_global_rmsd
            and self.global_rmsd_defer_threshold > 0
            and len(export_df) > self.global_rmsd_defer_threshold
        ):
            summary = {
                "scope": scope_name,
                "state": "deferred",
                "pose_count": int(len(export_df)),
                "pairwise_count": 0,
                "matrix_file": "",
                "inputs_file": str(inputs_file),
                "visualizations_dir": str(visualizations_dir) if visualizations_dir else "",
                "input_signature": input_signature,
                "mean_pairwise_rmsd": None,
                "median_pairwise_rmsd": None,
                "std_pairwise_rmsd": None,
                "min_pairwise_rmsd": None,
                "max_pairwise_rmsd": None,
                "reason": (
                    f"Deferred global best-pose RMSD because pose_count={len(export_df)} "
                    f"exceeds threshold={self.global_rmsd_defer_threshold}"
                ),
            }
            pd.DataFrame([summary]).to_csv(summary_csv, index=False)
            self._write_json(summary_json, summary)
            _persist_scope_state(input_signature, "deferred", summary)
            return summary

        try:
            from post_docking_analysis.enhanced_rmsd_analyzer import (
                analyze_conformational_diversity_enhanced,
                analyze_pose_clustering_enhanced,
                calculate_rmsd_matrix_from_pdbs,
                create_rmsd_visualizations_enhanced,
            )
            enhanced_available = True
        except Exception:
            enhanced_available = False

        pdb_files = [Path(str(value)) for value in export_df["pdb_file"].tolist()]
        matrix = np.zeros((len(pdb_files), len(pdb_files)), dtype=float)
        if len(pdb_files) > 1:
            if enhanced_available:
                requested_workers = self.rmsd_workers
                try:
                    matrix, _ = calculate_rmsd_matrix_from_pdbs(
                        pdb_files,
                        ligand_only=True,
                        num_workers=requested_workers,
                    )
                except Exception as exc:
                    should_retry_single = requested_workers != 1
                    if should_retry_single:
                        logger.warning(
                            "RMSD parallel run failed (%s). Retrying with single worker.",
                            exc,
                        )
                        try:
                            matrix, _ = calculate_rmsd_matrix_from_pdbs(
                                pdb_files,
                                ligand_only=True,
                                num_workers=1,
                            )
                        except Exception as retry_exc:
                            logger.warning(
                                "RMSD single-worker retry also failed (%s). "
                                "Continuing with empty RMSD matrix for this scope.",
                                retry_exc,
                            )
                            matrix = np.full((len(pdb_files), len(pdb_files)), np.nan, dtype=float)
                            np.fill_diagonal(matrix, 0.0)
                    else:
                        logger.warning(
                            "RMSD calculation failed with single worker (%s). "
                            "Continuing with empty RMSD matrix for this scope.",
                            exc,
                        )
                        matrix = np.full((len(pdb_files), len(pdb_files)), np.nan, dtype=float)
                        np.fill_diagonal(matrix, 0.0)
            else:
                matrix = np.full((len(pdb_files), len(pdb_files)), np.nan, dtype=float)
                np.fill_diagonal(matrix, 0.0)

        labels = export_df[label_column].astype(str).tolist() if label_column in export_df.columns else export_df["tag"].astype(str).tolist()
        matrix_df = pd.DataFrame(matrix, index=labels, columns=labels)
        matrix_file = scope_dir / matrix_filename
        matrix_df.to_csv(matrix_file)

        pair_values = self._pairwise_values(matrix)
        summary = {
            "scope": scope_name,
            "state": "completed",
            "pose_count": int(len(export_df)),
            "pairwise_count": int(pair_values.size),
            "matrix_file": str(matrix_file),
            "inputs_file": str(inputs_file),
            "visualizations_dir": str(visualizations_dir) if visualizations_dir else "",
            "input_signature": input_signature,
            **self._summarize_rmsd_values(pair_values),
        }
        pd.DataFrame([summary]).to_csv(summary_csv, index=False)
        self._write_json(summary_json, summary)
        _persist_scope_state(input_signature, "completed", summary)

        if enhanced_available and enable_enhanced_visualizations and visualizations_dir is not None:
            try:
                cluster_count = max(1, min(3, len(export_df)))
                clustering_results = analyze_pose_clustering_enhanced(
                    export_df,
                    matrix,
                    pdb_files,
                    n_clusters=cluster_count,
                )
                diversity_results = analyze_conformational_diversity_enhanced(
                    export_df,
                    matrix,
                )
                create_rmsd_visualizations_enhanced(
                    clustering_results,
                    diversity_results,
                    visualizations_dir,
                    dpi=300,
                )

                cluster_summary = clustering_results.get("cluster_summary", pd.DataFrame())
                if isinstance(cluster_summary, pd.DataFrame) and not cluster_summary.empty:
                    cluster_summary.to_csv(scope_dir / "cluster_summary.csv", index=False)

                cluster_centroids = clustering_results.get("cluster_centroids", pd.DataFrame())
                if isinstance(cluster_centroids, pd.DataFrame) and not cluster_centroids.empty:
                    cluster_centroids.to_csv(scope_dir / "cluster_centroids.csv", index=False)

                poses_with_clusters = clustering_results.get("poses_with_clusters", pd.DataFrame())
                if isinstance(poses_with_clusters, pd.DataFrame) and not poses_with_clusters.empty:
                    poses_with_clusters.to_csv(scope_dir / "poses_with_clusters.csv", index=False)

                diversity_metrics = diversity_results.get("diversity_metrics", pd.DataFrame())
                if isinstance(diversity_metrics, pd.DataFrame) and not diversity_metrics.empty:
                    diversity_metrics.to_csv(scope_dir / "diversity_metrics.csv", index=False)
            except Exception:
                self._save_rmsd_heatmap(matrix_df, scope_dir / heatmap_filename, heatmap_title)
        else:
            self._save_rmsd_heatmap(matrix_df, scope_dir / heatmap_filename, heatmap_title)

        return summary

    def _run_non_gnina_rmsd_bridge(
        self,
        engine_scores: pd.DataFrame,
        best_poses: pd.DataFrame,
        complexes_dir: Path,
        bridge_out: Path,
    ) -> Dict[str, object]:
        """Run a real RMSD stage for PDBQT-based engines using extracted complex PDBs."""
        rmsd_dir = bridge_out / "rmsd_analysis"
        rmsd_dir.mkdir(parents=True, exist_ok=True)
        per_complex_dir = rmsd_dir / "per_complex_all_poses"
        per_protein_dir = rmsd_dir / "per_protein_best_poses"
        global_best_dir = rmsd_dir / "global_best_poses"
        per_complex_dir.mkdir(exist_ok=True)
        per_protein_dir.mkdir(exist_ok=True)
        global_best_dir.mkdir(exist_ok=True)
        enabled_scopes = set(self.rmsd_scopes)
        enable_per_complex = "per_complex" in enabled_scopes
        enable_per_protein = "per_protein" in enabled_scopes
        enable_global = "global" in enabled_scopes

        per_complex_rows: List[Dict[str, object]] = []
        processed_complexes = 0

        if enable_per_complex:
            for tag, tag_df in engine_scores.groupby("tag", dropna=False):
                tag = str(tag or "").strip()
                if not tag:
                    continue
                tag_scope_dir = per_complex_dir / re.sub(r"[^A-Za-z0-9]+", "_", tag).strip("_")
                pose_structures_dir = tag_scope_dir / "pose_structures"
                pose_structures_dir.mkdir(parents=True, exist_ok=True)

                tag_df = tag_df.copy()
                tag_df["pose"] = pd.to_numeric(tag_df["pose"], errors="coerce").fillna(0).astype(int)
                tag_df["vina_affinity"] = pd.to_numeric(tag_df["vina_affinity"], errors="coerce")
                tag_df.sort_values(["pose"], inplace=True)

                pose_rows: List[Dict[str, object]] = []
                extraction_rows: List[Dict[str, object]] = []
                for _, row in tag_df.iterrows():
                    pose_number = int(row.get("pose", 0))
                    if pose_number <= 0:
                        continue
                    pose_label = f"mode_{pose_number}"
                    output_file = pose_structures_dir / f"{pose_label}.pdb"
                    status = self._write_pose_complex_pdb(
                        Path(str(row.get("pose_file") or "")),
                        str(row.get("protein") or ""),
                        str(row.get("ligand") or ""),
                        pose_number,
                        output_file,
                    )
                    extraction_rows.append(
                        {
                            "tag": tag,
                            "pose": pose_number,
                            "pose_label": pose_label,
                            "status": status,
                            "output_pdb": str(output_file),
                        }
                    )
                    if status != "extracted":
                        continue
                    pose_rows.append(
                        {
                            "tag": tag,
                            "pose": pose_number,
                            "pose_label": pose_label,
                            "protein": str(row.get("protein") or ""),
                            "protein_label": Path(str(row.get("protein") or "")).stem,
                            "ligand": str(row.get("ligand") or ""),
                            "vina_affinity": float(row["vina_affinity"]) if pd.notna(row.get("vina_affinity")) else None,
                            "pdb_file": output_file,
                        }
                    )

                if extraction_rows:
                    pd.DataFrame(extraction_rows).to_csv(tag_scope_dir / "pose_extraction_manifest.csv", index=False)
                if len(pose_rows) < 2:
                    continue

                pose_df = pd.DataFrame(pose_rows)
                scope_summary = self._run_rmsd_scope(
                    tag_scope_dir,
                    pose_df,
                    scope_name=f"per_complex_all_poses::{tag}",
                    heatmap_title=f"All-Pose RMSD: {Path(str(tag_df['protein'].iloc[0])).stem} | {str(tag_df['ligand'].iloc[0])}",
                    label_column="pose_label",
                    inputs_filename="all_pose_inputs.csv",
                    matrix_filename="all_pose_rmsd_matrix.csv",
                    summary_basename="all_pose_rmsd_summary",
                    heatmap_filename="all_pose_rmsd_heatmap.png",
                    enable_enhanced_visualizations=False,
                    visualizations_dir_name=None,
                )
                best_row = pose_df.loc[pose_df["vina_affinity"].idxmin()] if pose_df["vina_affinity"].notna().any() else pose_df.iloc[0]
                per_complex_rows.append(
                    {
                        "tag": tag,
                        "protein_name": Path(str(tag_df["protein"].iloc[0])).stem,
                        "protein_label": Path(str(tag_df["protein"].iloc[0])).stem,
                        "ligand": str(tag_df["ligand"].iloc[0]),
                        "pose_count": int(len(pose_df)),
                        "best_mode": int(best_row["pose"]),
                        "best_affinity": float(best_row["vina_affinity"]) if pd.notna(best_row["vina_affinity"]) else None,
                        "matrix_file": scope_summary.get("matrix_file", ""),
                        "inputs_file": scope_summary.get("inputs_file", ""),
                        "visualizations_dir": scope_summary.get("visualizations_dir", ""),
                        "pairwise_count": scope_summary.get("pairwise_count", 0),
                        "mean_pairwise_rmsd": scope_summary.get("mean_pairwise_rmsd"),
                        "median_pairwise_rmsd": scope_summary.get("median_pairwise_rmsd"),
                        "std_pairwise_rmsd": scope_summary.get("std_pairwise_rmsd"),
                        "min_pairwise_rmsd": scope_summary.get("min_pairwise_rmsd"),
                        "max_pairwise_rmsd": scope_summary.get("max_pairwise_rmsd"),
                    }
                )
                processed_complexes += 1

        per_complex_summary_file = per_complex_dir / "per_complex_rmsd_summary.csv"
        pd.DataFrame(per_complex_rows).to_csv(per_complex_summary_file, index=False)

        best_pose_rows: List[Dict[str, object]] = []
        best_pose_lookup = best_poses.copy()
        best_pose_lookup["pose"] = pd.to_numeric(best_pose_lookup["pose"], errors="coerce").fillna(0).astype(int)
        best_pose_lookup["vina_affinity"] = pd.to_numeric(best_pose_lookup["vina_affinity"], errors="coerce")
        for _, row in best_pose_lookup.iterrows():
            tag = str(row.get("tag") or "").strip()
            if not tag:
                continue
            pdb_file = complexes_dir / f"{tag}.pdb"
            if not pdb_file.exists():
                continue
            best_pose_rows.append(
                {
                    "tag": tag,
                    "protein_name": Path(str(row.get("protein") or "")).stem,
                    "protein_label": Path(str(row.get("protein") or "")).stem,
                    "ligand": str(row.get("ligand") or ""),
                    "vina_affinity": float(row["vina_affinity"]) if pd.notna(row.get("vina_affinity")) else None,
                    "pdb_file": pdb_file,
                }
            )

        best_pose_df = pd.DataFrame(best_pose_rows)
        per_protein_rows: List[Dict[str, object]] = []
        if enable_per_protein and not best_pose_df.empty:
            for protein_label, protein_df in best_pose_df.groupby("protein_label", dropna=False):
                protein_scope_dir = per_protein_dir / re.sub(r"[^A-Za-z0-9]+", "_", str(protein_label)).strip("_")
                protein_summary = self._run_rmsd_scope(
                    protein_scope_dir,
                    protein_df.sort_values(["vina_affinity", "ligand"], na_position="last").reset_index(drop=True),
                    scope_name=f"per_protein_best_poses::{protein_label}",
                    heatmap_title=f"Best-Pose RMSD: {protein_label}",
                )
                protein_summary.update(
                    {
                        "protein_name": str(protein_label),
                        "protein_label": str(protein_label),
                        "ligand_count": int(protein_df["ligand"].nunique()),
                    }
                )
                per_protein_rows.append(protein_summary)

        per_protein_summary_file = per_protein_dir / "per_protein_best_pose_rmsd_summary.csv"
        pd.DataFrame(per_protein_rows).to_csv(per_protein_summary_file, index=False)

        if enable_global:
            global_summary = self._run_rmsd_scope(
                global_best_dir,
                best_pose_df.sort_values(["vina_affinity", "protein_label", "ligand"], na_position="last").reset_index(drop=True),
                scope_name="global_best_poses",
                heatmap_title="Global Best-Pose RMSD Across Proteins",
            )
        else:
            global_summary = {
                "scope": "global_best_poses",
                "state": "skipped",
                "pose_count": 0,
                "pairwise_count": 0,
                "matrix_file": "",
                "inputs_file": "",
                "visualizations_dir": "",
                "input_signature": "",
                "mean_pairwise_rmsd": None,
                "median_pairwise_rmsd": None,
                "std_pairwise_rmsd": None,
                "min_pairwise_rmsd": None,
                "max_pairwise_rmsd": None,
                "reason": "Global RMSD scope disabled by configuration.",
            }
            pd.DataFrame([global_summary]).to_csv(global_best_dir / "best_pose_rmsd_summary.csv", index=False)
            self._write_json(global_best_dir / "best_pose_rmsd_summary.json", global_summary)
            self._write_json(
                global_best_dir / "best_pose_rmsd_summary_scope_state.json",
                {
                    "input_signature": "",
                    "state": "skipped",
                    "summary_file": str(global_best_dir / "best_pose_rmsd_summary.json"),
                    "scope": "global_best_poses",
                },
            )

        return {
            "state": "completed" if (processed_complexes > 0 or not best_pose_df.empty) else "skipped_or_failed",
            "processed_complexes": processed_complexes,
            "per_complex_summary": str(per_complex_summary_file),
            "per_protein_summary": str(per_protein_summary_file),
            "global_summary_csv": str(global_best_dir / "best_pose_rmsd_summary.csv"),
            "global_summary": global_summary,
        }

    def _run_non_gnina_downstream_bridge(
        self,
        favorite_engine: str,
        engine_scores: pd.DataFrame,
        output_dir: Path,
        notes_file: Path,
    ) -> None:
        bridge_out = output_dir / "deep_analysis"
        bridge_out.mkdir(parents=True, exist_ok=True)
        raw_dir = bridge_out / "raw_data"
        raw_dir.mkdir(parents=True, exist_ok=True)
        complexes_dir = bridge_out / "complexes"
        best_poses_pdb_dir = bridge_out / "best_poses_pdb"

        downstream_results = self._build_downstream_results(engine_scores)
        full_data = downstream_results["full_data"]
        best_poses = downstream_results["best_poses"]
        full_data[UNIFIED_COMPAT_COLUMNS].to_csv(raw_dir / "unified_all_scores.csv", index=False)
        best_poses[UNIFIED_COMPAT_COLUMNS].to_csv(raw_dir / "unified_best_poses.csv", index=False)

        simplified_scores = self._build_simplified_bridge_scores(engine_scores)
        simplified_scores.to_csv(bridge_out / "all_scores.csv", index=False)

        from post_docking_analysis.simplified_input_handler import find_receptor_files
        from post_docking_analysis.simplified_pipeline import SimplifiedPostDockingPipeline

        simplified_pipeline = SimplifiedPostDockingPipeline(
            sdf_folder=str(ensure_engine_layout(self.project_dir, favorite_engine)["poses"]),
            log_folder=str(ensure_engine_layout(self.project_dir, favorite_engine)["logs"]),
            receptors_folder=str(shared_receptors_dir(self.project_dir)),
            output_dir=str(bridge_out),
            pairlist_file=str(pairlist_path(self.project_dir)),
            run_rmsd=False,
            prompt_protein_names=self.prompt_protein_names,
            prompt_ligand_names=self.prompt_ligand_names,
            enable_poseview=self.enable_poseview,
            rmsd_scopes=self.rmsd_scopes,
            resume_rmsd=self.resume_rmsd,
            force_global_rmsd=self.force_global_rmsd,
            global_rmsd_defer_threshold=self.global_rmsd_defer_threshold,
            shared_mapping_dir=str(self.mapping_root),
            exclude_problematic_ligands=self.exclude_problematic_ligands,
            positive_affinity_threshold=self.positive_affinity_threshold,
            minimum_pose_count=self.minimum_pose_count,
            rmsd_workers=self.rmsd_workers,
        )
        simplified_pipeline.pairlist_df = self.pairlist_df.copy()
        simplified_pipeline.receptor_files = find_receptor_files(shared_receptors_dir(self.project_dir))
        simplified_pipeline.protein_name_map = build_protein_name_mapping(
            simplified_pipeline.receptor_files,
            simplified_pipeline.pairlist_df,
            overrides_file=simplified_pipeline.protein_name_override_file,
        )
        simplified_pipeline._prompt_for_protein_name_overrides()
        simplified_pipeline._persist_protein_name_mapping()
        simplified_pipeline.ligand_name_map = simplified_pipeline._build_ligand_name_map()
        simplified_pipeline._prompt_for_ligand_name_overrides()
        simplified_pipeline._persist_ligand_name_mapping()
        simplified_pipeline.scores_df = simplified_scores.copy()
        simplified_pipeline.complexes = self._build_simplified_bridge_complexes(best_poses)

        manifest_df = best_poses.copy()
        manifest_df["bridge_output_name"] = manifest_df["tag"].astype(str)
        extracted_count = self._extract_best_pose_complexes(manifest_df, complexes_dir)
        mirrored_best_pose_count = self._mirror_complexes_to_best_poses(complexes_dir, best_poses_pdb_dir)
        complex_index_file = self._write_complex_export_index(
            output_dir=bridge_out,
            complexes_dir=complexes_dir,
            best_poses_dir=best_poses_pdb_dir,
            favorite_engine=favorite_engine,
            source_table=manifest_df,
        )

        affinity_ok = simplified_pipeline._analyze_binding_affinity() if extracted_count else False
        polypharm_ok = simplified_pipeline._analyze_polypharmacology() if extracted_count else False
        extract_ok = simplified_pipeline._extract_poses() if extracted_count else False
        reports_ok = simplified_pipeline._generate_reports() if extracted_count else False
        viz_ok = False
        prolif_ok = False
        ligplot_ok = False
        poseview_ok = False
        pandamap_ok = False
        py3dmol_ok = False
        organize_ok = False
        fast_mode_note = ""
        if extracted_count:
            if self.fast_mode:
                fast_mode_note = (
                    "Fast mode active: skipped optional heavy visualization/interaction stages "
                    "(visualizations, ProLIF, LigPlot, PoseView, PandaMap, py3Dmol)."
                )
                logger.info("⚡ %s", fast_mode_note)
            else:
                viz_ok = simplified_pipeline._generate_visualizations()
                saved_complexes = simplified_pipeline.complexes
                try:
                    simplified_pipeline.complexes = []
                    prolif_ok = simplified_pipeline._generate_prolif_interaction_maps()
                finally:
                    simplified_pipeline.complexes = saved_complexes
                ligplot_ok = simplified_pipeline._generate_ligplot_diagrams()
                poseview_ok = simplified_pipeline._generate_poseview_diagrams()
                pandamap_ok = simplified_pipeline._generate_pandamap_analysis()
                py3dmol_ok = simplified_pipeline._generate_py3dmol_visualizations()
            organize_ok = simplified_pipeline._organize_output_artifacts()

        rmsd_result = (
            self._run_non_gnina_rmsd_bridge(engine_scores, best_poses, complexes_dir, bridge_out)
            if extracted_count
            else {"state": "skipped_no_best_pose_complexes", "processed_complexes": 0}
        )

        from post_docking_analysis.prolif_interaction_maps import PROLIF_AVAILABLE
        from post_docking_analysis.py3dmol_visualizer import PY3DMOL_AVAILABLE

        ligplot_root_configured = bool(simplified_pipeline._resolve_ligplus_root())
        poseview_enabled = bool(simplified_pipeline.poseview_config.get("enabled", False))
        if self.fast_mode:
            prolif_state = "skipped_fast_mode"
        else:
            prolif_state = (
                "completed"
                if isinstance(simplified_pipeline.results.get("prolif_interaction_maps"), dict)
                and simplified_pipeline.results.get("prolif_interaction_maps")
                else "missing_dependency"
                if not PROLIF_AVAILABLE
                else "failed_or_empty"
            )
        ligplot_summary = simplified_pipeline.results.get("ligplot_diagrams") or {}
        if self.fast_mode:
            ligplot_state = "skipped_fast_mode"
        else:
            ligplot_state = (
                "completed"
                if int(ligplot_summary.get("successful", 0)) > 0
                else "missing_configuration"
                if not ligplot_root_configured
                else "failed_or_empty"
            )
        poseview_summary = simplified_pipeline.results.get("poseview_summary") or {}
        if self.fast_mode:
            poseview_state = "skipped_fast_mode"
        else:
            poseview_state = (
                "completed"
                if int(poseview_summary.get("successful", 0)) > 0
                else "disabled"
                if not poseview_enabled
                else "failed_or_empty"
            )
        pandamap_summary = simplified_pipeline.results.get("pandamap_summary") or {}
        if self.fast_mode:
            pandamap_state = "skipped_fast_mode"
        else:
            pandamap_state = (
                "completed"
                if int(pandamap_summary.get("generated_2d_maps", 0)) > 0
                or int(pandamap_summary.get("generated_3d_visualizations", 0)) > 0
                else "failed_or_empty"
            )
        py3dmol_summary = simplified_pipeline.results.get("py3dmol_visualizations") or {}
        if self.fast_mode:
            py3dmol_state = "skipped_fast_mode"
        else:
            py3dmol_state = (
                "completed"
                if py3dmol_summary
                else "missing_dependency"
                if not PY3DMOL_AVAILABLE
                else "failed_or_empty"
            )
        if self.fast_mode:
            pymol_state = "skipped_fast_mode"
            quality_state = "skipped_fast_mode"
        else:
            # Retired legacy-route stages: structure-quality + PyMOL are no longer
            # executed through PostDockingAnalysisPipeline in bridge workflows.
            pymol_state = "skipped_retired_route"
            quality_state = "skipped_retired_route"
        simplified_visual_state = "skipped_fast_mode" if self.fast_mode else ("completed" if viz_ok else "skipped_or_failed")

        notes_lines = [
            f"Favorite engine: {favorite_engine}",
            f"Speed profile: {self.speed_profile}",
            f"Deep analysis output: {bridge_out}",
            "Continuation path: engine-aware structural bridge with simplified downstream reuse for PDBQT-based engines",
            f"Best-pose complexes extracted: {extracted_count}",
            f"Best-pose complexes mirrored to best_poses_pdb: {mirrored_best_pose_count}",
            f"Complex export index: {complex_index_file}",
            f"Simplified affinity analysis: {'completed' if affinity_ok else 'skipped_or_failed'}",
            f"Simplified polypharmacology: {'completed' if polypharm_ok else 'skipped_or_failed'}",
            f"Simplified pose extraction check: {'completed' if extract_ok else 'skipped_or_failed'}",
            f"Simplified reports: {'completed' if reports_ok else 'skipped_or_failed'}",
            f"Simplified visualizations: {simplified_visual_state}",
            f"RMSD analysis: {rmsd_result.get('state', 'failed_or_empty')}",
            f"RMSD per-complex sets processed: {int(rmsd_result.get('processed_complexes', 0) or 0)}",
            f"ProLIF interactions: {prolif_state}",
            f"LigPlot interactions: {ligplot_state}",
            f"PoseView interactions: {poseview_state}",
            f"PandaMap interactions: {pandamap_state}",
            f"py3Dmol visualizations: {py3dmol_state}",
            f"Artifact consolidation: {'completed' if organize_ok else 'skipped_or_failed'}",
            f"Structure-quality stage: {quality_state}",
            f"PyMOL visualization stage: {pymol_state}",
            fast_mode_note,
            "Unified score tables are mirrored under deep_analysis/raw_data/ for downstream compatibility.",
            "Non-GNINA RMSD uses PDBQT-derived complex PDB extractions for per-complex, per-protein-best, and global-best scopes.",
        ]
        notes_file.write_text(
            "\n".join([line for line in notes_lines if str(line).strip()]),
            encoding="utf-8",
        )

    def _run_downstream_bridge(self, favorite_engine: str, engine_scores: pd.DataFrame, output_dir: Path) -> None:
        notes_file = output_dir / "downstream_bridge.txt"
        if engine_scores.empty:
            notes_file.write_text(
                f"Favorite engine {favorite_engine} had no filtered scores to continue downstream.",
                encoding="utf-8",
            )
            return
        if favorite_engine != "gnina":
            self._run_non_gnina_downstream_bridge(favorite_engine, engine_scores, output_dir, notes_file)
            return

        engine_layout = ensure_engine_layout(self.project_dir, favorite_engine)
        if not list(engine_layout["poses"].glob("*.sdf")):
            notes_file.write_text(
                "GNINA favorite-engine continuation skipped because no SDF pose files were found.",
                encoding="utf-8",
            )
            return

        bridge_out = output_dir / "deep_analysis"
        from post_docking_analysis.simplified_pipeline import SimplifiedPostDockingPipeline

        pipeline = SimplifiedPostDockingPipeline(
            sdf_folder=str(engine_layout["poses"]),
            log_folder=str(engine_layout["logs"]),
            receptors_folder=str(shared_receptors_dir(self.project_dir)),
            output_dir=str(bridge_out),
            pairlist_file=str(pairlist_path(self.project_dir)),
            prompt_protein_names=self.prompt_protein_names,
            prompt_ligand_names=self.prompt_ligand_names,
            enable_poseview=self.enable_poseview,
            rmsd_scopes=self.rmsd_scopes,
            resume_rmsd=self.resume_rmsd,
            force_global_rmsd=self.force_global_rmsd,
            global_rmsd_defer_threshold=self.global_rmsd_defer_threshold,
            shared_mapping_dir=str(self.mapping_root),
            exclude_problematic_ligands=self.exclude_problematic_ligands,
            positive_affinity_threshold=self.positive_affinity_threshold,
            minimum_pose_count=self.minimum_pose_count,
            rmsd_workers=self.rmsd_workers,
        )
        success = pipeline.run()
        bridge_out = output_dir / "deep_analysis"
        complexes_dir = bridge_out / "complexes"
        best_poses_pdb_dir = bridge_out / "best_poses_pdb"
        if complexes_dir.exists():
            self._write_complex_export_index(
                output_dir=bridge_out,
                complexes_dir=complexes_dir,
                best_poses_dir=best_poses_pdb_dir if best_poses_pdb_dir.exists() else None,
                favorite_engine=favorite_engine,
                source_table=self._best_by_tag(engine_scores),
            )
        notes_file.write_text(
            f"GNINA downstream bridge {'completed' if success else 'failed'} at {bridge_out}",
            encoding="utf-8",
        )
