"""
Simplified Post-Docking Analysis Pipeline for GNINA.

Focuses on:
- 3-folder input structure (sdf_folder, log_folder, receptors_folder)
- Generate all_scores.csv from logs
- Create complex PDB files (receptor + ligand)
- Comparative benchmarking using PDB code matching + pairlist
- RMSD analysis
- Simplified visualizations (no PandaMap, no plugins)
"""
from pathlib import Path
from typing import Optional, Dict, List, Tuple, Any
from datetime import datetime, timezone
import subprocess
import pandas as pd
import logging
import shutil
import os
import sys
import json
import hashlib
import numpy as np
import re
import xml.etree.ElementTree as ET

from .generate_scores_csv import generate_all_scores_csv
from .complex_query import (
    filter_frame_by_complex_query,
    filter_paths_by_tags,
    tags_from_frame,
)
from .simplified_input_handler import (
    find_sdf_files,
    find_log_files,
    find_receptor_files,
    load_pairlist,
    match_poses_to_receptors,
    auto_detect_pairlist_file,
)
from .publication_pandamap import run_publication_pandamap_analysis
from .binding_affinity_analyzer import analyze_binding_affinities
from .hierarchical_analyzer import HierarchicalDockingAnalyzer
from .py3dmol_visualizer import (
    PY3DMOL_AVAILABLE,
    visualize_all_complexes,
    visualize_ligands_by_protein
)
from .prolif_interaction_maps import (
    PROLIF_AVAILABLE,
    ProLifInteractionMapper,
    create_interaction_maps_for_all_complexes,
)
from .plip_integration import PLIPAnalyzer
from .ligplot_integration import LigPlotRunner, export_ligplot_outputs
from .poseview_integration import run_poseview_analysis
from .ligand_naming import (
    build_ligand_name_mapping,
    ligand_mapping_to_dataframe,
    resolve_ligand_display_name,
)
from .protein_naming import (
    build_protein_name_mapping,
    mapping_to_dataframe,
    resolve_protein_display_name,
    extract_pdb_code,
    format_protein_label,
)


_VALID_RMSD_SCOPES = ("per_complex", "per_protein", "global")
_DEFAULT_RMSD_SCOPES = ("per_complex",)
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
_RUN_TRACKING_OPTIONAL_STATUS_VALUES = {
    "completed",
    "skipped_disabled",
    "skipped_missing_dependency",
    "failed_error",
}
_OPTIONAL_STEP_LABELS = {
    "poseview": "Generate PoseView diagrams",
    "pandamap": "Generate PandaMap analysis",
    "py3dmol": "Generate py3Dmol visualizations",
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
            return _VALID_RMSD_SCOPES
        if value in _VALID_RMSD_SCOPES and value not in normalized:
            normalized.append(value)
    return tuple(normalized) if normalized else _DEFAULT_RMSD_SCOPES


def _optional_feature_record(
    status: str,
    reason: str,
    *,
    step_status: str = "",
    error: str = "",
) -> Dict[str, str]:
    normalized_status = status if status in _RUN_TRACKING_OPTIONAL_STATUS_VALUES else "failed_error"
    return {
        "status": normalized_status,
        "reason": str(reason or ""),
        "step_status": str(step_status or ""),
        "error": str(error or ""),
    }


class SimplifiedPostDockingPipeline:
    """
    Simplified post-docking analysis pipeline for GNINA.
    
    Input Structure:
    - sdf_folder: Docking poses (SDF files)
    - log_folder: Docking logs (can be same as sdf_folder)
    - receptors_folder: Receptor PDBQT files
    - pairlist_file: Optional pairlist.csv for matching
    """
    
    def __init__(
        self,
        sdf_folder: str,
        log_folder: str,
        receptors_folder: str,
        output_dir: str,
        pairlist_file: Optional[str] = None,
        run_rmsd: bool = True,
        run_visualizations: bool = True,
        ligplus_root: Optional[str] = None,
        prompt_protein_names: bool = False,
        prompt_ligand_names: bool = False,
        enable_poseview: bool = False,
        complex_query: Optional[str] = None,
        rmsd_scopes: Optional[object] = None,
        resume_rmsd: bool = True,
        force_global_rmsd: bool = False,
        global_rmsd_defer_threshold: int = 150,
        shared_mapping_dir: Optional[str] = None,
        exclude_problematic_ligands: bool = False,
        positive_affinity_threshold: float = 0.0,
        minimum_pose_count: int = 1,
        rmsd_workers: int = 0,
        analysis_config: Optional[Dict[str, Any]] = None,
    ):
        """
        Initialize simplified pipeline.
        
        Parameters
        ----------
        sdf_folder : str
            Folder containing SDF pose files
        log_folder : str
            Folder containing log files (can be same as sdf_folder)
        receptors_folder : str
            Folder containing receptor PDBQT files
        output_dir : str
            Output directory for results
        pairlist_file : str, optional
            Path to pairlist.csv
        run_rmsd : bool
            Compatibility toggle (RMSD stage is now always enforced)
        run_visualizations : bool
            Compatibility toggle (core visualization/interaction stages are now always enforced)
        ligplus_root : str, optional
            Path to LigPlus root (optional). If omitted, pipeline checks
            LIGPLUS_ROOT / LIGPLUS_HOME.
        """
        self.sdf_folder = Path(sdf_folder)
        self.log_folder = Path(log_folder)
        self.receptors_folder = Path(receptors_folder)
        self.output_dir = Path(output_dir)
        self.pairlist_file = Path(pairlist_file) if pairlist_file else None
        self.requested_run_rmsd = bool(run_rmsd)
        self.requested_run_visualizations = bool(run_visualizations)
        self.run_rmsd = True
        self.run_visualizations = True
        self.ligplus_root = Path(ligplus_root).expanduser() if ligplus_root else None
        self.prompt_protein_names = prompt_protein_names
        self.prompt_ligand_names = prompt_ligand_names
        self.enable_poseview = enable_poseview
        self.complex_query = str(complex_query or "").strip()
        self.selected_query_tags = set()
        self.analysis_pairlist_file = self.pairlist_file
        self.shared_mapping_dir = (
            Path(shared_mapping_dir).expanduser().resolve() if shared_mapping_dir else None
        )
        if self.shared_mapping_dir is not None:
            self.shared_mapping_dir.mkdir(parents=True, exist_ok=True)
        self.rmsd_scopes = _normalize_rmsd_scopes(rmsd_scopes)
        self.resume_rmsd = bool(resume_rmsd)
        self.force_global_rmsd = bool(force_global_rmsd)
        self.exclude_problematic_ligands = bool(exclude_problematic_ligands)
        self.positive_affinity_threshold = float(positive_affinity_threshold)
        self.minimum_pose_count = max(int(minimum_pose_count or 1), 1)
        self.rmsd_workers = int(rmsd_workers or 0)
        try:
            self.global_rmsd_defer_threshold = max(int(global_rmsd_defer_threshold or 0), 0)
        except (TypeError, ValueError):
            self.global_rmsd_defer_threshold = 150
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.visualizations_root = self.output_dir / "visualizations"
        self.visual_2d_root = self.visualizations_root / "2d"
        self.visual_3d_root = self.visualizations_root / "3d"
        self.visual_interactions_2d_root = self.visual_2d_root / "interactions"
        self.visual_interactions_3d_root = self.visual_3d_root / "interactions"
        self.raw_data_root = self.output_dir / "raw_data"
        self.interactions_root = self.output_dir / "interactions"
        self.run_tracking_dir = self.output_dir / "run_tracking"
        self.run_manifest_file = self.run_tracking_dir / "run_manifest.json"

        # Setup logging
        log_file = self.output_dir / "simplified_pipeline.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        if not self.requested_run_rmsd:
            self.logger.warning("⚠️  Ignoring compatibility flag: RMSD stage is mandatory in clean interaction mode")
        if not self.requested_run_visualizations:
            self.logger.warning(
                "⚠️  Ignoring compatibility flag: core visualization/interaction stages are mandatory in clean interaction mode"
            )
        
        # Results storage
        self.results = {}
        self.complexes = []
        self.scores_df = None
        self.pairlist_df = pd.DataFrame()
        self.protein_name_map: Dict[str, Dict[str, str]] = {}
        self.ligand_name_map: Dict[str, str] = {}
        self.mapping_root = self.shared_mapping_dir or self.output_dir
        self.protein_name_override_file = self.mapping_root / "protein_name_overrides.csv"
        self.ligand_name_override_file = self.mapping_root / "ligand_name_overrides.csv"
        self.protein_name_mapping_file = self.mapping_root / "protein_name_mapping.csv"
        self.ligand_name_mapping_file = self.mapping_root / "ligand_name_mapping.csv"

        # Tuned defaults for visualization engines.
        self.pandamap_config = {
            "dpi": 350,
            "formats": ["pdf", "svg", "png"],
            "figure_width": 12,
            "figure_height": 9,
            "show_surface": True,
            "show_3d_cues": True,
            "color_by": "interaction",
            "size_scale": 1.35,
            "max_complexes": 30,
            "overwrite": False,
            "qc_min_non_generic_ligand_ratio": 0.60,
            "qc_require_reports_when_supported": True,
        }
        self.prolif_config = {
            "dpi": 350,
            "figsize": (14, 10),
            "max_complexes": 60,
            "write_html": True,
            # Docking-mode settings: analyze interactions across poses from each SDF.
            "poses_per_complex": None,
            "lignetwork_threshold": 0.3,
            "count_occurrences": False,
        }
        self.plip_config = {
            "max_complexes": 60,
            "max_workers": 4,
            "output_formats": ["png", "xml", "txt"],
            "generate_pymol_session": False,
            "timeout_seconds": 180,
        }
        self.ligplot_config = {
            "contact_type": "2",
            "no_abort": True,
            "hydrogenate": "auto",
            "strip_metals": "auto",
            "max_complexes": 40,
            "overwrite": False,
            "export_formats": ["png", "pdf"],
            "export_dpi": 350,
        }
        self.poseview_config = {
            "enabled": bool(enable_poseview),
            "max_complexes": 24,
            "output_formats": ["png", "svg", "pdf"],
            "poll_interval_seconds": 5,
            "timeout_seconds": 360,
            "max_retries": 2,
            "upload_poll_interval_seconds": 2,
            "upload_timeout_seconds": 180,
            "max_ligand_candidates": 8,
        }
        self.polypharm_config = {
            "enabled": True,
            "strong_threshold": -7.0,
            "moderate_threshold": -5.0,
            "near_reference_delta": 1.0,
            "min_targets_moderate": 2,
            "min_targets_high": 3,
            "exclude_reference_ligands": True,
            "top_ligands_plot_count": 20,
        }
        self.interaction_analytics_config = {
            "enabled": True,
            "chain_integrity_auto_repair": True,
        }
        self.analysis_config_overrides = analysis_config if isinstance(analysis_config, dict) else {}
        self._apply_analysis_config_overrides(self.analysis_config_overrides)

    @staticmethod
    def _utc_now_iso() -> str:
        return datetime.now(timezone.utc).isoformat()

    def _apply_analysis_config_overrides(self, overrides: Dict[str, Any]) -> None:
        """
        Apply config-file overrides to simplified analysis stage settings.
        """
        if not isinstance(overrides, dict) or not overrides:
            return

        visualization = overrides.get("visualization", {})
        if isinstance(visualization, dict):
            for key, value in (visualization.get("pandamap", {}) or {}).items():
                self.pandamap_config[key] = value
            for key, value in (visualization.get("prolif", {}) or {}).items():
                self.prolif_config[key] = value
            for key, value in (visualization.get("plip", {}) or {}).items():
                self.plip_config[key] = value
            for key, value in (visualization.get("ligplot", {}) or {}).items():
                self.ligplot_config[key] = value
            for key, value in (visualization.get("poseview", {}) or {}).items():
                self.poseview_config[key] = value

        polypharm = overrides.get("polypharmacology", {})
        if isinstance(polypharm, dict):
            for key, value in polypharm.items():
                self.polypharm_config[key] = value

        interaction_analytics = overrides.get("interaction_analytics", {})
        if isinstance(interaction_analytics, dict):
            for key, value in interaction_analytics.items():
                self.interaction_analytics_config[key] = value

    @staticmethod
    def _json_safe(value: Any) -> Any:
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, dict):
            return {str(key): SimplifiedPostDockingPipeline._json_safe(val) for key, val in value.items()}
        if isinstance(value, (list, tuple, set)):
            return [SimplifiedPostDockingPipeline._json_safe(item) for item in value]
        if isinstance(value, (str, int, float, bool)) or value is None:
            return value
        return str(value)

    def _classify_output_path(self, relative_path: Path) -> str:
        if not relative_path.parts:
            return "root"
        top_level = relative_path.parts[0].lower()
        known = {
            "best_poses",
            "complexes",
            "interactions",
            "raw_data",
            "reports",
            "run_tracking",
            "visualizations",
        }
        if top_level in known:
            return top_level
        extension = relative_path.suffix.lower()
        if extension in {".csv", ".tsv", ".json", ".txt", ".xlsx"}:
            return "tables_and_metadata"
        if extension in {".pdb", ".pdbqt", ".sdf", ".mol2"}:
            return "structures"
        if extension in {".png", ".svg", ".pdf", ".html"}:
            return "figures"
        if extension in {".log"}:
            return "logs"
        return "other"

    def _write_run_tracking_manifest(
        self,
        *,
        run_started_at: str,
        run_status: str,
        error_message: str,
        step_states: List[Dict[str, Any]],
    ) -> None:
        try:
            self.run_tracking_dir.mkdir(parents=True, exist_ok=True)

            step_rows: List[Dict[str, Any]] = []
            for state in step_states:
                step_rows.append(
                    {
                        "index": int(state.get("index", 0) or 0),
                        "step": str(state.get("step", "") or ""),
                        "required": bool(state.get("required", False)),
                        "status": str(state.get("status", "") or ""),
                        "started_at": str(state.get("started_at", "") or ""),
                        "ended_at": str(state.get("ended_at", "") or ""),
                        "details": str(state.get("details", "") or ""),
                        "error": str(state.get("error", "") or ""),
                    }
                )

            step_status_df = pd.DataFrame(
                step_rows,
                columns=_RUN_TRACKING_STEP_COLUMNS,
            )
            step_status_csv = self.run_tracking_dir / "step_status.csv"
            step_status_json = self.run_tracking_dir / "step_status.json"
            step_status_df.to_csv(step_status_csv, index=False)
            step_status_df.to_json(step_status_json, orient="records", indent=2)

            output_rows: List[Dict[str, Any]] = []
            for file_path in sorted(self.output_dir.rglob("*")):
                if not file_path.is_file():
                    continue
                relative = file_path.relative_to(self.output_dir)
                if relative.parts and relative.parts[0] == "run_tracking":
                    continue
                stat = file_path.stat()
                output_rows.append(
                    {
                        "relative_path": relative.as_posix(),
                        "category": self._classify_output_path(relative),
                        "extension": str(file_path.suffix.lower()),
                        "size_bytes": int(stat.st_size),
                        "modified_utc": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
                    }
                )
            outputs_df = pd.DataFrame(
                output_rows,
                columns=_RUN_TRACKING_OUTPUT_COLUMNS,
            )
            outputs_index_csv = self.run_tracking_dir / "outputs_index.csv"
            outputs_index_json = self.run_tracking_dir / "outputs_index.json"
            outputs_df.to_csv(outputs_index_csv, index=False)
            outputs_df.to_json(outputs_index_json, orient="records", indent=2)

            status_counts: Dict[str, int] = {}
            for row in step_rows:
                status = str(row.get("status", "") or "")
                if not status:
                    continue
                status_counts[status] = status_counts.get(status, 0) + 1

            optional_features = self._infer_optional_feature_statuses(step_rows)

            run_manifest = {
                "pipeline": "simplified_gnina",
                "input": {
                    "sdf_folder": str(self.sdf_folder),
                    "log_folder": str(self.log_folder),
                    "receptors_folder": str(self.receptors_folder),
                    "pairlist_file": str(self.pairlist_file) if self.pairlist_file else "",
                },
                "output_dir": str(self.output_dir),
                "run_started_at": run_started_at,
                "run_completed_at": self._utc_now_iso(),
                "run_status": str(run_status),
                "error_message": str(error_message or ""),
                "stage_contract": {
                    "requested_run_rmsd": bool(self.requested_run_rmsd),
                    "requested_run_visualizations": bool(self.requested_run_visualizations),
                    "effective_run_rmsd": bool(self.run_rmsd),
                    "effective_run_visualizations": bool(self.run_visualizations),
                    "requested_run_prolif": True,
                    "requested_run_ligplot": True,
                    "effective_run_prolif": True,
                    "effective_run_ligplot": True,
                },
                "rmsd": {
                    "scopes": list(self.rmsd_scopes),
                    "resume": bool(self.resume_rmsd),
                    "force_global": bool(self.force_global_rmsd),
                    "global_defer_threshold": int(self.global_rmsd_defer_threshold),
                    "workers": int(self.rmsd_workers),
                },
                "poseview_enabled": bool(self.poseview_config.get("enabled", False)),
                "step_status_file": str(step_status_csv),
                "outputs_index_file": str(outputs_index_csv),
                "steps_total": int(len(step_rows)),
                "step_status_counts": status_counts,
                "outputs_count": int(len(output_rows)),
                "optional_features": optional_features,
                "results": self._json_safe(self.results),
            }

            self.run_manifest_file.write_text(
                json.dumps(run_manifest, indent=2),
                encoding="utf-8",
            )
        except Exception as exc:
            self.logger.warning(f"⚠️  Failed to write run-tracking manifest: {exc}")

    def _infer_optional_feature_statuses(
        self,
        step_rows: List[Dict[str, Any]],
    ) -> Dict[str, Dict[str, str]]:
        step_index = {
            str(row.get("step", "") or ""): row
            for row in step_rows
        }
        complex_count = 0
        try:
            complex_count = int(len(self._get_complex_pdb_files()))
        except Exception:
            complex_count = 0

        poseview_enabled = bool(self.poseview_config.get("enabled", False))
        poseview_summary = self.results.get("poseview_summary")
        poseview_successful = 0
        if isinstance(poseview_summary, dict):
            poseview_successful = int(poseview_summary.get("successful", 0) or 0)

        pandamap_summary = self.results.get("pandamap_summary")
        pandamap_generated = 0
        if isinstance(pandamap_summary, dict):
            pandamap_generated = int(pandamap_summary.get("generated_2d_maps", 0) or 0) + int(
                pandamap_summary.get("generated_3d_visualizations", 0) or 0
            )

        py3dmol_payload = self.results.get("py3dmol_visualizations")
        py3dmol_has_outputs = bool(py3dmol_payload)

        poseview_row = step_index.get(_OPTIONAL_STEP_LABELS["poseview"], {})
        pandamap_row = step_index.get(_OPTIONAL_STEP_LABELS["pandamap"], {})
        py3dmol_row = step_index.get(_OPTIONAL_STEP_LABELS["py3dmol"], {})

        poseview_step_status = str(poseview_row.get("status", "") or "")
        poseview_error = str(poseview_row.get("error", "") or "")
        pandamap_step_status = str(pandamap_row.get("status", "") or "")
        pandamap_error = str(pandamap_row.get("error", "") or "")
        py3dmol_step_status = str(py3dmol_row.get("status", "") or "")
        py3dmol_error = str(py3dmol_row.get("error", "") or "")

        poseview_record = _optional_feature_record(
            "failed_error",
            "poseview_stage_failed",
            step_status=poseview_step_status,
            error=poseview_error,
        )
        if not poseview_enabled:
            poseview_record = _optional_feature_record(
                "skipped_disabled",
                "poseview_disabled",
                step_status=poseview_step_status,
                error=poseview_error,
            )
        elif complex_count <= 0:
            poseview_record = _optional_feature_record(
                "skipped_disabled",
                "no_complexes_available",
                step_status=poseview_step_status,
                error=poseview_error,
            )
        elif poseview_successful > 0:
            poseview_record = _optional_feature_record(
                "completed",
                "poseview_outputs_generated",
                step_status=poseview_step_status,
                error=poseview_error,
            )
        elif poseview_error:
            poseview_record = _optional_feature_record(
                "failed_error",
                "poseview_step_error",
                step_status=poseview_step_status,
                error=poseview_error,
            )

        pandamap_record = _optional_feature_record(
            "failed_error",
            "pandamap_stage_failed",
            step_status=pandamap_step_status,
            error=pandamap_error,
        )
        if complex_count <= 0:
            pandamap_record = _optional_feature_record(
                "skipped_disabled",
                "no_complexes_available",
                step_status=pandamap_step_status,
                error=pandamap_error,
            )
        elif pandamap_generated > 0:
            pandamap_record = _optional_feature_record(
                "completed",
                "pandamap_outputs_generated",
                step_status=pandamap_step_status,
                error=pandamap_error,
            )
        elif pandamap_step_status in {"completed", "failed"} and pandamap_error:
            pandamap_record = _optional_feature_record(
                "failed_error",
                "pandamap_step_error",
                step_status=pandamap_step_status,
                error=pandamap_error,
            )

        py3dmol_record = _optional_feature_record(
            "failed_error",
            "py3dmol_stage_failed",
            step_status=py3dmol_step_status,
            error=py3dmol_error,
        )
        if not PY3DMOL_AVAILABLE:
            py3dmol_record = _optional_feature_record(
                "skipped_missing_dependency",
                "py3dmol_not_installed",
                step_status=py3dmol_step_status,
                error=py3dmol_error,
            )
        elif complex_count <= 0:
            py3dmol_record = _optional_feature_record(
                "skipped_disabled",
                "no_complexes_available",
                step_status=py3dmol_step_status,
                error=py3dmol_error,
            )
        elif py3dmol_has_outputs:
            py3dmol_record = _optional_feature_record(
                "completed",
                "py3dmol_outputs_generated",
                step_status=py3dmol_step_status,
                error=py3dmol_error,
            )
        elif py3dmol_error:
            py3dmol_record = _optional_feature_record(
                "failed_error",
                "py3dmol_step_error",
                step_status=py3dmol_step_status,
                error=py3dmol_error,
            )

        return {
            "poseview": poseview_record,
            "pandamap": pandamap_record,
            "py3dmol": py3dmol_record,
        }
        
    def run(self) -> bool:
        """
        Run the simplified pipeline.
        
        Returns
        -------
        bool
            True if successful, False otherwise
        """
        self.logger.info("🚀 Starting Simplified Post-Docking Analysis Pipeline")
        self.logger.info(f"📂 SDF folder: {self.sdf_folder}")
        self.logger.info(f"📂 Log folder: {self.log_folder}")
        self.logger.info(f"📂 Receptors folder: {self.receptors_folder}")
        self.logger.info(f"📂 Output directory: {self.output_dir}")
        if self.run_rmsd:
            self.logger.info(
                "📂 RMSD scopes: %s | resume=%s | force_global=%s | defer_threshold=%s",
                ",".join(self.rmsd_scopes),
                self.resume_rmsd,
                self.force_global_rmsd,
                self.global_rmsd_defer_threshold,
            )
        run_started_at = self._utc_now_iso()
        run_status = "running"
        error_message = ""
        step_states: List[Dict[str, Any]] = []

        try:
            step_plan: List[Tuple[str, object, bool]] = [
                ("Find input files", self._find_input_files, True),
                ("Generate all_scores.csv", self._generate_scores_csv, True),
                ("Match poses to receptors", self._match_poses_to_receptors, True),
                ("Receptor chain-integrity QC + repair", self._receptor_chain_integrity_qc_and_repair, True),
                ("Create complex PDB files", self._create_complexes, True),
                ("Binding affinity analysis", self._analyze_binding_affinity, True),
                ("Polypharmacology analysis", self._analyze_polypharmacology, True),
                ("RMSD analysis", self._analyze_rmsd, True),
                ("Extract and organize poses", self._extract_poses, True),
                ("Generate reports", self._generate_reports, True),
                ("Generate summary visualizations", self._generate_visualizations, True),
                ("Generate PLIP interaction outputs", self._generate_plip_interaction_outputs, True),
                ("Generate ProLIF interaction maps", self._generate_prolif_interaction_maps, True),
                ("Generate LigPlot+ diagrams", self._generate_ligplot_diagrams, True),
                ("Normalize layered PLIP/ProLIF tables", self._normalize_interaction_tables, True),
                ("Build source-of-truth interaction analytics", self._build_source_of_truth_interaction_analytics, True),
                ("Generate PoseView diagrams", self._generate_poseview_diagrams, False),
                ("Generate PandaMap analysis", self._generate_pandamap_analysis, False),
                ("Generate py3Dmol visualizations", self._generate_py3dmol_visualizations, False),
                ("Consolidate output artifacts", self._organize_output_artifacts, True),
            ]

            total_steps = len(step_plan)
            for index, (step_label, step_func, required) in enumerate(step_plan, start=1):
                self.logger.info(f"📍 Step [{index}/{total_steps}] {step_label}")
                step_started_at = self._utc_now_iso()
                step_error = ""
                try:
                    ok = bool(step_func())
                except Exception as exc:
                    ok = False
                    step_error = str(exc)
                    self.logger.error(f"❌ Step exception [{index}/{total_steps}]: {step_label} -> {exc}", exc_info=True)

                status_text = "completed" if ok else ("failed" if required else "skipped_or_unavailable")
                status_icon = "✅" if ok else ("❌" if required else "⚠️")
                step_states.append(
                    {
                        "index": index,
                        "step": step_label,
                        "required": required,
                        "status": status_text,
                        "started_at": step_started_at,
                        "ended_at": self._utc_now_iso(),
                        "error": step_error,
                    }
                )
                self._write_run_tracking_manifest(
                    run_started_at=run_started_at,
                    run_status=run_status,
                    error_message=error_message,
                    step_states=step_states,
                )
                if required and not ok:
                    self.logger.error(f"❌ Step failed [{index}/{total_steps}]: {step_label}")
                    error_message = step_error or f"{step_label} returned an unsuccessful status"
                    run_status = "failed"
                    return False
                self.logger.info(f"{status_icon} Step [{index}/{total_steps}] {step_label}: {status_text}")

            run_status = "completed"
            self.logger.info("✅ Simplified pipeline completed successfully!")
            return True
            
        except Exception as e:
            error_message = str(e)
            self.logger.error(f"❌ Pipeline failed: {e}", exc_info=True)
            run_status = "failed"
            return False
        finally:
            self._write_run_tracking_manifest(
                run_started_at=run_started_at,
                run_status=run_status,
                error_message=error_message,
                step_states=step_states,
            )
    
    def _find_input_files(self) -> bool:
        """Find and validate input files."""
        self.logger.info("🔍 Finding input files...")
        
        self.sdf_files = find_sdf_files(self.sdf_folder)
        self.log_files = find_log_files(self.log_folder)
        self.receptor_files = find_receptor_files(self.receptors_folder)

        if self.pairlist_file is not None and not self.pairlist_file.exists():
            self.logger.warning(f"⚠️  Pairlist not found at provided path: {self.pairlist_file}")
            self.pairlist_file = None

        if self.pairlist_file is None:
            detected_pairlist = auto_detect_pairlist_file(
                self.sdf_folder,
                self.log_folder,
                self.receptors_folder,
            )
            if detected_pairlist is not None:
                self.pairlist_file = detected_pairlist
                self.logger.info(f"🔍 Auto-detected pairlist.csv: {self.pairlist_file}")

        self.pairlist_df = load_pairlist(self.pairlist_file)

        if self.complex_query and not self.pairlist_df.empty:
            original_pair_count = len(self.pairlist_df)
            self.pairlist_df = filter_frame_by_complex_query(self.pairlist_df, self.complex_query)
            if self.pairlist_df.empty:
                self.logger.error("❌ Complex query removed all pairlist rows")
                return False
            self.selected_query_tags = tags_from_frame(self.pairlist_df)
            self.sdf_files = filter_paths_by_tags(self.sdf_files, self.selected_query_tags)
            self.log_files = filter_paths_by_tags(self.log_files, self.selected_query_tags)
            selected_receptors = {
                Path(str(value)).stem.lower()
                for value in self.pairlist_df.get("receptor", pd.Series(dtype="object")).dropna().astype(str).tolist()
            }
            if selected_receptors:
                self.receptor_files = [
                    receptor_file
                    for receptor_file in self.receptor_files
                    if Path(receptor_file).stem.lower() in selected_receptors
                ]
            self.logger.info(
                f"🔎 Applied complex query to pairlist: {original_pair_count} -> {len(self.pairlist_df)} rows"
            )
            (self.output_dir / "complex_query.txt").write_text(self.complex_query + "\n", encoding="utf-8")
            filtered_pairlist_path = self.output_dir / "filtered_pairlist.csv"
            self.pairlist_df.to_csv(filtered_pairlist_path, index=False)
            self.analysis_pairlist_file = filtered_pairlist_path
            (self.output_dir / "complex_query_summary.txt").write_text(
                (
                    f"query={self.complex_query}\n"
                    f"original_pairlist_rows={original_pair_count}\n"
                    f"filtered_pairlist_rows={len(self.pairlist_df)}\n"
                    f"filtered_sdf_files={len(self.sdf_files)}\n"
                    f"filtered_log_files={len(self.log_files)}\n"
                    f"filtered_receptors={len(self.receptor_files)}\n"
                ),
                encoding="utf-8",
            )
        else:
            self.analysis_pairlist_file = self.pairlist_file

        self.protein_name_map = build_protein_name_mapping(
            self.receptor_files,
            self.pairlist_df,
            overrides_file=self.protein_name_override_file
        )
        self._prompt_for_protein_name_overrides()
        self._persist_protein_name_mapping()
        self.ligand_name_map = self._build_ligand_name_map()
        self._prompt_for_ligand_name_overrides()
        self._persist_ligand_name_mapping()

        empty_sdf_files = [f for f in self.sdf_files if f.exists() and f.stat().st_size == 0]
        if empty_sdf_files:
            self.logger.warning(
                f"⚠️  Skipping {len(empty_sdf_files)} empty SDF files "
                "(likely timed out or failed docking outputs)"
            )
            skipped_file = self.output_dir / "skipped_empty_sdf_files.txt"
            with open(skipped_file, "w") as handle:
                for sdf_path in empty_sdf_files:
                    handle.write(f"{sdf_path}\n")
            self.logger.info(f"  📝 Saved skipped empty SDF list: {skipped_file}")
            self.sdf_files = [f for f in self.sdf_files if f.stat().st_size > 0]
        
        self.logger.info(f"  Found {len(self.sdf_files)} SDF files")
        self.logger.info(f"  Found {len(self.log_files)} log files")
        self.logger.info(f"  Found {len(self.receptor_files)} receptor files")

        if not self.sdf_files:
            self.logger.error("❌ No valid (non-empty) SDF files found!")
            return False
        
        if not self.log_files:
            self.logger.error("❌ No log files found!")
            return False
        
        if not self.receptor_files:
            self.logger.warning("⚠️  No receptor files found!")
        
        return True

    def _prompt_for_protein_name_overrides(self) -> None:
        """
        Prompt for protein display-name overrides (one prompt per detected target).
        """
        if not self.prompt_protein_names:
            return
        if not sys.stdin.isatty():
            self.logger.warning("⚠️  --prompt-protein-names requested but no interactive TTY detected; skipping prompts")
            return

        code_map = self.protein_name_map.get("pdb_code_to_name", {})
        receptor_map = self.protein_name_map.get("receptor_to_name", {})
        if not code_map and not self.receptor_files:
            return

        self.logger.info(
            "📝 Protein naming prompts enabled: press Enter to keep each detected default name."
        )

        seen_codes = set()
        seen_receptors = set()
        prompt_targets: List[Tuple[str, str, str]] = []

        # Prefer one prompt per detected receptor file in this run.
        for receptor_file in self.receptor_files:
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
                    or resolve_protein_display_name(receptor_name, self.protein_name_map)
                ).strip()
                prompt_targets.append(("receptor", receptor_stem, current or receptor_stem))

        # Fallback when receptor files are unavailable but code map exists.
        if not prompt_targets:
            for code in sorted(code_map):
                current = str(code_map.get(code, "")).strip() or f"Protein {code}"
                prompt_targets.append(("code", code, current))

        for target_type, target_id, current in prompt_targets:
            if target_type == "code":
                prompt = f"  Protein name for PDB {target_id} [{current}]: "
            else:
                prompt = f"  Protein name for receptor {target_id} [{current}]: "
            try:
                updated = input(prompt).strip()
            except Exception:
                updated = ""
            final_name = updated or current

            if target_type == "code":
                code_map[target_id] = final_name
                # Re-apply receptor names from updated code map when possible.
                for receptor_key in list(receptor_map.keys()):
                    key_code = extract_pdb_code(receptor_key)
                    if key_code and key_code == target_id:
                        receptor_map[receptor_key] = final_name
            else:
                receptor_map[target_id] = final_name
                receptor_map[Path(target_id).stem] = final_name

    def _persist_protein_name_mapping(self) -> None:
        """Save auto-detected naming map and a user-editable override template."""
        map_df = mapping_to_dataframe(self.protein_name_map)
        if not map_df.empty:
            map_file = self.protein_name_mapping_file
            map_df.to_csv(map_file, index=False)
            output_map_file = self.output_dir / "protein_name_mapping.csv"
            if output_map_file != map_file:
                map_df.to_csv(output_map_file, index=False)
            self.logger.info(f"  🏷️  Saved protein name mapping: {map_file}")

        # Create a simple override template once; users can edit and rerun.
        if not self.protein_name_override_file.exists():
            rows = []
            for pdb_code, display_name in self.protein_name_map.get("pdb_code_to_name", {}).items():
                rows.append(
                    {
                        "pdb_code": pdb_code,
                        "receptor": "",
                        "display_name": display_name,
                        "notes": "Edit display_name if needed; keep pdb_code unchanged.",
                    }
                )
            pd.DataFrame(rows, columns=["pdb_code", "receptor", "display_name", "notes"]).to_csv(
                self.protein_name_override_file,
                index=False
            )
            self.logger.info(f"  📝 Created protein name override template: {self.protein_name_override_file}")

    def _build_ligand_name_map(self) -> Dict[str, str]:
        """Build ligand display-name mapping from pairlist and known defaults."""
        identifiers: List[str] = []
        if isinstance(self.pairlist_df, pd.DataFrame) and not self.pairlist_df.empty:
            for column in ("ligand", "ligand_name", "cocrystal_ligand_name"):
                if column in self.pairlist_df.columns:
                    identifiers.extend(self.pairlist_df[column].dropna().astype(str).tolist())
        return build_ligand_name_mapping(
            identifiers,
            pairlist_df=self.pairlist_df,
            overrides_file=self.ligand_name_override_file,
        )

    def _prompt_for_ligand_name_overrides(self) -> None:
        """Prompt for ligand display-name overrides when explicitly enabled."""
        if not self.prompt_ligand_names:
            return
        if not sys.stdin.isatty():
            self.logger.warning("⚠️  --prompt-ligand-names requested but no interactive TTY detected; skipping prompts")
            return
        if not self.ligand_name_map:
            return

        self.logger.info(
            "📝 Ligand naming prompts enabled: press Enter to keep each detected default name."
        )
        seen = set()
        prompt_targets: List[Tuple[str, str]] = []
        for source in sorted(self.ligand_name_map):
            stem = Path(str(source)).stem
            if source != stem:
                continue
            display = str(self.ligand_name_map.get(source, "")).strip()
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
            self.ligand_name_map[ligand_id] = final_name
            self.ligand_name_map[f"{ligand_id}.pdbqt"] = final_name

    def _persist_ligand_name_mapping(self) -> None:
        """Save ligand naming map and a user-editable override template."""
        map_df = ligand_mapping_to_dataframe(self.ligand_name_map)
        if not map_df.empty:
            map_file = self.ligand_name_mapping_file
            map_df.to_csv(map_file, index=False)
            output_map_file = self.output_dir / "ligand_name_mapping.csv"
            if output_map_file != map_file:
                map_df.to_csv(output_map_file, index=False)
            self.logger.info(f"  🏷️  Saved ligand name mapping: {map_file}")

        if not self.ligand_name_override_file.exists():
            rows = []
            seen_stems = set()
            for key, display_name in sorted(self.ligand_name_map.items()):
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
                self.ligand_name_override_file,
                index=False,
            )
            self.logger.info(f"  📝 Created ligand name override template: {self.ligand_name_override_file}")
    
    def _generate_scores_csv(self) -> bool:
        """Generate all_scores.csv from log files."""
        self.logger.info("📊 Generating all_scores.csv...")
        
        scores_csv = self.output_dir / "all_scores.csv"
        
        # Use log_folder as gnina_out_dir for generate_all_scores_csv
        success = generate_all_scores_csv(
            self.log_folder,
            scores_csv,
            self.analysis_pairlist_file,
            log_files=self.log_files,
        )
        
        if success and scores_csv.exists():
            self.scores_df = pd.read_csv(scores_csv)
            if self.complex_query:
                original_score_count = len(self.scores_df)
                if self.selected_query_tags:
                    tag_series = self.scores_df.get("tag", pd.Series([""] * len(self.scores_df), index=self.scores_df.index)).astype(str)
                    self.scores_df = self.scores_df.loc[tag_series.isin(self.selected_query_tags)].copy()
                else:
                    self.scores_df = filter_frame_by_complex_query(self.scores_df, self.complex_query)
                if self.scores_df.empty:
                    self.logger.error("❌ Complex query removed all score rows")
                    return False
                score_tags = tags_from_frame(self.scores_df)
                self.selected_query_tags = (
                    self.selected_query_tags.intersection(score_tags)
                    if self.selected_query_tags
                    else score_tags
                )
                if self.complex_query and not self.selected_query_tags:
                    self.logger.error("❌ Complex query matched pairlist rows, but none of those tags had usable score rows")
                    return False
                if self.selected_query_tags:
                    self.sdf_files = filter_paths_by_tags(self.sdf_files, self.selected_query_tags)
                    self.log_files = filter_paths_by_tags(self.log_files, self.selected_query_tags)
                    if not self.pairlist_df.empty:
                        pair_tags = tags_from_frame(self.pairlist_df)
                        keep_tags = pair_tags.intersection(self.selected_query_tags)
                        if not keep_tags:
                            self.logger.error("❌ Complex query selected score rows, but no canonical pairlist tags overlap after filtering")
                            return False
                        pair_tag_series = self.pairlist_df.apply(
                            lambda row: f"{row['receptor']}_{row['site_id']}_{row['ligand']}",
                            axis=1,
                        ).astype(str)
                        self.pairlist_df = self.pairlist_df.loc[pair_tag_series.isin(keep_tags)].copy()
                        if self.analysis_pairlist_file is not None:
                            self.pairlist_df.to_csv(self.analysis_pairlist_file, index=False)
                self.scores_df.to_csv(scores_csv, index=False)
                summary_lines = [
                    f"query={self.complex_query}",
                ]
                summary_path = self.output_dir / "complex_query_summary.txt"
                if summary_path.exists():
                    existing = {
                        line.split("=", 1)[0]: line.split("=", 1)[1]
                        for line in summary_path.read_text(encoding="utf-8").splitlines()
                        if "=" in line
                    }
                    for key in (
                        "original_pairlist_rows",
                        "filtered_pairlist_rows",
                        "filtered_sdf_files",
                        "filtered_log_files",
                        "filtered_receptors",
                    ):
                        if key in existing:
                            summary_lines.append(f"{key}={existing[key]}")
                summary_lines.extend(
                    [
                        f"filtered_score_rows={len(self.scores_df)}",
                        f"selected_tags={len(self.selected_query_tags)}",
                    ]
                )
                summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
                self.logger.info(
                    f"🔎 Applied complex query to scores: {original_score_count} -> {len(self.scores_df)} rows"
                )
            self.scores_df = self._exclude_problematic_ligands(scores_csv, self.scores_df)
            if self.scores_df.empty:
                self.logger.error(
                    "❌ No rows left in all_scores.csv after problematic-ligand filtering. "
                    "Disable exclusion or relax thresholds."
                )
                return False
            self.logger.info(f"✅ Generated all_scores.csv with {len(self.scores_df)} scores")
            return True
        else:
            self.logger.error("❌ Failed to generate all_scores.csv")
            return False

    @staticmethod
    def _pairlist_tag_series(frame: pd.DataFrame) -> pd.Series:
        if frame.empty:
            return pd.Series(dtype="object")
        if {"receptor", "site_id", "ligand"}.issubset(frame.columns):
            return frame.apply(
                lambda row: f"{row['receptor']}_{row['site_id']}_{row['ligand']}",
                axis=1,
            ).astype(str)
        if "tag" in frame.columns:
            return frame["tag"].astype(str)
        return pd.Series(dtype="object")

    def _exclude_problematic_ligands(self, scores_csv: Path, frame: pd.DataFrame) -> pd.DataFrame:
        if not self.exclude_problematic_ligands or frame.empty:
            return frame
        if "tag" not in frame.columns:
            self.logger.warning("⚠️  Problematic-ligand exclusion requested but scores have no 'tag' column")
            return frame

        work = frame.copy()
        work["tag"] = work["tag"].astype(str)
        work["vina_affinity"] = pd.to_numeric(work.get("vina_affinity"), errors="coerce")
        best_affinity = work.groupby("tag", dropna=False)["vina_affinity"].min()
        pose_counts = work.groupby("tag", dropna=False).size()

        score_tags = set(work["tag"].dropna().astype(str))
        pair_tags = set(self._pairlist_tag_series(self.pairlist_df).dropna().astype(str).tolist())
        all_tags = score_tags.union(pair_tags)

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
            return work

        keep_tags = all_tags.difference(excluded_tags)
        filtered = work.loc[work["tag"].isin(keep_tags)].copy()
        self.sdf_files = filter_paths_by_tags(self.sdf_files, keep_tags)
        self.log_files = filter_paths_by_tags(self.log_files, keep_tags)
        if not self.pairlist_df.empty:
            pair_tag_series = self._pairlist_tag_series(self.pairlist_df)
            self.pairlist_df = self.pairlist_df.loc[pair_tag_series.isin(keep_tags)].copy()
            if self.analysis_pairlist_file is not None:
                self.pairlist_df.to_csv(self.analysis_pairlist_file, index=False)

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
        exclusion_file = self.output_dir / "excluded_problematic_ligands.csv"
        exclusion_df.to_csv(exclusion_file, index=False)
        filtered.to_csv(scores_csv, index=False)

        if self.selected_query_tags:
            self.selected_query_tags = self.selected_query_tags.intersection(keep_tags)
        self.logger.info(
            f"🚫 Excluded {len(excluded_tags)} problematic ligands "
            f"(report: {exclusion_file})"
        )
        return filtered
    
    def _match_poses_to_receptors(self) -> bool:
        """Match SDF poses to receptors using pairlist or filename patterns."""
        self.logger.info("🔗 Matching poses to receptors...")

        self.complexes = match_poses_to_receptors(
            self.sdf_files,
            self.receptor_files,
            self.pairlist_df
        )

        # Attach resolved display names to each complex for downstream labels/manifests.
        for complex_info in self.complexes:
            receptor_name = str(complex_info.get("receptor_name") or complex_info.get("receptor_file") or "")
            display_name = resolve_protein_display_name(receptor_name, self.protein_name_map)
            complex_info["protein_display_name"] = display_name
            complex_info["pdb_code"] = extract_pdb_code(receptor_name)
            complex_info["protein_label"] = format_protein_label(receptor_name, self.protein_name_map)
        
        self.logger.info(f"✅ Matched {len(self.complexes)} complexes")
        return len(self.complexes) > 0

    @staticmethod
    def _infer_element_from_pdbqt(atom_name: str, atom_type: str) -> str:
        """
        Infer a valid element symbol from PDBQT atom name/type.
        """
        name = "".join(ch for ch in str(atom_name).strip() if ch.isalpha())
        adt = str(atom_type or "").strip()
        adt_upper = adt.upper()
        adt_lower = adt.lower()

        if adt_lower in {"cl", "br", "zn", "mg", "mn", "fe", "cu", "ni", "co", "na", "cd", "hg"}:
            return adt_lower.capitalize()
        if adt_lower == "ca":
            # Ambiguous in AutoDock types: "CA" can be aromatic carbon.
            if adt == "Ca":
                return "Ca"
            return "C"
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

        # Fallback from atom name.
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
    def _pdbqt_atom_to_pdb(cls, line: str, chain_id: str = "A") -> Optional[str]:
        """Convert a single ATOM/HETATM PDBQT record to a PDB-style fixed-width record."""
        if not line.startswith(("ATOM", "HETATM")):
            return None
        raw = line.rstrip("\n")
        if len(raw) < 80:
            raw = raw.ljust(80)
        # Keep ATOM semantics for receptor lines to avoid ligand misclassification.
        raw = "ATOM  " + raw[6:]
        raw = raw[:21] + chain_id + raw[22:]
        atom_name = raw[12:16]
        atom_type = raw.split()[-1] if raw.split() else ""
        element = cls._infer_element_from_pdbqt(atom_name, atom_type)
        raw = raw[:76] + f"{element:>2}" + raw[78:]
        return raw[:80]

    @staticmethod
    def _parse_atom_serial(pdb_line: str) -> Optional[int]:
        """Parse atom serial from a fixed-width PDB line."""
        try:
            token = str(pdb_line)[6:11].strip()
            return int(token) if token else None
        except Exception:
            return None

    @staticmethod
    def _infer_element_from_atom_name(atom_name: str) -> str:
        """Infer element from an atom name when element column is missing."""
        token = "".join(ch for ch in str(atom_name).strip() if ch.isalpha()).upper()
        if not token:
            return ""
        if token.startswith("CL"):
            return "Cl"
        if token.startswith("BR"):
            return "Br"
        if token.startswith("ZN"):
            return "Zn"
        if token.startswith("MG"):
            return "Mg"
        return token[0]

    def _read_receptor_lines_from_pdbqt(self, receptor_file: Path) -> List[str]:
        """
        Read receptor directly from PDBQT text to avoid OpenBabel aromatic-kekulization warnings.
        """
        receptor_lines: List[str] = []
        with open(receptor_file, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                converted = self._pdbqt_atom_to_pdb(line, chain_id="A")
                if converted:
                    receptor_lines.append(converted)
        return receptor_lines

    @staticmethod
    def _normalize_ligand_resname(raw_name: str) -> str:
        """
        Convert arbitrary ligand token into a 3-character PDB residue name.
        """
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
            candidates.extend(
                [
                    token[:3],
                    f"{token[0]}{token[-2:]}",
                    f"{token[:2]}{token[-1]}",
                    token[-3:],
                ]
            )
        elif len(token) == 2:
            candidates.extend([token + "X", f"{token[0]}X{token[1]}", f"X{token}"])
        else:
            candidates.extend([token + "XX", f"X{token}X", f"XX{token}"])

        normalized = []
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

    def _derive_ligand_identity(self, complex_info: Dict[str, object]) -> Tuple[str, str, str]:
        """
        Derive ligand (resname, chain, resnum) from pairlist-style names when available.
        """
        ligand_name = str(complex_info.get("ligand_name") or "")
        complex_name = str(complex_info.get("complex_name") or "")
        source = ligand_name or complex_name
        if not source:
            pose_file = complex_info.get("pose_file")
            if pose_file:
                source = Path(str(pose_file)).stem

        # Prefer explicit receptor ligand notation:
        # *_ligand_AQ4_A_999.pdbqt -> AQ4, A, 999
        ligand_match = re.search(r"_ligand_([A-Za-z0-9]{1,6})_([A-Za-z])_(\d+)", source)
        if ligand_match:
            resname = self._normalize_ligand_resname(ligand_match.group(1))
            chain = ligand_match.group(2).upper()
            resnum = ligand_match.group(3)
            return resname, chain, resnum

        stem = Path(source).stem
        token = stem.split("_")[-1] if stem else "LIG"
        if token.lower() in {"pdbqt", "poses", "top", "out"} and "_" in stem:
            for part in reversed(stem.split("_")):
                if part.lower() not in {"pdbqt", "poses", "top", "out"}:
                    token = part
                    break
        resname = self._normalize_ligand_resname(token)
        return resname, "B", "1"
    
    def _create_complexes(self) -> bool:
        """Create complex PDB files (receptor + ligand)."""
        self.logger.info("🧬 Creating complex PDB files...")
        
        complexes_dir = self.output_dir / "complexes"
        complexes_dir.mkdir(exist_ok=True)
        
        written = 0
        for complex_info in self.complexes:
            receptor_file = complex_info.get('receptor_file')
            pose_file = complex_info.get('pose_file')
            complex_name = complex_info.get('complex_name', 'unknown')
            
            if not receptor_file or not pose_file:
                self.logger.warning(f"  ⚠️  Missing files for {complex_name}")
                continue
            
            # Create complex PDB using pose_extractor logic
            try:
                output_pdb = complexes_dir / f"{complex_name}.pdb"
                
                # Use direct PDBQT parsing for receptor first; fallback to OpenBabel only if needed.
                from openbabel import pybel
                
                # Read receptor PDBQT
                receptor_lines: List[str] = []
                if Path(receptor_file).exists():
                    try:
                        receptor_lines = self._read_receptor_lines_from_pdbqt(Path(receptor_file))
                        if not receptor_lines:
                            self.logger.warning(
                                f"  ⚠️  No receptor ATOM/HETATM records found in {receptor_file}"
                            )
                            continue
                    except Exception as e:
                        self.logger.warning(f"  ⚠️  Could not read receptor {receptor_file}: {e}")
                        continue
                
                # Read ligand SDF
                ligand_lines = []
                ligand_conect_lines = []
                try:
                    lig_resname, lig_chain, lig_resnum = self._derive_ligand_identity(complex_info)
                    max_serial = 0
                    for rec in receptor_lines:
                        serial = self._parse_atom_serial(rec)
                        if serial is not None:
                            max_serial = max(max_serial, serial)
                    ligand_serial = max_serial
                    ligand_mol = next(pybel.readfile("sdf", str(pose_file)))
                    ligand_pdb = ligand_mol.write("pdb")
                    serial_map: Dict[int, int] = {}
                    raw_conect_lines: List[str] = []
                    for line in ligand_pdb.split('\n'):
                        if line.startswith('ATOM') or line.startswith('HETATM'):
                            line = line.ljust(80)
                            old_serial = self._parse_atom_serial(line)
                            ligand_serial += 1
                            new_line = f"HETATM{line[6:]}"
                            new_line = new_line[:6] + f"{ligand_serial:5d}" + new_line[11:]
                            new_line = new_line[:17] + f"{lig_resname:>3}" + new_line[20:]
                            new_line = new_line[:21] + lig_chain + new_line[22:]
                            new_line = new_line[:22] + str(lig_resnum)[-4:].rjust(4) + new_line[26:]
                            atom_name = new_line[12:16]
                            element = new_line[76:78].strip() or self._infer_element_from_atom_name(atom_name)
                            new_line = new_line[:76] + f"{element:>2}" + new_line[78:]
                            ligand_lines.append(new_line[:80])
                            if old_serial is not None:
                                serial_map[old_serial] = ligand_serial
                        elif line.startswith("CONECT"):
                            raw_conect_lines.append(line)

                    # Preserve ligand connectivity for tools that rely on bond topology.
                    for conect_line in raw_conect_lines:
                        parts = conect_line.split()
                        if len(parts) < 3:
                            continue
                        try:
                            src_old = int(parts[1])
                        except ValueError:
                            continue
                        src_new = serial_map.get(src_old)
                        if src_new is None:
                            continue

                        remapped_targets: List[int] = []
                        seen_targets = set()
                        for token in parts[2:]:
                            try:
                                dst_old = int(token)
                            except ValueError:
                                continue
                            dst_new = serial_map.get(dst_old)
                            if dst_new is None or dst_new == src_new:
                                continue
                            if dst_new in seen_targets:
                                continue
                            seen_targets.add(dst_new)
                            remapped_targets.append(dst_new)

                        if not remapped_targets:
                            continue

                        # PDB CONECT records carry up to 4 neighbors per line.
                        for i in range(0, len(remapped_targets), 4):
                            chunk = remapped_targets[i:i + 4]
                            record = "CONECT" + f"{src_new:5d}" + "".join(f"{val:5d}" for val in chunk)
                            ligand_conect_lines.append(record)
                except Exception as e:
                    self.logger.warning(f"  ⚠️  Could not read ligand {pose_file}: {e}")
                    continue

                if not ligand_lines:
                    self.logger.warning(f"  ⚠️  No ligand atoms found in {pose_file}")
                    continue
                
                # Combine receptor and ligand
                all_lines = receptor_lines + ["TER"] + ligand_lines + ["TER"] + ligand_conect_lines + ["END"]
                
                with open(output_pdb, 'w') as f:
                    f.write('\n'.join(all_lines))
                
                written += 1
                self.logger.debug(f"  ✅ Created complex: {complex_name}")
                
            except ImportError:
                self.logger.error("  ❌ OpenBabel not available for complex creation")
                return False
            except Exception as e:
                self.logger.warning(f"  ⚠️  Error creating complex {complex_name}: {e}")
                continue
        
        self.logger.info(f"✅ Created {written} complex PDB files")
        return written > 0

    @staticmethod
    def _chain_set_from_file(path: Path) -> set:
        chains: set = set()
        try:
            with path.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    if not line.startswith(("ATOM", "HETATM")):
                        continue
                    chain = line[21:22]
                    chains.add(chain if chain.strip() else "_BLANK_")
        except Exception:
            return set()
        return chains

    @staticmethod
    def _residue_spread_metrics(path: Path) -> Dict[str, object]:
        residues: Dict[Tuple[str, str, str], List[Tuple[float, float, float]]] = {}
        atom_rows = 0
        try:
            with path.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    if not line.startswith(("ATOM", "HETATM")):
                        continue
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except Exception:
                        continue
                    atom_rows += 1
                    key = (line[21:22], line[22:26].strip(), line[17:20].strip())
                    residues.setdefault(key, []).append((x, y, z))
        except Exception:
            return {
                "atom_rows": 0,
                "residue_keys": 0,
                "residue_spread_gt_12_count": 0,
                "residue_spread_gt_20_count": 0,
                "max_residue_spread_angstrom": 0.0,
            }

        def _point_distance(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
            dx = a[0] - b[0]
            dy = a[1] - b[1]
            dz = a[2] - b[2]
            return float((dx * dx + dy * dy + dz * dz) ** 0.5)

        def _spread(points: List[Tuple[float, float, float]]) -> float:
            if len(points) < 2:
                return 0.0
            max_d = 0.0
            for i, p in enumerate(points):
                for q in points[i + 1:]:
                    d = _point_distance(p, q)
                    if d > max_d:
                        max_d = d
            return max_d

        spread_values = [_spread(points) for points in residues.values()]
        spread_gt_12 = sum(1 for value in spread_values if value > 12.0)
        spread_gt_20 = sum(1 for value in spread_values if value > 20.0)
        max_spread = max(spread_values) if spread_values else 0.0
        return {
            "atom_rows": atom_rows,
            "residue_keys": len(residues),
            "residue_spread_gt_12_count": spread_gt_12,
            "residue_spread_gt_20_count": spread_gt_20,
            "max_residue_spread_angstrom": round(float(max_spread), 3),
        }

    def _candidate_cleaned_receptor_file(self, receptor_file: Path) -> Optional[Path]:
        stem = receptor_file.stem
        candidates = [
            receptor_file.with_suffix(".pdb"),
            self.receptors_folder.parent / "2-Raw_Protien" / f"{stem}_cleaned.pdb",
            self.receptors_folder.parent.parent / "2-Raw_Protien" / f"{stem}_cleaned.pdb",
            self.receptors_folder.parent.parent.parent / "2-Raw_Protien" / f"{stem}_cleaned.pdb",
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return None

    def _repair_receptor_with_openbabel(self, receptor_file: Path, cleaned_pdb: Path) -> Tuple[bool, str]:
        cmd = [
            "obabel",
            "-ipdb",
            str(cleaned_pdb),
            "-opdbqt",
            "-O",
            str(receptor_file),
            "-xr",
            "-xh",
            "-xn",
        ]
        try:
            completed = subprocess.run(cmd, capture_output=True, text=True, check=False)
        except FileNotFoundError:
            return False, "obabel_not_found"
        except Exception as exc:
            return False, f"repair_exec_error:{exc}"
        if int(completed.returncode) != 0:
            reason = (completed.stderr or completed.stdout or "").strip()
            return False, f"repair_failed:{reason[:300]}"
        chains = self._chain_set_from_file(receptor_file)
        if chains and chains != {"_BLANK_"}:
            return True, "repaired_with_openbabel"
        return False, "repair_produced_blank_chains"

    def _receptor_chain_integrity_qc_and_repair(self) -> bool:
        """Run receptor chain-integrity QC and optional OpenBabel-based repair."""
        self.logger.info("🧪 Running receptor chain-integrity QC...")
        report_rows: List[Dict[str, object]] = []
        repairs: List[Dict[str, object]] = []
        auto_repair = bool(self.interaction_analytics_config.get("chain_integrity_auto_repair", True))
        fatal_count = 0

        for receptor_file in self.receptor_files:
            receptor_path = Path(receptor_file)
            if not receptor_path.exists():
                continue
            cleaned_pdb = self._candidate_cleaned_receptor_file(receptor_path)
            receptor_chains = self._chain_set_from_file(receptor_path)
            cleaned_chains = self._chain_set_from_file(cleaned_pdb) if cleaned_pdb else set()
            spread = self._residue_spread_metrics(receptor_path)

            has_blank_only = receptor_chains == {"_BLANK_"}
            multichain_source = len(set(cleaned_chains) - {"_BLANK_"}) > 1
            merged_residue_risk = bool(has_blank_only and int(spread.get("residue_spread_gt_20_count", 0) or 0) > 0)

            if merged_residue_risk and multichain_source:
                status = "FAIL_MERGED_RESIDUES"
            elif merged_residue_risk:
                status = "WARN_MERGED_RESIDUES"
            elif has_blank_only:
                status = "WARN_BLANK_CHAIN_ONLY"
            else:
                status = "PASS"

            repaired = False
            repair_reason = ""
            if auto_repair and status in {"FAIL_MERGED_RESIDUES", "WARN_MERGED_RESIDUES", "WARN_BLANK_CHAIN_ONLY"} and cleaned_pdb:
                repaired, repair_reason = self._repair_receptor_with_openbabel(receptor_path, cleaned_pdb)
                repairs.append(
                    {
                        "receptor_pdbqt": str(receptor_path),
                        "cleaned_pdb": str(cleaned_pdb),
                        "repair_attempted": True,
                        "repair_success": bool(repaired),
                        "repair_reason": repair_reason,
                    }
                )
                if repaired:
                    receptor_chains = self._chain_set_from_file(receptor_path)
                    spread = self._residue_spread_metrics(receptor_path)
                    has_blank_only = receptor_chains == {"_BLANK_"}
                    merged_residue_risk = bool(has_blank_only and int(spread.get("residue_spread_gt_20_count", 0) or 0) > 0)
                    if has_blank_only and merged_residue_risk:
                        status = "FAIL_AFTER_REPAIR"
                    elif has_blank_only:
                        status = "WARN_BLANK_CHAIN_ONLY_AFTER_REPAIR"
                    else:
                        status = "PASS_REPAIRED"
            elif status != "PASS":
                repairs.append(
                    {
                        "receptor_pdbqt": str(receptor_path),
                        "cleaned_pdb": str(cleaned_pdb) if cleaned_pdb else "",
                        "repair_attempted": False,
                        "repair_success": False,
                        "repair_reason": "cleaned_pdb_missing" if cleaned_pdb is None else "auto_repair_disabled",
                    }
                )

            if str(status).startswith("FAIL"):
                fatal_count += 1

            report_rows.append(
                {
                    "receptor_pdbqt": str(receptor_path),
                    "cleaned_pdb": str(cleaned_pdb) if cleaned_pdb else "",
                    "atom_rows": int(spread.get("atom_rows", 0) or 0),
                    "residue_keys": int(spread.get("residue_keys", 0) or 0),
                    "receptor_chains": ",".join(sorted(receptor_chains)),
                    "cleaned_chains": ",".join(sorted(cleaned_chains)),
                    "residue_spread_gt_12_count": int(spread.get("residue_spread_gt_12_count", 0) or 0),
                    "residue_spread_gt_20_count": int(spread.get("residue_spread_gt_20_count", 0) or 0),
                    "max_residue_spread_angstrom": float(spread.get("max_residue_spread_angstrom", 0.0) or 0.0),
                    "status": status,
                    "auto_repair_enabled": bool(auto_repair),
                }
            )

        qc_dir = self.raw_data_root / "quality_control"
        qc_dir.mkdir(parents=True, exist_ok=True)
        report_csv = qc_dir / "receptor_chain_integrity_report.csv"
        repairs_csv = qc_dir / "receptor_chain_repair_actions.csv"
        pd.DataFrame(report_rows).to_csv(report_csv, index=False)
        pd.DataFrame(repairs).to_csv(repairs_csv, index=False)
        self.results["receptor_chain_integrity"] = {
            "report_csv": str(report_csv),
            "repairs_csv": str(repairs_csv),
            "rows": int(len(report_rows)),
            "fatal_count": int(fatal_count),
        }
        self.logger.info(f"✅ Receptor chain-integrity report written: {report_csv}")
        if fatal_count > 0:
            self.logger.error(f"❌ Receptor chain-integrity QC found {fatal_count} unrepaired fatal receptors")
            return False
        return True
    
    def _analyze_binding_affinity(self) -> bool:
        """
        Analyze binding affinities using hierarchical analysis.
        
        Hierarchy:
        1. Best pose per ligand-protein pair
        2. Best ligand per protein
        3. Cross-protein comparison (same ligands)
        4. Comparative (redocking) analysis
        """
        self.logger.info("📈 Analyzing binding affinities (hierarchical)...")
        
        if self.scores_df is None or self.scores_df.empty:
            self.logger.error("❌ No scores data available")
            return False
        
        analysis_pairlist = self.analysis_pairlist_file
        if not analysis_pairlist or not Path(analysis_pairlist).exists():
            self.logger.warning("⚠️  No pairlist.csv - falling back to basic analysis")
            return self._analyze_binding_affinity_basic()
        
        try:
            # Use hierarchical analyzer with pairlist mapping
            scores_csv = self.output_dir / "all_scores.csv"
            analysis_dir = self.output_dir / "analysis"
            analysis_dir.mkdir(exist_ok=True)
            
            # Run hierarchical analysis
            analyzer = HierarchicalDockingAnalyzer(
                str(scores_csv),
                str(analysis_pairlist),
                str(analysis_dir),
                protein_name_map=self.protein_name_map
            )
            
            hierarchical_results = analyzer.analyze_all()
            report = analyzer.generate_report()

            if "best_poses" in hierarchical_results:
                hierarchical_results["best_poses"] = self._apply_display_names(
                    hierarchical_results["best_poses"],
                    protein_col="protein"
                )
            
            # Print report
            print("\n" + report)
            
            # Store results
            self.results['hierarchical_analysis'] = hierarchical_results
            self.results['proteins'] = analyzer.proteins
            self.results['ligands'] = analyzer.ligands
            
            # Also organize best poses by affinity
            self._organize_poses_by_affinity(hierarchical_results['best_poses'])
            
            self.logger.info("✅ Hierarchical analysis completed")
            return True
            
        except Exception as e:
            self.logger.error(f"❌ Error in hierarchical analysis: {e}", exc_info=True)
            return False
    
    def _analyze_binding_affinity_basic(self) -> bool:
        """Basic binding affinity analysis without pairlist."""
        try:
            scores_csv = self.output_dir / "all_scores.csv"
            
            analysis_results = analyze_binding_affinities(
                str(scores_csv),
                comparative_benchmark="*",
                top_count=20
            )
            
            # Save results
            reports_dir = self.output_dir / "reports"
            reports_dir.mkdir(exist_ok=True)
            
            analysis_results['best_poses'].to_csv(
                reports_dir / "best_poses.csv", index=False
            )
            analysis_results['summary_stats'].to_csv(
                reports_dir / "summary_stats.csv", index=False
            )

            analysis_results['best_poses'] = self._apply_display_names(
                analysis_results.get('best_poses', pd.DataFrame()),
                protein_col="protein"
            )

            # Keep categorized best_poses tree available even in basic mode so
            # downstream PoseView stage can use the same directory contract.
            if isinstance(analysis_results.get('best_poses'), pd.DataFrame):
                self._organize_poses_by_affinity(analysis_results['best_poses'])
            
            self.results['affinity_analysis'] = analysis_results
            return True
            
        except Exception as e:
            self.logger.error(f"❌ Error in basic analysis: {e}", exc_info=True)
            return False

    @staticmethod
    def _is_reference_candidate_row(row: pd.Series) -> bool:
        """
        Detect whether a row should be treated as a reference/control entry.
        """
        site_id = str(row.get("site_id") or "").strip().lower()
        if site_id in {
            "reference",
            "comparative",
            "compartive",
            "lapi",
            "control",
            "benchmark",
            "known",
            "native",
            "redocking",
        }:
            return True

        text = " ".join(
            str(row.get(key) or "")
            for key in ("tag", "ligand", "ligand_name")
        ).lower()
        tokens = (
            "reference",
            "co-crystal",
            "cocrystal",
            "native",
            "control",
            "benchmark",
            "redocking",
        )
        return any(token in text for token in tokens)

    @staticmethod
    def _serialize_targets(values: List[str]) -> str:
        """Serialize target labels into deterministic semicolon-separated text."""
        cleaned = sorted({str(v).strip() for v in values if str(v).strip()})
        return ";".join(cleaned)

    def _analyze_polypharmacology(self) -> bool:
        """
        Ligand-centric polypharmacology analysis with reference-aware comparison.

        Outputs:
        - per ligand-protein pair table with delta-vs-reference
        - ligand summary (multi-target potential score + target lists)
        - affinity and delta matrices
        - optional visual summaries
        """
        if not self.polypharm_config.get("enabled", True):
            self.logger.info("⏭️  Polypharmacology analysis disabled")
            return True

        self.logger.info("🧭 Running polypharmacology analysis (reference-aware)...")

        best_poses = self._get_best_poses_for_visualization()
        if best_poses is None or best_poses.empty:
            self.logger.warning("⚠️  No best-pose data available for polypharmacology analysis")
            return True

        if "protein_label" not in best_poses.columns:
            best_poses = self._apply_display_names(best_poses, protein_col="protein")
        if "protein_label" not in best_poses.columns or "ligand" not in best_poses.columns or "vina_affinity" not in best_poses.columns:
            self.logger.warning("⚠️  Missing required columns for polypharmacology analysis")
            return True

        try:
            strong_threshold = float(self.polypharm_config.get("strong_threshold", -7.0))
            moderate_threshold = float(self.polypharm_config.get("moderate_threshold", -5.0))
            near_reference_delta = float(self.polypharm_config.get("near_reference_delta", 1.0))
            min_targets_moderate = int(self.polypharm_config.get("min_targets_moderate", 2))
            min_targets_high = int(self.polypharm_config.get("min_targets_high", 3))
            exclude_reference_ligands = bool(self.polypharm_config.get("exclude_reference_ligands", True))

            df = best_poses.copy()
            df["ligand"] = df["ligand"].astype(str)
            df["protein_label"] = df["protein_label"].astype(str)
            df["vina_affinity"] = pd.to_numeric(df["vina_affinity"], errors="coerce")
            df = df[np.isfinite(df["vina_affinity"])].copy()
            if df.empty:
                self.logger.warning("⚠️  No finite affinity values available for polypharmacology analysis")
                return True

            if "engine" in df and df["engine"].nunique() > 1:
                unavailable_dir = self.output_dir / "analysis" / "polypharmacology"
                unavailable_dir.mkdir(parents=True, exist_ok=True)
                status = {"status":"not_evaluable", "reason":"cross_engine_raw_score_comparison_requires_calibration"}
                (unavailable_dir / "status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
                self.results["polypharmacology_analysis"] = status
                return True
            # Ensure one best row per ligand-protein pair.
            idx = df.groupby(["ligand", "protein_label"])["vina_affinity"].idxmin()
            pair_df = df.loc[idx].copy().sort_values(["ligand", "protein_label"])

            pair_df["interpretation"] = "exploratory_uncalibrated_docking_scores_not_binding_or_selectivity"
            pair_df["reference_validation"] = "not_established_by_this_analysis"
            pair_df["is_reference_candidate"] = False
            pair_df["strong_binder"] = pair_df["vina_affinity"] <= strong_threshold
            pair_df["moderate_or_better_binder"] = pair_df["vina_affinity"] <= moderate_threshold
            pair_df["unfavorable_binder"] = pair_df["vina_affinity"] > 0.0
            pair_df["binding_category"] = pair_df["vina_affinity"].apply(self._affinity_category)
            pair_df["binding_category_label"] = pair_df["binding_category"].apply(self._format_category_label)

            reference_df = pair_df[pair_df["is_reference_candidate"]].copy()
            if not reference_df.empty:
                reference_by_protein = (
                    reference_df.groupby("protein_label", as_index=False)["vina_affinity"]
                    .min()
                    .rename(columns={"vina_affinity": "reference_affinity"})
                )
                pair_df = pair_df.merge(reference_by_protein, on="protein_label", how="left")
            else:
                reference_by_protein = pd.DataFrame(columns=["protein_label", "reference_affinity"])
                pair_df["reference_affinity"] = np.nan

            pair_df["delta_vs_reference"] = pair_df["vina_affinity"] - pair_df["reference_affinity"]

            within_label = f"within_{near_reference_delta:.1f}_kcal_of_reference"

            def _reference_status(row: pd.Series) -> str:
                ref = row.get("reference_affinity")
                delta = row.get("delta_vs_reference")
                if ref is None or not np.isfinite(ref):
                    return "no_reference_for_protein"
                if delta is None or not np.isfinite(delta):
                    return "no_reference_for_protein"
                if float(delta) <= 0.0:
                    return "better_or_equal_than_reference"
                if float(delta) <= near_reference_delta:
                    return within_label
                return "worse_than_reference"

            pair_df["reference_comparison"] = pair_df.apply(_reference_status, axis=1)
            pair_df["better_than_reference"] = pair_df["reference_comparison"] == "better_or_equal_than_reference"
            pair_df["near_reference"] = pair_df["reference_comparison"].isin(
                ["better_or_equal_than_reference", within_label]
            )
            pair_df["polypharm_hit"] = np.where(
                pair_df["reference_affinity"].notna(),
                pair_df["strong_binder"] & pair_df["near_reference"],
                pair_df["strong_binder"],
            )

            ligand_reference_flag = (
                pair_df.groupby("ligand")["is_reference_candidate"]
                .any()
                .rename("ligand_is_reference")
            )
            pair_df = pair_df.merge(ligand_reference_flag, on="ligand", how="left")

            candidate_df = pair_df.copy()
            if exclude_reference_ligands:
                candidate_df = candidate_df[~candidate_df["ligand_is_reference"]].copy()
            if candidate_df.empty:
                # Keep analysis alive when all entries are references/controls.
                candidate_df = pair_df.copy()

            total_proteins = int(pair_df["protein_label"].nunique())
            ligand_rows = []
            for ligand, ligand_df in candidate_df.groupby("ligand"):
                target_list = ligand_df["protein_label"].astype(str).tolist()
                strong_targets = ligand_df.loc[ligand_df["strong_binder"], "protein_label"].astype(str).tolist()
                moderate_targets = ligand_df.loc[ligand_df["moderate_or_better_binder"], "protein_label"].astype(str).tolist()
                poly_targets = ligand_df.loc[ligand_df["polypharm_hit"], "protein_label"].astype(str).tolist()
                better_targets = ligand_df.loc[ligand_df["better_than_reference"], "protein_label"].astype(str).tolist()
                near_targets = ligand_df.loc[ligand_df["near_reference"], "protein_label"].astype(str).tolist()
                unfavorable_targets = ligand_df.loc[ligand_df["unfavorable_binder"], "protein_label"].astype(str).tolist()
                reference_targets = ligand_df.loc[ligand_df["reference_affinity"].notna(), "protein_label"].astype(str).tolist()

                strong_count = len(set(strong_targets))
                moderate_count = len(set(moderate_targets))
                poly_count = len(set(poly_targets))
                better_count = len(set(better_targets))
                near_count = len(set(near_targets))
                unfavorable_count = len(set(unfavorable_targets))
                reference_count = len(set(reference_targets))

                mean_aff = float(ligand_df["vina_affinity"].mean())
                std_aff = float(ligand_df["vina_affinity"].std(ddof=0)) if len(ligand_df) > 1 else 0.0

                # Composite score balancing potency, cross-target reach, and reference comparison.
                poly_score = (
                    (2.0 * poly_count)
                    + (1.5 * better_count)
                    + (0.75 * near_count)
                    + (0.25 * moderate_count)
                    - (0.75 * unfavorable_count)
                    - (0.20 * std_aff)
                    + max(0.0, (-mean_aff - 5.0))
                )

                if poly_count >= min_targets_high:
                    level = "Exploratory docking hits on at least three targets"
                elif poly_count >= min_targets_moderate:
                    level = "Exploratory docking hits on at least two targets"
                elif strong_count >= min_targets_moderate:
                    level = "Exploratory multi-target score profile"
                elif strong_count >= 1:
                    level = "Exploratory single-target score profile"
                else:
                    level = "No docking hits under exploratory cutoffs"

                ligand_rows.append(
                    {
                        "ligand": ligand,
                        "polypharm_score": round(poly_score, 4),
                        "interpretation": "uncalibrated_docking_heuristic_not_potency_or_selectivity",
                        "reference_validation": "not_established_by_this_analysis",
                        "polypharm_level": level,
                        "targets_tested_count": len(set(target_list)),
                        "total_proteins_in_project": total_proteins,
                        "strong_targets_count": strong_count,
                        "moderate_or_better_targets_count": moderate_count,
                        "polypharm_targets_count": poly_count,
                        "better_than_reference_targets_count": better_count,
                        "near_reference_targets_count": near_count,
                        "reference_available_targets_count": reference_count,
                        "unfavorable_targets_count": unfavorable_count,
                        "mean_affinity": round(mean_aff, 4),
                        "min_affinity": round(float(ligand_df["vina_affinity"].min()), 4),
                        "max_affinity": round(float(ligand_df["vina_affinity"].max()), 4),
                        "std_affinity": round(std_aff, 4),
                        "polypharm_coverage": round(poly_count / total_proteins, 4) if total_proteins > 0 else np.nan,
                        "strong_target_proteins": self._serialize_targets(strong_targets),
                        "moderate_or_better_target_proteins": self._serialize_targets(moderate_targets),
                        "polypharm_target_proteins": self._serialize_targets(poly_targets),
                        "better_than_reference_proteins": self._serialize_targets(better_targets),
                        "near_reference_proteins": self._serialize_targets(near_targets),
                        "unfavorable_proteins": self._serialize_targets(unfavorable_targets),
                    }
                )

            summary_df = pd.DataFrame(ligand_rows)
            if not summary_df.empty:
                summary_df = summary_df.sort_values(
                    ["polypharm_score", "polypharm_targets_count", "strong_targets_count", "mean_affinity"],
                    ascending=[False, False, False, True],
                )

            analysis_dir = self.output_dir / "analysis"
            analysis_dir.mkdir(exist_ok=True)
            poly_dir = analysis_dir / "polypharmacology"
            poly_dir.mkdir(exist_ok=True)
            viz_dir = analysis_dir / "visualizations"
            viz_dir.mkdir(exist_ok=True)

            pair_file = poly_dir / "polypharmacology_pair_level.csv"
            pair_df.to_csv(pair_file, index=False)

            candidate_file = poly_dir / "polypharmacology_pair_level_non_reference.csv"
            candidate_df.to_csv(candidate_file, index=False)

            reference_file = poly_dir / "polypharmacology_reference_baseline_per_protein.csv"
            reference_by_protein.to_csv(reference_file, index=False)

            summary_file = poly_dir / "polypharmacology_ligand_summary.csv"
            summary_df.to_csv(summary_file, index=False)

            long_cols = [
                "ligand",
                "protein_label",
                "vina_affinity",
                "reference_affinity",
                "delta_vs_reference",
                "reference_comparison",
                "strong_binder",
                "polypharm_hit",
                "binding_category_label",
            ]
            comparison_long_df = candidate_df[[c for c in long_cols if c in candidate_df.columns]].copy()
            comparison_long_file = poly_dir / "polypharmacology_reference_comparison_long.csv"
            comparison_long_df.to_csv(comparison_long_file, index=False)

            hits_df = candidate_df[
                candidate_df["polypharm_hit"] | candidate_df["near_reference"] | candidate_df["better_than_reference"]
            ].copy()
            hits_file = poly_dir / "polypharmacology_reference_hits.csv"
            hits_df.to_csv(hits_file, index=False)

            affinity_matrix = candidate_df.pivot_table(
                index="ligand",
                columns="protein_label",
                values="vina_affinity",
                aggfunc="min",
            )
            if not affinity_matrix.empty:
                affinity_matrix = affinity_matrix.loc[affinity_matrix.mean(axis=1).sort_values().index]
            affinity_matrix_file = poly_dir / "polypharmacology_affinity_matrix.csv"
            affinity_matrix.to_csv(affinity_matrix_file)

            delta_matrix = candidate_df.pivot_table(
                index="ligand",
                columns="protein_label",
                values="delta_vs_reference",
                aggfunc="min",
            )
            delta_matrix_file = poly_dir / "polypharmacology_delta_vs_reference_matrix.csv"
            delta_matrix.to_csv(delta_matrix_file)

            # Visual outputs (optional).
            plotting_available = False
            if self.run_visualizations:
                try:
                    import matplotlib.pyplot as plt
                    import seaborn as sns
                    plotting_available = True
                except Exception:
                    plotting_available = False

            generated_figures: List[str] = []
            if plotting_available:
                if not affinity_matrix.empty:
                    fig, ax = plt.subplots(
                        figsize=(max(12, len(affinity_matrix.columns) * 0.85), max(8, len(affinity_matrix.index) * 0.4))
                    )
                    sns.heatmap(
                        affinity_matrix,
                        cmap="RdYlGn_r",
                        annot=True,
                        fmt=".2f",
                        linewidths=0.35,
                        cbar_kws={"label": "Best Affinity (kcal/mol)"},
                        ax=ax,
                    )
                    ax.set_title("Exploratory Docking Scores (Ligand x Protein; uncalibrated)")
                    ax.set_xlabel("Protein")
                    ax.set_ylabel("Ligand")
                    plt.tight_layout()
                    affinity_plot = viz_dir / "polypharmacology_affinity_heatmap.png"
                    plt.savefig(affinity_plot, dpi=300, bbox_inches="tight")
                    plt.close(fig)
                    generated_figures.append(str(affinity_plot))

                if not delta_matrix.empty and delta_matrix.notna().any().any():
                    fig, ax = plt.subplots(
                        figsize=(max(12, len(delta_matrix.columns) * 0.85), max(8, len(delta_matrix.index) * 0.4))
                    )
                    sns.heatmap(
                        delta_matrix,
                        cmap="RdYlGn_r",
                        center=0.0,
                        annot=True,
                        fmt=".2f",
                        linewidths=0.35,
                        cbar_kws={"label": "Delta vs Reference (kcal/mol)"},
                        ax=ax,
                    )
                    ax.set_title("Delta vs Reference by Ligand and Protein (negative is better)")
                    ax.set_xlabel("Protein")
                    ax.set_ylabel("Ligand")
                    plt.tight_layout()
                    delta_plot = viz_dir / "polypharmacology_reference_delta_heatmap.png"
                    plt.savefig(delta_plot, dpi=300, bbox_inches="tight")
                    plt.close(fig)
                    generated_figures.append(str(delta_plot))

                if not summary_df.empty:
                    top_n = int(self.polypharm_config.get("top_ligands_plot_count", 20))
                    top_df = summary_df.head(top_n).copy().sort_values("polypharm_score", ascending=True)
                    color_map = {
                        "High Polypharmacology Potential": "#2ca02c",
                        "Moderate Polypharmacology Potential": "#1f77b4",
                        "Potential Multi-Target Binder": "#ffbf00",
                        "Single-Target Dominant": "#ff7f0e",
                        "Low Polypharmacology Potential": "#7f7f7f",
                    }
                    colors = [color_map.get(v, "#7f7f7f") for v in top_df["polypharm_level"]]
                    fig, ax = plt.subplots(figsize=(12, max(7, len(top_df) * 0.45)))
                    ax.barh(top_df["ligand"], top_df["polypharm_score"], color=colors, edgecolor="black")
                    for idx_row, (_, row) in enumerate(top_df.iterrows()):
                        ax.text(
                            float(row["polypharm_score"]) + 0.05,
                            idx_row,
                            f"targets={int(row['polypharm_targets_count'])}, better_ref={int(row['better_than_reference_targets_count'])}",
                            va="center",
                            fontsize=8,
                        )
                    ax.set_xlabel("Uncalibrated docking heuristic")
                    ax.set_ylabel("Ligand")
                    ax.set_title("Exploratory Docking Heuristic (not binding or selectivity evidence)")
                    ax.grid(True, axis="x", alpha=0.25)
                    plt.tight_layout()
                    leaderboard_plot = viz_dir / "polypharmacology_leaderboard.png"
                    plt.savefig(leaderboard_plot, dpi=300, bbox_inches="tight")
                    plt.close(fig)
                    generated_figures.append(str(leaderboard_plot))

            self.results["polypharmacology_analysis"] = {
                "summary_file": str(summary_file),
                "pair_level_file": str(pair_file),
                "pair_level_non_reference_file": str(candidate_file),
                "reference_baseline_file": str(reference_file),
                "comparison_long_file": str(comparison_long_file),
                "reference_hits_file": str(hits_file),
                "affinity_matrix_file": str(affinity_matrix_file),
                "delta_matrix_file": str(delta_matrix_file),
                "generated_figures": generated_figures,
                "total_ligands_ranked": int(summary_df.shape[0]),
                "reference_proteins_count": int(reference_by_protein.shape[0]),
            }

            self.logger.info("✅ Polypharmacology analysis completed")
            self.logger.info(f"   Ligands ranked: {summary_df.shape[0]}")
            self.logger.info(f"   Proteins with reference baselines: {reference_by_protein.shape[0]}")
            if not summary_df.empty:
                top_row = summary_df.iloc[0]
                self.logger.info(
                    "   Top ligand: "
                    f"{top_row.get('ligand')} (score={float(top_row.get('polypharm_score', np.nan)):.2f}, "
                    f"polypharm targets={int(top_row.get('polypharm_targets_count', 0))})"
                )
            return True

        except Exception as exc:
            self.logger.warning(f"⚠️  Polypharmacology analysis failed: {exc}", exc_info=True)
            # Do not fail pipeline due to optional analysis block.
            return True

    def _apply_display_names(self, df: pd.DataFrame, protein_col: str = "protein") -> pd.DataFrame:
        """
        Add global protein/ligand naming columns using detected naming maps.
        """
        if df is None or df.empty:
            return df

        out_df = df.copy()
        source_series = None
        if protein_col in out_df.columns:
            source_series = out_df[protein_col].astype(str)
        elif "receptor" in out_df.columns:
            source_series = out_df["receptor"].astype(str)
        elif "tag" in out_df.columns:
            source_series = out_df["tag"].astype(str)

        if source_series is None:
            return out_df

        out_df["protein_display_name"] = source_series.apply(
            lambda value: resolve_protein_display_name(str(value), self.protein_name_map)
        )
        out_df["protein_label"] = source_series.apply(
            lambda value: format_protein_label(str(value), self.protein_name_map)
        )

        ligand_series = None
        if "ligand" in out_df.columns:
            ligand_series = out_df["ligand"].astype(str)
        elif "ligand_name" in out_df.columns:
            ligand_series = out_df["ligand_name"].astype(str)

        if ligand_series is not None:
            out_df["ligand_display_name"] = ligand_series.apply(
                lambda value: resolve_ligand_display_name(str(value), self.ligand_name_map)
            )
            out_df["ligand_label"] = out_df["ligand_display_name"]
        return out_df
    
    def _organize_poses_by_affinity(self, best_poses_df: pd.DataFrame):
        """Organize best poses by binding affinity strength."""
        best_poses_dir = self.output_dir / "best_poses"
        best_poses_dir.mkdir(exist_ok=True)
        
        strong_binders_dir = best_poses_dir / "strong_binders"
        moderate_binders_dir = best_poses_dir / "moderate_binders"
        weak_binders_dir = best_poses_dir / "weak_binders"
        
        for dir_path in [strong_binders_dir, moderate_binders_dir, weak_binders_dir]:
            dir_path.mkdir(exist_ok=True)
        
        # Strong: < -8.0, Moderate: -6.0 to -8.0, Weak: > -6.0
        strong_threshold = -8.0
        moderate_threshold = -6.0
        
        for _, row in best_poses_df.iterrows():
            tag = row.get('tag', '')
            affinity = row.get('vina_affinity', 0)
            
            # Try to find complex file
            complex_file = self.output_dir / "complexes" / f"{tag}.pdb"
            if not complex_file.exists():
                continue
            
            if affinity <= strong_threshold:
                target_dir = strong_binders_dir
            elif affinity <= moderate_threshold:
                target_dir = moderate_binders_dir
            else:
                target_dir = weak_binders_dir
            
            try:
                shutil.copy2(complex_file, target_dir / complex_file.name)
            except Exception:
                pass
        
        self.logger.info(f"✅ Organized poses by affinity")
    
    def _analyze_rmsd(self) -> bool:
        """
        Perform RMSD analysis aligned with docking hierarchy:
        1) all poses per protein-ligand combination
        2) per-protein aggregate from those combinations
        3) global aggregate across all proteins/ligands
        """
        self.logger.info("📐 Analyzing RMSD per protein-ligand combination...")

        if self.scores_df is None or self.scores_df.empty:
            self.logger.error("❌ No scores data for RMSD analysis")
            return False

        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolAlign
            from .pose_geometry import fixed_frame_rmsd
        except Exception:
            self.logger.error(
                "❌ RDKit not available for pose RMSD calculations; RMSD stage is mandatory"
            )
            return False

        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            plotting_available = True
        except Exception:
            plotting_available = False
        try:
            from .enhanced_rmsd_analyzer import (
                calculate_rmsd_matrix_from_pdbs,
                analyze_pose_clustering_enhanced,
                analyze_conformational_diversity_enhanced,
                create_rmsd_visualizations_enhanced,
            )
            enhanced_best_pose_rmsd = True
        except Exception as exc:
            enhanced_best_pose_rmsd = False
            self.logger.info(f"ℹ️  Enhanced best-pose RMSD helpers unavailable: {exc}")

        def _normalize_tag(tag: str) -> str:
            return self._normalize_tag(tag)

        def _pairwise_values(matrix: np.ndarray) -> np.ndarray:
            if matrix.size == 0 or matrix.ndim != 2:
                return np.array([], dtype=float)
            upper = matrix[np.triu_indices(matrix.shape[0], k=1)]
            return upper[np.isfinite(upper)]

        def _finite_float(value: object) -> Optional[float]:
            try:
                val = float(value)
            except (TypeError, ValueError):
                return None
            if not np.isfinite(val):
                return None
            return val

        def _summarize_values(values: np.ndarray) -> Dict[str, Optional[float]]:
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

        def _save_matrix_heatmap(matrix_df: pd.DataFrame, output_file: Path, title: str) -> None:
            if not plotting_available:
                return
            fig, ax = plt.subplots(figsize=(max(7, len(matrix_df.columns) * 0.7), max(6, len(matrix_df.index) * 0.7)))
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

        def _write_json(output_file: Path, payload: Dict[str, object]) -> None:
            with open(output_file, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2)

        def _scope_signature(export_df: pd.DataFrame, scope_name: str) -> str:
            signature_df = export_df.copy()
            for column in signature_df.columns:
                signature_df[column] = signature_df[column].map(
                    lambda value: "" if pd.isna(value) else str(value)
                )
            payload = json.dumps(
                {
                    "scope_name": scope_name,
                    "force_global_rmsd": self.force_global_rmsd,
                    "global_rmsd_defer_threshold": self.global_rmsd_defer_threshold,
                    "signature_csv": signature_df.to_csv(index=False),
                },
                sort_keys=True,
            )
            return hashlib.sha256(payload.encode("utf-8")).hexdigest()

        def _load_cached_scope_summary(
            scope_dir: Path,
            summary_basename: str,
            input_signature: str,
        ) -> Optional[Dict[str, object]]:
            if not self.resume_rmsd:
                return None
            state_file = scope_dir / f"{summary_basename}_scope_state.json"
            summary_json = scope_dir / f"{summary_basename}.json"
            if not state_file.exists() or not summary_json.exists():
                return None
            try:
                state_payload = json.loads(state_file.read_text(encoding="utf-8"))
                if str(state_payload.get("input_signature") or "") != input_signature:
                    return None
                summary = json.loads(summary_json.read_text(encoding="utf-8"))
                if not isinstance(summary, dict):
                    return None
                summary["resumed"] = True
                summary["state"] = str(summary.get("state") or state_payload.get("state") or "completed")
                return summary
            except Exception:
                return None

        def _persist_scope_state(
            scope_dir: Path,
            summary_basename: str,
            input_signature: str,
            state: str,
            summary: Dict[str, object],
        ) -> None:
            state_payload = {
                "input_signature": input_signature,
                "state": state,
                "summary_file": str(scope_dir / f"{summary_basename}.json"),
                "scope": str(summary.get("scope") or ""),
            }
            _write_json(scope_dir / f"{summary_basename}_scope_state.json", state_payload)

        def _write_scope_skip_summary(scope_dir: Path, scope_name: str, reason: str) -> Dict[str, object]:
            scope_dir.mkdir(parents=True, exist_ok=True)
            summary = {
                "scope": scope_name,
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
                "reason": reason,
            }
            pd.DataFrame([summary]).to_csv(scope_dir / "best_pose_rmsd_summary.csv", index=False)
            _write_json(scope_dir / "best_pose_rmsd_summary.json", summary)
            _persist_scope_state(scope_dir, "best_pose_rmsd_summary", "", "skipped", summary)
            return summary

        def _run_best_pose_scope(
            scope_dir: Path,
            pose_rows: pd.DataFrame,
            scope_name: str,
            heatmap_title: str,
        ) -> Dict[str, object]:
            scope_dir.mkdir(parents=True, exist_ok=True)
            export_df = pose_rows.copy()
            if "pdb_file" in export_df.columns:
                export_df["pdb_file"] = export_df["pdb_file"].astype(str)
            inputs_file = scope_dir / "best_pose_inputs.csv"
            export_df.to_csv(inputs_file, index=False)
            input_signature = _scope_signature(export_df, scope_name)
            cached_summary = _load_cached_scope_summary(scope_dir, "best_pose_rmsd_summary", input_signature)
            if cached_summary is not None:
                return cached_summary

            if export_df.empty:
                summary = {
                    "scope": scope_name,
                    "state": "skipped_empty",
                    "pose_count": 0,
                    "pairwise_count": 0,
                    "matrix_file": "",
                    "inputs_file": str(inputs_file),
                    "visualizations_dir": "",
                    "input_signature": input_signature,
                    "mean_pairwise_rmsd": None,
                    "median_pairwise_rmsd": None,
                    "std_pairwise_rmsd": None,
                    "min_pairwise_rmsd": None,
                    "max_pairwise_rmsd": None,
                }
                pd.DataFrame([summary]).to_csv(scope_dir / "best_pose_rmsd_summary.csv", index=False)
                _write_json(scope_dir / "best_pose_rmsd_summary.json", summary)
                _persist_scope_state(scope_dir, "best_pose_rmsd_summary", input_signature, "skipped_empty", summary)
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
                    "visualizations_dir": "",
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
                pd.DataFrame([summary]).to_csv(scope_dir / "best_pose_rmsd_summary.csv", index=False)
                _write_json(scope_dir / "best_pose_rmsd_summary.json", summary)
                _persist_scope_state(scope_dir, "best_pose_rmsd_summary", input_signature, "deferred", summary)
                return summary

            pdb_files = [Path(str(value)) for value in export_df["pdb_file"].tolist()]
            matrix = np.zeros((len(pdb_files), len(pdb_files)), dtype=float)
            if len(pdb_files) > 1:
                if enhanced_best_pose_rmsd:
                    matrix, _ = calculate_rmsd_matrix_from_pdbs(
                        pdb_files,
                        ligand_only=True,
                        num_workers=self.rmsd_workers,
                    )
                else:
                    matrix = np.full((len(pdb_files), len(pdb_files)), np.nan, dtype=float)
                    np.fill_diagonal(matrix, 0.0)

            labels = export_df["tag"].astype(str).tolist()
            matrix_df = pd.DataFrame(matrix, index=labels, columns=labels)
            matrix_file = scope_dir / "best_pose_rmsd_matrix.csv"
            matrix_df.to_csv(matrix_file)

            pair_values = _pairwise_values(matrix)
            summary = {
                "scope": scope_name,
                "state": "completed",
                "pose_count": int(len(export_df)),
                "pairwise_count": int(pair_values.size),
                "matrix_file": str(matrix_file),
                "inputs_file": str(inputs_file),
                "visualizations_dir": str(scope_dir / "visualizations"),
                "input_signature": input_signature,
                **_summarize_values(pair_values),
            }
            pd.DataFrame([summary]).to_csv(scope_dir / "best_pose_rmsd_summary.csv", index=False)
            _write_json(scope_dir / "best_pose_rmsd_summary.json", summary)
            _persist_scope_state(scope_dir, "best_pose_rmsd_summary", input_signature, "completed", summary)

            if enhanced_best_pose_rmsd:
                try:
                    vis_dir = scope_dir / "visualizations"
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
                        vis_dir,
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

                    overall_stats = diversity_results.get("overall_stats", {})
                    if isinstance(overall_stats, dict):
                        pd.DataFrame([overall_stats]).to_csv(
                            scope_dir / "diversity_overall_stats.csv",
                            index=False,
                        )
                except Exception as exc:
                    self.logger.warning(f"⚠️  Best-pose RMSD visuals unavailable for {scope_name}: {exc}")
                    if plotting_available:
                        _save_matrix_heatmap(matrix_df, scope_dir / "rmsd_heatmap.png", heatmap_title)
            elif plotting_available:
                _save_matrix_heatmap(matrix_df, scope_dir / "rmsd_heatmap.png", heatmap_title)

            return summary

        try:
            rmsd_dir = self.output_dir / "rmsd_analysis"
            rmsd_dir.mkdir(exist_ok=True)
            per_complex_dir = rmsd_dir / "per_complex_all_poses"
            per_protein_dir = rmsd_dir / "per_protein_best_poses"
            global_best_dir = rmsd_dir / "global_best_poses"
            per_complex_dir.mkdir(exist_ok=True)
            per_protein_dir.mkdir(exist_ok=True)
            global_best_dir.mkdir(exist_ok=True)

            # Map tag -> metadata and associated SDF file.
            tag_meta: Dict[str, Dict[str, str]] = {}
            for complex_info in self.complexes:
                tag = _normalize_tag(str(complex_info.get("complex_name") or ""))
                if not tag:
                    continue

                receptor_name = str(
                    complex_info.get("receptor_name")
                    or complex_info.get("receptor_file")
                    or ""
                )
                ligand_name = str(
                    complex_info.get("ligand_name")
                    or Path(str(complex_info.get("pose_file") or "")).stem
                    or "UnknownLigand"
                )
                protein_name = resolve_protein_display_name(receptor_name, self.protein_name_map)
                protein_label = format_protein_label(receptor_name, self.protein_name_map)
                sdf_path = Path(str(complex_info.get("pose_file") or ""))

                tag_meta[tag] = {
                    "protein_name": protein_name,
                    "protein_label": protein_label,
                    "ligand": ligand_name,
                    "receptor": receptor_name,
                    "sdf_file": str(sdf_path),
                }

            if not tag_meta:
                self.logger.warning("⚠️  No tag metadata available for RMSD analysis")
                return True

            scores_df = self.scores_df.copy()
            scores_df["tag_norm"] = scores_df["tag"].astype(str).apply(_normalize_tag)
            scores_df["mode_int"] = pd.to_numeric(scores_df.get("mode"), errors="coerce")
            enabled_scopes = set(self.rmsd_scopes)
            enable_per_complex = "per_complex" in enabled_scopes
            enable_per_protein = False
            enable_global = False
            for requested, directory in (("per_protein", per_protein_dir), ("global", global_best_dir)):
                if requested in enabled_scopes:
                    (directory / "not_evaluable.json").write_text(json.dumps({"status":"not_evaluable", "reason":"cross_ligand_or_receptor_RMSD_requires_explicit_atom_and_frame_mapping"}), encoding="utf-8")

            per_complex_rows: List[Dict[str, object]] = []
            processed_complexes = 0

            if enable_per_complex:
                for tag, meta in sorted(tag_meta.items()):
                    sdf_file = Path(meta["sdf_file"])
                    if not sdf_file.exists():
                        continue

                    supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=True, sanitize=True)
                    mols = list(supplier)  # Preserve original record indices, including invalid records.
                    if len(mols) < 2:
                        continue

                    tag_scores = scores_df[scores_df["tag_norm"] == tag].copy()
                    tag_scores = tag_scores.sort_values("mode_int", na_position="last")
                    pose_modes = list(range(1, len(mols) + 1))

                    n_poses = min(len(mols), len(pose_modes))
                    if n_poses < 2:
                        continue

                    mols = mols[:n_poses]
                    pose_modes = pose_modes[:n_poses]
                    mode_labels = [f"mode_{mode}" for mode in pose_modes]

                    rmsd_matrix = np.full((n_poses, n_poses), np.nan, dtype=float)
                    np.fill_diagonal(rmsd_matrix, 0.0)
                    for i in range(n_poses):
                        for j in range(i + 1, n_poses):
                            try:
                                rmsd = fixed_frame_rmsd(mols[i], mols[j])
                            except Exception:
                                rmsd = np.nan
                            rmsd_matrix[i, j] = rmsd
                            rmsd_matrix[j, i] = rmsd

                    pair_values = _pairwise_values(rmsd_matrix)
                    if len(pair_values) == 0:
                        continue

                    best_mode = pose_modes[0]
                    best_affinity = np.nan
                    if not tag_scores.empty and "vina_affinity" in tag_scores.columns:
                        best_row = tag_scores.loc[tag_scores["vina_affinity"].idxmin()]
                        try:
                            best_mode = int(best_row.get("mode_int"))
                        except Exception:
                            best_mode = pose_modes[0]
                        best_affinity = float(best_row.get("vina_affinity", np.nan))
                    if best_mode not in pose_modes:
                        best_mode = pose_modes[0]
                    best_idx = pose_modes.index(best_mode)

                    rmsd_to_best = []
                    for idx in range(n_poses):
                        if idx == best_idx:
                            continue
                        value = rmsd_matrix[best_idx, idx]
                        if np.isfinite(value):
                            rmsd_to_best.append(float(value))

                    tag_slug = self._slugify_name(tag)
                    complex_scope_dir = per_complex_dir / tag_slug
                    complex_scope_dir.mkdir(parents=True, exist_ok=True)
                    matrix_df = pd.DataFrame(rmsd_matrix, index=mode_labels, columns=mode_labels)
                    matrix_file = complex_scope_dir / "all_pose_rmsd_matrix.csv"
                    matrix_df.to_csv(matrix_file)

                    summary_row = {
                        "tag": tag,
                        "protein_name": meta["protein_name"],
                        "protein_label": meta["protein_label"],
                        "ligand": meta["ligand"],
                        "receptor": meta["receptor"],
                        "sdf_file": str(sdf_file),
                        "pose_count": n_poses,
                        "best_mode": int(best_mode),
                        "best_affinity": _finite_float(best_affinity),
                        "pairwise_count": int(len(pair_values)),
                        "mean_rmsd_to_best_pose": float(np.mean(rmsd_to_best)) if rmsd_to_best else None,
                        "median_rmsd_to_best_pose": float(np.median(rmsd_to_best)) if rmsd_to_best else None,
                        "affinity_category": self._affinity_category(best_affinity),
                        "affinity_category_label": self._format_category_label(self._affinity_category(best_affinity)),
                        **_summarize_values(pair_values),
                    }
                    pd.DataFrame([summary_row]).to_csv(
                        complex_scope_dir / "all_pose_rmsd_summary.csv",
                        index=False,
                    )
                    _write_json(complex_scope_dir / "all_pose_rmsd_summary.json", summary_row)

                    if plotting_available:
                        _save_matrix_heatmap(
                            matrix_df,
                            complex_scope_dir / "all_pose_rmsd_heatmap.png",
                            f"All-Pose RMSD: {meta['protein_label']} | {meta['ligand']}",
                        )

                    per_complex_rows.append(summary_row)
                    processed_complexes += 1
            else:
                self.logger.info("⏭️  Skipping per-complex RMSD scope")

            per_complex_summary_file = per_complex_dir / "per_complex_rmsd_summary.csv"
            if per_complex_rows:
                per_complex_df = pd.DataFrame(per_complex_rows).sort_values(
                    ["protein_label", "mean_pairwise_rmsd", "ligand"],
                    na_position="last",
                )
                per_complex_df.to_csv(per_complex_summary_file, index=False)
            else:
                pd.DataFrame(
                    [
                        {
                            "state": "skipped" if not enable_per_complex else "no_valid_pose_sets",
                            "reason": (
                                "per_complex scope disabled by configuration"
                                if not enable_per_complex
                                else "No valid pose sets were available for per-complex RMSD analysis"
                            ),
                        }
                    ]
                ).to_csv(per_complex_summary_file, index=False)
                if enable_per_complex:
                    self.logger.warning("⚠️  No valid pose sets for per-complex RMSD analysis")

            # Build best-pose metadata from the generated complex PDB files.
            complex_metadata = self._build_complex_metadata()
            best_pose_rows: List[Dict[str, object]] = []
            for pdb_file in self._get_complex_pdb_files():
                meta = complex_metadata.get(pdb_file.stem, {})
                best_pose_rows.append(
                    {
                        "tag": pdb_file.stem,
                        "protein_name": str(meta.get("protein_name") or ""),
                        "protein_display_name": str(meta.get("protein_name") or ""),
                        "protein_label": str(meta.get("protein_label") or "Unknown Protein"),
                        "ligand": str(meta.get("ligand_display") or self._infer_ligand_display_name(pdb_file.stem)),
                        "vina_affinity": _finite_float(meta.get("best_affinity")),
                        "affinity_category_label": str(meta.get("affinity_category_label") or "Uncategorized"),
                        "pdb_file": pdb_file,
                    }
                )

            best_pose_df = pd.DataFrame(
                best_pose_rows,
                columns=[
                    "tag",
                    "protein_name",
                    "protein_display_name",
                    "protein_label",
                    "ligand",
                    "vina_affinity",
                    "affinity_category_label",
                    "pdb_file",
                ],
            )
            per_protein_rows: List[Dict[str, object]] = []

            if enable_per_protein and not best_pose_df.empty:
                for protein_label, protein_df in best_pose_df.groupby("protein_label", dropna=False):
                    protein_scope_dir = per_protein_dir / self._slugify_name(str(protein_label) or "unknown_protein")
                    protein_scope_df = protein_df.sort_values(["vina_affinity", "ligand"], na_position="last").reset_index(drop=True)
                    protein_summary = _run_best_pose_scope(
                        protein_scope_dir,
                        protein_scope_df,
                        scope_name=f"per_protein_best_poses::{protein_label}",
                        heatmap_title=f"Best-Pose RMSD: {protein_label}",
                    )
                    protein_summary.update(
                        {
                            "protein_name": str(protein_scope_df["protein_name"].iloc[0]) if "protein_name" in protein_scope_df.columns else str(protein_label),
                            "protein_label": str(protein_label),
                            "ligand_count": int(protein_scope_df["ligand"].nunique()) if "ligand" in protein_scope_df.columns else int(len(protein_scope_df)),
                        }
                    )
                    per_protein_rows.append(protein_summary)
            elif not enable_per_protein:
                self.logger.info("⏭️  Skipping per-protein best-pose RMSD scope")

            per_protein_summary_file = per_protein_dir / "per_protein_best_pose_rmsd_summary.csv"
            per_protein_summary_df = pd.DataFrame(per_protein_rows)
            if not per_protein_summary_df.empty:
                per_protein_summary_df = per_protein_summary_df.sort_values(
                    ["mean_pairwise_rmsd", "protein_label"],
                    na_position="last",
                )
                per_protein_summary_df.to_csv(per_protein_summary_file, index=False)
            else:
                pd.DataFrame(
                    [
                        {
                            "state": "skipped" if not enable_per_protein else "no_best_pose_rows",
                            "reason": (
                                "per_protein scope disabled by configuration"
                                if not enable_per_protein
                                else "No best-pose rows were available for per-protein RMSD analysis"
                            ),
                        }
                    ]
                ).to_csv(per_protein_summary_file, index=False)

            if enable_global:
                global_summary = _run_best_pose_scope(
                    global_best_dir,
                    best_pose_df.sort_values(["vina_affinity", "protein_label", "ligand"], na_position="last").reset_index(drop=True),
                    scope_name="global_best_poses",
                    heatmap_title="Global Best-Pose RMSD Across Proteins",
                )
            else:
                self.logger.info("⏭️  Skipping global best-pose RMSD scope")
                global_summary = _write_scope_skip_summary(
                    global_best_dir,
                    "global_best_poses",
                    "global scope disabled by configuration",
                )
            global_summary_csv = global_best_dir / "best_pose_rmsd_summary.csv"
            global_summary_json = global_best_dir / "best_pose_rmsd_summary.json"

            self.results["rmsd_analysis"] = {
                "per_complex_all_poses": str(per_complex_dir),
                "per_complex_summary": str(per_complex_summary_file),
                "per_pair_summary": str(per_complex_summary_file),
                "per_protein_best_poses": str(per_protein_dir),
                "per_protein_summary": str(per_protein_summary_file),
                "global_best_poses": str(global_best_dir),
                "global_summary": global_summary,
                "global_summary_csv": str(global_summary_csv),
                "global_summary_json": str(global_summary_json),
            }

            self.logger.info("✅ RMSD analysis completed")
            self.logger.info(f"   Per-complex all-pose sets processed: {processed_complexes}")
            self.logger.info(f"   Per-protein best-pose groups: {len(per_protein_rows)}")
            if global_summary.get("state") == "deferred":
                self.logger.info(f"   Global best-pose RMSD deferred: {global_summary.get('reason', '')}")
            elif global_summary.get("mean_pairwise_rmsd") is not None:
                self.logger.info(
                    f"   Global best-pose mean pairwise RMSD: {float(global_summary['mean_pairwise_rmsd']):.2f} Å"
                )
            return True

        except Exception as e:
            self.logger.error(f"❌ RMSD analysis error: {e}", exc_info=True)
            return False
    
    def _extract_poses(self) -> bool:
        """Extract best poses and create complexes."""
        self.logger.info("📦 Extracting best poses...")
        
        # Complexes are already created in _create_complexes()
        # This step is mainly for organization and validation
        
        complexes_dir = self.output_dir / "complexes"
        if not complexes_dir.exists():
            self.logger.warning("⚠️  Complexes directory not found, creating...")
            complexes_dir.mkdir(exist_ok=True)
        
        complex_files = self._get_complex_pdb_files()
        self.logger.info(f"✅ Found {len(complex_files)} complex files")
        
        return len(complex_files) > 0
    
    def _generate_reports(self) -> bool:
        """Generate analysis reports."""
        self.logger.info("📄 Generating reports...")
        
        reports_dir = self.output_dir / "reports"
        reports_dir.mkdir(exist_ok=True)
        
        if self.scores_df is not None:
            # Best poses
            best_poses = self._best_poses_from_scores()
            
            best_poses.to_csv(reports_dir / "best_poses.csv", index=False)
            
            # Summary stats
            summary = self.scores_df.groupby('tag').agg({
                'vina_affinity': ['min', 'max', 'mean', 'std']
            }).round(3)
            summary.columns = ['_'.join(col).strip() for col in summary.columns]
            summary.to_csv(reports_dir / "summary_stats.csv")
            
            self.logger.info("✅ Reports generated")
            return True
        
        return False

    def _best_poses_from_scores(self) -> pd.DataFrame:
        """Select the best (lowest affinity) pose per tag from scores."""
        if self.scores_df is None or self.scores_df.empty:
            return pd.DataFrame()
        if not {"tag", "vina_affinity"}.issubset(self.scores_df.columns):
            return pd.DataFrame()

        scores_df = self.scores_df.copy()
        scores_df["tag_norm"] = scores_df["tag"].astype(str).apply(self._normalize_tag)
        scores_df["vina_affinity"] = pd.to_numeric(scores_df["vina_affinity"], errors="coerce")
        scores_df = scores_df[np.isfinite(scores_df["vina_affinity"])].copy()
        if scores_df.empty:
            return pd.DataFrame(columns=[col for col in self.scores_df.columns if col != "tag_norm"])

        best_idx = scores_df.groupby("tag_norm")["vina_affinity"].idxmin()
        return scores_df.loc[best_idx].drop(columns=["tag_norm"], errors="ignore").sort_values(
            ["vina_affinity", "tag"],
            na_position="last",
        )

    @staticmethod
    def _normalize_tag(tag: str) -> str:
        """Normalize tag/complex identifiers for safe joins."""
        return str(tag or "").replace(".log", "").strip()

    @staticmethod
    def _affinity_category(affinity: Optional[float]) -> str:
        """
        Categorize binding affinity into the strong/moderate/weak theme.
        """
        try:
            value = float(affinity)
        except (TypeError, ValueError):
            return "unknown"
        if not np.isfinite(value):
            return "unknown"
        if value <= -7.0:
            return "strong"
        if value <= -5.0:
            return "moderate"
        return "weak"

    @staticmethod
    def _format_category_label(category: str) -> str:
        labels = {
            "strong": "Strong Binder",
            "moderate": "Moderate Binder",
            "weak": "Weak Binder",
            "unknown": "Uncategorized",
        }
        return labels.get(str(category), str(category).capitalize())

    @staticmethod
    def _infer_ligand_display_name(complex_name: str) -> str:
        """
        Infer human-readable ligand name from complex naming patterns.
        """
        text = str(complex_name or "")
        ligand_match = re.search(r"_ligand_([A-Za-z0-9]{1,20})_([A-Za-z])_(\d+)", text)
        if ligand_match:
            return ligand_match.group(1)

        stem = Path(text).stem
        parts = [p for p in stem.split("_") if p]
        for marker in ("Series", "Reference", "Comparative", "Compartive"):
            if marker in parts:
                idx = parts.index(marker)
                if idx + 1 < len(parts):
                    return parts[idx + 1]

        for token in reversed(parts):
            if token.lower() in {"pdbqt", "poses", "top", "out", "cleaned"}:
                continue
            return token
        return "UnknownLigand"

    def _pairlist_metadata_frame(self) -> pd.DataFrame:
        """Normalize pairlist metadata for safe joins onto score-derived tables."""
        if self.pairlist_df is None or self.pairlist_df.empty:
            return pd.DataFrame()

        pair_df = self.pairlist_df.copy()
        required = {"receptor", "site_id", "ligand"}
        if not required.issubset(pair_df.columns):
            return pd.DataFrame()

        pair_df["tag_norm"] = (
            pair_df["receptor"].astype(str)
            + "_"
            + pair_df["site_id"].astype(str)
            + "_"
            + pair_df["ligand"].astype(str)
        ).apply(self._normalize_tag)

        if "ligand" in pair_df.columns:
            pair_df["ligand_display_name"] = pair_df["ligand"].astype(str).apply(
                lambda value: resolve_ligand_display_name(value, self.ligand_name_map)
            )

        keep_cols = [
            "tag_norm",
            "pair_source",
            "selection_mode",
            "is_cocrystal_benchmark",
            "cocrystal_ligand_name",
            "ligand_display_name",
        ]
        existing = [col for col in keep_cols if col in pair_df.columns]
        return pair_df[existing].drop_duplicates("tag_norm", keep="first")

    def _build_complex_metadata(self) -> Dict[str, Dict[str, object]]:
        """
        Build per-complex labels (protein, ligand, affinity, category) for
        consistent visualization naming and titles across tools.
        """
        metadata: Dict[str, Dict[str, object]] = {}

        tag_best_affinity: Dict[str, float] = {}
        best_pose_df = self._best_poses_from_scores()
        if not best_pose_df.empty:
            for _, row in best_pose_df.iterrows():
                tag_best_affinity[self._normalize_tag(str(row.get("tag") or ""))] = float(row["vina_affinity"])

        for complex_info in self.complexes:
            complex_name = str(complex_info.get("complex_name") or "").strip()
            if not complex_name:
                continue
            tag_norm = self._normalize_tag(complex_name)

            receptor_name = str(
                complex_info.get("receptor_name")
                or complex_info.get("receptor_file")
                or ""
            )
            protein_label = format_protein_label(receptor_name or complex_name, self.protein_name_map)
            protein_name = resolve_protein_display_name(receptor_name or complex_name, self.protein_name_map)
            pdb_code = extract_pdb_code(receptor_name) or extract_pdb_code(complex_name) or ""

            ligand_display = str(complex_info.get("ligand_name") or "").strip()
            if not ligand_display or ligand_display.lower().endswith(".pdbqt"):
                ligand_display = self._infer_ligand_display_name(complex_name)

            best_affinity = tag_best_affinity.get(tag_norm)
            category = self._affinity_category(best_affinity)

            metadata[complex_name] = {
                "tag_norm": tag_norm,
                "protein_name": protein_name,
                "protein_label": protein_label,
                "pdb_code": pdb_code,
                "ligand_display": ligand_display,
                "best_affinity": best_affinity,
                "affinity_category": category,
                "affinity_category_label": self._format_category_label(category),
            }

        return metadata

    @staticmethod
    def _best_pose_per_protein_ligand(df: pd.DataFrame) -> pd.DataFrame:
        """Collapse to one best pose per protein-ligand pair when columns are available."""
        if df is None or df.empty:
            return pd.DataFrame()
        if {'protein', 'ligand'}.issubset(df.columns):
            return df.loc[df.groupby(['protein', 'ligand'])['vina_affinity'].idxmin()].copy()
        return df.copy()

    def _best_poses_with_complex_metadata(self) -> pd.DataFrame:
        """
        Enrich score-derived best poses with receptor/ligand metadata from the
        matched complex manifest. This keeps score-based best-pose selection in
        one place while downstream stages reuse consistent labels.
        """
        best_poses = self._best_poses_from_scores()
        if best_poses.empty:
            return pd.DataFrame()

        manifest_rows: List[Dict[str, str]] = []
        for complex_info in self.complexes:
            complex_name = str(complex_info.get("complex_name") or "").strip()
            if not complex_name:
                continue

            receptor_name = str(
                complex_info.get("receptor_name")
                or complex_info.get("receptor_file")
                or ""
            ).strip()
            receptor_stem = Path(receptor_name).stem if receptor_name else ""
            ligand_identifier = str(complex_info.get("ligand_name") or "").strip()
            ligand_name = resolve_ligand_display_name(
                ligand_identifier or complex_name,
                self.ligand_name_map,
            )

            manifest_rows.append(
                {
                    "tag_norm": self._normalize_tag(complex_name),
                    "complex_name": complex_name,
                    "protein": receptor_stem or receptor_name,
                    "receptor": receptor_name or receptor_stem,
                    "ligand": ligand_name,
                    "ligand_identifier": ligand_identifier or ligand_name,
                }
            )

        if not manifest_rows:
            return best_poses

        manifest_df = pd.DataFrame(manifest_rows).drop_duplicates("tag_norm", keep="first")
        best_poses = best_poses.copy()
        best_poses["tag_norm"] = best_poses["tag"].astype(str).apply(self._normalize_tag)
        enriched = best_poses.merge(manifest_df, on="tag_norm", how="left")
        pair_meta = self._pairlist_metadata_frame()
        if not pair_meta.empty:
            enriched = enriched.merge(pair_meta, on="tag_norm", how="left", suffixes=("", "_pair"))
            if "ligand_display_name_pair" in enriched.columns and "ligand_display_name" not in enriched.columns:
                enriched = enriched.rename(columns={"ligand_display_name_pair": "ligand_display_name"})
        enriched = enriched.drop(columns=["tag_norm"], errors="ignore")
        return self._apply_display_names(enriched, protein_col="protein")

    def _get_complex_pdb_files(self) -> List[Path]:
        """Return available complex PDB files under the pipeline output."""
        complexes_dir = self.output_dir / "complexes"
        if not complexes_dir.exists():
            return []
        return sorted(path for path in complexes_dir.glob("*.pdb") if not path.name.endswith(".receptor.pdb"))

    def _get_best_poses_for_visualization(self) -> pd.DataFrame:
        """Return best poses from current results with display names applied."""
        best_poses = None
        if 'affinity_analysis' in self.results:
            best_poses = self.results['affinity_analysis'].get('best_poses')
        elif 'hierarchical_analysis' in self.results:
            best_poses = self.results['hierarchical_analysis'].get('best_poses')
        else:
            best_poses = self._best_poses_with_complex_metadata()

        if best_poses is None or best_poses.empty:
            return pd.DataFrame()
        best_poses = self._best_pose_per_protein_ligand(best_poses)
        return self._apply_display_names(best_poses, protein_col="protein")

    def _reference_mask(self, df: pd.DataFrame) -> pd.Series:
        """Return a boolean mask for benchmark/reference-like rows."""
        if df is None or df.empty:
            return pd.Series(dtype=bool)

        ref_mask = pd.Series(False, index=df.index)
        if "is_cocrystal_benchmark" in df.columns:
            ref_mask |= df["is_cocrystal_benchmark"].fillna(False).astype(bool)
        if "is_reference_candidate" in df.columns:
            ref_mask |= df["is_reference_candidate"].fillna(False).astype(bool)
        ref_mask |= df.apply(self._is_reference_candidate_row, axis=1)
        return ref_mask

    def _select_nonreference_shared_ligands(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Select shared non-reference ligands for cross-protein comparisons.

        Falls back to all non-reference ligands when no ligand is present on
        more than one protein target.
        """
        if df is None or df.empty:
            return pd.DataFrame()

        pool = df[~self._reference_mask(df)].copy()
        if pool.empty:
            return pd.DataFrame()

        ligand_counts = pool.groupby("ligand")["protein"].nunique()
        shared_ligands = ligand_counts[ligand_counts >= 2].index.tolist()
        if shared_ligands:
            pool = pool[pool["ligand"].isin(shared_ligands)].copy()
        return pool
    
    def _generate_visualizations(self) -> bool:
        """Generate enhanced visualizations with 2D plots and heatmaps."""
        self.logger.info("📊 Generating enhanced visualizations...")
        
        viz_dir = self.visual_2d_root / "summary"
        viz_dir.mkdir(parents=True, exist_ok=True)
        
        if self.scores_df is None or self.scores_df.empty:
            self.logger.warning("⚠️  No scores data for visualizations")
            return False
        
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns

            best_poses = self._get_best_poses_for_visualization()

            def _save_placeholder_plot(output_file: Path, title: str, message: str) -> None:
                fig, ax = plt.subplots(figsize=(10, 5))
                ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=11)
                ax.set_title(title, fontsize=13, fontweight="bold")
                ax.axis("off")
                plt.tight_layout()
                plt.savefig(output_file, dpi=300, bbox_inches="tight")
                plt.close(fig)
            
            # Set publication style
            plt.style.use('default')
            sns.set_palette("husl")
            
            # 1. Affinity distribution histogram across best protein-ligand pairs.
            fig, ax = plt.subplots(figsize=(12, 7))
            distribution_df = best_poses.copy()
            distribution_df["vina_affinity"] = pd.to_numeric(
                distribution_df["vina_affinity"], errors="coerce"
            )
            distribution_df = distribution_df[np.isfinite(distribution_df["vina_affinity"])].copy()
            favorable = distribution_df[distribution_df['vina_affinity'] <= 0]['vina_affinity']
            unfavorable = distribution_df[distribution_df['vina_affinity'] > 0]['vina_affinity']
            bins = 32
            if not favorable.empty:
                ax.hist(
                    favorable,
                    bins=bins,
                    alpha=0.75,
                    edgecolor='black',
                    color='skyblue',
                    label='Favorable (<= 0)',
                )
            if not unfavorable.empty:
                ax.hist(
                    unfavorable,
                    bins=bins,
                    alpha=0.80,
                    edgecolor='black',
                    color='salmon',
                    label='Unfavorable (> 0)',
                )
            mean_aff = distribution_df['vina_affinity'].mean()
            median_aff = distribution_df['vina_affinity'].median()
            positive_count = int((distribution_df['vina_affinity'] > 0).sum())
            total_count = int(len(distribution_df))
            positive_pct = (100.0 * positive_count / total_count) if total_count else 0.0
            ax.axvline(mean_aff, color='red', linestyle='--', label=f'Mean: {mean_aff:.2f}')
            ax.axvline(median_aff, color='green', linestyle='--', label=f'Median: {median_aff:.2f}')
            ax.axvline(0.0, color='darkorange', linestyle=':', linewidth=2, label='Energetic split (0)')
            ax.text(
                0.02,
                0.95,
                f"Unfavorable (>0): {positive_count}/{total_count} ({positive_pct:.1f}%)",
                transform=ax.transAxes,
                ha='left',
                va='top',
                fontsize=10,
                bbox=dict(facecolor='white', alpha=0.85, edgecolor='gray'),
            )
            ax.set_xlabel('Vina Affinity (kcal/mol)', fontsize=12)
            ax.set_ylabel('Frequency', fontsize=12)
            ax.set_title('Distribution of Best-Pose Affinities (Protein-Ligand Pairs)', fontsize=14, fontweight='bold')
            ax.legend()
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(viz_dir / "affinity_distribution.png", dpi=300, bbox_inches='tight')
            plt.close()

            candidate_df = self._select_nonreference_shared_ligands(best_poses)
            if candidate_df.empty:
                candidate_df = best_poses.copy()
            benchmark_df = best_poses[self._reference_mask(best_poses)].copy()

            # 2. Top non-reference performer per protein with optional benchmark markers.
            if not best_poses.empty and {'protein', 'ligand'}.issubset(best_poses.columns):
                per_protein_best = candidate_df.loc[
                    candidate_df.groupby('protein')['vina_affinity'].idxmin()
                ].copy()
                per_protein_best = per_protein_best.sort_values('vina_affinity')

                benchmark_by_protein = {}
                if not benchmark_df.empty:
                    benchmark_by_protein = (
                        benchmark_df.groupby("protein")["vina_affinity"].min().to_dict()
                    )

                labels = []
                for _, row in per_protein_best.iterrows():
                    protein_label = str(
                        row.get('protein_label')
                        or format_protein_label(str(row.get('protein') or row.get('receptor') or ""), self.protein_name_map)
                    )
                    ligand_name = str(
                        row.get('ligand_display_name')
                        or row.get('ligand_label')
                        or row.get('ligand')
                        or row.get('ligand_name')
                        or row.get('tag')
                        or ""
                    )
                    labels.append(f"{protein_label} | {ligand_name}")

                fig, ax = plt.subplots(figsize=(13, max(6, len(per_protein_best) * 0.45)))
                ax.barh(range(len(per_protein_best)), per_protein_best['vina_affinity'].values, color='lightcoral', edgecolor='black')
                ax.set_yticks(range(len(per_protein_best)))
                ax.set_yticklabels(labels, fontsize=9)
                for idx, (_, row) in enumerate(per_protein_best.iterrows()):
                    benchmark_aff = benchmark_by_protein.get(str(row.get("protein")))
                    if benchmark_aff is not None and np.isfinite(float(benchmark_aff)):
                        ax.scatter(
                            float(benchmark_aff),
                            idx,
                            marker="D",
                            s=34,
                            color="black",
                            zorder=5,
                        )
                ax.set_xlabel('Vina Affinity (kcal/mol)', fontsize=12)
                ax.set_title('Top Performer per Protein', fontsize=14, fontweight='bold')
                ax.axvline(0.0, color='darkorange', linestyle=':', alpha=0.8, label='Unfavorable (>0)')
                if benchmark_by_protein:
                    ax.scatter([], [], marker="D", s=34, color="black", label="Cocrystal benchmark")
                ax.grid(True, alpha=0.3, axis='x')
                ax.legend()
                plt.tight_layout()
                plt.savefig(viz_dir / "top_performers.png", dpi=300, bbox_inches='tight')
                plt.close()

            # 3. Affinity heatmap (protein × ligand) with global labels.
            if not best_poses.empty and {'protein_label', 'ligand'}.issubset(best_poses.columns):
                heatmap_df = best_poses.copy()
                heatmap_df["ligand_plot_label"] = heatmap_df.apply(
                    lambda row: str(
                        row.get("ligand_display_name")
                        or row.get("ligand_label")
                        or row.get("ligand")
                        or ""
                    ),
                    axis=1,
                )
                pivot_data = heatmap_df.pivot_table(
                    values='vina_affinity',
                    index='protein_label',
                    columns='ligand_plot_label',
                    aggfunc='min'
                )

                if len(pivot_data) > 1 and len(pivot_data.columns) > 1:
                    fig, ax = plt.subplots(figsize=(max(12, len(pivot_data.columns) * 0.8), max(8, len(pivot_data) * 0.7)))
                    sns.heatmap(
                        pivot_data,
                        annot=True,
                        fmt='.2f',
                        cmap='RdYlGn_r',
                        cbar_kws={'label': 'Binding Affinity (kcal/mol)'},
                        ax=ax,
                        linewidths=0.5,
                    )
                    ax.set_title('Binding Affinity Heatmap (Protein × Ligand)', fontsize=14, fontweight='bold')
                    ax.set_xlabel('Ligand', fontsize=12)
                    ax.set_ylabel('Protein', fontsize=12)
                    plt.tight_layout()
                    plt.savefig(viz_dir / "affinity_heatmap.png", dpi=300, bbox_inches='tight')
                    plt.close()

            # 4. Clustered affinity distribution per protein with best performer marker.
            if not best_poses.empty and {'protein_label', 'ligand'}.issubset(best_poses.columns):
                clustered = best_poses.copy()
                fig, ax = plt.subplots(figsize=(max(12, len(clustered['protein_label'].unique()) * 1.2), 7))
                sns.boxplot(data=clustered, x='protein_label', y='vina_affinity', ax=ax, color='lightblue')
                sns.stripplot(
                    data=clustered,
                    x='protein_label',
                    y='vina_affinity',
                    ax=ax,
                    color='black',
                    alpha=0.55,
                    jitter=0.20,
                    size=4,
                )
                best_rows = clustered.loc[clustered.groupby('protein_label')['vina_affinity'].idxmin()]
                for _, row in best_rows.iterrows():
                    protein = row['protein_label']
                    xpos = list(clustered['protein_label'].unique()).index(protein)
                    ax.scatter(xpos, row['vina_affinity'], marker='*', s=160, color='gold', edgecolor='black', zorder=6)
                    ax.text(
                        xpos + 0.03,
                        row['vina_affinity'] + 0.12,
                        str(
                            row.get('ligand_display_name')
                            or row.get('ligand_label')
                            or row.get('ligand')
                            or ""
                        ),
                        fontsize=8,
                        ha='left',
                        va='bottom',
                    )
                ax.axhline(-7.0, color='red', linestyle='--', alpha=0.55, label='Strong binder threshold (-7)')
                ax.axhline(0.0, color='darkorange', linestyle=':', alpha=0.85, label='Unfavorable (>0)')
                ax.set_xlabel('Protein', fontsize=12)
                ax.set_ylabel('Best-Pose Affinity by Ligand (kcal/mol)', fontsize=12)
                ax.set_title('Best Binding Affinity by Protein (Ligand Clusters)', fontsize=14, fontweight='bold')
                ax.legend()
                plt.xticks(rotation=40, ha='right')
                plt.tight_layout()
                plt.savefig(viz_dir / "affinity_by_protein.png", dpi=300, bbox_inches='tight')
                plt.close()

            # 5. Strong/moderate/weak theme across proteins.
            if not best_poses.empty and {'protein_label', 'vina_affinity'}.issubset(best_poses.columns):
                category_df = best_poses.copy()
                category_df["affinity_category"] = category_df["vina_affinity"].apply(self._affinity_category)
                category_df["affinity_category_label"] = category_df["affinity_category"].apply(self._format_category_label)
                order = ["Strong Binder", "Moderate Binder", "Weak Binder", "Uncategorized"]
                grouped = (
                    category_df.groupby(["protein_label", "affinity_category_label"])
                    .size()
                    .reset_index(name="count")
                )
                if not grouped.empty:
                    grouped["affinity_category_label"] = pd.Categorical(
                        grouped["affinity_category_label"],
                        categories=order,
                        ordered=True,
                    )
                    pivot_counts = grouped.pivot_table(
                        values="count",
                        index="protein_label",
                        columns="affinity_category_label",
                        aggfunc="sum",
                        fill_value=0,
                    )
                    fig, ax = plt.subplots(figsize=(max(12, len(pivot_counts) * 1.15), 7))
                    pivot_counts = pivot_counts[[c for c in order if c in pivot_counts.columns]]
                    pivot_counts.plot(
                        kind="bar",
                        stacked=True,
                        ax=ax,
                        color=["#2ca02c", "#ffbf00", "#d62728", "#7f7f7f"][: len(pivot_counts.columns)],
                        edgecolor="black",
                    )
                    ax.set_xlabel("Protein")
                    ax.set_ylabel("Number of Ligands (Best Pose per Pair)")
                    ax.set_title("Binding Category Profile by Protein (Strong / Moderate / Weak)")
                    ax.legend(title="Category")
                    plt.xticks(rotation=40, ha="right")
                    plt.tight_layout()
                    plt.savefig(viz_dir / "binding_category_by_protein.png", dpi=300, bbox_inches="tight")
                    plt.close(fig)

            # 6. Cross-protein comparison for shared non-reference ligands.
            if not best_poses.empty and {'protein_label', 'ligand', 'vina_affinity'}.issubset(best_poses.columns):
                cross_df = self._select_nonreference_shared_ligands(best_poses)

                if cross_df.empty:
                    _save_placeholder_plot(
                        viz_dir / "cross_protein_comparison.png",
                        "Cross-Protein Comparison",
                        "No shared non-reference ligands were detected for cross-protein comparison.",
                    )
                else:
                    cross_df = cross_df.loc[
                        cross_df.groupby(["protein_label", "ligand"])["vina_affinity"].idxmin()
                    ].copy()
                    cross_df["ligand_plot_label"] = cross_df.apply(
                        lambda row: str(
                            row.get("ligand_display_name")
                            or row.get("ligand_label")
                            or row.get("ligand")
                            or ""
                        ),
                        axis=1,
                    )
                    matrix = cross_df.pivot_table(
                        values="vina_affinity",
                        index="ligand_plot_label",
                        columns="protein_label",
                        aggfunc="min",
                    )
                    matrix = matrix.dropna(how="all", axis=0).dropna(how="all", axis=1)

                    if matrix.empty:
                        _save_placeholder_plot(
                            viz_dir / "cross_protein_comparison.png",
                            "Cross-Protein Comparison",
                            "Shared non-reference ligands were detected, but no matrix could be formed.",
                        )
                    else:
                        ligand_rank = matrix.mean(axis=1, skipna=True).sort_values()
                        matrix = matrix.loc[ligand_rank.head(min(16, len(ligand_rank))).index]
                        protein_order = matrix.mean(axis=0, skipna=True).sort_values().index.tolist()
                        matrix = matrix[protein_order]
                        top_profile_ligands = (
                            matrix.mean(axis=1, skipna=True).sort_values().head(min(6, len(matrix))).index.tolist()
                        )

                        fig, (ax_heat, ax_profile) = plt.subplots(
                            1,
                            2,
                            figsize=(max(15, len(protein_order) * 1.1 + 8), max(6, len(matrix) * 0.42 + 2)),
                            gridspec_kw={"width_ratios": [1.55, 1.0]},
                        )

                        sns.heatmap(
                            matrix,
                            ax=ax_heat,
                            cmap="RdYlGn_r",
                            annot=True,
                            fmt=".2f",
                            linewidths=0.3,
                            linecolor="white",
                            cbar_kws={"label": "Vina Affinity (kcal/mol)"},
                        )
                        ax_heat.set_xlabel("Protein Target")
                        ax_heat.set_ylabel("Ligand (Top by Mean Affinity)")
                        ax_heat.set_title("Cross-Protein Affinity Heatmap (Shared Non-Reference Ligands)")

                        x_positions = list(range(len(protein_order)))
                        for ligand in top_profile_ligands:
                            row = pd.to_numeric(matrix.loc[ligand].reindex(protein_order), errors="coerce")
                            y_values = row.to_numpy(dtype=float)
                            if int(np.isfinite(y_values).sum()) == 0:
                                continue
                            ax_profile.plot(
                                x_positions,
                                y_values,
                                marker="o",
                                linewidth=2,
                                markersize=5,
                                label=str(ligand),
                            )
                        ax_profile.set_xticks(x_positions)
                        ax_profile.set_xticklabels(protein_order, rotation=45, ha="right")
                        ax_profile.set_xlabel("Protein Target")
                        ax_profile.set_ylabel("Vina Affinity (kcal/mol)")
                        ax_profile.set_title("Top Shared Ligand Profiles")
                        ax_profile.axhline(y=-7.0, color="red", linestyle="--", alpha=0.4)
                        ax_profile.axhline(y=0.0, color="darkorange", linestyle=":", alpha=0.7)
                        if top_profile_ligands:
                            ax_profile.legend(title="Ligand", fontsize=8, title_fontsize=9, loc="best")

                        plt.tight_layout()
                        plt.savefig(viz_dir / "cross_protein_comparison.png", dpi=300, bbox_inches="tight")
                        plt.close(fig)
            else:
                _save_placeholder_plot(
                    viz_dir / "cross_protein_comparison.png",
                    "Cross-Protein Comparison",
                    "Insufficient data columns for cross-protein comparison.",
                )

            # 7. Cocrystal/reference benchmark view with best screening context.
            if not best_poses.empty and {'protein_label', 'ligand', 'vina_affinity'}.issubset(best_poses.columns):
                comp_df = best_poses.copy()
                comp_df = comp_df[self._reference_mask(comp_df)].copy()

                if comp_df.empty:
                    _save_placeholder_plot(
                        viz_dir / "comparative_redocking.png",
                        "Comparative / Reference Results",
                        "No cocrystal/reference benchmark rows were detected in this run.",
                    )
                else:
                    comp_df = comp_df.loc[
                        comp_df.groupby(["protein_label", "ligand"])["vina_affinity"].idxmin()
                    ].copy()
                    comp_df["source_class"] = comp_df.apply(
                        lambda row: "Cocrystal Benchmark" if bool(row.get("is_cocrystal_benchmark")) else "Reference/Control",
                        axis=1,
                    )
                    comp_df = comp_df.sort_values("vina_affinity")
                    comp_df["ligand_plot_label"] = comp_df.apply(
                        lambda row: str(
                            row.get("ligand_display_name")
                            or row.get("ligand_label")
                            or row.get("ligand")
                            or ""
                        ),
                        axis=1,
                    )
                    comp_df["plot_label"] = (
                        comp_df["protein_label"].astype(str) + " | " + comp_df["ligand_plot_label"].astype(str)
                    )

                    # Per-protein best screening binder marker for context.
                    best_series_by_protein = {}
                    if not candidate_df.empty:
                        best_series_by_protein = candidate_df.groupby("protein_label")["vina_affinity"].min().to_dict()

                    palette = {
                        "Cocrystal Benchmark": "#4c72b0",
                        "Reference/Control": "#55a868",
                    }
                    colors = comp_df["source_class"].map(palette).fillna("#999999")

                    fig, ax = plt.subplots(figsize=(max(11, len(comp_df) * 0.42), max(6, len(comp_df) * 0.32)))
                    bars = ax.barh(
                        comp_df["plot_label"],
                        comp_df["vina_affinity"],
                        color=colors,
                        edgecolor="black",
                    )
                    for bar, (_, row) in zip(bars, comp_df.iterrows()):
                        val = float(row["vina_affinity"])
                        ax.text(
                            val + (0.14 if val >= 0 else -0.14),
                            bar.get_y() + bar.get_height() / 2,
                            f"{val:.2f}",
                            va="center",
                            ha="left" if val >= 0 else "right",
                            fontsize=8,
                        )
                        series_best = best_series_by_protein.get(str(row["protein_label"]))
                        if series_best is not None and np.isfinite(float(series_best)):
                            ax.scatter(
                                float(series_best),
                                bar.get_y() + bar.get_height() / 2,
                                marker="D",
                                s=28,
                                color="black",
                                zorder=5,
                            )

                    import matplotlib.patches as mpatches

                    legend_handles = []
                    present_classes = [c for c in palette.keys() if c in set(comp_df["source_class"].tolist())]
                    for cls_name in present_classes:
                        legend_handles.append(mpatches.Patch(color=palette[cls_name], label=cls_name))
                    if best_series_by_protein:
                        legend_handles.append(
                            mpatches.Patch(color="black", label="◆ Best Screening Ligand per Protein (marker)")
                        )

                    ax.set_xlabel("Vina Affinity (kcal/mol)")
                    ax.set_ylabel("Protein | Ligand")
                    ax.set_title("Cocrystal / Reference Benchmark Overview")
                    ax.axvline(x=-7.0, color='red', linestyle='--', alpha=0.45)
                    ax.axvline(x=0.0, color='darkorange', linestyle=':', alpha=0.6)
                    if legend_handles:
                        ax.legend(handles=legend_handles, loc="best", fontsize=8, title="Entry Type")
                    plt.tight_layout()
                    plt.savefig(viz_dir / "comparative_redocking.png", dpi=300, bbox_inches='tight')
                    plt.close(fig)
            else:
                _save_placeholder_plot(
                    viz_dir / "comparative_redocking.png",
                    "Comparative / Reference Results",
                    "Insufficient data columns for comparative/reference visualization.",
                )
            
            self.logger.info("✅ Enhanced visualizations generated")
            return True
            
        except ImportError:
            self.logger.warning("⚠️  Matplotlib/Seaborn not available for visualizations")
            return False
        except Exception as e:
            self.logger.warning(f"⚠️  Error generating visualizations: {e}", exc_info=True)
            return False
    
    def _generate_py3dmol_visualizations(self) -> bool:
        """Generate py3Dmol 3D visualizations for complexes."""
        self.logger.info("🌐 Generating py3Dmol 3D visualizations...")

        if not PY3DMOL_AVAILABLE:
            self.logger.warning("⚠️  py3Dmol is not installed in current environment - skipping 3D visualizations")
            return False
        
        complex_pdb_files = self._get_complex_pdb_files()
        if not complex_pdb_files:
            self.logger.warning("⚠️  No complexes found for py3Dmol visualization")
            return False
        complexes_dir = self.output_dir / "complexes"
        
        viz_3d_dir = self.visual_3d_root / "py3dmol"
        viz_3d_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            # Visualize individual complexes
            individual_dir = viz_3d_dir / "individual_complexes"
            created_individual = visualize_all_complexes(
                complexes_dir,
                individual_dir,
                width=800,
                height=600
            )
            
            self.logger.info(f"  ✅ Created {len(created_individual)} individual 3D visualizations")
            
            # Visualize aggregated ligands by protein
            best_pose_df = self._get_best_poses_for_visualization()
            if not best_pose_df.empty:
                aggregated_dir = viz_3d_dir / "aggregated_by_protein"
                created_aggregated = visualize_ligands_by_protein(
                    complexes_dir,
                    self.receptors_folder,
                    best_pose_df,
                    aggregated_dir,
                    width=1000,
                    height=800
                )
                
                self.logger.info(f"  ✅ Created {len(created_aggregated)} aggregated visualizations")
                self.results['py3dmol_visualizations'] = {
                    'individual': created_individual,
                    'aggregated': created_aggregated
                }
            
            return True
            
        except Exception as e:
            self.logger.warning(f"⚠️  py3Dmol visualization error: {e}")
            return False

    def _infer_source_layer_for_complex(self, complex_name: str) -> str:
        normalized = self._normalize_tag(complex_name)
        if "_ligand_" in str(normalized).lower():
            return "redocked_reference_best_pose"
        if isinstance(self.pairlist_df, pd.DataFrame) and not self.pairlist_df.empty:
            pair_df = self._pairlist_metadata_frame()
            if not pair_df.empty and "tag_norm" in pair_df.columns:
                matches = pair_df[pair_df["tag_norm"].astype(str) == normalized]
                if not matches.empty:
                    row = matches.iloc[0]
                    cocrystal_flag = str(row.get("is_cocrystal_benchmark", "")).strip().lower()
                    if cocrystal_flag in {"true", "1", "yes"}:
                        return "redocked_reference_best_pose"
                    pair_source = str(row.get("pair_source", "")).strip().lower()
                    selection_mode = str(row.get("selection_mode", "")).strip().lower()
                    if "comparative" in pair_source or "reference" in pair_source or "benchmark" in pair_source:
                        return "redocked_reference_best_pose"
                    if selection_mode in {"comparative", "reference", "benchmark", "native"}:
                        return "redocked_reference_best_pose"
        return "novel_best_pose"

    def _generate_plip_interaction_outputs(self) -> bool:
        """Run PLIP over complexes and persist XML/diagram artifacts + manifest."""
        self.logger.info("🔬 Generating PLIP interaction outputs...")
        complex_files = self._get_complex_pdb_files()
        if not complex_files:
            self.logger.error("❌ No complexes found for PLIP analysis")
            return False

        output_root = self.interactions_root / "plip"
        output_root.mkdir(parents=True, exist_ok=True)
        analyzer = PLIPAnalyzer(config=dict(self.plip_config))
        if not analyzer.plip_available:
            self.logger.error("❌ PLIP is not installed; mandatory PLIP stage cannot run")
            return False

        max_complexes = self.plip_config.get("max_complexes")
        selected_files = complex_files
        if isinstance(max_complexes, int) and max_complexes > 0:
            selected_files = complex_files[:max_complexes]

        records: List[Dict[str, object]] = []
        successful = 0
        for pdb_file in selected_files:
            result = analyzer.analyze_complex(pdb_file, output_root, complex_name=pdb_file.stem)
            if result.success:
                successful += 1
            records.append(
                {
                    "complex_name": str(result.complex_name),
                    "pdb_file": str(result.pdb_file),
                    "success": bool(result.success),
                    "png_diagram": str(result.png_diagram) if result.png_diagram else "",
                    "pymol_session": str(result.pymol_session) if result.pymol_session else "",
                    "xml_report": str(result.xml_report) if result.xml_report else "",
                    "txt_report": str(result.txt_report) if result.txt_report else "",
                    "binding_site_count": int(result.binding_site_count or 0),
                    "interactions": json.dumps(result.interactions or {}, sort_keys=True),
                    "error_message": str(result.error_message or ""),
                    "source_layer": self._infer_source_layer_for_complex(str(result.complex_name)),
                }
            )

        manifest_csv = output_root / "plip_job_manifest.csv"
        pd.DataFrame(records).to_csv(manifest_csv, index=False)
        summary = {
            "processed": int(len(records)),
            "successful": int(successful),
            "failed": int(max(0, len(records) - successful)),
            "output_dir": str(output_root),
            "manifest_csv": str(manifest_csv),
        }
        summary_file = output_root / "plip_stage_summary.json"
        summary_file.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        self.results["plip_interaction_outputs"] = {
            **summary,
            "summary_file": str(summary_file),
            "records": records,
        }
        if successful <= 0:
            self.logger.error("❌ PLIP stage produced zero successful outputs")
            return False
        self.logger.info(f"✅ PLIP processed {len(records)} complexes ({successful} successful)")
        return True

    @staticmethod
    def _residue_uid(pdb_id: str, chain: str, resname: str, resnum: str) -> str:
        return f"{str(pdb_id).upper()}:{str(chain or 'UNK').upper()}:{str(resname or '').upper()}:{str(resnum or '')}"

    def _parse_plip_xml_to_rows(
        self,
        xml_path: Path,
        complex_name: str,
        source_layer: str,
        metadata: Dict[str, Dict[str, object]],
    ) -> List[Dict[str, object]]:
        rows: List[Dict[str, object]] = []
        if not xml_path.exists():
            return rows
        meta = metadata.get(complex_name, {})
        pdb_id = str(meta.get("pdb_code") or "") or str(meta.get("protein_name") or "")
        if not pdb_id:
            pdb_id = complex_name[:4].upper()
        protein_name = str(meta.get("protein_name") or "")
        protein_label = str(meta.get("protein_label") or "")
        ligand_name = str(meta.get("ligand_display") or self._infer_ligand_display_name(complex_name))
        ligand_class = "redocked_reference" if str(source_layer).startswith("redocked_reference") else "novel"
        try:
            root = ET.parse(xml_path).getroot()
        except Exception:
            return rows
        for site in root.findall(".//bindingsite"):
            interactions = site.find("interactions")
            if interactions is None:
                continue
            for section in list(interactions):
                interaction_type = str(section.tag or "").replace("_interactions", "")
                for item in list(section):
                    residue_name = str(item.findtext("restype") or "").strip().upper()
                    residue_number = str(item.findtext("resnr") or "").strip()
                    chain = str(item.findtext("reschain") or "UNK").strip() or "UNK"
                    if not residue_name or not residue_number:
                        continue
                    rows.append(
                        {
                            "complex_name": complex_name,
                            "target_id": str(pdb_id).upper(),
                            "pdb_id": str(pdb_id).upper(),
                            "protein_name": protein_name,
                            "protein_label": protein_label,
                            "ligand_id": ligand_name,
                            "ligand_name": ligand_name,
                            "ligand_class": ligand_class,
                            "pose_id": "best_pose",
                            "pose_rank": 1,
                            "source_layer": source_layer,
                            "tool": "PLIP",
                            "chain": chain,
                            "residue_name": residue_name,
                            "residue_number": residue_number,
                            "residue_uid": self._residue_uid(str(pdb_id), chain, residue_name, residue_number),
                            "interaction_type": interaction_type,
                            "distance": str(item.findtext("dist") or item.findtext("dist_h-a") or item.findtext("dist_d-a") or ""),
                            "angle": str(item.findtext("angle") or ""),
                            "notes": str(xml_path),
                        }
                    )
        return rows

    def _normalize_interaction_tables(self) -> bool:
        """Emit layered + normalized interaction tables from PLIP and ProLIF outputs."""
        self.logger.info("🧩 Normalizing layered PLIP/ProLIF interaction tables...")
        normalized_root = self.interactions_root / "normalized"
        layered_root = self.interactions_root / "layered"
        normalized_root.mkdir(parents=True, exist_ok=True)
        layered_root.mkdir(parents=True, exist_ok=True)

        metadata = self._build_complex_metadata()
        plip_rows: List[Dict[str, object]] = []
        plip_result = self.results.get("plip_interaction_outputs") or {}
        plip_records = list(plip_result.get("records") or [])
        for record in plip_records:
            if not bool(record.get("success", False)):
                continue
            complex_name = str(record.get("complex_name") or "")
            xml_report = Path(str(record.get("xml_report") or ""))
            source_layer = str(record.get("source_layer") or self._infer_source_layer_for_complex(complex_name))
            if not complex_name or not xml_report.exists():
                continue
            plip_rows.extend(self._parse_plip_xml_to_rows(xml_report, complex_name, source_layer, metadata))

        prolif_rows: List[Dict[str, object]] = []
        prolif_frequency_files = sorted((self.interactions_root / "prolif" / "across_poses" / "tables").glob("*_frequency.csv"))
        for frequency_file in prolif_frequency_files:
            complex_name = frequency_file.stem.replace("_frequency", "")
            source_layer = self._infer_source_layer_for_complex(complex_name)
            meta = metadata.get(complex_name, {})
            pdb_id = str(meta.get("pdb_code") or "") or complex_name[:4].upper()
            protein_name = str(meta.get("protein_name") or "")
            protein_label = str(meta.get("protein_label") or "")
            ligand_name = str(meta.get("ligand_display") or self._infer_ligand_display_name(complex_name))
            ligand_class = "redocked_reference" if str(source_layer).startswith("redocked_reference") else "novel"
            try:
                freq_df = pd.read_csv(frequency_file)
            except Exception:
                continue
            for _, row in freq_df.iterrows():
                residue_name = str(row.get("protein", "") or "")
                interaction_type = str(row.get("interaction", "") or "")
                frequency = row.get("frequency", "")
                prolif_rows.append(
                    {
                        "complex_name": complex_name,
                        "target_id": str(pdb_id).upper(),
                        "pdb_id": str(pdb_id).upper(),
                        "protein_name": protein_name,
                        "protein_label": protein_label,
                        "ligand_id": ligand_name,
                        "ligand_name": ligand_name,
                        "ligand_class": ligand_class,
                        "pose_id": "across_poses",
                        "pose_rank": "",
                        "source_layer": source_layer,
                        "tool": "ProLIF",
                        "chain": "",
                        "residue_name": residue_name,
                        "residue_number": "",
                        "residue_uid": f"{str(pdb_id).upper()}::{residue_name}",
                        "interaction_type": interaction_type,
                        "distance": "",
                        "angle": "",
                        "notes": str(frequency_file),
                        "frequency": frequency,
                    }
                )

        plip_csv = normalized_root / "plip_interactions_long.csv"
        prolif_csv = normalized_root / "prolif_interactions_frequency.csv"
        merged_csv = normalized_root / "interaction_long_table.csv"
        pd.DataFrame(plip_rows).to_csv(plip_csv, index=False)
        pd.DataFrame(prolif_rows).to_csv(prolif_csv, index=False)
        pd.DataFrame(plip_rows + prolif_rows).to_csv(merged_csv, index=False)

        raw_rows = [row for row in plip_rows if str(row.get("source_layer", "")) == "raw_cocrystal"]
        redocked_rows = [row for row in plip_rows if str(row.get("source_layer", "")).startswith("redocked_reference")]
        novel_rows = [
            row
            for row in plip_rows
            if str(row.get("source_layer", "")) not in {"raw_cocrystal"} and not str(row.get("source_layer", "")).startswith("redocked_reference")
        ]
        pd.DataFrame(raw_rows).to_csv(layered_root / "raw_cocrystal_plip_interactions.csv", index=False)
        pd.DataFrame(redocked_rows).to_csv(layered_root / "redocked_reference_plip_interactions.csv", index=False)
        pd.DataFrame(novel_rows).to_csv(layered_root / "novel_plip_interactions.csv", index=False)
        layered_summary_rows = [
            {"layer": "raw_cocrystal", "rows": len(raw_rows)},
            {"layer": "redocked_reference", "rows": len(redocked_rows)},
            {"layer": "novel", "rows": len(novel_rows)},
            {"layer": "all", "rows": len(plip_rows)},
        ]
        layered_summary_csv = layered_root / "layered_plip_summary.csv"
        pd.DataFrame(layered_summary_rows).to_csv(layered_summary_csv, index=False)

        self.results["normalized_interactions"] = {
            "plip_csv": str(plip_csv),
            "prolif_csv": str(prolif_csv),
            "merged_csv": str(merged_csv),
            "layered_summary_csv": str(layered_summary_csv),
            "plip_rows": int(len(plip_rows)),
            "prolif_rows": int(len(prolif_rows)),
        }
        if len(plip_rows) <= 0:
            self.logger.error("❌ Normalization stage found zero PLIP interaction rows")
            return False
        self.logger.info(f"✅ Normalized PLIP rows: {len(plip_rows)} | ProLIF rows: {len(prolif_rows)}")
        return True

    def _build_source_of_truth_interaction_analytics(self) -> bool:
        """Build reference-anchored interaction agreement table from normalized PLIP rows."""
        self.logger.info("📘 Building source-of-truth interaction analytics...")
        if not bool(self.interaction_analytics_config.get("enabled", True)):
            self.logger.info("⏭️  Source-of-truth interaction analytics disabled")
            return True

        normalized = self.results.get("normalized_interactions") or {}
        plip_csv = Path(str(normalized.get("plip_csv") or self.interactions_root / "normalized" / "plip_interactions_long.csv"))
        if not plip_csv.exists():
            self.logger.error(f"❌ Missing normalized PLIP table: {plip_csv}")
            return False

        try:
            plip_df = pd.read_csv(plip_csv)
        except Exception as exc:
            self.logger.error(f"❌ Could not read normalized PLIP table: {exc}")
            return False
        if plip_df.empty:
            self.logger.error("❌ Source-of-truth stage cannot run with empty PLIP table")
            return False

        group_cols = ["complex_name", "target_id", "ligand_id", "source_layer"]
        event_rows: List[Dict[str, object]] = []
        grouped = plip_df.groupby(group_cols, dropna=False)
        event_map: Dict[str, set] = {}
        target_of_complex: Dict[str, str] = {}
        ligand_of_complex: Dict[str, str] = {}
        layer_of_complex: Dict[str, str] = {}
        for keys, block in grouped:
            complex_name = str(keys[0])
            target_id = str(keys[1])
            ligand_id = str(keys[2])
            source_layer = str(keys[3])
            events = {
                f"{str(row.get('interaction_type', ''))}|{str(row.get('residue_uid', ''))}"
                for _, row in block.iterrows()
            }
            event_map[complex_name] = events
            target_of_complex[complex_name] = target_id
            ligand_of_complex[complex_name] = ligand_id
            layer_of_complex[complex_name] = source_layer

        reference_by_target: Dict[str, str] = {}
        for complex_name, target_id in target_of_complex.items():
            layer = str(layer_of_complex.get(complex_name, ""))
            ligand_id = str(ligand_of_complex.get(complex_name, ""))
            if layer.startswith("redocked_reference") or "_ligand_" in ligand_id.lower():
                current = reference_by_target.get(target_id)
                if current is None or len(event_map.get(complex_name, set())) > len(event_map.get(current, set())):
                    reference_by_target[target_id] = complex_name

        best_pose_df = self._best_poses_with_complex_metadata()
        affinity_map: Dict[str, float] = {}
        if isinstance(best_pose_df, pd.DataFrame) and not best_pose_df.empty:
            for _, row in best_pose_df.iterrows():
                tag = self._normalize_tag(str(row.get("tag") or ""))
                try:
                    affinity_map[tag] = float(row.get("vina_affinity"))
                except Exception:
                    continue

        for complex_name, events in event_map.items():
            target_id = target_of_complex.get(complex_name, "")
            reference_complex = reference_by_target.get(target_id, "")
            ref_events = event_map.get(reference_complex, set()) if reference_complex else set()
            matched_events = events.intersection(ref_events) if ref_events else set()
            match_pct = (100.0 * len(matched_events) / len(ref_events)) if ref_events else np.nan
            event_rows.append(
                {
                    "complex_name": complex_name,
                    "target_id": target_id,
                    "ligand_id": ligand_of_complex.get(complex_name, ""),
                    "source_layer": layer_of_complex.get(complex_name, ""),
                    "reference_complex_name": reference_complex,
                    "reference_event_count": int(len(ref_events)),
                    "complex_event_count": int(len(events)),
                    "matched_event_count_vs_reference": int(len(matched_events)),
                    "interaction_event_match_percentage_vs_reference": float(match_pct) if np.isfinite(match_pct) else "",
                    "vina_affinity": affinity_map.get(self._normalize_tag(complex_name), ""),
                }
            )

        source_root = self.interactions_root / "source_of_truth"
        source_root.mkdir(parents=True, exist_ok=True)
        source_csv = source_root / "source_of_truth_interactions.csv"
        summary_txt = source_root / "source_of_truth_summary.txt"
        source_df = pd.DataFrame(event_rows).sort_values(
            ["target_id", "source_layer", "ligand_id", "complex_name"],
            na_position="last",
        )
        source_df.to_csv(source_csv, index=False)
        summary_txt.write_text(
            "\n".join(
                [
                    f"rows: {len(source_df)}",
                    f"targets: {source_df['target_id'].nunique() if not source_df.empty else 0}",
                    f"ligands: {source_df['ligand_id'].nunique() if not source_df.empty else 0}",
                    f"reference_targets: {len(reference_by_target)}",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        self.results["source_of_truth"] = {
            "source_csv": str(source_csv),
            "summary_txt": str(summary_txt),
            "rows": int(len(source_df)),
            "reference_targets": int(len(reference_by_target)),
        }
        if source_df.empty:
            self.logger.error("❌ Source-of-truth stage produced no rows")
            return False
        self.logger.info(f"✅ Source-of-truth table written: {source_csv}")
        return True
    
    def _generate_prolif_interaction_maps(self) -> bool:
        """Generate ProLIF 2D interaction maps for complexes."""
        self.logger.info("📊 Generating ProLIF 2D interaction maps...")

        if not PROLIF_AVAILABLE:
            self.logger.error("❌ ProLIF is not installed; mandatory interaction stage cannot run")
            return False
        
        complex_pdb_files = self._get_complex_pdb_files()
        if not complex_pdb_files:
            self.logger.error("❌ No complexes found for ProLIF interaction maps")
            return False
        complexes_dir = self.output_dir / "complexes"
        
        interaction_maps_dir = self.interactions_root / "prolif"
        across_root = interaction_maps_dir / "across_poses"
        best_root = interaction_maps_dir / "best_pose"
        across_visuals_dir = across_root / "visuals"
        across_tables_dir = across_root / "tables"
        best_visuals_dir = best_root / "visuals"
        best_tables_dir = best_root / "tables"
        summaries_dir = interaction_maps_dir / "summaries"
        across_visuals_dir.mkdir(parents=True, exist_ok=True)
        across_tables_dir.mkdir(parents=True, exist_ok=True)
        best_visuals_dir.mkdir(parents=True, exist_ok=True)
        best_tables_dir.mkdir(parents=True, exist_ok=True)
        summaries_dir.mkdir(parents=True, exist_ok=True)
        prolif_visual_mirror = self.visual_interactions_2d_root / "prolif"
        prolif_visual_mirror.mkdir(parents=True, exist_ok=True)
        
        try:
            complex_metadata = self._build_complex_metadata()
            max_complexes = self.prolif_config.get("max_complexes")
            selected_pdb_files = complex_pdb_files
            if isinstance(max_complexes, int) and max_complexes > 0:
                selected_pdb_files = complex_pdb_files[:max_complexes]
            created_maps = create_interaction_maps_for_all_complexes(
                complexes_dir,
                across_visuals_dir,
                ligand_resname="UNK",
                dpi=self.prolif_config.get("dpi", 300),
                figsize=self.prolif_config.get("figsize", (12, 8)),
                max_complexes=max_complexes,
                overwrite=False,
                write_html=bool(self.prolif_config.get("write_html", True)),
                complex_manifest=self.complexes,
                complex_metadata=complex_metadata,
                poses_per_complex=self.prolif_config.get("poses_per_complex"),
                lignetwork_threshold=float(self.prolif_config.get("lignetwork_threshold", 0.3)),
                count_occurrences=bool(self.prolif_config.get("count_occurrences", False)),
            )
            moved_across_tables = self._move_files_by_suffix(
                across_visuals_dir,
                across_tables_dir,
                {".csv"},
            )

            # Also generate best-pose single-complex maps for cleaner per-ligand interpretation.
            mapper = ProLifInteractionMapper()
            best_pose_maps: Dict[str, str] = {}
            for pdb_file in selected_pdb_files:
                meta = complex_metadata.get(pdb_file.stem, {})
                protein_slug = self._slugify_name(str(meta.get("protein_label") or "Protein"))
                ligand_slug = self._slugify_name(str(meta.get("ligand_display") or self._infer_ligand_display_name(pdb_file.stem)))
                out_name = f"{protein_slug}__{ligand_slug}__{pdb_file.stem}__best_pose_2d.png"
                output_png = best_visuals_dir / out_name[:220]
                ligand_hint, _, _ = self._derive_ligand_identity({"complex_name": pdb_file.stem})
                display_title = (
                    f"ProLIF Best-Pose Map: {meta.get('protein_label') or pdb_file.stem} | "
                    f"Ligand {meta.get('ligand_display') or ligand_hint}"
                )
                if mapper.create_interaction_map(
                    pdb_file,
                    output_png,
                    ligand_resname=ligand_hint,
                    dpi=self.prolif_config.get("dpi", 300),
                    figsize=self.prolif_config.get("figsize", (12, 8)),
                    write_html=bool(self.prolif_config.get("write_html", True)),
                    display_title=display_title,
                ):
                    best_pose_maps[pdb_file.stem] = str(output_png)
            moved_best_tables = self._move_files_by_suffix(
                best_visuals_dir,
                best_tables_dir,
                {".csv"},
            )

            visual_ext = {".png", ".svg", ".pdf", ".html"}
            mirrored_across = self._copy_visual_assets(
                across_visuals_dir,
                prolif_visual_mirror / "across_poses",
                visual_ext,
            )
            mirrored_best = self._copy_visual_assets(
                best_visuals_dir,
                prolif_visual_mirror / "best_pose",
                visual_ext,
            )

            # Build a compact summary table from per-complex frequency files.
            frequency_rows: List[pd.DataFrame] = []
            for frequency_file in sorted(across_tables_dir.glob("*_frequency.csv")):
                try:
                    freq_df = pd.read_csv(frequency_file)
                except Exception:
                    continue
                if freq_df.empty or "frequency" not in freq_df.columns:
                    continue
                complex_stem = frequency_file.stem.replace("_frequency", "")
                freq_df["complex_stem"] = complex_stem
                frequency_rows.append(freq_df)

            summary_paths: Dict[str, str] = {}
            if frequency_rows:
                all_freq_df = pd.concat(frequency_rows, ignore_index=True)
                all_freq_csv = summaries_dir / "across_pose_interaction_frequency.csv"
                all_freq_df.to_csv(all_freq_csv, index=False)
                summary_paths["across_pose_interaction_frequency_csv"] = str(all_freq_csv)

                grouped = (
                    all_freq_df.groupby(["protein", "interaction"], dropna=False)["frequency"]
                    .mean()
                    .reset_index()
                    .sort_values("frequency", ascending=False)
                )
                grouped_csv = summaries_dir / "interaction_frequency_by_residue_type.csv"
                grouped.to_csv(grouped_csv, index=False)
                summary_paths["interaction_frequency_by_residue_type_csv"] = str(grouped_csv)

            stage_summary = {
                "across_pose_maps": len(created_maps),
                "best_pose_maps": len(best_pose_maps),
                "across_poses_visuals_dir": str(across_visuals_dir),
                "across_poses_tables_dir": str(across_tables_dir),
                "best_pose_visuals_dir": str(best_visuals_dir),
                "best_pose_tables_dir": str(best_tables_dir),
                "summaries_dir": str(summaries_dir),
                "visual_mirror_dir": str(prolif_visual_mirror),
                "moved_across_table_files": moved_across_tables,
                "moved_best_pose_table_files": moved_best_tables,
                "mirrored_across_visual_assets": mirrored_across,
                "mirrored_best_pose_visual_assets": mirrored_best,
                **summary_paths,
            }
            summary_file = summaries_dir / "prolif_stage_summary.json"
            with open(summary_file, "w", encoding="utf-8") as handle:
                json.dump(stage_summary, handle, indent=2)

            self.logger.info(f"  ✅ Created {len(created_maps)} across-pose ProLIF maps")
            self.logger.info(f"  ✅ Created {len(best_pose_maps)} best-pose ProLIF maps")
            self.results["prolif_interaction_maps"] = {
                "across_pose_maps": {key: str(path) for key, path in created_maps.items()},
                "best_pose_maps": best_pose_maps,
                "summary_file": str(summary_file),
                **stage_summary,
            }

            if not created_maps and not best_pose_maps:
                self.logger.error("❌ ProLIF stage produced no interaction maps")
                return False
            
            return True
            
        except Exception as e:
            self.logger.warning(f"⚠️  ProLIF interaction map error: {e}")
            return False

    def _resolve_ligplus_root(self) -> Optional[Path]:
        """Resolve LigPlus root from explicit arg or environment."""
        candidates: List[Path] = []
        if self.ligplus_root:
            candidates.append(self.ligplus_root)

        for env_key in ("LIGPLUS_ROOT", "LIGPLUS_HOME"):
            raw_value = os.environ.get(env_key, "").strip()
            if raw_value:
                candidates.append(Path(raw_value).expanduser())

        expanded_candidates: List[Path] = []
        for candidate in candidates:
            expanded_candidates.append(candidate)
            try:
                if candidate.exists() and candidate.is_dir():
                    for child in sorted(candidate.iterdir()):
                        if child.is_dir():
                            expanded_candidates.append(child)
            except Exception:
                continue

        seen: set = set()
        for candidate in expanded_candidates:
            try:
                resolved = candidate.resolve()
            except Exception:
                resolved = candidate
            resolved_key = str(resolved)
            if resolved_key in seen:
                continue
            seen.add(resolved_key)

            if not resolved.exists() or not resolved.is_dir():
                continue

            has_expected_layout = (
                (resolved / "lib" / "params").exists()
                and (
                    (resolved / "lib" / "exe").exists()
                    or (resolved / "lib" / "exe_mac").exists()
                    or (resolved / "lib" / "exe_mac64").exists()
                    or (resolved / "lib" / "exe_linux").exists()
                    or (resolved / "lib" / "exe_linux64").exists()
                    or (resolved / "lib" / "exe_win").exists()
                    or (resolved / "lib" / "exe_win32").exists()
                )
            )

            if has_expected_layout:
                return resolved

        return None

    def _generate_ligplot_diagrams(self) -> bool:
        """Generate LigPlot+ interaction diagrams for complexes."""
        self.logger.info("🧾 Generating LigPlot+ 2D interaction diagrams...")

        complex_files = self._get_complex_pdb_files()
        if not complex_files:
            self.logger.error("❌ No complexes found for LigPlot analysis")
            return False

        ligplus_root = self._resolve_ligplus_root()
        if ligplus_root is None:
            self.logger.error(
                "❌ LigPlus root not configured for mandatory LigPlot stage "
                "(set LIGPLUS_ROOT/LIGPLUS_HOME or pass --ligplus-root)"
            )
            return False

        try:
            runner = LigPlotRunner(
                ligplus_root=ligplus_root,
                contact_type=str(self.ligplot_config.get("contact_type", "2")),
                no_abort=bool(self.ligplot_config.get("no_abort", True)),
                hydrogenate=str(self.ligplot_config.get("hydrogenate", "auto")),
                strip_metals=str(self.ligplot_config.get("strip_metals", "auto")),
            )
            runner.validate()
        except Exception as exc:
            self.logger.error(f"❌ LigPlot setup unavailable for mandatory stage: {exc}")
            return False

        output_root = self.interactions_root / "ligplot"
        output_root.mkdir(parents=True, exist_ok=True)
        complex_metadata = self._build_complex_metadata()

        max_complexes = self.ligplot_config.get("max_complexes")
        if isinstance(max_complexes, int) and max_complexes > 0:
            complex_files = complex_files[:max_complexes]

        overwrite = bool(self.ligplot_config.get("overwrite", False))
        export_formats = self.ligplot_config.get("export_formats", ["png", "pdf"])
        if isinstance(export_formats, str):
            export_formats = [export_formats]
        export_dpi = int(self.ligplot_config.get("export_dpi", 350))

        processed = 0
        successful = 0
        exported = 0
        failed = 0
        failure_reasons: Dict[str, int] = {}

        for pdb_file in complex_files:
            processed += 1
            meta = complex_metadata.get(pdb_file.stem, {})
            protein_slug = self._slugify_name(str(meta.get("protein_label") or "Protein"))
            ligand_slug = self._slugify_name(str(meta.get("ligand_display") or self._infer_ligand_display_name(pdb_file.stem)))
            category_slug = self._slugify_name(str(meta.get("affinity_category_label") or "Uncategorized"))
            target_name = f"{protein_slug}__{ligand_slug}__{category_slug}__{pdb_file.stem}"
            if len(target_name) > 180:
                target_name = f"{target_name[:175]}__id"
            target_dir = output_root / target_name
            target_dir.mkdir(parents=True, exist_ok=True)

            # Hint ligand identity from complex naming to reduce wrong-residue picks.
            lig_hint, _, _ = self._derive_ligand_identity({"complex_name": pdb_file.stem})
            ligand_hint = None if lig_hint in {"UNK", "LIG"} else lig_hint

            ok = runner.run_single(
                pdb_file,
                target_dir,
                ligand_resname_hint=ligand_hint,
                overwrite=overwrite,
            )
            if not ok:
                failed += 1
                run_report = getattr(runner, "last_run_report", {}) or {}
                reason = str(run_report.get("failure_reason") or "unknown failure")
                failure_reasons[reason] = failure_reasons.get(reason, 0) + 1
                continue

            successful += 1
            if export_formats and "none" not in {str(v).lower() for v in export_formats}:
                try:
                    if export_ligplot_outputs(
                        target_dir,
                        export_formats,
                        overwrite=overwrite,
                        dpi=export_dpi,
                    ):
                        exported += 1
                except Exception as exc:
                    self.logger.warning(f"⚠️  LigPlot export warning for {pdb_file.name}: {exc}")

        self.logger.info(f"  ✅ LigPlot processed {processed} complexes")
        self.logger.info(f"  ✅ LigPlot successful runs: {successful}")
        self.logger.info(f"  ✅ LigPlot exported diagrams: {exported}")
        if failed:
            self.logger.warning(f"  ⚠️  LigPlot failed runs: {failed}")

        stage_summary = {
            "ligplus_root": str(ligplus_root),
            "processed": processed,
            "successful": successful,
            "failed": failed,
            "exported": exported,
            "failure_reasons": failure_reasons,
            "output_dir": str(output_root),
        }
        summary_file = output_root / "ligplot_stage_summary.json"
        with open(summary_file, "w", encoding="utf-8") as handle:
            json.dump(stage_summary, handle, indent=2)

        ligplot_visual_mirror = self.visual_interactions_2d_root / "ligplot"
        ligplot_visual_mirror.mkdir(parents=True, exist_ok=True)
        mirrored_ligplot_assets = self._copy_visual_assets(
            output_root,
            ligplot_visual_mirror,
            {".png", ".svg", ".pdf", ".html"},
        )

        if failure_reasons:
            failure_rows = [
                {"reason": reason, "count": count}
                for reason, count in sorted(failure_reasons.items(), key=lambda item: item[1], reverse=True)
            ]
            pd.DataFrame(failure_rows).to_csv(output_root / "ligplot_failure_reasons.csv", index=False)

        self.results["ligplot_diagrams"] = {
            "ligplus_root": str(ligplus_root),
            "processed": processed,
            "successful": successful,
            "failed": failed,
            "exported": exported,
            "failure_reasons": failure_reasons,
            "output_dir": str(output_root),
            "summary_file": str(summary_file),
            "visual_mirror_dir": str(ligplot_visual_mirror),
            "mirrored_visual_assets": mirrored_ligplot_assets,
        }
        if processed > 0 and successful == 0:
            self.logger.error("❌ LigPlot stage produced zero successful diagrams")
            return False
        return True

    def _generate_poseview_diagrams(self) -> bool:
        """Generate PoseView API interaction diagrams from best-poses or complexes."""
        if not self.poseview_config.get("enabled", False):
            self.logger.info("⏭️  PoseView stage disabled (enable with --enable-poseview)")
            return True

        self.logger.info("🧭 Generating PoseView interaction diagrams...")
        best_poses_dir = self.output_dir / "best_poses"
        complexes_dir = self.output_dir / "complexes"

        poseview_input_dir: Optional[Path] = None
        if best_poses_dir.exists():
            poseview_input_dir = best_poses_dir
        elif complexes_dir.exists():
            # Fallback when run used basic affinity mode (no categorized best_poses tree).
            poseview_input_dir = complexes_dir
            self.logger.info("ℹ️  best_poses directory not found; using complexes directory for PoseView")
        else:
            self.logger.warning("⚠️  Neither best_poses nor complexes directory found - skipping PoseView")
            return True

        poseview_dir = self.interactions_root / "poseview"
        poseview_output_dir = poseview_dir / "2d_diagrams"
        poseview_dir.mkdir(parents=True, exist_ok=True)
        poseview_output_dir.mkdir(parents=True, exist_ok=True)

        try:
            config = {
                "output_formats": self.poseview_config.get("output_formats", ["png", "svg", "pdf"]),
                "poll_interval_seconds": int(self.poseview_config.get("poll_interval_seconds", 5)),
                "timeout_seconds": int(self.poseview_config.get("timeout_seconds", 360)),
                "max_retries": int(self.poseview_config.get("max_retries", 2)),
                "max_complexes": int(self.poseview_config.get("max_complexes", 24)),
                "upload_poll_interval_seconds": int(self.poseview_config.get("upload_poll_interval_seconds", 2)),
                "upload_timeout_seconds": int(self.poseview_config.get("upload_timeout_seconds", 180)),
                "max_ligand_candidates": int(self.poseview_config.get("max_ligand_candidates", 8)),
            }
            summary = run_poseview_analysis(poseview_input_dir, poseview_output_dir, config=config)
            summary_file = poseview_dir / "poseview_stage_summary.json"
            with open(summary_file, "w", encoding="utf-8") as handle:
                json.dump(summary, handle, indent=2)

            poseview_visual_mirror = self.visual_interactions_2d_root / "poseview"
            poseview_visual_mirror.mkdir(parents=True, exist_ok=True)
            mirrored_assets = self._copy_visual_assets(
                poseview_output_dir,
                poseview_visual_mirror,
                {".png", ".svg", ".pdf", ".html"},
            )

            self.results["poseview_summary"] = summary
            self.results["poseview_summary"]["summary_file"] = str(summary_file)
            self.results["poseview_summary"]["visual_mirror_dir"] = str(poseview_visual_mirror)
            self.results["poseview_summary"]["mirrored_visual_assets"] = mirrored_assets
            analyzed = int(summary.get("total_analyzed", 0))
            successful = int(summary.get("successful", 0))
            failed = int(summary.get("failed", max(0, analyzed - successful)))
            self.logger.info(f"  ✅ PoseView processed {analyzed} entries ({successful} successful)")
            self.logger.info(f"  ✅ PoseView visual mirror assets: {mirrored_assets}")
            if failed:
                self.logger.warning(f"  ⚠️  PoseView failed entries: {failed}")
                failure_reasons = summary.get("failure_reasons") or {}
                if isinstance(failure_reasons, dict) and failure_reasons:
                    top_reasons = sorted(
                        ((str(reason), int(count)) for reason, count in failure_reasons.items()),
                        key=lambda item: item[1],
                        reverse=True,
                    )
                    for reason, count in top_reasons[:5]:
                        self.logger.warning(f"     - {reason} ({count})")
            return True
        except Exception as exc:
            self.logger.warning(f"⚠️  PoseView stage failed: {exc}")
            return True
    
    def _generate_pandamap_analysis(self) -> bool:
        """Generate publication-quality PandaMap interaction analysis."""
        self.logger.info("🐼 Generating publication-quality PandaMap analysis...")
        
        complex_pdb_files = self._get_complex_pdb_files()
        if not complex_pdb_files:
            self.logger.warning("⚠️  No complexes found for PandaMap analysis")
            return False
        complexes_dir = self.output_dir / "complexes"
        
        pandamap_dir = self.interactions_root / "pandamap"
        pandamap_dir.mkdir(parents=True, exist_ok=True)
        self._seed_pandamap_cache_from_previous_session(pandamap_dir)
        
        try:
            pandamap_cfg = dict(self.pandamap_config)
            pandamap_cfg["protein_name_map"] = self.protein_name_map
            pandamap_cfg["complex_metadata"] = self._build_complex_metadata()
            summary = run_publication_pandamap_analysis(
                complexes_dir=complexes_dir,
                output_dir=pandamap_dir,
                ligand_name="UNK",
                conda_env="pandamap",
                config=pandamap_cfg,
                max_complexes=int(self.pandamap_config.get("max_complexes", 30))
            )
            
            if summary:
                maps_2d_dir = None
                for candidate in (pandamap_dir / "maps_2d", pandamap_dir / "2d_interaction_maps"):
                    if candidate.exists():
                        maps_2d_dir = candidate
                        break
                maps_3d_dir = None
                for candidate in (pandamap_dir / "maps_3d", pandamap_dir / "3d_visualizations"):
                    if candidate.exists():
                        maps_3d_dir = candidate
                        break

                pandamap_visual_2d = self.visual_interactions_2d_root / "pandamap"
                pandamap_visual_3d = self.visual_interactions_3d_root / "pandamap"
                pandamap_visual_2d.mkdir(parents=True, exist_ok=True)
                pandamap_visual_3d.mkdir(parents=True, exist_ok=True)
                mirrored_2d = self._copy_visual_assets(
                    maps_2d_dir,
                    pandamap_visual_2d,
                    {".png", ".svg", ".pdf", ".html"},
                ) if maps_2d_dir else 0
                mirrored_3d = self._copy_visual_assets(
                    maps_3d_dir,
                    pandamap_visual_3d,
                    {".png", ".svg", ".pdf", ".html"},
                ) if maps_3d_dir else 0

                summary["maps_2d_dir"] = str(maps_2d_dir) if maps_2d_dir else ""
                summary["maps_3d_dir"] = str(maps_3d_dir) if maps_3d_dir else ""
                summary["visual_mirror_2d_dir"] = str(pandamap_visual_2d)
                summary["visual_mirror_3d_dir"] = str(pandamap_visual_3d)
                summary["mirrored_2d_assets"] = mirrored_2d
                summary["mirrored_3d_assets"] = mirrored_3d
                qc_ok, qc_summary = self._evaluate_pandamap_quality_contract(summary)
                summary["quality_control"] = qc_summary

                self.logger.info(f"✅ PandaMap analysis completed")
                self.logger.info(f"   📊 Generated {summary.get('generated_2d_maps', 0)} 2D maps")
                self.logger.info(f"   🌐 Generated {summary.get('generated_3d_visualizations', 0)} 3D visualizations")
                self.logger.info(f"   ✅ PandaMap mirrored assets: 2D={mirrored_2d}, 3D={mirrored_3d}")
                if not qc_ok:
                    summary["quality_gate_failed"] = True
                    self.results["pandamap_summary"] = summary
                    self.logger.error(
                        "❌ PandaMap QC gate failed: %s",
                        "; ".join(qc_summary.get("failure_reasons", []) or ["unknown quality gate error"]),
                    )
                    return False
                self.results['pandamap_summary'] = summary
                return True
            else:
                self.logger.warning("⚠️  PandaMap analysis returned no results")
                return False
                
        except Exception as e:
            self.logger.warning(f"⚠️  PandaMap analysis error: {e}")
            # Don't fail the pipeline if PandaMap is not available
            return True

    def _evaluate_pandamap_quality_contract(self, summary: Dict[str, Any]) -> Tuple[bool, Dict[str, Any]]:
        """
        Evaluate PandaMap chemistry/readability quality gates.

        The stage fails when ligand identity confidence is too low or when
        report support exists but no report outputs are produced.
        """
        analysis_results = summary.get("analysis_results", [])
        if not isinstance(analysis_results, list):
            analysis_results = []

        total_entries = len(analysis_results)
        fallback_entries = 0
        report_entries = 0
        for row in analysis_results:
            if not isinstance(row, dict):
                continue
            if str(row.get("ligand_resolution_source", "")).strip() == "generic_fallback":
                fallback_entries += 1
            if str(row.get("report_file", "")).strip():
                report_entries += 1

        non_generic_ratio = 0.0
        if total_entries > 0:
            non_generic_ratio = (total_entries - fallback_entries) / float(total_entries)
        min_non_generic_ratio = float(self.pandamap_config.get("qc_min_non_generic_ligand_ratio", 0.60))

        runtime_capabilities = summary.get("pandamap_runtime_capabilities", {})
        if not isinstance(runtime_capabilities, dict):
            runtime_capabilities = {}
        option_matrix = runtime_capabilities.get("option_matrix", {})
        if not isinstance(option_matrix, dict):
            option_matrix = {}
        report_supported = bool(option_matrix.get("text_report", False))
        require_reports = bool(self.pandamap_config.get("qc_require_reports_when_supported", True))

        failures: List[str] = []
        if total_entries > 0 and non_generic_ratio < min_non_generic_ratio:
            failures.append(
                f"non_generic_ligand_ratio={non_generic_ratio:.2f} below minimum {min_non_generic_ratio:.2f}"
            )
        if require_reports and report_supported and total_entries > 0 and report_entries == 0:
            failures.append("installed PandaMap supports text reports but no report_file outputs were generated")

        qc_summary = {
            "status": "passed" if not failures else "failed",
            "checks": {
                "total_entries": total_entries,
                "generic_fallback_entries": fallback_entries,
                "non_generic_ligand_ratio": round(non_generic_ratio, 4),
                "minimum_non_generic_ligand_ratio": round(min_non_generic_ratio, 4),
                "report_supported": report_supported,
                "report_required_when_supported": require_reports,
                "entries_with_report_file": report_entries,
            },
            "failure_reasons": failures,
        }
        return len(failures) == 0, qc_summary

    def _resolve_session_context(self) -> Optional[Tuple[Path, Path, Path]]:
        """
        Resolve (sessions_root, session_root, relative_output_path) for sessionized outputs.
        Returns None when output_dir is not under analysis/sessions/<session_id>/...
        """
        sessions_root: Optional[Path] = None
        for candidate in [self.output_dir, *self.output_dir.parents]:
            if candidate.name == "sessions":
                sessions_root = candidate
                break
        if sessions_root is None:
            return None
        try:
            relative = self.output_dir.relative_to(sessions_root)
        except Exception:
            return None
        if not relative.parts:
            return None
        session_root = sessions_root / relative.parts[0]
        relative_output = Path(*relative.parts[1:]) if len(relative.parts) > 1 else Path(".")
        return sessions_root, session_root, relative_output

    @staticmethod
    def _directory_has_assets(directory: Path, suffixes: set[str]) -> bool:
        if not directory.exists():
            return False
        for file_path in directory.rglob("*"):
            if file_path.is_file() and file_path.suffix.lower() in suffixes:
                return True
        return False

    def _seed_pandamap_cache_from_previous_session(self, pandamap_dir: Path) -> None:
        """
        If current session has no PandaMap outputs yet, copy assets from the most
        recent previous session (same relative output scope) to avoid re-rendering.
        """
        overwrite = bool(self.pandamap_config.get("overwrite", False))
        if overwrite:
            return
        asset_suffixes = {".pdf", ".svg", ".png", ".html", ".json", ".csv", ".txt"}
        if self._directory_has_assets(pandamap_dir, asset_suffixes):
            return

        context = self._resolve_session_context()
        if context is None:
            return
        sessions_root, session_root, relative_output = context

        prior_sessions = sorted(
            [path for path in sessions_root.iterdir() if path.is_dir() and path != session_root],
            key=lambda path: path.stat().st_mtime,
            reverse=True,
        )
        for prior_session in prior_sessions:
            source_root = prior_session / relative_output if str(relative_output) != "." else prior_session
            source_pandamap = source_root / "interactions" / "pandamap"
            if not self._directory_has_assets(source_pandamap, asset_suffixes):
                continue
            copied = self._copy_visual_assets(source_pandamap, pandamap_dir, asset_suffixes)
            if copied <= 0:
                continue
            self.logger.info(
                "↪️ Seeded PandaMap cache from previous session: %s (copied %s files)",
                source_pandamap,
                copied,
            )
            return

    @staticmethod
    def _slugify_name(value: str) -> str:
        clean = []
        for ch in str(value):
            if ch.isalnum():
                clean.append(ch)
            elif clean and clean[-1] != "_":
                clean.append("_")
        return "".join(clean).strip("_") or "unknown"

    def _alias_visual_filename(self, filename: str) -> str:
        """
        Prefix visualization filenames with resolved protein display names when
        identifiers include PDB-like codes.
        """
        path = Path(filename)
        stem = path.stem
        suffix = path.suffix
        pdb_code = extract_pdb_code(stem)
        if not pdb_code:
            return filename

        display_name = format_protein_label(stem, self.protein_name_map)
        display_slug = self._slugify_name(display_name)
        stem_lower = stem.lower()
        if display_slug.lower() in stem_lower:
            return filename
        return f"{display_slug}__{stem}{suffix}"

    @staticmethod
    def _copy_with_collision(source: Path, destination: Path) -> Path:
        """Copy file while preserving existing outputs via deterministic suffixes."""
        destination.parent.mkdir(parents=True, exist_ok=True)
        if not destination.exists():
            shutil.copy2(source, destination)
            return destination

        counter = 2
        while True:
            candidate = destination.with_name(f"{destination.stem}__{counter}{destination.suffix}")
            if not candidate.exists():
                shutil.copy2(source, candidate)
                return candidate
            counter += 1

    @staticmethod
    def _move_files_by_suffix(source_dir: Path, target_dir: Path, suffixes: set[str]) -> int:
        """Move files with matching suffixes from source to target directory."""
        moved = 0
        if not source_dir.exists():
            return moved
        target_dir.mkdir(parents=True, exist_ok=True)
        for file_path in sorted(source_dir.glob("*")):
            if not file_path.is_file():
                continue
            if file_path.suffix.lower() not in suffixes:
                continue
            destination = target_dir / file_path.name
            if destination.exists():
                counter = 2
                while True:
                    candidate = destination.with_name(
                        f"{destination.stem}__{counter}{destination.suffix}"
                    )
                    if not candidate.exists():
                        destination = candidate
                        break
                    counter += 1
            try:
                file_path.rename(destination)
                moved += 1
            except Exception:
                continue
        return moved

    @staticmethod
    def _copy_visual_assets(source_dir: Path, target_dir: Path, suffixes: set[str]) -> int:
        """Copy visual assets from source subtree into target subtree."""
        copied = 0
        if not source_dir.exists():
            return copied
        for file_path in sorted(source_dir.rglob("*")):
            if not file_path.is_file():
                continue
            if file_path.suffix.lower() not in suffixes:
                continue
            relative = file_path.relative_to(source_dir)
            destination = target_dir / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists():
                continue
            try:
                shutil.copy2(file_path, destination)
                copied += 1
            except Exception:
                continue
        return copied

    def _write_start_here_file(self, visual_manifest_path: Path, raw_inventory_path: Path) -> None:
        """Write a compact navigation guide for end users."""
        start_here = self.output_dir / "START_HERE.md"
        lines = [
            "# Post-Docking Output Guide",
            "",
            "## Quick Navigation",
            f"- `visualizations/2d/summary`: primary 2D affinity and ranking figures",
            f"- `visualizations/2d/interactions`: 2D interaction diagrams (ProLIF/PoseView/LigPlot/PandaMap)",
            f"- `visualizations/3d`: interactive and static 3D visualizations (py3Dmol/PyMOL/PandaMap)",
            f"- `interactions/`: raw interaction-tool outputs and per-tool metadata",
            f"- `analysis/`, `reports/`, `rmsd_analysis/`: processed analysis tables and reports",
            f"- `{raw_inventory_path.relative_to(self.output_dir)}`: raw artifact inventory",
            f"- `{visual_manifest_path.relative_to(self.output_dir)}`: visualization manifest",
            "",
            "## Notes",
            "- `raw_data/` contains indexed copies of key tables/logs/metadata for streamlined review.",
            "- Empty directories are automatically pruned at the end of each run.",
            "",
        ]
        start_here.write_text("\n".join(lines), encoding="utf-8")

    @staticmethod
    def _prune_empty_directories(root: Path) -> int:
        """Remove empty directories bottom-up and return number removed."""
        removed = 0
        all_dirs = sorted(
            (path for path in root.rglob("*") if path.is_dir()),
            key=lambda path: len(path.parts),
            reverse=True,
        )
        for directory in all_dirs:
            try:
                if directory == root:
                    continue
                if any(directory.iterdir()):
                    continue
                directory.rmdir()
                removed += 1
            except Exception:
                continue
        return removed

    def _organize_output_artifacts(self) -> bool:
        """
        Consolidate outputs into a predictable navigation/index structure.
        """
        self.logger.info("🗂️  Consolidating outputs into canonical visualization/raw-data indexes ...")

        try:
            self.visualizations_root.mkdir(parents=True, exist_ok=True)
            self.visual_2d_root.mkdir(parents=True, exist_ok=True)
            self.visual_3d_root.mkdir(parents=True, exist_ok=True)
            self.visual_interactions_2d_root.mkdir(parents=True, exist_ok=True)
            self.visual_interactions_3d_root.mkdir(parents=True, exist_ok=True)
            self.raw_data_root.mkdir(parents=True, exist_ok=True)

            visual_ext = {".png", ".svg", ".pdf", ".html", ".pse", ".jpg", ".jpeg", ".gif"}
            raw_ext = {
                ".csv",
                ".json",
                ".txt",
                ".xlsx",
                ".log",
                ".md",
                ".pml",
                ".drw",
                ".sum",
                ".ps",
                ".pdb",
                ".pdbqt",
                ".sdf",
                ".mol2",
            }

            # Backfill visualization mirrors from raw interaction outputs when needed.
            mirror_jobs = [
                (self.interactions_root / "prolif", self.visual_interactions_2d_root / "prolif"),
                (self.interactions_root / "poseview", self.visual_interactions_2d_root / "poseview"),
                (self.interactions_root / "ligplot", self.visual_interactions_2d_root / "ligplot"),
                (self.interactions_root / "pandamap" / "maps_2d", self.visual_interactions_2d_root / "pandamap"),
                (self.interactions_root / "pandamap" / "2d_interaction_maps", self.visual_interactions_2d_root / "pandamap"),
                (self.interactions_root / "pandamap" / "maps_3d", self.visual_interactions_3d_root / "pandamap"),
                (self.interactions_root / "pandamap" / "3d_visualizations", self.visual_interactions_3d_root / "pandamap"),
                (self.output_dir / "pymol_visualizations", self.visual_3d_root / "pymol"),
                (self.output_dir / "3d_visualizations", self.visual_3d_root / "py3dmol"),
            ]
            mirrored_assets = 0
            for src_dir, dst_dir in mirror_jobs:
                mirrored_assets += self._copy_visual_assets(src_dir, dst_dir, visual_ext)

            manifest_rows = []
            for visual_file in sorted(self.visualizations_root.rglob("*")):
                if not visual_file.is_file():
                    continue
                if visual_file.suffix.lower() not in visual_ext:
                    continue
                rel_visual = visual_file.relative_to(self.output_dir)
                rel_parts = rel_visual.parts
                if len(rel_parts) >= 3:
                    visual_type = f"{rel_parts[1]}/{rel_parts[2]}"
                elif len(rel_parts) >= 2:
                    visual_type = rel_parts[1]
                else:
                    visual_type = "general"
                source_name = visual_file.stem
                manifest_rows.append(
                    {
                        "visualization_file": str(rel_visual),
                        "visualization_type": visual_type,
                        "source_file": str(rel_visual),
                        "pdb_code": extract_pdb_code(source_name) or "",
                        "protein_display_name": format_protein_label(source_name, self.protein_name_map),
                    }
                )

            raw_bucket_dirs = {
                "tables": self.raw_data_root / "tables",
                "logs": self.raw_data_root / "logs",
                "structures": self.raw_data_root / "structures",
                "metadata": self.raw_data_root / "metadata",
                "tool_files": self.raw_data_root / "tool_files",
            }
            for bucket_dir in raw_bucket_dirs.values():
                bucket_dir.mkdir(parents=True, exist_ok=True)

            raw_inventory_rows = []
            for file_path in self.output_dir.rglob("*"):
                if not file_path.is_file():
                    continue
                if file_path.suffix.lower() not in raw_ext:
                    continue
                if self.raw_data_root in file_path.parents:
                    continue
                if self.visualizations_root in file_path.parents:
                    continue

                rel_path = file_path.relative_to(self.output_dir)
                ext = file_path.suffix.lower()
                if ext in {".csv", ".xlsx"}:
                    bucket = "tables"
                elif ext in {".log"}:
                    bucket = "logs"
                elif ext in {".pdb", ".pdbqt", ".sdf", ".mol2"}:
                    bucket = "structures"
                elif ext in {".pml", ".drw", ".sum", ".ps"}:
                    bucket = "tool_files"
                else:
                    bucket = "metadata"

                destination = self._copy_with_collision(
                    file_path,
                    raw_bucket_dirs[bucket] / rel_path,
                )
                raw_inventory_rows.append(
                    {
                        "bucket": bucket,
                        "source_file": str(rel_path),
                        "indexed_copy": str(destination.relative_to(self.output_dir)),
                    }
                )

            manifest_df = pd.DataFrame(
                manifest_rows,
                columns=[
                    "visualization_file",
                    "visualization_type",
                    "source_file",
                    "pdb_code",
                    "protein_display_name",
                ],
            )
            manifest_csv = self.raw_data_root / "visualization_manifest.csv"
            manifest_json = self.raw_data_root / "visualization_manifest.json"
            manifest_df.to_csv(manifest_csv, index=False)
            manifest_df.to_json(manifest_json, orient="records", indent=2)

            raw_inventory_df = pd.DataFrame(
                raw_inventory_rows,
                columns=["bucket", "source_file", "indexed_copy"],
            )
            raw_inventory_csv = self.raw_data_root / "raw_artifact_inventory.csv"
            raw_inventory_json = self.raw_data_root / "raw_artifact_inventory.json"
            raw_inventory_df.to_csv(raw_inventory_csv, index=False)
            raw_inventory_df.to_json(raw_inventory_json, orient="records", indent=2)

            self._write_start_here_file(manifest_csv, raw_inventory_csv)

            self.logger.info(f"✅ Visual assets indexed: {len(manifest_rows)}")
            self.logger.info(f"✅ Raw artifacts indexed: {len(raw_inventory_rows)}")
            self.logger.info(f"✅ Mirrored visual assets during consolidation: {mirrored_assets}")
            self.logger.info(f"✅ Visualization manifest: {manifest_csv}")
            self.logger.info(f"✅ Raw artifact inventory: {raw_inventory_csv}")

            pruned_count = self._prune_empty_directories(self.output_dir)
            if pruned_count:
                self.logger.info(f"🧹 Pruned {pruned_count} empty directories")

            self.results["visualization_manifest"] = str(manifest_csv)
            self.results["raw_artifact_inventory"] = str(raw_inventory_csv)
            return True
        except Exception as exc:
            self.logger.error(f"❌ Artifact organization failed: {exc}", exc_info=True)
            return False
