from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Literal, Optional


WorkflowTarget = Literal[
    "workflow.init",
    "workflow.init_project",
    "workflow.clone_checkpoint",
    "workflow.clone_maturation",
    "pdb.collect",
    "pdb.collect.append",
    "pdb.fetch",
    "pdb.run",
    "pdb.batch",
    "pdb.prepare_protein",
    "pdb.prepare_ligand",
    "pdb.prepare_both",
    "prep.pairlist",
    "prep.project",
    "dock.run",
    "dock.deploy",
    "dock.sync",
    "dock.submit",
    "dock.engine.gnina",
    "dock.engine.vina",
    "dock.engine.smina",
    "dock.engine.autodock4",
    "analyze.comparative",
    "analyze.favorite_engine",
    "analyze.stage.hierarchical",
    "analyze.stage.polypharmacology",
    "analyze.stage.rmsd",
    "analyze.stage.reports",
    "analyze.stage.visualizations",
    "analyze.stage.structure_quality",
    "analyze.interactions.pandamap",
    "analyze.interactions.prolif",
    "analyze.interactions.ligplot",
    "analyze.interactions.poseview",
    "analyze.interactions.clean",
    "analyze.visuals.py3dmol",
    "analyze.visuals.pymol",
]

WorkflowStepStatus = Literal[
    "not_started",
    "in_progress",
    "completed",
    "validated",
    "failed",
    "skipped",
    "needs_review",
    "partial",
    "blocked",
]


BackgroundTaskStatus = Literal[
    "queued",
    "running",
    "completed",
    "failed",
    "cancelled",
]


@dataclass
class WorkflowContext:
    project_root: str = ""
    analysis_output_dir: str = ""
    favorite_engine: str = ""
    selected_engines: List[str] = field(default_factory=list)
    enabled_panels: List[str] = field(default_factory=list)
    pair_mode: str = ""
    active_target: str = ""
    layout_profile: str = ""
    raw_proteins_dir: str = ""
    raw_ligands_dir: str = ""
    raw_ligands_sdf_dir: str = ""
    prepared_proteins_dir: str = ""
    prepared_ligands_dir: str = ""
    docking_root: str = ""
    post_docking_root: str = ""

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class WorkflowStepRecord:
    status: WorkflowStepStatus
    started_at: str
    finished_at: str
    inputs: Dict[str, object] = field(default_factory=dict)
    outputs: Dict[str, object] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class WorkflowStepResult:
    target: str
    status: WorkflowStepStatus
    root: str
    inputs: Dict[str, object] = field(default_factory=dict)
    outputs: Dict[str, object] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class BackgroundTaskRecord:
    task_id: str
    label: str
    target: str
    status: BackgroundTaskStatus
    created_at: str
    updated_at: str
    root: str
    progress: int = 0
    note: str = ""
    result_status: str = ""
    error: str = ""

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class DockForgeFeatureFlags:
    """Feature toggles used to stage new workflow capabilities safely."""

    enable_timeline: bool = True
    enable_background_tasks: bool = True
    enable_checkpoint_revise: bool = True
    enable_execution_environment_adapters: bool = True
    enable_post_docking_first_class: bool = True
    enable_sqlite_dual_write: bool = False

    def to_dict(self) -> Dict[str, bool]:
        return {
            "enable_timeline": bool(self.enable_timeline),
            "enable_background_tasks": bool(self.enable_background_tasks),
            "enable_checkpoint_revise": bool(self.enable_checkpoint_revise),
            "enable_execution_environment_adapters": bool(self.enable_execution_environment_adapters),
            "enable_post_docking_first_class": bool(self.enable_post_docking_first_class),
            "enable_sqlite_dual_write": bool(self.enable_sqlite_dual_write),
        }


@dataclass
class CheckpointMetadataRecord:
    """Lineage metadata for Checkpoint & Revise workspace snapshots."""

    checkpoint_id: str
    source_project_dir: str
    target_project_dir: str
    created_at: str
    layout_profile: str
    marker_file: str = ""
    status: str = "created"
    note: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)
