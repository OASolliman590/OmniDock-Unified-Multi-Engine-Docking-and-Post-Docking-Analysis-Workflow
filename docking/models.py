from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Dict, List, Literal, Optional


DockingEngine = Literal["gnina", "vina", "smina", "autodock4"]
AnalysisMode = Literal["single_engine", "comparative_all_engines", "favorite_engine_continue"]
EnsembleAggregationStrategy = Literal["best_score", "mean_score", "boltzmann_weighted"]
LigandPreparationProfile = Literal[
    "openbabel_only",
    "meeko_only",
    "autodocktools_only",
    "openbabel_meeko",
    "openbabel_meeko_autodock",
    "openbabel_autodocktools",
    "engine_aware_full",
]

LIGAND_PREPARATION_PROFILES: List[str] = [
    "openbabel_only",
    "meeko_only",
    "autodocktools_only",
    "openbabel_meeko",
    "openbabel_meeko_autodock",
    "openbabel_autodocktools",
    "engine_aware_full",
]

VINA_FAMILY_ENGINES = {"gnina", "vina", "smina"}
AD4_ENGINES = {"autodock4"}

DEFAULT_PREPARATION_PH = 7.4
MIN_PREPARATION_PH = 0.0
MAX_PREPARATION_PH = 14.0


@dataclass
class PairlistRow:
    receptor: str
    site_id: str
    ligand: str
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    reference_source: str = ""
    reference_frame_id: str = ""
    receptor_frame_id: str = ""
    reference_pdb_id: str = ""
    reference_pose_file: str = ""
    reference_ligand_file: str = ""

    @property
    def tag(self) -> str:
        return f"{self.receptor}_{self.site_id}_{self.ligand}"

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class SiteRecord:
    pdb_id: str
    selected_ligand: str = ""
    center_x: float = 0.0
    center_y: float = 0.0
    center_z: float = 0.0

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class ProjectManifest:
    project_name: str
    project_root: str
    created_at: str
    asset_mode: str
    engines: List[str]
    pairlist_file: str
    receptors_dir: str
    ligands_dir: str
    metadata_dir: str
    layout_profile: str = "canonical"
    docking_root: str = ""
    post_docking_root: str = ""
    raw_proteins_dir: str = ""
    raw_ligands_dir: str = ""
    raw_ligands_sdf_dir: str = ""
    prepared_proteins_dir: str = ""
    prepared_ligands_dir: str = ""
    enabled_panels: List[str] = field(default_factory=list)
    pairlist_mode: str = ""
    has_cocrystal_benchmark_rows: bool = False
    source_paths: Dict[str, str] = field(default_factory=dict)
    favorite_engine: str = ""
    engine_settings: Dict[str, Dict[str, object]] = field(default_factory=dict)
    pair_curation_state_file: str = ""
    latest_pair_round: str = ""
    deployment_root: str = ""
    latest_rerun_manifest_file: str = ""
    top_pose_selection_policy: str = "best_affinity"
    top_pose_global_aggregation: str = "best_target"
    hpc_profile: Dict[str, str] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class EngineJobResult:
    tag: str
    command: List[str]
    pose_file: str
    log_file: str
    status: str
    returncode: Optional[int] = None
    error: str = ""
    engine: str = ""
    fingerprint: str = ""
    completion_file: str = ""
    input_files: List[Dict[str, str]] = field(default_factory=list)
    executables: List[str] = field(default_factory=list)
    skip_completed: bool = False

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class LigandPreparationCompatibility:
    requested_profile: str
    effective_profile: str
    selected_engines: List[str]
    is_valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class EnsembleReceptorConfig:
    """
    Extension-point model for future ensemble receptor docking support.

    Notes
    -----
    Current implementation is a contract scaffold only. Engine runners may
    inspect this payload in runtime config but are not required to execute
    multi-conformation docking yet.
    """

    enabled: bool = False
    strategy: EnsembleAggregationStrategy = "best_score"
    conformer_paths: List[str] = field(default_factory=list)
    max_conformers: int = 0
    selection_note: str = ""

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


def normalize_ensemble_aggregation_strategy(raw: object) -> str:
    token = str(raw or "").strip().lower()
    if token in {"best", "best_score", "min"}:
        return "best_score"
    if token in {"mean", "average", "mean_score"}:
        return "mean_score"
    if token in {"boltzmann", "boltzmann_weighted"}:
        return "boltzmann_weighted"
    return "best_score"


def normalize_engine_names(engines: Optional[List[str]]) -> List[str]:
    normalized: List[str] = []
    for engine in engines or []:
        token = str(engine or "").strip().lower()
        if token and token not in normalized:
            normalized.append(token)
    return normalized


def normalize_ligand_preparation_profile(raw: str) -> str:
    token = str(raw or "").strip().lower()
    alias_map = {
        "": "engine_aware_full",
        "auto": "engine_aware_full",
        "default": "engine_aware_full",
        "obabel": "openbabel_only",
        "openbabel": "openbabel_only",
        "meeko": "meeko_only",
        "autodocktools": "autodocktools_only",
        "adt": "autodocktools_only",
    }
    normalized = alias_map.get(token, token)
    if normalized in LIGAND_PREPARATION_PROFILES:
        return normalized
    return "engine_aware_full"


def resolve_effective_ligand_preparation_profile(profile: str, selected_engines: Optional[List[str]] = None) -> str:
    requested = normalize_ligand_preparation_profile(profile)
    engines = normalize_engine_names(selected_engines)
    if requested != "engine_aware_full":
        return requested
    if "autodock4" in engines:
        return "openbabel_meeko_autodock"
    return "openbabel_meeko"


def validate_ligand_preparation_profile(
    profile: str,
    selected_engines: Optional[List[str]] = None,
) -> LigandPreparationCompatibility:
    engines = normalize_engine_names(selected_engines)
    requested_profile = normalize_ligand_preparation_profile(profile)
    effective_profile = resolve_effective_ligand_preparation_profile(requested_profile, engines)
    errors: List[str] = []
    warnings: List[str] = []

    unknown = [engine for engine in engines if engine not in (VINA_FAMILY_ENGINES | AD4_ENGINES)]
    if unknown:
        warnings.append(
            "Unknown engine names were ignored for preparation compatibility: " + ", ".join(sorted(unknown))
        )
        engines = [engine for engine in engines if engine not in unknown]

    includes_ad4 = "autodock4" in engines
    if includes_ad4 and effective_profile in {"openbabel_only", "meeko_only", "openbabel_meeko"}:
        errors.append(
            "AutoDock4 selected, but the ligand preparation profile does not include AutoDockTools compatibility. "
            "Use engine_aware_full, openbabel_meeko_autodock, openbabel_autodocktools, or autodocktools_only."
        )

    if not engines:
        warnings.append(
            "No engines were selected; engine_aware_full resolves to Open Babel -> Meeko by default."
        )

    if "autodock4" not in engines and requested_profile == "autodocktools_only":
        warnings.append(
            "AutoDockTools-only preparation selected without AutoDock4 engine. "
            "This may be unnecessary for Vina/Smina/GNINA-only runs."
        )

    return LigandPreparationCompatibility(
        requested_profile=requested_profile,
        effective_profile=effective_profile,
        selected_engines=engines,
        is_valid=not errors,
        errors=errors,
        warnings=warnings,
    )


def validate_preparation_ph(raw: object) -> tuple[bool, float, str]:
    """
    Validate preparation protonation pH for Open Babel enrichment stages.

    Returns
    -------
    tuple[bool, float, str]
        (is_valid, normalized_ph, error_message)
    """
    try:
        value = float(raw)
    except Exception:
        return False, float(DEFAULT_PREPARATION_PH), "pH must be numeric."

    if value < float(MIN_PREPARATION_PH) or value > float(MAX_PREPARATION_PH):
        return (
            False,
            float(DEFAULT_PREPARATION_PH),
            f"pH must be between {MIN_PREPARATION_PH:.1f} and {MAX_PREPARATION_PH:.1f}.",
        )
    return True, float(value), ""
