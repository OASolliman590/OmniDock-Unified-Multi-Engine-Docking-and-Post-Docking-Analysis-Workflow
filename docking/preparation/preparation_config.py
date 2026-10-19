"""Preparation configuration for the AutoDock preparation pipeline (Spec 038 R2).

Moved verbatim from autodock_preparation.py; autodock_preparation re-exports these names.
"""

from dataclasses import dataclass
from typing import List, Optional

from docking.models import resolve_effective_ligand_preparation_profile


def profile_requires_autodocktools(profile: str, selected_engines: Optional[List[str]] = None) -> bool:
    effective = resolve_effective_ligand_preparation_profile(profile, selected_engines)
    if effective in {"autodocktools_only", "openbabel_autodocktools", "openbabel_meeko_autodock"}:
        return True
    return False


@dataclass
class PreparationConfig:
    """Configuration for AutoDock preparation"""
    ligands_input: str
    receptors_input: str
    ligands_output: str
    receptors_output: str
    # Spec 034 R3: no silent pH or force field. Each is an explicit value with its source
    # ("user_entered" from CLI/prompt, "config_file" from an explicit configuration value).
    force_field: Optional[str] = None
    force_field_source: Optional[str] = None
    ph: Optional[float] = None
    ph_source: Optional[str] = None
    # Spec 034 R2c: per-ligand protonation policy and the explicit state map (CSV) when required.
    protonation_policy: str = "ph_model"
    protonation_state_map: str = ""
    # Strict receptor preparation refuses permissive bad-residue deletion (Spec 033).
    allow_bad_res: bool = False
    default_altloc: str = "A"
    # Spec 033 R2c (Scientific Lead, option A): PDB2PQR protonation at pH/force field.
    receptor_use_pdb2pqr: bool = True
    validate_outputs: bool = True
    min_file_size_kb: int = 1
    ligand_preparation_backend: str = "engine_aware_full"  # backward-compatible field
    ligand_preparation_profile: str = "engine_aware_full"
    selected_engines: List[str] = None
    autodocktools_prepare_ligand4: str = ""
    autodocktools_prepare_receptor4: str = ""
    autodocktools_python: str = ""
    # Spec 037 R1: the project whose layered protonation policy (receptor and ligand overrides) applies.
    project_dir: str = ""
    # Spec 037 R4a: prepare every entity again, even when its existing output is reusable.
    force: bool = False
