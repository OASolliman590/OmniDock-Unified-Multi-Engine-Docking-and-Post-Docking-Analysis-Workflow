"""
Docked ligand microspecies for the best-pose table (Spec 036 R1c / Spec 034 R2d).

Source order, per best-pose row:
1. The pose's Meeko ``REMARK SMILES`` (read with ``atom_mapping.read_pdbqt_lineage_smiles``).
   This is the state that was actually docked.
2. Otherwise the ligand preparation provenance ``preparation_steps/<ligand>.json``, field
   ``protonation.pdbqt_state`` (the measured state of the docking input).
3. Otherwise explicit empty values with ``docked_state_source = unavailable``.

Nothing is re-protonated here. Net charge and the canonical microspecies SMILES are measured
with RDKit, using the same canonical form as ``docking.preparation.ligand_preparation``
(hydrogens removed, stereochemistry removed, charges kept).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from post_docking_analysis.atom_mapping import LineageError, read_pdbqt_lineage_smiles


DOCKED_STATE_COLUMNS = (
    "docked_microspecies_smiles",
    "docked_net_charge",
    "protonation_policy",
    "docked_state_source",
)
SOURCE_POSE_REMARK = "pose_remark_smiles"
SOURCE_PREPARATION = "preparation_provenance"
SOURCE_UNAVAILABLE = "unavailable"

_PREPARATION_STEP_DIRS = (
    "prepared_ligands/preparation_steps",
    "3-Preparation/4-Prepared_Ligand/preparation_steps",
    "ligands/preparation_steps",
)


def _unavailable() -> Dict[str, object]:
    return {
        "docked_microspecies_smiles": "",
        "docked_net_charge": "",
        "protonation_policy": "",
        "docked_state_source": SOURCE_UNAVAILABLE,
    }


def _measure_smiles(smiles: str) -> Optional[Tuple[str, int]]:
    """Return (canonical microspecies SMILES, net formal charge) or None when RDKit cannot read it."""
    try:
        from rdkit import Chem  # type: ignore
        from docking.preparation.ligand_preparation import canonical_microspecies
    except ImportError:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return canonical_microspecies(mol), int(Chem.GetFormalCharge(mol))


def _provenance_path(project_dir: Path, ligand: str) -> Optional[Path]:
    stem = Path(str(ligand or "")).stem
    if not stem:
        return None
    for relative in _PREPARATION_STEP_DIRS:
        candidate = Path(project_dir) / relative / f"{stem}.json"
        if candidate.is_file():
            return candidate
    return None


def _load_provenance(project_dir: Path, ligand: str, cache: Dict[str, Optional[dict]]) -> Optional[dict]:
    key = str(ligand or "")
    if key in cache:
        return cache[key]
    path = _provenance_path(project_dir, key)
    payload: Optional[dict] = None
    if path is not None:
        try:
            loaded = json.loads(path.read_text(encoding="utf-8"))
            payload = loaded if isinstance(loaded, dict) else None
        except (OSError, ValueError):
            payload = None
    cache[key] = payload
    return payload


def _provenance_policy(payload: Optional[dict]) -> str:
    if not payload:
        return ""
    return str(payload.get("protonation_policy") or "").strip()


def _provenance_state(payload: Optional[dict]) -> Optional[Tuple[str, int]]:
    if not payload:
        return None
    block = payload.get("protonation")
    if not isinstance(block, dict):
        return None
    state = block.get("pdbqt_state")
    if not isinstance(state, dict):
        return None
    smiles = str(state.get("microspecies_smiles") or "").strip()
    charge = state.get("net_formal_charge")
    if not smiles or charge is None:
        return None
    try:
        return smiles, int(charge)
    except (TypeError, ValueError):
        return None


def read_docked_microspecies(
    project_dir: Path,
    *,
    ligand: str,
    pose_file: object,
    pose_index: object,
    cache: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    """Return the four docked-state fields for one best-pose row."""
    cache = cache if cache is not None else {}
    provenance_cache: Dict[str, Optional[dict]] = cache.setdefault("provenance", {})
    pose_cache: Dict[Tuple[str, str], Optional[str]] = cache.setdefault("pose_smiles", {})
    provenance = _load_provenance(Path(project_dir), str(ligand or ""), provenance_cache)
    policy = _provenance_policy(provenance)

    pose_text = str(pose_file or "").strip()
    if pose_text.lower().endswith(".pdbqt") and Path(pose_text).is_file():
        try:
            model = int(float(pose_index)) if str(pose_index).strip() not in {"", "nan", "None"} else 1
        except (TypeError, ValueError):
            model = 1
        key = (pose_text, str(model))
        if key not in pose_cache:
            try:
                pose_cache[key] = read_pdbqt_lineage_smiles(Path(pose_text), model)
            except (LineageError, OSError, ValueError):
                pose_cache[key] = None
        remark = pose_cache[key]
        if remark:
            measured = _measure_smiles(remark)
            if measured is not None:
                smiles, charge = measured
                return {
                    "docked_microspecies_smiles": smiles,
                    "docked_net_charge": charge,
                    "protonation_policy": policy,
                    "docked_state_source": SOURCE_POSE_REMARK,
                }

    state = _provenance_state(provenance)
    if state is not None:
        smiles, charge = state
        return {
            "docked_microspecies_smiles": smiles,
            "docked_net_charge": charge,
            "protonation_policy": policy,
            "docked_state_source": SOURCE_PREPARATION,
        }
    return _unavailable()


def annotate_docked_microspecies(frame: pd.DataFrame, project_dir: Optional[Path]) -> pd.DataFrame:
    """Add the four docked-state columns to a best-pose frame (one row per engine/tag)."""
    working = frame.copy()
    if working.empty or project_dir is None:
        for column in DOCKED_STATE_COLUMNS:
            working[column] = "" if column != "docked_net_charge" else None
        if not working.empty:
            working["docked_state_source"] = SOURCE_UNAVAILABLE
        return working
    cache: Dict[str, object] = {}
    records: List[Dict[str, object]] = []
    for _, row in working.iterrows():
        records.append(
            read_docked_microspecies(
                Path(project_dir),
                ligand=str(row.get("ligand", "") or ""),
                pose_file=row.get("pose_file", ""),
                pose_index=row.get("pose", 1),
                cache=cache,
            )
        )
    annotated = pd.DataFrame(records, index=working.index, columns=list(DOCKED_STATE_COLUMNS))
    for column in DOCKED_STATE_COLUMNS:
        working[column] = annotated[column]
    return working
