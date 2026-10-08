"""Approved Spec 031 binding-site center contract."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Dict, Optional, Sequence

import numpy as np

try:
    import Bio
    from Bio.PDB import MMCIFParser, PDBParser
except ImportError:  # pragma: no cover
    Bio = None
    MMCIFParser = PDBParser = None


CENTER_METHOD = "binding_site_center_v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _element(atom) -> str:
    value = str(getattr(atom, "element", "") or "").strip().upper()
    if value:
        return value
    return "".join(char for char in str(atom.get_name()) if char.isalpha())[:1].upper()


def _altloc_key(atom) -> tuple:
    occupancy = atom.get_occupancy()
    occupancy_value = float(occupancy) if occupancy is not None else 0.0
    altloc = str(atom.get_altloc() or "").strip()
    return (-occupancy_value, 0 if altloc == "A" else 1, altloc)


def _selected_residue_atoms(residue) -> list:
    variants: Dict[str, list] = {}
    for atom in residue.get_unpacked_list():
        variants.setdefault(str(atom.get_name()), []).append(atom)
    return [sorted(atoms, key=_altloc_key)[0] for _, atoms in sorted(variants.items())]


def resolve_binding_site_center(
    *,
    structure_file: Optional[Path] = None,
    accession: str = "",
    ligand_name: str = "",
    chain_id: str = "",
    residue_number: Optional[int] = None,
    insertion_code: str = "",
    model_number: Optional[int] = None,
    explicit_center: Optional[Sequence[float]] = None,
    explicit_center_source: str = "",
) -> Dict[str, object]:
    """Resolve an explicit user center or selected holo-ligand heavy-atom centroid."""
    if explicit_center is not None:
        center = np.asarray(list(explicit_center), dtype=float)
        if center.shape != (3,) or not np.isfinite(center).all():
            return {"status": "failed", "reason": "invalid_explicit_center", "center_method": CENTER_METHOD}
        return {
            "status": "completed",
            "reason": "explicit_user_coordinates",
            "binding_site_center": center,
            "ligand_centroid": None,
            "center_method": CENTER_METHOD,
            "center_units": "angstrom",
            "coordinate_frame": "source_structure_cartesian",
            "explicit_center_source": str(explicit_center_source),
            "compatibility_aliases_written": True,
            "overall_center": center,
            "ligand_center": center,
        }
    if structure_file is None or not ligand_name or not chain_id or residue_number is None:
        return {
            "status": "skipped_missing_configuration",
            "reason": "explicit_ligand_selection_or_user_center_required",
            "center_method": CENTER_METHOD,
        }
    if Bio is None:
        return {"status": "skipped_missing_dependency", "reason": "biopython", "center_method": CENTER_METHOD}
    source = Path(structure_file).expanduser().resolve()
    if not source.is_file():
        return {"status": "failed", "reason": "missing_structure_file", "center_method": CENTER_METHOD}
    parser = MMCIFParser(QUIET=True) if source.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    structure = parser.get_structure(str(accession or source.stem), str(source))
    models = list(structure.get_models())
    compatibility_default = False
    if model_number is None:
        if len(models) == 1:
            selected_model = models[0]
            selected_model_number = 1
        else:
            selected_model = models[0]
            selected_model_number = 1
            compatibility_default = True
    else:
        requested = int(model_number)
        matches = [model for ordinal, model in enumerate(models, start=1) if ordinal == requested]
        if not matches:
            return {"status": "failed", "reason": "model_not_found", "center_method": CENTER_METHOD}
        selected_model = matches[0]
        selected_model_number = requested
    if chain_id not in selected_model:
        return {"status": "failed", "reason": "chain_not_found", "center_method": CENTER_METHOD}
    icode = str(insertion_code or "").strip()
    matches = [
        residue for residue in selected_model[chain_id]
        if residue.get_resname().strip().upper() == str(ligand_name).strip().upper()
        and int(residue.get_id()[1]) == int(residue_number)
        and str(residue.get_id()[2] or "").strip() == icode
    ]
    if len(matches) != 1:
        return {
            "status": "failed",
            "reason": "selected_ligand_not_found" if not matches else "selected_ligand_ambiguous",
            "center_method": CENTER_METHOD,
        }
    atoms = [atom for atom in _selected_residue_atoms(matches[0]) if _element(atom) != "H"]
    if not atoms:
        return {"status": "failed", "reason": "selected_ligand_has_no_heavy_atoms", "center_method": CENTER_METHOD}
    center = np.mean(np.asarray([atom.get_coord() for atom in atoms], dtype=float), axis=0)
    return {
        "status": "completed",
        "reason": "selected_holo_ligand_heavy_atom_centroid",
        "binding_site_center": center,
        "ligand_centroid": center,
        "overall_center": center,
        "ligand_center": center,
        "contact_residue_centroid": None,
        "center_method": CENTER_METHOD,
        "center_units": "angstrom",
        "coordinate_frame": "source_structure_cartesian",
        "source_structure_path": source.as_posix(),
        "source_structure_sha256": _sha256(source),
        "accession": str(accession or source.stem),
        "model_number": selected_model_number,
        "model_selection_reason": "compatibility_default_model_1" if compatibility_default else "explicit_or_single_model",
        "chain_id": str(chain_id),
        "residue_name": str(ligand_name).strip().upper(),
        "residue_number": int(residue_number),
        "insertion_code": icode,
        "heavy_atom_count": len(atoms),
        "hydrogen_policy": "excluded",
        "altloc_policy": "highest_occupancy_then_A_then_lexical",
        "parser": "biopython",
        "parser_version": str(getattr(Bio, "__version__", "")),
        "compatibility_aliases_written": True,
    }
