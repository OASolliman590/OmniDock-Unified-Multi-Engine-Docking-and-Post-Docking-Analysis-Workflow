"""
Validation helpers for exported protein+ligand complex PDB files.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple


STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "ASH", "CYM", "CYX", "GLH", "HID", "HIE", "HIP", "LYN", "MSE",
}


def validate_complex_pdb_structure(
    complex_pdb: Path,
    *,
    receptor_reference: Optional[Path] = None,
) -> Dict[str, object]:
    """
    Validate a protein+pose complex PDB export for downstream rescoring/visualization.

    Validation checks are intentionally lightweight and format-oriented:
    - file exists and is non-empty
    - contains parseable ATOM/HETATM coordinate records
    - contains both receptor atoms (ATOM) and ligand atoms (HETATM)
    - atom serials are unique
    - PDB terminator (`END`) present (warn-only if missing)
    - optional receptor atom-count mismatch warning against the reference receptor
    """
    pdb_path = Path(complex_pdb).expanduser().resolve()
    errors: List[str] = []
    warnings: List[str] = []

    result: Dict[str, object] = {
        "complex_pdb": str(pdb_path),
        "is_valid": False,
        "errors": errors,
        "warnings": warnings,
        "atom_count": 0,
        "receptor_atom_count": 0,
        "ligand_atom_count": 0,
        "has_end_record": False,
    }

    if not pdb_path.exists():
        errors.append("complex file is missing")
        return result
    if not pdb_path.is_file():
        errors.append("complex path is not a file")
        return result
    if pdb_path.stat().st_size <= 0:
        errors.append("complex file is empty")
        return result

    lines = pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    atom_lines: List[str] = []
    receptor_atom_count = 0
    ligand_atom_count = 0
    serials: List[int] = []
    parse_failures = 0

    for line in lines:
        if line.startswith("ATOM"):
            atom_lines.append(line)
            receptor_atom_count += 1
        elif line.startswith("HETATM"):
            atom_lines.append(line)
            ligand_atom_count += 1

    for line in atom_lines:
        try:
            serial = int(line[6:11].strip())
            serials.append(serial)
        except Exception:
            parse_failures += 1
        try:
            float(line[30:38].strip())
            float(line[38:46].strip())
            float(line[46:54].strip())
        except Exception:
            parse_failures += 1

    has_end = any(line.strip() == "END" for line in lines)
    if not has_end:
        warnings.append("complex file missing END record")

    if not atom_lines:
        errors.append("complex has no ATOM/HETATM coordinate records")
    if receptor_atom_count == 0:
        errors.append("complex has no receptor ATOM records")
    if ligand_atom_count == 0:
        errors.append("complex has no ligand HETATM records")
    if parse_failures > 0:
        errors.append(f"complex contains {parse_failures} malformed atom/coordinate fields")

    if serials and len(serials) != len(set(serials)):
        errors.append("complex contains duplicate atom serial numbers")

    if receptor_reference is not None:
        receptor_path = Path(receptor_reference).expanduser().resolve()
        if receptor_path.exists() and receptor_path.is_file():
            ref_lines = receptor_path.read_text(encoding="utf-8", errors="ignore").splitlines()
            ref_atoms = sum(1 for line in ref_lines if line.startswith("ATOM"))
            if ref_atoms > 0 and ref_atoms != receptor_atom_count:
                warnings.append(
                    f"receptor atom count mismatch vs reference ({receptor_atom_count} != {ref_atoms})"
                )

    result.update(
        {
            "is_valid": not errors,
            "atom_count": len(atom_lines),
            "receptor_atom_count": receptor_atom_count,
            "ligand_atom_count": ligand_atom_count,
            "has_end_record": has_end,
        }
    )
    return result


def _strict_pdb_atom(line: str) -> Tuple[int, str, str, str, str, int, float, float, float, str]:
    """Parse the fixed columns consumed by the MD-export profile."""
    if len(line) < 78:
        raise ValueError("atom_record_shorter_than_78_columns")
    record = line[0:6].strip()
    if record not in {"ATOM", "HETATM"}:
        raise ValueError("unsupported_atom_record")
    serial = int(line[6:11].strip())
    atom_name = line[12:16].strip()
    altloc = line[16]
    resname = line[17:20].strip()
    chain = line[21]
    resseq = int(line[22:26].strip())
    x = float(line[30:38].strip())
    y = float(line[38:46].strip())
    z = float(line[46:54].strip())
    if not all(math.isfinite(value) for value in (x, y, z)):
        raise ValueError("non_finite_coordinates")
    element = line[76:78].strip()
    if not atom_name:
        raise ValueError("missing_atom_name")
    if not resname:
        raise ValueError("missing_residue_name")
    if chain == " ":
        raise ValueError("missing_chain_id_column_22")
    if not element:
        raise ValueError("missing_element_columns_77_78")
    return serial, atom_name, altloc, resname, chain, resseq, x, y, z, element


def validate_md_receptor_pdb(receptor_pdb: Path) -> Dict[str, object]:
    """Validate an already-prepared, protein-only receptor without cleaning it."""
    path = Path(receptor_pdb).expanduser().resolve()
    errors: List[str] = []
    atom_count = 0
    amino_acid_atom_count = 0
    amino_acid_residues: Set[Tuple[str, int, str]] = set()
    serials: Set[int] = set()
    if not path.is_file():
        return {
            "status": "skipped_missing_configuration",
            "is_valid": False,
            "errors": [f"receptor_file_missing:{path}"],
            "atom_count": 0,
            "amino_acid_atom_count": 0,
            "amino_acid_residue_count": 0,
        }
    if path.suffix.lower() != ".pdb":
        return {
            "status": "not_comparable",
            "is_valid": False,
            "errors": [f"consumer_profile_v1_requires_prepared_pdb:{path.suffix.lower()}"],
            "atom_count": 0,
            "amino_acid_atom_count": 0,
            "amino_acid_residue_count": 0,
        }
    for line_number, line in enumerate(path.read_text(encoding="utf-8", errors="strict").splitlines(), start=1):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_count += 1
        try:
            serial, _, altloc, resname, chain, resseq, *_ = _strict_pdb_atom(line)
        except Exception as exc:
            errors.append(f"line_{line_number}:{exc}")
            continue
        if serial in serials:
            errors.append(f"line_{line_number}:duplicate_atom_serial:{serial}")
        serials.add(serial)
        if altloc != " ":
            errors.append(f"line_{line_number}:unresolved_altloc:{altloc}")
        if line.startswith("ATOM  ") and resname.upper() in STANDARD_AMINO_ACIDS:
            amino_acid_atom_count += 1
            amino_acid_residues.add((chain, resseq, resname.upper()))
    if atom_count == 0:
        errors.append("receptor_has_no_atoms")
    if amino_acid_atom_count == 0:
        errors.append("protein_profile_has_no_recognized_amino_acid_atoms")
    return {
        "status": "completed" if not errors else "failed",
        "is_valid": not errors,
        "errors": errors,
        "atom_count": atom_count,
        "amino_acid_atom_count": amino_acid_atom_count,
        "amino_acid_residue_count": len(amino_acid_residues),
        "polymer_class": "protein" if amino_acid_atom_count else "not_protein",
    }


def validate_md_system_pdb(
    system_pdb: Path,
    *,
    expected_receptor_atoms: int,
    expected_ligand_atoms: int,
    ligand_chain: str,
    ligand_resname: str = "LIG",
    expected_ligand_bonds: Optional[Sequence[Tuple[int, int]]] = None,
) -> Dict[str, object]:
    """Strict fixed-column/count/name/connectivity validation for a combined system."""
    path = Path(system_pdb).expanduser().resolve()
    errors: List[str] = []
    serials: Set[int] = set()
    ligand_serials: List[int] = []
    ligand_names: List[str] = []
    receptor_atoms = 0
    ligand_atoms = 0
    conect_edges: Set[Tuple[int, int]] = set()
    lines = path.read_text(encoding="utf-8", errors="strict").splitlines() if path.is_file() else []
    if not lines:
        errors.append("system_pdb_missing_or_empty")
    for line_number, line in enumerate(lines, start=1):
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                serial, atom_name, altloc, resname, chain, _, *_rest, element = _strict_pdb_atom(line)
            except Exception as exc:
                errors.append(f"line_{line_number}:{exc}")
                continue
            if serial in serials:
                errors.append(f"line_{line_number}:duplicate_atom_serial:{serial}")
            serials.add(serial)
            if altloc != " ":
                errors.append(f"line_{line_number}:nonblank_altloc:{altloc}")
            if chain == ligand_chain and resname == ligand_resname:
                ligand_atoms += 1
                ligand_serials.append(serial)
                ligand_names.append(atom_name)
                expected_element = "".join(ch for ch in atom_name if ch.isalpha()).upper()
                if expected_element and not expected_element.startswith(element.upper()):
                    errors.append(f"line_{line_number}:atom_name_element_mismatch:{atom_name}:{element}")
            else:
                receptor_atoms += 1
        elif line.startswith("CONECT"):
            try:
                values = [int(line[index : index + 5].strip()) for index in range(6, len(line), 5) if line[index : index + 5].strip()]
            except ValueError:
                errors.append(f"line_{line_number}:invalid_conect")
                continue
            if values:
                for target in values[1:]:
                    conect_edges.add(tuple(sorted((values[0], target))))
    if receptor_atoms != int(expected_receptor_atoms):
        errors.append(f"receptor_atom_count_mismatch:{receptor_atoms}!={expected_receptor_atoms}")
    if ligand_atoms != int(expected_ligand_atoms):
        errors.append(f"ligand_atom_count_mismatch:{ligand_atoms}!={expected_ligand_atoms}")
    if len(ligand_names) != len(set(ligand_names)):
        errors.append("duplicate_ligand_atom_names")
    if not any(line.startswith("TER") for line in lines):
        errors.append("missing_TER_record")
    if not any(line.strip() == "END" for line in lines):
        errors.append("missing_END_record")

    if expected_ligand_bonds is not None and len(ligand_serials) == int(expected_ligand_atoms):
        serial_by_index = {index + 1: serial for index, serial in enumerate(ligand_serials)}
        required_edges = {
            tuple(sorted((serial_by_index[int(left)], serial_by_index[int(right)])))
            for left, right in expected_ligand_bonds
            if int(left) in serial_by_index and int(right) in serial_by_index
        }
        missing = sorted(required_edges - conect_edges)
        if missing:
            errors.append(f"conect_missing_bonds:{len(missing)}")
    return {
        "status": "completed" if not errors else "failed",
        "is_valid": not errors,
        "errors": errors,
        "atom_count": receptor_atoms + ligand_atoms,
        "receptor_atom_count": receptor_atoms,
        "ligand_atom_count": ligand_atoms,
        "ligand_unique_name_count": len(set(ligand_names)),
        "conect_edge_count": len(conect_edges),
    }

