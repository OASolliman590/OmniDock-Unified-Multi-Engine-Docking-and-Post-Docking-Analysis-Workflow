"""
Validation helpers for exported protein+ligand complex PDB files.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional


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

