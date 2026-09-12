"""Exact pose records and graph-aware, receptor-frame docking RMSD.

No atom-count truncation, coordinate fitting, inferred bonds, or record fallback.
RDKit is required for graph identity and symmetry. Coordinate-only files are
unsupported unless Meeko's authoritative SMILES/atom-index mapping is present.
"""
from __future__ import annotations
import hashlib
import re
from pathlib import Path
import numpy as np

class PoseGeometryError(ValueError):
    pass

def content_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def selected_record(path: Path, pose_index: int = 1) -> str:
    try:
        index = int(pose_index)
    except (ValueError, TypeError):
        raise PoseGeometryError("invalid_pose_index")
    if index < 1 or float(pose_index) != index:
        raise PoseGeometryError("invalid_pose_index")
    text = Path(path).read_text(encoding="utf-8", errors="strict")
    if Path(path).suffix.lower() in {".sdf", ".mol"}:
        records = []
        for position, block in enumerate(text.split("$$$$")):
            if not block.strip():
                continue
            # Remove only the record-separator newline; an empty molecule title
            # is a legitimate first header line and must remain present.
            if position and block.startswith("\r\n"):
                block = block[2:]
            elif position and block.startswith("\n"):
                block = block[1:]
            records.append(block)
    elif re.search(r"^MODEL\s", text, flags=re.M):
        records = []
        current = None
        for line in text.splitlines():
            if line.startswith("MODEL"):
                if current is not None:
                    raise PoseGeometryError("unterminated_model")
                current = []
            elif line.startswith("ENDMDL"):
                if current is not None:
                    records.append("\n".join(current))
                    current = None
            elif current is not None:
                current.append(line)
        if current is not None:
            raise PoseGeometryError("unterminated_model")
    else:
        records = [text]
    if index > len(records):
        raise PoseGeometryError("pose_index_out_of_range")
    return records[index - 1]

def load_pose_molecule(path: Path, pose_index: int = 1):
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise PoseGeometryError("rdkit_required_for_chemical_mapping") from exc
    path = Path(path)
    record = selected_record(path, pose_index)
    if path.suffix.lower() in {".sdf", ".mol"}:
        molecule = Chem.MolFromMolBlock(record, sanitize=True, removeHs=False, strictParsing=True)
        if molecule is None:
            raise PoseGeometryError("invalid_molecule_graph")
    elif path.suffix.lower() == ".pdbqt":
        # Meeko's atom mapping numbers are 1-based SMILES index / PDBQT serial.
        header = path.read_text(encoding="utf-8")
        smiles_match = re.search(r"^REMARK SMILES (?!IDX|H PARENT)(\S+)", record, re.M)
        if not smiles_match:
            smiles_match = re.search(r"^REMARK SMILES (?!IDX|H PARENT)(\S+)", header, re.M)
        mapping_lines = re.findall(r"^REMARK SMILES IDX (.+)$", record, re.M)
        if not mapping_lines:
            mapping_lines = re.findall(r"^REMARK SMILES IDX (.+)$", header.split("MODEL", 1)[0], re.M)
        if not smiles_match or not mapping_lines:
            raise PoseGeometryError("missing_authoritative_atom_mapping")
        molecule = Chem.MolFromSmiles(smiles_match.group(1))
        if molecule is None:
            raise PoseGeometryError("invalid_smiles")
        atoms, atom_types = {}, {}
        for line in record.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                serial = int(line[6:11])
                if serial in atoms:
                    raise PoseGeometryError("duplicate_atom_serial")
                atoms[serial] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                atom_types[serial] = line.split()[-1]
        numbers = [int(value) for line in mapping_lines for value in line.split()]
        if len(numbers) % 2:
            raise PoseGeometryError("invalid_atom_mapping")
        mapping = dict(zip(numbers[::2], numbers[1::2]))
        if len(mapping) != len(numbers) // 2 or len(set(mapping.values())) != len(mapping):
            raise PoseGeometryError("ambiguous_atom_mapping")
        if set(mapping) != set(range(1, molecule.GetNumAtoms() + 1)):
            raise PoseGeometryError("incomplete_atom_mapping")
        # Meeko closure pseudoatoms have no chemical counterpart in the original
        # SMILES. Their paired CG atoms are ordinary mapped carbons.
        ignored_types = {"H", "HD", "HS", "G0", "G1", "G2", "G3"}
        heavy_serials = {serial for serial, atom_type in atom_types.items() if atom_type not in ignored_types}
        if heavy_serials != set(mapping.values()):
            raise PoseGeometryError("unmapped_heavy_atoms")
        conformer = Chem.Conformer(molecule.GetNumAtoms())
        for index, serial in mapping.items():
            if serial not in atoms:
                raise PoseGeometryError("mapped_atom_missing")
            atomic_symbol = {"A":"C", "NA":"N", "OA":"O", "SA":"S", "CG0":"C", "CG1":"C", "CG2":"C", "CG3":"C"}.get(atom_types[serial], atom_types[serial])
            if atomic_symbol.lower() != molecule.GetAtomWithIdx(index - 1).GetSymbol().lower():
                raise PoseGeometryError("mapped_element_mismatch")
            conformer.SetAtomPosition(index - 1, atoms[serial])
        molecule.AddConformer(conformer)
    else:
        raise PoseGeometryError("coordinate_only_format_requires_chemical_template")
    molecule = Chem.RemoveHs(molecule)
    if molecule.GetNumAtoms() == 0 or molecule.GetNumConformers() != 1:
        raise PoseGeometryError("missing_pose_coordinates")
    if not np.isfinite(molecule.GetConformer().GetPositions()).all():
        raise PoseGeometryError("nonfinite_coordinates")
    return molecule

def fixed_frame_rmsd(molecule_a, molecule_b) -> float:
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    if Chem.MolToSmiles(molecule_a, isomericSmiles=True) != Chem.MolToSmiles(molecule_b, isomericSmiles=True):
        raise PoseGeometryError("chemical_identity_mismatch")
    matches = molecule_b.GetSubstructMatches(molecule_a, uniquify=False, useChirality=True, maxMatches=100001)
    if not matches or len(matches) > 100000:
        raise PoseGeometryError("ambiguous_or_excessive_atom_mappings")
    mappings = [list(enumerate(match)) for match in matches]
    # Explicit stereochemistry-aware graph mappings; CalcRMS never fits coordinates.
    return float(rdMolAlign.CalcRMS(molecule_a, molecule_b, map=mappings))

def pose_rmsd(path_a, path_b, *, pose_a=1, pose_b=1) -> float:
    return fixed_frame_rmsd(load_pose_molecule(path_a, pose_a), load_pose_molecule(path_b, pose_b))
