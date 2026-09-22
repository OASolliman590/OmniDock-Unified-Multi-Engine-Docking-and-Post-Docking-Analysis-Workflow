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


_MODEL_LINE = re.compile(r"^MODEL(?:\s|$)")
_END_MODEL_LINE = re.compile(r"^ENDMDL(?:\s|$)")


def _model_identifier(line: str) -> int:
    """Return the positive one-based MODEL identifier used as pose ordinal."""
    tokens = line.split()
    if len(tokens) != 2:
        raise PoseGeometryError("invalid_model_identifier")
    try:
        identifier = int(tokens[1])
    except (TypeError, ValueError):
        raise PoseGeometryError("invalid_model_identifier")
    if identifier < 1:
        raise PoseGeometryError("invalid_model_identifier")
    return identifier


def _split_sdf_records(text: str) -> list[str]:
    """Split SDF records while retaining empty and malformed record ordinals.

    A terminal separator's whitespace tail is not a record. Every block before
    a separator, including an empty block, remains addressable by its 1-based
    pose index. An unterminated non-empty tail is retained so a requested pose
    is reported as malformed rather than silently re-numbered.
    """
    lines = text.splitlines(keepends=True)
    separator_positions = [
        position for position, line in enumerate(lines) if line.strip() == "$$$$"
    ]
    if not separator_positions:
        return [text]

    records: list[str] = []
    start = 0
    for position in separator_positions:
        records.append("".join(lines[start:position]))
        start = position + 1
    trailing = "".join(lines[start:])
    if trailing.strip():
        records.append(trailing)
    return records


def _split_pdbqt_models(text: str) -> tuple[list[str], str, bool]:
    """Return model records and the pre-model global header.

    The returned model records contain only lines between their own MODEL and
    ENDMDL boundaries. Positive one-based MODEL identifiers must be sequential
    because callers address records by the same ordinal. Once a model has
    started, non-whitespace content outside a model is rejected so metadata from
    a later/trailing scope cannot be mistaken for global metadata.
    """
    lines = text.splitlines()
    has_models = any(_MODEL_LINE.match(line) for line in lines)
    if not has_models:
        return [text], text, False

    records: list[str] = []
    global_lines: list[str] = []
    current: list[str] | None = None
    started_models = False
    for line in lines:
        if _MODEL_LINE.match(line):
            if current is not None:
                raise PoseGeometryError("nested_model")
            identifier = _model_identifier(line)
            expected_identifier = len(records) + 1
            if identifier != expected_identifier:
                raise PoseGeometryError("nonsequential_model_identifier")
            current = []
            started_models = True
            continue
        if _END_MODEL_LINE.match(line):
            if current is None:
                raise PoseGeometryError("orphan_endmdl")
            records.append("\n".join(current))
            current = None
            continue
        if current is not None:
            current.append(line)
        elif not started_models:
            if line.strip() and not line.lstrip().startswith("REMARK"):
                raise PoseGeometryError("pdbqt_content_outside_model")
            global_lines.append(line)
        elif line.strip() and line.strip() != "END":
            raise PoseGeometryError("pdbqt_content_outside_model")

    if current is not None:
        raise PoseGeometryError("unterminated_model")
    if not records:
        raise PoseGeometryError("no_models")
    return records, "\n".join(global_lines), True


def _pdbqt_metadata(scope: str) -> tuple[bool, str | None, list[int], str | None]:
    """Parse one metadata scope without borrowing from another scope.

    ``REMARK SMILES IDX`` can be wrapped over multiple lines by Meeko, so all
    such lines are concatenated. Repeated SMILES declarations and malformed or
    partial declarations are errors; a missing pair is reported by the caller.
    """
    metadata_seen = False
    smiles_values: list[str] = []
    mapping_values: list[int] = []
    error: str | None = None
    for line in scope.splitlines():
        tokens = line.split()
        if len(tokens) < 2 or tokens[0] != "REMARK" or tokens[1] != "SMILES":
            continue
        # Meeko has emitted both ``REMARK H PARENT`` and the older
        # ``REMARK SMILES H PARENT`` form. Hydrogen-parent metadata is not the
        # heavy-atom graph mapping used here, so preserve and ignore it.
        if len(tokens) >= 4 and tokens[2] == "H" and tokens[3] == "PARENT":
            continue
        metadata_seen = True
        if len(tokens) >= 3 and tokens[2] == "IDX":
            if len(tokens) == 3:
                error = error or "invalid_atom_mapping"
                continue
            try:
                values = [int(value) for value in tokens[3:]]
            except (TypeError, ValueError):
                error = error or "invalid_atom_mapping"
                continue
            if len(values) % 2:
                error = error or "invalid_atom_mapping"
                continue
            mapping_values.extend(values)
            continue
        if len(tokens) != 3:
            error = error or "invalid_smiles_metadata"
            continue
        smiles_values.append(tokens[2])

    if len(smiles_values) > 1:
        error = error or "duplicate_smiles_metadata"
    if mapping_values:
        pairs = list(zip(mapping_values[::2], mapping_values[1::2]))
        if len({pair[0] for pair in pairs}) != len(pairs) or len({pair[1] for pair in pairs}) != len(pairs):
            error = error or "ambiguous_atom_mapping"
    if len(smiles_values) == 0 and mapping_values:
        error = error or "partial_atom_mapping_metadata"
    if len(smiles_values) > 0 and not mapping_values:
        error = error or "partial_atom_mapping_metadata"
    if metadata_seen and not mapping_values and not smiles_values:
        error = error or "partial_atom_mapping_metadata"
    if error:
        return metadata_seen, None, mapping_values, error
    if not metadata_seen:
        return False, None, [], None
    return True, smiles_values[0], mapping_values, None


def _resolve_pdbqt_metadata(text: str, selected: str) -> tuple[str, list[int]]:
    """Resolve chemistry metadata from the selected model or global header.

    A complete selected-model scope takes precedence. A complete global scope
    is used only when the selected model contains no metadata at all. The two
    scopes are never combined.
    """
    _, global_scope, has_models = _split_pdbqt_models(text)
    selected_seen, selected_smiles, selected_mapping, selected_error = _pdbqt_metadata(selected)
    global_seen, global_smiles, global_mapping, global_error = _pdbqt_metadata(global_scope)
    if selected_seen:
        if selected_error:
            raise PoseGeometryError(selected_error)
        if not selected_smiles or not selected_mapping:
            raise PoseGeometryError("partial_atom_mapping_metadata")
        if global_error:
            raise PoseGeometryError(global_error)
        if global_seen:
            selected_pairs = dict(zip(selected_mapping[::2], selected_mapping[1::2]))
            global_pairs = dict(zip(global_mapping[::2], global_mapping[1::2]))
            if selected_smiles != global_smiles or selected_pairs != global_pairs:
                raise PoseGeometryError("conflicting_metadata_scopes")
        return selected_smiles, selected_mapping

    if has_models:
        if global_error:
            raise PoseGeometryError(global_error)
        if not global_seen or not global_smiles or not global_mapping:
            raise PoseGeometryError("missing_authoritative_atom_mapping")
        return global_smiles, global_mapping

    if global_error:
        raise PoseGeometryError(global_error)
    if not global_seen or not global_smiles or not global_mapping:
        raise PoseGeometryError("missing_authoritative_atom_mapping")
    return global_smiles, global_mapping

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
        records = _split_sdf_records(text)
    elif any(_MODEL_LINE.match(line) for line in text.splitlines()):
        records, _, _ = _split_pdbqt_models(text)
    else:
        records = [text]
    if index > len(records):
        raise PoseGeometryError("pose_index_out_of_range")
    selected = records[index - 1]
    if Path(path).suffix.lower() in {".sdf", ".mol"} and not selected.strip():
        raise PoseGeometryError("empty_pose_record")
    return selected

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
        text = path.read_text(encoding="utf-8", errors="strict")
        smiles, numbers = _resolve_pdbqt_metadata(text, record)
        molecule = Chem.MolFromSmiles(smiles)
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
