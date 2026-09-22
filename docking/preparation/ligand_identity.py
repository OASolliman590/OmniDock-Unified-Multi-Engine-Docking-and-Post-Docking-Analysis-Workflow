"""Conservative source/normalised ligand identity contracts.

The preparation backends may add hydrogens and may change formal charges only
when an explicitly selected pH procedure declares that operation.  They may
not replace the heavy atom graph, aromaticity, isotopes, radicals, bond
orders, or defined stereochemistry.  RDKit is loaded lazily because chemistry
is an optional dependency for the rest of OmniDock.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


class LigandIdentityError(ValueError):
    """Raised when ligand identity cannot be established without guessing."""


def sha256_file(path: Path) -> str:
    """Hash exact artifact bytes for durable provenance."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_rdkit() -> Any:
    try:
        from rdkit import Chem
    except ImportError as exc:  # pragma: no cover - dependency/host boundary
        raise ImportError("RDKit is required for the ligand chemistry contract") from exc
    return Chem


def ensure_single_component(molecule: Any, *, context: str = "ligand") -> int:
    """Reject disconnected ordinary ligand input before any backend runs."""

    if molecule is None:
        raise LigandIdentityError(f"{context} chemistry could not be parsed")
    Chem = _require_rdkit()
    try:
        fragments = Chem.GetMolFrags(molecule, asMols=False, sanitizeFrags=False)
        count = len(fragments)
    except Exception as exc:
        raise LigandIdentityError(f"{context} components could not be evaluated") from exc
    if count != 1:
        raise LigandIdentityError(
            f"Ordinary ligand preparation requires exactly one connected component; "
            f"found {count} components in {context}. Multi-component/covalent input "
            "requires an explicit supported preparation mode."
        )
    return count


def _is_hydrogen(atom: Any) -> bool:
    return int(atom.GetAtomicNum()) == 1


def heavy_atom_indices(molecule: Any) -> List[int]:
    return [int(atom.GetIdx()) for atom in molecule.GetAtoms() if not _is_hydrogen(atom)]


def _total_hydrogens(atom: Any) -> int:
    # ``includeNeighbors=True`` counts explicit hydrogen atoms attached to the
    # heavy atom as well as implicit hydrogens.  Omitting it would make an
    # explicit-H source appear to gain/lose a proton merely through indexing.
    try:
        return int(atom.GetTotalNumHs(includeNeighbors=True))
    except TypeError:  # pragma: no cover - old RDKit compatibility
        return int(atom.GetTotalNumHs())


def _coordinate(molecule: Any, index: int) -> Optional[Tuple[float, float, float]]:
    try:
        if molecule.GetNumConformers() == 0:
            return None
        point = molecule.GetConformer().GetAtomPosition(int(index))
        values = (float(point.x), float(point.y), float(point.z))
        return values if all(math.isfinite(value) for value in values) else None
    except Exception:
        return None


def _is_3d(molecule: Any) -> bool:
    if molecule is None or molecule.GetNumConformers() == 0:
        return False
    try:
        if molecule.GetConformer().Is3D():
            return True
    except Exception:
        pass
    positions = molecule.GetConformer().GetPositions()
    return bool(len(positions)) and max(abs(float(position[2])) for position in positions) > 1e-3


def require_finite_3d_coordinates(molecule: Any, *, context: str = "ligand") -> None:
    """Require one finite, explicitly 3D conformer for a prepared artifact.

    A graph identity ledger can be useful for a 2D molecule, but a docking
    input also needs actual coordinates.  Keep this check separate from
    ``build_identity_ledger`` so graph-only callers can still validate an
    input before a backend is selected.
    """

    if molecule is None or molecule.GetNumConformers() == 0 or not _is_3d(molecule):
        raise LigandIdentityError(f"{context} must contain a finite 3D conformer")
    try:
        conformer = molecule.GetConformer()
        positions = conformer.GetPositions()
        if len(positions) != molecule.GetNumAtoms():
            raise LigandIdentityError(f"{context} has incomplete 3D coordinates")
        if not all(
            math.isfinite(float(value))
            for position in positions
            for value in (position[0], position[1], position[2])
        ):
            raise LigandIdentityError(f"{context} has nonfinite 3D coordinates")
    except LigandIdentityError:
        raise
    except Exception as exc:
        raise LigandIdentityError(f"{context} coordinates could not be evaluated") from exc


def _heavy_copy(molecule: Any) -> Tuple[Any, Dict[int, int]]:
    """Return an H-free copy and original-index -> copy-index mapping."""

    Chem = _require_rdkit()
    original_heavy_indices = heavy_atom_indices(molecule)
    copy = Chem.Mol(molecule)
    for atom in copy.GetAtoms():
        atom.SetIntProp("_omnidock_original_index", int(atom.GetIdx()))
    try:
        heavy = Chem.RemoveHs(copy, sanitize=False)
    except TypeError:  # pragma: no cover - old RDKit compatibility
        heavy = Chem.RemoveHs(copy)
    original_to_heavy: Dict[int, int] = {}
    for atom in heavy.GetAtoms():
        if atom.HasProp("_omnidock_original_index"):
            original_to_heavy[int(atom.GetIntProp("_omnidock_original_index"))] = int(atom.GetIdx())
    if len(original_to_heavy) != len(original_heavy_indices):
        # RDKit normally preserves atom properties through RemoveHs.  Keep a
        # conservative positional fallback for older builds where it drops
        # private properties; RemoveHs preserves the relative heavy-atom order.
        if heavy.GetNumAtoms() != len(original_heavy_indices):
            raise LigandIdentityError("RDKit could not retain heavy-atom indices while removing hydrogens")
        original_to_heavy = {
            original_index: heavy_index
            for heavy_index, original_index in enumerate(original_heavy_indices)
        }
    return heavy, original_to_heavy


def _neutral_graph_copy(molecule: Any) -> Tuple[Any, Dict[int, int]]:
    """Create an H-free graph query while ignoring only formal charges.

    Formal charge is checked independently below, so neutralising the query
    allows RDKit's supported ``useChirality`` matching to establish a mapping
    across an allowed protonation change.  The copy is never sanitised after
    neutralisation; its graph and stereo tags remain authoritative.
    """

    heavy, original_to_heavy = _heavy_copy(molecule)
    for atom in heavy.GetAtoms():
        atom.SetFormalCharge(0)
    return heavy, original_to_heavy


def _mapping_candidates(
    source: Any,
    normalized: Any,
    *,
    max_mappings: int = 10000,
) -> List[Dict[int, int]]:
    source_indices = heavy_atom_indices(source)
    normalized_indices = heavy_atom_indices(normalized)
    if len(source_indices) != len(normalized_indices):
        return []
    source_graph, source_to_graph = _neutral_graph_copy(source)
    normalized_graph, normalized_to_graph = _neutral_graph_copy(normalized)
    if len(source_graph.GetAtoms()) != len(normalized_graph.GetAtoms()):
        return []
    try:
        matches = source_graph.GetSubstructMatches(
            normalized_graph,
            uniquify=False,
            useChirality=True,
            maxMatches=max_mappings + 1,
        )
    except Exception as exc:
        raise LigandIdentityError(
            "RDKit could not establish a stereochemistry-aware heavy-atom mapping"
        ) from exc
    if len(matches) >= max_mappings + 1:
        raise LigandIdentityError(
            f"heavy-atom mapping is ambiguous beyond the safety limit of {max_mappings} mappings"
        )
    if not matches:
        return []
    graph_to_source = {graph_index: source_index for source_index, graph_index in source_to_graph.items()}
    graph_to_normalized = {
        graph_index: normalized_index for normalized_index, graph_index in normalized_to_graph.items()
    }
    mappings: List[Dict[int, int]] = []
    for match in matches:
        if len(match) != len(source_graph.GetAtoms()):
            continue
        try:
            # ``GetSubstructMatches`` is called on the source target with the
            # normalized molecule as query: each tuple position is a query
            # (normalized) index and its value is a source target index.
            mappings.append(
                {
                    graph_to_source[source_graph_index]: graph_to_normalized[normalized_graph_index]
                    for normalized_graph_index, source_graph_index in enumerate(match)
                }
            )
        except KeyError as exc:
            raise LigandIdentityError("RDKit returned an incomplete heavy-atom mapping") from exc
    return mappings


def _choose_mapping(source: Any, normalized: Any, mappings: Sequence[Mapping[int, int]]) -> Dict[int, int]:
    source_indices = heavy_atom_indices(source)
    normalized_coordinates = {index: _coordinate(normalized, index) for index in heavy_atom_indices(normalized)}
    source_coordinates = {index: _coordinate(source, index) for index in source_indices}

    def distance(mapping: Mapping[int, int]) -> float:
        total = 0.0
        for source_index, normalized_index in mapping.items():
            left = source_coordinates[source_index]
            right = normalized_coordinates[normalized_index]
            if left is not None and right is not None:
                total += sum((a - b) ** 2 for a, b in zip(left, right))
        return total

    return dict(min(mappings, key=lambda mapping: (distance(mapping), tuple(sorted(mapping.items())))))


def _stereo_defined(atom: Any) -> bool:
    return str(atom.GetChiralTag()).strip() not in {"", "CHI_UNSPECIFIED"}


def _stereo_state(molecule: Any, labels: Mapping[int, int]) -> Tuple[str, str]:
    """Return mapping-aware atom/bond stereo signatures.

    RDKit's ``useChirality=True`` match is the acceptance check.  These
    signatures make an undefined-to-defined transition explicit and provide a
    stable audit record; they are deliberately not used as a CIP R/S shortcut.
    """

    source_indices = heavy_atom_indices(molecule)
    atom_states = []
    for index in source_indices:
        atom = molecule.GetAtomWithIdx(index)
        atom_states.append((labels.get(index, index), _stereo_defined(atom)))
    bond_states = []
    for source_index in source_indices:
        for neighbor in molecule.GetAtomWithIdx(source_index).GetNeighbors():
            neighbor_index = int(neighbor.GetIdx())
            if _is_hydrogen(neighbor) or neighbor_index <= source_index:
                continue
            bond = molecule.GetBondBetweenAtoms(source_index, neighbor_index)
            stereo = str(bond.GetStereo()).strip()
            left = labels.get(source_index, source_index)
            right = labels.get(neighbor_index, neighbor_index)
            bond_states.append(
                (
                    min(left, right),
                    max(left, right),
                    stereo not in {"", "STEREONONE"},
                )
            )
    atom_states.sort(key=lambda item: item[0])
    bond_states.sort(key=lambda item: (item[0], item[1]))
    return repr(atom_states), repr(bond_states)


def _heavy_bond_descriptor(bond: Any) -> Optional[Tuple[float, bool, bool]]:
    if bond is None:
        return None
    return (
        float(bond.GetBondTypeAsDouble()),
        bool(bond.GetIsAromatic()),
        bool(bond.GetIsConjugated()),
    )


def _isomeric_smiles(molecule: Any) -> str:
    Chem = _require_rdkit()
    try:
        heavy, _ = _heavy_copy(molecule)
        return str(Chem.MolToSmiles(heavy, isomericSmiles=True))
    except Exception as exc:
        raise LigandIdentityError("RDKit could not serialize the ligand isomeric graph") from exc


def build_identity_ledger(
    source: Any,
    normalized: Any,
    *,
    allow_hydrogen_changes: bool = True,
    allow_formal_charge_changes: bool = False,
    preserve_coordinates: bool = True,
    coordinate_precision: int = 3,
) -> Dict[str, object]:
    """Compare molecules, allowing only explicitly selected H/charge changes."""

    ensure_single_component(source, context="source ligand")
    ensure_single_component(normalized, context="normalized ligand")
    source_indices = heavy_atom_indices(source)
    normalized_indices = heavy_atom_indices(normalized)
    mappings = _mapping_candidates(source, normalized)
    differences: Dict[str, List[Dict[str, object]]] = {"hydrogens": [], "formal_charges": []}
    if not mappings:
        reason = (
            f"heavy-atom count changed ({len(source_indices)} -> {len(normalized_indices)})"
            if len(source_indices) != len(normalized_indices)
            else "normalized ligand heavy-atom graph or defined stereochemistry changed"
        )
        return {
            "is_valid": False,
            "errors": [reason],
            "warnings": [],
            "source_heavy_atom_count": len(source_indices),
            "normalized_heavy_atom_count": len(normalized_indices),
            "atom_map": [],
            "differences": differences,
            "source_isomeric_smiles": _isomeric_smiles(source),
            "normalized_isomeric_smiles": _isomeric_smiles(normalized),
            "source_is_3d": _is_3d(source),
            "normalized_is_3d": _is_3d(normalized),
            "mapping_count": 0,
        }
    mapping = _choose_mapping(source, normalized, mappings)
    errors: List[str] = []
    warnings: List[str] = []
    atom_map: List[Dict[str, object]] = []
    for source_index in sorted(mapping):
        normalized_index = mapping[source_index]
        source_atom = source.GetAtomWithIdx(source_index)
        normalized_atom = normalized.GetAtomWithIdx(normalized_index)
        source_h = _total_hydrogens(source_atom)
        normalized_h = _total_hydrogens(normalized_atom)
        source_charge = int(source_atom.GetFormalCharge())
        normalized_charge = int(normalized_atom.GetFormalCharge())
        delta_h = normalized_h - source_h
        delta_charge = normalized_charge - source_charge
        if (
            allow_formal_charge_changes
            and (delta_h or delta_charge)
            and delta_charge != delta_h
        ):
            errors.append(
                f"formal charge/proton delta is not a paired protonation change for mapped heavy atom "
                f"{source_index}->{normalized_index} (delta charge {delta_charge}, delta H {delta_h})"
            )
        if delta_h:
            differences["hydrogens"].append(
                {
                    "source_index": source_index,
                    "normalized_index": normalized_index,
                    "source": source_h,
                    "normalized": normalized_h,
                    "delta": delta_h,
                }
            )
            if not allow_hydrogen_changes:
                errors.append(f"hydrogen count changed for mapped heavy atom {source_index}->{normalized_index}")
        if delta_charge:
            differences["formal_charges"].append(
                {
                    "source_index": source_index,
                    "normalized_index": normalized_index,
                    "source": source_charge,
                    "normalized": normalized_charge,
                    "delta": delta_charge,
                }
            )
            if not allow_formal_charge_changes:
                errors.append(f"formal charge changed for mapped heavy atom {source_index}->{normalized_index}")
        atom_map.append(
            {
                "source_index": source_index,
                "normalized_index": normalized_index,
                "element": str(source_atom.GetSymbol()),
                "isotope": int(source_atom.GetIsotope()),
                "radical_electrons": int(source_atom.GetNumRadicalElectrons()),
                "aromatic": bool(source_atom.GetIsAromatic()),
                "source_formal_charge": source_charge,
                "normalized_formal_charge": normalized_charge,
                "source_hydrogens": source_h,
                "normalized_hydrogens": normalized_h,
                "source_coordinates": _coordinate(source, source_index),
                "normalized_coordinates": _coordinate(normalized, normalized_index),
            }
        )

    # RDKit's graph matcher supplies the candidate mapping, but explicit
    # element/isotope/radical/aromatic/bond checks keep the contract robust to
    # query defaults that may treat an unspecified isotope as a wildcard.
    for source_index, normalized_index in mapping.items():
        source_atom = source.GetAtomWithIdx(source_index)
        normalized_atom = normalized.GetAtomWithIdx(normalized_index)
        if (
            int(source_atom.GetAtomicNum()) != int(normalized_atom.GetAtomicNum())
            or str(source_atom.GetSymbol()) != str(normalized_atom.GetSymbol())
            or int(source_atom.GetIsotope()) != int(normalized_atom.GetIsotope())
            or int(source_atom.GetNumRadicalElectrons()) != int(normalized_atom.GetNumRadicalElectrons())
            or bool(source_atom.GetIsAromatic()) != bool(normalized_atom.GetIsAromatic())
        ):
            errors.append(f"mapped heavy-atom element/isotope/radical/aromatic identity changed for {source_index}->{normalized_index}")
        for other_source, other_normalized in mapping.items():
            if other_source <= source_index:
                continue
            source_bond = source.GetBondBetweenAtoms(source_index, other_source)
            normalized_bond = normalized.GetBondBetweenAtoms(normalized_index, other_normalized)
            if _heavy_bond_descriptor(source_bond) != _heavy_bond_descriptor(normalized_bond):
                errors.append(
                    f"mapped heavy-atom connectivity or bond order changed for {source_index}-{other_source}"
                )

    # The RDKit match above checks chirality and stereobonds in the graph's
    # actual neighbour ordering.  Explicit state comparison below catches a
    # missing/added stereo tag even when the query atom was unspecified.
    source_stereo = _stereo_state(source, {index: index for index in source_indices})
    normalized_stereo = _stereo_state(normalized, {value: key for key, value in mapping.items()})
    if source_stereo != normalized_stereo:
        errors.append("defined stereochemistry or stereobond state changed in the mapped ligand graph")

    source_is_3d = _is_3d(source)
    normalized_is_3d = _is_3d(normalized)
    if preserve_coordinates and source_is_3d:
        if not normalized_is_3d:
            errors.append("normalization removed the source 3D conformer")
        else:
            for source_index, normalized_index in mapping.items():
                source_position = _coordinate(source, source_index)
                normalized_position = _coordinate(normalized, normalized_index)
                if source_position is None or normalized_position is None:
                    errors.append(f"missing finite coordinates for mapped heavy atom {source_index}->{normalized_index}")
                    continue
                if tuple(round(value, coordinate_precision) for value in source_position) != tuple(
                    round(value, coordinate_precision) for value in normalized_position
                ):
                    errors.append(
                        f"source coordinates changed beyond {coordinate_precision}-decimal precision for "
                        f"mapped heavy atom {source_index}->{normalized_index}"
                    )

    if differences["formal_charges"] and allow_formal_charge_changes and not any(
        "paired protonation" in error for error in errors
    ):
        warnings.append("formal-charge changes accepted under the declared pH normalization procedure")
    if differences["hydrogens"] and allow_hydrogen_changes:
        warnings.append("hydrogen-count changes accepted under the declared normalization procedure")
    return {
        "is_valid": not errors,
        "errors": errors,
        "warnings": warnings,
        "source_heavy_atom_count": len(source_indices),
        "normalized_heavy_atom_count": len(normalized_indices),
        "atom_map": atom_map,
        "differences": differences,
        "source_isomeric_smiles": _isomeric_smiles(source),
        "normalized_isomeric_smiles": _isomeric_smiles(normalized),
        "source_is_3d": source_is_3d,
        "normalized_is_3d": normalized_is_3d,
        "mapping_count": len(mappings),
        "mapping_scope": "rdkit_stereochemistry_aware;heavy_atoms_only",
    }


def validate_mapped_pdbqt(
    input_molecule: Any,
    output_file: Path,
    *,
    coordinate_precision: int = 3,
) -> Dict[str, object]:
    """Validate a PDBQT only when its authoritative Meeko mapping is usable."""

    try:
        require_finite_3d_coordinates(input_molecule, context="prepared input ligand")
    except LigandIdentityError as exc:
        return {
            "status": "invalid",
            "is_valid": False,
            "scope": "input_coordinates_required",
            "errors": [str(exc)],
            "warnings": [],
        }

    try:
        from post_docking_analysis.pose_geometry import load_pose_molecule

        output_molecule = load_pose_molecule(Path(output_file), pose_index=1)
    except Exception as exc:
        reason = str(exc)
        absent_mapping = any(
            token in reason
            for token in (
                "missing_authoritative_atom_mapping",
                "coordinate_only_format_requires_chemical_template",
            )
        )
        return {
            "status": "limited" if absent_mapping else "invalid",
            "is_valid": False,
            "scope": (
                "pdbqt_syntax_only;authoritative_atom_mapping_unavailable"
                if absent_mapping
                else "pdbqt_authoritative_mapping_malformed"
            ),
            "errors": [
                f"authoritative PDBQT atom mapping {'unavailable' if absent_mapping else 'invalid'}: {exc}"
            ],
            "warnings": ["PDBQT output was not certified for exact chemistry because no authoritative mapping was available"],
        }
    ledger = build_identity_ledger(
        input_molecule,
        output_molecule,
        allow_hydrogen_changes=False,
        allow_formal_charge_changes=False,
        preserve_coordinates=True,
        coordinate_precision=coordinate_precision,
    )
    if not ledger["is_valid"]:
        return {
            "status": "invalid",
            "is_valid": False,
            "scope": "meeko_authoritative_mapping;heavy_atom_graph_and_coordinates_exact",
            "errors": list(ledger["errors"]),
            "warnings": list(ledger["warnings"]),
            "identity_ledger": ledger,
        }
    return {
        "status": "exact",
        "is_valid": True,
        "scope": "meeko_authoritative_mapping;heavy_atom_graph_and_coordinates_exact",
        "errors": [],
        "warnings": list(ledger["warnings"]),
        "identity_ledger": ledger,
    }
