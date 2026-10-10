"""Chemistry backend contract for strict MD-input export.

Open Babel is intentionally imported lazily.  The ordinary post-docking
environment can run without it; the md-inputs stage then reports
``skipped_missing_dependency`` and emits no molecular files.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Protocol, Sequence, Tuple


@dataclass
class MDAtom:
    index: int
    element: str
    x: float
    y: float
    z: float
    formal_charge: int = 0
    name: str = ""


@dataclass
class MDBond:
    begin: int
    end: int
    order: str


@dataclass
class PreparedLigand:
    atoms: List[MDAtom]
    bonds: List[MDBond]
    net_charge: int
    backend: str
    backend_version: str
    protonated: bool
    pH: Optional[float]  # None for a docked microspecies (no pH model was run)
    implicit_hydrogen_count: int
    heavy_coordinate_max_delta: float
    native: object = field(default=None, repr=False)


class ChemistryBackend(Protocol):
    name: str
    version: str

    def prepare(
        self,
        source_file: Path,
        *,
        source_pose: int,
        pH: float,
        protonate: bool,
        pose_heavy_coordinates: Optional[Sequence[Sequence[float]]] = None,
        pose_heavy_elements: Optional[Sequence[str]] = None,
        pose_to_topology: Optional[Sequence[int]] = None,
    ) -> PreparedLigand: ...

    def prepare_docked_state(
        self,
        source_file: Path,
        *,
        source_pose: int,
        expected_net_charge: int,
        add_hydrogens: bool,
        pose_heavy_coordinates: Optional[Sequence[Sequence[float]]] = None,
        pose_heavy_elements: Optional[Sequence[str]] = None,
        pose_to_smiles: Optional[Sequence[int]] = None,
        microspecies_smiles: Optional[str] = None,
    ) -> PreparedLigand: ...

    def write(self, ligand: PreparedLigand, output_dir: Path, formats: Sequence[str]) -> Dict[str, Path]: ...


def _normalize_element(value: object) -> str:
    letters = "".join(ch for ch in str(value or "").strip() if ch.isalpha())
    if not letters:
        return ""
    if len(letters) > 1 and letters[1].islower():
        return letters[:2].title()
    return letters[:1].upper()


def assign_unique_atom_names(atoms: Sequence[MDAtom]) -> List[str]:
    """Assign deterministic, residue-unique atom names of at most four chars."""
    counts: Dict[str, int] = {}
    names: List[str] = []
    for atom in atoms:
        element = _normalize_element(atom.element)
        if not element:
            raise ValueError("missing_atom_element")
        counts[element] = counts.get(element, 0) + 1
        name = f"{element}{counts[element]}"
        if len(name) > 4:
            raise ValueError(f"atom_name_overflow:{element}:{counts[element]}")
        atom.name = name
        names.append(name)
    if len(set(names)) != len(names):
        raise ValueError("duplicate_atom_names")
    return names


def openbabel_capability() -> Dict[str, object]:
    try:
        from openbabel import openbabel as ob  # type: ignore
    except Exception as exc:
        return {
            "status": "skipped_missing_dependency",
            "backend": "openbabel",
            "backend_version": "",
            "reason": f"openbabel_import_failed:{exc.__class__.__name__}",
        }
    version = ""
    try:
        version = str(ob.OBReleaseVersion())
    except Exception:
        version = "unknown"
    return {
        "status": "completed",
        "backend": "openbabel",
        "backend_version": version,
        "reason": "",
    }


def get_openbabel_backend() -> Optional["OpenBabelChemistryBackend"]:
    capability = openbabel_capability()
    if capability["status"] != "completed":
        return None
    return OpenBabelChemistryBackend()


def _selected_sdf_record(path: Path, pose_number: int) -> str:
    records = [record for record in path.read_text(encoding="utf-8", errors="strict").split("$$$$") if record.strip()]
    if pose_number < 1 or pose_number > len(records):
        raise ValueError(f"pose_index_out_of_range:{pose_number}>{len(records)}")
    return records[pose_number - 1].lstrip("\r\n") + "\n$$$$\n"


class OpenBabelChemistryBackend:
    """Open Babel 3 adapter used by the approved predicted-protonation path."""

    name = "openbabel"

    def __init__(self) -> None:
        from openbabel import openbabel as ob  # type: ignore

        self.ob = ob
        try:
            self.version = str(ob.OBReleaseVersion())
        except Exception:
            self.version = "unknown"

    def _read_source(self, source_file: Path, source_pose: int):
        source = Path(source_file).expanduser().resolve()
        suffix = source.suffix.lower()
        conversion = self.ob.OBConversion()
        mol = self.ob.OBMol()
        if suffix in {".sdf", ".mol"}:
            if not conversion.SetInFormat("sdf" if suffix == ".sdf" else "mol"):
                raise RuntimeError(f"openbabel_input_format_unavailable:{suffix}")
            if suffix == ".sdf":
                ok = conversion.ReadString(mol, _selected_sdf_record(source, source_pose))
            else:
                if source_pose != 1:
                    raise ValueError("mol_topology_supports_pose_1_only")
                ok = conversion.ReadFile(mol, str(source))
        elif suffix == ".mol2":
            if source_pose != 1:
                raise ValueError("mol2_topology_supports_pose_1_only")
            if not conversion.SetInFormat("mol2"):
                raise RuntimeError("openbabel_input_format_unavailable:mol2")
            ok = conversion.ReadFile(mol, str(source))
        else:
            raise ValueError(f"unsupported_connectivity_format:{suffix}")
        if not ok or mol.NumAtoms() <= 0 or mol.NumBonds() <= 0:
            raise ValueError("connectivity_source_parse_failed")
        return mol

    def prepare(
        self,
        source_file: Path,
        *,
        source_pose: int,
        pH: float,
        protonate: bool,
        pose_heavy_coordinates: Optional[Sequence[Sequence[float]]] = None,
        pose_heavy_elements: Optional[Sequence[str]] = None,
        pose_to_topology: Optional[Sequence[int]] = None,
    ) -> PreparedLigand:
        mol = self._read_source(source_file, int(source_pose))
        heavy_atoms = [atom for atom in self.ob.OBMolAtomIter(mol) if int(atom.GetAtomicNum()) != 1]

        if pose_heavy_coordinates is not None:
            coords = [tuple(float(value) for value in row) for row in pose_heavy_coordinates]
            elements = [_normalize_element(value) for value in (pose_heavy_elements or [])]
            if len(coords) != len(heavy_atoms) or len(elements) != len(heavy_atoms):
                raise ValueError("topology_pose_heavy_atom_count_mismatch")
            mapping = list(pose_to_topology or [])
            if not mapping:
                raise ValueError("missing_pose_to_topology_atom_map")
            if sorted(mapping) != list(range(len(heavy_atoms))):
                raise ValueError("invalid_pose_to_topology_atom_map")
            topology_elements = [_normalize_element(self.ob.GetSymbol(atom.GetAtomicNum())) for atom in heavy_atoms]
            for pose_index, topology_index in enumerate(mapping):
                if elements[pose_index] != topology_elements[topology_index]:
                    raise ValueError("topology_atom_map_element_mismatch")
                x, y, z = coords[pose_index]
                heavy_atoms[topology_index].SetVector(x, y, z)

        before = [
            (float(atom.GetX()), float(atom.GetY()), float(atom.GetZ()))
            for atom in self.ob.OBMolAtomIter(mol)
            if int(atom.GetAtomicNum()) != 1
        ]

        if protonate:
            mol.DeleteHydrogens()
            if not bool(mol.CorrectForPH(float(pH))):
                raise ValueError("openbabel_ph_correction_failed")
            if not bool(mol.AddHydrogens()):
                # AddHydrogens returns false when no atoms were added; that is valid only
                # when the molecule is already valence-complete after pH correction.
                pass

        after_heavy = [atom for atom in self.ob.OBMolAtomIter(mol) if int(atom.GetAtomicNum()) != 1]
        maximum_delta = 0.0
        if len(after_heavy) != len(before):
            raise ValueError("heavy_atom_identity_changed_during_protonation")
        for atom, old in zip(after_heavy, before):
            delta = max(abs(float(atom.GetX()) - old[0]), abs(float(atom.GetY()) - old[1]), abs(float(atom.GetZ()) - old[2]))
            maximum_delta = max(maximum_delta, delta)
        if maximum_delta > 1.0e-4:
            raise ValueError(f"heavy_atom_coordinates_changed:{maximum_delta:.8f}")
        return self._package(mol, protonated=bool(protonate), pH=float(pH), maximum_delta=maximum_delta)

    def prepare_docked_state(
        self,
        source_file: Path,
        *,
        source_pose: int,
        expected_net_charge: int,
        add_hydrogens: bool,
        pose_heavy_coordinates: Optional[Sequence[Sequence[float]]] = None,
        pose_heavy_elements: Optional[Sequence[str]] = None,
        pose_to_smiles: Optional[Sequence[int]] = None,
        microspecies_smiles: Optional[str] = None,
    ) -> PreparedLigand:
        """Spec 036 R1b: build the MD ligand from the docked microspecies, with no pH model.

        With ``microspecies_smiles`` (Meeko lineage), the heavy atoms are taken from the pose
        coordinates through ``pose_to_smiles`` (pose heavy index -> 0-based SMILES heavy index)
        and the formal charges come from the SMILES. Without it, the docked SDF record is used
        as supplied (its charges and hydrogens are the docked state). Hydrogens are added at
        the pose coordinates without ``-p``; the net charge must equal ``expected_net_charge``.
        """
        if microspecies_smiles is None:
            mol = self._read_source(Path(source_file), int(source_pose))
            reference: List[Tuple[float, float, float]] = []
            before_h = sum(1 for atom in self.ob.OBMolAtomIter(mol) if int(atom.GetAtomicNum()) == 1)
        else:
            mol, reference = self._molecule_from_microspecies(
                microspecies_smiles,
                pose_heavy_coordinates=pose_heavy_coordinates,
                pose_heavy_elements=pose_heavy_elements,
                pose_to_smiles=pose_to_smiles,
            )
            before_h = sum(1 for atom in self.ob.OBMolAtomIter(mol) if int(atom.GetAtomicNum()) == 1)
        implicit_before = sum(int(atom.GetImplicitHCount()) for atom in self.ob.OBMolAtomIter(mol))
        if add_hydrogens:
            mol.SetDimension(3)
            mol.AddHydrogens()
        after_h = sum(1 for atom in self.ob.OBMolAtomIter(mol) if int(atom.GetAtomicNum()) == 1)
        if after_h - before_h != (implicit_before if add_hydrogens else 0):
            raise ValueError(f"hydrogen_count_inconsistent_with_docked_state:{after_h - before_h}!={implicit_before}")
        net = int(mol.GetTotalCharge())
        if net != int(expected_net_charge):
            raise ValueError(f"md_state_differs_from_docked_state:{net}!={int(expected_net_charge)}")
        maximum_delta = 0.0
        if reference:
            heavy_atoms = [atom for atom in self.ob.OBMolAtomIter(mol) if int(atom.GetAtomicNum()) != 1]
            for index, atom in enumerate(heavy_atoms):
                x, y, z = reference[index]
                delta = max(abs(float(atom.GetX()) - x), abs(float(atom.GetY()) - y), abs(float(atom.GetZ()) - z))
                maximum_delta = max(maximum_delta, delta)
            if maximum_delta > 1.0e-4:
                raise ValueError(f"heavy_atom_coordinates_changed:{maximum_delta:.8f}")
        return self._package(mol, protonated=bool(add_hydrogens), pH=None, maximum_delta=maximum_delta)

    def _molecule_from_microspecies(
        self,
        smiles: str,
        *,
        pose_heavy_coordinates: Optional[Sequence[Sequence[float]]],
        pose_heavy_elements: Optional[Sequence[str]],
        pose_to_smiles: Optional[Sequence[int]],
    ):
        if pose_heavy_coordinates is None or pose_heavy_elements is None or pose_to_smiles is None:
            raise ValueError("docked_state_requires_pose_coordinates_and_lineage")
        conversion = self.ob.OBConversion()
        if not conversion.SetInFormat("smi"):
            raise RuntimeError("openbabel_smiles_format_unavailable")
        mol = self.ob.OBMol()
        if not conversion.ReadString(mol, smiles) or mol.NumAtoms() <= 0:
            raise ValueError("microspecies_smiles_unreadable")
        if any(int(atom.GetAtomicNum()) == 1 for atom in self.ob.OBMolAtomIter(mol)):
            raise ValueError("microspecies_smiles_must_be_heavy_atoms_only")
        mapping = [int(item) for item in pose_to_smiles]
        if len(mapping) != mol.NumAtoms() or sorted(mapping) != list(range(mol.NumAtoms())):
            raise ValueError("pose_to_smiles_not_permutation")
        coords = [tuple(float(value) for value in row) for row in pose_heavy_coordinates]
        elements = [_normalize_element(value) for value in pose_heavy_elements]
        if len(coords) != len(mapping) or len(elements) != len(mapping):
            raise ValueError("docked_state_heavy_atom_count_mismatch")
        reference: List[Tuple[float, float, float]] = [(0.0, 0.0, 0.0)] * len(mapping)
        for pose_index, smiles_index in enumerate(mapping):
            atom = mol.GetAtom(smiles_index + 1)
            if _normalize_element(self.ob.GetSymbol(atom.GetAtomicNum())) != elements[pose_index]:
                raise ValueError("docked_state_atom_map_element_mismatch")
            x, y, z = coords[pose_index]
            atom.SetVector(x, y, z)
            reference[smiles_index] = (x, y, z)
        # Reference is kept in SMILES order for the coordinate check, which walks the OBMol in order.
        ordered_reference = [reference[index] for index in range(mol.NumAtoms())]
        return mol, ordered_reference

    def _package(self, mol, *, protonated: bool, pH: Optional[float], maximum_delta: float) -> PreparedLigand:
        """Measure a finished OBMol and return the MD-facing ligand record (shared by every path)."""
        atoms: List[MDAtom] = []
        implicit_hydrogen_count = 0
        for ordinal, atom in enumerate(self.ob.OBMolAtomIter(mol), start=1):
            implicit_hydrogen_count += int(atom.GetImplicitHCount())
            atoms.append(
                MDAtom(
                    index=ordinal,
                    element=_normalize_element(self.ob.GetSymbol(atom.GetAtomicNum())),
                    x=float(atom.GetX()),
                    y=float(atom.GetY()),
                    z=float(atom.GetZ()),
                    formal_charge=int(atom.GetFormalCharge()),
                )
            )
        if implicit_hydrogen_count != 0:
            raise ValueError(f"implicit_hydrogens_remain:{implicit_hydrogen_count}")
        if protonated and not any(atom.element.upper() == "H" for atom in atoms):
            raise ValueError("protonation_added_no_explicit_hydrogens")

        bonds: List[MDBond] = []
        for bond in self.ob.OBMolBondIter(mol):
            order = "ar" if bool(bond.IsAromatic()) else str(int(bond.GetBondOrder()))
            bonds.append(MDBond(int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx()), order))

        assign_unique_atom_names(atoms)
        residue = mol.NewResidue()
        residue.SetName("LIG")
        residue.SetNum(1)
        for atom_record, ob_atom in zip(atoms, self.ob.OBMolAtomIter(mol)):
            residue.AddAtom(ob_atom)
            residue.SetAtomID(ob_atom, atom_record.name)

        return PreparedLigand(
            atoms=atoms,
            bonds=bonds,
            net_charge=int(mol.GetTotalCharge()),
            backend=self.name,
            backend_version=self.version,
            protonated=bool(protonated),
            pH=pH,
            implicit_hydrogen_count=implicit_hydrogen_count,
            heavy_coordinate_max_delta=maximum_delta,
            native=mol,
        )

    def write(self, ligand: PreparedLigand, output_dir: Path, formats: Sequence[str]) -> Dict[str, Path]:
        output_root = Path(output_dir).expanduser().resolve()
        output_root.mkdir(parents=True, exist_ok=True)
        requested = {str(value).strip().lower() for value in formats if str(value).strip()}
        requested.add("mol2")
        written: Dict[str, Path] = {}
        for fmt in sorted(requested):
            if fmt not in {"mol2", "sdf", "pdb"}:
                raise ValueError(f"unsupported_ligand_output_format:{fmt}")
            conversion = self.ob.OBConversion()
            out_format = "sdf" if fmt == "sdf" else fmt
            if not conversion.SetOutFormat(out_format):
                raise RuntimeError(f"openbabel_output_format_unavailable:{fmt}")
            output = output_root / f"ligand.{fmt}"
            if not conversion.WriteFile(ligand.native, str(output)):
                raise RuntimeError(f"openbabel_write_failed:{fmt}")
            conversion.CloseOutFile()
            if not output.exists() or output.stat().st_size <= 0:
                raise RuntimeError(f"empty_ligand_output:{fmt}")
            if fmt == "mol2":
                lines = output.read_text(encoding="utf-8", errors="strict").splitlines()
                section = ""
                atom_index = 0
                rewritten: List[str] = []
                for line in lines:
                    if line.startswith("@<TRIPOS>"):
                        section = line.strip().upper()
                        rewritten.append(line)
                        continue
                    if section == "@<TRIPOS>ATOM" and line.strip():
                        tokens = line.split()
                        if len(tokens) < 6 or atom_index >= len(ligand.atoms):
                            raise RuntimeError("invalid_openbabel_mol2_atom_section")
                        tokens[1] = ligand.atoms[atom_index].name
                        atom_index += 1
                        rewritten.append(" ".join(tokens))
                    else:
                        rewritten.append(line)
                if atom_index != len(ligand.atoms):
                    raise RuntimeError("openbabel_mol2_atom_count_mismatch")
                output.write_text("\n".join(rewritten) + "\n", encoding="utf-8")
            verifier = self.ob.OBConversion()
            verified = self.ob.OBMol()
            if not verifier.SetInFormat(out_format) or not verifier.ReadFile(verified, str(output)):
                raise RuntimeError(f"ligand_output_parseback_failed:{fmt}")
            if verified.NumAtoms() != len(ligand.atoms) or verified.NumBonds() != len(ligand.bonds):
                raise RuntimeError(
                    f"ligand_output_graph_count_mismatch:{fmt}:"
                    f"{verified.NumAtoms()}/{verified.NumBonds()}!={len(ligand.atoms)}/{len(ligand.bonds)}"
                )
            written[fmt] = output
        return written


# ---------------------------------------------------------------- Spec 036 R1b: docked microspecies


def _rdkit_chem():
    try:
        from rdkit import Chem  # type: ignore
    except ImportError as exc:
        raise RuntimeError("rdkit_unavailable") from exc
    return Chem


def docked_microspecies_from_pdbqt(pose_file: Path, pose: int) -> Dict[str, object]:
    """Docked microspecies of one Vina/Smina pose from its Meeko ``REMARK SMILES`` lineage.

    ``pose_to_smiles[i]`` is the 0-based SMILES heavy-atom index of pose heavy atom ``i``, in the
    same heavy-atom order that ``load_pdbqt_lineage_pose`` uses for the coordinates. Raises
    ``ValueError`` with a machine-readable reason when the lineage cannot be used.
    """
    from docking.preparation.ligand_preparation import canonical_microspecies
    from post_docking_analysis.atom_mapping import _lineage_remarks, _pdbqt_atoms, _pdbqt_model_lines

    lines = _pdbqt_model_lines(Path(pose_file), int(pose))
    smiles, flat_pairs = _lineage_remarks(lines)
    heavy_serials = [serial for serial, _, _, is_heavy in _pdbqt_atoms(lines) if is_heavy]
    position = {serial: index for index, serial in enumerate(heavy_serials)}
    pose_to_smiles: List[int] = [-1] * len(heavy_serials)
    for index in range(0, len(flat_pairs), 2):
        smiles_index_1based, serial = flat_pairs[index], flat_pairs[index + 1]
        if serial not in position:
            raise ValueError("lineage_idx_serial_not_a_heavy_atom")
        pose_to_smiles[position[serial]] = smiles_index_1based - 1
    if any(item < 0 for item in pose_to_smiles):
        raise ValueError("incomplete_lineage_idx")
    Chem = _rdkit_chem()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("invalid_lineage_smiles")
    if mol.GetNumAtoms() != len(heavy_serials):
        raise ValueError("lineage_smiles_heavy_atom_count_mismatch")
    return {
        "source": "meeko_remark_smiles",
        "smiles": smiles,
        "canonical_smiles": canonical_microspecies(mol),
        "net_formal_charge": int(Chem.GetFormalCharge(mol)),
        "formal_charges": [int(atom.GetFormalCharge()) for atom in mol.GetAtoms()],
        "pose_to_smiles": pose_to_smiles,
        "heavy_atom_count": len(heavy_serials),
    }


def docked_microspecies_from_sdf(pose_file: Path, pose: int) -> Dict[str, object]:
    """Docked microspecies of one GNINA SDF pose record: its own charges and hydrogens."""
    from docking.preparation.ligand_preparation import canonical_microspecies

    records = [record for record in Path(pose_file).read_text(encoding="utf-8", errors="strict").split("$$$$") if record.strip()]
    if pose < 1 or pose > len(records):
        raise ValueError(f"pose_index_out_of_range:{pose}>{len(records)}")
    Chem = _rdkit_chem()
    mol = Chem.MolFromMolBlock(records[pose - 1].lstrip("\r\n"), sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError("pose_sdf_record_unreadable")
    canonical = ""
    sanitized = False
    try:
        Chem.SanitizeMol(mol)
        canonical = canonical_microspecies(mol)
        sanitized = True
    except Exception:
        canonical = ""
    return {
        "source": "native_pose_sdf",
        "smiles": canonical,
        "canonical_smiles": canonical,
        "sanitized": sanitized,
        "net_formal_charge": int(Chem.GetFormalCharge(mol)),
        "formal_charges": [int(atom.GetFormalCharge()) for atom in mol.GetAtoms()],
        "pose_to_smiles": None,
        "heavy_atom_count": sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1),
    }
