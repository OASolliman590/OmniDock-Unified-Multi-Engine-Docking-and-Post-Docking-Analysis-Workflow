"""Chemistry backend contract for strict MD-input export.

Open Babel is intentionally imported lazily.  The ordinary post-docking
environment can run without it; the md-inputs stage then reports
``skipped_missing_dependency`` and emits no molecular files.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Protocol, Sequence


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
    pH: float
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
        if protonate and not any(atom.element.upper() == "H" for atom in atoms):
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
            protonated=bool(protonate),
            pH=float(pH),
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
