"""Chemistry-required ligand preparation contracts.

These tests intentionally import RDKit at collection time.  Local filtered
verification explicitly excludes this chemistry-required module when the
Windows host blocks RDKit's optional DLL; the chemistry CI environment must
execute it.
"""
from __future__ import annotations

import json
from pathlib import Path
import subprocess

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from docking.preparation import ligand_preparation as preparation
from docking.preparation.ligand_identity import (
    LigandIdentityError,
    build_identity_ledger,
    ensure_single_component,
    require_finite_3d_coordinates,
    validate_mapped_pdbqt,
)
from docking.preparation.structure_contract import coordinates as fixed_coordinates


def _embedded(smiles: str, seed: int = 19):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(molecule, randomSeed=seed) == 0
    return molecule


def _copy_heavy_coordinates(source, target):
    source_conformer = source.GetConformer()
    target_conformer = target.GetConformer()
    for index, atom in enumerate(target.GetAtoms()):
        if atom.GetAtomicNum() != 1:
            target_conformer.SetAtomPosition(index, source_conformer.GetAtomPosition(index))


def _write_sdf(path: Path, molecule) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(path))
    try:
        writer.write(molecule)
    finally:
        writer.close()
    return path


def _write_mapped_pdbqt(
    path: Path,
    molecule,
    *,
    smiles: str | None = None,
    atom_types: tuple[str, ...] | None = None,
    shift_atom: int | None = None,
    shift: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> Path:
    """Write a small fixed-column Meeko-style mapped PDBQT fixture."""

    heavy = Chem.RemoveHs(molecule)
    smiles = smiles or Chem.MolToSmiles(heavy, canonical=False, isomericSmiles=True)
    atom_types = atom_types or tuple(
        {"C": "C", "N": "NA", "O": "OA", "S": "SA"}.get(atom.GetSymbol(), atom.GetSymbol())
        for atom in heavy.GetAtoms()
    )
    assert len(atom_types) == heavy.GetNumAtoms()
    conformer = molecule.GetConformer()
    lines = [
        f"REMARK SMILES {smiles}",
        "REMARK SMILES IDX " + " ".join(
            f"{index} {index}" for index in range(1, heavy.GetNumAtoms() + 1)
        ),
        "ROOT",
    ]
    for index, atom in enumerate(heavy.GetAtoms()):
        point = conformer.GetAtomPosition(index)
        xyz = [float(point.x), float(point.y), float(point.z)]
        if index == shift_atom:
            xyz = [value + delta for value, delta in zip(xyz, shift)]
        base = (
            f"{'HETATM':<6}{index + 1:5d} {('A' + str(index + 1))[:4]:>4} "
            f"{'LIG':>3} {'A':1}{1:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.00:6.2f}{0.00:6.2f}"
            f"          {atom.GetSymbol():>2}{'':>2}"
        )
        # Exercise the same fixed-column coordinate parser used for PDB/PDBQT
        # records before appending the AutoDock atom type field.
        fixed_coordinates(base)
        lines.append(base[:66] + f"    {0.000:6.3f} {atom_types[index]}")
    lines.extend(["ENDROOT", "TORSDOF 0"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_chiral_protonation_only_change_is_mapped_and_coordinates_are_preserved():
    source = _embedded("C[C@H](N)CC")
    normalized = _embedded("C[C@H]([NH3+])CC")
    _copy_heavy_coordinates(source, normalized)

    assert "@" in Chem.MolToSmiles(Chem.RemoveHs(source), isomericSmiles=True)

    ledger = build_identity_ledger(
        source,
        normalized,
        allow_hydrogen_changes=True,
        allow_formal_charge_changes=True,
        preserve_coordinates=True,
    )

    assert ledger["is_valid"] is True
    assert ledger["mapping_count"] >= 1
    assert ledger["differences"]["formal_charges"]
    assert ledger["differences"]["hydrogens"]
    assert all(item["delta"] == 1 for item in ledger["differences"]["formal_charges"])
    assert all(item["delta"] == 1 for item in ledger["differences"]["hydrogens"])


def test_mapped_pdbqt_exact_contract_accepts_explicit_h_prepared_input(tmp_path):
    prepared = _embedded("CCO")
    output = _write_mapped_pdbqt(tmp_path / "mapped.pdbqt", prepared)
    result = validate_mapped_pdbqt(prepared, output)
    assert result["status"] == "exact"
    assert result["is_valid"] is True


@pytest.mark.parametrize(
    "smiles,shift_atom",
    [("CC[O-]", None), ("CCO", 0)],
)
def test_mapped_pdbqt_changed_state_or_coordinates_is_rejected(tmp_path, smiles, shift_atom):
    prepared = _embedded("CCO")
    output = _write_mapped_pdbqt(
        tmp_path / ("changed_state.pdbqt" if shift_atom is None else "shifted.pdbqt"),
        prepared,
        smiles=smiles,
        shift_atom=shift_atom,
        shift=(0.02, 0.0, 0.0),
    )
    result = validate_mapped_pdbqt(prepared, output)
    assert result["status"] == "invalid"
    assert result["is_valid"] is False
    assert result["errors"]


def test_mapped_pdbqt_requires_finite_3d_prepared_input(tmp_path):
    source = Chem.MolFromSmiles("CCO")
    output = _write_mapped_pdbqt(tmp_path / "mapped_2d_input.pdbqt", _embedded("CCO"))
    result = validate_mapped_pdbqt(source, output)
    assert result["status"] == "invalid"
    assert result["scope"] == "input_coordinates_required"
    assert any("3D" in str(error) for error in result["errors"])


def test_rdkit_normalization_creates_finite_3d_for_2d_source(tmp_path):
    source_molecule = Chem.MolFromSmiles("CCO")
    source = _write_sdf(tmp_path / "source_2d.sdf", source_molecule)
    output = tmp_path / "normalized.sdf"
    result = preparation._normalize_with_rdkit(source, output)
    normalized = preparation._load_rdkit_molecule(output)
    require_finite_3d_coordinates(normalized, context="test normalized ligand")
    assert result["generated_3d"] is True
    assert result["coordinate_provenance"] == "rdkit_generated_seeded"


def test_failed_backend_does_not_publish_or_replace_preexisting_artifacts(tmp_path, monkeypatch):
    source_molecule = _embedded("CCO")
    source = _write_sdf(tmp_path / "source.sdf", source_molecule)
    output = tmp_path / "prepared.pdbqt"
    normalized = output.with_suffix(".normalized.sdf")
    output.write_bytes(b"prior pdbqt artifact")
    normalized.write_bytes(b"prior normalized artifact")

    def staged_normalizer(input_path, output_path, *, protonation_ph):
        _write_sdf(output_path, source_molecule)
        return {
            "normalization_backend": "openbabel",
            "source_file": str(input_path),
            "normalized_sdf": str(output_path),
            "protonation_ph": protonation_ph,
            "protonation_applied": True,
            "normalization_provenance": {
                "procedure": "test_openbabel_pH_protonation",
                "allowed_changes": ["hydrogens", "formal_charge"],
            },
        }

    monkeypatch.setenv("PDBWIZARD_LIGAND_PREP_PROFILE", "openbabel_meeko")
    monkeypatch.setattr(preparation, "normalize_ligand_to_sdf", staged_normalizer)

    def failing_backend(*args, **kwargs):
        raise subprocess.CalledProcessError(1, ["mk_prepare_ligand.py"])

    monkeypatch.setattr(preparation, "_prepare_with_meeko", failing_backend)
    with pytest.raises(RuntimeError, match="OpenBabel->Meeko"):
        preparation.prepare_ligand_for_vina_family(source, output)
    assert output.read_bytes() == b"prior pdbqt artifact"
    assert normalized.read_bytes() == b"prior normalized artifact"


@pytest.mark.parametrize("mutation", ["bond_order", "isotope", "radical"])
def test_heavy_atom_graph_mutations_are_rejected(mutation):
    source = _embedded("CCOC")
    normalized = Chem.Mol(source)
    rw_mol = Chem.RWMol(normalized)
    if mutation == "bond_order":
        rw_mol.GetBondWithIdx(0).SetBondType(Chem.BondType.DOUBLE)
    elif mutation == "isotope":
        rw_mol.GetAtomWithIdx(0).SetIsotope(13)
    else:
        rw_mol.GetAtomWithIdx(0).SetNumRadicalElectrons(1)
    normalized = rw_mol.GetMol()
    ledger = build_identity_ledger(source, normalized)
    assert ledger["is_valid"] is False
    assert ledger["errors"]


def test_stereochemical_inversion_is_rejected_even_when_graph_size_matches():
    source = _embedded("C[C@H](O)F")
    inverted = _embedded("C[C@@H](O)F")
    _copy_heavy_coordinates(source, inverted)
    ledger = build_identity_ledger(source, inverted)
    assert ledger["is_valid"] is False
    assert any("stere" in str(error).lower() for error in ledger["errors"])


def test_non_self_inverse_atom_permutation_is_mapped_in_the_correct_direction():
    source = _embedded("C1CCOC1")
    # A five-cycle permutation catches source/query mapping inversion while
    # retaining exactly the same heavy graph and coordinates.
    normalized = Chem.RenumberAtoms(source, [1, 2, 3, 4, 0, *range(5, source.GetNumAtoms())])
    ledger = build_identity_ledger(source, normalized, preserve_coordinates=True)
    assert ledger["is_valid"] is True
    assert {item["source_index"] for item in ledger["atom_map"]} == {0, 1, 2, 3, 4}


def test_undefined_source_stereo_cannot_be_assigned_by_normalization():
    source = _embedded("CC(O)F")
    assigned = _embedded("C[C@H](O)F")
    _copy_heavy_coordinates(source, assigned)
    ledger = build_identity_ledger(source, assigned)
    assert ledger["is_valid"] is False
    assert ledger["errors"]


def test_disconnected_ordinary_ligand_is_rejected_before_backend(tmp_path, monkeypatch):
    disconnected = Chem.AddHs(Chem.MolFromSmiles("CC.O"))
    source = _write_sdf(tmp_path / "disconnected.sdf", disconnected)
    calls = []

    def backend(*args, **kwargs):
        calls.append((args, kwargs))
        raise AssertionError("backend must not run for disconnected ordinary input")

    monkeypatch.setattr(preparation, "_prepare_with_meeko", backend)
    monkeypatch.setenv("PDBWIZARD_LIGAND_PREP_PROFILE", "meeko_only")
    with pytest.raises(ValueError, match="connected component"):
        preparation.prepare_ligand_for_vina_family(source, tmp_path / "out.pdbqt")
    assert calls == []


def test_source_coordinates_changed_by_normalizer_are_rejected(tmp_path):
    source = _embedded("CCO")
    normalized = Chem.Mol(source)
    normalized.GetConformer().SetAtomPosition(0, (0.025, 0.0, 0.0))
    ledger = build_identity_ledger(source, normalized, preserve_coordinates=True)
    assert ledger["is_valid"] is False
    assert any("coordinates" in str(error) for error in ledger["errors"])


def test_charge_change_without_paired_proton_is_rejected():
    source = _embedded("CCO")
    normalized = Chem.Mol(source)
    normalized.GetAtomWithIdx(2).SetFormalCharge(-1)
    ledger = build_identity_ledger(
        source,
        normalized,
        allow_hydrogen_changes=True,
        allow_formal_charge_changes=True,
    )
    assert ledger["is_valid"] is False
    assert any("paired protonation" in str(error) for error in ledger["errors"])


def test_absent_and_malformed_pdbqt_mapping_have_distinct_statuses(tmp_path):
    molecule = _embedded("CCO")
    absent = tmp_path / "absent.pdbqt"
    absent.write_text(
        "ROOT\nHETATM    1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00      0.000 C\n"
        "ENDROOT\nTORSDOF 0\n",
        encoding="utf-8",
    )
    absent_result = validate_mapped_pdbqt(molecule, absent)
    assert absent_result["status"] == "limited"

    malformed = tmp_path / "malformed.pdbqt"
    malformed.write_text(
        "REMARK SMILES CCO\nREMARK SMILES IDX 1 bad\n" + absent.read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    malformed_result = validate_mapped_pdbqt(molecule, malformed)
    assert malformed_result["status"] == "invalid"
