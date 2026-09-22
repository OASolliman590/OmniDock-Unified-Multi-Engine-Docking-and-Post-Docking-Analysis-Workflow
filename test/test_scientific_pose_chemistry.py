"""Chemistry-aware pose metadata contracts requiring the declared RDKit extra."""

from pathlib import Path

import pytest

from rdkit import Chem
from post_docking_analysis.pose_geometry import PoseGeometryError, load_pose_molecule


def _atom(serial: int, atom_type: str, x: float = 0.0) -> str:
    return (
        f"HETATM{serial:5d}  C{serial:<2d} LIG A   1    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00     0.000 {atom_type}\n"
    )


def _model(number: int, *, smiles: str | None = None, mapping: str = "", atom_types: tuple[str, str] = ("C", "C"), x: float = 0.0) -> str:
    metadata = ""
    if smiles is not None:
        metadata += f"REMARK SMILES {smiles}\n"
    if mapping:
        metadata += f"REMARK SMILES IDX {mapping}\n"
    return (
        f"MODEL {number}\n"
        + metadata
        + _atom(1, atom_types[0], x)
        + _atom(2, atom_types[1], x + 1.5)
        + "ENDMDL\n"
    )


def test_global_metadata_is_used_for_each_coherent_model(tmp_path: Path) -> None:
    path = tmp_path / "global.pdbqt"
    path.write_text(
        "REMARK SMILES CC\n"
        "REMARK SMILES H PARENT 1 0\n"
        "REMARK SMILES IDX 1 1\n"
        "REMARK SMILES IDX 2 2\n"
        + _model(1, atom_types=("C", "C"))
        + _model(2, atom_types=("C", "C"), x=5.0),
        encoding="utf-8",
    )

    assert load_pose_molecule(path, 1).GetNumAtoms() == 2
    assert load_pose_molecule(path, 2).GetNumAtoms() == 2


def test_per_model_metadata_is_selected_without_borrowing_another_model(tmp_path: Path) -> None:
    path = tmp_path / "per_model.pdbqt"
    path.write_text(
        _model(1, smiles="CC", mapping="1 1 2 2", atom_types=("C", "C"))
        + _model(2, smiles="CO", mapping="1 1 2 2", atom_types=("C", "OA")),
        encoding="utf-8",
    )

    assert Chem.MolToSmiles(load_pose_molecule(path, 1)) == "CC"
    assert Chem.MolToSmiles(load_pose_molecule(path, 2)) == "CO"


def test_selected_model_partial_metadata_cannot_mix_with_global_mapping(tmp_path: Path) -> None:
    path = tmp_path / "partial.pdbqt"
    path.write_text(
        "REMARK SMILES CC\nREMARK SMILES IDX 1 1 2 2\n"
        + _model(1, atom_types=("C", "C"))
        + _model(2, smiles="CO", atom_types=("C", "OA")),
        encoding="utf-8",
    )

    with pytest.raises(PoseGeometryError, match="partial_atom_mapping_metadata"):
        load_pose_molecule(path, 2)


def test_selected_model_cannot_borrow_complete_mapping_from_another_model(tmp_path: Path) -> None:
    path = tmp_path / "cross_model.pdbqt"
    path.write_text(
        _model(1, smiles="CC", mapping="1 1 2 2", atom_types=("C", "C"))
        + _model(2, atom_types=("C", "C")),
        encoding="utf-8",
    )

    with pytest.raises(PoseGeometryError, match="missing_authoritative_atom_mapping"):
        load_pose_molecule(path, 2)


def test_conflicting_complete_global_and_selected_metadata_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "conflict.pdbqt"
    path.write_text(
        "REMARK SMILES CC\nREMARK SMILES IDX 1 1 2 2\n"
        + _model(1, smiles="CO", mapping="1 1 2 2", atom_types=("C", "OA")),
        encoding="utf-8",
    )

    with pytest.raises(PoseGeometryError, match="conflicting_metadata_scopes"):
        load_pose_molecule(path, 1)


def test_duplicate_selected_smiles_metadata_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.pdbqt"
    path.write_text(
        "MODEL 1\nREMARK SMILES CC\nREMARK SMILES CC\nREMARK SMILES IDX 1 1 2 2\n"
        + _atom(1, "C")
        + _atom(2, "C", 1.5)
        + "ENDMDL\n",
        encoding="utf-8",
    )

    with pytest.raises(PoseGeometryError, match="duplicate_smiles_metadata"):
        load_pose_molecule(path, 1)


def test_duplicate_selected_mapping_metadata_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "duplicate_mapping.pdbqt"
    path.write_text(
        "MODEL 1\nREMARK SMILES CC\nREMARK SMILES IDX 1 1 2 2 1 1\n"
        + _atom(1, "C")
        + _atom(2, "C", 1.5)
        + "ENDMDL\n",
        encoding="utf-8",
    )

    with pytest.raises(PoseGeometryError, match="ambiguous_atom_mapping"):
        load_pose_molecule(path, 1)


def test_wrapped_mapping_and_glue_atoms_remain_valid(tmp_path: Path) -> None:
    path = tmp_path / "wrapped_glue.pdbqt"
    path.write_text(
        "MODEL 1\nREMARK SMILES C1CCC1\n"
        "REMARK SMILES IDX 1 1 2 2\nREMARK SMILES IDX 3 3 4 4\n"
        + _atom(1, "CG0")
        + _atom(2, "CG0", 1.5)
        + _atom(3, "CG0", 3.0)
        + _atom(4, "CG0", 4.5)
        + _atom(5, "G0", 4.5)
        + "ENDMDL\n",
        encoding="utf-8",
    )

    molecule = load_pose_molecule(path, 1)
    assert molecule.GetNumAtoms() == 4
    assert Chem.MolToSmiles(molecule) == "C1CCC1"
