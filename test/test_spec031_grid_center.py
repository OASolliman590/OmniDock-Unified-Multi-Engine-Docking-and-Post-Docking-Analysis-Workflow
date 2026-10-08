from pathlib import Path

import numpy as np

from docking.preparation.binding_site_center import CENTER_METHOD, resolve_binding_site_center


def _atom(serial: int, name: str, altloc: str, resname: str, chain: str, resnum: int,
          x: float, occupancy: float, element: str) -> str:
    return (
        f"HETATM{serial:5d} {name:>4}{altloc:1}{resname:>3} {chain:1}{resnum:4d}    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{occupancy:6.2f}{20.0:6.2f}          {element:>2}\n"
    )


def test_holo_center_is_selected_ligand_heavy_atom_centroid_with_altloc_policy(tmp_path: Path) -> None:
    pdb = tmp_path / "holo.pdb"
    pdb.write_text(
        _atom(1, "C1", "A", "LIG", "A", 101, 0.0, 0.40, "C")
        + _atom(2, "C1", "B", "LIG", "A", 101, 10.0, 0.60, "C")
        + _atom(3, "C2", "A", "LIG", "A", 101, 2.0, 0.50, "C")
        + _atom(4, "C2", "B", "LIG", "A", 101, 20.0, 0.50, "C")
        + _atom(5, "O1", " ", "LIG", "A", 101, 4.0, 1.00, "O")
        + _atom(6, "H1", " ", "LIG", "A", 101, 100.0, 1.00, "H")
        + _atom(7, "C1", " ", "LIG", "B", 202, 999.0, 1.00, "C")
        + "END\n",
        encoding="utf-8",
    )
    result = resolve_binding_site_center(
        structure_file=pdb, accession="TEST", ligand_name="LIG", chain_id="A", residue_number=101
    )
    assert result["status"] == "completed"
    assert result["center_method"] == CENTER_METHOD
    np.testing.assert_allclose(result["binding_site_center"], [16.0 / 3.0, 0.0, 0.0])
    np.testing.assert_allclose(result["overall_center"], result["binding_site_center"])
    np.testing.assert_allclose(result["ligand_center"], result["binding_site_center"])
    assert result["heavy_atom_count"] == 3
    assert result["source_structure_sha256"]


def test_no_ligand_selection_does_not_choose_first_hetero_residue(tmp_path: Path) -> None:
    pdb = tmp_path / "two_ligands.pdb"
    pdb.write_text(
        _atom(1, "C1", " ", "AAA", "A", 1, 1.0, 1.0, "C")
        + _atom(2, "C1", " ", "BBB", "A", 2, 2.0, 1.0, "C")
        + "END\n",
        encoding="utf-8",
    )
    result = resolve_binding_site_center(structure_file=pdb)
    assert result["status"] == "skipped_missing_configuration"
    assert "binding_site_center" not in result


def test_apo_requires_explicit_coordinates_and_records_source() -> None:
    missing = resolve_binding_site_center()
    assert missing["status"] == "skipped_missing_configuration"
    explicit = resolve_binding_site_center(
        explicit_center=(1.0, 2.0, 3.0), explicit_center_source="Scientific Lead box record"
    )
    assert explicit["status"] == "completed"
    assert explicit["reason"] == "explicit_user_coordinates"
    assert explicit["explicit_center_source"] == "Scientific Lead box record"
    np.testing.assert_allclose(explicit["binding_site_center"], [1.0, 2.0, 3.0])

