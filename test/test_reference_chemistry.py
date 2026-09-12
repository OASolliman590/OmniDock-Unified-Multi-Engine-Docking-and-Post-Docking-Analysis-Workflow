import json

import pytest

Chem = pytest.importorskip("rdkit.Chem")
from workflow.reference_chemistry import retrieve_reference_sdf


def reference_fixture(tmp_path, *, shift=0.0):
    molecule = Chem.MolFromSmiles("CCO")
    conf = Chem.Conformer(3)
    positions = [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (2.9, 0.2, 0.0)]
    lines = []
    for i, (atom, position) in enumerate(zip(molecule.GetAtoms(), positions), 1):
        conf.SetAtomPosition(i - 1, (position[0] + shift, position[1], position[2]))
        element = atom.GetSymbol()
        name = f"{element}{i}"
        x, y, z = position
        lines.append(f"HETATM{i:5d} {name:>4} LIG A 101    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 10.00          {element:>2}  ")
    molecule.AddConformer(conf)
    pdb = tmp_path / "1abc_ligand_LIG_A_101.pdb"
    pdb.write_text("\n".join(lines) + "\nEND\n")
    pdb.with_suffix(".pdb.preparation.json").write_text(json.dumps({"reference_source": "cocrystal", "reference_pdb_id": "1ABC", "reference_frame_id": "native-frame"}))
    sdf = Chem.MolToMolBlock(molecule) + "\n$$$$\n"
    return pdb, sdf


def test_authoritative_instance_retains_native_coordinates_and_graph(tmp_path):
    pdb, sdf = reference_fixture(tmp_path)
    output = tmp_path / "reference.sdf"
    urls = []

    def fetch(url):
        urls.append(url)
        return sdf

    result = retrieve_reference_sdf(pdb, "1ABC", output, fetch_text=fetch)
    molecule = Chem.SDMolSupplier(str(output), removeHs=False)[0]
    assert Chem.MolToSmiles(molecule) == "CCO"
    assert molecule.GetConformer().GetAtomPosition(1).x == 1.5
    assert result["reference_frame_id"] == "native-frame"
    assert result["chemistry_status"] == "verified_instance_graph"
    assert len(result["sdf_atom_map"]) == 3
    assert "auth_asym_id=A&auth_seq_id=101" in urls[0]


def test_instance_coordinates_are_never_aligned_to_hide_mismatch(tmp_path):
    pdb, sdf = reference_fixture(tmp_path, shift=100)
    output = tmp_path / "reference.sdf"
    with pytest.raises(ValueError, match="native coordinates"):
        retrieve_reference_sdf(pdb, "1ABC", output, fetch_text=lambda _: sdf)
    assert not output.exists()


def test_generated_conformer_cannot_become_experimental_reference(tmp_path):
    pdb, sdf = reference_fixture(tmp_path)
    pdb.with_suffix(".pdb.preparation.json").write_text("{}")
    with pytest.raises(ValueError, match="provenance"):
        retrieve_reference_sdf(pdb, "1ABC", tmp_path / "reference.sdf", fetch_text=lambda _: sdf)


def test_multiple_instance_records_are_rejected(tmp_path):
    pdb, sdf = reference_fixture(tmp_path)
    with pytest.raises(ValueError, match="one valid molecular graph"):
        retrieve_reference_sdf(pdb, "1ABC", tmp_path / "reference.sdf", fetch_text=lambda _: sdf + sdf)
