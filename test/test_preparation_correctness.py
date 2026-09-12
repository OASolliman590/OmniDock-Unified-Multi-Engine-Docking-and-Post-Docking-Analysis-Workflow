"""Scientific preparation regressions using real parsers and small offline fixtures."""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
import subprocess

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from core_pipeline import MolecularDockingPipeline, extract_residue_level_coordinates
from cli_pipeline import run_single_pdb_cli
from docking.preparation import ligand_preparation as preparation
from docking.preparation import receptor_preparation
from docking.preparation.asset_identity import pair_reference_metadata
from docking.preparation.excel_sites import load_site_catalog_from_summary
from docking.preparation.ligand_quality import sanitize_ligand_pdb_file, sanitize_prepared_ligand_pdbqt, validate_prepared_ligand_pdbqt
from docking.preparation.project_builder import DockingPreparationConfig, DockingProjectBuilder
from docking.preparation.receptor_quality import validate_receptor_file
from docking.preparation.structure_contract import write_selection


def atom(serial, name="C1", *, x=0.0, y=0.0, z=0.0, element="C", record="HETATM", resname="LIG", chain="L", resid=1, alt="", occupancy=1.0, charge=""):
    return f"{record:<6}{serial:5d} {name:>4}{alt:1}{resname:>3} {chain:1}{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{0.0:6.2f}          {element:>2}{charge:>2}"


def pdbqt(serial, kind="C", **kwargs):
    return atom(serial, **kwargs)[:66] + f"    {0.0:6.3f} {kind}"


def ligand_file(path, types=("C",)):
    path.write_text("ROOT\n" + "\n".join(pdbqt(i + 1, kind) for i, kind in enumerate(types)) + "\nENDROOT\nTORSDOF 0\n", encoding="utf-8")
    return path


def sdf(path):
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=17)
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    return path


def test_extraction_preserves_identity_charge_connectivity_and_shared_frame(tmp_path):
    source = tmp_path / "1ABC.pdb"
    source.write_text("\n".join([atom(5, "N1", element="N", charge="1+", occupancy=.6), atom(8, "C2", x=1.4), atom(12, "CA", x=4, record="ATOM", resname="ALA", chain="A"), "CONECT    5    8    8", "CONECT    8    5    5", "END"]) + "\n")
    pipeline = MolecularDockingPipeline(str(tmp_path))
    extracted = Path(pipeline.save_hetatm_as_pdb(str(source), "LIG", "L", 1, pdb_id="1ABC"))
    text = extracted.read_text()
    assert "CONECT    5    8    8" in text
    assert "1+" in text and "N1" in text and "  0.60" in text
    receptor = Path(pipeline.clean_pdb(str(source), [], pdb_id="1ABC", remove_instances=[("L", "1", "", "LIG")]))
    assert "LIG" not in receptor.read_text()
    ligand_meta = json.loads(extracted.with_suffix(".pdb.preparation.json").read_text())
    receptor_meta = json.loads(receptor.with_suffix(".pdb.preparation.json").read_text())
    assert ligand_meta["reference_frame_id"] == receptor_meta["receptor_frame_id"]
    assert ligand_meta["atom_map"][0]["source_serial"] == 5
    assert ligand_meta["reference_pdb_id"] == "1ABC"


def test_altloc_requires_one_complete_conformer_and_keeps_occupancy(tmp_path):
    source, target = tmp_path / "alt.pdb", tmp_path / "selected.pdb"
    source.write_text("\n".join([atom(1, alt="A", x=0, occupancy=.6), atom(2, alt="B", x=10, occupancy=.4), atom(3, name="C2", alt="B", x=11)]) + "\n")
    with pytest.raises(ValueError, match="Incomplete alternate"):
        sanitize_ligand_pdb_file(source, target)
    assert not target.exists()
    source.write_text(source.read_text() + atom(4, name="C2", alt="A", x=1, occupancy=.6) + "\n")
    sanitize_ligand_pdb_file(source, target)
    selected = [line for line in target.read_text().splitlines() if line.startswith("HETATM")]
    assert [float(line[30:38]) for line in selected] == [0.0, 1.0]
    assert all(line[16] == " " and float(line[54:60]) == .6 for line in selected)


def test_multimodel_is_rejected_instead_of_merging(tmp_path):
    source = tmp_path / "models.pdb"
    source.write_text("MODEL        1\n" + atom(1) + "\nENDMDL\nMODEL        2\n" + atom(1, x=10) + "\nENDMDL\n")
    with pytest.raises(ValueError, match="Multiple models"):
        write_selection(source, tmp_path / "out.pdb")


def test_default_cli_removes_only_selected_ligand_instance(tmp_path, monkeypatch):
    source = tmp_path / "1ABC.pdb"
    source.write_text("\n".join([atom(1), atom(2, x=8, chain="M", resid=2), atom(3, record="ATOM", resname="ALA", chain="A", x=2), atom(4, element="O", resname="HOH", x=30)]) + "\n")
    pipeline = MolecularDockingPipeline(str(tmp_path))
    monkeypatch.setattr(pipeline, "fetch_pdb", lambda _id: str(source))
    assert run_single_pdb_cli(pipeline, "1ABC", {"ligand_selection": {"1abc": "LIG"}})
    lines = (tmp_path / "1ABC_cleaned.pdb").read_text().splitlines()
    assert not any(line.startswith("HETATM") and line[21] == "L" for line in lines)
    assert any(line.startswith("HETATM") and line[21] == "M" for line in lines)
    assert not any("HOH" in line for line in lines)


def test_site_envelope_contacts_are_order_independent_and_chain_specific(tmp_path):
    source = tmp_path / "site.pdb"
    ligand = [atom(1, x=-10), atom(2, name="C2", x=10)]
    protein = [atom(3, name="CA", x=11, record="ATOM", resname="ALA", chain="A"), atom(4, name="CB", x=12, record="ATOM", resname="ALA", chain="A"), atom(5, name="CA", x=-11, record="ATOM", resname="ALA", chain="B")]
    source.write_text("\n".join(ligand + protein) + "\nEND\n")
    first = extract_residue_level_coordinates(str(source), "LIG", "L", 1)
    source.write_text("\n".join(ligand + list(reversed(protein))) + "\nEND\n")
    second = extract_residue_level_coordinates(str(source), "LIG", "L", 1)
    assert first["num_interacting_residues"] == 2
    assert first["num_interacting_atoms"] == 3
    np.testing.assert_array_equal(first["overall_center"], second["overall_center"])
    assert first["size_x"] == 30
    assert set(first["residue_averages"]) == {"A:1:ALA", "B:1:ALA"}


def test_unsupported_pocket_metrics_are_not_fabricated(tmp_path):
    source = tmp_path / "empty.pdb"
    source.write_text("END\n")
    pipeline = MolecularDockingPipeline(str(tmp_path))
    values = pipeline.analyze_pocket_properties(str(source), np.zeros(3))
    assert values["pocket_volume_A3"] is None
    assert values["druggability_score"] is None
    assert values["druggability_interpretation"] == "Not evaluated"


@pytest.mark.parametrize("coordinate", [float("nan"), float("inf")])
def test_nonfinite_receptor_is_rejected(tmp_path, coordinate):
    receptor = tmp_path / "receptor.pdb"
    receptor.write_text("\n".join(atom(i, x=coordinate) for i in range(1, 101)))
    assert "invalid_coordinates" in {issue.reason for issue in validate_receptor_file(receptor)}


def test_macrocycle_closure_types_are_preserved_and_engine_qualified(tmp_path):
    source = ligand_file(tmp_path / "ring.pdbqt", ("CG0", "CG0", "G0", "G0"))
    before = source.read_bytes()
    assert validate_prepared_ligand_pdbqt(source, ["vina"]) == []
    sanitize_prepared_ligand_pdbqt(source, source)
    assert source.read_bytes() == before
    assert "unsupported_macrocycle_protocol" in {issue.reason for issue in validate_prepared_ligand_pdbqt(source, ["autodock4"])}


def test_real_meeko_macrocycle_output_survives_validation_and_copy(tmp_path):
    meeko = pytest.importorskip("meeko", exc_type=ImportError)
    MoleculePreparation, PDBQTWriterLegacy = meeko.MoleculePreparation, meeko.PDBQTWriterLegacy
    molecule = Chem.AddHs(Chem.MolFromSmiles("C1CCCCCCCCCCC1"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=19) == 0
    setup = MoleculePreparation().prepare(molecule)[0]
    text, okay, error = PDBQTWriterLegacy.write_string(setup)
    assert okay, error
    assert "CG0" in text and "G0" in text
    source = tmp_path / "actual_macrocycle.pdbqt"
    source.write_text(text, encoding="utf-8")
    assert validate_prepared_ligand_pdbqt(source, ["vina"]) == []
    target = tmp_path / "validated_macrocycle.pdbqt"
    sanitize_prepared_ligand_pdbqt(source, target)
    assert source.read_bytes() == target.read_bytes()


def test_unknown_type_repair_cannot_delete_chemistry(tmp_path):
    source = ligand_file(tmp_path / "unknown.pdbqt", ("Unknown",))
    before = source.read_bytes()
    with pytest.raises(ValueError, match="regeneration"):
        sanitize_prepared_ligand_pdbqt(source, source)
    assert source.read_bytes() == before


def test_pdb_cannot_be_treated_as_authoritative_ligand_chemistry(tmp_path):
    source = tmp_path / "ligand.pdb"
    source.write_text(atom(1))
    with pytest.raises(ValueError, match="authoritative SDF"):
        preparation.prepare_ligand_for_vina_family(source, tmp_path / "out.pdbqt")


def test_ph_failure_does_not_fall_back_to_add_hydrogens(tmp_path, monkeypatch):
    source = sdf(tmp_path / "ethanol.sdf")
    calls = []
    monkeypatch.setattr(preparation.shutil, "which", lambda command: command)
    def command(args):
        calls.append(args)
        if "-p" in args:
            raise subprocess.CalledProcessError(1, args, stderr="failed protonation")
        Path(args[args.index("-O") + 1]).write_bytes(Path(args[1]).read_bytes())
        return subprocess.CompletedProcess(args, 0, "", "1 molecule converted")
    monkeypatch.setattr(preparation, "_run_command", command)
    with pytest.raises(RuntimeError, match="pH-dependent"):
        preparation.normalize_ligand_to_sdf(source, tmp_path / "out.sdf", protonation_ph=5.0)
    assert not any("-h" in args for args in calls)
    assert not any("--gen3d" in args for args in calls)


def test_missing_openbabel_does_not_substitute_untreated_rdkit_state(tmp_path, monkeypatch):
    source = sdf(tmp_path / "ethanol.sdf")
    monkeypatch.setattr(preparation.shutil, "which", lambda _command: None)
    with pytest.raises(FileNotFoundError, match="Open Babel"):
        preparation.normalize_ligand_to_sdf(source, tmp_path / "out.sdf")


def test_receptor_backend_cannot_silently_drop_heavy_atoms(tmp_path, monkeypatch):
    source, output = tmp_path / "receptor.pdb", tmp_path / "receptor.pdbqt"
    source.write_text(atom(1, x=0) + "\n" + atom(2, x=2) + "\n")
    def run(args, **kwargs):
        Path(args[args.index("-p") + 1]).write_text(pdbqt(1) + "\n")
        return subprocess.CompletedProcess(args, 0, "dropped residue", "diagnostic")
    monkeypatch.setattr(receptor_preparation.subprocess, "run", run)
    with pytest.raises(ValueError, match="conservation failed"):
        receptor_preparation.prepare_receptor(source, output, {"preparation": {}})
    assert not output.exists()
    assert "diagnostic" in output.with_suffix(".pdbqt.preparation.log").read_text()


def test_strict_receptor_propagates_original_coordinate_frame(tmp_path, monkeypatch):
    source, output = tmp_path / "receptor.pdb", tmp_path / "receptor.pdbqt"
    source.write_text(atom(1) + "\n")
    source.with_suffix(".pdb.preparation.json").write_text(json.dumps({"receptor_frame_id": "coordinates:original"}))
    def run(args, **kwargs):
        Path(args[args.index("-p") + 1]).write_text(pdbqt(1) + "\n")
        return subprocess.CompletedProcess(args, 0, "", "")
    monkeypatch.setattr(receptor_preparation.subprocess, "run", run)
    report = receptor_preparation.prepare_receptor(source, output, {"preparation": {"ph": 6.5}})
    assert report["receptor_frame_id"] == "coordinates:original"
    assert report["protonation_status"] == "template_selected_not_ph_titrated"


def test_builder_selects_cleaned_receptor_and_rejects_ambiguous_conformations(tmp_path):
    for name in ["1ABC.pdbqt", "1ABC_cleaned.pdbqt"]:
        (tmp_path / name).write_text("prepared")
    builder = DockingProjectBuilder(SimpleNamespace())
    index = builder._index_assets(tmp_path, {".pdbqt"})
    assert builder._resolve_asset("1ABC", index, "receptor").name == "1ABC_cleaned.pdbqt"
    (tmp_path / "1ABC_conformation2.pdbqt").write_text("prepared")
    index = builder._index_assets(tmp_path, {".pdbqt"})
    with pytest.raises(ValueError, match="Ambiguous"):
        builder._resolve_asset("1ABC", index, "receptor")
    assert builder._resolve_asset("1ABC_cleaned.pdbqt", index, "receptor").name == "1ABC_cleaned.pdbqt"


def test_materialization_rejects_stale_copied_input(tmp_path):
    source, destination = tmp_path / "source", tmp_path / "copied"
    source.write_text("new chemistry")
    destination.write_text("old chemistry")
    builder = DockingProjectBuilder(SimpleNamespace(asset_mode="copy"))
    with pytest.raises(ValueError, match="Stale"):
        builder._materialize_asset(source, destination)


def test_reference_metadata_does_not_mislabel_cross_docking(tmp_path):
    receptor, ligand = tmp_path / "rec.pdbqt", tmp_path / "lig.pdbqt"
    receptor.with_suffix(".pdbqt.preparation.json").write_text(json.dumps({"receptor_frame_id": "frame2"}))
    ligand.with_suffix(".pdbqt.preparation.json").write_text(json.dumps({"reference_frame_id": "frame1", "reference_source": "experimental:1ABC", "reference_pose_file": "native.sdf"}))
    metadata = pair_reference_metadata(receptor, ligand)
    assert metadata["receptor_frame_id"] == "frame2"
    assert metadata["reference_pose_file"] == ""


def test_spreadsheet_rejects_duplicate_sites_and_nonfinite_centers(tmp_path):
    book = tmp_path / "sites.xlsx"
    rows = [{"PDB_ID": "1ABC", "Property": f"active_site_center_{axis}", "Value": 0} for axis in "xyz"]
    with pd.ExcelWriter(book) as writer:
        pd.DataFrame(rows + [rows[0]]).to_excel(writer, sheet_name="Summary", index=False)
    with pytest.raises(ValueError, match="Duplicate"):
        load_site_catalog_from_summary(book)
    rows[0]["Value"] = "inf"
    with pd.ExcelWriter(book) as writer:
        pd.DataFrame(rows).to_excel(writer, sheet_name="Summary", index=False)
    with pytest.raises(ValueError, match="Missing active site"):
        load_site_catalog_from_summary(book)


def test_project_build_preserves_verified_native_reference_and_frame(tmp_path):
    receptors, ligands = tmp_path / "receptors", tmp_path / "ligands"
    receptors.mkdir()
    ligands.mkdir()
    receptor = receptors / "1ABC_cleaned.pdbqt"
    receptor.write_text(pdbqt(1, record="ATOM", name="CA", resname="ALA") + "\n")
    ligand = ligand_file(ligands / "1ABC_ligand_LIG.pdbqt")
    native = sdf(tmp_path / "1ABC_native.sdf")
    reference = {
        "reference_source": "experimental:1ABC",
        "reference_pdb_id": "1ABC",
        "reference_pose_file": str(native),
        "reference_frame_id": "coordinates:original",
        "chemistry_status": "verified_instance_graph",
    }
    ligand.with_suffix(".pdbqt.preparation.json").write_text(json.dumps(reference))
    receptor.with_suffix(".pdbqt.preparation.json").write_text(json.dumps({"receptor_frame_id": "coordinates:original"}))
    pairlist = tmp_path / "input.csv"
    pd.DataFrame([dict(receptor=receptor.name, site_id="native", ligand=ligand.name,
                       center_x=0, center_y=0, center_z=0, size_x=20, size_y=20, size_z=20)]).to_csv(pairlist, index=False)
    result = DockingProjectBuilder(DockingPreparationConfig(
        prepared_proteins=receptors, prepared_ligands=ligands,
        pair_intent=None, excel_path=None, output_dir=tmp_path / "project",
        pairlist_file=pairlist, asset_mode="copy", engines=["vina"],
    )).build()
    row = pd.read_csv(result["pairlist_file"]).iloc[0]
    assert row["reference_source"] == "experimental:1ABC"
    assert row["reference_pose_file"] == str(native)
    assert row["reference_frame_id"] == row["receptor_frame_id"] == "coordinates:original"
    assert row["reference_pdb_id"] == "1ABC"
    materialized = Path(result["docking_root"]) / "ligands" / (ligand.name + ".preparation.json")
    assert json.loads(materialized.read_text())["chemistry_status"] == "verified_instance_graph"
