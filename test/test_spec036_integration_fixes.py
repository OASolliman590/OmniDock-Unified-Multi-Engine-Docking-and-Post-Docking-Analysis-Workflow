"""Spec 036 integration blockers found by the real 1IEP re-dock (regression tests, one group per blocker).

1. `prep project --pairlist-file` keeps every pairlist column; the redocking reference lineage is read from
   metadata/pair_intent.csv when the pairlist lacks it, and is keyed by the R6 pair tag.
2. `prep project` keeps the project protonation_policy and the engine settings.
3. The pair-intent snapshot is always metadata/pair_intent.csv.
4. (md-inputs layout resolution) lives in test_spec036_protonation_policy.py (canonical and legacy).
5. One lineage comparison entry point with the aromatic normalisation; charge-blind parent mapping for redocking.
6. `main.py dock run` accepts and forwards --replicates and --energy-range.
7. num_modes: basic mode refuses an explicit --num-modes; the effective source is recorded.
8. Replicate jobs record wall_time_s and duration_s.
9. `analyze md-inputs` prints the status and the reason when blocked.
"""

from __future__ import annotations

import json
import stat
from pathlib import Path

import pandas as pd
import pytest

from docking.models import PairlistRow
from docking.parameter_schema import resolve_parameter_schema, validate_parameter_schema
from docking.preparation.project_builder import DockingPreparationConfig, DockingProjectBuilder
from docking.project_layout import (
    PAIRLIST_COLUMNS,
    load_protonation_policy,
    set_protonation_policy,
)
from docking.runners.vina import VinaRunner
from post_docking_analysis.atom_mapping import (
    IN_PLACE_METHOD,
    compare_lineage_graph_poses,
    compare_pose_pair,
    load_pdbqt_lineage_pose,
    load_sdf_graph_pose,
)
from post_docking_analysis.reference_policy import classify_reference_row
from post_docking_analysis.redocking_validation import _compute_pose_rmsd, run_redocking_validation
from post_docking_analysis.replicates import pose_reproducibility_table
from workflow import execution

from test_spec036_docking_params import _dock_args, _dock_project, _hexagon_atoms, _pdbqt_text, _run_manifest
from test_spec036_protonation_policy import _auto_project, _prepared_dirs_for

REPO_ROOT = Path(__file__).resolve().parents[1]
MODEL1_PDBQT = REPO_ROOT / "test" / "fixtures" / "spec034" / "sti_vina_model1.pdbqt"
CRYSTAL_SDF = REPO_ROOT / "test" / "fixtures" / "spec032" / "phase2_1iep_inputs" / "STI_1IEP_A201.sdf"
PLUS_ONE_SMILES = "Cc1ccc(NC(=O)c2ccc(CN3CC[NH+](C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"
NEUTRAL_SMILES = "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"
RECEPTOR = "1IEP_A_protein.pdbqt"
LIGAND = "STI_A_201.pdbqt"
PAIR_TAG = "1IEP_A_protein_site_1_STI_A_201.pdbqt"
LEGACY_TAG = "1IEP_A_protein.pdbqt_site_1_STI_A_201.pdbqt"


def _plus_one_copy(tmp_path: Path) -> Path:
    """The committed STI model 1 with its lineage SMILES in the +1 N-methylpiperazinium form."""
    text = MODEL1_PDBQT.read_text(encoding="utf-8")
    old = f"REMARK SMILES {NEUTRAL_SMILES}"
    assert text.count(old) == 1
    target = tmp_path / "sti_vina_model1_plus1.pdbqt"
    target.write_text(text.replace(old, f"REMARK SMILES {PLUS_ONE_SMILES}"), encoding="utf-8")
    return target


# ----------------------------------------------------------------------------- blocker 1: pairlist columns


_CARRIED_ROW = {
    "receptor": RECEPTOR,
    "site_id": "site_1",
    "ligand": LIGAND,
    "center_x": 15.6,
    "center_y": 53.4,
    "center_z": 15.5,
    "size_x": 19.4,
    "size_y": 19.4,
    "size_z": 19.4,
    "protein_display_name": "A",
    "ligand_display_name": "STI A 201",
    "pdb_id": "1IEP",
    "pair_source": "cocrystal",
    "is_cocrystal_benchmark": True,
    "cocrystal_ligand_name": "STI_A_201",
    "cocrystal_ligand_display_name": "STI A 201",
    "selection_mode": "cocrystal_only",
    "box_method": "rg_scaled_v1",
    "ligand_rg_angstrom": 6.69,
    "edge_angstrom": 19.4,
    "box_containment_status": "contained",
    "box_warnings": "",
    "reference_ligand_file": "/refs/STI_A_201.sdf",
    "reference_receptor_file": "/refs/1IEP.pdb",
    "reference_source_structure_id": "1IEP",
}


def _valid_ligand_pdbqt() -> str:
    """A prepared-ligand PDBQT with the ROOT and TORSDOF records that the materialisation check requires."""
    body = _pdbqt_text(_hexagon_atoms()).splitlines()
    return "\n".join(["ROOT", *body, "ENDROOT", "TORSDOF 0", ""])


def _prep_inputs(tmp_path: Path, *, pairlist_row: dict, intent_rows: list, intent_name: str = "pair_intent.csv"):
    proteins = tmp_path / "prepared_proteins"
    ligands = tmp_path / "prepared_ligands"
    proteins.mkdir(parents=True)
    ligands.mkdir(parents=True)
    (proteins / RECEPTOR).write_text("ATOM\n", encoding="utf-8")
    (ligands / LIGAND).write_text(_valid_ligand_pdbqt(), encoding="utf-8")
    pairlist = tmp_path / "pairlist_in.csv"
    pd.DataFrame([pairlist_row]).to_csv(pairlist, index=False)
    intent = tmp_path / intent_name
    pd.DataFrame(intent_rows).to_csv(intent, index=False)
    return proteins, ligands, pairlist, intent


def _run_prep_project(tmp_path: Path, project: Path, *, pairlist_row: dict, intent_rows: list, intent_name: str = "pair_intent.csv"):
    proteins, ligands, pairlist, intent = _prep_inputs(
        tmp_path / "inputs", pairlist_row=pairlist_row, intent_rows=intent_rows, intent_name=intent_name
    )
    config = DockingPreparationConfig(
        prepared_proteins=proteins,
        prepared_ligands=ligands,
        pair_intent=intent,
        excel_path=None,
        output_dir=project,
        pairlist_file=pairlist,
        engines=["vina"],
        layout_profile="canonical",
        asset_mode="copy",
    )
    DockingProjectBuilder(config).build()
    return project


def test_prep_project_keeps_every_existing_pairlist_column(tmp_path):
    project = tmp_path / "project"
    execution.run_workflow_init(str(project), layout_profile="canonical")
    intent_rows = [{"receptor": RECEPTOR, "ligand": LIGAND, "site_id": "site_1", "is_cocrystal_benchmark": True}]
    _run_prep_project(tmp_path, project, pairlist_row=dict(_CARRIED_ROW), intent_rows=intent_rows)

    written = pd.read_csv(project / "pairlist.csv")
    source_columns = list(_CARRIED_ROW)
    assert set(source_columns) <= set(written.columns), set(source_columns) - set(written.columns)
    assert set(PAIRLIST_COLUMNS) <= set(written.columns)
    row = written.iloc[0]
    assert str(row["pdb_id"]) == "1IEP"
    assert str(row["cocrystal_ligand_name"]) == "STI_A_201"
    assert str(row["is_cocrystal_benchmark"]).lower() == "true"
    assert str(row["reference_ligand_file"]) == "/refs/STI_A_201.sdf"
    assert str(row["reference_source_structure_id"]) == "1IEP"
    # No float NaN is written for blank cells.
    assert "nan" not in (project / "pairlist.csv").read_text(encoding="utf-8").lower()


def test_prep_project_pairlist_union_keeps_blank_metadata_blank(tmp_path):
    project = tmp_path / "project"
    execution.run_workflow_init(str(project), layout_profile="canonical")
    stripped = {key: value for key, value in _CARRIED_ROW.items() if key in PAIRLIST_COLUMNS}
    stripped["is_cocrystal_benchmark"] = ""
    _run_prep_project(tmp_path, project, pairlist_row=stripped, intent_rows=[{"receptor": RECEPTOR, "ligand": LIGAND}])
    written = pd.read_csv(project / "pairlist.csv")
    assert list(written.columns) == list(PAIRLIST_COLUMNS)
    assert pd.isna(written.loc[0, "is_cocrystal_benchmark"])


# ----------------------------------------------------------------------------- blocker 2: protonation policy


def test_prep_project_keeps_protonation_policy_and_engine_settings(tmp_path):
    project = tmp_path / "project"
    execution.run_workflow_init(str(project), layout_profile="canonical")
    state_map = tmp_path / "ligand_state_map.csv"
    state_map.write_text("ligand,net_charge\nSTI_A_201,1\n", encoding="utf-8")
    block = set_protonation_policy(
        project,
        receptor_ph=7.4,
        receptor_force_field="AMBER",
        ligand_policy="explicit_state",
        ligand_state_map=state_map,
    )
    manifest_file = project / "project_manifest.json"
    payload = json.loads(manifest_file.read_text(encoding="utf-8"))
    payload["engine_settings"] = {"vina": {"version": "1.2.5", "binary": "/bin/vina"}}
    payload["created_at"] = "2026-01-01T00:00:00+00:00"
    manifest_file.write_text(json.dumps(payload), encoding="utf-8")

    intent_rows = [{"receptor": RECEPTOR, "ligand": LIGAND}]
    _run_prep_project(tmp_path, project, pairlist_row=dict(_CARRIED_ROW), intent_rows=intent_rows)

    after = load_protonation_policy(project)
    assert after is not None and after["ligand_policy"] == "explicit_state"
    assert after["receptor_ph"] == block["receptor_ph"] == 7.4
    assert after["ligand_state_map"] == block["ligand_state_map"]
    kept = json.loads(manifest_file.read_text(encoding="utf-8"))
    assert kept["engine_settings"]["vina"]["version"] == "1.2.5"
    assert kept["created_at"] == "2026-01-01T00:00:00+00:00"


# ----------------------------------------------------------------------------- blocker 3: snapshot name


def test_pair_intent_snapshot_is_always_named_pair_intent_csv(tmp_path):
    project = tmp_path / "project"
    execution.run_workflow_init(str(project), layout_profile="canonical")
    intent_rows = [{"receptor": RECEPTOR, "ligand": LIGAND, "is_cocrystal_benchmark": True}]
    _run_prep_project(
        tmp_path, project, pairlist_row=dict(_CARRIED_ROW), intent_rows=intent_rows, intent_name="my_reviewed_intent.csv"
    )
    snapshot = project / "metadata" / "pair_intent.csv"
    assert snapshot.is_file()
    assert not (project / "metadata" / "my_reviewed_intent.csv").exists()
    assert pd.read_csv(snapshot)["is_cocrystal_benchmark"].tolist() == [True]


# ----------------------------------------------------------------------------- blocker 1: reference lineage


def _pipeline_for(project: Path):
    from post_docking_analysis.multi_engine_pipeline_impl import MultiEngineAnalysisPipeline

    pipeline = MultiEngineAnalysisPipeline.__new__(MultiEngineAnalysisPipeline)
    pipeline.project_dir = project
    pipeline.manifest = {"layout_profile": "canonical"}
    return pipeline


def _pairlist_row(**overrides):
    row = {
        "receptor": RECEPTOR,
        "site_id": "site_1",
        "ligand": LIGAND,
        "center_x": 15.6,
        "center_y": 53.4,
        "center_z": 15.5,
        "size_x": 19.4,
        "size_y": 19.4,
        "size_z": 19.4,
    }
    row.update(overrides)
    return row


def _intent_row(**overrides):
    row = {
        "receptor": RECEPTOR,
        "site_id": "site_1",
        "ligand": LIGAND,
        "pdb_id": "1IEP",
        "pair_source": "cocrystal",
        "is_cocrystal_benchmark": True,
        "cocrystal_ligand_name": "STI_A_201",
        "reference_ligand_file": str(CRYSTAL_SDF),
        "reference_receptor_file": "/refs/1IEP.pdb",
        "reference_source_structure_id": "1IEP",
    }
    row.update(overrides)
    return row


def test_pair_metadata_is_keyed_by_the_r6_pair_tag_and_falls_back_to_pair_intent(tmp_path):
    project = tmp_path / "project"
    (project / "metadata").mkdir(parents=True)
    # The pairlist lost the metadata (the pre-fix `prep project` behaviour); pair_intent.csv has it.
    pd.DataFrame([_pairlist_row()]).to_csv(project / "pairlist.csv", index=False)
    pd.DataFrame([_intent_row()]).to_csv(project / "metadata" / "pair_intent.csv", index=False)

    meta = _pipeline_for(project)._pair_metadata_frame().set_index("tag")
    assert PAIR_TAG in meta.index and LEGACY_TAG in meta.index
    row = meta.loc[PAIR_TAG]
    assert bool(row["is_cocrystal_benchmark"]) is True
    assert row["cocrystal_ligand_name"] == "STI_A_201"
    assert row["reference_ligand_file"] == str(CRYSTAL_SDF)
    assert row["reference_receptor_file"] == "/refs/1IEP.pdb"
    assert row["reference_source_structure_id"] == "1IEP"


def test_cocrystal_benchmark_pair_is_a_reference_row_when_pairlist_lacks_lineage(tmp_path):
    project = tmp_path / "project"
    (project / "metadata").mkdir(parents=True)
    pd.DataFrame([_pairlist_row()]).to_csv(project / "pairlist.csv", index=False)
    pd.DataFrame([_intent_row()]).to_csv(project / "metadata" / "pair_intent.csv", index=False)

    pipeline = _pipeline_for(project)
    pose_row = {"protein": RECEPTOR, "ligand": LIGAND, "tag": PAIR_TAG, "engine": "vina", "pose_file": "/x.pdbqt"}
    merged = pipeline._merge_pair_metadata(pd.DataFrame([pose_row]), pipeline._pair_metadata_frame())
    record = merged.iloc[0].to_dict()
    classification = classify_reference_row(record)
    assert classification.is_reference is True
    assert classification.classification == "reference"
    assert classification.validation_anchor_eligible is True


def test_pairlist_without_intent_still_yields_non_reference(tmp_path):
    project = tmp_path / "project"
    project.mkdir()
    pd.DataFrame([_pairlist_row()]).to_csv(project / "pairlist.csv", index=False)
    pipeline = _pipeline_for(project)
    merged = pipeline._merge_pair_metadata(
        pd.DataFrame([{"protein": RECEPTOR, "ligand": LIGAND, "tag": PAIR_TAG, "engine": "vina"}]),
        pipeline._pair_metadata_frame(),
    )
    assert classify_reference_row(merged.iloc[0].to_dict()).classification == "non_reference"


# ----------------------------------------------------------------------------- blocker 5: lineage comparison


def test_committed_sti_model1_in_place_rmsd_through_the_shared_entry_point():
    docked = load_pdbqt_lineage_pose(MODEL1_PDBQT, 1)
    crystal = load_sdf_graph_pose(CRYSTAL_SDF, 1)
    # The crystal SDF is Kekule and the lineage is aromatic: the shared normalisation makes them comparable.
    result = compare_pose_pair(docked, crystal, rmsd_frame=IN_PLACE_METHOD)
    assert result.comparable, result.reason
    assert result.rmsd_angstrom == pytest.approx(0.857, abs=0.005)
    assert result.charge_blind_parent_mapping is False
    assert result.to_dict()["charge_blind_parent_mapping"] is False


def test_plus_one_copy_matches_the_crystal_only_with_charge_blind_parent_mapping(tmp_path):
    plus_one = _plus_one_copy(tmp_path)

    # Without the convention the +1 charge on the piperazine nitrogen is a node mismatch.
    strict = compare_lineage_graph_poses(
        load_pdbqt_lineage_pose(plus_one, 1), load_sdf_graph_pose(CRYSTAL_SDF, 1), rmsd_frame=IN_PLACE_METHOD
    )
    assert strict.comparable is False
    assert strict.reason == "graphs_not_isomorphic"

    rmsd, error, details = _compute_pose_rmsd(plus_one, CRYSTAL_SDF, docked_pose_index=1, reference_pose_index=1)
    assert error == ""
    assert rmsd == pytest.approx(0.857, abs=0.005)
    assert details["charge_blind_parent_mapping"] is True
    assert details["rmsd_frame"] == IN_PLACE_METHOD


def test_neutral_fixture_redocking_is_also_recorded_as_charge_blind(tmp_path):
    rmsd, error, details = _compute_pose_rmsd(MODEL1_PDBQT, CRYSTAL_SDF, docked_pose_index=1, reference_pose_index=1)
    assert error == ""
    assert rmsd == pytest.approx(0.857, abs=0.005)
    assert details["charge_blind_parent_mapping"] is True


def test_redocking_validation_row_records_charge_blind_mapping(tmp_path):
    plus_one = _plus_one_copy(tmp_path)
    best = pd.DataFrame(
        [
            {
                "protein": RECEPTOR,
                "ligand": LIGAND,
                "tag": PAIR_TAG,
                "engine": "vina",
                "pose_file": str(plus_one),
                "pose": 1,
                "affinity_kcal_mol": -13.2,
                "is_cocrystal_benchmark": True,
                "cocrystal_ligand_name": "STI_A_201",
                "pdb_id": "1IEP",
                "reference_ligand_file": str(CRYSTAL_SDF),
                "reference_receptor_file": "/refs/1IEP.pdb",
                "reference_source_structure_id": "1IEP",
            }
        ]
    )
    outputs = run_redocking_validation(project_dir=tmp_path, best_by_engine=best, output_dir=tmp_path / "reports")
    row = outputs["validation_df"].iloc[0]
    assert row["redocking_classification"] == "pass"
    assert bool(row["charge_blind_parent_mapping"]) is True
    assert float(row["redocking_rmsd_angstrom"]) == pytest.approx(0.857, abs=0.005)


def test_pose_reproducibility_compares_through_the_same_normalised_entry_point(tmp_path):
    # Two replicate top poses: a Meeko lineage PDBQT (aromatic) and the Kekule crystal SDF.
    normalized = pd.DataFrame(
        [
            {"engine": "vina", "protein": RECEPTOR, "tag": PAIR_TAG, "replicate_id": 1, "seed": 1,
             "pose": 1, "pose_file": str(MODEL1_PDBQT), "affinity_kcal_mol": -13.2},
            {"engine": "vina", "protein": RECEPTOR, "tag": PAIR_TAG, "replicate_id": 2, "seed": 2,
             "pose": 1, "pose_file": str(CRYSTAL_SDF), "affinity_kcal_mol": -13.0},
        ]
    )
    table = pose_reproducibility_table(normalized)
    row = table.iloc[0]
    assert row["pose_reproducibility_status"] == "completed"
    assert float(row["pose_reproducibility_max_rmsd_angstrom"]) == pytest.approx(0.857, abs=0.005)


# ----------------------------------------------------------------------------- blocker 6: dock run flags


def test_workflow_dock_run_accepts_and_forwards_replicates_and_energy_range():
    from workflow.cli import _build_dock_argv, build_parser

    parser = build_parser()
    args = parser.parse_args(
        [
            "dock",
            "run",
            "--project-dir",
            "/tmp/none",
            "--engines",
            "vina",
            "--seed",
            "7",
            "--replicates",
            "2",
            "--energy-range",
            "4.5",
        ]
    )
    assert args.replicates == 2
    assert args.energy_range == pytest.approx(4.5)
    argv = _build_dock_argv(args)
    assert argv[argv.index("--replicates") + 1] == "2"
    assert argv[argv.index("--energy-range") + 1] == "4.5"


def test_workflow_dock_run_without_the_flags_leaves_the_docking_defaults(tmp_path):
    from workflow.cli import _build_dock_argv, build_parser

    args = build_parser().parse_args(["dock", "run", "--project-dir", str(tmp_path), "--seed", "7"])
    argv = _build_dock_argv(args)
    assert "--replicates" not in argv and "--energy-range" not in argv and "--num-modes" not in argv


# ----------------------------------------------------------------------------- blocker 7: num_modes


def test_basic_preset_sets_num_modes_and_records_the_preset_source():
    schema = resolve_parameter_schema(
        mode="basic", preset="exhaustive", exhaustiveness=32, num_modes=20, seed=1,
        box_scale=1.0, box_padding=0.0, runtime_by_engine={"vina": {}},
    )
    assert schema.common.num_modes == 40
    assert schema.common.num_modes_source == "basic_preset:exhaustive"
    assert validate_parameter_schema(schema, ["vina"])[0] == []


def test_basic_mode_refuses_an_explicit_num_modes_with_a_clear_message():
    schema = resolve_parameter_schema(
        mode="basic", preset="exhaustive", exhaustiveness=32, num_modes=5, seed=1,
        box_scale=1.0, box_padding=0.0, runtime_by_engine={"vina": {}}, num_modes_explicit=True,
    )
    errors, _warnings = validate_parameter_schema(schema, ["vina"])
    assert any("--num-modes" in error and "--parameter-mode advanced" in error for error in errors)


def test_advanced_mode_honours_an_explicit_num_modes():
    schema = resolve_parameter_schema(
        mode="advanced", preset="exhaustive", exhaustiveness=32, num_modes=5, seed=1,
        box_scale=1.0, box_padding=0.0, runtime_by_engine={"vina": {}}, num_modes_explicit=True,
    )
    assert schema.common.num_modes == 5
    assert schema.common.num_modes_source == "user_explicit"
    assert validate_parameter_schema(schema, ["vina"])[0] == []


def test_dock_run_cli_refuses_explicit_num_modes_in_basic_mode(tmp_path, monkeypatch, capsys):
    import docking.cli as docking_cli

    root = _dock_project(tmp_path, monkeypatch)
    assert docking_cli.dock_main(_dock_args(root, "--seed", "100", "--num-modes", "5")) == 1
    assert "--num-modes` is not used in basic mode" in capsys.readouterr().out


def test_dock_run_records_the_effective_num_modes_source(tmp_path, monkeypatch, capsys):
    import docking.cli as docking_cli

    root = _dock_project(tmp_path, monkeypatch)
    assert docking_cli.dock_main(_dock_args(root, "--seed", "100")) == 0
    out = capsys.readouterr().out
    assert "num_modes=40 (source: basic_preset:exhaustive)" in out
    manifest = _run_manifest(root)
    assert manifest["effective"]["num_modes"] == 40
    assert manifest["effective"]["num_modes_source"] == "basic_preset:exhaustive"
    assert manifest["jobs"][0]["num_modes"] == 40


# ----------------------------------------------------------------------------- blocker 8: wall time


def _fake_vina(tmp_path: Path) -> Path:
    fake = tmp_path / "bin" / "fake_vina"
    fake.parent.mkdir(exist_ok=True)
    fake.write_text(
        "#!/bin/bash\n"
        "out=\"\"\n"
        "while [ $# -gt 0 ]; do if [ \"$1\" = \"--out\" ]; then out=\"$2\"; fi; shift; done\n"
        "cat > \"$out\" <<'EOF'\n"
        "MODEL 1\nREMARK VINA RESULT:     -7.5      0.000      0.000\nENDMDL\n"
        "EOF\n",
        encoding="utf-8",
    )
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return fake


def test_replicate_jobs_record_wall_time_and_duration(tmp_path):
    fake = _fake_vina(tmp_path)
    row = PairlistRow(
        receptor=RECEPTOR, site_id="site_1", ligand=LIGAND,
        center_x=0.0, center_y=0.0, center_z=0.0, size_x=4.0, size_y=4.0, size_z=4.0,
    )
    runner = VinaRunner(
        tmp_path / "proj",
        runtime={"binary": str(fake), "exhaustiveness": 8, "num_modes": 1, "seed": 9, "replicates": 2,
                 "energy_range": 3.0},
    )
    payload = runner.run([row], dry_run=False)
    assert len(payload["jobs"]) == 2
    for job in payload["jobs"]:
        assert job["status"] == "completed"
        assert isinstance(job["wall_time_s"], float) and job["wall_time_s"] >= 0.0
        assert job["duration_s"] == job["wall_time_s"]
    on_disk = json.loads((tmp_path / "proj" / "engines" / "vina" / "run_manifest.json").read_text(encoding="utf-8"))
    assert all(job["wall_time_s"] is not None and job["duration_s"] is not None for job in on_disk["jobs"])


def test_dry_run_jobs_have_no_wall_time(tmp_path):
    row = PairlistRow(
        receptor=RECEPTOR, site_id="site_1", ligand=LIGAND,
        center_x=0.0, center_y=0.0, center_z=0.0, size_x=4.0, size_y=4.0, size_z=4.0,
    )
    runner = VinaRunner(tmp_path / "proj", runtime={"binary": "vina", "seed": 9, "replicates": 1})
    payload = runner.run([row], dry_run=True)
    assert payload["jobs"][0]["wall_time_s"] is None


# ----------------------------------------------------------------------------- blocker 9: md-inputs status


def test_md_inputs_cli_prints_the_status_and_the_reason_when_blocked(tmp_path, capsys):
    from workflow.cli import main

    project = _auto_project(tmp_path, "canonical")
    (_prepared_dirs_for(project)[0] / "prot.pdbqt.preparation.json").unlink()
    exit_code = main(
        ["analyze", "md-inputs", "--project-dir", str(project), "--engine", "vina", "--auto-maps"]
    )
    out = capsys.readouterr().out
    assert exit_code == 1
    assert "MD-input status: blocked" in out
    assert "MD-input reason: md_input_map_unavailable:receptor_preparation_record_missing" in out
