"""Spec 037 R4c: ``dock deploy`` applies the dock run parameter schema.

Rules tested (each on a real canonical project, QC gates disabled, no submission):
- the seed is required;
- replicates must be at least 1;
- energy_range is passed explicitly, must be positive, and is recorded with its value;
- basic mode refuses an explicit --num-modes and an explicit --exhaustiveness;
- advanced mode accepts them and records them as user_explicit;
- the effective values and their sources are recorded in the deployment bundle manifest.
"""

from __future__ import annotations

import json
from pathlib import Path

import argparse

import pandas as pd
import pytest

import docking.cli as docking_cli
from docking.project_layout import bootstrap_project_layout, pairlist_path


PAIR_COLUMNS = [
    "receptor", "site_id", "ligand", "center_x", "center_y", "center_z", "size_x", "size_y", "size_z",
]


def _project(root: Path) -> Path:
    bootstrap_project_layout(root, engines=["vina"], layout_profile="canonical")
    pd.DataFrame([{
        "receptor": "1IEP.pdbqt", "site_id": "site_1", "ligand": "STI.pdbqt",
        "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
        "size_x": 20.0, "size_y": 20.0, "size_z": 20.0,
    }], columns=PAIR_COLUMNS).to_csv(pairlist_path(root), index=False)
    return root


def _deploy(project: Path, *extra: str) -> int:
    return docking_cli.deploy_main([
        "--project-dir", str(project),
        "--engines", "vina",
        "--slurm-partition", "test-partition",
        "--no-ligand-qc-gate",
        "--no-receptor-qc-gate",
        *extra,
    ])


def _bundle_manifest(project: Path) -> dict:
    path = next(project.rglob("deployment_manifest.json"))
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture()
def project(tmp_path):
    return _project(tmp_path / "project")


def test_deploy_requires_a_seed_and_writes_no_bundle(project, capsys):
    status = _deploy(project)
    assert status == 1
    assert "seed" in capsys.readouterr().out
    assert not list(project.rglob("deployment_manifest.json"))


def test_deploy_refuses_zero_replicates(project, capsys):
    status = _deploy(project, "--seed", "7", "--replicates", "0")
    assert status == 1
    assert "`replicates` must be >= 1" in capsys.readouterr().out


def test_deploy_passes_energy_range_explicitly_and_records_it(project):
    assert _deploy(project, "--seed", "7") == 0
    manifest = _bundle_manifest(project)
    assert manifest["engine_effective_parameters"]["vina"]["energy_range"] == pytest.approx(3.0)
    runtime = json.loads(next(project.rglob("manifest/deployment_manifest.json")).read_text(encoding="utf-8"))
    assert runtime["runtime"]["energy_range"] == pytest.approx(3.0)


def test_deploy_records_an_explicit_energy_range(tmp_path):
    project = _project(tmp_path / "project_energy")
    assert _deploy(project, "--seed", "7", "--energy-range", "4.5") == 0
    assert _bundle_manifest(project)["engine_effective_parameters"]["vina"]["energy_range"] == pytest.approx(4.5)


def test_deploy_refuses_a_non_positive_energy_range(project, capsys):
    status = _deploy(project, "--seed", "7", "--energy-range", "0")
    assert status == 1
    assert "`energy_range` must be > 0" in capsys.readouterr().out


def test_basic_mode_refuses_an_explicit_num_modes(project, capsys):
    status = _deploy(project, "--seed", "7", "--num-modes", "10")
    assert status == 1
    assert "`--num-modes` is not used in basic mode" in capsys.readouterr().out
    assert not list(project.rglob("deployment_manifest.json"))


def test_basic_mode_refuses_an_explicit_exhaustiveness(project, capsys):
    status = _deploy(project, "--seed", "7", "--exhaustiveness", "8")
    assert status == 1
    assert "`--exhaustiveness` is not used in basic mode" in capsys.readouterr().out


def test_basic_default_records_the_preset_values_and_sources(project):
    assert _deploy(project, "--seed", "7") == 0
    effective = _bundle_manifest(project)["engine_effective_parameters"]["vina"]
    assert effective["exhaustiveness"] == 32
    assert effective["exhaustiveness_source"] == "basic_preset:exhaustive"
    assert effective["num_modes"] == 40
    assert effective["num_modes_source"] == "basic_preset:exhaustive"
    assert effective["parameter_mode"] == "basic"
    assert effective["parameter_preset"] == "exhaustive"


def test_advanced_mode_accepts_and_records_explicit_values(project):
    status = _deploy(project, "--seed", "7", "--parameter-mode", "advanced", "--exhaustiveness", "8", "--num-modes", "12")
    assert status == 0
    effective = _bundle_manifest(project)["engine_effective_parameters"]["vina"]
    assert (effective["exhaustiveness"], effective["exhaustiveness_source"]) == (8, "user_explicit")
    assert (effective["num_modes"], effective["num_modes_source"]) == (12, "user_explicit")


def test_bundle_records_the_full_parameter_schema(project):
    assert _deploy(project, "--seed", "99", "--replicates", "2") == 0
    schema = _bundle_manifest(project)["parameter_schema"]
    assert schema["common"]["seed"] == 99
    assert schema["common"]["replicates"] == 2
    assert schema["mode"] == "basic"


def test_dock_run_refuses_the_same_explicit_num_modes_as_deploy(project, capsys):
    status = docking_cli.dock_main([
        "--project-dir", str(project),
        "--engines", "vina",
        "--seed", "7",
        "--num-modes", "10",
        "--dry-run",
        "--no-ligand-qc-gate",
        "--no-receptor-qc-gate",
    ])
    assert status == 1
    assert "`--num-modes` is not used in basic mode" in capsys.readouterr().out


def test_workflow_deploy_forwards_seed_and_parameter_flags():
    from workflow.cli import _add_deploy_args, _build_deploy_argv

    parser = argparse.ArgumentParser()
    _add_deploy_args(parser)
    args = parser.parse_args([
        "--project-dir", "project", "--seed", "5", "--replicates", "2", "--energy-range", "4",
        "--parameter-mode", "advanced", "--parameter-preset", "balanced",
    ])
    argv = _build_deploy_argv(args)
    assert argv[argv.index("--seed") + 1] == "5"
    assert argv[argv.index("--replicates") + 1] == "2"
    assert float(argv[argv.index("--energy-range") + 1]) == pytest.approx(4.0)
    assert argv[argv.index("--parameter-mode") + 1] == "advanced"
    assert argv[argv.index("--parameter-preset") + 1] == "balanced"
    # Not given, so not forwarded: basic mode then takes its values from the preset.
    assert "--exhaustiveness" not in argv and "--num-modes" not in argv
