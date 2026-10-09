"""Offline scoring propagation contracts; no engines or schedulers are executed."""

import json
import shlex
import subprocess
from dataclasses import asdict
from pathlib import Path

import pandas as pd
import pytest

from docking.cli import deploy_main, dock_main
from docking.engine_registry import build_runner
from docking.models import PairlistRow
from docking.project_layout import (
    bootstrap_project_layout,
    load_manifest,
    pairlist_path,
    shared_ligands_dir,
    shared_receptors_dir,
)
from workflow.cli import main as workflow_main

pytestmark = pytest.mark.filterwarnings("ignore:Created compatibility symlink:RuntimeWarning")


@pytest.fixture(autouse=True)
def prohibit_processes(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("Offline command planning must not launch a process")

    monkeypatch.setattr(subprocess, "run", forbidden)
    monkeypatch.setattr(subprocess, "Popen", forbidden)


@pytest.fixture
def row():
    return PairlistRow("fixture_r.pdbqt", "fixture_site", "fixture_l.pdbqt", 1, 2, 3, 20, 22, 24)


@pytest.fixture
def project(tmp_path, row):
    root = tmp_path / "project"
    bootstrap_project_layout(root, engines=["gnina", "vina", "smina"], layout_profile="canonical")
    pd.DataFrame([asdict(row)]).to_csv(pairlist_path(root), index=False)
    (shared_receptors_dir(root) / row.receptor).write_text("REMARK synthetic fixture\n")
    (shared_ligands_dir(root) / row.ligand).write_text("REMARK synthetic fixture\n")
    return root


def command_for(engine, project, row, runtime):
    runner = build_runner(engine, project, runtime)
    return runner.build_command(row, Path("fixture.pose"), Path("fixture.log"))


def assert_option(command, flag, value):
    assert command.count(flag) == 1
    assert command[command.index(flag) + 1] == str(value)


@pytest.mark.parametrize("engine,scoring", [("gnina", "vina"), ("vina", "vinardo"), ("smina", "vinardo")])
def test_runner_forwards_explicit_scoring_once(engine, scoring, project, row):
    command = command_for(engine, project, row, {"scoring": scoring})
    assert_option(command, "--scoring", scoring)


@pytest.mark.parametrize("image", [None, "fixture.sif"])
def test_gnina_forwards_cpu_with_gpu_device_zero(image, project, row):
    command = command_for("gnina", project, row, {"image": image, "use_gpu": True, "device": 0, "cpu": 8, "seed": 0})
    assert_option(command, "--device", 0)
    assert_option(command, "--cpu", 8)
    assert_option(command, "--seed", 0)
    if image:
        assert command[:4] == ["apptainer", "exec", "--nv", "fixture.sif"]


@pytest.mark.parametrize("scoring", [None, "", "   "])
@pytest.mark.parametrize("engine", ["gnina", "vina", "smina"])
def test_omitted_or_empty_scoring_preserves_existing_command(engine, scoring, project, row):
    command = command_for(engine, project, row, {"scoring": scoring})
    receptor = str(shared_receptors_dir(project) / row.receptor)
    ligand = str(shared_ligands_dir(project) / row.ligand)
    expected = [engine, "-r" if engine == "gnina" else "--receptor", receptor,
                "-l" if engine == "gnina" else "--ligand", ligand,
                "--center_x", "1", "--center_y", "2", "--center_z", "3",
                "--size_x", "20", "--size_y", "22", "--size_z", "24",
                "--exhaustiveness", "16", "--num_modes", "20"]
    if engine == "gnina":
        expected += ["--cnn_scoring", "rescore", "--out", "fixture.pose", "--log", "fixture.log"]
    else:
        expected += ["--out", "fixture.pose"]
    assert command == expected


def test_gnina_cpu_only_runtime_still_omits_device_and_container_nv(project, row):
    command = command_for("gnina", project, row, {"image": "fixture.sif", "use_gpu": False, "device": 0, "cpu": 8})
    assert_option(command, "--cpu", 8)
    assert "--device" not in command
    assert "--nv" not in command


def invoke(entry, operation, argv):
    if entry == "workflow":
        return workflow_main(["dock", operation, *argv])
    return (deploy_main if operation == "deploy" else dock_main)(argv)


def common_argv(project):
    # QC gates are explicitly disabled only for these non-molecular fixtures.
    return ["--project-dir", str(project), "--engines", "gnina,vina,smina",
            "--no-ligand-qc-gate", "--no-receptor-qc-gate",
            "--gnina-device", "0", "--gnina-cpu", "8", "--vina-cpu", "4",
            "--smina-cpu", "4", "--exhaustiveness", "32", "--num-modes", "20"]


@pytest.mark.parametrize("entry", ["workflow", "docking"])
@pytest.mark.parametrize("explicit", [False, True])
@pytest.mark.parametrize("seed", [0, 101])
def test_local_cli_writes_scoring_and_resources_to_real_run_manifest(entry, explicit, seed, project):
    argv = common_argv(project) + ["--dry-run", "--parameter-mode", "advanced",
                                    "--execution-environment", "local_gpu", "--seed", str(seed)]
    if explicit:
        argv += ["--gnina-scoring", "vina", "--vina-scoring", "vinardo", "--smina-scoring", "vinardo"]
    assert invoke(entry, "run", argv) == 0
    for engine, scoring in [("gnina", "vina"), ("vina", "vinardo"), ("smina", "vinardo")]:
        runner = build_runner(engine, project, {})
        payload = json.loads((runner.layout["root"] / "run_manifest.json").read_text())
        assert payload["runtime"].get("scoring") == (scoring if explicit else None)
        assert len(payload["jobs"]) == 1
        job = payload["jobs"][0]
        assert job["status"] == "dry_run"
        command = job["command"]
        if explicit:
            assert_option(command, "--scoring", scoring)
        else:
            assert "--scoring" not in command
        assert_option(command, "--exhaustiveness", 32)
        assert_option(command, "--num_modes", 20)
        assert_option(command, "--seed", seed)
        assert_option(command, "--size_y", "22.0")
        assert_option(command, "--cpu", 8 if engine == "gnina" else 4)
        if engine == "gnina":
            assert_option(command, "--device", 0)
            assert_option(command, "--cnn_scoring", "rescore")


@pytest.mark.parametrize("entry", ["workflow", "docking"])
@pytest.mark.parametrize("scheduler", ["slurm", "condor"])
@pytest.mark.parametrize("source", ["omitted", "profile", "cli_override"])
def test_deploy_cli_freezes_scoring_with_profile_precedence(entry, scheduler, source, project, tmp_path):
    profile = {"name": "offline-fixture", "scheduler": scheduler,
               "slurm": {"partition": "fixture-partition"}, "engines": {}}
    for engine in ["gnina", "vina", "smina"]:
        runtime = {"seed": 0}
        if source != "omitted":
            runtime["scoring"] = "vina"
        profile["engines"][engine] = {"runtime": runtime}
    profile_file = tmp_path / "profile.json"
    profile_file.write_text(json.dumps(profile))
    argv = common_argv(project) + ["--hpc-profile-file", str(profile_file)]
    if source == "cli_override":
        argv += ["--gnina-scoring", "vinardo", "--vina-scoring", "vinardo", "--smina-scoring", "vinardo"]
    assert invoke(entry, "deploy", argv) == 0
    stage = Path(load_manifest(project)["deployment_root"])
    for engine in ["gnina", "vina", "smina"]:
        payload = json.loads((stage / engine / "manifest" / "deployment_manifest.json").read_text())
        expected_scoring = {"omitted": None, "profile": "vina", "cli_override": "vinardo"}[source]
        assert payload["runtime"].get("scoring") == expected_scoring
        assert len(payload["jobs"]) == 1
        command = payload["jobs"][0]["command"]
        if expected_scoring:
            assert_option(command, "--scoring", expected_scoring)
        else:
            assert "--scoring" not in command
        assert_option(command, "--seed", 0)
        assert_option(command, "--exhaustiveness", 32)
        assert_option(command, "--num_modes", 20)
        assert_option(command, "--cpu", 8 if engine == "gnina" else 4)
        if engine == "gnina":
            assert_option(command, "--device", 0)
        # The scheduler's executable manifest must carry the same frozen argv.
        command_line = Path(payload["job_manifest"]).read_text().strip().split("\t")[2]
        assert shlex.split(command_line) == command


@pytest.mark.parametrize("entry", ["workflow", "docking"])
@pytest.mark.parametrize("operation", ["run", "deploy"])
@pytest.mark.parametrize("option", [["--gnina-scoring"], ["--vina-scoring"], ["--gnina-scorng", "vina"]])
def test_cli_rejects_missing_scoring_values_and_unknown_options(entry, operation, option, project):
    with pytest.raises(SystemExit) as error:
        invoke(entry, operation, ["--project-dir", str(project), *option])
    assert error.value.code == 2


def test_explicit_empirical_scoring_preserves_named_gnina_metrics(project, row):
    runner = build_runner("gnina", project, {"scoring": "vina"})
    (runner.layout["poses"] / f"{row.tag}.sdf").write_text("$$$$\n")
    (runner.layout["logs"] / f"{row.tag}.log").write_text(
        "mode | affinity | intramol | CNN pose score | CNN affinity\n"
        "1 -7.00 0.00 0.90 4.00\n2 -8.00 0.00 0.20 8.00\n"
    )
    scores = runner.collect_normalized_scores([row])
    assert scores["score_name_primary"].tolist() == ["cnn_affinity", "cnn_affinity"]
    assert scores["score_primary"].tolist() == [4.0, 8.0]
    assert scores["score_name_secondary"].tolist() == ["cnn_score", "cnn_score"]
    assert scores["score_secondary"].tolist() == [0.9, 0.2]
    named_scores = pd.read_csv(runner.layout["scores"] / "all_scores.csv")
    assert named_scores["cnn_score"].tolist() == [0.9, 0.2]
    assert named_scores["cnn_affinity"].tolist() == [4.0, 8.0]
