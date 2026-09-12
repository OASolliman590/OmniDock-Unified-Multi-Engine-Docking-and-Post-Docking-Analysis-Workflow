"""Regression coverage for scientific execution contracts (no docking install)."""
from dataclasses import replace
import json
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

import pandas as pd
import pytest

from docking.cli import dock_main, deploy_main, _pairlist_rows_from_df
from docking.deployment import generate_slurm_deployment, generate_condor_deployment
from docking.models import PairlistRow
from docking.parameter_schema import CommonDockingParameters, DockingParameterSchema, transform_pairlist_rows, validate_parameter_schema
from docking.preflight import run_docking_preflight
from docking.project_layout import ensure_project_layout, save_manifest, pairlist_path, ensure_post_docking_compat_shim
from docking.runners.vina import VinaRunner
from docking.runners.gnina import GninaRunner
from docking.runners.autodock4 import AutoDock4Runner
from docking.runners.job_contract import execute_job, validate_pose_file, parse_sdf_scores


ATOM = "ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\n"
PDBQT = "MODEL 1\nREMARK VINA RESULT: -8.0 0.0 0.0\n" + ATOM + "ENDMDL\n"
SDF = "ligand\nfixture\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> <minimizedAffinity>\n-8.0\n\n> <CNNaffinity>\n7.0\n\n> <CNNscore>\n0.9\n\n$$$$\n"


@pytest.fixture
def project(tmp_path):
    root = tmp_path / "project with spaces"
    layout = ensure_project_layout(root)
    (layout["receptors"] / "R.pdbqt").write_text(ATOM)
    (layout["ligands"] / "L.pdbqt").write_text("ROOT\n" + ATOM + "ENDROOT\nTORSDOF 0\n")
    pair = PairlistRow("R.pdbqt", "site_1", "L.pdbqt", 0, 0, 0, 20, 20, 20)
    pd.DataFrame([pair.to_dict()]).to_csv(pairlist_path(root), index=False)
    save_manifest(root, {"engines": ["vina"], "pairlist_file": str(pairlist_path(root))})
    return root, pair


def output_process(pose_text=PDBQT):
    def run(command, **kwargs):
        Path(command[command.index("--out") + 1]).write_text(pose_text, encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "engine output", "")
    return run


def runner(project, **runtime):
    return VinaRunner(project[0], {"binary": sys.executable, "seed": 42, **runtime})


def test_resume_requires_exact_inputs_protocol_and_validated_output(project):
    root, pair = project
    engine = runner(project)
    with patch("docking.runners.job_contract.subprocess.run", side_effect=output_process()):
        result = engine.run([pair])
    assert result["jobs"][0]["status"] == "completed"
    assert engine.plan_jobs([pair], skip_completed=True)[0].status == "skipped"
    assert runner(project, seed=43).plan_jobs([pair], skip_completed=True)[0].status == "planned"
    assert engine.plan_jobs([replace(pair, size_x=24)], skip_completed=True)[0].status == "planned"
    (root / "ligands" / "L.pdbqt").write_text(ATOM + "REMARK changed\n")
    assert engine.plan_jobs([pair], skip_completed=True)[0].status == "planned"


def test_empty_legacy_pose_does_not_resume(project):
    engine = runner(project)
    job = engine.plan_jobs([project[1]])[0]
    Path(job.pose_file).touch()
    assert engine.plan_jobs([project[1]], skip_completed=True)[0].status == "planned"


def test_failed_rerun_retires_old_pose_and_clears_scores(project):
    engine = runner(project)
    with patch("docking.runners.job_contract.subprocess.run", side_effect=output_process()):
        engine.run([project[1]])
    with patch("docking.runners.job_contract.subprocess.run", return_value=subprocess.CompletedProcess([], 1, "", "failure")):
        result = engine.run([project[1]])
    assert result["jobs"][0]["status"] == "failed"
    assert not Path(result["jobs"][0]["pose_file"]).exists()
    assert pd.read_csv(engine.layout["scores"] / "normalized_scores.csv").empty
    assert list((engine.layout["poses"] / ".attempts").rglob("*.pdbqt"))
    completion = json.loads(Path(result["jobs"][0]["completion_file"]).read_text())
    assert completion["status"] == "failed"


@pytest.mark.parametrize("mode", ["missing", "empty", "truncated", "startup"])
def test_invalid_output_and_startup_failure_never_complete(project, mode):
    engine = runner(project)
    def process(command, **kwargs):
        if mode == "startup":
            raise FileNotFoundError("engine vanished")
        if mode != "missing":
            Path(command[command.index("--out") + 1]).write_text("" if mode == "empty" else PDBQT.replace("ENDMDL\n", ""))
        return subprocess.CompletedProcess(command, 0, "", "")
    with patch("docking.runners.job_contract.subprocess.run", side_effect=process):
        result = engine.run([project[1]])
    assert result["jobs"][0]["status"] == "failed"
    assert pd.read_csv(engine.layout["scores"] / "normalized_scores.csv").empty


def test_changed_input_after_plan_does_not_execute(project):
    engine = runner(project)
    job = engine.plan_jobs([project[1]])[0].to_dict()
    Path(job["input_files"][0]["path"]).write_text("changed")
    with patch("docking.runners.job_contract.subprocess.run") as process:
        result = execute_job(job)
    assert result["status"] == "failed"
    process.assert_not_called()


def test_input_mutation_during_execution_invalidates_result(project):
    engine = runner(project)
    def process(command, **kwargs):
        result = output_process()(command, **kwargs)
        (project[0] / "ligands" / "L.pdbqt").write_text("mutated during run")
        return result
    with patch("docking.runners.job_contract.subprocess.run", side_effect=process):
        result = engine.run([project[1]])
    assert result["jobs"][0]["status"] == "failed"
    assert "changed during execution" in result["jobs"][0]["error"]


def test_existing_job_lock_prevents_artifact_retirement(project):
    engine = runner(project)
    job = engine.plan_jobs([project[1]])[0].to_dict()
    pose = Path(job["pose_file"])
    pose.write_text(PDBQT)
    Path(job["completion_file"] + ".lock").write_text("active attempt")
    with patch("docking.runners.job_contract.subprocess.run") as process:
        result = execute_job(job)
    assert result["status"] == "failed" and "job lock" in result["error"]
    assert pose.read_text() == PDBQT
    process.assert_not_called()


def test_interruption_clears_old_score_table(project):
    engine = runner(project)
    with patch("docking.runners.job_contract.subprocess.run", side_effect=output_process()):
        engine.run([project[1]])
    with patch("docking.runners.job_contract.subprocess.run", side_effect=KeyboardInterrupt):
        with pytest.raises(KeyboardInterrupt):
            engine.run([project[1]])
    assert pd.read_csv(engine.layout["scores"] / "normalized_scores.csv").empty
    assert not list(engine.layout["poses"].glob("*.pdbqt"))


def test_cli_returns_failure_for_failed_engine(project):
    with patch("docking.cli.run_docking_preflight", return_value={"ok": True}), patch("docking.runners.job_contract.subprocess.run", return_value=subprocess.CompletedProcess([], 1, "", "failure")):
        assert dock_main(["--project-dir", str(project[0]), "--engines", "vina"]) == 1


@pytest.mark.parametrize("factory,options", [(generate_slurm_deployment, {"slurm_options": {"partition": "test"}}), (generate_condor_deployment, {})])
def test_portable_deployment_uses_shared_executor_and_posix_paths(project, factory, options):
    root, pair = project
    result = factory(project_root=root, engines=["vina"], pairlist_rows=[pair], runtime_by_engine={"vina": {"seed": 42}}, round_id="r1", docking_mode="screen", remote_project_root="/cluster/project", **options)
    assert result["execution_project_root"] == "/cluster/project"
    payload = json.loads(next((root / "deployments" / "r1" / "screen" / "vina" / "manifest").glob("deployment_manifest.json")).read_text())
    job = payload["jobs"][0]
    assert job["completion_file"].startswith("/cluster/project/")
    assert all("\\" not in token for token in job["command"])
    assert all(item["path"].startswith("/cluster/project/") for item in job["input_files"])
    tsv = next((root / "deployments" / "r1" / "screen" / "vina" / "manifest").glob("*.tsv")).read_text()
    assert "execute_job.py" in tsv


def test_deploy_preserves_reviewed_local_protocol(project):
    root, _ = project
    assert dock_main(["--project-dir", str(root), "--engines", "vina", "--dry-run", "--seed", "42", "--box-scale", "2", "--box-padding", "1", "--exhaustiveness", "28", "--no-ligand-qc-gate", "--no-receptor-qc-gate"]) == 0
    assert deploy_main(["--project-dir", str(root), "--engines", "vina", "--slurm-partition", "test", "--no-ligand-qc-gate", "--no-receptor-qc-gate"]) == 0
    payload = json.loads(next((root / "deployments").rglob("vina/manifest/deployment_manifest.json")).read_text())
    command = payload["jobs"][0]["command"]
    assert command[command.index("--seed") + 1] == "42"
    assert command[command.index("--size_x") + 1] == "42.0"
    assert command[command.index("--exhaustiveness") + 1] == "28"


def test_gnina_resource_flags_independent(project):
    root, pair = project
    cpu = GninaRunner(root, {"use_gpu": False, "cpu": 4}).build_command(pair, root / "a.sdf", root / "a.log")
    gpu = GninaRunner(root, {"use_gpu": True, "device": 0, "cpu": 4}).build_command(pair, root / "a.sdf", root / "a.log")
    assert "--no_gpu" in cpu
    assert "--cpu" in gpu and "--device" in gpu and "--no_gpu" not in gpu


def test_gnina_scores_come_from_poses_not_runner_log(project):
    root, pair = project
    engine = GninaRunner(root, {"binary": sys.executable, "use_gpu": False})
    (engine.layout["logs"] / (pair.tag + ".runner.log")).write_text("mode | affinity\n1 -99.00 0.00 0.99 12.00\n")
    with patch("docking.runners.job_contract.subprocess.run", side_effect=output_process(SDF)):
        result = engine.run([pair])
    assert result["jobs"][0]["status"] == "completed"
    scores = pd.read_csv(engine.layout["scores"] / "normalized_scores.csv")
    assert scores.tag.tolist() == [pair.tag]
    assert scores.affinity_kcal_mol.tolist() == [-8.0]


def test_ad4_grid_and_parameter_files_are_portable_and_consistent(project):
    root, pair = project
    custom = root / "custom parameters.dat"
    custom.write_text("# test parameter data\n")
    engine = AutoDock4Runner(root, {"parameter_file": str(custom)})
    job = engine.plan_jobs([pair])[0]
    directory = engine.layout["root"] / "pairs" / pair.tag
    gpf, dpf = [(directory / name).read_text() for name in ("grid.gpf", "docking.dpf")]
    assert "npts 54 54 54" in gpf
    assert "parameter_file parameters.dat" in gpf and "parameter_file parameters.dat" in dpf
    assert str(root) not in gpf + dpf
    assert "receptor receptor.pdbqt" in gpf and "move ligand.pdbqt" in dpf
    assert "cd " in job.command[-1]
    assert (directory / "parameters.dat").read_bytes() == custom.read_bytes()


def test_ad4_adt_branch_supplies_same_custom_parameters(project):
    root, pair = project
    custom = root / "parameters.dat"
    custom.write_text("# parameters")
    engine = AutoDock4Runner(root, {"parameter_file": str(custom), "autodocktools_python": "python2", "autodocktools_prepare_gpf4": "/tools/prepare_gpf4.py", "autodocktools_prepare_dpf4": "/tools/prepare_dpf4.py"})
    command = engine.build_command(pair, root / "output.dlg", root / "output.log")[-1]
    assert command.count("parameter_file=parameters.dat") == 2
    assert "npts=54,54,54" in command


@pytest.mark.parametrize("kind", ["ligand", "receptor"])
def test_required_qc_exceptions_block(project, kind):
    root, pair = project
    with patch(f"docking.preflight.audit_project_{kind}s", side_effect=RuntimeError("QC unavailable")):
        result = run_docking_preflight(project_dir=root, engines=["vina"], pairlist_rows=[pair], runtime_by_engine={}, dry_run=True, **{f"enable_{kind}_qc": True})
    assert not result["ok"]


@pytest.mark.parametrize("attribute,value", [("center_x", float("nan")), ("size_x", -1), ("size_y", float("inf")), ("size_z", 0)])
def test_invalid_geometry_blocks(project, attribute, value):
    root, pair = project
    bad = replace(pair, **{attribute: value})
    result = run_docking_preflight(project_dir=root, engines=["vina"], pairlist_rows=[bad], runtime_by_engine={}, dry_run=True)
    assert not result["ok"]
    with pytest.raises(ValueError):
        runner(project).plan_jobs([bad])


def test_nonfinite_box_schema_rejected():
    errors, _ = validate_parameter_schema(DockingParameterSchema(common=CommonDockingParameters(box_scale=float("nan"), box_padding=float("inf"))), ["vina"])
    assert len(errors) == 2


def test_reference_coordinate_provenance_survives_conversion_and_scaling(project):
    pair = replace(project[1], reference_source="crystal", reference_frame_id="source-frame", receptor_frame_id="source-frame", reference_pdb_id="1ABC", reference_pose_file="reference.sdf")
    parsed = _pairlist_rows_from_df(pd.DataFrame([pair.to_dict()]))[0]
    scaled = transform_pairlist_rows([parsed], DockingParameterSchema(common=CommonDockingParameters(box_scale=2)))[0]
    assert scaled.reference_frame_id == "source-frame" and scaled.receptor_frame_id == "source-frame"
    assert scaled.reference_pose_file == "reference.sdf" and scaled.size_x == 40


def test_symlink_privilege_failure_has_usable_fallback(tmp_path):
    with patch.object(Path, "symlink_to", side_effect=OSError("privilege unavailable")):
        result = ensure_post_docking_compat_shim(tmp_path)
    assert result["status"] == "directory_fallback"
    assert result["canonical_root"].is_dir()


def test_multimodel_sdf_and_ad4_completion_validation(tmp_path):
    sdf = tmp_path / "two.sdf"
    sdf.write_text(SDF + SDF)
    assert len(parse_sdf_scores(sdf)) == 2
    dlg = tmp_path / "docking.dlg"
    docked = "MODEL 1\nUSER Estimated Free Energy of Binding = -8.0 kcal/mol\n" + ATOM + "ENDMDL\n"
    text = "\n".join("DOCKED: " + line for line in docked.splitlines()) + "\n"
    dlg.write_text(text)
    with pytest.raises(ValueError, match="Successful Completion"):
        validate_pose_file(dlg, "autodock4")
    dlg.write_text(text + "autodock4: Successful Completion.\n")
    assert validate_pose_file(dlg, "autodock4") == 1


def test_bundled_executor_runs_without_repository_imports(project):
    root, pair = project
    generate_slurm_deployment(project_root=root, engines=["vina"], pairlist_rows=[pair], runtime_by_engine={"vina": {"binary": sys.executable}}, round_id="standalone", docking_mode="screen", slurm_options={"partition": "test"})
    manifest = next((root / "deployments" / "standalone").rglob("vina/manifest/deployment_manifest.json"))
    payload = json.loads(manifest.read_text())
    job = payload["jobs"][0]
    job["command"] = [sys.executable, "-c", f"from pathlib import Path; Path({job['pose_file']!r}).write_text({PDBQT!r})"]
    manifest.write_text(json.dumps(payload))
    helper = manifest.parent.parent / "job_scripts" / "execute_job.py"
    completed = subprocess.run([sys.executable, str(helper), str(manifest), "0"], cwd=root, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert json.loads(Path(job["completion_file"]).read_text())["status"] == "completed"


def test_wrapped_or_partially_resolved_engines_cannot_certify_resume(project):
    from docking.runners.job_contract import execution_identity, completion_valid
    root, pair = project
    wrapped = VinaRunner(root, {"conda_env": "docking", "binary": "vina"})
    job = wrapped.plan_jobs([pair])[0].to_dict()
    assert job["executables"] == []
    with patch("docking.runners.job_contract.subprocess.run", side_effect=output_process()):
        assert execute_job(job)["status"] == "completed"
    assert not completion_valid(job)
    assert execution_identity({"executables": [sys.executable, "nonexistent-omnidock-engine-12345"]}) == {}


def test_attempt_provenance_does_not_repeat_long_pose_names(tmp_path):
    from docking.runners.job_contract import file_hash
    # Leave space for the actual pose/completion filenames but make the former
    # `<pose name>.<UUID>/result.json.<UUID>.tmp` layout exceed 260 characters.
    output_dir = tmp_path / ("x" * max(1, 140 - len(str(tmp_path))))
    output_dir.mkdir()
    pose = output_dir / ("long_protein_ligand_site_name" * 2 + ".pdbqt")
    job = {"engine": "vina", "command": ["fake", "--out", str(pose)],
           "pose_file": str(pose), "log_file": str(output_dir / "job.log"),
           "completion_file": str(output_dir / "completion.json"),
           "fingerprint": "fixture", "input_files": [], "executables": [sys.executable]}
    with patch("docking.runners.job_contract.subprocess.run", side_effect=output_process()):
        assert execute_job(job)["status"] == "completed"
        assert execute_job(job)["status"] == "completed"
    attempts = list((output_dir / ".attempts").iterdir())
    assert len(attempts) == 2 and all(len(path.name) == 32 for path in attempts)
    archived_pose = next((output_dir / ".attempts").rglob(pose.name))
    assert file_hash(archived_pose) == file_hash(pose)
    assert all((attempt / "result.json").is_file() for attempt in attempts)
    assert not list(output_dir.rglob("*.tmp"))
