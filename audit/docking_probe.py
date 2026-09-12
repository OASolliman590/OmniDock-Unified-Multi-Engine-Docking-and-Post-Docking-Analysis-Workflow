"""Read-only audit of source behavior; generated fixtures live under audit/probe-data."""
from __future__ import annotations

import contextlib
import io
import json
import subprocess
import sys
import warnings
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pandas as pd
from docking.models import PairlistRow
from docking.project_layout import ensure_project_layout, save_manifest, pairlist_path
from docking.runners.vina import VinaRunner
from docking.runners.gnina import GninaRunner
from docking.runners.autodock4 import AutoDock4Runner
from docking.deployment import generate_slurm_deployment
from docking.cli import dock_main, deploy_main
from docking.preflight import run_docking_preflight
from docking.parameter_schema import CommonDockingParameters, DockingParameterSchema, validate_parameter_schema
from post_docking_analysis.docking_parser import parse_vina_pdbqt


def main():
    root = ROOT / "audit" / "probe-data" / "docking"
    root.mkdir(parents=True, exist_ok=True)
    # Pre-existing compatibility directory avoids Windows symlink privilege needs.
    (root / "5-Post_Docking_Analysis").mkdir(exist_ok=True)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    layout = ensure_project_layout(root)
    atom = "ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\n"
    (layout["receptors"] / "R.pdbqt").write_text(atom, encoding="utf-8")
    (layout["ligands"] / "L.pdbqt").write_text("ROOT\n" + atom + "ENDROOT\nTORSDOF 0\n", encoding="utf-8")
    pair = PairlistRow("R.pdbqt", "site_1", "L.pdbqt", 0, 0, 0, 20, 20, 20)
    pd.DataFrame([pair.to_dict()]).to_csv(pairlist_path(root), index=False)
    save_manifest(root, {"engines": ["vina"], "pairlist_file": str(pairlist_path(root))})
    results = {}

    vina = VinaRunner(root, {"seed": 42})
    planned = vina.plan_jobs([pair])[0]
    Path(planned.pose_file).write_text("", encoding="utf-8")
    results["empty_pose_status"] = vina.plan_jobs([pair], skip_completed=True)[0].status
    Path(planned.pose_file).write_text("MODEL 1\nREMARK VINA RESULT: -9.0 0.0 0.0\nENDMDL\n", encoding="utf-8")
    Path(planned.log_file).write_text("OLD LOG", encoding="utf-8")
    with patch("docking.runners.base.subprocess.run", return_value=subprocess.CompletedProcess([], 1, "", "engine failed")):
        failed = vina.run([pair])
    results["failed_run"] = {
        "status": failed["jobs"][0]["status"],
        "retained_pair_log": Path(planned.log_file).read_text(encoding="utf-8"),
        "normalized_scores": pd.read_csv(vina.layout["scores"] / "normalized_scores.csv").to_dict("records"),
    }
    with patch("docking.cli.run_docking_preflight", return_value={"ok": True}), patch(
        "docking.runners.base.subprocess.run", return_value=subprocess.CompletedProcess([], 1, "", "engine failed")
    ), contextlib.redirect_stdout(io.StringIO()):
        results["failed_cli_exit_code"] = dock_main(["--project-dir", str(root), "--engines", "vina"])

    cpu = GninaRunner(root, {"use_gpu": False, "cpu": 4}).build_command(pair, root / "cpu.sdf", root / "cpu.log")
    gpu = GninaRunner(root, {"use_gpu": True, "device": 0, "cpu": 4}).build_command(pair, root / "gpu.sdf", root / "gpu.log")
    results["gnina_cpu_command"] = cpu
    results["gnina_gpu_command"] = gpu
    gnina = GninaRunner(root, {"use_gpu": False, "cpu": 4})
    table = "mode | affinity | intramol | CNN | CNN\n1 -8.00 0.00 0.90 7.00\n"
    with patch("docking.runners.base.subprocess.run", return_value=subprocess.CompletedProcess([], 0, table, "")), contextlib.redirect_stdout(io.StringIO()):
        gnina_result = gnina.run([pair])
    results["gnina_success_without_pose"] = {
        "status": gnina_result["jobs"][0]["status"],
        "pose_exists": Path(gnina_result["jobs"][0]["pose_file"]).exists(),
        "normalized_rows": pd.read_csv(gnina.layout["scores"] / "normalized_scores.csv").to_dict("records"),
    }
    with patch("docking.preflight.audit_project_receptors", side_effect=RuntimeError("synthetic QC failure")):
        results["qc_exception"] = run_docking_preflight(
            project_dir=root, engines=["vina"], pairlist_rows=[pair], runtime_by_engine={},
            enable_receptor_qc=True, dry_run=True,
        )
    bad_pair = PairlistRow("R.pdbqt", "site_1", "L.pdbqt", float("nan"), 0, 0, -20, float("inf"), 0)
    results["invalid_box_preflight"] = run_docking_preflight(
        project_dir=root, engines=["vina"], pairlist_rows=[bad_pair], runtime_by_engine={}, dry_run=True,
    )
    results["nan_schema_validation"] = validate_parameter_schema(
        DockingParameterSchema(common=CommonDockingParameters(box_scale=float("nan"), box_padding=float("inf"))), ["vina"]
    )
    smina_fixture = root / "smina_missing_rmsd.pdbqt"
    smina_fixture.write_text("MODEL 1\nREMARK minimizedAffinity -8.0\nENDMDL\n", encoding="utf-8")
    results["smina_missing_rmsd"] = parse_vina_pdbqt(smina_fixture).to_dict("records")

    ad4 = AutoDock4Runner(root, {"parameter_file": "custom-atom-parameters.dat"})
    ad_job = ad4.plan_jobs([pair])[0]
    ad_pair = ad4.layout["root"] / "pairs" / pair.tag
    results["ad4"] = {
        "npts_for_20A": ad4._grid_points(20, 0.375),
        "gpf": (ad_pair / "grid.gpf").read_text(encoding="utf-8"),
        "dpf": (ad_pair / "docking.dpf").read_text(encoding="utf-8"),
    }
    deployment = generate_slurm_deployment(
        project_root=root, pairlist_rows=[pair], engines=["autodock4"],
        runtime_by_engine={"autodock4": {"parameter_file": "custom-atom-parameters.dat"}},
        slurm_options={"partition": "audit"},
        round_id="audit_ad4", docking_mode="screen", remote_project_root=str(root.parent / "remote-docking"),
    )
    results["ad4_deployment"] = deployment
    results["ad4_config_still_has_local_paths"] = str(root) in (ad_pair / "grid.gpf").read_text(encoding="utf-8")

    with contextlib.redirect_stdout(io.StringIO()):
        results["dry_run_exit_code"] = dock_main([
            "--project-dir", str(root), "--engines", "vina", "--dry-run",
            "--parameter-mode", "advanced", "--seed", "42", "--box-scale", "2",
            "--no-ligand-qc-gate", "--no-receptor-qc-gate",
        ])
        local_manifest = json.loads((vina.layout["root"] / "run_manifest.json").read_text(encoding="utf-8"))
        results["local_parameter_command"] = local_manifest["jobs"][0]["command"]
        results["deploy_exit_code"] = deploy_main([
            "--project-dir", str(root), "--engines", "vina", "--round", "audit_parity",
            "--no-ligand-qc-gate", "--no-receptor-qc-gate",
            "--slurm-partition", "audit",
        ])
    results["parity_deployment_files"] = {
        str(p.relative_to(root)): p.read_text(encoding="utf-8")
        for p in root.rglob("job_array.tsv") if "audit_parity" in str(p)
    }
    out = ROOT / "audit" / "docking-probe-results.json"
    out.write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
    print(out)
    print(json.dumps({key: results[key] for key in ["empty_pose_status", "failed_cli_exit_code", "ad4_config_still_has_local_paths", "dry_run_exit_code", "deploy_exit_code"]}, indent=2))


if __name__ == "__main__":
    main()
