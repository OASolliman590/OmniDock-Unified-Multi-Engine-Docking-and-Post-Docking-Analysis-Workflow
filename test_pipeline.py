#!/usr/bin/env python3
"""End-to-end smoke test for the unified docking pipeline CLI."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional

ENGINES: List[str] = ["gnina", "vina", "smina", "autodock4"]
REPO_ROOT = Path(__file__).resolve().parent

RECEPTOR_PDBQT = """\
ATOM      1  N   ALA A   1      -0.500   0.000   0.000  1.00  0.00           N
ATOM      2  C   ALA A   1       0.500   0.000   0.000  1.00  0.00           C
END
"""

LIGAND_PDBQT = """\
ROOT
ATOM      1  C   LIG A   1       0.200   0.100   0.000  1.00  0.00           C
ENDROOT
TORSDOF 0
"""


def _run(command: Iterable[str], cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    result = subprocess.run(
        list(command),
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "Command failed\n"
            f"cmd: {' '.join(command)}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
    return result


def _write_pairlist(path: Path) -> None:
    rows = [
        {
            "receptor": "R1.pdbqt",
            "site_id": "site_1",
            "ligand": "L1.pdbqt",
            "center_x": 0.0,
            "center_y": 0.0,
            "center_z": 0.0,
            "size_x": 20.0,
            "size_y": 20.0,
            "size_z": 20.0,
        }
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "receptor",
                "site_id",
                "ligand",
                "center_x",
                "center_y",
                "center_z",
                "size_x",
                "size_y",
                "size_z",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def _assert_file(path: Path, label: str) -> None:
    if not path.exists():
        raise RuntimeError(f"Missing {label}: {path}")


def _assert_manifest_job_count(project_dir: Path, engine: str, expected: int = 1) -> None:
    candidates = [
        project_dir / "4-Docking" / f"{engine}_out" / "run_manifest.json",
        project_dir / "4-Docking" / "run_manifest.json",
    ]
    manifest_path = next((candidate for candidate in candidates if candidate.exists()), None)
    if manifest_path is None:
        raise RuntimeError(f"Missing {engine} run manifest in: {', '.join(str(path) for path in candidates)}")
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    if str(payload.get("engine", "")).lower() != engine:
        raise RuntimeError(
            f"Expected {engine} run manifest, got engine={payload.get('engine')} at {manifest_path}"
        )
    jobs = payload.get("jobs", [])
    if len(jobs) != expected:
        raise RuntimeError(f"{engine} expected {expected} job(s), found {len(jobs)}")
    statuses = sorted({str(job.get("status", "")).lower() for job in jobs})
    if statuses != ["dry_run"]:
        raise RuntimeError(f"{engine} expected dry_run status, got: {statuses}")


def _assert_autodock4_pair_assets(project_dir: Path) -> None:
    pair_dir = project_dir / "4-Docking" / "autodock4_out" / "pairs" / "R1.pdbqt_site_1_L1.pdbqt"
    _assert_file(pair_dir / "grid.gpf", "AutoDock4 GPF")
    _assert_file(pair_dir / "docking.dpf", "AutoDock4 DPF")


def run_smoke_test(keep_temp: bool = False) -> Dict[str, str]:
    temp_root = Path(tempfile.mkdtemp(prefix="pdbwiz_smoke_all_engines_"))
    project_dir = temp_root / "project"
    prepared_proteins = temp_root / "prepared_proteins"
    prepared_ligands = temp_root / "prepared_ligands"
    pairlist = temp_root / "pairlist.csv"

    prepared_proteins.mkdir(parents=True, exist_ok=True)
    prepared_ligands.mkdir(parents=True, exist_ok=True)

    (prepared_proteins / "R1.pdbqt").write_text(RECEPTOR_PDBQT, encoding="utf-8")
    (prepared_ligands / "L1.pdbqt").write_text(LIGAND_PDBQT, encoding="utf-8")
    _write_pairlist(pairlist)

    print("1) Initializing project layout...")
    _run(
        [
            "python",
            "main.py",
            "workflow",
            "init",
            "--project-dir",
            str(project_dir),
            "--layout-profile",
            "docking_legacy",
            "--engines",
            ",".join(ENGINES),
        ],
        cwd=REPO_ROOT,
    )

    print("2) Materializing docking project from pairlist...")
    _run(
        [
            "python",
            "main.py",
            "prep",
            "project",
            "--project-dir",
            str(project_dir),
            "--prepared-proteins",
            str(prepared_proteins),
            "--prepared-ligands",
            str(prepared_ligands),
            "--pairlist-file",
            str(pairlist),
            "--layout-profile",
            "docking_legacy",
            "--engines",
            ",".join(ENGINES),
        ],
        cwd=REPO_ROOT,
    )

    print("3) Running all engines in dry-run mode...")
    _run(
        [
            "python",
            "main.py",
            "dock",
            "run",
            "--project-dir",
            str(project_dir),
            "--engines",
            ",".join(ENGINES),
            "--dry-run",
            "--no-ligand-qc-gate",
            "--no-receptor-qc-gate",
        ],
        cwd=REPO_ROOT,
    )

    print("4) Verifying run manifests and generated AD4 parameter files...")
    for engine in ENGINES:
        _assert_manifest_job_count(project_dir, engine, expected=1)
    _assert_autodock4_pair_assets(project_dir)

    artifacts = {
        "temp_root": str(temp_root),
        "project_dir": str(project_dir),
        "pairlist": str(project_dir / "4-Docking" / "pairlist.csv"),
    }
    if not keep_temp:
        shutil.rmtree(temp_root, ignore_errors=True)
        artifacts["temp_root"] = "(deleted)"
        artifacts["project_dir"] = "(deleted)"
    return artifacts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a current-architecture smoke test for all docking engines.")
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep the temporary project on disk for inspection.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    print("🧪 Running unified all-engines smoke test")
    print("=" * 60)
    try:
        artifacts = run_smoke_test(keep_temp=args.keep_temp)
    except Exception as exc:
        print(f"❌ Smoke test failed: {exc}")
        return 1

    print("✅ Smoke test passed")
    print(f"   temp_root: {artifacts['temp_root']}")
    print(f"   project_dir: {artifacts['project_dir']}")
    print(f"   pairlist: {artifacts['pairlist']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
