from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

import docking.cli as docking_cli
import docking.deployment as deployment
from docking.models import PairlistRow
from workflow.cli import _add_deploy_args, _build_deploy_argv


PAIR_COLUMNS = [
    "receptor",
    "site_id",
    "ligand",
    "center_x",
    "center_y",
    "center_z",
    "size_x",
    "size_y",
    "size_z",
]


def _pair_row(index: int) -> dict[str, object]:
    return {
        "receptor": f"R{index}.pdbqt",
        "site_id": "site_1",
        "ligand": f"L{index}.pdbqt",
        "center_x": float(index),
        "center_y": float(index + 1),
        "center_z": float(index + 2),
        "size_x": 20.0,
        "size_y": 20.0,
        "size_z": 20.0,
    }


def _write_project(project_root: Path, canonical_rows: int = 3) -> tuple[Path, Path]:
    project_root.mkdir(parents=True)
    canonical_pairlist = project_root / "pairlist.csv"
    pd.DataFrame([_pair_row(index) for index in range(1, canonical_rows + 1)], columns=PAIR_COLUMNS).to_csv(
        canonical_pairlist,
        index=False,
    )
    curation_state = project_root / "metadata" / "pair_curation_state.json"
    curation_state.parent.mkdir(parents=True)
    curation_state.write_text(json.dumps({"rounds": ["round_001"]}, indent=2), encoding="utf-8")
    (project_root / "project_manifest.json").write_text(
        json.dumps(
            {
                "engines": ["vina"],
                "pairlist_file": str(canonical_pairlist),
                "latest_pair_round": "round_001",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return canonical_pairlist, curation_state


def test_deploy_pairlist_file_uses_exact_rows_and_preserves_canonical_state(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    project_root = tmp_path / "project"
    canonical_pairlist, curation_state = _write_project(project_root)
    exact_pairlist = tmp_path / "exact_pairs.csv"
    pd.DataFrame([_pair_row(7), _pair_row(9)], columns=PAIR_COLUMNS).to_csv(exact_pairlist, index=False)

    canonical_before = canonical_pairlist.read_bytes()
    curation_before = curation_state.read_bytes()
    captured: dict[str, object] = {}

    def _capture_deployment(**kwargs):
        captured.update(kwargs)
        return {
            "generation_stage_root": str(project_root / "deployments" / "round_001" / "screen"),
            "stage_root": str(project_root / "deployments" / "round_001" / "screen"),
            "engine_jobs": {"vina": {"planned": len(kwargs["pairlist_rows"]), "skipped": 0}},
        }

    monkeypatch.setattr(docking_cli, "generate_slurm_deployment", _capture_deployment)

    exit_code = docking_cli.deploy_main(
        [
            "--project-dir",
            str(project_root),
            "--seed",
            "12345",
            "--pairlist-file",
            str(exact_pairlist),
            "--engines",
            "vina",
            "--no-ligand-qc-gate",
            "--no-receptor-qc-gate",
        ]
    )

    assert exit_code == 0
    rows = captured["pairlist_rows"]
    # Spec 036 R6: tags use the receptor file stem (no .pdbqt).
    assert [row.tag for row in rows] == [
        "R7_site_1_L7.pdbqt",
        "R9_site_1_L9.pdbqt",
    ]
    assert captured["pair_source_file"] == str(exact_pairlist.resolve())
    assert captured["pair_source_kind"] == "exact_pairlist"
    assert captured["pair_source_sha256"] == hashlib.sha256(exact_pairlist.read_bytes()).hexdigest()
    assert captured["rerun_manifest_file"] == ""
    assert canonical_pairlist.read_bytes() == canonical_before
    assert curation_state.read_bytes() == curation_before


def test_pairlist_file_and_rerun_manifest_are_mutually_exclusive(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as exc_info:
        docking_cli.deploy_main(
            [
                "--project-dir",
                str(tmp_path / "project"),
                "--pairlist-file",
                str(tmp_path / "pairs.csv"),
                "--from-rerun-manifest",
                str(tmp_path / "rerun.csv"),
            ]
        )
    assert exc_info.value.code == 2


def test_workflow_deploy_forwards_exact_pairlist_and_enforces_exclusion() -> None:
    parser = argparse.ArgumentParser()
    _add_deploy_args(parser)
    args = parser.parse_args(["--project-dir", "project", "--pairlist-file", "subset.csv"])

    argv = _build_deploy_argv(args)

    assert argv[argv.index("--pairlist-file") + 1] == "subset.csv"
    assert "--from-rerun-manifest" not in argv
    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "--project-dir",
                "project",
                "--pairlist-file",
                "subset.csv",
                "--from-rerun-manifest",
                "rerun.csv",
            ]
        )


@pytest.mark.parametrize("scheduler", ["slurm", "condor"])
def test_deployment_manifests_include_pair_source_provenance(
    scheduler: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class _Runner:
        @staticmethod
        def plan_jobs(_rows, *, skip_completed=False):
            return []

    monkeypatch.setattr(deployment, "build_runner", lambda *_args, **_kwargs: _Runner())
    project_root = tmp_path / "project"
    pair_source = tmp_path / "exact.csv"
    pair_source.write_text("receptor,site_id,ligand\nR1,site_1,L1\n", encoding="utf-8")
    digest = hashlib.sha256(pair_source.read_bytes()).hexdigest()
    row = PairlistRow(
        receptor="R1.pdbqt",
        site_id="site_1",
        ligand="L1.pdbqt",
        center_x=1.0,
        center_y=2.0,
        center_z=3.0,
        size_x=20.0,
        size_y=20.0,
        size_z=20.0,
    )

    common = {
        "project_root": project_root,
        "engines": ["vina"],
        "pairlist_rows": [row],
        "runtime_by_engine": {"vina": {}},
        "round_id": "round_exact",
        "docking_mode": "screen",
        "pair_source_file": str(pair_source),
        "pair_source_sha256": digest,
        "pair_source_kind": "exact_pairlist",
    }
    if scheduler == "slurm":
        summary = deployment.generate_slurm_deployment(
            **common,
            slurm_options={"partition": "test"},
        )
    else:
        summary = deployment.generate_condor_deployment(**common)

    assert summary["pair_count"] == 1
    assert summary["pair_source_sha256"] == digest
    assert summary["pair_source_kind"] == "exact_pairlist"
    project_manifest = json.loads(
        (Path(summary["generation_stage_root"]) / "deployment_manifest.json").read_text(encoding="utf-8")
    )
    engine_manifest = json.loads(
        (
            Path(summary["generation_stage_root"])
            / "vina"
            / "manifest"
            / "deployment_manifest.json"
        ).read_text(encoding="utf-8")
    )
    for manifest in (project_manifest, engine_manifest):
        assert manifest["pair_source_file"] == str(pair_source.resolve())
        assert manifest["pair_source_sha256"] == digest
        assert manifest["pair_source_kind"] == "exact_pairlist"


def test_console_entry_point_targets_packaged_workflow_cli() -> None:
    setup_text = (Path(__file__).resolve().parents[1] / "setup.py").read_text(encoding="utf-8")
    assert "pdb-prepare-wizard=workflow.cli:main" in setup_text
    assert "pdb-prepare-wizard=main:main" not in setup_text

    from workflow.cli import main

    assert callable(main)


def test_direct_runtime_imports_are_declared() -> None:
    root = Path(__file__).resolve().parents[1]

    def _names(path: Path) -> set[str]:
        return {
            line.split(";", 1)[0].split("#", 1)[0].strip().split(">=", 1)[0].lower()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }

    root_requirements = _names(root / "requirements.txt")
    post_requirements = _names(root / "post_docking_analysis" / "requirements.txt")
    assert {"pyyaml", "scipy", "scikit-learn"}.issubset(root_requirements)
    assert {"pyyaml", "scipy", "scikit-learn"}.issubset(post_requirements)
