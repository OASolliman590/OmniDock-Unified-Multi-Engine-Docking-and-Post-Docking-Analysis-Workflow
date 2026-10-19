"""Spec 037 R4a: prepare once.

An existing prepared output is reused (``reused_unchanged``) only when its provenance matches the
current resolved policy hash, the request fingerprint, and the input and output SHA-256. The tests count
backend invocations (the PDB2PQR, Meeko and ligand-conversion calls), so "prepared exactly once" is
shown by the counts and not only by the status field.

Receptor tests use a fake Meeko (and fake PDB2PQR and Open Babel for the PDB2PQR path) on PATH.
Ligand tests run real RDKit and Open Babel preparation on the STI fixture (about 1 s per ligand).
The shell test runs ``prep_autodock_enhanced.sh`` end to end with counting wrappers on PATH.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from docking.preparation import ligand_preparation as lp
from docking.preparation import prepared_reuse
from docking.preparation import receptor_preparation as rp
from docking.preparation.receptor_preparation import prepare_receptor
from docking.project_layout import set_protonation_override, set_protonation_policy
from test_spec033_receptor_prep import GOOD_PDBQT_LINES, SOURCE_PDB, _config, _write
from test_spec036_protonation_policy import _write_manifest

REPO_ROOT = Path(__file__).resolve().parents[1]
STI_SDF = REPO_ROOT / "test" / "fixtures" / "spec032" / "phase2_1iep_inputs" / "STI_1IEP_A201.sdf"
SHELL_SCRIPT = REPO_ROOT / "prep_autodock_enhanced.sh"

try:
    import rdkit  # noqa: F401

    HAS_RDKIT = True
except ImportError:  # pragma: no cover - environment dependent
    HAS_RDKIT = False

HAS_OBABEL = shutil.which("obabel") is not None
HAS_JQ = shutil.which("jq") is not None


# ----------------------------------------------------------------------------- receptor helpers


def _counting_receptor_bin(tmp_path: Path, monkeypatch) -> Path:
    """PATH with fake PDB2PQR, Open Babel and Meeko. Each invocation appends one line to a count log.

    Meeko copies GOOD_PDBQT_LINES to its -p output; PDB2PQR and Open Babel copy their input to the output.
    """
    bin_dir = tmp_path / "counting_bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    log = tmp_path / "backend_invocations.log"
    log.write_text("", encoding="utf-8")
    good = tmp_path / "good_meeko_output.pdbqt"
    good.write_text("\n".join(GOOD_PDBQT_LINES) + "\n", encoding="utf-8")
    # Shell builtins only: the receptor tests restrict PATH to this folder. --help probes are not counted.
    copy = (
        "copy_file() { while IFS= read -r line || [ -n \"$line\" ]; do printf '%s\\n' \"$line\"; done < \"$1\" > \"$2\"; }\n"
    )
    scripts = {
        "mk_prepare_receptor.py": (
            "[ \"$1\" = \"--help\" ] && exit 0\n"
            "out=\"\"\n"
            "while [ $# -gt 0 ]; do\n"
            "  if [ \"$1\" = \"-p\" ]; then out=\"$2\"; fi\n"
            "  shift\n"
            "done\n"
            "echo meeko >> \"$BACKEND_LOG\"\n"
            "copy_file \"$GOOD_MEEKO_OUTPUT\" \"$out\"\n"
        ),
        "pdb2pqr30": (
            "n=$#\n"
            "eval src=\\${$((n-1))}\n"
            "eval dst=\\${$n}\n"
            "echo pdb2pqr >> \"$BACKEND_LOG\"\n"
            "copy_file \"$src\" \"$dst\"\n"
        ),
        "obabel": (
            "in=\"$1\"; out=\"\"\n"
            "while [ $# -gt 0 ]; do\n"
            "  if [ \"$1\" = \"-O\" ]; then out=\"$2\"; fi\n"
            "  shift\n"
            "done\n"
            "echo obabel >> \"$BACKEND_LOG\"\n"
            "copy_file \"$in\" \"$out\"\n"
        ),
    }
    for name, body in scripts.items():
        script = bin_dir / name
        script.write_text("#!/bin/sh\n" + copy + body, encoding="utf-8")
        script.chmod(0o755)
    monkeypatch.setenv("BACKEND_LOG", str(log))
    monkeypatch.setenv("GOOD_MEEKO_OUTPUT", str(good))
    monkeypatch.setenv("PATH", str(bin_dir))
    return log


def _invocations(log: Path, backend: str = "meeko") -> int:
    return sum(1 for line in log.read_text(encoding="utf-8").splitlines() if line.strip() == backend)


def _receptor_setup(tmp_path: Path) -> tuple:
    source = _write(tmp_path / "raw" / "1IEP.pdb", SOURCE_PDB)
    destination = tmp_path / "prepared" / "1IEP.pdbqt"
    return source, destination


# ----------------------------------------------------------------------------- receptor: prepare once


def test_receptor_second_run_reuses_and_meeko_runs_once(tmp_path, monkeypatch):
    log = _counting_receptor_bin(tmp_path, monkeypatch)
    source, destination = _receptor_setup(tmp_path)
    config = _config(receptor_use_pdb2pqr=False)

    first = prepare_receptor(source, destination, config)
    second = prepare_receptor(source, destination, config)

    assert first["preparation_status"] == "prepared"
    assert second["status"] == "reused_unchanged"
    assert second["preparation_status"] == "reused_unchanged"
    assert _invocations(log) == 1, "Meeko must run exactly once for an unchanged receptor"


def test_receptor_changed_input_is_prepared_again(tmp_path, monkeypatch):
    log = _counting_receptor_bin(tmp_path, monkeypatch)
    source, destination = _receptor_setup(tmp_path)
    config = _config(receptor_use_pdb2pqr=False)
    prepare_receptor(source, destination, config)

    source.write_text(SOURCE_PDB + "\nREMARK edited input\n", encoding="utf-8")
    again = prepare_receptor(source, destination, config)

    assert again["preparation_status"] == "prepared"
    assert _invocations(log) == 2
    assert again["input_sha256"] == prepared_reuse.sha256_file(source)


def test_receptor_tampered_output_is_prepared_again(tmp_path, monkeypatch):
    log = _counting_receptor_bin(tmp_path, monkeypatch)
    source, destination = _receptor_setup(tmp_path)
    config = _config(receptor_use_pdb2pqr=False)
    prepare_receptor(source, destination, config)

    destination.write_text(destination.read_text(encoding="utf-8") + "REMARK tampered\n", encoding="utf-8")
    again = prepare_receptor(source, destination, config)

    assert again["preparation_status"] == "prepared"
    assert _invocations(log) == 2


def test_receptor_force_prepares_again_even_when_reusable(tmp_path, monkeypatch):
    log = _counting_receptor_bin(tmp_path, monkeypatch)
    source, destination = _receptor_setup(tmp_path)
    config = _config(receptor_use_pdb2pqr=False)
    prepare_receptor(source, destination, config)

    forced = prepare_receptor(source, destination, config, force=True)
    via_config = prepare_receptor(source, destination, _config(receptor_use_pdb2pqr=False, force=True))

    assert forced["preparation_status"] == "prepared"
    assert via_config["preparation_status"] == "prepared"
    assert _invocations(log) == 3


def test_receptor_changed_request_pH_is_prepared_again(tmp_path, monkeypatch):
    log = _counting_receptor_bin(tmp_path, monkeypatch)
    source, destination = _receptor_setup(tmp_path)
    prepare_receptor(source, destination, _config(ph=7.4))
    assert _invocations(log, "pdb2pqr") == 1

    again = prepare_receptor(source, destination, _config(ph=6.0))

    assert again["preparation_status"] == "prepared"
    assert _invocations(log, "pdb2pqr") == 2
    assert again["requested_ph"] == 6.0


def test_receptor_changed_project_policy_is_prepared_again(tmp_path, monkeypatch):
    log = _counting_receptor_bin(tmp_path, monkeypatch)
    project = tmp_path / "project"
    _write_manifest(project)
    set_protonation_policy(project, receptor_ph=7.4, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.4)
    source, destination = _receptor_setup(tmp_path)
    config = _config(ph=7.4, project_dir=str(project))
    first = prepare_receptor(source, destination, config)
    assert first["policy_level"] == "project"
    assert prepare_receptor(source, destination, config)["status"] == "reused_unchanged"

    # The override changes the resolved policy hash for this receptor: the record no longer matches.
    set_protonation_override(project, kind="receptor", name="1IEP", receptor_ph=6.0)
    again = prepare_receptor(source, destination, config)

    assert again["policy_level"] == "receptor"
    assert again["preparation_status"] == "prepared"
    assert again["policy_hash"] != first["policy_hash"]
    assert _invocations(log, "pdb2pqr") == 2


def test_receptor_policy_hash_mismatch_alone_is_refused_for_reuse():
    record = {"policy_hash": "old", "request_fingerprint": "same", "input_sha256": "i", "output_sha256": "o"}
    reason = prepared_reuse.reuse_block_reason(
        record, output=Path(__file__), input_file=Path(__file__), policy_hash="new", fingerprint="same", force=False
    )
    assert reason == "policy_changed"


def test_record_written_before_prepare_once_is_prepared_again(tmp_path, monkeypatch):
    log = _counting_receptor_bin(tmp_path, monkeypatch)
    source, destination = _receptor_setup(tmp_path)
    prepare_receptor(source, destination, _config(receptor_use_pdb2pqr=False))
    provenance = destination.parent / (destination.name + ".preparation.json")
    payload = json.loads(provenance.read_text(encoding="utf-8"))
    payload.pop("request_fingerprint")
    provenance.write_text(json.dumps(payload), encoding="utf-8")

    again = prepare_receptor(source, destination, _config(receptor_use_pdb2pqr=False))

    assert again["preparation_status"] == "prepared"
    assert _invocations(log) == 2


# ----------------------------------------------------------------------------- ligand: prepare once


def _ligand_env(monkeypatch, report_dir: Path, ph: str = "7.4") -> None:
    monkeypatch.delenv("PDBWIZARD_PROJECT_DIR", raising=False)
    monkeypatch.setenv("PDBWIZARD_LIGAND_PREP_PROFILE", "engine_aware_full")
    monkeypatch.setenv("PDBWIZARD_SELECTED_ENGINES", "vina")
    monkeypatch.setenv("PDBWIZARD_LIGAND_PREP_PROTONATION_POLICY", "ph_model")
    monkeypatch.setenv("PDBWIZARD_LIGAND_PREP_PH", ph)
    monkeypatch.setenv("PDBWIZARD_LIGAND_PREP_PH_SOURCE", "user_entered")
    monkeypatch.setenv("PDBWIZARD_LIGAND_PREP_REPORT_DIR", str(report_dir))
    monkeypatch.delenv("PDBWIZARD_LIGAND_PREP_STATE_MAP", raising=False)


@pytest.fixture
def ligand_backend_counter(monkeypatch):
    """Count the ligand conversion backend (Meeko PDBQT conversion) without changing its behaviour."""
    calls = {"count": 0}
    original = lp._prepare_with_meeko

    def counted(input_sdf, output_pdbqt):
        calls["count"] += 1
        return original(input_sdf, output_pdbqt)

    monkeypatch.setattr(lp, "_prepare_with_meeko", counted)
    return calls


@pytest.mark.skipif(not HAS_RDKIT, reason="skipped_missing_dependency: needs RDKit")
def test_ligand_prepared_once_then_reused_then_changed_input_refused(tmp_path, monkeypatch, ligand_backend_counter):
    raw = tmp_path / "raw" / "STI.sdf"
    raw.parent.mkdir(parents=True)
    shutil.copy(STI_SDF, raw)
    report_dir = tmp_path / "steps"
    output = tmp_path / "prepared" / "STI.pdbqt"
    _ligand_env(monkeypatch, report_dir)

    first = lp.prepare_ligand_for_vina_family(raw, output)
    second = lp.prepare_ligand_for_vina_family(raw, output)
    assert first["preparation_status"] == "prepared"
    assert second["preparation_status"] == "reused_unchanged"
    assert second["status"] == "reused_unchanged"
    assert ligand_backend_counter["count"] == 1, "ligand conversion must run exactly once for an unchanged input"

    raw.write_bytes(raw.read_bytes() + b"\n")
    third = lp.prepare_ligand_for_vina_family(raw, output)
    assert third["preparation_status"] == "prepared"
    assert third["reuse_block_reason"] == "input_changed"
    assert ligand_backend_counter["count"] == 2


def test_ligand_tampered_output_and_force_are_prepared_again(tmp_path, monkeypatch, ligand_backend_counter):
    raw = tmp_path / "raw" / "STI.sdf"
    raw.parent.mkdir(parents=True)
    shutil.copy(STI_SDF, raw)
    report_dir = tmp_path / "steps"
    output = tmp_path / "prepared" / "STI.pdbqt"
    _ligand_env(monkeypatch, report_dir)
    lp.prepare_ligand_for_vina_family(raw, output)

    output.write_text(output.read_text(encoding="utf-8") + "REMARK tampered\n", encoding="utf-8")
    tampered = lp.prepare_ligand_for_vina_family(raw, output)
    assert tampered["preparation_status"] == "prepared"
    assert tampered["reuse_block_reason"] == "output_changed"
    assert ligand_backend_counter["count"] == 2

    forced = lp.prepare_ligand_for_vina_family(raw, output, force=True)
    assert forced["preparation_status"] == "prepared"
    assert forced["reuse_block_reason"] == "forced"
    assert ligand_backend_counter["count"] == 3


def test_ligand_changed_pH_without_project_policy_is_prepared_again(tmp_path, monkeypatch, ligand_backend_counter):
    raw = tmp_path / "raw" / "STI.sdf"
    raw.parent.mkdir(parents=True)
    shutil.copy(STI_SDF, raw)
    report_dir = tmp_path / "steps"
    output = tmp_path / "prepared" / "STI.pdbqt"
    _ligand_env(monkeypatch, report_dir, ph="7.4")
    lp.prepare_ligand_for_vina_family(raw, output)

    _ligand_env(monkeypatch, report_dir, ph="6.0")
    changed = lp.prepare_ligand_for_vina_family(raw, output)

    assert changed["preparation_status"] == "prepared"
    assert changed["reuse_block_reason"] == "request_changed"
    assert ligand_backend_counter["count"] == 2


def test_ligand_changed_project_override_is_prepared_again(tmp_path, monkeypatch, ligand_backend_counter):
    project = tmp_path / "project"
    _write_manifest(project)
    set_protonation_policy(project, receptor_ph=7.4, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.4)
    raw = tmp_path / "raw" / "NEUT.sdf"
    raw.parent.mkdir(parents=True)
    shutil.copy(STI_SDF, raw)
    report_dir = tmp_path / "steps"
    output = tmp_path / "prepared" / "NEUT.pdbqt"
    _ligand_env(monkeypatch, report_dir, ph="7.4")
    monkeypatch.setenv("PDBWIZARD_PROJECT_DIR", str(project))
    first = lp.prepare_ligand_for_vina_family(raw, output)
    assert first["policy_level"] == "project"
    assert lp.prepare_ligand_for_vina_family(raw, output)["preparation_status"] == "reused_unchanged"

    set_protonation_override(project, kind="ligand", name="NEUT", ligand_ph=6.0)
    changed = lp.prepare_ligand_for_vina_family(raw, output)

    assert changed["policy_level"] == "ligand"
    assert changed["policy_hash"] != first["policy_hash"]
    assert changed["preparation_status"] == "prepared"
    assert ligand_backend_counter["count"] == 2


# ----------------------------------------------------------------------------- shell loop: end to end


def _real_tool(name: str) -> str:
    found = shutil.which(name)
    assert found, f"{name} is required for the shell test"
    return found


@pytest.mark.skipif(not (HAS_RDKIT and HAS_OBABEL and HAS_JQ), reason="skipped_missing_dependency: RDKit, Open Babel and jq")
def test_shell_loop_prepares_once_reuses_then_reprepares_tampered_and_forced(tmp_path, monkeypatch):
    real_ligand_tool = _real_tool("mk_prepare_ligand.py")
    project = tmp_path / "shell_project"
    raw_lig = project / "raw_ligands"
    raw_rec = project / "raw_receptors"
    raw_lig.mkdir(parents=True)
    raw_rec.mkdir(parents=True)
    shutil.copy(STI_SDF, raw_lig / "STI.sdf")
    _write(raw_rec / "1IEP.pdb", SOURCE_PDB)
    prepared_lig = project / "ligands_prep"
    prepared_rec = project / "receptors_prep"
    config_path = project / ".meta" / "config.json"
    config_path.parent.mkdir(parents=True)

    def write_config(force: bool) -> None:
        config = {
            "input": {
                "ligands": {"path": str(raw_lig), "formats": ["sdf"], "in_same_folder": False},
                "receptors": {"path": str(raw_rec), "formats": ["pdb"], "in_same_folder": False},
            },
            "output": {"ligands": str(prepared_lig), "receptors": str(prepared_rec), "logs": "logs"},
            "preparation": {
                "force_field": "AMBER", "force_field_source": "user_entered",
                "ph": 7.4, "ph_source": "user_entered",
                "protonation_policy": "ph_model", "protonation_state_map": "",
                "allow_bad_res": False, "receptor_use_pdb2pqr": False,
                "ligand_preparation_profile": "engine_aware_full",
                "ligand_preparation_backend": "engine_aware_full",
                "selected_engines": ["vina"], "force": force,
            },
            "quality_control": {"validate_outputs": True, "check_file_sizes": True, "min_file_size_kb": 1},
        }
        config_path.write_text(json.dumps(config), encoding="utf-8")

    bin_dir = tmp_path / "shell_bin"
    bin_dir.mkdir()
    log = tmp_path / "shell_backend.log"
    log.write_text("", encoding="utf-8")
    good = tmp_path / "good_meeko_output.pdbqt"
    # The shell loop rejects outputs under 1 KB (validate_file), so the fixture is padded with REMARK lines.
    padding = [f"REMARK padding line {index:04d} for the 1 KB output check" for index in range(24)]
    good.write_text("\n".join(padding + GOOD_PDBQT_LINES) + "\n", encoding="utf-8")
    fake_receptor = bin_dir / "mk_prepare_receptor.py"
    fake_receptor.write_text(
        "#!/bin/sh\n"
        "[ \"$1\" = \"--help\" ] && exit 0\n"
        "out=\"\"\n"
        "while [ $# -gt 0 ]; do\n"
        "  if [ \"$1\" = \"-p\" ]; then out=\"$2\"; fi\n"
        "  shift\n"
        "done\n"
        f"echo receptor >> '{log}'\n"
        f"cp '{good}' \"$out\"\n",
        encoding="utf-8",
    )
    counting_ligand = bin_dir / "mk_prepare_ligand.py"
    counting_ligand.write_text(
        "#!/bin/sh\n"
        "[ \"$1\" = \"--help\" ] && exec '" + real_ligand_tool + "' \"$@\"\n"
        f"echo ligand >> '{log}'\n"
        "exec '" + real_ligand_tool + "' \"$@\"\n",
        encoding="utf-8",
    )
    for tool in (fake_receptor, counting_ligand):
        tool.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("PYTHONPATH", str(REPO_ROOT))

    def run_loop(force: bool = False) -> subprocess.CompletedProcess:
        write_config(force)
        return subprocess.run(
            ["bash", str(SHELL_SCRIPT), str(config_path)],
            capture_output=True,
            text=True,
            cwd=str(tmp_path),
            timeout=600,
        )

    def counts() -> tuple:
        lines = log.read_text(encoding="utf-8").splitlines()
        return lines.count("ligand"), lines.count("receptor")

    first = run_loop()
    assert first.returncode == 0, first.stdout[-2000:] + first.stderr[-2000:]
    assert counts() == (1, 1)
    assert (prepared_lig / "STI.pdbqt").is_file() and (prepared_rec / "1IEP.pdbqt").is_file()

    second = run_loop()
    assert second.returncode == 0, second.stdout[-2000:]
    assert counts() == (1, 1), "an unchanged run must not invoke any backend"
    ligand_skips = (prepared_lig / "ligand_preparation_skipped.csv").read_text(encoding="utf-8")
    receptor_skips = (prepared_rec / "receptor_preparation_skipped.csv").read_text(encoding="utf-8")
    assert '"skipped","reused_unchanged"' in ligand_skips
    assert '"skipped","reused_unchanged"' in receptor_skips
    summary = (prepared_rec / "preparation_summary.txt").read_text(encoding="utf-8")
    assert "Reused unchanged" in summary and "Ligands: 1" in summary and "Receptors: 1" in summary

    (prepared_rec / "1IEP.pdbqt").write_text(
        (prepared_rec / "1IEP.pdbqt").read_text(encoding="utf-8") + "REMARK tampered\n", encoding="utf-8"
    )
    third = run_loop()
    assert third.returncode == 0, third.stdout[-2000:]
    assert counts() == (1, 2), "only the tampered receptor is prepared again"

    forced = run_loop(force=True)
    assert forced.returncode == 0, forced.stdout[-2000:]
    assert counts() == (2, 3), "--force (config preparation.force) prepares every entity again"
