"""Spec 033 R2/R3: strict receptor preparation and project-local preparation config.

Unit tests use synthetic files and fake backend executables placed on PATH.
Dependency-backed tests run the real PDB2PQR, Open Babel and Meeko chain on the
1IEP chain-A protein-only receptor and skip with a reason when those tools are
not installed.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from docking.preparation.receptor_preparation import (
    ReceptorPreparationError,
    polar_hydrogen_gate,
    prepare_receptor,
    receptor_backend_capabilities,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_1IEP = REPO_ROOT / "test" / "fixtures" / "spec032" / "phase2_1iep_inputs" / "1IEP.pdb"

SOURCE_PDB = "\n".join(
    [
        "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N",
        "ATOM      2  CA  GLY A   1       1.450   0.000   0.000  1.00  0.00           C",
        "ATOM      3  C   GLY A   1       2.000   1.400   0.000  1.00  0.00           C",
        "ATOM      4  O   GLY A   1       1.300   2.400   0.000  1.00  0.00           O",
        "END",
    ]
) + "\n"


def _pdbqt_line(serial: int, name: str, xyz, atom_type: str) -> str:
    x, y, z = xyz
    return (
        f"ATOM  {serial:5d} {name:<4s} GLY A   1    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00    {0.000:6.3f} {atom_type}"
    )


GOOD_PDBQT_LINES = [
    _pdbqt_line(1, "N", (0.0, 0.0, 0.0), "N"),
    _pdbqt_line(2, "CA", (1.45, 0.0, 0.0), "C"),
    _pdbqt_line(3, "C", (2.0, 1.4, 0.0), "C"),
    _pdbqt_line(4, "O", (1.3, 2.4, 0.0), "OA"),
    _pdbqt_line(5, "HN", (0.0, 0.8, 0.6), "HD"),  # 1.0 A from N: a polar hydrogen on the backbone N
]


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _config(**preparation) -> dict:
    base = {"force_field": "AMBER", "ph": 7.4, "allow_bad_res": False, "ligand_preparation_profile": "engine_aware_full"}
    base.update(preparation)
    return {"preparation": base}


def _fake_meeko_bin(tmp_path: Path, output_lines, monkeypatch) -> Path:
    """Put a fake mk_prepare_receptor.py first on PATH. It copies FAKE_MEEKO_OUTPUT to the -p target.

    The script uses only shell builtins, so PATH can be restricted to this folder.
    """
    bin_dir = tmp_path / "fake_bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    fake_output = _write(tmp_path / "fake_meeko_output.pdbqt", "\n".join(output_lines) + "\n")
    script = bin_dir / "mk_prepare_receptor.py"
    script.write_text(
        "#!/bin/sh\n"
        "out=\"\"\n"
        "while [ $# -gt 0 ]; do\n"
        "  if [ \"$1\" = \"-p\" ]; then out=\"$2\"; fi\n"
        "  shift\n"
        "done\n"
        "while IFS= read -r line || [ -n \"$line\" ]; do printf '%s\\n' \"$line\"; done < \"$FAKE_MEEKO_OUTPUT\" > \"$out\"\n",
        encoding="utf-8",
    )
    script.chmod(0o755)
    monkeypatch.setenv("FAKE_MEEKO_OUTPUT", str(fake_output))
    monkeypatch.setenv("PATH", str(bin_dir))
    return bin_dir


def _have(executable: str) -> bool:
    return shutil.which(executable) is not None


REQUIRES_1IEP_TOOLS = pytest.mark.skipif(
    not (_have("pdb2pqr30") and _have("obabel") and _have("mk_prepare_receptor.py")),
    reason="needs pdb2pqr30, obabel and mk_prepare_receptor.py on PATH (conda-forge pdb2pqr and meeko)",
)


# ---------------------------------------------------------------- polar-H gate (R2d)


def test_polar_hydrogen_gate_rejects_receptor_without_hd(tmp_path):
    pdbqt = _write(tmp_path / "no_hd.pdbqt", "\n".join(GOOD_PDBQT_LINES[:4]) + "\n")
    gate = polar_hydrogen_gate(pdbqt)
    assert gate["passed"] is False
    assert gate["hd_count"] == 0
    assert gate["hd_near_nitrogen_count"] == 0


def test_polar_hydrogen_gate_accepts_hd_near_backbone_nitrogen(tmp_path):
    pdbqt = _write(tmp_path / "with_hd.pdbqt", "\n".join(GOOD_PDBQT_LINES) + "\n")
    gate = polar_hydrogen_gate(pdbqt)
    assert gate["passed"] is True
    assert gate["hd_count"] == 1
    assert gate["hd_near_nitrogen_count"] == 1


def test_polar_hydrogen_gate_rejects_hd_far_from_any_nitrogen(tmp_path):
    far_hd = _pdbqt_line(5, "HX", (5.0, 5.0, 5.0), "HD")
    pdbqt = _write(tmp_path / "far_hd.pdbqt", "\n".join(GOOD_PDBQT_LINES[:4] + [far_hd]) + "\n")
    gate = polar_hydrogen_gate(pdbqt)
    assert gate["passed"] is False
    assert gate["hd_count"] == 1
    assert gate["hd_near_nitrogen_count"] == 0


def test_receptor_without_hd_fails_gate_with_explicit_reason(tmp_path, monkeypatch):
    source = _write(tmp_path / "src" / "receptor.pdb", SOURCE_PDB)
    destination = tmp_path / "out" / "receptor.pdbqt"
    _fake_meeko_bin(tmp_path, GOOD_PDBQT_LINES[:4], monkeypatch)
    with pytest.raises(ReceptorPreparationError) as excinfo:
        prepare_receptor(source, destination, _config(receptor_use_pdb2pqr=False))
    assert excinfo.value.reason == "polar_hydrogen_gate_failed"
    assert not destination.exists()


# ---------------------------------------------------------------- heavy-atom conservation (R2a)


def test_synthetic_heavy_atom_loss_is_detected_and_named(tmp_path, monkeypatch):
    source = _write(tmp_path / "src" / "receptor.pdb", SOURCE_PDB)
    destination = tmp_path / "out" / "receptor.pdbqt"
    _fake_meeko_bin(tmp_path, GOOD_PDBQT_LINES[:3] + [GOOD_PDBQT_LINES[4]], monkeypatch)  # drops O
    with pytest.raises(ReceptorPreparationError) as excinfo:
        prepare_receptor(source, destination, _config(receptor_use_pdb2pqr=False))
    assert excinfo.value.reason == "heavy_atom_conservation_failed"
    assert "GLY 1:O" in excinfo.value.details
    assert not destination.exists()


def test_synthetic_heavy_atom_move_is_detected(tmp_path, monkeypatch):
    moved_o = _pdbqt_line(4, "O", (1.3, 2.9, 0.0), "OA")  # 0.5 A away from the source O
    source = _write(tmp_path / "src" / "receptor.pdb", SOURCE_PDB)
    destination = tmp_path / "out" / "receptor.pdbqt"
    _fake_meeko_bin(tmp_path, GOOD_PDBQT_LINES[:3] + [moved_o, GOOD_PDBQT_LINES[4]], monkeypatch)
    with pytest.raises(ReceptorPreparationError) as excinfo:
        prepare_receptor(source, destination, _config(receptor_use_pdb2pqr=False))
    assert excinfo.value.reason == "heavy_atom_conservation_failed"
    assert not destination.exists()


# ---------------------------------------------------------------- backends on PATH (R2b, R2c)


def test_missing_meeko_on_path_fails_explicitly_without_output(tmp_path, monkeypatch):
    empty_bin = tmp_path / "empty_bin"
    empty_bin.mkdir()
    monkeypatch.setenv("PATH", str(empty_bin))
    source = _write(tmp_path / "src" / "receptor.pdb", SOURCE_PDB)
    destination = tmp_path / "out" / "receptor.pdbqt"
    with pytest.raises(ReceptorPreparationError) as excinfo:
        prepare_receptor(source, destination, _config(receptor_use_pdb2pqr=False))
    assert excinfo.value.reason == "backend_unavailable"
    assert "mk_prepare_receptor.py" in excinfo.value.details
    assert not destination.exists()


def test_missing_pdb2pqr_on_path_fails_by_default_without_output(tmp_path, monkeypatch):
    # Meeko is present, PDB2PQR is not. The default configuration (no receptor_use_pdb2pqr key)
    # requires PDB2PQR, so the receptor must fail rather than silently skip protonation.
    _fake_meeko_bin(tmp_path, GOOD_PDBQT_LINES, monkeypatch)
    source = _write(tmp_path / "src" / "receptor.pdb", SOURCE_PDB)
    destination = tmp_path / "out" / "receptor.pdbqt"
    with pytest.raises(ReceptorPreparationError) as excinfo:
        prepare_receptor(source, destination, _config())
    assert excinfo.value.reason == "backend_unavailable"
    assert "pdb2pqr30" in excinfo.value.details
    assert not destination.exists()


def test_explicit_false_uses_template_protonation_and_records_provenance(tmp_path, monkeypatch):
    _fake_meeko_bin(tmp_path, GOOD_PDBQT_LINES, monkeypatch)
    source = _write(tmp_path / "src" / "receptor.pdb", SOURCE_PDB)
    destination = tmp_path / "out" / "receptor.pdbqt"
    payload = prepare_receptor(source, destination, _config(receptor_use_pdb2pqr=False))
    assert payload["protonation_status"] == "template_selected_not_ph_titrated"
    assert payload["hd_count"] == 1 and payload["hd_near_nitrogen_count"] == 1
    assert payload["heavy_atom_conservation"] == "passed"
    provenance = json.loads((tmp_path / "out" / "receptor.pdbqt.preparation.json").read_text(encoding="utf-8"))
    assert provenance["input_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert provenance["output_sha256"] == hashlib.sha256(destination.read_bytes()).hexdigest()
    assert set(provenance["tool_versions"]) == {"pdb2pqr", "meeko"}
    assert provenance["backend"] == "meeko"
    assert (tmp_path / "out" / "receptor.pdbqt.preparation.log").exists()


@pytest.mark.parametrize(
    "preparation, reason",
    [
        ({"receptor_use_pdb2pqr": "false"}, "invalid_configuration"),  # a string is not a disable switch
        ({"receptor_use_pdb2pqr": 0}, "invalid_configuration"),
        ({"allow_bad_res": True, "receptor_use_pdb2pqr": False}, "allow_bad_res_unsupported"),
    ],
)
def test_ambiguous_or_permissive_preparation_settings_are_refused(tmp_path, monkeypatch, preparation, reason):
    _fake_meeko_bin(tmp_path, GOOD_PDBQT_LINES, monkeypatch)
    source = _write(tmp_path / "src" / "receptor.pdb", SOURCE_PDB)
    destination = tmp_path / "out" / "receptor.pdbqt"
    with pytest.raises(ReceptorPreparationError) as excinfo:
        prepare_receptor(source, destination, _config(**preparation))
    assert excinfo.value.reason == reason
    assert not destination.exists()


def test_failed_run_removes_stale_receptor_output(tmp_path, monkeypatch):
    destination = _write(tmp_path / "out" / "receptor.pdbqt", "\n".join(GOOD_PDBQT_LINES) + "\n")
    source = _write(tmp_path / "src" / "receptor.pdb", SOURCE_PDB)
    empty_bin = tmp_path / "empty_bin"
    empty_bin.mkdir()
    monkeypatch.setenv("PATH", str(empty_bin))
    with pytest.raises(ReceptorPreparationError):
        prepare_receptor(source, destination, _config(receptor_use_pdb2pqr=False))
    assert not destination.exists()


def test_capability_check_reports_none_when_backends_absent(tmp_path, monkeypatch):
    empty_bin = tmp_path / "empty_bin"
    empty_bin.mkdir()
    monkeypatch.setenv("PATH", str(empty_bin))
    assert receptor_backend_capabilities() == {"pdb2pqr": None, "meeko": None}


def test_capability_check_reports_versions_when_present(tmp_path, monkeypatch):
    _fake_meeko_bin(tmp_path, GOOD_PDBQT_LINES, monkeypatch)
    capabilities = receptor_backend_capabilities()
    assert set(capabilities) == {"pdb2pqr", "meeko"}
    assert capabilities["pdb2pqr"] is None  # pdb2pqr30 is not on this restricted PATH
    assert capabilities["meeko"] is not None


# ---------------------------------------------------------------- project-local config (R3)


def test_config_creation_writes_nothing_into_working_directory(tmp_path, monkeypatch):
    from autodock_preparation import AutoDockPreparationPipeline, PreparationConfig

    working_dir = tmp_path / "caller_cwd"
    working_dir.mkdir()
    project = tmp_path / "project" / "1-Preparation"
    monkeypatch.chdir(working_dir)
    before = sorted(p.name for p in working_dir.iterdir())

    pipeline = AutoDockPreparationPipeline(
        PreparationConfig(
            ligands_input=str(tmp_path / "ligands_raw"),
            receptors_input=str(tmp_path / "receptors_raw"),
            ligands_output=str(project / "ligands" / "prepared"),
            receptors_output=str(project / "receptors" / "prepared"),
        )
    )
    config_path = Path(pipeline.create_config_file())

    assert sorted(p.name for p in working_dir.iterdir()) == before == []
    assert config_path.is_file()
    assert config_path.resolve().is_relative_to(project.resolve())
    assert config_path.parent.name == "autodock_preparation"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    assert config["output"]["logs"] == "logs"  # relative: logs follow the config into the project
    assert config["preparation"]["receptor_use_pdb2pqr"] is True
    assert config["preparation"]["allow_bad_res"] is False


def test_shell_requires_explicit_config_and_writes_nothing_to_cwd(tmp_path):
    working_dir = tmp_path / "caller_cwd"
    working_dir.mkdir()
    result = subprocess.run(
        ["bash", str(REPO_ROOT / "prep_autodock_enhanced.sh")],
        cwd=working_dir,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "configuration path is required" in result.stderr
    assert list(working_dir.iterdir()) == []


# ---------------------------------------------------------------- shell receptor route (R2b)


@pytest.mark.skipif(
    not (_have("jq") and _have("obabel")),
    reason="needs jq and obabel (the shell wrapper requires both); the restricted PATH is built inside the test",
)
def test_shell_records_backend_unavailable_receptor_failure_and_no_pdbqt(tmp_path):
    from autodock_preparation import AutoDockPreparationPipeline, PreparationConfig

    # Restricted PATH: only jq and obabel (required by the shell) plus standard utilities.
    tool_bin = tmp_path / "tool_bin"
    tool_bin.mkdir()
    for tool in ("jq", "obabel"):
        (tool_bin / tool).symlink_to(shutil.which(tool))
    path = os.pathsep.join([str(tool_bin), "/usr/bin", "/bin"])
    if any(shutil.which(name, path=path) for name in ("pdb2pqr30", "mk_prepare_receptor.py")):
        pytest.skip("PDB2PQR or Meeko is installed in a system folder; restricted PATH cannot be isolated")

    raw = tmp_path / "receptors_raw"
    _write(raw / "rec_chainA.pdb", SOURCE_PDB)
    ligands_raw = tmp_path / "ligands_raw"
    ligands_raw.mkdir()
    receptors_out = tmp_path / "project" / "receptors" / "prepared"
    ligands_out = tmp_path / "project" / "ligands" / "prepared"

    pipeline = AutoDockPreparationPipeline(
        PreparationConfig(
            ligands_input=str(ligands_raw),
            receptors_input=str(raw),
            ligands_output=str(ligands_out),
            receptors_output=str(receptors_out),
        )
    )
    config_path = pipeline.create_config_file()
    env = dict(os.environ, PATH=path, PDBWIZARD_PYTHON=sys.executable)
    result = subprocess.run(
        ["bash", str(REPO_ROOT / "prep_autodock_enhanced.sh"), config_path],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        env=env,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    failures = (receptors_out / "receptor_preparation_failures.csv").read_text(encoding="utf-8").splitlines()
    assert len(failures) == 2, failures  # header + one failed receptor
    assert '"failed","backend_unavailable"' in failures[1]
    assert list(receptors_out.glob("*.pdbqt")) == []
    assert "receptor failure(s)" in result.stdout + result.stderr
    # Logs stay inside the project, next to the configuration. The working directory stays empty.
    assert list(Path(config_path).parent.glob("logs/autodock_prep_*.log"))
    assert sorted(p.name for p in tmp_path.iterdir()) == ["ligands_raw", "project", "receptors_raw", "tool_bin"]


# ---------------------------------------------------------------- 1IEP chain A (dependency-backed)


def _chain_a_protein(tmp_path: Path) -> Path:
    """Protein-only chain A of the 1IEP fixture: ATOM records only, no HETATM, waters or altlocs."""
    atoms = [
        line
        for line in FIXTURE_1IEP.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM  ") and line[21:22] == "A"
    ]
    assert len(atoms) == 2229
    return _write(tmp_path / "1IEP_chainA_protein.pdb", "\n".join(atoms + ["TER", "END"]) + "\n")


@REQUIRES_1IEP_TOOLS
def test_1iep_chain_a_pdb2pqr_default_is_refused_on_heavy_atom_change(tmp_path):
    # Observed on pdb2pqr 3.4.1 (AMBER, pH 7.4): PDB2PQR adds a terminal carboxylate oxygen at
    # GLN 498 and moves His 295/375/396, Asn 414, Gln 252 and Thr 272 side-chain atoms. The strict
    # conservation gate must fail closed with a precise reason, not weaken the check.
    source = _chain_a_protein(tmp_path)
    destination = tmp_path / "out" / "1IEP_chainA_receptor.pdbqt"
    with pytest.raises(ReceptorPreparationError) as excinfo:
        prepare_receptor(source, destination, _config())
    assert excinfo.value.reason == "heavy_atom_conservation_failed"
    assert "GLN 498" in excinfo.value.details
    assert not destination.exists()


@REQUIRES_1IEP_TOOLS
@pytest.mark.xfail(
    strict=True,
    raises=ReceptorPreparationError,
    reason=(
        "Spec 033 acceptance 3 is blocked by a Scientific Lead decision: PDB2PQR changes heavy atoms "
        "of 1IEP chain A (C-terminal O at GLN 498; His/Asn/Gln/Thr side-chain flips). The strict "
        "conservation gate refuses it. See test_1iep_chain_a_pdb2pqr_default_is_refused_on_heavy_atom_change."
    ),
)
def test_1iep_chain_a_pdb2pqr_receptor_is_prepared_at_ph_7_4(tmp_path):
    source = _chain_a_protein(tmp_path)
    destination = tmp_path / "out" / "1IEP_chainA_receptor.pdbqt"
    payload = prepare_receptor(source, destination, _config())
    assert payload["protonation_status"] == "pdb2pqr_applied"
    assert payload["requested_ph"] == 7.4
    assert payload["hd_count"] > 0
    assert payload["heavy_atom_conservation"] == "passed"


@REQUIRES_1IEP_TOOLS
def test_1iep_chain_a_template_receptor_has_polar_hydrogens_and_provenance(tmp_path):
    source = _chain_a_protein(tmp_path)
    destination = tmp_path / "out" / "1IEP_chainA_receptor.pdbqt"
    payload = prepare_receptor(source, destination, _config(receptor_use_pdb2pqr=False))
    assert payload["protonation_status"] == "template_selected_not_ph_titrated"
    assert payload["heavy_atom_count"] == 2229
    assert payload["heavy_atom_conservation"] == "passed"
    assert payload["hd_count"] > 0
    assert payload["hd_near_nitrogen_count"] > 0
    assert payload["input_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert payload["output_sha256"] == hashlib.sha256(destination.read_bytes()).hexdigest()
    assert payload["tool_versions"]["pdb2pqr"] is not None
    assert payload["tool_versions"]["meeko"] is not None
