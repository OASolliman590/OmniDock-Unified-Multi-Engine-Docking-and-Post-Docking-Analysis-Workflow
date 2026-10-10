"""Spec 034 R2/R3: ligand protonation policy, measured provenance, and explicit receptor pH entry.

Dependency-backed tests use the Phase 2 STI (imatinib, 1IEP A201) fixture with the real Open Babel,
RDKit and Meeko tools and skip with a reason when a tool is not installed. Pure-logic tests run
everywhere.
"""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
STI_SDF = REPO_ROOT / "test" / "fixtures" / "spec032" / "phase2_1iep_inputs" / "STI_1IEP_A201.sdf"
FIXTURE_1IEP = REPO_ROOT / "test" / "fixtures" / "spec032" / "phase2_1iep_inputs" / "1IEP.pdb"
# Imatinib with the N-methylpiperazine nitrogen protonated (+1 microspecies).
STI_PLUS1_SMILES = "Cc1ccc(NC(=O)c2ccc(CN3CC[NH+](C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"
STI_NEUTRAL_SMILES = "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"

HAS_OBABEL = shutil.which("obabel") is not None
HAS_MEEKO = shutil.which("mk_prepare_ligand.py") is not None
HAS_PDB2PQR = shutil.which("pdb2pqr30") is not None and shutil.which("mk_prepare_receptor.py") is not None
try:  # RDKit is optional for the module import but required for every measurement.
    import rdkit  # noqa: F401

    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

REQUIRES_OBABEL = pytest.mark.skipif(not (HAS_OBABEL and HAS_RDKIT), reason="needs Open Babel and RDKit")
REQUIRES_MEEKO = pytest.mark.skipif(not (HAS_OBABEL and HAS_RDKIT and HAS_MEEKO), reason="needs Open Babel, RDKit and Meeko (mk_prepare_ligand.py)")
REQUIRES_RDKIT = pytest.mark.skipif(not HAS_RDKIT, reason="needs RDKit")
REQUIRES_PDB2PQR = pytest.mark.skipif(not HAS_PDB2PQR, reason="needs PDB2PQR (pdb2pqr30) and Meeko (mk_prepare_receptor.py)")

PROTONATION_ENV = (
    "PDBWIZARD_LIGAND_PREP_PH",
    "PDBWIZARD_LIGAND_PREP_PH_SOURCE",
    "PDBWIZARD_LIGAND_PREP_PROTONATION_POLICY",
    "PDBWIZARD_LIGAND_PREP_STATE_MAP",
    "PDBWIZARD_LIGAND_PREP_PROFILE",
    "PDBWIZARD_LIGAND_PREP_BACKEND",
    "PDBWIZARD_SELECTED_ENGINES",
    "PDBWIZARD_LIGAND_PREP_REPORT_DIR",
)


@pytest.fixture
def clean_env(monkeypatch):
    for name in PROTONATION_ENV:
        monkeypatch.delenv(name, raising=False)
    return monkeypatch


def _write_state_map(path: Path, rows: str) -> Path:
    path.write_text("ligand,smiles,net_charge\n" + rows, encoding="utf-8")
    return path


def _measure(path: Path) -> dict:
    from docking.preparation.ligand_preparation import measure_ligand_file

    return measure_ligand_file(path)


# ---------------------------------------------------------------- R2a/R2b: the L2 root cause


@REQUIRES_OBABEL
def test_pre_fix_sequence_leaves_sti_neutral(tmp_path):
    """Regression pin on Open Babel: gen3d then '-h -p 7.4' (the old order) leaves STI neutral."""
    gen3d = tmp_path / "gen3d.sdf"
    old = tmp_path / "old_hp.sdf"
    subprocess.run(["obabel", str(STI_SDF), "-O", str(gen3d), "--gen3d"], check=True, capture_output=True)
    subprocess.run(["obabel", str(gen3d), "-O", str(old), "-h", "-p", "7.4"], check=True, capture_output=True)
    state = _measure(old)
    assert state["net_formal_charge"] == 0
    assert state["charged_atoms"] == []
    assert state["hydrogen_count"] > 0  # hydrogens were added, but the pH model did nothing


@REQUIRES_OBABEL
def test_fixed_ph_model_path_measures_a_charged_state(clean_env, tmp_path):
    from docking.preparation import ligand_preparation as lp

    commands: list[list[str]] = []
    original_run = lp._run_command

    def _recording_run(command, env=None, cwd=None):
        commands.append(list(command))
        return original_run(command, env=env, cwd=cwd)

    clean_env.setattr(lp, "_run_command", _recording_run)
    output = tmp_path / "sti_ph_model.sdf"
    normalization = lp._normalize_with_openbabel(STI_SDF, output, protonation_ph=7.4)

    ph_commands = [command for command in commands if "-p" in command]
    assert len(ph_commands) == 1, "exactly one pH-model step"
    assert "-d" in ph_commands[0] and "-h" not in ph_commands[0], "pH model runs on the hydrogen-free molecule"
    assert not any("-h" in command and "-p" in command for command in commands), "no pH-and-H combined step"

    state = normalization["output_state"]
    assert normalization["protonation_step_ran"] is True
    assert normalization["protonation_applied"] is True
    assert state["net_formal_charge"] == 2
    charged = sorted((atom["element"], atom["charge"]) for atom in state["charged_atoms"])
    assert charged == [("N", 1), ("N", 1)]
    assert state["hydrogen_count"] > 0


@REQUIRES_MEEKO
def test_ph_model_docking_input_keeps_the_measured_charge(clean_env, tmp_path):
    from docking.preparation import ligand_preparation as lp

    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PROFILE", "openbabel_meeko")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PH", "7.4")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PH_SOURCE", "user_entered")
    destination = tmp_path / "sti_ph_model.pdbqt"
    summary = lp.prepare_ligand_for_vina_family(STI_SDF, destination)

    protonation = summary["protonation"]
    assert summary["protonation_applied"] is True
    assert protonation["human_review_required"] is True
    assert protonation["ph_source"] == "user_entered"
    assert protonation["prepared_state"]["net_formal_charge"] == 2
    assert protonation["pdbqt_state"]["net_formal_charge"] == 2
    assert protonation["pdbqt_state_check"] == "passed"
    assert summary["output_contract_validation"]["is_valid"] is True


def test_provenance_never_claims_protonation_without_a_measurement(clean_env, monkeypatch, tmp_path):
    from docking.preparation import ligand_preparation as lp

    if not (HAS_OBABEL and HAS_RDKIT):
        pytest.skip("needs Open Babel and RDKit to run the pH step")

    def _unmeasurable(path):
        raise ValueError("measurement unavailable in this test")

    monkeypatch.setattr(lp, "measure_ligand_file", _unmeasurable)
    output = tmp_path / "unmeasured.sdf"
    normalization = lp._normalize_with_openbabel(STI_SDF, output, protonation_ph=7.4)
    assert normalization["protonation_step_ran"] is True
    assert normalization["output_state"] is None
    assert normalization["protonation_applied"] is False, "no measurement, so no protonation claim"
    assert "measurement unavailable" in normalization["output_state_error"]

    # The contract refuses a summary that claims protonation without a measured prepared state.
    summary = {
        "requested_profile": "openbabel_meeko",
        "effective_profile": "openbabel_meeko",
        "preparation_method": "openbabel_then_meeko",
        "protonation_ph": 7.4,
        "normalization": {"normalization_backend": "openbabel"},
        "output_file": str(tmp_path / "missing.pdbqt"),
        "protonation": {"protonation_applied": True, "ph_model_run": True, "prepared_state": None},
    }
    contract = lp.validate_ligand_preparation_output_contract(summary)
    assert contract["is_valid"] is False
    assert any("without a measured prepared state" in error for error in contract["errors"])


def test_as_input_profile_never_claims_a_pH_step(clean_env, tmp_path):
    from docking.preparation import ligand_preparation as lp

    if not HAS_RDKIT:
        pytest.skip("needs RDKit")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PH", "7.4")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PH_SOURCE", "user_entered")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PROTONATION_POLICY", "as_input")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PROFILE", "openbabel_meeko")
    if not HAS_MEEKO:
        pytest.skip("needs Meeko for the PDBQT step")
    summary = lp.prepare_ligand_for_vina_family(STI_SDF, tmp_path / "as_input.pdbqt")
    assert summary["protonation_applied"] is False
    assert summary["protonation"]["ph_model_run"] is False


# ---------------------------------------------------------------- R2c: explicit_state


@REQUIRES_MEEKO
def test_explicit_state_smiles_prepares_the_requested_microspecies(clean_env, tmp_path):
    from docking.preparation import ligand_preparation as lp

    state_map = _write_state_map(tmp_path / "states.csv", f"STI_1IEP_A201,{STI_PLUS1_SMILES},\n")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PROFILE", "openbabel_meeko")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PROTONATION_POLICY", "explicit_state")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_STATE_MAP", str(state_map))
    destination = tmp_path / "sti_plus1.pdbqt"
    summary = lp.prepare_ligand_for_vina_family(STI_SDF, destination)

    protonation = summary["protonation"]
    assert protonation["state_source"] == "explicit_state_smiles"
    assert protonation["explicit_state"]["smiles"] == STI_PLUS1_SMILES
    assert protonation["prepared_state"]["net_formal_charge"] == 1
    assert protonation["pdbqt_state"]["net_formal_charge"] == 1
    assert protonation["pdbqt_state_check"] == "passed"
    assert protonation["pdbqt_canonical_match"] is True
    # The PDBQT's own REMARK SMILES is the same molecule as the given SMILES, charges included.
    remark = lp.measure_pdbqt_remark_smiles(destination)
    assert remark["microspecies_smiles"] == lp.canonical_smiles_from_text(STI_PLUS1_SMILES)
    assert remark["net_formal_charge"] == 1
    # No pH step ran, so no protonation claim is made for an explicit state.
    assert summary["protonation_applied"] is False
    assert protonation["human_review_required"] is False


@REQUIRES_RDKIT
def test_explicit_state_uses_input_heavy_atom_coordinates(clean_env, tmp_path):
    from docking.preparation import ligand_preparation as lp

    out = tmp_path / "explicit.sdf"
    result = lp.build_explicit_state_sdf(STI_SDF, STI_PLUS1_SMILES, out)
    assert result["coordinates_source"] == "input_heavy_atoms"

    def _heavy_coordinates(path: Path):
        from rdkit import Chem

        mol = Chem.RemoveHs(Chem.SDMolSupplier(str(path), removeHs=False)[0])
        conformer = mol.GetConformer()
        return sorted(tuple(round(v, 3) for v in conformer.GetAtomPosition(i)) for i in range(mol.GetNumAtoms()))

    assert _heavy_coordinates(out) == _heavy_coordinates(STI_SDF)


@REQUIRES_RDKIT
def test_explicit_state_rejects_a_smiles_for_a_different_molecule(clean_env, tmp_path):
    from docking.preparation import ligand_preparation as lp

    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PROFILE", "openbabel_meeko")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PROTONATION_POLICY", "explicit_state")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_STATE_MAP", str(_write_state_map(tmp_path / "bad.csv", "STI_1IEP_A201,CCO,\n")))
    with pytest.raises(lp.ProtonationStateError) as info:
        lp.prepare_ligand_for_vina_family(STI_SDF, tmp_path / "bad.pdbqt")
    assert info.value.reason == "explicit_state_mismatch"
    assert not (tmp_path / "bad.pdbqt").exists()


@REQUIRES_OBABEL
@pytest.mark.parametrize("net_charge", [0, 1])
def test_explicit_state_net_charge_mismatch_fails(clean_env, tmp_path, net_charge):
    from docking.preparation import ligand_preparation as lp

    # The pH model gives +2 for STI, so an approved net charge of 0 or +1 must be refused.
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PROFILE", "openbabel_meeko")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PH", "7.4")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PH_SOURCE", "user_entered")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PROTONATION_POLICY", "explicit_state")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_STATE_MAP", str(_write_state_map(tmp_path / "net.csv", f"STI_1IEP_A201,,{net_charge}\n")))
    destination = tmp_path / "net_mismatch.pdbqt"
    with pytest.raises(lp.ProtonationStateError) as info:
        lp.prepare_ligand_for_vina_family(STI_SDF, destination)
    assert info.value.reason == "explicit_state_net_charge_mismatch"
    assert not destination.exists()


@REQUIRES_OBABEL
def test_explicit_state_net_charge_match_prepares(clean_env, tmp_path):
    from docking.preparation import ligand_preparation as lp

    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PROFILE", "openbabel_only")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PH", "7.4")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PH_SOURCE", "user_entered")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PROTONATION_POLICY", "explicit_state")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_STATE_MAP", str(_write_state_map(tmp_path / "net2.csv", "STI_1IEP_A201,,2\n")))
    summary = lp.prepare_ligand_for_vina_family(STI_SDF, tmp_path / "net2.pdbqt")
    assert summary["protonation"]["prepared_state"]["net_formal_charge"] == 2


def test_state_map_rows_must_be_consistent(tmp_path):
    from docking.preparation import ligand_preparation as lp

    if not HAS_RDKIT:
        pytest.skip("needs RDKit")
    bad = _write_state_map(tmp_path / "disagree.csv", f"STI,{STI_PLUS1_SMILES},0\n")
    with pytest.raises(lp.ProtonationStateError) as info:
        lp.load_protonation_state_map(bad)
    assert info.value.reason == "state_map_invalid"


# ---------------------------------------------------------------- R2c: as_input


@REQUIRES_RDKIT
def test_as_input_keeps_the_input_charges_and_adds_missing_hydrogens(clean_env, tmp_path):
    from docking.preparation import ligand_preparation as lp

    # Input with explicit hydrogens and a +1 charge (built by the explicit path, then re-read).
    plus1 = tmp_path / "plus1_input.sdf"
    lp.build_explicit_state_sdf(STI_SDF, STI_PLUS1_SMILES, plus1)
    kept = lp.prepare_as_input_sdf(plus1, tmp_path / "plus1_kept.sdf")
    assert kept["hydrogens_added"] is False
    assert kept["output_state"]["net_formal_charge"] == 1
    assert kept["protonation_applied"] is False

    added = lp.prepare_as_input_sdf(STI_SDF, tmp_path / "sti_as_input.sdf")
    assert added["hydrogens_added"] is True
    assert added["output_state"]["net_formal_charge"] == 0
    assert added["protonation_step_ran"] is False


@REQUIRES_RDKIT
def test_as_input_policy_is_recorded_in_provenance(clean_env, tmp_path):
    from docking.preparation import ligand_preparation as lp

    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PROTONATION_POLICY", "as_input")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PROFILE", "openbabel_meeko")
    if not HAS_MEEKO:
        pytest.skip("needs Meeko for the PDBQT step")
    summary = lp.prepare_ligand_for_vina_family(STI_SDF, tmp_path / "policy.pdbqt")
    assert summary["protonation_policy"] == "as_input"
    assert summary["protonation"]["state_source"] == "as_input"
    assert summary["protonation"]["human_review_required"] is False


# ---------------------------------------------------------------- R2c: ligand CLI and policy plumbing


def test_ligand_cli_requires_ph_for_the_default_policy(clean_env, capsys, tmp_path):
    from docking.preparation import ligand_preparation as lp

    code = lp.main(["--input", str(STI_SDF), "--output", str(tmp_path / "x.pdbqt")])
    assert code == 2
    assert "--ph is required" in capsys.readouterr().err


def test_ligand_pH_without_source_is_refused(clean_env):
    from docking.preparation import ligand_preparation as lp

    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PH", "7.4")
    with pytest.raises(lp.ProtonationStateError) as info:
        lp._selected_ph_with_source()
    assert info.value.reason == "ph_source_missing"


def test_state_map_rejects_unknown_columns_and_duplicates(tmp_path):
    from docking.preparation import ligand_preparation as lp

    no_ligand = tmp_path / "no_ligand.csv"
    no_ligand.write_text("name,net_charge\nSTI,1\n", encoding="utf-8")
    with pytest.raises(lp.ProtonationStateError) as info:
        lp.load_protonation_state_map(no_ligand)
    assert info.value.reason == "state_map_invalid"

    duplicate = _write_state_map(tmp_path / "dup.csv", "STI,,1\nSTI,,1\n")
    with pytest.raises(lp.ProtonationStateError) as info:
        lp.load_protonation_state_map(duplicate)
    assert "duplicate" in info.value.details


# ---------------------------------------------------------------- R3: explicit receptor pH


def test_prepare_protein_without_ph_exits_nonzero_with_message(capsys, tmp_path):
    from workflow.cli import main

    with pytest.raises(SystemExit) as info:
        main(
            [
                "pdb",
                "prepare-protein",
                "--receptors-input",
                str(tmp_path / "raw"),
                "--receptors-output",
                str(tmp_path / "out"),
                "--force-field",
                "AMBER",
            ]
        )
    assert info.value.code != 0
    assert "--ph" in capsys.readouterr().err


def test_prepare_protein_without_force_field_exits_nonzero(capsys, tmp_path):
    from workflow.cli import main

    with pytest.raises(SystemExit) as info:
        main(
            [
                "pdb",
                "prepare-protein",
                "--receptors-input",
                str(tmp_path / "raw"),
                "--receptors-output",
                str(tmp_path / "out"),
                "--ph",
                "7.4",
            ]
        )
    assert info.value.code != 0
    assert "--force-field" in capsys.readouterr().err


def test_prepare_protein_rejects_a_force_field_outside_pdb2pqr(capsys, tmp_path):
    from workflow.cli import main

    with pytest.raises(SystemExit) as info:
        main(["pdb", "prepare-protein", "--receptors-input", "x", "--receptors-output", "y", "--ph", "7.4", "--force-field", "OPLS"])
    assert info.value.code == 2
    assert "invalid choice" in capsys.readouterr().err


def test_receptor_refuses_missing_ph_source(tmp_path):
    from docking.preparation.receptor_preparation import ReceptorPreparationError, prepare_receptor

    source = tmp_path / "rec.pdb"
    source.write_text("ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n", encoding="utf-8")
    config = {"preparation": {"force_field": "AMBER", "force_field_source": "user_entered", "ph": 7.4}}
    with pytest.raises(ReceptorPreparationError) as info:
        prepare_receptor(source, tmp_path / "out.pdbqt", config)
    assert info.value.reason == "ph_source_missing"


def test_receptor_refuses_missing_ph(tmp_path):
    from docking.preparation.receptor_preparation import ReceptorPreparationError, prepare_receptor

    source = tmp_path / "rec.pdb"
    source.write_text("ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n", encoding="utf-8")
    config = {"preparation": {"force_field": "AMBER", "force_field_source": "user_entered", "ph_source": "user_entered"}}
    with pytest.raises(ReceptorPreparationError) as info:
        prepare_receptor(source, tmp_path / "out.pdbqt", config)
    assert info.value.reason == "ph_required"


def _chain_a_protein(path: Path) -> Path:
    """Protein-only chain A of the 1IEP fixture (PDB2PQR removes HETATM records)."""
    lines = [
        line
        for line in FIXTURE_1IEP.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM  ") and line[21:22] == "A"
    ]
    path.write_text("\n".join(lines + ["END"]) + "\n", encoding="utf-8")
    return path


@REQUIRES_PDB2PQR
@pytest.mark.parametrize("ph_source", ["user_entered", "config_file"])
def test_receptor_provenance_records_ph_and_force_field_source(tmp_path, ph_source):
    from docking.preparation.receptor_preparation import prepare_receptor

    source = _chain_a_protein(tmp_path / "1IEP_chainA.pdb")
    config = {
        "preparation": {
            "force_field": "AMBER",
            "force_field_source": ph_source,
            "ph": 7.4,
            "ph_source": ph_source,
            "allow_bad_res": False,
            "ligand_preparation_profile": "engine_aware_full",
        }
    }
    destination = tmp_path / "out" / "1IEP_chainA.pdbqt"
    payload = prepare_receptor(source, destination, config)
    provenance = json.loads((destination.parent / (destination.name + ".preparation.json")).read_text(encoding="utf-8"))
    for record in (payload, provenance):
        assert record["ph_source"] == ph_source
        assert record["force_field_source"] == ph_source
        assert record["requested_ph"] == 7.4
        assert record["requested_force_field"] == "AMBER"
        assert record["protonation_status"] == "pdb2pqr_applied"


@REQUIRES_OBABEL
def test_pdbqt_charge_change_fails_as_lost_in_conversion(clean_env, monkeypatch, tmp_path):
    """Spec 034 R2b: a docking PDBQT whose REMARK SMILES changes the net charge is refused."""
    from docking.preparation import ligand_preparation as lp

    neutral_sti = STI_NEUTRAL_SMILES

    def _neutralising_meeko(input_sdf, output_pdbqt):
        # Stand-in for a conversion that drops the protonation state.
        output_pdbqt.write_text(f"REMARK SMILES {neutral_sti}\nREMARK SMILES IDX 1 1\n", encoding="utf-8")

    monkeypatch.setattr(lp, "_prepare_with_meeko", _neutralising_meeko)
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PROFILE", "openbabel_meeko")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PH", "7.4")
    clean_env.setenv("PDBWIZARD_LIGAND_PREP_PH_SOURCE", "user_entered")
    destination = tmp_path / "lost.pdbqt"
    with pytest.raises(lp.ProtonationStateError) as info:
        lp.prepare_ligand_for_vina_family(STI_SDF, destination)
    assert info.value.reason == "protonation_state_lost_in_conversion"
    assert not destination.exists(), "a PDBQT that lost its state must not remain as docking input"


# ---------------------------------------------------------------- R3 / R2c: interactive prompts


class _Answer:
    def __init__(self, value):
        self.value = value

    def ask(self):
        return self.value


class _FakeQuestionary:
    """Accepts every default (as pressing Enter does) and records each prompt."""

    def __init__(self):
        self.prompts: list[tuple[str, str, object]] = []

    def Choice(self, title, value=None, checked=False):  # noqa: N802 - mirrors questionary's API
        return value

    def select(self, message, choices=None, default=None):
        self.prompts.append(("select", message, default))
        return _Answer(default)

    def text(self, message, default=""):
        self.prompts.append(("text", message, default))
        return _Answer(default)

    def confirm(self, message, default=True):
        self.prompts.append(("confirm", message, default))
        return _Answer(default)


def test_interactive_prompts_suggest_7_4_and_amber_and_confirm_measured_charges(monkeypatch, tmp_path):
    import types

    import workflow.interactive as interactive

    captured: dict = {}

    def _fake_run(target, *args, **kwargs):
        captured.update(kwargs, target=target)
        states = [
            {
                "ligand": "STI_1IEP_A201",
                "protonation_policy": "ph_model",
                "protonation_applied": True,
                "net_formal_charge": 2,
                "charged_atoms": "N29+1; N32+1",
                "microspecies_smiles": "Cc1ccc(NC(=O)c2ccc(C[NH+]3CC[NH+](C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
                "pdbqt_state_check": "passed",
                "human_review_required": True,
                "is_valid": True,
            }
        ]
        return types.SimpleNamespace(status="completed", outputs={"ligand_protonation_states": states}, notes=[])

    monkeypatch.setattr(interactive, "run_autodock_prepare", _fake_run)
    monkeypatch.setattr(interactive, "load_manifest", lambda root: {"engines": []})
    layout = {
        "raw_proteins": tmp_path / "raw_proteins",
        "raw_ligands": tmp_path / "raw_ligands",
        "prepared_proteins": tmp_path / "prepared_proteins",
        "prepared_ligands": tmp_path / "prepared_ligands",
    }
    monkeypatch.setattr(interactive, "_project_layout", lambda root: layout)
    fake = _FakeQuestionary()
    printed: list[str] = []
    monkeypatch.setattr("builtins.print", lambda *args, **kwargs: printed.append(" ".join(str(a) for a in args)))

    interactive._prepare_assets(fake, tmp_path)

    ph_prompts = [p for p in fake.prompts if p[0] == "text" and p[1].startswith("pH (suggested 7.4")]
    assert ph_prompts and ph_prompts[0][2] == "7.4", "pH prompt suggests 7.4 and needs confirmation"
    force_prompts = [p for p in fake.prompts if p[0] == "select" and p[1] == "PDB2PQR force field"]
    assert force_prompts and force_prompts[0][2] == "AMBER", "force field is chosen from the PDB2PQR list, AMBER suggested"
    assert any(p[0] == "confirm" and p[1].startswith("Accept the measured protonation states") for p in fake.prompts)
    assert captured["ph"] == 7.4 and captured["ph_source"] == "user_entered"
    assert captured["force_field"] == "AMBER" and captured["force_field_source"] == "user_entered"
    assert captured["protonation_policy"] == "ph_model"
    assert any("net charge +2" in line for line in printed), "the measured charge is shown per ligand"
