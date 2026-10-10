"""Spec 037 R1 (layered protonation policy) and R2 (ligand Rg at preparation).

Synthetic project: 2 receptors (1IEP, 2ABC) x 2 ligands.
  - STI: an ``explicit_state`` ligand (a ``ligands.STI`` override, +1 N-methylpiperazinium microspecies).
  - NEUT: a ``ph_model`` ligand under the project default (pH 7.4).
  - 2ABC: a receptor override at pH 6.0 (1IEP keeps the project default 7.4).

Ligand preparation is real (RDKit and Open Babel, about 1 s per ligand). Receptor preparation needs
PDB2PQR, which is not installed here, so the receptor preparation records are written as stubs that use
the same resolver. The receptor module's pH selection is tested through its PDB2PQR command line.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from openpyxl import Workbook

from docking.box_policy import RG_SOURCE_COMPUTED, RG_SOURCE_PREPARATION, compute_ligand_box
from docking.preparation import ligand_preparation as lp
from docking.preparation import pairlist_builder
from docking.preparation import receptor_preparation as rp
from docking.project_layout import (
    LAYOUT_CANONICAL,
    PAIRLIST_COLUMNS,
    ProtonationPolicyConflict,
    ProtonationPolicyError,
    check_preparation_policy_current,
    manifest_path,
    load_protonation_policy,
    resolve_entity_protonation,
    set_protonation_override,
    set_protonation_policy,
)
from test_spec032_md_inputs import _write_protein
from test_spec036_docking_params import pair_patches  # noqa: F401  (monkeypatch fixture for alias files)
from test_spec036_protonation_policy import _write_manifest
from workflow.cli import main

REPO_ROOT = Path(__file__).resolve().parents[1]
STI_SDF = REPO_ROOT / "test" / "fixtures" / "spec032" / "phase2_1iep_inputs" / "STI_1IEP_A201.sdf"
PLUS1_SMILES = "Cc1ccc(NC(=O)c2ccc(CN3CC[NH+](C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"
RECEPTOR_1 = "1IEP.pdbqt"
RECEPTOR_2 = "2ABC.pdbqt"
LIGAND_EXPLICIT = "STI.pdbqt"
LIGAND_PH = "NEUT.pdbqt"

try:
    import rdkit  # noqa: F401

    HAS_RDKIT = True
except ImportError:  # pragma: no cover - environment dependent
    HAS_RDKIT = False

pytestmark = pytest.mark.skipif(not HAS_RDKIT, reason="skipped_missing_dependency: needs RDKit")


# ----------------------------------------------------------------------------- helpers


def _make_project(tmp_path: Path) -> Path:
    """Project with the flat default (ph_model, pH 7.4) and the explicit_state override for STI."""
    project = tmp_path / "project"
    _write_manifest(project, layout_profile="canonical")
    set_protonation_policy(project, receptor_ph=7.4, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.4)
    set_protonation_override(project, kind="ligand", name="STI", ligand_policy="explicit_state", smiles=PLUS1_SMILES)
    return project


def _prepared_dirs(project: Path) -> tuple:
    return project / "prepared_proteins", project / "prepared_ligands"


def _prepare_ligands(project: Path, monkeypatch, names=("STI", "NEUT")) -> dict:
    """Prepare each ligand exactly once with the pipeline's environment, and write its step report."""
    _, prepared = _prepared_dirs(project)
    prepared.mkdir(parents=True, exist_ok=True)
    raw = project / "raw_ligands"
    raw.mkdir(parents=True, exist_ok=True)
    # The pipeline passes the project default policy and pH through the environment (Spec 034 R2c).
    monkeypatch.setenv("PDBWIZARD_PROJECT_DIR", str(project))
    monkeypatch.setenv("PDBWIZARD_LIGAND_PREP_PROTONATION_POLICY", "ph_model")
    monkeypatch.setenv("PDBWIZARD_LIGAND_PREP_PH", "7.4")
    monkeypatch.setenv("PDBWIZARD_LIGAND_PREP_PH_SOURCE", "user_entered")
    monkeypatch.setenv("PDBWIZARD_LIGAND_PREP_PROFILE", "engine_aware_full")
    monkeypatch.setenv("PDBWIZARD_SELECTED_ENGINES", "vina")
    monkeypatch.delenv("PDBWIZARD_LIGAND_PREP_STATE_MAP", raising=False)
    summaries: dict = {}
    for name in names:
        source = raw / f"{name}.sdf"
        shutil.copy(STI_SDF, source)
        summaries[name] = lp.prepare_ligand_for_vina_family(source, prepared / f"{name}.pdbqt")
        lp.write_ligand_preparation_step_report(summaries[name], report_dir=prepared / "preparation_steps")
    return summaries


def _write_receptor_stub(project: Path, receptor_file: str) -> None:
    """Prepared receptor and its preparation record, stamped with the resolved policy (see module docstring)."""
    prepared, _ = _prepared_dirs(project)
    prepared.mkdir(parents=True, exist_ok=True)
    (prepared / receptor_file).write_text("ATOM\n", encoding="utf-8")
    resolution = resolve_entity_protonation(project, "receptor", receptor_file)
    record = {
        "input_file": str(project / "raw_receptors" / f"{Path(receptor_file).stem}.pdb"),
        "output_file": str(prepared / receptor_file),
        "policy_level": resolution["policy_level"],
        "policy_hash": resolution["policy_hash"],
        "policy_hash_version": resolution["policy_hash_version"],
        "requested_ph": resolution["receptor_ph"],
        "requested_force_field": resolution["receptor_force_field"],
        "is_valid": True,
    }
    (prepared / f"{receptor_file}.preparation.json").write_text(json.dumps(record), encoding="utf-8")


def _write_site_excel(path: Path) -> Path:
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "Summary"
    sheet.append(["PDB_ID", "Property", "Value"])
    for pdb_id, centre in (("1IEP", (0.0, 0.0, 0.0)), ("2ABC", (1.0, 1.0, 1.0))):
        for prop, value in (
            ("active_site_center_x", centre[0]),
            ("active_site_center_y", centre[1]),
            ("active_site_center_z", centre[2]),
            ("binding_site_center_method", "binding_site_center_v1"),
        ):
            sheet.append([pdb_id, prop, value])
    workbook.save(path)
    return path


def _build_pairlist(project: Path, tmp_path: Path):
    excel = _write_site_excel(tmp_path / "sites.xlsx")
    prepared_proteins, prepared_ligands = _prepared_dirs(project)
    summary = pairlist_builder.build_pairlists(
        project_root=project,
        prepared_proteins=prepared_proteins,
        prepared_ligands=prepared_ligands,
        excel_path=excel,
        mode="curated_per_protein",
        curated_mapping={
            RECEPTOR_1: [LIGAND_EXPLICIT, LIGAND_PH],
            RECEPTOR_2: [LIGAND_EXPLICIT, LIGAND_PH],
        },
        layout_profile=LAYOUT_CANONICAL,
    )
    return summary, pd.read_csv(summary["pairlist_file"])


def _pair(frame: pd.DataFrame, receptor: str, ligand: str) -> pd.Series:
    rows = frame[(frame["receptor"] == receptor) & (frame["ligand"] == ligand)]
    assert len(rows) == 1, f"expected one pair {receptor} x {ligand}, found {len(rows)}"
    return rows.iloc[0]


def _independent_rg(path: Path) -> float:
    """Heavy-atom unweighted Rg of a PDBQT, computed here without the code under test."""
    coords = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        symbol = line[77:79].strip() if len(line) >= 79 and line[77:79].strip() else line.split()[-1]
        if symbol.split(".")[0].upper() in {"H", "HD", "HS", "D"}:
            continue
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    points = np.asarray(coords, dtype=float)
    centre = points.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum((points - centre) ** 2, axis=1))))


def _record(project: Path, name: str) -> dict:
    _, prepared = _prepared_dirs(project)
    return json.loads((prepared / "preparation_steps" / f"{name}.json").read_text(encoding="utf-8"))


def _write_record(project: Path, name: str, payload: dict) -> None:
    _, prepared = _prepared_dirs(project)
    (prepared / "preparation_steps" / f"{name}.json").write_text(json.dumps(payload), encoding="utf-8")


@pytest.fixture
def env_project(tmp_path, monkeypatch):
    """Project with both ligands prepared and both receptors recorded under the current policy."""
    project = _make_project(tmp_path)
    summaries = _prepare_ligands(project, monkeypatch)
    _write_receptor_stub(project, RECEPTOR_1)
    _write_receptor_stub(project, RECEPTOR_2)
    return project, summaries


# ----------------------------------------------------------------------------- R1a: resolver levels


def test_resolver_reports_project_receptor_and_ligand_levels(tmp_path):
    project = _make_project(tmp_path)
    set_protonation_override(project, kind="receptor", name="2ABC", receptor_ph=6.0)

    inherited = resolve_entity_protonation(project, "receptor", RECEPTOR_1)
    assert inherited["policy_level"] == "project"
    assert inherited["receptor_ph"] == 7.4
    assert inherited["receptor_force_field"] == "AMBER"
    assert set(inherited["field_levels"].values()) == {"project"}

    overridden = resolve_entity_protonation(project, "receptor", RECEPTOR_2)
    assert overridden["policy_level"] == "receptor"
    assert overridden["receptor_ph"] == 6.0
    # The force field is not overridden, so it is inherited from the project default.
    assert overridden["field_levels"] == {"receptor_ph": "receptor", "receptor_force_field": "project"}
    assert overridden["receptor_force_field"] == "AMBER"

    explicit = resolve_entity_protonation(project, "ligand", LIGAND_EXPLICIT)
    assert explicit["policy_level"] == "ligand"
    assert explicit["ligand_policy"] == "explicit_state"
    assert explicit["explicit_state"]["smiles"] == PLUS1_SMILES
    assert explicit["ligand_ph"] is None

    default = resolve_entity_protonation(project, "ligand", LIGAND_PH)
    assert default["policy_level"] == "project"
    assert default["ligand_policy"] == "ph_model"
    assert default["ligand_ph"] == 7.4


def test_flat_spec036_block_is_the_project_default_and_overrides_do_not_change_it(tmp_path):
    project = tmp_path / "flat"
    _write_manifest(project, layout_profile="canonical")
    set_protonation_policy(project, receptor_ph=7.4, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.4)
    flat = dict(load_protonation_policy(project))
    resolution = resolve_entity_protonation(project, "ligand", "ANY.pdbqt")
    assert resolution["policy_level"] == "project"
    assert resolution["policy_hash"]
    set_protonation_override(project, kind="ligand", name="ANY", ligand_ph=6.5)
    after = load_protonation_policy(project)
    assert {key: after[key] for key in flat} == flat
    assert resolve_entity_protonation(project, "ligand", "ANY.pdbqt")["policy_level"] == "ligand"


def test_no_stored_policy_resolves_to_none_without_a_hash(tmp_path):
    project = tmp_path / "bare"
    _write_manifest(project, layout_profile="canonical")
    resolution = resolve_entity_protonation(project, "receptor", RECEPTOR_1)
    assert resolution["policy_level"] == "none"
    assert resolution["policy_hash"] is None


def test_override_needs_the_project_default_and_refuses_a_different_value_without_replace(tmp_path):
    bare = tmp_path / "bare"
    _write_manifest(bare, layout_profile="canonical")
    with pytest.raises(ProtonationPolicyError, match="protonation_default_missing"):
        set_protonation_override(bare, kind="receptor", name="1IEP", receptor_ph=6.0)

    project = _make_project(tmp_path)
    set_protonation_override(project, kind="receptor", name="2ABC", receptor_ph=6.0)
    with pytest.raises(ProtonationPolicyConflict, match="protonation_policy_conflict:receptors.2ABC.receptor_ph"):
        set_protonation_override(project, kind="receptor", name="2ABC", receptor_ph=6.5)
    set_protonation_override(project, kind="receptor", name="2ABC", receptor_ph=6.5, replace=True)
    assert resolve_entity_protonation(project, "receptor", RECEPTOR_2)["receptor_ph"] == 6.5


def test_state_map_csv_with_many_ligands_stores_only_the_named_row(tmp_path):
    project = tmp_path / "csv_project"
    _write_manifest(project, layout_profile="canonical")
    set_protonation_policy(project, receptor_ph=7.4, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.4)
    csv_path = tmp_path / "states.csv"
    csv_path.write_text(
        "ligand,smiles,net_charge\n"
        f"STI,{PLUS1_SMILES},\n"
        "NEUT,,0\n",
        encoding="utf-8",
    )
    code = main([
        "workflow", "protonation-policy", "--project-dir", str(project),
        "--ligand", "STI", "--ligand-policy", "explicit_state", "--ligand-state-map", str(csv_path),
    ])
    assert code == 0
    stored = load_protonation_policy(project)["ligands"]
    assert set(stored) == {"STI"}
    assert stored["STI"]["smiles"] == PLUS1_SMILES
    assert stored["STI"]["state_map"]["sha256"]


def test_cli_sets_receptor_override_then_refuses_a_conflicting_value(tmp_path, capsys):
    project = _make_project(tmp_path)
    assert main(["workflow", "protonation-policy", "--project-dir", str(project), "--receptor", "2ABC", "--receptor-ph", "6.0"]) == 0
    assert resolve_entity_protonation(project, "receptor", RECEPTOR_2)["policy_level"] == "receptor"
    code = main(["workflow", "protonation-policy", "--project-dir", str(project), "--receptor", "2ABC", "--receptor-ph", "5.5"])
    captured = capsys.readouterr()
    assert code == 2
    assert "protonation_policy_conflict:receptors.2ABC.receptor_ph" in (captured.out + captured.err)


def test_cli_flat_default_without_its_required_values_is_refused(tmp_path, capsys):
    project = tmp_path / "flat_missing"
    _write_manifest(project, layout_profile="canonical")
    code = main(["workflow", "protonation-policy", "--project-dir", str(project), "--receptor-ph", "7.4"])
    assert code == 2
    assert "--receptor-force-field" in capsys.readouterr().err


# ----------------------------------------------------------------------------- R1b: prepare once, reuse everywhere


def test_each_ligand_is_prepared_once_and_explicit_state_microspecies_is_shared(env_project, pair_patches, tmp_path):
    project, summaries = env_project
    # One preparation per ligand file. Both receptors then reuse the same prepared file.
    assert set(summaries) == {"STI", "NEUT"}
    summary, frame = _build_pairlist(project, tmp_path)
    assert len(frame) == 4

    explicit_pairs = frame[frame["ligand"] == LIGAND_EXPLICIT]
    assert set(explicit_pairs["receptor"]) == {RECEPTOR_1, RECEPTOR_2}
    expected = lp.canonical_smiles_from_text(PLUS1_SMILES)
    record = _record(project, "STI")
    assert record["protonation"]["prepared_state"]["microspecies_smiles"] == expected
    assert record["policy_level"] == "ligand"
    assert record["policy_hash"] == resolve_entity_protonation(project, "ligand", "STI")["policy_hash"]
    # The same prepared file serves both pairs: one sha256, one Rg.
    assert explicit_pairs["ligand_rg_angstrom"].nunique() == 1
    assert set(explicit_pairs["pair_protonation_status"]) == {"explicit_state_no_ph"}
    assert explicit_pairs["ligand_ph"].isna().all()


def test_receptor_override_ph_mismatch_is_flagged_only_on_that_receptors_ph_model_pairs(env_project, pair_patches, tmp_path):
    project, _ = env_project
    set_protonation_override(project, kind="receptor", name="2ABC", receptor_ph=6.0)
    _write_receptor_stub(project, RECEPTOR_2)  # the receptor's policy changed, so its record is re-written
    _, frame = _build_pairlist(project, tmp_path)

    coherent = _pair(frame, RECEPTOR_1, LIGAND_PH)
    assert coherent["pair_protonation_status"] == "coherent"
    assert coherent["receptor_policy_level"] == "project" and coherent["ligand_policy_level"] == "project"
    assert float(coherent["receptor_ph"]) == 7.4 and float(coherent["ligand_ph"]) == 7.4

    mismatch = _pair(frame, RECEPTOR_2, LIGAND_PH)
    assert mismatch["pair_protonation_status"] == "pair_ph_mismatch"
    assert mismatch["receptor_policy_level"] == "receptor"
    assert float(mismatch["receptor_ph"]) == 6.0 and float(mismatch["ligand_ph"]) == 7.4

    for receptor in (RECEPTOR_1, RECEPTOR_2):
        explicit = _pair(frame, receptor, LIGAND_EXPLICIT)
        assert explicit["pair_protonation_status"] == "explicit_state_no_ph"
        assert pd.isna(explicit["ligand_ph"])
    assert set(frame["pair_protonation_status"]) == {"coherent", "pair_ph_mismatch", "explicit_state_no_ph"}


def test_pairlist_columns_carry_the_pair_marks_and_manifest_records_the_counts(env_project, pair_patches, tmp_path):
    project, _ = env_project
    for column in ("receptor_policy_level", "ligand_policy_level", "receptor_ph", "ligand_ph", "pair_protonation_status", "rg_source"):
        assert column in PAIRLIST_COLUMNS
    summary, frame = _build_pairlist(project, tmp_path)
    assert list(frame.columns)[: len(PAIRLIST_COLUMNS)] == PAIRLIST_COLUMNS
    manifest = json.loads(manifest_path(project, LAYOUT_CANONICAL).read_text(encoding="utf-8"))
    assert manifest["pair_protonation"]["status_counts"] == summary["pair_protonation_status_counts"]
    assert sum(summary["pair_protonation_status_counts"].values()) == 4


# ----------------------------------------------------------------------------- R1d: a changed policy refuses stale preparations


def test_changed_ligand_default_refuses_the_stale_ligand_preparation(env_project, pair_patches, tmp_path):
    project, _ = env_project
    set_protonation_policy(
        project, receptor_ph=7.4, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.0, replace=True
    )
    with pytest.raises(ProtonationPolicyError) as refused:
        _build_pairlist(project, tmp_path)
    assert refused.value.reason == "stale_preparation_policy_changed"
    assert "ligand:NEUT" in refused.value.details


def test_changed_receptor_default_refuses_the_stale_receptor_record(env_project, pair_patches, tmp_path):
    project, _ = env_project
    set_protonation_policy(
        project, receptor_ph=6.5, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.4, replace=True
    )
    with pytest.raises(ProtonationPolicyError) as refused:
        _build_pairlist(project, tmp_path)
    assert refused.value.reason == "stale_preparation_policy_changed"
    assert "receptor:1IEP" in refused.value.details


def test_a_preparation_record_without_a_policy_hash_is_refused_under_a_stored_policy(env_project, pair_patches, tmp_path):
    project, _ = env_project
    legacy = _record(project, "NEUT")
    legacy.pop("policy_hash")
    _write_record(project, "NEUT", legacy)
    with pytest.raises(ProtonationPolicyError) as refused:
        _build_pairlist(project, tmp_path)
    assert refused.value.reason == "preparation_policy_unrecorded"


def test_check_is_silent_for_an_entity_whose_record_matches_and_for_no_record(tmp_path):
    project = _make_project(tmp_path)
    resolution = resolve_entity_protonation(project, "ligand", LIGAND_PH)
    check_preparation_policy_current({"policy_hash": resolution["policy_hash"]}, resolution)
    check_preparation_policy_current(None, resolution)
    with pytest.raises(ProtonationPolicyError, match="stale_preparation_policy_changed"):
        check_preparation_policy_current({"policy_hash": "0" * 64}, resolution)


# ----------------------------------------------------------------------------- R1: receptor module uses the override


def _capture_pdb2pqr(monkeypatch) -> list:
    """Stop the receptor run at the PDB2PQR command and record that command line."""
    captured: list = []
    monkeypatch.setattr(rp.shutil, "which", lambda name: f"/fake/bin/{name}")

    def fake_run(command, capture_output=False, text=False, env=None, **_kwargs):
        if str(command[0]).endswith(rp.PDB2PQR_EXECUTABLE):
            captured.append(list(command))
        return subprocess.CompletedProcess(command, 1, "", "stopped before PDB2PQR in the test")

    monkeypatch.setattr(rp.subprocess, "run", fake_run)
    return captured


def _receptor_config(project: Path, ph: float) -> dict:
    return {
        "preparation": {
            "ph": ph,
            "ph_source": "user_entered",
            "force_field": "AMBER",
            "force_field_source": "user_entered",
            "project_dir": str(project),
        }
    }


def test_receptor_module_gives_pdb2pqr_the_override_pH(tmp_path, monkeypatch):
    project = _make_project(tmp_path)
    set_protonation_override(project, kind="receptor", name="2ABC", receptor_ph=6.0)
    source = tmp_path / "raw" / "2ABC.pdb"
    _write_protein(source)
    captured = _capture_pdb2pqr(monkeypatch)
    # The config carries the project default (7.4); the override (6.0) for 2ABC must win.
    with pytest.raises(rp.ReceptorPreparationError):
        rp.prepare_receptor(source, tmp_path / "out" / "2ABC.pdbqt", _receptor_config(project, 7.4))
    assert captured, "PDB2PQR was not reached"
    assert captured[0][captured[0].index("--with-ph") + 1] == "6.0"


def test_receptor_module_refuses_a_config_pH_that_differs_from_the_project_default(tmp_path, monkeypatch):
    project = _make_project(tmp_path)
    source = tmp_path / "raw" / "1IEP.pdb"
    _write_protein(source)
    captured = _capture_pdb2pqr(monkeypatch)
    with pytest.raises(rp.ReceptorPreparationError) as refused:
        rp.prepare_receptor(source, tmp_path / "out" / "1IEP.pdbqt", _receptor_config(project, 7.0))
    assert refused.value.reason == "protonation_policy_conflict"
    assert captured == []


# ----------------------------------------------------------------------------- R2: Rg at ligand preparation


def test_ligand_provenance_records_rg_method_conformer_and_the_unverified_note(env_project):
    project, _ = env_project
    _, prepared = _prepared_dirs(project)
    record = _record(project, "NEUT")
    assert record["rg_method"] == "heavy_atom_unweighted_v1"
    assert record["rg_conformer"] == "prepared"
    assert record["radius_of_gyration_angstrom"] == pytest.approx(_independent_rg(prepared / "NEUT.pdbqt"), abs=1e-6)
    assert record["output_sha256"]
    assert "unverified assumption" in record["rg_note"] and "Spec 037 R2c" in record["rg_note"]
    assert record["rg_definition_status"] == "unverified_assumption_spec037_r2c"


def test_box_reads_rg_from_ligand_preparation_and_says_so(env_project):
    project, _ = env_project
    _, prepared = _prepared_dirs(project)
    box = compute_ligand_box(prepared / "STI.pdbqt")
    record = _record(project, "STI")
    assert box.rg_source == RG_SOURCE_PREPARATION == "ligand_preparation"
    assert box.ligand_rg_angstrom == pytest.approx(record["radius_of_gyration_angstrom"], abs=1e-12)
    assert box.edge_angstrom == pytest.approx(2.9 * record["radius_of_gyration_angstrom"])


def test_box_uses_the_recorded_rg_not_a_recomputation(env_project):
    project, _ = env_project
    _, prepared = _prepared_dirs(project)
    record = _record(project, "STI")
    record["radius_of_gyration_angstrom"] = 5.5  # the file hash is unchanged, so the record is trusted
    _write_record(project, "STI", record)
    box = compute_ligand_box(prepared / "STI.pdbqt")
    assert box.rg_source == "ligand_preparation"
    assert box.ligand_rg_angstrom == pytest.approx(5.5)


def test_box_falls_back_to_computing_rg_with_a_labelled_note(env_project):
    project, _ = env_project
    _, prepared = _prepared_dirs(project)
    shutil.rmtree(prepared / "preparation_steps")
    box = compute_ligand_box(prepared / "STI.pdbqt")
    assert box.rg_source == RG_SOURCE_COMPUTED == "computed_at_box_step"
    assert box.ligand_rg_angstrom == pytest.approx(_independent_rg(prepared / "STI.pdbqt"), abs=1e-6)
    assert any("Rg computed at the box step" in note for note in box.warnings)


def test_changed_conformer_is_not_trusted_for_rg(env_project):
    project, _ = env_project
    _, prepared = _prepared_dirs(project)
    with (prepared / "NEUT.pdbqt").open("a", encoding="utf-8") as handle:
        handle.write("REMARK edited after preparation\n")
    box = compute_ligand_box(prepared / "NEUT.pdbqt")
    assert box.rg_source == RG_SOURCE_COMPUTED
    assert any("output_sha256 differs" in note for note in box.warnings)


def test_pairlist_rows_record_where_the_box_rg_came_from(env_project, pair_patches, tmp_path):
    project, _ = env_project
    _, frame = _build_pairlist(project, tmp_path)
    explicit = _pair(frame, RECEPTOR_1, LIGAND_EXPLICIT)
    assert explicit["rg_source"] == "ligand_preparation"
    assert float(explicit["ligand_rg_angstrom"]) == pytest.approx(_record(project, "STI")["radius_of_gyration_angstrom"], abs=1e-9)


def test_state_given_without_an_explicit_policy_is_refused_not_silently_ignored(tmp_path):
    project = _make_project(tmp_path)  # project default is ph_model, so a state under it would be ignored
    with pytest.raises(ProtonationPolicyError, match="explicit_state_fields_need_explicit_policy"):
        set_protonation_override(project, kind="ligand", name="NEUT", smiles=PLUS1_SMILES)
    assert "NEUT" not in (load_protonation_policy(project).get("ligands") or {})
