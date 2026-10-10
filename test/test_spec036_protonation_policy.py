"""Spec 036 R1a/R1b and the R6 md-input map generation.

- R1a: one project protonation policy (manifest block), explicit values checked against it.
- R1b: the MD export uses the docked microspecies. Vina/Smina lineage on the committed fixture
  ``test/fixtures/spec034/sti_vina_model1.pdbqt`` (neutral) and a tmp copy whose REMARK SMILES is the
  +1 N-methylpiperazinium microspecies. The committed fixture is never modified.
- R6 subset: ``resolve_md_input_maps`` generates tags, receptor and topology maps from the project.

The MD export tests need Open Babel and RDKit and skip with a reason otherwise. Policy, CLI and map
tests are dependency-light.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd
import pytest

from docking.project_layout import (
    ProtonationPolicyConflict,
    check_protonation_policy_conflicts,
    load_protonation_policy,
    protonation_state_map_path,
    resolve_protonation_inputs,
    set_protonation_policy,
)
from post_docking_analysis.md_chemistry import (
    OpenBabelChemistryBackend,
    openbabel_capability,
)
from post_docking_analysis.md_inputs import (
    MDInputsRequest,
    MDMapError,
    check_request_protonation_against_project,
    export_md_inputs,
    resolve_md_input_maps,
)
from test_spec032_md_inputs import FakeChemistryBackend, _write_protein

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_PDBQT = REPO_ROOT / "test" / "fixtures" / "spec034" / "sti_vina_model1.pdbqt"
CRYSTAL_SDF = REPO_ROOT / "test" / "fixtures" / "spec032" / "phase2_1iep_inputs" / "STI_1IEP_A201.sdf"
NEUTRAL_SMILES = "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"
# Imatinib with the distal N-methylpiperazine nitrogen protonated (+1 microspecies).
PLUS1_SMILES = "Cc1ccc(NC(=O)c2ccc(CN3CC[NH+](C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"
NATIVE_SCORE = -17.895  # REMARK INTER + INTRA of the committed MODEL 1

try:
    import rdkit  # noqa: F401

    HAS_RDKIT = True
except ImportError:  # pragma: no cover - environment dependent
    HAS_RDKIT = False

REQUIRES_OBABEL = pytest.mark.skipif(
    openbabel_capability()["status"] != "completed" or not HAS_RDKIT,
    reason="skipped_missing_dependency: needs Open Babel and RDKit",
)
REQUIRES_RDKIT = pytest.mark.skipif(not HAS_RDKIT, reason="skipped_missing_dependency: needs RDKit")


def _canonical(smiles: str) -> str:
    from docking.preparation.ligand_preparation import canonical_smiles_from_text

    return canonical_smiles_from_text(smiles)


def _pose_copy(directory: Path, *, smiles: str | None = None, name: str = "sti_vina_model1.pdbqt") -> Path:
    """Copy of the committed fixture with a native score line (absent in the fixture) and optional SMILES swap."""
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / name
    lines = []
    for line in FIXTURE_PDBQT.read_text(encoding="utf-8").splitlines():
        if smiles is not None and line.startswith("REMARK SMILES ") and not line.startswith("REMARK SMILES IDX"):
            line = f"REMARK SMILES {smiles}"
        lines.append(line)
        if line.startswith("MODEL 1"):
            # The committed MODEL 1 has no REMARK VINA RESULT; the native-score gate needs one.
            lines.append(f"REMARK VINA RESULT: {NATIVE_SCORE:.3f}      0.000      0.000")
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def _write_manifest(project: Path, layout_profile: str = "canonical") -> None:
    project.mkdir(parents=True, exist_ok=True)
    (project / "project_manifest.json").write_text(
        json.dumps(
            {
                "project_name": "spec036",
                "project_root": str(project),
                "layout_profile": layout_profile,
                "engines": ["vina"],
                "favorite_engine": "vina",
            }
        ),
        encoding="utf-8",
    )


def _prepared_dirs_for(project: Path) -> tuple:
    """Prepared receptor and ligand directories of the project, as its recorded layout profile puts them."""
    profile = json.loads((project / "project_manifest.json").read_text(encoding="utf-8"))["layout_profile"]
    if profile == "docking_legacy":
        return project / "3-Preparation" / "2-Prepared_Protiens", project / "3-Preparation" / "4-Prepared_Ligand"
    return project / "prepared_proteins", project / "prepared_ligands"


def _make_docked_project(
    tmp_path: Path,
    *,
    pose_smiles: str | None = None,
    policy: str = "ph_model",
    prepared_smiles: str = NEUTRAL_SMILES,
    ph: float | None = 7.4,
    write_provenance: bool = True,
):
    """A Vina project with one selected pose, a receptor map and an explicit topology map."""
    project = tmp_path / "project"
    _write_manifest(project)
    pose = _pose_copy(tmp_path / "poses", smiles=pose_smiles)
    selection = project / "4-Working" / "scores" / "unified"
    selection.mkdir(parents=True)
    pd.DataFrame(
        [
            {
                "engine": "vina",
                "tag": "prot_site_sti",
                "protein": "prot",
                "ligand": "sti",
                "pose": 1,
                "pose_file": str(pose),
                "affinity_kcal_mol": NATIVE_SCORE,
            }
        ]
    ).to_csv(selection / "best_pose_per_tag_by_engine.csv", index=False)
    tags = project / "tags.csv"
    pd.DataFrame([{"tag": "prot_site_sti"}]).to_csv(tags, index=False)
    receptor = project / "prepared" / "prot.pdb"
    _write_protein(receptor)
    receptor_map = project / "receptor_map.csv"
    pd.DataFrame([{"tag": "prot_site_sti", "receptor_file": str(receptor)}]).to_csv(receptor_map, index=False)
    topology_map = project / "topology_map.csv"
    pd.DataFrame(
        [{"tag": "prot_site_sti", "topology_file": str(CRYSTAL_SDF), "topology_pose": 1, "atom_map": ""}]
    ).to_csv(topology_map, index=False)
    if write_provenance:
        # The docking input carries the pose's REMARK SMILES: the neutral fixture unless a +1 copy is used.
        _write_step_report(
            project,
            policy=policy,
            prepared_smiles=prepared_smiles,
            ph=ph,
            pdbqt_smiles=pose_smiles or NEUTRAL_SMILES,
        )
    return project, pose, tags, receptor_map, topology_map


def _write_step_report(
    project: Path,
    *,
    policy: str,
    prepared_smiles: str,
    ph: float | None,
    pdbqt_smiles: str | None,
    name: str = "sti",
    input_file: Path = CRYSTAL_SDF,
) -> Path:
    """Ligand preparation step report (Spec 034 R2b) for the docked ligand."""
    from rdkit import Chem

    prepared = _canonical(prepared_smiles)
    net = int(Chem.GetFormalCharge(Chem.MolFromSmiles(prepared_smiles)))
    pdbqt_block = None
    if pdbqt_smiles is not None:
        pdbqt_canonical = _canonical(pdbqt_smiles)
        pdbqt_block = {
            "smiles": pdbqt_smiles,
            "net_formal_charge": int(Chem.GetFormalCharge(Chem.MolFromSmiles(pdbqt_smiles))),
            "microspecies_smiles": pdbqt_canonical,
            "stereo_compared": False,
        }
    report_dir = _prepared_dirs_for(project)[1] / "preparation_steps"
    report_dir.mkdir(parents=True, exist_ok=True)
    report = report_dir / f"{name}.json"
    report.write_text(
        json.dumps(
            {
                "input_file": str(input_file),
                "output_file": str(_prepared_dirs_for(project)[1] / f"{name}.pdbqt"),
                "protonation": {
                    "policy": policy,
                    "ph": ph if policy == "ph_model" else None,
                    "ph_model_run": policy == "ph_model",
                    "prepared_state": {
                        "net_formal_charge": net,
                        "microspecies_smiles": prepared,
                        "charged_atoms": [],
                    },
                    "pdbqt_state": pdbqt_block,
                },
            }
        ),
        encoding="utf-8",
    )
    return report


def _export(project, tags, receptor_map, topology_map, *, backend=None, **overrides):
    payload = {
        "project_dir": str(project),
        "engine": "vina",
        "tags_file": str(tags),
        "receptor_map": str(receptor_map),
        "topology_map": str(topology_map),
        "ligand_formats": ["mol2"],
        **overrides,
    }
    request = MDInputsRequest.from_dict(payload)
    return export_md_inputs(request, chemistry_backend=backend or OpenBabelChemistryBackend())


def _provenance(result) -> dict:
    return json.loads(Path(result["rows"][0]["provenance_file"]).read_text(encoding="utf-8"))


# ---------------------------------------------------------------- R1a: policy round trip and conflicts


def _write_state_map(path: Path, smiles: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"ligand,smiles\nsti,{smiles}\n", encoding="utf-8")
    return path


def test_protonation_policy_round_trip_stores_hashed_relative_state_map(tmp_path):
    project = tmp_path / "project"
    _write_manifest(project)
    state_map = _write_state_map(project / "states" / "sti_states.csv", PLUS1_SMILES)

    assert load_protonation_policy(project) is None
    stored = set_protonation_policy(
        project,
        receptor_ph=7.4,
        receptor_force_field="amber",
        ligand_policy="explicit_state",
        ligand_state_map=state_map,
    )
    loaded = load_protonation_policy(project)
    assert loaded == stored
    assert loaded["receptor_ph"] == 7.4
    assert loaded["receptor_force_field"] == "AMBER"
    assert loaded["ligand_policy"] == "explicit_state"
    assert loaded["ligand_ph"] is None
    assert loaded["source"] == "user_entered"
    assert loaded["set_at"].endswith("Z") or "T" in loaded["set_at"]
    assert loaded["ligand_state_map"]["path"] == "states/sti_states.csv"
    assert len(loaded["ligand_state_map"]["sha256"]) == 64
    assert protonation_state_map_path(project, loaded) == state_map.resolve()

    # The manifest keeps every other key it already had.
    manifest = json.loads((project / "project_manifest.json").read_text(encoding="utf-8"))
    assert manifest["engines"] == ["vina"]
    assert manifest["protonation_policy"]["ligand_policy"] == "explicit_state"


def test_stored_state_map_that_changed_fails_closed(tmp_path):
    project = tmp_path / "project"
    _write_manifest(project)
    state_map = _write_state_map(project / "sti_states.csv", PLUS1_SMILES)
    set_protonation_policy(
        project, receptor_ph=7.4, receptor_force_field="AMBER", ligand_policy="explicit_state", ligand_state_map=state_map
    )
    state_map.write_text("ligand,smiles\nsti," + NEUTRAL_SMILES + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="protonation_state_map_hash_mismatch"):
        protonation_state_map_path(project, load_protonation_policy(project))


def test_conflicting_explicit_value_fails_with_both_values_shown(tmp_path):
    project = tmp_path / "project"
    _write_manifest(project)
    set_protonation_policy(project, receptor_ph=7.4, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.4)

    with pytest.raises(ProtonationPolicyConflict) as caught:
        resolve_protonation_inputs(project, receptor_ph=6.0)
    error = caught.value
    assert error.reason == "protonation_policy_conflict"
    assert error.field == "receptor_ph"
    assert "7.4" in str(error) and "6.0" in str(error)
    assert str(error).startswith("protonation_policy_conflict:receptor_ph")

    with pytest.raises(ProtonationPolicyConflict, match="protonation_policy_conflict:receptor_force_field"):
        resolve_protonation_inputs(project, receptor_force_field="CHARMM")
    with pytest.raises(ProtonationPolicyConflict, match="protonation_policy_conflict:ligand_policy"):
        resolve_protonation_inputs(project, ligand_policy="as_input")


def test_matching_explicit_values_pass_and_unset_values_inherit(tmp_path):
    project = tmp_path / "project"
    _write_manifest(project)
    set_protonation_policy(project, receptor_ph=7.4, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.0)

    matched = resolve_protonation_inputs(project, receptor_ph=7.4, receptor_force_field="amber", ligand_policy="PH_MODEL")
    assert matched["receptor_ph"] == 7.4
    assert matched["receptor_ph_source"] == "user_entered"
    assert matched["receptor_force_field"] == "amber"

    inherited = resolve_protonation_inputs(project)
    assert inherited["receptor_ph"] == 7.4
    assert inherited["receptor_ph_source"] == "user_entered"
    assert inherited["receptor_force_field"] == "AMBER"
    assert inherited["ligand_policy"] == "ph_model"
    assert inherited["ligand_ph"] == 7.0
    assert check_protonation_policy_conflicts(load_protonation_policy(project), {"ligand_ph": 7.0}) == []


def test_without_project_policy_explicit_values_are_passed_through(tmp_path):
    project = tmp_path / "project"
    _write_manifest(project)
    resolved = resolve_protonation_inputs(project, receptor_ph=6.5)
    assert resolved["project_policy"] is False
    assert resolved["receptor_ph"] == 6.5
    assert resolved["ligand_policy"] is None


def test_set_policy_refuses_a_different_stored_value_without_replace(tmp_path):
    project = tmp_path / "project"
    _write_manifest(project)
    set_protonation_policy(project, receptor_ph=7.4, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.4)
    with pytest.raises(ProtonationPolicyConflict):
        set_protonation_policy(project, receptor_ph=6.0, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.4)
    replaced = set_protonation_policy(
        project, receptor_ph=6.0, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.4, replace=True
    )
    assert replaced["receptor_ph"] == 6.0


def test_cli_protonation_policy_command_sets_then_enforces(tmp_path, capsys):
    from workflow.cli import main

    project = tmp_path / "project"
    _write_manifest(project)
    code = main(
        [
            "workflow",
            "protonation-policy",
            "--project-dir",
            str(project),
            "--receptor-ph",
            "7.4",
            "--receptor-force-field",
            "AMBER",
            "--ligand-policy",
            "ph_model",
            "--ligand-ph",
            "7.4",
        ]
    )
    assert code == 0
    assert load_protonation_policy(project)["receptor_ph"] == 7.4

    code = main(
        [
            "workflow",
            "protonation-policy",
            "--project-dir",
            str(project),
            "--receptor-ph",
            "6.0",
            "--receptor-force-field",
            "AMBER",
            "--ligand-policy",
            "ph_model",
            "--ligand-ph",
            "7.4",
        ]
    )
    captured = capsys.readouterr()
    assert code != 0
    assert "protonation_policy_conflict:receptor_ph" in (captured.out + captured.err)
    assert load_protonation_policy(project)["receptor_ph"] == 7.4


def test_cli_prepare_protein_conflict_is_refused_before_any_work(tmp_path, capsys):
    from workflow.cli import main

    project = tmp_path / "project"
    _write_manifest(project)
    set_protonation_policy(project, receptor_ph=7.4, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.4)
    with pytest.raises(SystemExit) as caught:
        main(
            [
                "pdb",
                "prepare-protein",
                "--project-dir",
                str(project),
                "--receptors-input",
                str(tmp_path / "raw"),
                "--receptors-output",
                str(tmp_path / "out"),
                "--ph",
                "6.0",
                "--force-field",
                "AMBER",
            ]
        )
    assert caught.value.code == 2
    assert "protonation_policy_conflict:receptor_ph" in capsys.readouterr().err


# ---------------------------------------------------------------- R1b: CLI contract for md-inputs


def test_md_inputs_cli_default_is_docked_state_and_ph_alone_is_refused(tmp_path, capsys):
    from workflow.cli import build_parser, main

    parser = build_parser()
    base = ["analyze", "md-inputs", "--project-dir", str(tmp_path), "--engine", "vina"]
    args = parser.parse_args(base)
    assert args.ph is None and args.reprotonate_at_ph is None and args.protonation_policy is None
    assert args.tags_file is None and args.receptor_map is None and args.auto_maps is False

    args = parser.parse_args(base + ["--reprotonate-at-ph", "7.4"])
    assert args.reprotonate_at_ph == 7.4

    args = parser.parse_args(base + ["--protonation-policy", "openbabel_predicted", "--ph", "7.4"])
    assert args.protonation_policy == "openbabel_predicted" and args.ph == 7.4

    # Refused in main() before any work: a bare pH and the legacy flag without a pH.
    with pytest.raises(SystemExit):
        main(base + ["--ph", "7.4"])
    assert "--reprotonate-at-ph" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        main(base + ["--protonation-policy", "openbabel_predicted"])
    with pytest.raises(SystemExit):
        main(base + ["--auto-maps", "--tags-file", "tags.csv"])


def test_request_override_is_recorded_and_legacy_flag_maps_to_it(tmp_path):
    base = {"project_dir": str(tmp_path), "engine": "vina", "tags_file": "", "receptor_map": ""}
    tags = tmp_path / "t.csv"
    tags.write_text("tag\nx\n", encoding="utf-8")
    receptor = tmp_path / "r.csv"
    receptor.write_text("tag,receptor_file\nx,y\n", encoding="utf-8")
    payload = {**base, "tags_file": str(tags), "receptor_map": str(receptor)}

    default = MDInputsRequest.from_dict(payload)
    assert default.protonation_policy == "docked_state"
    assert default.pH is None and not default.reprotonation_override

    override = MDInputsRequest.from_dict({**payload, "reprotonate_at_ph": 6.5})
    assert override.reprotonation_override and override.pH == 6.5
    assert override.protonation_policy == "openbabel_predicted"
    assert override.override_source == "reprotonate_at_ph"

    legacy = MDInputsRequest.from_dict({**payload, "protonation_policy": "openbabel_predicted", "pH": 7.4})
    assert legacy.pH == 7.4 and legacy.override_source == "legacy_openbabel_predicted"

    with pytest.raises(ValueError, match="conflicting_reprotonation_pH"):
        MDInputsRequest.from_dict(
            {**payload, "protonation_policy": "openbabel_predicted", "pH": 7.4, "reprotonate_at_ph": 6.0}
        )
    with pytest.raises(ValueError, match="bare pH is ambiguous"):
        MDInputsRequest.from_dict({**payload, "pH": 7.4})


def test_override_against_explicit_state_project_policy_conflicts(tmp_path):
    project = tmp_path / "project"
    _write_manifest(project)
    state_map = _write_state_map(project / "sti_states.csv", PLUS1_SMILES)
    set_protonation_policy(
        project, receptor_ph=7.4, receptor_force_field="AMBER", ligand_policy="explicit_state", ligand_state_map=state_map
    )
    tags = project / "t.csv"
    tags.write_text("tag\nx\n", encoding="utf-8")
    receptor = project / "r.csv"
    receptor.write_text("tag,receptor_file\nx,y\n", encoding="utf-8")
    request = MDInputsRequest.from_dict(
        {
            "project_dir": str(project),
            "engine": "vina",
            "tags_file": str(tags),
            "receptor_map": str(receptor),
            "reprotonate_at_ph": 7.4,
        }
    )
    with pytest.raises(ProtonationPolicyConflict, match="protonation_policy_conflict:ligand_policy"):
        check_request_protonation_against_project(project, request)


# ---------------------------------------------------------------- R1b: MD export from the docked microspecies


@REQUIRES_OBABEL
def test_docked_neutral_state_with_ph_model_provenance_is_predicted_and_neutral(tmp_path):
    project, _, tags, receptor_map, topology_map = _make_docked_project(tmp_path, policy="ph_model", ph=7.4)
    result = _export(project, tags, receptor_map, topology_map)
    assert result["status"] == "completed", result["rows"][0]
    provenance = _provenance(result)
    g6 = provenance["gates"]["G6"]
    assert provenance["net_charge"] == 0
    assert provenance["charge_authority"] == "predicted"
    assert provenance["scientific_review"] == "human_review_required"
    assert g6["charge_authority"] == "predicted"
    assert g6["docked_state"]["net_formal_charge"] == 0
    assert g6["docked_state"]["source"] == "meeko_remark_smiles"
    assert g6["heavy_coordinate_max_delta"] == pytest.approx(0.0, abs=1e-9)
    assert g6["implicit_hydrogen_count"] == 0
    assert g6["reprotonation_override"] is None
    assert g6["pH_source"] == "preparation_provenance"
    assert g6["pH"] == 7.4
    assert provenance["gates"]["G3"]["atom_map_source"] == "meeko_smiles_idx_lineage_v1"
    assert provenance["gates"]["G3"]["lineage"]["mapping_charge_policy"] == "charge_blind_protonation_state_v1"


@REQUIRES_OBABEL
def test_docked_neutral_state_with_explicit_state_provenance_is_docked_state(tmp_path):
    project, _, tags, receptor_map, topology_map = _make_docked_project(
        tmp_path, policy="explicit_state", ph=None
    )
    result = _export(project, tags, receptor_map, topology_map)
    assert result["status"] == "completed", result["rows"][0]
    provenance = _provenance(result)
    assert provenance["net_charge"] == 0
    assert provenance["charge_authority"] == "docked_state"
    assert provenance["scientific_review"] == "completed"
    assert provenance["gates"]["G6"]["heavy_coordinate_max_delta"] == pytest.approx(0.0, abs=1e-9)


@REQUIRES_OBABEL
def test_docked_plus_one_state_gives_plus_one_md_ligand_with_zero_heavy_delta(tmp_path):
    pose_smiles = PLUS1_SMILES
    project, _, tags, receptor_map, topology_map = _make_docked_project(
        tmp_path, pose_smiles=pose_smiles, policy="explicit_state", prepared_smiles=pose_smiles, ph=None
    )
    result = _export(project, tags, receptor_map, topology_map)
    assert result["status"] == "completed", result["rows"][0]
    provenance = _provenance(result)
    assert provenance["net_charge"] == 1
    assert provenance["charge_authority"] == "docked_state"
    g6 = provenance["gates"]["G6"]
    assert g6["net_charge"] == 1
    assert g6["docked_state"]["net_formal_charge"] == 1
    assert g6["docked_state"]["microspecies_smiles"] == _canonical(PLUS1_SMILES)
    assert g6["heavy_coordinate_max_delta"] == pytest.approx(0.0, abs=1e-9)
    # The MD ligand carries the protonated amine: 32 hydrogens for the +1 form (31 for the neutral form).
    ligand_mol2 = Path(result["rows"][0]["output_dir"]) / "ligand.mol2"
    atoms = ligand_mol2.read_text(encoding="utf-8").split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>BOND")[0]
    hydrogens = [line for line in atoms.splitlines() if len(line.split()) >= 6 and line.split()[5].startswith("H")]
    assert len(hydrogens) == 32


@REQUIRES_OBABEL
def test_docked_state_provenance_mismatch_fails_g6(tmp_path):
    # The pose says +1 (REMARK SMILES) but the preparation record says neutral.
    project, _, tags, receptor_map, topology_map = _make_docked_project(
        tmp_path, pose_smiles=PLUS1_SMILES, policy="ph_model", prepared_smiles=NEUTRAL_SMILES, ph=7.4
    )
    result = _export(project, tags, receptor_map, topology_map)
    row = result["rows"][0]
    assert row["status"] == "failed"
    assert row["gates"]["G6"]["status"] == "failed"
    assert "docked_state_differs_from_preparation_provenance" in row["gates"]["G6"]["reason"]
    assert result["status"] == "failed"


@REQUIRES_OBABEL
def test_md_ligand_charge_that_differs_from_docked_state_fails_g6(tmp_path):
    project, _, tags, receptor_map, topology_map = _make_docked_project(
        tmp_path, pose_smiles=PLUS1_SMILES, policy="explicit_state", prepared_smiles=PLUS1_SMILES, ph=None
    )

    class NeutralisingBackend(FakeChemistryBackend):
        """Returns a neutral ligand whatever the docked state is: the export must refuse it."""

        def prepare_docked_state(self, source_file, **kwargs):
            ligand = self.prepare(source_file, source_pose=1, pH=None, protonate=True)
            return ligand

    result = _export(project, tags, receptor_map, topology_map, backend=NeutralisingBackend())
    row = result["rows"][0]
    assert row["status"] == "failed"
    assert row["gates"]["G6"]["status"] == "failed"
    assert row["gates"]["G6"]["reason"].startswith("md_state_differs_from_docked_state:0!=1")


@REQUIRES_OBABEL
def test_reprotonation_override_is_recorded_and_keeps_predicted_path(tmp_path):
    project, _, tags, receptor_map, topology_map = _make_docked_project(tmp_path, policy="ph_model", ph=7.4)
    result = _export(project, tags, receptor_map, topology_map, reprotonate_at_ph=7.4)
    assert result["status"] == "completed", result["rows"][0]
    provenance = _provenance(result)
    override = provenance["reprotonation_override"]
    assert override["pH"] == 7.4
    assert override["source"] == "reprotonate_at_ph"
    assert override["docked_net_charge"] == 0
    # The Open Babel pH model over-protonates STI (Spec 034 L2: the +2 dication at pH 7.4). The override
    # reproduces that predicted state and records that it differs from the docked microspecies.
    assert override["md_net_charge"] == 2
    assert override["charge_changed_from_docked"] is True
    assert provenance["charge_authority"] == "predicted"
    assert provenance["gates"]["G6"]["pH_source"] == "reprotonate_at_ph"
    assert provenance["gates"]["G6"]["pH"] == 7.4


@REQUIRES_OBABEL
def test_legacy_openbabel_predicted_flag_is_recorded_as_override(tmp_path):
    project, _, tags, receptor_map, topology_map = _make_docked_project(tmp_path, policy="ph_model", ph=7.4)
    result = export_md_inputs(
        MDInputsRequest.from_dict(
            {
                "project_dir": str(project),
                "engine": "vina",
                "tags_file": str(tags),
                "receptor_map": str(receptor_map),
                "topology_map": str(topology_map),
                "protonation_policy": "openbabel_predicted",
                "pH": 7.4,
                "ligand_formats": ["mol2"],
            }
        ),
        chemistry_backend=OpenBabelChemistryBackend(),
    )
    provenance = _provenance(result)
    assert provenance["reprotonation_override"]["source"] == "legacy_openbabel_predicted"
    assert provenance["reprotonation_override"]["pH"] == 7.4
    assert provenance["charge_authority"] == "predicted"


# ---------------------------------------------------------------- R6 subset: generated md-input maps


def _auto_project(tmp_path: Path, layout_profile: str = "canonical"):
    project = tmp_path / "auto_project"
    _write_manifest(project, layout_profile)
    pose = _pose_copy(tmp_path / "poses")
    selection = project / "4-Working" / "scores" / "unified"
    selection.mkdir(parents=True)
    pd.DataFrame(
        [
            {
                "engine": "vina",
                "tag": "prot_site_sti",
                "protein": "prot",
                "ligand": "sti",
                "pose": 1,
                "pose_file": str(pose),
                "affinity_kcal_mol": NATIVE_SCORE,
            }
        ]
    ).to_csv(selection / "best_pose_per_tag_by_engine.csv", index=False)
    receptors = _prepared_dirs_for(project)[0]
    receptors.mkdir(parents=True)
    chain_pdb = project / "0-Input" / "receptors" / "prot_chain_A.pdb"
    _write_protein(chain_pdb)
    (receptors / "prot.pdbqt").write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00    -0.300 N\n", encoding="utf-8")
    (receptors / "prot.pdbqt.preparation.json").write_text(
        json.dumps({"input_file": str(chain_pdb), "is_valid": True}), encoding="utf-8"
    )
    _write_step_report(project, policy="ph_model", prepared_smiles=NEUTRAL_SMILES, ph=7.4, pdbqt_smiles=NEUTRAL_SMILES)
    from docking.project_layout import PAIRLIST_COLUMNS

    pd.DataFrame(columns=PAIRLIST_COLUMNS).to_csv(project / "pairlist.csv", index=False)
    return project


def test_auto_maps_are_generated_from_the_project_and_hashed(tmp_path):
    project = _auto_project(tmp_path)
    resolved = resolve_md_input_maps(project, "vina", auto_maps=True)
    maps_dir = project / ".meta" / "md_inputs_maps"
    assert Path(resolved["tags_file"]).parent == maps_dir
    tags = pd.read_csv(resolved["tags_file"], dtype=str)
    assert tags["tag"].tolist() == ["prot_site_sti"]
    receptor_map = pd.read_csv(resolved["receptor_map"], dtype=str)
    assert receptor_map.loc[0, "tag"] == "prot_site_sti"
    assert Path(receptor_map.loc[0, "receptor_file"]) == (project / "0-Input" / "receptors" / "prot_chain_A.pdb").resolve()
    topology_map = pd.read_csv(resolved["topology_map"], dtype=str, keep_default_na=False)
    assert topology_map.loc[0, "topology_file"] == str(CRYSTAL_SDF.resolve())
    assert topology_map.loc[0, "topology_pose"] == "1"
    assert topology_map.loc[0, "atom_map"] == ""

    manifest = json.loads(Path(resolved["manifest_file"]).read_text(encoding="utf-8"))
    from post_docking_analysis.md_inputs import sha256_file

    for name in ("tags", "receptor_map", "topology_map"):
        record = manifest["maps"][name]
        assert record["sha256"] == sha256_file(Path(record["path"]))
    assert manifest["sources"]["receptors"][0]["receptor_source_sha256"] == sha256_file(
        project / "0-Input" / "receptors" / "prot_chain_A.pdb"
    )
    assert manifest["sources"]["ligands"][0]["topology_source_sha256"] == sha256_file(CRYSTAL_SDF)


def test_gnina_does_not_generate_a_topology_map(tmp_path):
    project = _auto_project(tmp_path)
    selection = project / "4-Working" / "scores" / "unified" / "best_pose_per_tag_by_engine.csv"
    frame = pd.read_csv(selection)
    frame["engine"] = "gnina"
    frame.to_csv(selection, index=False)
    resolved = resolve_md_input_maps(project, "gnina", auto_maps=True)
    assert resolved["topology_map"] == ""
    assert resolved["sources"]["topology_map"] == "not_required_for_gnina"


def test_explicit_map_files_are_kept_when_only_missing_maps_are_generated(tmp_path):
    project = _auto_project(tmp_path)
    tags = project / "explicit_tags.csv"
    pd.DataFrame([{"tag": "prot_site_sti"}]).to_csv(tags, index=False)
    resolved = resolve_md_input_maps(project, "vina", tags_file=tags)
    assert resolved["tags_file"] == str(tags)
    assert resolved["sources"]["tags"] == "explicit"
    assert resolved["sources"]["receptor_map"] == "auto_generated"
    assert resolved["sources"]["topology_map"] == "auto_generated"


def test_auto_maps_refuse_explicit_files_as_ambiguous(tmp_path):
    project = _auto_project(tmp_path)
    tags = project / "explicit_tags.csv"
    pd.DataFrame([{"tag": "prot_site_sti"}]).to_csv(tags, index=False)
    with pytest.raises(MDMapError) as caught:
        resolve_md_input_maps(project, "vina", tags_file=tags, auto_maps=True)
    assert caught.value.reason == "auto_maps_conflicts_with_explicit_map_file"


def test_missing_ligand_provenance_fails_with_explicit_reason(tmp_path):
    project = _auto_project(tmp_path)
    shutil.rmtree(_prepared_dirs_for(project)[1] / "preparation_steps")
    with pytest.raises(MDMapError) as caught:
        resolve_md_input_maps(project, "vina", auto_maps=True)
    assert caught.value.reason == "ligand_provenance_missing:sti"


def test_ambiguous_ligand_provenance_fails_instead_of_guessing(tmp_path):
    project = _auto_project(tmp_path)
    report_dir = _prepared_dirs_for(project)[1] / "preparation_steps"
    shutil.copy2(report_dir / "sti.json", report_dir / "sti_duplicate.json")
    with pytest.raises(MDMapError) as caught:
        resolve_md_input_maps(project, "vina", auto_maps=True)
    assert caught.value.reason == "ligand_provenance_ambiguous:sti:2"


def test_missing_receptor_lineage_fails_with_explicit_reason(tmp_path):
    project = _auto_project(tmp_path)
    (_prepared_dirs_for(project)[0] / "prot.pdbqt.preparation.json").unlink()
    with pytest.raises(MDMapError) as caught:
        resolve_md_input_maps(project, "vina", auto_maps=True)
    assert caught.value.reason.startswith("receptor_preparation_record_missing")


def test_generated_maps_drive_a_full_docked_export(tmp_path):
    project = _auto_project(tmp_path)
    resolved = resolve_md_input_maps(project, "vina", auto_maps=True)
    if not (openbabel_capability()["status"] == "completed" and HAS_RDKIT):
        pytest.skip("skipped_missing_dependency: needs Open Babel and RDKit")
    request = MDInputsRequest.from_dict(
        {
            "project_dir": str(project),
            "engine": "vina",
            "tags_file": resolved["tags_file"],
            "receptor_map": resolved["receptor_map"],
            "topology_map": resolved["topology_map"],
            "ligand_formats": ["mol2"],
        }
    )
    result = export_md_inputs(request, chemistry_backend=OpenBabelChemistryBackend())
    assert result["status"] == "completed", result["rows"][0]
    assert _provenance(result)["charge_authority"] == "predicted"


# ---------------------------------------------------------------- R1a: CLI inheritance and interactive once-per-project


def test_prepare_commands_inherit_the_project_policy_when_values_are_not_given(tmp_path):
    from workflow.cli import _resolve_prepare_protonation, build_parser

    project = tmp_path / "project"
    _write_manifest(project)
    set_protonation_policy(project, receptor_ph=7.4, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.0)
    parser = build_parser()
    args = parser.parse_args(
        ["pdb", "prepare-both", "--project-dir", str(project), "--ph", "7.4"]
    )
    # The receptor and ligand pH must agree for one run; here they do not, so the run is refused.
    with pytest.raises(SystemExit):
        _resolve_prepare_protonation(args, parser)

    args = parser.parse_args(["pdb", "prepare-protein", "--project-dir", str(project)])
    resolved = _resolve_prepare_protonation(args, parser)
    assert resolved["ph"] == 7.4
    assert resolved["ph_source"] == "user_entered"
    assert resolved["force_field"] == "AMBER"
    assert resolved["protonation_policy"] == "ph_model"

    args = parser.parse_args(["pdb", "prepare-ligand", "--project-dir", str(project)])
    resolved = _resolve_prepare_protonation(args, parser)
    assert resolved["ph"] == 7.0

    # Without a project, the Spec 034 rule holds: the value must be given explicitly.
    args = parser.parse_args(["pdb", "prepare-protein"])
    with pytest.raises(SystemExit):
        _resolve_prepare_protonation(args, parser)


def test_interactive_policy_is_asked_once_and_stored(tmp_path, monkeypatch):
    import workflow.interactive as interactive

    project = tmp_path / "project"
    _write_manifest(project)
    state_map = _write_state_map(project / "states.csv", PLUS1_SMILES)
    text_answers = iter(["7.4", str(state_map)])  # receptor pH, then the state map (explicit_state asks no ligand pH)
    select_answers = iter(["AMBER", "explicit_state"])
    asked: list[str] = []

    def fake_text(questionary, message, default=None, required=False):
        asked.append(message)
        return next(text_answers)

    def fake_select(questionary, message, choices, default=None):
        asked.append(message)
        return next(select_answers)

    monkeypatch.setattr(interactive, "_text", fake_text)
    monkeypatch.setattr(interactive, "_select", fake_select)
    monkeypatch.setattr(interactive, "_header", lambda title: None)

    block = interactive._ask_protonation_policy(None, project)
    assert block is not None
    assert block["receptor_ph"] == 7.4
    assert block["receptor_force_field"] == "AMBER"
    assert block["ligand_policy"] == "explicit_state"
    assert load_protonation_policy(project) == block
    assert any(question.startswith("pH (suggested 7.4") and "PDB2PQR receptor" in question for question in asked)
    assert "Ligand protonation policy" in asked and "Path to the protonation state map (CSV)" in asked


@REQUIRES_OBABEL
def test_cli_md_inputs_with_auto_maps_runs_the_dag_and_completes(tmp_path):
    from workflow.cli import main

    project = _auto_project(tmp_path)
    code = main(["analyze", "md-inputs", "--project-dir", str(project), "--engine", "vina", "--auto-maps"])
    manifest = json.loads((project / "5-Analysis" / "md_inputs" / "md_inputs_manifest.json").read_text(encoding="utf-8"))
    assert code == 0, manifest["rows"]
    assert manifest["status"] == "completed"
    assert manifest["counts"]["completed"] == 1
    assert manifest["rows"][0]["gates"]["G6"]["charge_authority"] == "predicted"
    assert (project / ".meta" / "md_inputs_maps" / "auto_maps_manifest.json").is_file()


@pytest.mark.parametrize("layout_profile", ["canonical", "docking_legacy"])
def test_auto_maps_resolve_receptor_lineage_and_ligand_provenance_through_layout(tmp_path, layout_profile):
    """Spec 036 R6: the auto maps find the receptor lineage and the ligand provenance through the layout."""
    project = _auto_project(tmp_path, layout_profile)
    resolved = resolve_md_input_maps(project, "vina", auto_maps=True)
    receptor_map = pd.read_csv(resolved["receptor_map"], dtype=str)
    assert Path(receptor_map.loc[0, "receptor_file"]).name == "prot_chain_A.pdb"
    topology_map = pd.read_csv(resolved["topology_map"], dtype=str)
    assert Path(topology_map.loc[0, "topology_file"]).resolve() == CRYSTAL_SDF.resolve()


def test_auto_maps_do_not_read_numbered_directories_of_a_canonical_project(tmp_path):
    """A canonical project is read from prepared_proteins only; a numbered 1-Preparation copy is ignored."""
    project = _auto_project(tmp_path, "canonical")
    shutil.rmtree(project / "prepared_proteins")
    stale = project / "1-Preparation" / "receptors" / "prepared"
    stale.mkdir(parents=True)
    (stale / "prot.pdbqt").write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00    -0.300 N\n", encoding="utf-8")
    with pytest.raises(MDMapError) as caught:
        resolve_md_input_maps(project, "vina", auto_maps=True)
    assert caught.value.reason == "receptor_lineage_missing:prot"


def test_manifest_prepared_dir_outside_the_project_is_not_used(tmp_path):
    """A copied project whose manifest points at another project's prepared directory must not read it."""
    project = _auto_project(tmp_path / "copy", "canonical")
    _auto_project(tmp_path / "original", "canonical")
    shutil.rmtree(project / "prepared_proteins")
    manifest_path = project / "project_manifest.json"
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["prepared_proteins_dir"] = str(tmp_path / "original" / "auto_project" / "prepared_proteins")
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(MDMapError) as caught:
        resolve_md_input_maps(project, "vina", auto_maps=True)
    assert caught.value.reason == "receptor_lineage_missing:prot"
