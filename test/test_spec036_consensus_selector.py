"""Spec 036 R2, R3, R4 and R1c / R6 contract tests.

- R3: removed consensus modes fail with a message naming consensus v2; no silent fallback;
  single-engine ranking is explicit; retired legacy RMSD is gone.
- R4: one pose selector. GNINA picks the highest cnn_score (not the lowest cnn_affinity), and
  every writer agrees on a mixed-engine fixture.
- R2: per-target redocking status and its attachment to tables.
- R1c: docked microspecies from the pose REMARK SMILES, then provenance, then unavailable.
- R3d: enhanced RMSD plots carry the diagnostic_not_spec031 label.
"""
from __future__ import annotations

import inspect
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from post_docking_analysis import consensus as consensus_module  # noqa: E402
from post_docking_analysis import geometric_consensus  # noqa: E402
from post_docking_analysis.consensus import (  # noqa: E402
    CONSENSUS_V2_METHOD,
    SINGLE_ENGINE_NATIVE_METHOD,
    build_consensus_rankings,
    normalize_consensus_mode,
    select_representative_poses_v2,
)
from post_docking_analysis.docked_microspecies import (  # noqa: E402
    DOCKED_STATE_COLUMNS,
    annotate_docked_microspecies,
    read_docked_microspecies,
)
from post_docking_analysis.enhanced_rmsd_analyzer import DIAGNOSTIC_LABEL, _stamp_diagnostic_label  # noqa: E402
from post_docking_analysis.multi_engine_pipeline_impl import MultiEngineAnalysisPipeline  # noqa: E402
from post_docking_analysis.pose_selection import select_best_pose_rows  # noqa: E402
from post_docking_analysis.redocking_validation import (  # noqa: E402
    REDOCKING_STATUS_COLUMNS,
    attach_redocking_status,
    build_target_redocking_status,
)
from post_docking_analysis.top_pose_selector import build_top_pose_atlas  # noqa: E402


REMOVED_MODES = ["dockbox_geometric", "weighted_hybrid", "strict_consensus", "favorite_guardrails"]


def _gnina_worst_pose_fixture() -> pd.DataFrame:
    """Pose 2 has the lowest cnn_affinity (the bug's pick). Pose 1 has the highest cnn_score (correct)."""
    return pd.DataFrame(
        [
            {"engine": "gnina", "protein": "p", "tag": "t", "ligand": "L", "site_id": "s", "pose": 1,
             "affinity_kcal_mol": -7.0, "cnn_score": 0.95, "cnn_affinity": 6.0},
            {"engine": "gnina", "protein": "p", "tag": "t", "ligand": "L", "site_id": "s", "pose": 2,
             "affinity_kcal_mol": -9.5, "cnn_score": 0.40, "cnn_affinity": 4.5},
            {"engine": "gnina", "protein": "p", "tag": "t", "ligand": "L", "site_id": "s", "pose": 3,
             "affinity_kcal_mol": -8.0, "cnn_score": 0.70, "cnn_affinity": 7.0},
        ]
    )


def _mixed_engine_fixture() -> pd.DataFrame:
    rows = []
    for tag, ligand, base in (("p_site_1_A_one", "A", 0), ("p_site_1_B_two", "B", 1)):
        rows.append({"engine": "gnina", "protein": "p", "tag": tag, "ligand": ligand, "site_id": "site_1",
                     "pose": 1, "affinity_kcal_mol": -7.0 - base, "cnn_score": 0.2 + 0.1 * base,
                     "cnn_affinity": 5.0 + base})
        rows.append({"engine": "gnina", "protein": "p", "tag": tag, "ligand": ligand, "site_id": "site_1",
                     "pose": 2, "affinity_kcal_mol": -9.5 + base, "cnn_score": 0.9 - 0.4 * base,
                     "cnn_affinity": 4.5 - base})
        rows.append({"engine": "vina", "protein": "p", "tag": tag, "ligand": ligand, "site_id": "site_1",
                     "pose": 1, "affinity_kcal_mol": -8.0 - base})
        rows.append({"engine": "vina", "protein": "p", "tag": tag, "ligand": ligand, "site_id": "site_1",
                     "pose": 2, "affinity_kcal_mol": -9.0 - base})
    return pd.DataFrame(rows)


def _selected_pose(frame: pd.DataFrame, engine: str, tag: str) -> int:
    rows = frame[(frame["engine"] == engine) & (frame["tag"] == tag)]
    return int(rows.iloc[0]["pose"])


# ---------------------------------------------------------------- R4: one selector

def test_gnina_selector_picks_highest_cnn_score_not_lowest_cnn_affinity() -> None:
    selected = select_best_pose_rows(_gnina_worst_pose_fixture())
    assert len(selected) == 1
    assert int(selected.iloc[0]["pose"]) == 1
    assert selected.iloc[0]["v2_ranking_metric_name"] == "cnn_affinity"
    assert selected.iloc[0]["v2_pose_selection_status"] == "completed"


def test_gnina_worst_pose_bug_is_gone_in_every_writer() -> None:
    fixture = _gnina_worst_pose_fixture()
    obj = object.__new__(MultiEngineAnalysisPipeline)
    obj.best_pose_selection_metric = "auto"
    single_writer, metric = MultiEngineAnalysisPipeline._select_best_pose_rows(
        obj, fixture, requested_metric="cnn_affinity"
    )
    assert int(single_writer.iloc[0]["pose"]) == 1
    assert metric == "cnn_score"
    v2 = select_representative_poses_v2(fixture)
    assert int(v2.iloc[0]["pose"]) == 1
    by_tag = MultiEngineAnalysisPipeline._best_by_tag(fixture)
    assert int(by_tag.iloc[0]["pose"]) == 1


def test_top_pose_atlas_uses_the_gnina_selected_pose() -> None:
    fixture = _gnina_worst_pose_fixture()
    atlas = build_top_pose_atlas(best_by_engine=fixture, consensus_df=None, selection_policy="best_affinity")
    per_protein = atlas["top_pose_per_ligand_per_protein"]
    assert len(per_protein) == 1
    assert int(per_protein.iloc[0]["pose"]) == 1


def test_all_writers_agree_on_mixed_engine_fixture() -> None:
    fixture = _mixed_engine_fixture()
    selected = select_best_pose_rows(fixture)
    v2 = select_representative_poses_v2(fixture)
    obj = object.__new__(MultiEngineAnalysisPipeline)
    obj.best_pose_selection_metric = "auto"
    gnina_only = fixture[fixture["engine"] == "gnina"]
    single_writer, _ = MultiEngineAnalysisPipeline._select_best_pose_rows(obj, gnina_only)

    expected = {}
    for engine in ("gnina", "vina"):
        for tag in ("p_site_1_A_one", "p_site_1_B_two"):
            expected[(engine, tag)] = _selected_pose(selected, engine, tag)
    for key in expected:
        assert _selected_pose(v2, *key) == expected[key], key
    for tag in ("p_site_1_A_one", "p_site_1_B_two"):
        assert _selected_pose(single_writer, "gnina", tag) == expected[("gnina", tag)], tag
    # GNINA: highest cnn_score is pose 1 for A (0.2) vs pose 2 (0.9); the selector gives pose 2.
    assert expected[("gnina", "p_site_1_A_one")] == 2
    # Vina: lowest affinity is pose 2 for both tags.
    assert expected[("vina", "p_site_1_A_one")] == 2 and expected[("vina", "p_site_1_B_two")] == 2


def test_selector_keeps_auditable_row_when_gnina_cnn_score_missing() -> None:
    frame = pd.DataFrame(
        [{"engine": "gnina", "protein": "p", "tag": "t", "ligand": "L", "site_id": "s", "pose": 1,
          "affinity_kcal_mol": -9.0, "cnn_score": np.nan, "cnn_affinity": np.nan}]
    )
    selected = select_best_pose_rows(frame)
    assert selected.iloc[0]["v2_pose_selection_status"] == "missing_required_cnn_score"
    assert pd.isna(selected.iloc[0]["v2_ranking_metric"])


# ---------------------------------------------------------------- R3: consensus modes

def test_default_consensus_mode_is_v2_everywhere() -> None:
    assert normalize_consensus_mode("") == CONSENSUS_V2_METHOD
    assert normalize_consensus_mode(None) == CONSENSUS_V2_METHOD
    default = inspect.signature(build_consensus_rankings).parameters["consensus_mode"].default
    assert default == CONSENSUS_V2_METHOD
    init_params = inspect.signature(MultiEngineAnalysisPipeline.__init__).parameters
    assert init_params["consensus_mode"].default == CONSENSUS_V2_METHOD


@pytest.mark.parametrize("mode", REMOVED_MODES)
def test_removed_mode_raises_and_names_v2(mode: str) -> None:
    with pytest.raises(ValueError) as excinfo:
        normalize_consensus_mode(mode)
    assert "consensus_rank_geometry_qc_v2" in str(excinfo.value)
    assert mode in str(excinfo.value)


def test_unknown_mode_raises_instead_of_silent_legacy_fallback() -> None:
    with pytest.raises(ValueError) as excinfo:
        normalize_consensus_mode("made_up_mode")
    assert "consensus_rank_geometry_qc_v2" in str(excinfo.value)


def test_removed_mode_raises_in_build_consensus_rankings() -> None:
    with pytest.raises(ValueError, match="consensus_rank_geometry_qc_v2"):
        build_consensus_rankings(_gnina_worst_pose_fixture(), consensus_mode="dockbox_geometric")


def test_legacy_consensus_helpers_are_retired() -> None:
    assert not hasattr(geometric_consensus, "compute_geometric_consensus")
    assert not hasattr(geometric_consensus, "_pairwise_rmsd")
    assert not hasattr(geometric_consensus, "_centroid_sorted")
    assert not hasattr(consensus_module, "compute_geometric_consensus")
    assert set(consensus_module.SUPPORTED_CONSENSUS_MODES) == {CONSENSUS_V2_METHOD, SINGLE_ENGINE_NATIVE_METHOD}


def _single_vina_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"engine": "vina", "protein": "P1", "ligand": "A", "site_id": "site_1", "tag": "P1_site_1_A",
             "pose": 1, "affinity_kcal_mol": -10.0},
            {"engine": "vina", "protein": "P1", "ligand": "B", "site_id": "site_1", "tag": "P1_site_1_B",
             "pose": 1, "affinity_kcal_mol": -8.0},
            {"engine": "vina", "protein": "P1", "ligand": "C", "site_id": "site_1", "tag": "P1_site_1_C",
             "pose": 1, "affinity_kcal_mol": -6.0},
        ]
    )


def test_single_engine_native_ranks_by_native_score_and_is_labelled() -> None:
    ranked = build_consensus_rankings(_single_vina_frame(), consensus_mode=SINGLE_ENGINE_NATIVE_METHOD)
    assert list(ranked["tag"]) == ["P1_site_1_A", "P1_site_1_B", "P1_site_1_C"]
    assert set(ranked["consensus_mode"]) == {SINGLE_ENGINE_NATIVE_METHOD}
    assert set(ranked["consensus_status"]) == {"single_engine"}
    assert ranked["agreement_fraction"].isna().all()
    assert bool(ranked["single_engine"].all())


def test_single_engine_native_refuses_multi_engine_input() -> None:
    frame = pd.concat([_single_vina_frame(), _gnina_worst_pose_fixture()], ignore_index=True)
    with pytest.raises(ValueError, match="exactly one engine"):
        build_consensus_rankings(frame, consensus_mode=SINGLE_ENGINE_NATIVE_METHOD)


def test_v2_multi_engine_output_is_unchanged_by_the_shared_selector() -> None:
    fixture = _mixed_engine_fixture()
    ranked = build_consensus_rankings(fixture, consensus_mode=CONSENSUS_V2_METHOD)
    assert set(ranked["consensus_mode"]) == {CONSENSUS_V2_METHOD}
    assert set(ranked["consensus_status"]) <= {"completed", "single_engine", "consensus_incomplete"}


def test_enhanced_rmsd_plots_carry_the_diagnostic_label() -> None:
    assert DIAGNOSTIC_LABEL == "diagnostic_not_spec031"
    fig, ax = plt.subplots()
    _stamp_diagnostic_label()
    texts = [text.get_text() for text in fig.texts]
    plt.close(fig)
    assert any(DIAGNOSTIC_LABEL in text for text in texts)


# ---------------------------------------------------------------- R2: redocking status

def _validation_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"protein": "P_pass", "tag": "ref_good", "engine": "vina", "redocking_classification": "pass",
             "redocking_rmsd_angstrom": 1.2, "validation_reason": ""},
            {"protein": "P_pass", "tag": "ref_bad", "engine": "vina", "redocking_classification": "fail",
             "redocking_rmsd_angstrom": 6.0, "validation_reason": ""},
            {"protein": "P_fail", "tag": "ref_far", "engine": "vina", "redocking_classification": "fail",
             "redocking_rmsd_angstrom": 5.0, "validation_reason": ""},
            {"protein": "P_warn", "tag": "ref_warn", "engine": "vina", "redocking_classification": "warn",
             "redocking_rmsd_angstrom": 2.8, "validation_reason": ""},
            {"protein": "P_nc", "tag": "ref_nc", "engine": "vina", "redocking_classification": "not_evaluable",
             "redocking_rmsd_angstrom": np.nan, "validation_reason": "not_comparable:graphs_not_isomorphic"},
            {"protein": "P_missing", "tag": "ref_missing", "engine": "vina",
             "redocking_classification": "not_evaluable", "redocking_rmsd_angstrom": np.nan,
             "validation_reason": "missing_docked_pose_file"},
        ]
    )


def test_redocking_status_per_target() -> None:
    status = build_target_redocking_status(_validation_frame(), targets=["P_pass", "P_fail", "P_warn", "P_nc", "P_missing", "P_none"])
    by = status.set_index("protein")
    assert by.loc["P_pass", "redocking_validation_status"] == "validated"
    assert by.loc["P_pass", "redocking_reference"] == "ref_good;engine=vina"
    assert "1.200A" in by.loc["P_pass", "redocking_validation_reason"]
    assert by.loc["P_fail", "redocking_validation_status"] == "failed"
    assert "above_3.5A" in by.loc["P_fail", "redocking_validation_reason"]
    assert by.loc["P_warn", "redocking_validation_status"] == "failed"
    assert "warn_band" in by.loc["P_warn", "redocking_validation_reason"]
    assert by.loc["P_nc", "redocking_validation_status"] == "not_comparable"
    assert by.loc["P_nc", "redocking_validation_reason"] == "not_comparable:graphs_not_isomorphic"
    assert by.loc["P_missing", "redocking_validation_status"] == "not_evaluated"
    assert by.loc["P_missing", "redocking_validation_reason"] == "missing_docked_pose_file"
    assert by.loc["P_none", "redocking_validation_status"] == "not_evaluated"
    assert by.loc["P_none", "redocking_validation_reason"] == "no_reference_rows_for_target"


def test_redocking_status_attaches_to_tables_and_marks_unvalidated_targets() -> None:
    status = build_target_redocking_status(_validation_frame())
    table = pd.DataFrame({"protein": ["P_pass", "P_fail", "P_unknown"], "tag": ["a", "b", "c"]})
    attached = attach_redocking_status(table, status)
    for column in REDOCKING_STATUS_COLUMNS:
        assert column in attached.columns
    assert list(attached["redocking_validation_status"]) == ["validated", "failed", "not_evaluated"]
    assert attached.loc[2, "redocking_validation_reason"] == "no_reference_rows_for_target"
    assert attached.loc[0, "redocking_reference"] == "ref_good;engine=vina"


def test_redocking_attach_replaces_existing_columns() -> None:
    status = build_target_redocking_status(_validation_frame())
    table = pd.DataFrame({"protein": ["P_pass"], "redocking_validation_status": ["stale"]})
    attached = attach_redocking_status(table, status)
    assert list(attached["redocking_validation_status"]) == ["validated"]


def test_redocking_status_empty_validation_is_not_evaluated() -> None:
    status = build_target_redocking_status(pd.DataFrame(), targets=["P1"])
    assert list(status["redocking_validation_status"]) == ["not_evaluated"]


# ---------------------------------------------------------------- R1c: docked microspecies

_CHARGED_SMILES = "CC[NH3+]"  # net +1, docked microspecies for the fixture


def _write_pdbqt(path: Path, *, remark_smiles: str | None) -> None:
    # Meeko remarks sit inside each MODEL block, as Vina copies them into every output model.
    lines = ["MODEL 1", "REMARK VINA RESULT:    -7.0      0.000      0.000"]
    if remark_smiles:
        lines.append(f"REMARK SMILES {remark_smiles}")
        lines.append("REMARK SMILES IDX 1 1 2 2 3 3")
    lines.append("ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00    +0.000 C")
    lines.append("ENDMDL")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_provenance(project: Path, stem: str, *, smiles: str, charge: int, policy: str) -> None:
    folder = project / "prepared_ligands" / "preparation_steps"
    folder.mkdir(parents=True, exist_ok=True)
    payload = {
        "protonation_policy": policy,
        "protonation": {"pdbqt_state": {"microspecies_smiles": smiles, "net_formal_charge": charge}},
    }
    (folder / f"{stem}.json").write_text(json.dumps(payload), encoding="utf-8")


def test_docked_microspecies_from_pose_remark_smiles(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    project = tmp_path / "project"
    project.mkdir()
    pose = tmp_path / "pose.pdbqt"
    _write_pdbqt(pose, remark_smiles=_CHARGED_SMILES)
    result = read_docked_microspecies(project, ligand="LIG.pdbqt", pose_file=str(pose), pose_index=1)
    assert result["docked_microspecies_smiles"] == _CHARGED_SMILES
    assert int(result["docked_net_charge"]) == 1
    assert result["docked_state_source"] == "pose_remark_smiles"


def test_docked_microspecies_falls_back_to_preparation_provenance(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    project = tmp_path / "project"
    _write_provenance(project, "LIG", smiles=_CHARGED_SMILES, charge=1, policy="explicit_state")
    pose = tmp_path / "pose.pdbqt"
    _write_pdbqt(pose, remark_smiles=None)
    result = read_docked_microspecies(project, ligand="LIG.pdbqt", pose_file=str(pose), pose_index=1)
    assert result["docked_microspecies_smiles"] == _CHARGED_SMILES
    assert int(result["docked_net_charge"]) == 1
    assert result["protonation_policy"] == "explicit_state"
    assert result["docked_state_source"] == "preparation_provenance"


def test_docked_microspecies_unavailable_gives_explicit_empty_values(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    project = tmp_path / "project"
    project.mkdir()
    pose = tmp_path / "pose.pdbqt"
    _write_pdbqt(pose, remark_smiles=None)
    result = read_docked_microspecies(project, ligand="LIG.pdbqt", pose_file=str(pose), pose_index=1)
    assert result == {
        "docked_microspecies_smiles": "",
        "docked_net_charge": "",
        "protonation_policy": "",
        "docked_state_source": "unavailable",
    }


def test_annotate_docked_microspecies_adds_columns_to_best_pose_table(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    project = tmp_path / "project"
    project.mkdir()
    pose = tmp_path / "pose.pdbqt"
    _write_pdbqt(pose, remark_smiles=_CHARGED_SMILES)
    frame = pd.DataFrame([
        {"engine": "vina", "tag": "t1", "ligand": "LIG.pdbqt", "pose_file": str(pose), "pose": 1},
        {"engine": "vina", "tag": "t2", "ligand": "MISSING.pdbqt", "pose_file": "", "pose": 1},
    ])
    annotated = annotate_docked_microspecies(frame, project)
    for column in DOCKED_STATE_COLUMNS:
        assert column in annotated.columns
    assert annotated.loc[0, "docked_state_source"] == "pose_remark_smiles"
    assert annotated.loc[1, "docked_state_source"] == "unavailable"
    assert annotated.loc[1, "docked_microspecies_smiles"] == ""
