"""Spec 036 integration follow-ups (T008): manifest bootstrap, pairlist columns, GNINA pose export,
best_engine_per_complex, replicate pooling and pose reproducibility (R5b), box default in the webui.

Replicate fixtures are tmp copies of ``test/fixtures/spec034/sti_vina_model1.pdbqt``. Each copy gets
its own ``REMARK VINA RESULT`` affinity and a translated x coordinate. The committed fixture is never changed.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from docking.models import PairlistRow
from docking.preparation.pairlist_builder import _write_csv
from docking.project_layout import (
    PAIRLIST_COLUMNS,
    bootstrap_project_layout,
    ensure_engine_layout,
    ensure_pairlist_stub,
    load_protonation_policy,
    manifest_path,
    pairlist_path,
    set_protonation_policy,
)
from docking.runners.base import build_pair_index, replicate_job_tag, split_replicate_stem
from post_docking_analysis.generate_scores_csv import _resolve_tag_from_pairlist, load_pairlist_mapping
from post_docking_analysis.engine_detector import _pose_pair_count
from post_docking_analysis.multi_engine_pipeline_impl import MultiEngineAnalysisPipeline
from post_docking_analysis.pose_extractor import extract_best_poses_from_gnina
from post_docking_analysis.pose_selection import select_best_pose_rows
from post_docking_analysis.replicates import (
    load_replicate_seeds,
    pose_reproducibility_table,
    split_pose_stem,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
STI_FIXTURE = REPO_ROOT / "test" / "fixtures" / "spec034" / "sti_vina_model1.pdbqt"
RECEPTOR_FILE = "1IEP_A_protein.pdbqt"
LIGAND_FILE = "STI_A_201.pdbqt"
SITE = "site_1"
PAIR_TAG = f"1IEP_A_protein_{SITE}_{LIGAND_FILE}"  # receptor-stem tag (Spec 036 R6)
LEGACY_TAG = f"{RECEPTOR_FILE}_{SITE}_{LIGAND_FILE}"
BASE_SEED = 11


# ----------------------------------------------------------------------------- helpers

def _shifted_vina_pose(text: str, shift_x: float, affinity: float) -> str:
    """One Vina MODEL: the fixture coordinates moved by ``shift_x`` on x, with a REMARK VINA RESULT line."""
    lines = text.splitlines()
    assert lines[0] == "MODEL 1"
    out = ["MODEL 1", f"REMARK VINA RESULT:    {affinity:.3f}      0.000      0.000"]
    for line in lines[1:]:
        if line.startswith(("ATOM", "HETATM")):
            x = float(line[30:38]) + shift_x
            line = f"{line[:30]}{x:8.3f}{line[38:]}"
        out.append(line)
    return "\n".join(out) + "\n"


def _write_replicate_poses(poses_dir: Path, specs) -> None:
    """specs: [(replicate_id, shift_x, affinity)]. Writes ``<pair>__repNN.pdbqt`` into poses_dir."""
    text = STI_FIXTURE.read_text(encoding="utf-8")
    poses_dir.mkdir(parents=True, exist_ok=True)
    for replicate_id, shift_x, affinity in specs:
        name = f"{replicate_job_tag(PAIR_TAG, replicate_id)}.pdbqt"
        (poses_dir / name).write_text(_shifted_vina_pose(text, shift_x, affinity), encoding="utf-8")


def _vina_project(root: Path, replicate_specs) -> tuple[Path, Path, dict]:
    """A canonical vina project with one pair and the given replicate poses plus a run manifest."""
    bootstrap_project_layout(root, engines=["vina"], project_name="spec036-followups",
                             favorite_engine="vina", layout_profile="canonical")
    pd.DataFrame([{
        "receptor": RECEPTOR_FILE, "site_id": SITE, "ligand": LIGAND_FILE,
        "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
        "size_x": 20.0, "size_y": 20.0, "size_z": 20.0,
        "protein_display_name": "1IEP_A_protein", "ligand_display_name": "STI_A_201",
    }]).to_csv(pairlist_path(root), index=False)
    layout = ensure_engine_layout(root, "vina", layout_profile="canonical")
    _write_replicate_poses(layout["poses"], replicate_specs)
    seeds = [BASE_SEED + index for index in range(len(replicate_specs))]
    manifest = {
        "engine": "vina",
        "effective": {"engine": "vina", "replicates": len(seeds), "seeds": seeds},
        "jobs": [
            {"replicate_id": rid, "seed": BASE_SEED + rid - 1, "pair_tag": PAIR_TAG,
             "tag": replicate_job_tag(PAIR_TAG, rid), "status": "completed"}
            for rid, _, _ in replicate_specs
        ],
    }
    (layout["root"] / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return root, layout["poses"], manifest


def _pipeline(project_dir: Path, output_dir: Path) -> MultiEngineAnalysisPipeline:
    return MultiEngineAnalysisPipeline(
        project_dir=str(project_dir),
        output_dir=str(output_dir),
        analysis_mode="comparative_all_engines",
        analysis_scope="comparison_only",
        normalization_method="per_engine_rank",
        favorite_engine="vina",
        hit_class_policy="target_percentile",
        hit_class_strong_percentile=0.10,
        hit_class_moderate_percentile=0.40,
        top_pose_selection_policy="best_affinity",
        top_pose_global_aggregation="best_target",
    )


def _gnina_sdf_record(name: str) -> str:
    return (f"{name}\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n")


def _gnina_extract(tmp_path: Path, rows: list[dict], sdf_tags: dict[str, int], criterion=None) -> pd.DataFrame:
    """Run the GNINA pose extractor on synthetic all_scores rows. sdf_tags: source tag -> pose count."""
    inp = tmp_path / "input"
    gdir = inp / "gnina_out"
    gdir.mkdir(parents=True, exist_ok=True)
    for source_tag, count in sdf_tags.items():
        (gdir / f"{source_tag}_top.sdf").write_text(
            "".join(_gnina_sdf_record(f"{source_tag}_{i}") for i in range(1, count + 1)), encoding="utf-8")
    scores = tmp_path / "all_scores.csv"
    pd.DataFrame(rows).to_csv(scores, index=False)
    out = tmp_path / "out"
    extract_best_poses_from_gnina(inp, out, config={}, gnina_dir=gdir, receptors_dir=inp / "receptors",
                                  scores_csv=scores, best_pose_criterion=criterion)
    return pd.read_csv(out / "pose_extraction_manifest.csv")


# ----------------------------------------------------------------------------- item 1: bootstrap keeps policy

def test_bootstrap_again_keeps_protonation_policy_and_unknown_keys(tmp_path):
    root = tmp_path / "project"
    bootstrap_project_layout(root, engines=["vina"], layout_profile="canonical")
    set_protonation_policy(root, receptor_ph=7.4, receptor_force_field="amber",
                           ligand_policy="ph_model", ligand_ph=7.4)
    path = manifest_path(root, "canonical")
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["custom_unknown_key"] = {"kept": True}
    path.write_text(json.dumps(payload), encoding="utf-8")
    before = load_protonation_policy(root)
    assert before is not None

    bootstrap_project_layout(root, engines=["vina", "smina"], layout_profile="canonical")

    assert load_protonation_policy(root) == before
    after = json.loads(path.read_text(encoding="utf-8"))
    assert after["custom_unknown_key"] == {"kept": True}
    assert sorted(after["engines"]) == ["smina", "vina"]


def test_workflow_init_on_existing_project_keeps_protonation_policy(tmp_path):
    from workflow import execution

    root = tmp_path / "project"
    bootstrap_project_layout(root, engines=["vina"], layout_profile="canonical")
    set_protonation_policy(root, receptor_ph=7.4, receptor_force_field="amber",
                           ligand_policy="as_input")
    policy = load_protonation_policy(root)
    result = execution.run_workflow_init(str(root), engines=["vina"], layout_profile="canonical")
    assert result.status == "completed"
    assert load_protonation_policy(root) == policy


# ----------------------------------------------------------------------------- item 2: one pairlist column set

def test_pairlist_header_has_box_columns_in_every_writer(tmp_path):
    stub = ensure_pairlist_stub(tmp_path / "project")
    header = stub.read_text(encoding="utf-8").splitlines()[0].split(",")
    assert header == PAIRLIST_COLUMNS
    for column in ("box_method", "ligand_rg_angstrom", "edge_angstrom", "box_containment_status", "box_warnings"):
        assert column in PAIRLIST_COLUMNS


def test_pairlist_rows_with_box_fields_write_without_error(tmp_path):
    row = PairlistRow(receptor=RECEPTOR_FILE, site_id=SITE, ligand=LIGAND_FILE, center_x=1.0, center_y=2.0,
                      center_z=3.0, size_x=18.0, size_y=18.0, size_z=18.0, box_method="rg_scaled_v1",
                      ligand_rg_angstrom=3.1, edge_angstrom=18.0, box_containment_status="contained",
                      box_warnings="")
    out = tmp_path / "pairlist.csv"
    _write_csv(out, [row.to_dict()], PAIRLIST_COLUMNS)
    frame = pd.read_csv(out)
    assert list(frame.columns) == PAIRLIST_COLUMNS
    assert frame.loc[0, "box_method"] == "rg_scaled_v1"


def test_old_pairlist_without_box_columns_still_loads(tmp_path):
    from docking.cli import _pairlist_rows_from_df
    from docking.preparation.project_builder import _recorded_box_provenance

    old = pd.DataFrame([{
        "receptor": RECEPTOR_FILE, "site_id": SITE, "ligand": LIGAND_FILE,
        "center_x": 1.0, "center_y": 2.0, "center_z": 3.0, "size_x": 20.0, "size_y": 20.0, "size_z": 20.0,
    }])
    rows = _pairlist_rows_from_df(old)
    assert rows[0].box_method is None or rows[0].box_method == ""
    assert rows[0].tag == PAIR_TAG
    provenance = _recorded_box_provenance(old.iloc[0])
    assert provenance["box_method"] == ""
    assert provenance["edge_angstrom"] is None


# ----------------------------------------------------------------------------- item 3: GNINA worst-pose regression

def _gnina_worst_pose_rows(tag: str) -> list[dict]:
    """Mirrors test_spec036_consensus_selector._gnina_worst_pose_fixture: pose 2 has the lowest cnn_affinity."""
    return [
        {"tag": tag, "mode": 1, "vina_affinity": -7.0, "cnn_affinity": 6.0, "cnn_score": 0.95},
        {"tag": tag, "mode": 2, "vina_affinity": -9.5, "cnn_affinity": 4.5, "cnn_score": 0.40},
        {"tag": tag, "mode": 3, "vina_affinity": -8.0, "cnn_affinity": 7.0, "cnn_score": 0.70},
    ]


@pytest.mark.parametrize("criterion", [None, "affinity", "cnn", "cnn_affinity"])
def test_gnina_pose_export_selects_highest_cnn_score_for_every_criterion(tmp_path, criterion):
    manifest = _gnina_extract(tmp_path, _gnina_worst_pose_rows("R1_site_1_L1"),
                              {"R1_site_1_L1": 3}, criterion=criterion)
    assert len(manifest) == 1
    assert int(manifest.loc[0, "selected_pose"]) == 1
    assert manifest.loc[0, "selection_key"] == "cnn_score:higher_is_better"
    assert manifest.loc[0, "selection_rule"] == "consensus_rank_geometry_qc_v2/pose_selection"


def test_gnina_pose_export_pools_replicates_and_records_the_source(tmp_path):
    rows = [
        {"tag": "R1_site_1_L1__rep01", "mode": 1, "vina_affinity": -8.0, "cnn_affinity": 5.0, "cnn_score": 0.60},
        {"tag": "R1_site_1_L1__rep01", "mode": 2, "vina_affinity": -9.0, "cnn_affinity": 4.0, "cnn_score": 0.20},
        {"tag": "R1_site_1_L1__rep02", "mode": 1, "vina_affinity": -7.5, "cnn_affinity": 6.5, "cnn_score": 0.90},
        {"tag": "R1_site_1_L1__rep02", "mode": 2, "vina_affinity": -9.4, "cnn_affinity": 3.0, "cnn_score": 0.30},
    ]
    manifest = _gnina_extract(tmp_path, rows, {"R1_site_1_L1__rep01": 2, "R1_site_1_L1__rep02": 2})
    assert len(manifest) == 1
    row = manifest.iloc[0]
    assert row["tag"] == "R1_site_1_L1"
    assert row["source_tag"] == "R1_site_1_L1__rep02"
    assert int(row["source_replicate_id"]) == 2
    assert int(row["selected_pose"]) == 1
    assert float(row["selected_score"]) == pytest.approx(6.5)


def test_gnina_row_without_cnn_score_is_recorded_not_replaced(tmp_path):
    rows = [{"tag": "R1_site_1_L1", "mode": 1, "vina_affinity": -9.0, "cnn_affinity": np.nan, "cnn_score": np.nan},
            {"tag": "R1_site_1_L1", "mode": 2, "vina_affinity": -10.0, "cnn_affinity": np.nan, "cnn_score": np.nan}]
    manifest = _gnina_extract(tmp_path, rows, {"R1_site_1_L1": 2})
    assert manifest.loc[0, "status"] == "missing_required_cnn_score"
    assert not list((tmp_path / "out").rglob("*.pdb")), "no pose may be exported from an incomplete v2 selection"


# ----------------------------------------------------------------------------- item 4: no cross-engine winner

def _two_engine_best(project: Path) -> pd.DataFrame:
    frame = pd.DataFrame([
        {"engine": "gnina", "protein": RECEPTOR_FILE, "tag": PAIR_TAG, "ligand": LIGAND_FILE, "site_id": SITE,
         "pose": 1, "affinity_kcal_mol": -8.0, "cnn_score": 0.5, "cnn_affinity": 5.0, "replicate_id": 1,
         "seed": 1, "pose_file": "x.sdf"},
        {"engine": "vina", "protein": RECEPTOR_FILE, "tag": PAIR_TAG, "ligand": LIGAND_FILE, "site_id": SITE,
         "pose": 2, "affinity_kcal_mol": -7.8, "replicate_id": 2, "seed": 12, "pose_file": "y.pdbqt"},
        {"engine": "vina", "protein": RECEPTOR_FILE, "tag": "1IEP_A_protein_site_1_OTHER.pdbqt", "ligand": "OTHER.pdbqt",
         "site_id": SITE, "pose": 1, "affinity_kcal_mol": -6.0, "replicate_id": 1, "seed": 11, "pose_file": "z"},
    ])
    return frame


def test_best_engine_per_complex_reports_native_scores_side_by_side_without_a_winner():
    table = MultiEngineAnalysisPipeline._per_engine_side_by_side(_two_engine_best(None))
    assert "affinity_kcal_mol" not in table.columns
    assert "engine" not in table.columns
    shared = table.set_index("tag").loc[PAIR_TAG]
    assert shared["comparison_status"] == "not_comparable_across_engines"
    assert shared["gnina_affinity_kcal_mol"] == pytest.approx(-8.0)
    assert shared["vina_affinity_kcal_mol"] == pytest.approx(-7.8)
    assert shared["gnina_cnn_score"] == pytest.approx(0.5)
    assert shared["vina_replicate_id"] == 2
    single = table.set_index("tag").loc["1IEP_A_protein_site_1_OTHER.pdbqt"]
    assert single["comparison_status"] == "single_engine"
    assert single["engines_missing"] == "gnina"


def test_comparative_run_writes_side_by_side_best_engine_per_complex(tmp_path):
    """Full comparative run on a two-engine project (vina and smina poses of one pair)."""
    root = tmp_path / "project"
    bootstrap_project_layout(root, engines=["vina", "smina"], layout_profile="canonical")
    pd.DataFrame([{
        "receptor": RECEPTOR_FILE, "site_id": SITE, "ligand": LIGAND_FILE,
        "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
        "size_x": 20.0, "size_y": 20.0, "size_z": 20.0,
    }]).to_csv(pairlist_path(root), index=False)
    text = STI_FIXTURE.read_text(encoding="utf-8")
    for engine, affinity in (("vina", -9.0), ("smina", -6.0)):
        layout = ensure_engine_layout(root, engine, layout_profile="canonical")
        (layout["poses"] / f"{PAIR_TAG}.pdbqt").write_text(_shifted_vina_pose(text, 0.0, affinity), encoding="utf-8")
        (layout["logs"] / f"{PAIR_TAG}.log").write_text("REMARK\n", encoding="utf-8")
    output = root / "analysis"
    assert _pipeline(root, output).run() is True
    csv_path = next(output.rglob("best_engine_per_complex.csv"))
    frame = pd.read_csv(csv_path)
    assert len(frame) == 1
    assert frame.loc[0, "comparison_status"] == "not_comparable_across_engines"
    assert frame.loc[0, "vina_affinity_kcal_mol"] == pytest.approx(-9.0)
    assert frame.loc[0, "smina_affinity_kcal_mol"] == pytest.approx(-6.0)
    assert "affinity_kcal_mol" not in frame.columns


# ----------------------------------------------------------------------------- item 5: replicates

@pytest.fixture()
def three_replicates(tmp_path):
    """Replicate 1 and 2 differ by 0.5 A in place, replicate 3 is 4 A away. Affinities: -9.0, -8.5, -9.5."""
    return _vina_project(tmp_path / "project", [(1, 0.0, -9.0), (2, 0.5, -8.5), (3, 4.0, -9.5)])


def test_replicate_stems_and_receptor_stem_tags_parse():
    assert split_replicate_stem(f"{PAIR_TAG}__rep01") == (PAIR_TAG, 1)
    assert split_pose_stem(f"{PAIR_TAG}__rep03") == (PAIR_TAG, 3)
    assert split_replicate_stem(LEGACY_TAG) == (LEGACY_TAG, None)
    assert replicate_job_tag(PAIR_TAG, 2) == f"{PAIR_TAG}__rep02"
    # The receptor stem (no .pdbqt) with the ligand file name (with .pdbqt) is the pair tag.
    row = PairlistRow(receptor=RECEPTOR_FILE, site_id=SITE, ligand=LIGAND_FILE, center_x=0, center_y=0,
                      center_z=0, size_x=1, size_y=1, size_z=1)
    assert row.tag == PAIR_TAG
    assert row.legacy_tag == LEGACY_TAG
    index = build_pair_index([row])
    assert index[PAIR_TAG] is row and index[LEGACY_TAG] is row


def test_replicate_seeds_come_from_the_run_manifest(three_replicates):
    root, poses, _ = three_replicates
    seeds = load_replicate_seeds(root / "engines" / "vina" / "run_manifest.json"
                                 if (root / "engines").exists() else
                                 next(root.rglob("run_manifest.json")))
    assert seeds == {1: 11, 2: 12, 3: 13}


def test_score_import_keeps_every_replicate_pose_under_the_pair_tag(three_replicates):
    root, poses, _ = three_replicates
    pipeline = _pipeline(root, root / "analysis")
    layout = {"root": next(root.rglob("run_manifest.json")).parent, "poses": poses,
              "logs": poses.parent / "logs", "scores": poses.parent / "scores"}
    frame = pipeline._build_normalized_frame("vina", layout, pipeline._pair_index())
    assert len(frame) == 3
    assert set(frame["tag"]) == {PAIR_TAG}
    assert sorted(frame["replicate_id"].astype(int)) == [1, 2, 3]
    assert sorted(frame["seed"].astype(int)) == [11, 12, 13]
    assert set(frame["pose_file"]) == {str(poses / f"{replicate_job_tag(PAIR_TAG, r)}.pdbqt") for r in (1, 2, 3)}


def test_pooled_selection_picks_the_best_replicate_and_records_it(three_replicates):
    root, poses, _ = three_replicates
    pipeline = _pipeline(root, root / "analysis")
    layout = {"root": next(root.rglob("run_manifest.json")).parent, "poses": poses,
              "logs": poses.parent / "logs", "scores": poses.parent / "scores"}
    normalized = pipeline._build_normalized_frame("vina", layout, pipeline._pair_index())
    best = pipeline._select_best_pose_table(normalized)
    assert len(best) == 1
    row = best.iloc[0]
    assert int(row["replicate_id"]) == 3
    assert int(row["seed"]) == 13
    assert row["pose_file"].endswith(f"{replicate_job_tag(PAIR_TAG, 3)}.pdbqt")
    assert int(row["replicates_pooled"]) == 3
    assert row["affinity_kcal_mol"] == pytest.approx(-9.5, abs=1e-6)


def test_selection_does_not_depend_on_replicate_row_order(three_replicates):
    root, poses, _ = three_replicates
    pipeline = _pipeline(root, root / "analysis")
    layout = {"root": next(root.rglob("run_manifest.json")).parent, "poses": poses,
              "logs": poses.parent / "logs", "scores": poses.parent / "scores"}
    normalized = pipeline._build_normalized_frame("vina", layout, pipeline._pair_index())
    forward = select_best_pose_rows(normalized)
    backward = select_best_pose_rows(normalized.iloc[::-1].reset_index(drop=True))
    assert forward.iloc[0]["pose_file"] == backward.iloc[0]["pose_file"]


def test_three_replicates_give_pairwise_in_place_rmsd_and_fraction_within_two_angstrom(three_replicates):
    root, poses, _ = three_replicates
    pipeline = _pipeline(root, root / "analysis")
    layout = {"root": next(root.rglob("run_manifest.json")).parent, "poses": poses,
              "logs": poses.parent / "logs", "scores": poses.parent / "scores"}
    normalized = pipeline._build_normalized_frame("vina", layout, pipeline._pair_index())
    table = pose_reproducibility_table(normalized)
    assert len(table) == 1
    row = table.iloc[0]
    assert row["pose_reproducibility_status"] == "completed"
    assert int(row["pose_reproducibility_replicates"]) == 3
    assert int(row["pose_reproducibility_pairs"]) == 3
    assert int(row["pose_reproducibility_pairs_comparable"]) == 3
    # In-place shifts: 0.5 (rep1-rep2), 4.0 (rep1-rep3), 3.5 (rep2-rep3).
    assert row["pose_reproducibility_max_rmsd_angstrom"] == pytest.approx(4.0, abs=1e-3)
    assert row["pose_reproducibility_median_rmsd_angstrom"] == pytest.approx(3.5, abs=1e-3)
    assert row["pose_reproducibility_fraction_within_2A"] == pytest.approx(1 / 3)
    assert row["pose_reproducibility_method"] == "spec031_in_place_heavy_atom_rmsd_v1"


def test_best_pose_table_carries_the_reproducibility_metric(three_replicates):
    root, poses, _ = three_replicates
    pipeline = _pipeline(root, root / "analysis")
    layout = {"root": next(root.rglob("run_manifest.json")).parent, "poses": poses,
              "logs": poses.parent / "logs", "scores": poses.parent / "scores"}
    normalized = pipeline._build_normalized_frame("vina", layout, pipeline._pair_index())
    best = pipeline._select_best_pose_table(normalized)
    assert best.iloc[0]["pose_reproducibility_status"] == "completed"
    assert best.iloc[0]["pose_reproducibility_max_rmsd_angstrom"] == pytest.approx(4.0, abs=1e-3)
    banner = "\n".join(pipeline._pose_reproducibility_banner_lines(best))
    assert "within_2A=0.33" in banner


def test_one_replicate_is_single_replicate(tmp_path):
    root, poses, _ = _vina_project(tmp_path / "project", [(1, 0.0, -9.0)])
    pipeline = _pipeline(root, root / "analysis")
    layout = {"root": next(root.rglob("run_manifest.json")).parent, "poses": poses,
              "logs": poses.parent / "logs", "scores": poses.parent / "scores"}
    normalized = pipeline._build_normalized_frame("vina", layout, pipeline._pair_index())
    row = pose_reproducibility_table(normalized).iloc[0]
    assert row["pose_reproducibility_status"] == "single_replicate"
    assert int(row["pose_reproducibility_replicates"]) == 1
    assert math.isnan(row["pose_reproducibility_max_rmsd_angstrom"])


def test_legacy_single_pose_names_still_resolve_to_their_pair(tmp_path):
    root = tmp_path / "project"
    bootstrap_project_layout(root, engines=["vina"], layout_profile="canonical")
    pd.DataFrame([{"receptor": RECEPTOR_FILE, "site_id": SITE, "ligand": LIGAND_FILE,
                   "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
                   "size_x": 20.0, "size_y": 20.0, "size_z": 20.0}]).to_csv(pairlist_path(root), index=False)
    layout = ensure_engine_layout(root, "vina", layout_profile="canonical")
    text = STI_FIXTURE.read_text(encoding="utf-8")
    (layout["poses"] / f"{LEGACY_TAG}.pdbqt").write_text(_shifted_vina_pose(text, 0.0, -8.0), encoding="utf-8")
    pipeline = _pipeline(root, root / "analysis")
    normalized = pipeline._build_normalized_frame(
        "vina", {"root": layout["root"], "poses": layout["poses"], "logs": layout["logs"], "scores": layout["scores"]},
        pipeline._pair_index())
    assert set(normalized["tag"]) == {LEGACY_TAG}
    assert normalized["replicate_id"].isna().all()


def test_detector_counts_pairs_not_replicate_files(three_replicates):
    _, poses, _ = three_replicates
    files = sorted(poses.glob("*.pdbqt"))
    assert len(files) == 3
    assert _pose_pair_count(files) == 1


def test_gnina_log_resolver_keeps_the_replicate_suffix():
    mapping = load_pairlist_mapping(_write_pairlist_for_mapping())
    assert _resolve_tag_from_pairlist(f"{PAIR_TAG}__rep02", mapping) == (f"{PAIR_TAG}__rep02", True)
    assert _resolve_tag_from_pairlist(LEGACY_TAG, mapping) == (LEGACY_TAG, True)


def _write_pairlist_for_mapping() -> Path:
    import tempfile

    handle = tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False, encoding="utf-8")
    handle.write("receptor,site_id,ligand\n" f"{RECEPTOR_FILE},{SITE},{LIGAND_FILE}\n")
    handle.close()
    return Path(handle.name)


# ----------------------------------------------------------------------------- item 6: webui box default

def test_webui_box_size_is_blank_by_default_and_omitted_from_the_command():
    from webui.targets import build_argv, get_target

    field = next(f for f in get_target("prep.pairlist").fields if f.name == "default_box_size")
    assert field.default is None
    argv = build_argv("prep.pairlist", {"mode": "cocrystal_only"}, "/tmp/project")
    assert "--default-box-size" not in argv
    argv_fixed = build_argv("prep.pairlist", {"mode": "cocrystal_only", "default_box_size": 16}, "/tmp/project")
    assert argv_fixed[argv_fixed.index("--default-box-size") + 1] == "16"
