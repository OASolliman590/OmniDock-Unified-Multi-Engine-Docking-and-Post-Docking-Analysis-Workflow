from pathlib import Path

import pandas as pd

from post_docking_analysis.redocking_validation import run_redocking_validation
from post_docking_analysis.reference_policy import (
    LEGACY_METHOD,
    REFERENCE_METHOD,
    classify_reference_frame,
    classify_reference_row,
    reference_mask,
)


def _write_sdf(path: Path, atoms, bonds) -> None:
    lines = ["ref", "  DockForge", "", f"{len(atoms):>3}{len(bonds):>3}  0  0  0  0            999 V2000"]
    for x, y, z, element in atoms:
        lines.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {element:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
    for left, right, order in bonds:
        lines.append(f"{left:>3}{right:>3}{order:>3}  0  0  0  0")
    lines.extend(["M  END", "$$$$", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def test_explicit_flag_requires_identity_or_lineage() -> None:
    incomplete = classify_reference_row({"is_cocrystal_benchmark": True, "ligand": "anything"})
    assert incomplete.classification == "unclassified"
    assert not incomplete.is_reference
    explicit = classify_reference_row({
        "is_cocrystal_benchmark": True,
        "cocrystal_ligand_name": "ATP_A_101",
    })
    assert explicit.is_reference
    assert explicit.method == REFERENCE_METHOD
    assert not explicit.validation_anchor_eligible


def test_names_and_false_like_values_never_create_scientific_reference() -> None:
    frame = pd.DataFrame([
        {"ligand": "novel_ligand_1", "tag": "reference_fake", "is_cocrystal_benchmark": False},
        {"ligand": "1ABC_ligand_ATP_A_1", "is_cocrystal_benchmark": "False"},
        {"ligand": "native_control", "is_cocrystal_benchmark": 0},
        {"ligand": "benchmark", "is_cocrystal_benchmark": None},
    ])
    assert not reference_mask(frame).any()


def test_legacy_mode_emits_candidate_but_never_anchor() -> None:
    result = classify_reference_row({"tag": "reference_old"}, legacy_mode=True)
    assert result.method == LEGACY_METHOD
    assert result.legacy_candidate
    assert not result.is_reference
    assert not result.validation_anchor_eligible


def test_frame_outputs_are_versioned_and_distinguishable() -> None:
    frame = classify_reference_frame(pd.DataFrame([
        {"is_cocrystal_benchmark": True, "cocrystal_ligand_name": "ATP"},
        {"tag": "reference_old"},
    ]), legacy_mode=True)
    assert list(frame["reference_classification_method"]) == [REFERENCE_METHOD, LEGACY_METHOD]
    assert list(frame["legacy_reference_candidate"]) == [False, True]


def test_redocking_uses_graph_mapping_and_explicit_source_correspondence(tmp_path: Path) -> None:
    docked = tmp_path / "docked.sdf"
    reference = tmp_path / "reference.sdf"
    _write_sdf(docked, [(0, 0, 0, "C"), (1, 0, 0, "O"), (0, 1, 0, "N")], [(1, 2, 1), (1, 3, 1)])
    _write_sdf(reference, [(0, 1, 0, "N"), (0, 0, 0, "C"), (1, 0, 0, "O")], [(2, 3, 1), (2, 1, 1)])
    rows = pd.DataFrame([{
        "protein": "1ABC", "pdb_id": "1ABC", "ligand": "ATP", "tag": "1ABC__ATP",
        "engine": "vina", "pose": 1, "affinity_kcal_mol": -7.0,
        "pose_file": str(docked), "reference_pose_file": str(reference),
        "is_cocrystal_benchmark": True, "cocrystal_ligand_name": "ATP",
    }])
    bundle = run_redocking_validation(project_dir=tmp_path, best_by_engine=rows, output_dir=tmp_path / "reports")
    row = bundle["validation_df"].iloc[0]
    assert row["redocking_classification"] == "pass"
    assert row["mapping_method"] == "spec031-atom-mapping-v1"
    assert bundle["reference_baselines_df"].iloc[0]["anchor_status"] == "eligible"

