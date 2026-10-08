from pathlib import Path

import numpy as np
import pandas as pd

from post_docking_analysis.geometric_consensus import _extract_pdb_like_pose
from post_docking_analysis.redocking_validation import run_redocking_validation


def _atom_line(index: int, x: float, y: float, z: float, element: str = "C") -> str:
    return (
        f"HETATM{index:5d} {element}{index:<3d} LIG A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
    )


def _write_models(
    path: Path,
    models: list[list[tuple[float, float, float]]],
    *,
    close_final_model: bool = True,
) -> None:
    lines: list[str] = []
    for model_index, coords in enumerate(models, start=1):
        lines.append(f"MODEL     {model_index:4d}")
        lines.extend(
            _atom_line(atom_index, x, y, z)
            for atom_index, (x, y, z) in enumerate(coords, start=1)
        )
        if model_index < len(models) or close_final_model:
            lines.append("ENDMDL")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_sdf(path: Path, coords: list[tuple[float, float, float]]) -> None:
    lines = ["LIG", "  DockForge", "", f"{len(coords):>3}{2:>3}  0  0  0  0            999 V2000"]
    lines.extend(
        f"{x:10.4f}{y:10.4f}{z:10.4f} C   0  0  0  0  0  0  0  0  0  0  0  0"
        for x, y, z in coords
    )
    lines.extend(["  1  2  1  0  0  0  0", "  1  3  1  0  0  0  0", "M  END", "$$$$", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def test_extract_pdbqt_pose_uses_one_based_model_ordinal(tmp_path: Path) -> None:
    pose_file = tmp_path / "three_models.pdbqt"
    models = [
        [(1.0, 0.0, 0.0), (1.0, 1.0, 0.0)],
        [(2.0, 0.0, 0.0), (2.0, 2.0, 0.0)],
        [(3.0, 0.0, 0.0), (3.0, 3.0, 0.0)],
    ]
    _write_models(pose_file, models)

    coords, elements, error = _extract_pdb_like_pose(pose_file, pose_index=2)

    assert error == ""
    assert elements == ["C", "C"]
    np.testing.assert_allclose(coords, np.asarray(models[1], dtype=float))


def test_extract_pdbqt_pose_accepts_final_model_closed_by_eof(tmp_path: Path) -> None:
    pose_file = tmp_path / "unterminated_final_model.pdbqt"
    final_coords = [(4.0, 0.0, 0.0), (4.0, 2.0, 0.0)]
    _write_models(
        pose_file,
        [[(1.0, 0.0, 0.0), (1.0, 1.0, 0.0)], final_coords],
        close_final_model=False,
    )

    coords, _, error = _extract_pdb_like_pose(pose_file, pose_index=2)

    assert error == ""
    np.testing.assert_allclose(coords, np.asarray(final_coords, dtype=float))


def test_extract_pdbqt_pose_rejects_out_of_range_without_fallback(tmp_path: Path) -> None:
    pose_file = tmp_path / "two_models.pdbqt"
    _write_models(
        pose_file,
        [[(1.0, 0.0, 0.0)], [(2.0, 0.0, 0.0)]],
    )

    coords, elements, error = _extract_pdb_like_pose(pose_file, pose_index=3)

    assert coords.shape == (0, 3)
    assert elements == []
    assert error == "pose_index_out_of_range:3>2"


def test_extract_pdbqt_pose_rejects_non_positive_or_fractional_indices(tmp_path: Path) -> None:
    pose_file = tmp_path / "one_model.pdbqt"
    _write_models(pose_file, [[(1.0, 0.0, 0.0)]])

    for invalid_index in (0, -1, 1.5):
        coords, elements, error = _extract_pdb_like_pose(pose_file, pose_index=invalid_index)
        assert coords.shape == (0, 3)
        assert elements == []
        assert error == "invalid_pose_index"


def test_extract_single_pose_file_rejects_second_pose(tmp_path: Path) -> None:
    pose_file = tmp_path / "single_pose.pdb"
    pose_file.write_text(_atom_line(1, 1.0, 2.0, 3.0) + "\n", encoding="utf-8")

    first_coords, _, first_error = _extract_pdb_like_pose(pose_file, pose_index=1)
    second_coords, _, second_error = _extract_pdb_like_pose(pose_file, pose_index=2)

    assert first_error == ""
    np.testing.assert_allclose(first_coords, np.asarray([[1.0, 2.0, 3.0]]))
    assert second_coords.shape == (0, 3)
    assert second_error == "pose_index_out_of_range:2>1"


def test_extract_pdbqt_pose_reports_malformed_model_boundaries(tmp_path: Path) -> None:
    nested_file = tmp_path / "nested_model.pdbqt"
    nested_file.write_text(
        "MODEL 1\n" + _atom_line(1, 0.0, 0.0, 0.0) + "\nMODEL 2\n",
        encoding="utf-8",
    )
    unexpected_end_file = tmp_path / "unexpected_end.pdbqt"
    unexpected_end_file.write_text("ENDMDL\n", encoding="utf-8")
    mixed_file = tmp_path / "atom_before_model.pdbqt"
    mixed_file.write_text(
        _atom_line(1, 0.0, 0.0, 0.0) + "\nMODEL 1\nENDMDL\n",
        encoding="utf-8",
    )

    nested_coords, _, nested_error = _extract_pdb_like_pose(nested_file, pose_index=1)
    end_coords, _, end_error = _extract_pdb_like_pose(unexpected_end_file, pose_index=1)
    mixed_coords, _, mixed_error = _extract_pdb_like_pose(mixed_file, pose_index=1)

    assert nested_coords.shape == (0, 3)
    assert nested_error == "malformed_model:nested_model"
    assert end_coords.shape == (0, 3)
    assert end_error == "malformed_model:unexpected_endmdl"
    assert mixed_coords.shape == (0, 3)
    assert mixed_error == "malformed_model:atom_outside_model"


def test_redocking_validation_uses_selected_best_row_pose(tmp_path: Path) -> None:
    docked_pose = tmp_path / "docked_multi_pose.pdbqt"
    reference_pose = tmp_path / "reference_pose.sdf"
    docked_topology = tmp_path / "docked_topology.sdf"
    reference_coords = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0)]
    different_coords = [(0.0, 0.0, 0.0), (8.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    _write_models(docked_pose, [different_coords, reference_coords])
    _write_sdf(reference_pose, reference_coords)
    _write_sdf(docked_topology, reference_coords)

    outputs = run_redocking_validation(
        project_dir=tmp_path,
        best_by_engine=pd.DataFrame(
            [
                {
                    "protein": "P1",
                    "ligand": "REF",
                    "tag": "P1_site_1_REF",
                    "engine": "vina",
                    "pose": 2,
                    "pose_file": str(docked_pose),
                    "topology_file": str(docked_topology),
                    "reference_pose_file": str(reference_pose),
                    "affinity_kcal_mol": -8.0,
                    "is_cocrystal_benchmark": True,
                    "cocrystal_ligand_name": "REF",
                    "pdb_id": "P1",
                }
            ]
        ),
        output_dir=tmp_path / "reports",
    )

    validation = outputs["validation_df"]
    assert len(validation) == 1
    assert int(validation.iloc[0]["docked_pose_index"]) == 2
    assert validation.iloc[0]["redocking_classification"] == "pass"
    assert float(validation.iloc[0]["redocking_rmsd_angstrom"]) < 1e-9
