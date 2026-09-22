"""Reference-gate input contract tests independent of chemistry backends."""

from pathlib import Path

import pandas as pd
import pytest

from post_docking_analysis.redocking_validation import run_redocking_validation


def _reference_row(*, pose: int, pose_file: str) -> dict[str, object]:
    return {
        "engine": "vina",
        "tag": "P_reference_L",
        "protein": "P",
        "ligand": "L",
        "site_id": "reference",
        "pose": pose,
        "pose_file": pose_file,
        "affinity_kcal_mol": -8.0,
        "is_cocrystal_benchmark": True,
        "reference_pose_file": "reference.sdf",
        "reference_source": "cocrystal:1ABC",
        "reference_pdb_id": "1ABC",
        "reference_frame_id": "frame",
        "receptor_frame_id": "frame",
    }


def test_duplicate_reference_rows_are_rejected_before_denominator_changes(tmp_path: Path) -> None:
    first = _reference_row(pose=1, pose_file="missing-first.sdf")
    second = _reference_row(pose=2, pose_file="missing-second.sdf")

    with pytest.raises(ValueError, match="duplicate_reference_input_rows"):
        run_redocking_validation(
            project_dir=tmp_path,
            best_by_engine=pd.DataFrame([first, second]),
            output_dir=tmp_path / "validation",
        )
