"""Parser-only pose record contracts that do not require RDKit."""

from pathlib import Path

import pytest

from post_docking_analysis.pose_geometry import PoseGeometryError, selected_record


def test_sdf_empty_record_keeps_ordinal_and_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "poses.sdf"
    path.write_text("$$$$\nsecond record\n$$$$\n", encoding="utf-8")

    with pytest.raises(PoseGeometryError, match="empty_pose_record"):
        selected_record(path, 1)
    assert selected_record(path, 2).strip() == "second record"
    with pytest.raises(PoseGeometryError, match="pose_index_out_of_range"):
        selected_record(path, 3)


def test_sdf_interior_empty_record_does_not_shift_later_pose(tmp_path: Path) -> None:
    path = tmp_path / "poses.sdf"
    path.write_text("first record\n$$$$\n$$$$\nthird record\n$$$$\n", encoding="utf-8")

    assert selected_record(path, 1).strip() == "first record"
    with pytest.raises(PoseGeometryError, match="empty_pose_record"):
        selected_record(path, 2)
    assert selected_record(path, 3).strip() == "third record"


def test_sdf_separator_must_be_a_standalone_line(tmp_path: Path) -> None:
    path = tmp_path / "poses.sdf"
    path.write_text("title contains $$$$ as data\n", encoding="utf-8")

    assert selected_record(path, 1) == "title contains $$$$ as data\n"


def test_sdf_unterminated_nonempty_tail_keeps_its_ordinal(tmp_path: Path) -> None:
    path = tmp_path / "poses.sdf"
    path.write_text("first\n$$$$\nunterminated", encoding="utf-8")

    assert selected_record(path, 2).strip() == "unterminated"


def test_model_records_are_isolated_and_global_header_is_not_a_model(tmp_path: Path) -> None:
    path = tmp_path / "poses.pdbqt"
    path.write_text(
        "REMARK GLOBAL HEADER\n"
        "MODEL 1\n"
        "ATOM model-one\n"
        "ENDMDL\n"
        "MODEL 2\n"
        "ATOM model-two\n"
        "ENDMDL\nEND\n",
        encoding="utf-8",
    )

    assert selected_record(path, 1).splitlines() == ["ATOM model-one"]
    assert selected_record(path, 2).splitlines() == ["ATOM model-two"]


@pytest.mark.parametrize(
    ("text", "reason"),
    [
        ("MODEL 1\nMODEL 2\nENDMDL\n", "nested_model"),
        ("MODEL 1\nATOM x\n", "unterminated_model"),
        ("ENDMDL\nMODEL 1\nENDMDL\n", "orphan_endmdl"),
        ("ATOM outside\nMODEL 1\nENDMDL\n", "pdbqt_content_outside_model"),
        ("MODEL 1\nENDMDL\nREMARK after model\n", "pdbqt_content_outside_model"),
        ("MODEL 2\nENDMDL\n", "nonsequential_model_identifier"),
        ("MODEL -1\nENDMDL\n", "invalid_model_identifier"),
        ("MODEL x\nENDMDL\n", "invalid_model_identifier"),
    ],
)
def test_model_boundary_errors_are_fail_closed(tmp_path: Path, text: str, reason: str) -> None:
    path = tmp_path / "malformed.pdbqt"
    path.write_text(text, encoding="utf-8")

    with pytest.raises(PoseGeometryError, match=reason):
        selected_record(path, 1)
