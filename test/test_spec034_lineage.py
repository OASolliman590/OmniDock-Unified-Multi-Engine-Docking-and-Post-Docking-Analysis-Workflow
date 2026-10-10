"""Spec 034 R1: Meeko SMILES/IDX lineage source (``meeko_smiles_idx_lineage_v1``).

Covers the PDBQT lineage reader, the bond-order convention shared with the Spec 031
Kekule SDF reference, the redocking reuse (R1d) and the Spec 032 atom_map derivation (R1c).

Real-data fixtures:
- ``test/fixtures/spec034/sti_vina_model1.pdbqt`` is MODEL 1 of the retained Phase 2 Vina
  output for 1IEP/STI (public-derived redocking pose).  Scores are not asserted.
- ``test/fixtures/spec032/phase2_1iep_inputs/STI_1IEP_A201.sdf`` is the crystal reference.
- The full four-pose check reads a Phase 2 project directory named by the environment
  variable ``SPEC034_PHASE2_PROJECT`` and skips with a reason when it is not set.
"""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from post_docking_analysis.atom_mapping import (
    LINEAGE_METHOD,
    LineageError,
    compare_graph_poses,
    compare_lineage_graph_poses,
    load_pdbqt_lineage_pose,
    load_sdf_graph_pose,
)
from post_docking_analysis.md_inputs import MDInputsRequest, export_md_inputs
from post_docking_analysis.redocking_validation import _compute_pose_rmsd, run_redocking_validation
from test_spec031_pose_selection import _write_models, _write_sdf
from test_spec032_md_inputs import (
    FakeChemistryBackend,
    _pdb_atom,
    _project,
    _sdf_record,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = REPO_ROOT / "test" / "fixtures"
MODEL1_PDBQT = FIXTURE_DIR / "spec034" / "sti_vina_model1.pdbqt"
CRYSTAL_SDF = FIXTURE_DIR / "spec032" / "phase2_1iep_inputs" / "STI_1IEP_A201.sdf"

# Synthetic phenol: SMILES atom 1 (O) .. 7; PDBQT serials are deliberately not in SMILES order.
PHENOL_SMILES = "Oc1ccccc1"
PHENOL_IDX_PAIRS = [(1, 3), (2, 1), (3, 6), (4, 2), (5, 5), (6, 4), (7, 7)]  # (smiles_1based, serial)
PHENOL_HYDROGEN_SERIAL = 8


def _phenol_coordinates() -> np.ndarray:
    rng = np.random.default_rng(34)
    return rng.uniform(-3.0, 3.0, size=(7, 3))  # SMILES order


def _write_phenol_pdbqt(
    path: Path,
    *,
    idx_pairs=None,
    models: int = 0,
    shift: tuple = (0.0, 0.0, 0.0),
) -> Path:
    coords = _phenol_coordinates() + np.asarray(shift, dtype=float)
    smiles_of_serial = {serial: smiles for smiles, serial in PHENOL_IDX_PAIRS}
    pairs = PHENOL_IDX_PAIRS if idx_pairs is None else idx_pairs
    atom_lines = []
    for serial in range(1, 8):
        smiles_index = smiles_of_serial[serial]
        x, y, z = coords[smiles_index - 1]
        element_type = "OA" if smiles_index == 1 else "A"
        atom_lines.append(_pdbqt_atom(serial, f"A{serial}", x, y, z, element_type))
    x, y, z = coords[0]
    atom_lines.append(_pdbqt_atom(PHENOL_HYDROGEN_SERIAL, "H", x + 0.9, y, z, "HD"))
    idx_text = " ".join(f"{smiles} {serial}" for smiles, serial in pairs)
    header = [
        "REMARK synthetic Meeko-style lineage fixture (Spec 034 unit test)",
        f"REMARK SMILES {PHENOL_SMILES}",
        f"REMARK SMILES IDX {idx_text}",
    ]
    body = header + atom_lines + ["TORSDOF 0"]
    if models:
        chunks = []
        for number in range(1, models + 1):
            chunks.append(f"MODEL {number}")
            chunks.extend(body)
            chunks.append("ENDMDL")
        text = "\n".join(chunks) + "\n"
    else:
        text = "\n".join(body) + "\n"
    path.write_text(text, encoding="utf-8")
    return path


def _pdbqt_atom(serial: int, name: str, x: float, y: float, z: float, ad_type: str) -> str:
    return (
        f"ATOM  {serial:5d} {name:<4s} LIG A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    {0.0:6.3f} {ad_type:<2s}"
    )


def _write_phenol_kekule_sdf(path: Path) -> Path:
    """Kekule reference (bond orders 1/2 only) with the same coordinates in SMILES order."""
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    mol = Chem.MolFromSmiles(PHENOL_SMILES)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for index, (x, y, z) in enumerate(_phenol_coordinates()):
        conformer.SetAtomPosition(index, Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conformer, assignId=True)
    block = Chem.MolToMolBlock(mol, kekulize=False)
    path.write_text(block + "$$$$\n", encoding="utf-8")
    return path


def _expected_phenol_mapping() -> list:
    """Pose heavy position i (PDBQT serial i) -> reference (SDF) 0-based index."""
    smiles_of_serial = {serial: smiles for smiles, serial in PHENOL_IDX_PAIRS}
    return [smiles_of_serial[serial] - 1 for serial in range(1, 8)]


# ---------------------------------------------------------------------------
# Parser on synthetic Meeko-style PDBQT
# ---------------------------------------------------------------------------


def test_synthetic_lineage_builds_heavy_atom_graph_in_pdbqt_order(tmp_path):
    pdbqt = _write_phenol_pdbqt(tmp_path / "phenol.pdbqt")
    pose = load_pdbqt_lineage_pose(pdbqt, 1)
    assert pose.lineage_source == LINEAGE_METHOD
    # PDBQT serial 1..7 are SMILES atoms 2,4,1,6,5,3,7: the O is the third heavy record.
    assert pose.elements == ("C", "C", "O", "C", "C", "C", "C")
    assert pose.coords.shape == (7, 3)
    np.testing.assert_allclose(pose.coords, _phenol_coordinates()[_expected_phenol_mapping()], atol=2e-3)
    assert pose.graph.number_of_nodes() == 7
    assert pose.graph.number_of_edges() == 7
    labels = sorted(data["order"] for _, _, data in pose.graph.edges(data=True))
    assert labels == ["1"] + ["aromatic"] * 6
    assert all(data["formal_charge"] == 0 for _, data in pose.graph.nodes(data=True))


def test_synthetic_lineage_matches_kekule_reference_under_one_convention(tmp_path):
    pdbqt = _write_phenol_pdbqt(tmp_path / "phenol.pdbqt")
    reference_sdf = _write_phenol_kekule_sdf(tmp_path / "phenol_kekule.sdf")
    pose = load_pdbqt_lineage_pose(pdbqt, 1)
    reference = load_sdf_graph_pose(reference_sdf, 1)
    # Kekule phenol: O-C single, three ring singles, three ring doubles.
    assert sorted(data["order"] for _, _, data in reference.graph.edges(data=True)) == ["1"] * 4 + ["2"] * 3

    raw = compare_graph_poses(pose, reference)
    assert raw.status == "not_comparable"
    assert raw.reason == "graphs_not_isomorphic"

    result = compare_lineage_graph_poses(pose, reference)
    assert result.status == "comparable"
    assert result.reason == "graph_isomorphism_complete"
    assert result.mapped_heavy_atoms == 7
    assert result.mapping_coverage == pytest.approx(1.0)
    # Coordinates are written to 3 decimals in the PDBQT, so the true-mapping RMSD is ~1e-3 A.
    assert result.rmsd_angstrom == pytest.approx(0.0, abs=2e-3)
    assert list(result.selected_mapping) == _expected_phenol_mapping()
    assert result.lineage_source == LINEAGE_METHOD
    assert result.reference_lineage_source == "explicit_sdf_topology"
    assert result.to_dict()["lineage_source"] == LINEAGE_METHOD


def test_missing_idx_pair_raises_incomplete_lineage_and_is_not_comparable(tmp_path):
    pairs_without_last = PHENOL_IDX_PAIRS[:-1]
    pdbqt = _write_phenol_pdbqt(tmp_path / "phenol_missing.pdbqt", idx_pairs=pairs_without_last)
    with pytest.raises(LineageError) as excinfo:
        load_pdbqt_lineage_pose(pdbqt, 1)
    assert excinfo.value.reason == "incomplete_lineage_idx"


def test_no_lineage_remarks_is_reported_as_absent_not_as_failure(tmp_path):
    pdbqt = tmp_path / "bare.pdbqt"
    pdbqt.write_text(
        "MODEL 1\nREMARK VINA RESULT: -8.500 0.000 0.000\n"
        + _pdbqt_atom(1, "C1", 0.0, 0.0, 0.0, "C")
        + "\nENDMDL\n",
        encoding="utf-8",
    )
    with pytest.raises(LineageError) as excinfo:
        load_pdbqt_lineage_pose(pdbqt, 1)
    assert excinfo.value.reason == "no_lineage_remarks"


def test_model_index_selects_the_requested_vina_model(tmp_path):
    pdbqt = _write_phenol_pdbqt(tmp_path / "phenol_models.pdbqt", models=2)
    second = load_pdbqt_lineage_pose(pdbqt, 2)
    first = load_pdbqt_lineage_pose(pdbqt, 1)
    assert second.topology_pose_index == 2
    assert first.graph.number_of_nodes() == second.graph.number_of_nodes() == 7
    with pytest.raises(LineageError) as excinfo:
        load_pdbqt_lineage_pose(pdbqt, 3)
    assert excinfo.value.reason.startswith("model_not_found")


# ---------------------------------------------------------------------------
# Real Phase 2 data: committed MODEL 1 fixture and the crystal SDF
# ---------------------------------------------------------------------------


def test_committed_vina_model1_gives_complete_lineage_mapping_against_crystal_sdf(tmp_path):
    pdbqt = tmp_path / MODEL1_PDBQT.name
    reference = tmp_path / CRYSTAL_SDF.name
    shutil.copy2(MODEL1_PDBQT, pdbqt)
    shutil.copy2(CRYSTAL_SDF, reference)

    pose = load_pdbqt_lineage_pose(pdbqt, 1)
    assert pose.graph.number_of_nodes() == 37
    result = compare_lineage_graph_poses(pose, load_sdf_graph_pose(reference, 1))

    assert result.status != "not_comparable"
    assert result.status == "comparable"
    assert result.mapped_heavy_atoms == 37
    assert result.total_heavy_atoms_a == result.total_heavy_atoms_b == 37
    assert result.mapping_coverage == pytest.approx(1.0)
    assert result.lineage_source == LINEAGE_METHOD
    assert result.valid_mapping_count >= 1
    assert result.rmsd_angstrom is not None and np.isfinite(result.rmsd_angstrom)


def test_redocking_uses_lineage_pose_for_pdbqt_with_remarks(tmp_path):
    pdbqt = tmp_path / MODEL1_PDBQT.name
    reference = tmp_path / CRYSTAL_SDF.name
    shutil.copy2(MODEL1_PDBQT, pdbqt)
    shutil.copy2(CRYSTAL_SDF, reference)

    rmsd, error, details = _compute_pose_rmsd(pdbqt, reference, docked_pose_index=1, reference_pose_index=1)

    assert error == ""
    assert rmsd is not None and np.isfinite(rmsd)
    assert details["status"] == "comparable"
    assert details["lineage_source"] == LINEAGE_METHOD


def test_redocking_reports_incomplete_lineage_as_not_comparable(tmp_path):
    pdbqt = _write_phenol_pdbqt(tmp_path / "phenol_missing.pdbqt", idx_pairs=PHENOL_IDX_PAIRS[:-1])
    reference = _write_phenol_kekule_sdf(tmp_path / "phenol_kekule.sdf")

    rmsd, error, details = _compute_pose_rmsd(pdbqt, reference, docked_pose_index=1, reference_pose_index=1)

    assert rmsd is None
    assert "incomplete_lineage_idx" in error
    assert details["status"] == "not_comparable"


@pytest.mark.skipif(
    not os.environ.get("SPEC034_PHASE2_PROJECT"),
    reason="Phase 2 scratch project not available; set SPEC034_PHASE2_PROJECT to its directory",
)
def test_phase2_scratch_all_four_vina_poses_map_completely(tmp_path):
    project = Path(os.environ["SPEC034_PHASE2_PROJECT"])
    poses = sorted((project / "engines" / "vina" / "poses").glob("*.pdbqt"))
    if not poses:
        pytest.skip(f"no Vina pose PDBQT under {project}")
    pdbqt = tmp_path / poses[0].name
    reference = tmp_path / CRYSTAL_SDF.name
    for source in poses:
        shutil.copy2(source, tmp_path / source.name)
    shutil.copy2(CRYSTAL_SDF, reference)
    reference_pose = load_sdf_graph_pose(reference, 1)
    model_count = sum(1 for line in pdbqt.read_text(encoding="utf-8").splitlines() if line.startswith("MODEL"))
    assert model_count == 4
    for model in range(1, model_count + 1):
        result = compare_lineage_graph_poses(load_pdbqt_lineage_pose(pdbqt, model), reference_pose)
        assert result.status == "comparable", (model, result.reason)
        assert result.mapping_coverage == pytest.approx(1.0)
        assert result.lineage_source == LINEAGE_METHOD


# ---------------------------------------------------------------------------
# Spec 032 atom_map derivation (R1c)
# ---------------------------------------------------------------------------


class RecordingBackend(FakeChemistryBackend):
    """Fake backend that records the pose_to_topology permutation it was given."""

    def __init__(self) -> None:
        super().__init__()
        self.pose_to_topology_calls = []

    def prepare(self, source_file, **kwargs):
        self.pose_to_topology_calls.append(kwargs.get("pose_to_topology"))
        return super().prepare(source_file, **kwargs)


def _co_pose_pdbqt(
    path: Path,
    *,
    idx_line: str | None = "REMARK SMILES IDX 1 2 2 1",
    with_smiles: bool = True,
) -> Path:
    # SMILES "CO": atom 1 is C (PDBQT serial 2), atom 2 is O (PDBQT serial 1).
    lines = ["MODEL 1", "REMARK VINA RESULT: -8.500 0.000 0.000"]
    if with_smiles:
        lines.append("REMARK SMILES CO")
    if idx_line is not None:
        lines.append(idx_line)
    lines += [
        _pdb_atom(1, "O1", "LIG", "B", 1, 21.2, 0.0, 0.0, "O", record="HETATM") + "  0.000 OA",
        _pdb_atom(2, "C1", "LIG", "B", 1, 20.0, 0.0, 0.0, "C", record="HETATM") + "  0.000 C",
        "ENDMDL",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _vina_lineage_request(
    tmp_path: Path,
    *,
    idx_line: str | None,
    atom_map: str,
    with_smiles: bool = True,
):
    project, _, tags, _, receptor_map = _project(tmp_path)
    pose = project / "3-Docking" / "vina" / "poses" / "prot_site_lig.pdbqt"
    pose.parent.mkdir(parents=True)
    _co_pose_pdbqt(pose, idx_line=idx_line, with_smiles=with_smiles)
    selection = project / "4-Working" / "scores" / "unified" / "best_pose_per_tag_by_engine.csv"
    frame = pd.read_csv(selection)
    frame.loc[0, ["engine", "pose", "pose_file", "affinity_kcal_mol"]] = ["vina", 1, str(pose), -8.5]
    frame.to_csv(selection, index=False)
    topology = project / "topologies" / "co.sdf"
    topology.parent.mkdir(parents=True)
    topology.write_text(_sdf_record("topology", 0.0, 0.0), encoding="utf-8")
    topology_map = project / "topology_map.csv"
    pd.DataFrame(
        [{"tag": "prot_site_lig", "topology_file": str(topology), "topology_pose": 1, "atom_map": atom_map}]
    ).to_csv(topology_map, index=False)
    request = MDInputsRequest.from_dict(
        {
            "project_dir": str(project),
            "engine": "vina",
            "tags_file": str(tags),
            "receptor_map": str(receptor_map),
            "topology_map": str(topology_map),
            "pH": 7.4,
            "protonation_policy": "openbabel_predicted",
        }
    )
    return request


def test_vina_row_with_empty_atom_map_derives_permutation_from_lineage(tmp_path):
    request = _vina_lineage_request(tmp_path, idx_line="REMARK SMILES IDX 1 2 2 1", atom_map="")
    backend = RecordingBackend()
    result = export_md_inputs(request, chemistry_backend=backend)

    row = result["rows"][0]
    assert row["status"] == "completed", row
    provenance = json.loads(Path(row["provenance_file"]).read_text(encoding="utf-8"))
    assert provenance["gates"]["G3"]["atom_map"] == [1, 0]
    assert provenance["gates"]["G3"]["atom_map_source"] == LINEAGE_METHOD
    assert provenance["atom_map_source"] == LINEAGE_METHOD
    assert backend.pose_to_topology_calls == [[1, 0]]


def test_explicit_atom_map_wins_over_lineage(tmp_path):
    request = _vina_lineage_request(tmp_path, idx_line="REMARK SMILES IDX 1 2 2 1", atom_map="[0, 1]")
    backend = RecordingBackend()
    result = export_md_inputs(request, chemistry_backend=backend)

    row = result["rows"][0]
    assert row["status"] == "completed", row
    provenance = json.loads(Path(row["provenance_file"]).read_text(encoding="utf-8"))
    assert provenance["gates"]["G3"]["atom_map"] == [0, 1]
    assert provenance["gates"]["G3"]["atom_map_source"] == "explicit_topology_map"
    assert provenance["atom_map_source"] == "explicit_topology_map"
    assert backend.pose_to_topology_calls == [[0, 1]]


def test_empty_atom_map_without_lineage_stays_not_comparable(tmp_path):
    request = _vina_lineage_request(tmp_path, idx_line=None, atom_map="", with_smiles=False)
    result = export_md_inputs(request, chemistry_backend=RecordingBackend())

    row = result["rows"][0]
    assert row["status"] == "not_comparable"
    assert row["reason"] == "missing_pose_to_topology_atom_map:no_lineage_remarks"
    assert row["gates"]["G3"]["status"] == "not_comparable"


def test_empty_atom_map_with_incomplete_lineage_stays_not_comparable(tmp_path):
    request = _vina_lineage_request(tmp_path, idx_line="REMARK SMILES IDX 1 2", atom_map="")
    result = export_md_inputs(request, chemistry_backend=RecordingBackend())

    row = result["rows"][0]
    assert row["status"] == "not_comparable"
    assert row["reason"] == "incomplete_lineage_idx"
    assert row["gates"]["G3"]["status"] == "not_comparable"


# ---------------------------------------------------------------------------
# Spec 034 T002a: redocking primary is in-place heavy-atom RMSD (in_place_v1)
# ---------------------------------------------------------------------------


def test_rigid_5A_translation_in_place_is_5_kabsch_is_0_and_decision_fails(tmp_path):
    # Explicit-topology path (no lineage remarks): three carbons, path graph C2-C1-C3.
    reference_coords = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0)]
    shifted_coords = [(x + 5.0, y, z) for x, y, z in reference_coords]
    docked_pose = tmp_path / "docked_shifted.pdbqt"
    reference_pose = tmp_path / "reference_pose.sdf"
    topology = tmp_path / "topology.sdf"
    _write_models(docked_pose, [shifted_coords])
    _write_sdf(reference_pose, reference_coords)
    _write_sdf(topology, reference_coords)

    rmsd, error, details = _compute_pose_rmsd(
        docked_pose,
        reference_pose,
        docked_pose_index=1,
        reference_pose_index=1,
        docked_topology_file=topology,
    )
    assert error == ""
    assert rmsd == pytest.approx(5.0, abs=1e-6)
    assert details["rmsd_frame"] == "in_place_v1"
    assert details["kabsch_rmsd_angstrom"] == pytest.approx(0.0, abs=1e-6)

    outputs = run_redocking_validation(
        project_dir=tmp_path,
        best_by_engine=pd.DataFrame(
            [
                {
                    "protein": "P1",
                    "ligand": "REF",
                    "tag": "P1_site_1_REF",
                    "engine": "vina",
                    "pose": 1,
                    "pose_file": str(docked_pose),
                    "topology_file": str(topology),
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
    row = validation.iloc[0]
    assert float(row["redocking_rmsd_angstrom"]) == pytest.approx(5.0, abs=1e-6)
    assert row["redocking_classification"] == "fail"
    assert row["redocking_rmsd_frame"] == "in_place_v1"
    assert float(row["kabsch_secondary_rmsd_angstrom"]) == pytest.approx(0.0, abs=1e-6)
    assert row["rmsd_secondary_method"] == "kabsch_secondary"
    summary = Path(outputs["summary_file"]).read_text(encoding="utf-8")
    assert "in_place_v1" in summary
    assert "decision_rmsd" in summary and "kabsch_secondary" in summary


def test_rigid_5A_translation_lineage_path_in_place_is_5_kabsch_is_0(tmp_path):
    shifted = _write_phenol_pdbqt(tmp_path / "phenol_shifted.pdbqt", shift=(5.0, 0.0, 0.0))
    reference = _write_phenol_kekule_sdf(tmp_path / "phenol_kekule.sdf")

    rmsd, error, details = _compute_pose_rmsd(shifted, reference, docked_pose_index=1, reference_pose_index=1)

    assert error == ""
    assert rmsd == pytest.approx(5.0, abs=2e-3)
    assert details["lineage_source"] == LINEAGE_METHOD
    assert details["rmsd_frame"] == "in_place_v1"
    assert details["kabsch_rmsd_angstrom"] == pytest.approx(0.0, abs=2e-3)


def test_committed_sti_model1_in_place_rmsd_against_crystal(tmp_path):
    pdbqt = tmp_path / MODEL1_PDBQT.name
    reference = tmp_path / CRYSTAL_SDF.name
    shutil.copy2(MODEL1_PDBQT, pdbqt)
    shutil.copy2(CRYSTAL_SDF, reference)

    rmsd, error, details = _compute_pose_rmsd(pdbqt, reference, docked_pose_index=1, reference_pose_index=1)

    assert error == ""
    assert details["status"] == "comparable"
    assert details["rmsd_frame"] == "in_place_v1"
    assert rmsd == pytest.approx(0.857, abs=0.005)
    assert details["kabsch_rmsd_angstrom"] == pytest.approx(1.116, abs=0.005)
