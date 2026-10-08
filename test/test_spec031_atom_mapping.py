from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from post_docking_analysis.atom_mapping import (
    MAPPING_METHOD,
    attach_coordinates_to_topology,
    backend_capability,
    compare_graph_poses,
    load_sdf_graph_pose,
)
from post_docking_analysis.geometric_consensus import compute_geometric_consensus_v2


def _write_sdf(path: Path, atoms: list[tuple[float, float, float, str]], bonds: list[tuple[int, int, int]]) -> None:
    lines = [
        "spec031",
        "  DockForge",
        "",
        f"{len(atoms):>3}{len(bonds):>3}  0  0  0  0            999 V2000",
    ]
    for x, y, z, element in atoms:
        lines.append(
            f"{x:10.4f}{y:10.4f}{z:10.4f} {element:<3} 0  0  0  0  0  0  0  0  0  0  0  0"
        )
    for atom_a, atom_b, order in bonds:
        lines.append(f"{atom_a:>3}{atom_b:>3}{order:>3}  0  0  0  0")
    lines.extend(["M  END", "$$$$", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def test_backend_capability_is_explicit() -> None:
    capability = backend_capability()
    assert capability["status"] == "completed"
    assert capability["mapping_method"] == MAPPING_METHOD
    assert capability["supports_edge_labels"] is True


def test_permuted_atoms_map_by_bond_graph(tmp_path: Path) -> None:
    first = tmp_path / "first.sdf"
    second = tmp_path / "second.sdf"
    atoms = [(0.0, 0.0, 0.0, "C"), (1.2, 0.0, 0.0, "N"), (2.2, 0.7, 0.0, "O")]
    _write_sdf(first, atoms, [(1, 2, 1), (2, 3, 2)])
    permuted = [atoms[2], atoms[0], atoms[1]]
    _write_sdf(second, permuted, [(2, 3, 1), (3, 1, 2)])

    result = compare_graph_poses(load_sdf_graph_pose(first), load_sdf_graph_pose(second))

    assert result.comparable
    assert result.rmsd_angstrom is not None and result.rmsd_angstrom < 1e-8
    assert result.mapped_heavy_atoms == 3
    assert result.mapping_coverage == 1.0


def test_symmetry_uses_minimum_valid_mapping(tmp_path: Path) -> None:
    reference = tmp_path / "reference.sdf"
    pose = tmp_path / "pose.sdf"
    _write_sdf(
        reference,
        [(0.0, 0.0, 0.0, "C"), (-1.0, 0.0, 0.0, "O"), (1.0, 0.0, 0.0, "O")],
        [(1, 2, 1), (1, 3, 1)],
    )
    _write_sdf(
        pose,
        [(0.0, 0.0, 0.0, "C"), (1.0, 0.0, 0.0, "O"), (-1.0, 0.0, 0.0, "O")],
        [(1, 2, 1), (1, 3, 1)],
    )

    result = compare_graph_poses(load_sdf_graph_pose(reference), load_sdf_graph_pose(pose))

    assert result.comparable
    assert result.valid_mapping_count == 2
    assert result.rmsd_angstrom is not None and result.rmsd_angstrom < 1e-8


def test_same_formula_constitutional_isomer_is_not_comparable(tmp_path: Path) -> None:
    chain = tmp_path / "chain.sdf"
    branched = tmp_path / "branched.sdf"
    atoms = [
        (0.0, 0.0, 0.0, "C"),
        (1.0, 0.0, 0.0, "C"),
        (2.0, 0.0, 0.0, "C"),
        (3.0, 0.0, 0.0, "C"),
    ]
    _write_sdf(chain, atoms, [(1, 2, 1), (2, 3, 1), (3, 4, 1)])
    _write_sdf(branched, atoms, [(1, 2, 1), (1, 3, 1), (1, 4, 1)])

    result = compare_graph_poses(load_sdf_graph_pose(chain), load_sdf_graph_pose(branched))

    assert result.status == "not_comparable"
    assert result.reason == "graphs_not_isomorphic"
    assert result.rmsd_angstrom is None


def test_atom_count_mismatch_has_no_numeric_rmsd(tmp_path: Path) -> None:
    first = tmp_path / "first.sdf"
    second = tmp_path / "second.sdf"
    _write_sdf(first, [(0.0, 0.0, 0.0, "C"), (1.0, 0.0, 0.0, "O")], [(1, 2, 1)])
    _write_sdf(
        second,
        [(0.0, 0.0, 0.0, "C"), (1.0, 0.0, 0.0, "O"), (2.0, 0.0, 0.0, "N")],
        [(1, 2, 1), (2, 3, 1)],
    )

    result = compare_graph_poses(load_sdf_graph_pose(first), load_sdf_graph_pose(second))

    assert result.status == "not_comparable"
    assert result.reason == "atom_count_mismatch"
    assert result.rmsd_angstrom is None


def test_explicit_topology_can_attach_docked_coordinates(tmp_path: Path) -> None:
    topology = tmp_path / "ligand.sdf"
    atoms = [(0.0, 0.0, 0.0, "C"), (1.0, 0.0, 0.0, "N")]
    _write_sdf(topology, atoms, [(1, 2, 1)])
    coords = np.asarray([[5.0, 2.0, 1.0], [6.0, 2.0, 1.0]], dtype=float)

    attached = attach_coordinates_to_topology(coords, ["C", "N"], topology)
    result = compare_graph_poses(load_sdf_graph_pose(topology), attached)

    assert result.comparable
    assert result.rmsd_angstrom is not None and result.rmsd_angstrom < 1e-8


def test_bond_order_mismatch_is_not_comparable(tmp_path: Path) -> None:
    single = tmp_path / "single.sdf"
    double = tmp_path / "double.sdf"
    atoms = [(0.0, 0.0, 0.0, "C"), (1.2, 0.0, 0.0, "O")]
    _write_sdf(single, atoms, [(1, 2, 1)])
    _write_sdf(double, atoms, [(1, 2, 2)])

    result = compare_graph_poses(load_sdf_graph_pose(single), load_sdf_graph_pose(double))

    assert result.status == "not_comparable"
    assert result.reason == "graphs_not_isomorphic"


def test_v2_geometric_consensus_never_soft_aligns_incompatible_graphs(tmp_path: Path) -> None:
    first = tmp_path / "first.sdf"
    second = tmp_path / "second.sdf"
    atoms = [(0.0, 0.0, 0.0, "C"), (1.2, 0.0, 0.0, "O")]
    _write_sdf(first, atoms, [(1, 2, 1)])
    _write_sdf(second, atoms, [(1, 2, 2)])
    frame = pd.DataFrame([
        {"engine": "vina", "protein": "p", "tag": "t", "pose_file": str(first), "pose": 1},
        {"engine": "smina", "protein": "p", "tag": "t", "pose_file": str(second), "pose": 1},
    ])
    result = compute_geometric_consensus_v2(frame).iloc[0]
    assert result["geometry_status"] == "not_comparable"
    assert pd.isna(result["geometric_pairwise_rmsd_mean"])
    assert "soft_alignment" not in result["geometric_pairwise_rmsd_json"]
