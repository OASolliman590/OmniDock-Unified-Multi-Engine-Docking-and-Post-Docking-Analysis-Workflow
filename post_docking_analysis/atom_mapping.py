"""Scientifically conservative graph/symmetry-aware heavy-atom RMSD.

Spec 031 method ``spec031-atom-mapping-v1`` requires complete heavy-atom
graph isomorphism.  This module never falls back to element sorting,
nearest-neighbour assignment, or common-substructure truncation.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np

try:
    import networkx as nx
except ImportError:  # pragma: no cover - exercised through capability status
    nx = None  # type: ignore[assignment]


MAPPING_METHOD = "spec031-atom-mapping-v1"
ALIGNMENT_METHOD = "kabsch"
DEFAULT_MAX_ISOMORPHISMS = 100_000


_SDF_CHARGE_CODES = {
    0: 0,
    1: 3,
    2: 2,
    3: 1,
    5: -1,
    6: -2,
    7: -3,
}


@dataclass
class GraphPose:
    coords: np.ndarray
    elements: Tuple[str, ...]
    graph: object
    topology_path: str
    topology_sha256: str
    topology_pose_index: int


@dataclass
class AtomMappingResult:
    status: str
    reason: str
    rmsd_angstrom: Optional[float] = None
    mapping_method: str = MAPPING_METHOD
    alignment_method: str = ALIGNMENT_METHOD
    mapped_heavy_atoms: int = 0
    total_heavy_atoms_a: int = 0
    total_heavy_atoms_b: int = 0
    mapping_coverage: float = 0.0
    valid_mapping_count: int = 0
    selected_mapping: Optional[Tuple[int, ...]] = None
    backend: str = "networkx"
    backend_version: str = ""
    topology_a: str = ""
    topology_b: str = ""
    topology_sha256_a: str = ""
    topology_sha256_b: str = ""

    @property
    def comparable(self) -> bool:
        return self.status == "comparable" and self.rmsd_angstrom is not None

    def to_dict(self) -> Dict[str, object]:
        return {
            "status": self.status,
            "reason": self.reason,
            "rmsd_angstrom": self.rmsd_angstrom,
            "mapping_method": self.mapping_method,
            "alignment_method": self.alignment_method,
            "mapped_heavy_atoms": self.mapped_heavy_atoms,
            "total_heavy_atoms_a": self.total_heavy_atoms_a,
            "total_heavy_atoms_b": self.total_heavy_atoms_b,
            "mapping_coverage": self.mapping_coverage,
            "valid_mapping_count": self.valid_mapping_count,
            "selected_mapping": list(self.selected_mapping) if self.selected_mapping is not None else None,
            "backend": self.backend,
            "backend_version": self.backend_version,
            "topology_a": self.topology_a,
            "topology_b": self.topology_b,
            "topology_sha256_a": self.topology_sha256_a,
            "topology_sha256_b": self.topology_sha256_b,
        }


def backend_capability() -> Dict[str, object]:
    available = nx is not None
    return {
        "status": "completed" if available else "skipped_missing_dependency",
        "backend": "networkx",
        "backend_version": str(getattr(nx, "__version__", "")) if available else "",
        "mapping_method": MAPPING_METHOD,
        "supports_node_labels": bool(available),
        "supports_edge_labels": bool(available),
        "supports_isomorphism_enumeration": bool(available),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _normalize_element(value: object) -> str:
    letters = "".join(ch for ch in str(value or "").strip() if ch.isalpha())
    if not letters:
        return ""
    if len(letters) >= 2 and letters[1].islower():
        return letters[:2].upper()
    return letters[:1].upper()


def _sdf_blocks(path: Path) -> List[List[str]]:
    text = path.read_text(encoding="utf-8", errors="strict")
    return [block.splitlines() for block in text.split("$$$$") if block.strip()]


def _parse_counts(line: str) -> Tuple[int, int]:
    try:
        return int(line[0:3].strip()), int(line[3:6].strip())
    except (TypeError, ValueError):
        parts = str(line).split()
        if len(parts) < 2:
            raise ValueError("invalid_sdf_counts")
        return int(parts[0]), int(parts[1])


def _parse_sdf_atom(line: str) -> Tuple[float, float, float, str, int, str]:
    try:
        x = float(line[0:10].strip())
        y = float(line[10:20].strip())
        z = float(line[20:30].strip())
        element = _normalize_element(line[31:34].strip())
        charge_code = int((line[36:39] or "0").strip() or 0)
        stereo = str((line[39:42] or "0").strip() or "0")
    except (TypeError, ValueError):
        parts = str(line).split()
        if len(parts) < 4:
            raise ValueError("invalid_sdf_atom")
        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
        element = _normalize_element(parts[3])
        charge_code = int(parts[5]) if len(parts) > 5 and parts[5].lstrip("+-").isdigit() else 0
        stereo = parts[6] if len(parts) > 6 else "0"
    if not element:
        raise ValueError("missing_sdf_element")
    return x, y, z, element, _SDF_CHARGE_CODES.get(charge_code, 0), stereo


def _parse_sdf_bond(line: str) -> Tuple[int, int, str, str]:
    try:
        atom_a = int(line[0:3].strip()) - 1
        atom_b = int(line[3:6].strip()) - 1
        raw_order = int(line[6:9].strip())
        stereo = str((line[9:12] or "0").strip() or "0")
    except (TypeError, ValueError):
        parts = str(line).split()
        if len(parts) < 3:
            raise ValueError("invalid_sdf_bond")
        atom_a, atom_b = int(parts[0]) - 1, int(parts[1]) - 1
        raw_order = int(parts[2])
        stereo = parts[3] if len(parts) > 3 else "0"
    order = "aromatic" if raw_order == 4 else str(raw_order)
    return atom_a, atom_b, order, stereo


def load_sdf_graph_pose(path: Path, pose_index: int = 1) -> GraphPose:
    """Load one SDF block as a heavy-atom, bond-labelled graph pose."""
    if nx is None:
        raise RuntimeError("skipped_missing_dependency:networkx")
    source = Path(path).expanduser().resolve()
    if source.suffix.lower() not in {".sdf", ".mol"}:
        raise ValueError("unsupported_topology_format")
    if int(pose_index) < 1:
        raise ValueError("invalid_pose_index")
    blocks = _sdf_blocks(source)
    if int(pose_index) > len(blocks):
        raise ValueError(f"pose_index_out_of_range:{pose_index}>{len(blocks)}")
    lines = blocks[int(pose_index) - 1]
    if len(lines) < 4:
        raise ValueError("invalid_sdf_header")
    atom_count, bond_count = _parse_counts(lines[3])
    if atom_count <= 0 or len(lines) < 4 + atom_count + bond_count:
        raise ValueError("truncated_sdf_mol_block")

    raw_atoms = [_parse_sdf_atom(line) for line in lines[4 : 4 + atom_count]]
    charge_overrides: Dict[int, int] = {}
    property_start = 4 + atom_count + bond_count
    for line in lines[property_start:]:
        if not line.startswith("M  CHG"):
            continue
        parts = line.split()
        try:
            pair_count = int(parts[2])
            for offset in range(pair_count):
                atom_index = int(parts[3 + 2 * offset]) - 1
                charge_overrides[atom_index] = int(parts[4 + 2 * offset])
        except (IndexError, ValueError):
            raise ValueError("invalid_sdf_charge_record")

    heavy_raw_indices = [index for index, atom in enumerate(raw_atoms) if atom[3] != "H"]
    if not heavy_raw_indices:
        raise ValueError("no_heavy_atoms")
    raw_to_heavy = {raw: heavy for heavy, raw in enumerate(heavy_raw_indices)}
    graph = nx.Graph()
    coords: List[List[float]] = []
    elements: List[str] = []
    for heavy_index, raw_index in enumerate(heavy_raw_indices):
        x, y, z, element, formal_charge, stereo = raw_atoms[raw_index]
        formal_charge = charge_overrides.get(raw_index, formal_charge)
        coords.append([x, y, z])
        elements.append(element)
        graph.add_node(
            heavy_index,
            element=element,
            formal_charge=int(formal_charge),
            stereo=str(stereo),
        )
    for line in lines[4 + atom_count : 4 + atom_count + bond_count]:
        atom_a, atom_b, order, stereo = _parse_sdf_bond(line)
        if atom_a in raw_to_heavy and atom_b in raw_to_heavy:
            graph.add_edge(
                raw_to_heavy[atom_a],
                raw_to_heavy[atom_b],
                order=order,
                stereo=str(stereo),
            )
    if graph.number_of_nodes() != len(coords):
        raise ValueError("topology_node_count_mismatch")
    return GraphPose(
        coords=np.asarray(coords, dtype=float),
        elements=tuple(elements),
        graph=graph,
        topology_path=source.as_posix(),
        topology_sha256=_sha256(source),
        topology_pose_index=int(pose_index),
    )


def attach_coordinates_to_topology(
    coords: np.ndarray,
    elements: Sequence[str],
    topology_file: Path,
    *,
    topology_pose_index: int = 1,
) -> GraphPose:
    """Attach pose coordinates to an explicitly declared ordered topology."""
    topology = load_sdf_graph_pose(topology_file, topology_pose_index)
    normalized_elements = tuple(_normalize_element(value) for value in elements)
    pose_coords = np.asarray(coords, dtype=float)
    if pose_coords.shape != topology.coords.shape:
        raise ValueError("atom_count_mismatch")
    if normalized_elements != topology.elements:
        raise ValueError("topology_element_order_mismatch")
    topology.coords = pose_coords
    return topology


def _node_match(left: Dict[str, object], right: Dict[str, object]) -> bool:
    return (
        left.get("element") == right.get("element")
        and int(left.get("formal_charge", 0)) == int(right.get("formal_charge", 0))
        and str(left.get("stereo", "0")) == str(right.get("stereo", "0"))
    )


def _edge_match(left: Dict[str, object], right: Dict[str, object]) -> bool:
    return left.get("order") == right.get("order") and str(left.get("stereo", "0")) == str(
        right.get("stereo", "0")
    )


def kabsch_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    left = np.asarray(coords_a, dtype=float)
    right = np.asarray(coords_b, dtype=float)
    if left.shape != right.shape or left.ndim != 2 or left.shape[1] != 3 or left.shape[0] == 0:
        return float("nan")
    left_centered = left - np.mean(left, axis=0)
    right_centered = right - np.mean(right, axis=0)
    covariance = left_centered.T @ right_centered
    u, _, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T
    aligned = left_centered @ rotation
    delta = aligned - right_centered
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def compare_graph_poses(
    pose_a: GraphPose,
    pose_b: GraphPose,
    *,
    max_isomorphisms: int = DEFAULT_MAX_ISOMORPHISMS,
) -> AtomMappingResult:
    count_a = int(pose_a.coords.shape[0])
    count_b = int(pose_b.coords.shape[0])
    base = {
        "total_heavy_atoms_a": count_a,
        "total_heavy_atoms_b": count_b,
        "backend_version": str(getattr(nx, "__version__", "")) if nx is not None else "",
        "topology_a": pose_a.topology_path,
        "topology_b": pose_b.topology_path,
        "topology_sha256_a": pose_a.topology_sha256,
        "topology_sha256_b": pose_b.topology_sha256,
    }
    if nx is None:
        return AtomMappingResult(status="not_comparable", reason="skipped_missing_dependency:networkx", **base)
    if count_a != count_b:
        return AtomMappingResult(status="not_comparable", reason="atom_count_mismatch", **base)
    if count_a <= 0:
        return AtomMappingResult(status="not_comparable", reason="no_heavy_atoms", **base)

    shared_lineage = (
        bool(pose_a.topology_sha256)
        and pose_a.topology_sha256 == pose_b.topology_sha256
        and pose_a.topology_pose_index == pose_b.topology_pose_index
        and pose_a.elements == pose_b.elements
        and list(pose_a.graph.nodes(data=True)) == list(pose_b.graph.nodes(data=True))
        and list(pose_a.graph.edges(data=True)) == list(pose_b.graph.edges(data=True))
    )
    if shared_lineage:
        rmsd = kabsch_rmsd(pose_a.coords, pose_b.coords)
        if np.isfinite(rmsd):
            return AtomMappingResult(
                status="comparable",
                reason="shared_topology_lineage_exact_order",
                rmsd_angstrom=float(rmsd),
                mapped_heavy_atoms=count_a,
                mapping_coverage=1.0,
                valid_mapping_count=1,
                selected_mapping=tuple(range(count_a)),
                **base,
            )

    matcher = nx.algorithms.isomorphism.GraphMatcher(
        pose_a.graph,
        pose_b.graph,
        node_match=_node_match,
        edge_match=_edge_match,
    )
    minimum = float("inf")
    selected: Optional[Tuple[int, ...]] = None
    mapping_count = 0
    for mapping in matcher.isomorphisms_iter():
        mapping_count += 1
        if mapping_count > int(max_isomorphisms):
            return AtomMappingResult(
                status="not_comparable",
                reason=f"mapping_limit_exceeded:{max_isomorphisms}",
                valid_mapping_count=mapping_count - 1,
                **base,
            )
        order = tuple(int(mapping[index]) for index in range(count_a))
        candidate = kabsch_rmsd(pose_a.coords, pose_b.coords[list(order), :])
        if not np.isfinite(candidate):
            continue
        if candidate < minimum - 1e-12 or (
            abs(candidate - minimum) <= 1e-12 and (selected is None or order < selected)
        ):
            minimum = float(candidate)
            selected = order

    if mapping_count == 0:
        return AtomMappingResult(status="not_comparable", reason="graphs_not_isomorphic", **base)
    if selected is None or not np.isfinite(minimum):
        return AtomMappingResult(
            status="not_comparable",
            reason="rmsd_not_finite",
            valid_mapping_count=mapping_count,
            **base,
        )
    return AtomMappingResult(
        status="comparable",
        reason="graph_isomorphism_complete",
        rmsd_angstrom=minimum,
        mapped_heavy_atoms=count_a,
        mapping_coverage=1.0,
        valid_mapping_count=mapping_count,
        selected_mapping=selected,
        **base,
    )


def not_comparable_result(
    reason: str,
    *,
    total_a: int = 0,
    total_b: int = 0,
) -> AtomMappingResult:
    return AtomMappingResult(
        status="not_comparable",
        reason=str(reason),
        rmsd_angstrom=None,
        total_heavy_atoms_a=int(total_a),
        total_heavy_atoms_b=int(total_b),
        backend_version=str(getattr(nx, "__version__", "")) if nx is not None else "",
    )
