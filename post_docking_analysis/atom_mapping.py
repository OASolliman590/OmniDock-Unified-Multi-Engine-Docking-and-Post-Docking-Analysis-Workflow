"""Scientifically conservative graph/symmetry-aware heavy-atom RMSD.

Spec 031 method ``spec031-atom-mapping-v1`` requires complete heavy-atom
graph isomorphism.  This module never falls back to element sorting,
nearest-neighbour assignment, or common-substructure truncation.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, replace
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

# RMSD frames.  "kabsch" is the accepted Spec 031 value (superposed).  "in_place_v1" is the
# Spec 034 T002a redocking primary (no superposition, same valid mappings).
KABSCH_METHOD = "kabsch"
IN_PLACE_METHOD = "in_place_v1"

# Spec 034 R1: lineage source for Meeko-prepared PDBQT poses. Explicit SDF
# topologies (the Spec 031 source) are labelled EXPLICIT_TOPOLOGY_SOURCE.
LINEAGE_METHOD = "meeko_smiles_idx_lineage_v1"
EXPLICIT_TOPOLOGY_SOURCE = "explicit_sdf_topology"
NO_LINEAGE_REASON = "no_lineage_remarks"

# AutoDock atom type -> element symbol (shared by PDBQT readers).
AUTODOCK_ELEMENTS = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N",
    "OA": "O", "OS": "O", "O": "O", "S": "S", "SA": "S",
    "P": "P", "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
    "H": "H", "HD": "H", "HS": "H", "MG": "Mg", "MN": "Mn",
    "ZN": "Zn", "CA": "Ca", "FE": "Fe", "CU": "Cu",
}


class LineageError(ValueError):
    """A PDBQT lineage pose cannot be built; ``reason`` is machine-readable."""

    def __init__(self, reason: str) -> None:
        super().__init__(str(reason))
        self.reason = str(reason)


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
    lineage_source: str = EXPLICIT_TOPOLOGY_SOURCE


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
    lineage_source: str = ""
    reference_lineage_source: str = ""
    rmsd_frame: str = ""
    kabsch_rmsd_angstrom: Optional[float] = None

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
            "lineage_source": self.lineage_source,
            "reference_lineage_source": self.reference_lineage_source,
            "rmsd_frame": self.rmsd_frame,
            "kabsch_rmsd_angstrom": self.kabsch_rmsd_angstrom,
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


def in_place_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    """Heavy-atom RMSD with no superposition (both poses already share the receptor frame)."""
    left = np.asarray(coords_a, dtype=float)
    right = np.asarray(coords_b, dtype=float)
    if left.shape != right.shape or left.ndim != 2 or left.shape[1] != 3 or left.shape[0] == 0:
        return float("nan")
    delta = left - right
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def _pair_base_fields(pose_a: GraphPose, pose_b: GraphPose, *, rmsd_frame: str, alignment_method: str) -> Dict[str, object]:
    return {
        "total_heavy_atoms_a": int(pose_a.coords.shape[0]),
        "total_heavy_atoms_b": int(pose_b.coords.shape[0]),
        "backend_version": str(getattr(nx, "__version__", "")) if nx is not None else "",
        "topology_a": pose_a.topology_path,
        "topology_b": pose_b.topology_path,
        "topology_sha256_a": pose_a.topology_sha256,
        "topology_sha256_b": pose_b.topology_sha256,
        "lineage_source": pose_a.lineage_source,
        "reference_lineage_source": pose_b.lineage_source,
        "rmsd_frame": rmsd_frame,
        "alignment_method": alignment_method,
    }


def _shared_lineage(pose_a: GraphPose, pose_b: GraphPose) -> bool:
    return (
        bool(pose_a.topology_sha256)
        and pose_a.topology_sha256 == pose_b.topology_sha256
        and pose_a.topology_pose_index == pose_b.topology_pose_index
        and pose_a.elements == pose_b.elements
        and list(pose_a.graph.nodes(data=True)) == list(pose_b.graph.nodes(data=True))
        and list(pose_a.graph.edges(data=True)) == list(pose_b.graph.edges(data=True))
    )


@dataclass
class _MappingSearch:
    """Minima over the same valid symmetry mappings, for both RMSD frames."""

    mapping_count: int = 0
    limit_exceeded: bool = False
    kabsch_min: float = float("inf")
    kabsch_order: Optional[Tuple[int, ...]] = None
    in_place_min: float = float("inf")
    in_place_order: Optional[Tuple[int, ...]] = None


def _offer_candidate(
    search: _MappingSearch,
    frame: str,
    candidate: float,
    order: Tuple[int, ...],
) -> None:
    if not np.isfinite(candidate):
        return
    if frame == KABSCH_METHOD:
        current, selected = search.kabsch_min, search.kabsch_order
    else:
        current, selected = search.in_place_min, search.in_place_order
    if candidate < current - 1e-12 or (abs(candidate - current) <= 1e-12 and (selected is None or order < selected)):
        if frame == KABSCH_METHOD:
            search.kabsch_min, search.kabsch_order = float(candidate), order
        else:
            search.in_place_min, search.in_place_order = float(candidate), order


def _search_mappings(pose_a: GraphPose, pose_b: GraphPose, max_isomorphisms: int) -> _MappingSearch:
    """Enumerate every bond-labelled, element- and charge-matched isomorphism once."""
    count_a = int(pose_a.coords.shape[0])
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        pose_a.graph,
        pose_b.graph,
        node_match=_node_match,
        edge_match=_edge_match,
    )
    search = _MappingSearch()
    for mapping in matcher.isomorphisms_iter():
        search.mapping_count += 1
        if search.mapping_count > int(max_isomorphisms):
            search.mapping_count -= 1
            search.limit_exceeded = True
            break
        order = tuple(int(mapping[index]) for index in range(count_a))
        permuted = pose_b.coords[list(order), :]
        _offer_candidate(search, KABSCH_METHOD, kabsch_rmsd(pose_a.coords, permuted), order)
        _offer_candidate(search, IN_PLACE_METHOD, in_place_rmsd(pose_a.coords, permuted), order)
    return search


def _precheck(pose_a: GraphPose, pose_b: GraphPose, base: Dict[str, object]) -> Optional[AtomMappingResult]:
    count_a = int(pose_a.coords.shape[0])
    count_b = int(pose_b.coords.shape[0])
    if nx is None:
        return AtomMappingResult(status="not_comparable", reason="skipped_missing_dependency:networkx", **base)
    if count_a != count_b:
        return AtomMappingResult(status="not_comparable", reason="atom_count_mismatch", **base)
    if count_a <= 0:
        return AtomMappingResult(status="not_comparable", reason="no_heavy_atoms", **base)
    return None


def compare_graph_poses(
    pose_a: GraphPose,
    pose_b: GraphPose,
    *,
    max_isomorphisms: int = DEFAULT_MAX_ISOMORPHISMS,
) -> AtomMappingResult:
    """Spec 031 comparison: minimum Kabsch-superposed heavy-atom RMSD over valid mappings."""
    count_a = int(pose_a.coords.shape[0])
    base = _pair_base_fields(pose_a, pose_b, rmsd_frame=KABSCH_METHOD, alignment_method=ALIGNMENT_METHOD)
    early = _precheck(pose_a, pose_b, base)
    if early is not None:
        return early

    if _shared_lineage(pose_a, pose_b):
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
                kabsch_rmsd_angstrom=float(rmsd),
                **base,
            )

    search = _search_mappings(pose_a, pose_b, max_isomorphisms)
    if search.limit_exceeded:
        return AtomMappingResult(
            status="not_comparable",
            reason=f"mapping_limit_exceeded:{max_isomorphisms}",
            valid_mapping_count=search.mapping_count,
            **base,
        )
    if search.mapping_count == 0:
        return AtomMappingResult(status="not_comparable", reason="graphs_not_isomorphic", **base)
    if search.kabsch_order is None or not np.isfinite(search.kabsch_min):
        return AtomMappingResult(
            status="not_comparable",
            reason="rmsd_not_finite",
            valid_mapping_count=search.mapping_count,
            **base,
        )
    return AtomMappingResult(
        status="comparable",
        reason="graph_isomorphism_complete",
        rmsd_angstrom=search.kabsch_min,
        mapped_heavy_atoms=count_a,
        mapping_coverage=1.0,
        valid_mapping_count=search.mapping_count,
        selected_mapping=search.kabsch_order,
        kabsch_rmsd_angstrom=search.kabsch_min,
        **base,
    )


def compare_graph_poses_in_place(
    pose_a: GraphPose,
    pose_b: GraphPose,
    *,
    max_isomorphisms: int = DEFAULT_MAX_ISOMORPHISMS,
) -> AtomMappingResult:
    """Spec 034 T002a redocking primary: minimum in-place heavy-atom RMSD over the same
    valid mappings used by ``compare_graph_poses``.  The Kabsch minimum over that mapping
    set is reported as ``kabsch_rmsd_angstrom`` (secondary, not used for classification).
    """
    count_a = int(pose_a.coords.shape[0])
    base = _pair_base_fields(pose_a, pose_b, rmsd_frame=IN_PLACE_METHOD, alignment_method="in_place")
    early = _precheck(pose_a, pose_b, base)
    if early is not None:
        return early

    if _shared_lineage(pose_a, pose_b):
        in_place = in_place_rmsd(pose_a.coords, pose_b.coords)
        if np.isfinite(in_place):
            kabsch = kabsch_rmsd(pose_a.coords, pose_b.coords)
            return AtomMappingResult(
                status="comparable",
                reason="shared_topology_lineage_exact_order",
                rmsd_angstrom=float(in_place),
                mapped_heavy_atoms=count_a,
                mapping_coverage=1.0,
                valid_mapping_count=1,
                selected_mapping=tuple(range(count_a)),
                kabsch_rmsd_angstrom=float(kabsch) if np.isfinite(kabsch) else None,
                **base,
            )

    search = _search_mappings(pose_a, pose_b, max_isomorphisms)
    if search.limit_exceeded:
        return AtomMappingResult(
            status="not_comparable",
            reason=f"mapping_limit_exceeded:{max_isomorphisms}",
            valid_mapping_count=search.mapping_count,
            **base,
        )
    if search.mapping_count == 0:
        return AtomMappingResult(status="not_comparable", reason="graphs_not_isomorphic", **base)
    if search.in_place_order is None or not np.isfinite(search.in_place_min):
        return AtomMappingResult(
            status="not_comparable",
            reason="rmsd_not_finite",
            valid_mapping_count=search.mapping_count,
            **base,
        )
    kabsch_secondary = search.kabsch_min if np.isfinite(search.kabsch_min) else None
    return AtomMappingResult(
        status="comparable",
        reason="graph_isomorphism_complete",
        rmsd_angstrom=search.in_place_min,
        mapped_heavy_atoms=count_a,
        mapping_coverage=1.0,
        valid_mapping_count=search.mapping_count,
        selected_mapping=search.in_place_order,
        kabsch_rmsd_angstrom=kabsch_secondary,
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


# ---------------------------------------------------------------------------
# Spec 034 R1: Meeko SMILES / SMILES IDX lineage source (meeko_smiles_idx_lineage_v1)
#
# Meeko writes ``REMARK SMILES <smiles>`` and ``REMARK SMILES IDX`` pairs
# ``<smiles_index_1based> <pdbqt_serial_1based>`` (meeko.writer.remark_index_map).
# Vina copies these remarks into every output MODEL.  The graph built here has one
# node per PDBQT heavy atom, in PDBQT file order, so node i carries the coordinates of
# the i-th heavy ATOM record of the selected model.
# ---------------------------------------------------------------------------


def _load_rdkit_chem():
    try:
        from rdkit import Chem  # type: ignore
    except ImportError as exc:
        raise RuntimeError("skipped_missing_dependency:rdkit") from exc
    return Chem


def _pdbqt_model_lines(path: Path, model_index: int) -> List[str]:
    try:
        text_lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    except (OSError, UnicodeDecodeError) as exc:
        raise LineageError("unreadable_pdbqt") from exc
    if not any(line.startswith("MODEL") for line in text_lines):
        if int(model_index) != 1:
            raise LineageError(f"model_not_found:{model_index}")
        return text_lines
    selected: List[str] = []
    capture = False
    found = False
    for line in text_lines:
        if line.startswith("MODEL"):
            try:
                number = int(line.split()[1])
            except (IndexError, ValueError) as exc:
                raise LineageError("invalid_pdbqt_model_record") from exc
            capture = number == int(model_index)
            found = found or capture
            continue
        if line.startswith("ENDMDL"):
            if capture:
                break
            capture = False
            continue
        if capture:
            selected.append(line)
    if not found:
        raise LineageError(f"model_not_found:{model_index}")
    return selected


def _pdbqt_atoms(lines: Sequence[str]) -> List[Tuple[int, str, Tuple[float, float, float], bool]]:
    """Return (serial, element, xyz, is_heavy) for every ATOM/HETATM record."""
    atoms: List[Tuple[int, str, Tuple[float, float, float], bool]] = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        parts = line.split()
        try:
            serial = int(line[6:11])
        except ValueError:
            try:
                serial = int(parts[1])
            except (IndexError, ValueError) as exc:
                raise LineageError("invalid_pdbqt_atom_serial") from exc
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            try:
                xyz = (float(parts[6]), float(parts[7]), float(parts[8]))
            except (IndexError, ValueError) as exc:
                raise LineageError("invalid_pdbqt_coordinates") from exc
        token = parts[-1].upper() if parts else ""
        element = AUTODOCK_ELEMENTS.get(token)
        if not element:
            raw_element = line[76:78].strip().upper() if len(line) >= 78 else ""
            element = AUTODOCK_ELEMENTS.get(raw_element)
        if not element:
            raise LineageError(f"unknown_autodock_atom_type:{token}")
        atoms.append((serial, element, xyz, element.upper() != "H"))
    return atoms


def _lineage_remarks(lines: Sequence[str]) -> Tuple[str, List[int]]:
    smiles: List[str] = []
    flat_pairs: List[int] = []
    idx_seen = False
    for line in lines:
        if line.startswith("REMARK SMILES IDX"):
            try:
                values = [int(token) for token in line.split()[3:]]
            except ValueError as exc:
                raise LineageError("invalid_lineage_idx") from exc
            if len(values) % 2:
                raise LineageError("invalid_lineage_idx")
            flat_pairs.extend(values)
            idx_seen = True
        elif line.startswith("REMARK SMILES "):
            parts = line.split()
            if len(parts) < 3:
                raise LineageError("invalid_lineage_smiles")
            smiles.append(parts[2])
    if not smiles and not idx_seen:
        raise LineageError(NO_LINEAGE_REASON)
    if len(smiles) > 1:
        raise LineageError("ambiguous_lineage_smiles")
    if not smiles:
        raise LineageError("missing_lineage_smiles")
    if not idx_seen:
        raise LineageError("incomplete_lineage_idx")
    return smiles[0], flat_pairs


def read_pdbqt_lineage_smiles(path: Path, model_index: int = 1) -> Optional[str]:
    """Return the Meeko ``REMARK SMILES`` text of one PDBQT model (Spec 036 R1c).

    Returns None when the model carries no SMILES remark. Raises ``LineageError`` when the
    file or model cannot be read. Only the SMILES text is read here; the strict lineage
    checks stay in ``load_pdbqt_lineage_pose``.
    """
    source = Path(path).expanduser().resolve()
    if source.suffix.lower() != ".pdbqt":
        raise LineageError("unsupported_pose_format")
    if int(model_index) < 1:
        raise LineageError("invalid_pose_index")
    for line in _pdbqt_model_lines(source, int(model_index)):
        if line.startswith("REMARK SMILES ") and not line.startswith("REMARK SMILES IDX"):
            text = line[len("REMARK SMILES "):].strip()
            return text or None
    return None


def load_pdbqt_lineage_pose(path: Path, model_index: int = 1) -> GraphPose:
    """Build a heavy-atom, bond-labelled pose graph from Meeko lineage remarks.

    Raises ``LineageError`` whose ``reason`` is machine-readable.  ``no_lineage_remarks``
    means the file carries no SMILES/IDX remarks; every other reason is a lineage that
    is present but cannot be used, and must not fall back silently.
    """
    if nx is None:
        raise LineageError("skipped_missing_dependency:networkx")
    source = Path(path).expanduser().resolve()
    if source.suffix.lower() != ".pdbqt":
        raise LineageError("unsupported_pose_format")
    if int(model_index) < 1:
        raise LineageError("invalid_pose_index")
    lines = _pdbqt_model_lines(source, int(model_index))
    smiles_text, flat_pairs = _lineage_remarks(lines)
    try:
        Chem = _load_rdkit_chem()
    except RuntimeError as exc:
        raise LineageError(str(exc)) from exc
    mol = Chem.MolFromSmiles(smiles_text)
    if mol is None:
        raise LineageError("invalid_lineage_smiles")
    for atom in mol.GetAtoms():
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            raise LineageError("lineage_stereo_unsupported")
    for bond in mol.GetBonds():
        if bond.GetStereo() != Chem.BondStereo.STEREONONE:
            raise LineageError("lineage_stereo_unsupported")

    atoms = _pdbqt_atoms(lines)
    heavy = [(serial, element, xyz) for serial, element, xyz, is_heavy in atoms if is_heavy]
    hydrogen_serials = {serial for serial, _, _, is_heavy in atoms if not is_heavy}
    if not heavy:
        raise LineageError("no_heavy_atoms")
    serial_to_position: Dict[int, int] = {}
    for position, (serial, _, _) in enumerate(heavy):
        if serial in serial_to_position:
            raise LineageError("ambiguous_pdbqt_serial")
        serial_to_position[serial] = position
    pose_element = {serial: element for serial, element, _ in heavy}

    smiles_heavy = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1}
    smiles_atom_count = mol.GetNumAtoms()
    smiles_to_serial: Dict[int, int] = {}
    serials_seen: set = set()
    pairs = [(flat_pairs[index], flat_pairs[index + 1]) for index in range(0, len(flat_pairs), 2)]
    for smiles_index_1based, serial in pairs:
        if not 1 <= smiles_index_1based <= smiles_atom_count:
            raise LineageError("invalid_lineage_idx")
        smiles_index = smiles_index_1based - 1
        if smiles_index not in smiles_heavy or serial in hydrogen_serials:
            raise LineageError("lineage_idx_references_hydrogen")
        if serial not in serial_to_position:
            raise LineageError("invalid_lineage_idx")
        if smiles_index in smiles_to_serial or serial in serials_seen:
            raise LineageError("ambiguous_lineage_idx")
        smiles_to_serial[smiles_index] = serial
        serials_seen.add(serial)
    if set(smiles_to_serial) != smiles_heavy or serials_seen != set(serial_to_position):
        raise LineageError("incomplete_lineage_idx")
    for smiles_index, serial in smiles_to_serial.items():
        if mol.GetAtomWithIdx(smiles_index).GetSymbol() != pose_element[serial]:
            raise LineageError("lineage_idx_element_mismatch")

    graph = nx.Graph()
    smiles_to_node: Dict[int, int] = {}
    for smiles_index, serial in smiles_to_serial.items():
        node = serial_to_position[serial]
        smiles_to_node[smiles_index] = node
        atom = mol.GetAtomWithIdx(smiles_index)
        graph.add_node(
            node,
            element=pose_element[serial],
            formal_charge=int(atom.GetFormalCharge()),
            stereo="0",
            total_h=int(atom.GetTotalNumHs(includeNeighbors=True)),
        )
    for bond in mol.GetBonds():
        begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if begin not in smiles_to_node or end not in smiles_to_node:
            continue
        bond_type = bond.GetBondType()
        if bond_type == Chem.BondType.AROMATIC:
            order = "aromatic"
        elif bond_type == Chem.BondType.SINGLE:
            order = "1"
        elif bond_type == Chem.BondType.DOUBLE:
            order = "2"
        elif bond_type == Chem.BondType.TRIPLE:
            order = "3"
        else:
            raise LineageError(f"unsupported_lineage_bond_type:{bond_type}")
        graph.add_edge(smiles_to_node[begin], smiles_to_node[end], order=order, stereo="0")

    coords = np.asarray([xyz for _, _, xyz in heavy], dtype=float)
    return GraphPose(
        coords=coords,
        elements=tuple(element for _, element, _ in heavy),
        graph=graph,
        topology_path=source.as_posix(),
        topology_sha256=_sha256(source),
        topology_pose_index=int(model_index),
        lineage_source=LINEAGE_METHOD,
    )


def _aromatic_convention_graph(graph) -> "nx.Graph":
    """Return ``graph`` with bond labels in one convention: aromatic rings as ``aromatic``.

    Aromatic labels are Kekulise-then-perceive: Kekule forms of one ring are not unique,
    so a Kekule assignment copied from another source is not a stable label.  RDKit
    aromaticity perception is Kekule-independent and is applied to both sides of a
    lineage comparison.  Element, formal charge and stereo labels are unchanged.
    """
    Chem = _load_rdkit_chem()
    rw = Chem.RWMol()
    nodes = sorted(graph.nodes())
    atom_index: Dict[object, int] = {}
    for node in nodes:
        data = graph.nodes[node]
        atom = Chem.Atom(str(data["element"]))
        atom.SetFormalCharge(int(data.get("formal_charge", 0)))
        if "total_h" in data:
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(int(data["total_h"]))
        atom_index[node] = rw.AddAtom(atom)
    aromatic_nodes: set = set()
    bond_types = {"1": Chem.BondType.SINGLE, "2": Chem.BondType.DOUBLE, "3": Chem.BondType.TRIPLE}
    for begin, end, data in graph.edges(data=True):
        order = str(data.get("order"))
        if order == "aromatic":
            bond_type = Chem.BondType.AROMATIC
            aromatic_nodes.update((begin, end))
        elif order in bond_types:
            bond_type = bond_types[order]
        else:
            raise ValueError(f"unsupported_bond_order:{order}")
        rw.AddBond(atom_index[begin], atom_index[end], bond_type)
        if order == "aromatic":
            rw.GetBondBetweenAtoms(atom_index[begin], atom_index[end]).SetIsAromatic(True)
    for node in aromatic_nodes:
        rw.GetAtomWithIdx(atom_index[node]).SetIsAromatic(True)
    mol = rw.GetMol()
    if aromatic_nodes:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    Chem.SanitizeMol(mol)

    inverse = {index: node for node, index in atom_index.items()}
    normalized = nx.Graph()
    for node in nodes:
        normalized.add_node(node, **dict(graph.nodes[node]))
    for bond in mol.GetBonds():
        begin = inverse[bond.GetBeginAtomIdx()]
        end = inverse[bond.GetEndAtomIdx()]
        if bond.GetIsAromatic():
            order = "aromatic"
        elif bond.GetBondType() == Chem.BondType.SINGLE:
            order = "1"
        elif bond.GetBondType() == Chem.BondType.DOUBLE:
            order = "2"
        elif bond.GetBondType() == Chem.BondType.TRIPLE:
            order = "3"
        else:
            raise ValueError(f"unsupported_normalized_bond:{bond.GetBondType()}")
        stereo = graph.edges[begin, end].get("stereo", "0")
        normalized.add_edge(begin, end, order=order, stereo=str(stereo))
    return normalized


def _normalize_pose_convention(pose: GraphPose) -> GraphPose:
    return replace(pose, graph=_aromatic_convention_graph(pose.graph))


def compare_lineage_graph_poses(
    pose_a: GraphPose,
    pose_b: GraphPose,
    *,
    max_isomorphisms: int = DEFAULT_MAX_ISOMORPHISMS,
    rmsd_frame: str = KABSCH_METHOD,
) -> AtomMappingResult:
    """Compare a lineage pose with a reference graph under one bond-order convention.

    Both graphs go through ``_aromatic_convention_graph`` before the mapping comparison.
    ``rmsd_frame`` selects ``compare_graph_poses`` (``kabsch``, the default, used by the
    Spec 032 atom_map derivation) or ``compare_graph_poses_in_place`` (``in_place_v1``,
    the Spec 034 T002a redocking primary).  Failure to normalise is ``not_comparable``;
    nothing is loosened.  At least one side must carry the lineage source.
    """
    if rmsd_frame not in (KABSCH_METHOD, IN_PLACE_METHOD):
        raise ValueError(f"unknown_rmsd_frame:{rmsd_frame}")
    if LINEAGE_METHOD not in (pose_a.lineage_source, pose_b.lineage_source):
        raise ValueError("lineage_pose_required")
    compare = compare_graph_poses_in_place if rmsd_frame == IN_PLACE_METHOD else compare_graph_poses
    alignment = "in_place" if rmsd_frame == IN_PLACE_METHOD else ALIGNMENT_METHOD
    try:
        normalized_a = _normalize_pose_convention(pose_a)
        normalized_b = _normalize_pose_convention(pose_b)
    except RuntimeError as exc:
        reason = str(exc) if str(exc).startswith("skipped_missing_dependency") else f"bond_convention_failed:{exc}"
        return AtomMappingResult(
            status="not_comparable",
            reason=reason,
            **_pair_base_fields(pose_a, pose_b, rmsd_frame=rmsd_frame, alignment_method=alignment),
        )
    except Exception as exc:  # RDKit sanitisation errors do not share a base class
        return AtomMappingResult(
            status="not_comparable",
            reason=f"bond_convention_failed:{exc.__class__.__name__}",
            **_pair_base_fields(pose_a, pose_b, rmsd_frame=rmsd_frame, alignment_method=alignment),
        )
    return compare(normalized_a, normalized_b, max_isomorphisms=max_isomorphisms)
