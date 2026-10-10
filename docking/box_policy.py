"""Spec 036 R5a docking-box policy.

Default box method ``rg_scaled_v1``: cubic edge = 2.9 x radius of gyration of the
docked ligand's heavy atoms (Feinstein & Brylinski, J Cheminform 2015, 7:18).
The Rg is conformation-dependent and is computed from the input conformer
(the prepared ligand file). A fixed edge is allowed only as an explicit user
value (``user_fixed``).

Scientific choice recorded in the output: ``RG_DEFINITION`` is the unweighted
(geometric) Rg over heavy atoms. It is written into every pair record so it can
be audited or changed in one place.

Spec 037 R2b: the Rg is read first from the ligand preparation provenance
(``preparation_steps/<ligand>.json``, computed once on the prepared conformer and
checked against the file's SHA-256). When that record is missing or does not match
the file, the Rg is computed here with ``rg_source = computed_at_box_step`` and a
note is added to the box warnings. The definition itself is an unverified assumption
(Spec 037 R2c; the methods section of Feinstein & Brylinski 2015 was not retrieved).
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


BOX_METHOD_RG = "rg_scaled_v1"
BOX_METHOD_USER_FIXED = "user_fixed"
RG_EDGE_FACTOR = 2.9
RG_DEFINITION = "heavy_atom_unweighted_v1"
RG_SOURCE_PREPARATION = "ligand_preparation"
RG_SOURCE_COMPUTED = "computed_at_box_step"
RG_DEFINITION_STATUS = "unverified_assumption_spec037_r2c"
PREPARATION_STEPS_DIRNAME = "preparation_steps"
EDGE_WARNING_ANGSTROM = 30.0
CONTAINMENT_NOT_EVALUATED = "not_evaluated_no_reference"
CONTAINMENT_CONTAINED = "contained"
CONTAINMENT_NOT_CONTAINED = "not_contained"

_HYDROGEN_TYPES = {"H", "HD", "HS", "D"}


class BoxError(ValueError):
    """Raised when a box cannot be defined without guessing."""


@dataclass
class PairBox:
    box_method: str
    ligand_rg_angstrom: Optional[float]
    edge_angstrom: float
    rg_definition: str = RG_DEFINITION
    ligand_heavy_atom_count: int = 0
    ligand_file: str = ""
    warnings: List[str] = field(default_factory=list)
    rg_source: str = ""


def _is_heavy_type(token: str) -> bool:
    symbol = str(token or "").strip().upper()
    if not symbol:
        return False
    return symbol.split(".")[0] not in _HYDROGEN_TYPES


def _pdb_like_heavy_coords(path: Path, *, pdbqt: bool) -> np.ndarray:
    coords: List[List[float]] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        models_seen = 0
        for raw_line in handle:
            if raw_line.startswith("MODEL"):
                models_seen += 1
                if models_seen > 1:
                    break
                continue
            if raw_line.startswith("ENDMDL"):
                break
            if not raw_line.startswith(("ATOM", "HETATM")):
                continue
            try:
                x = float(raw_line[30:38])
                y = float(raw_line[38:46])
                z = float(raw_line[46:54])
            except ValueError:
                continue
            if pdbqt:
                if len(raw_line) >= 79 and raw_line[77:79].strip():
                    symbol = raw_line[77:79].strip()
                else:
                    symbol = raw_line.split()[-1]
            else:
                symbol = raw_line[76:78].strip() if len(raw_line) >= 78 and raw_line[76:78].strip() else raw_line[12:16].strip()[:1]
            if _is_heavy_type(symbol):
                coords.append([x, y, z])
    return np.asarray(coords, dtype=float).reshape(-1, 3)


def _sdf_heavy_coords(path: Path) -> np.ndarray:
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.read().splitlines()
    if len(lines) < 4:
        raise BoxError(f"SDF file too short to hold a molecule block: {path}")
    try:
        natoms = int(lines[3][0:3])
    except ValueError as exc:
        raise BoxError(f"SDF counts line unreadable in {path}") from exc
    coords: List[List[float]] = []
    for line in lines[4 : 4 + natoms]:
        tokens = line.split()
        if len(tokens) < 4:
            continue
        if _is_heavy_type(tokens[3]):
            coords.append([float(tokens[0]), float(tokens[1]), float(tokens[2])])
    return np.asarray(coords, dtype=float).reshape(-1, 3)


def _mol2_heavy_coords(path: Path) -> np.ndarray:
    coords: List[List[float]] = []
    in_atoms = False
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if stripped.startswith("@<TRIPOS>"):
                in_atoms = stripped == "@<TRIPOS>ATOM"
                continue
            if not in_atoms or not stripped:
                continue
            tokens = stripped.split()
            if len(tokens) < 6:
                continue
            if _is_heavy_type(tokens[5]):
                coords.append([float(tokens[2]), float(tokens[3]), float(tokens[4])])
    return np.asarray(coords, dtype=float).reshape(-1, 3)


def heavy_atom_coordinates(path: Path) -> np.ndarray:
    """Heavy-atom Cartesian coordinates (N x 3) of the first molecule in a file."""
    path = Path(path)
    if not path.is_file():
        raise BoxError(f"Structure file not found: {path}")
    suffix = path.suffix.lower()
    if suffix == ".pdbqt":
        coords = _pdb_like_heavy_coords(path, pdbqt=True)
    elif suffix == ".pdb":
        coords = _pdb_like_heavy_coords(path, pdbqt=False)
    elif suffix == ".sdf":
        coords = _sdf_heavy_coords(path)
    elif suffix == ".mol2":
        coords = _mol2_heavy_coords(path)
    else:
        raise BoxError(f"Unsupported ligand format for box calculation: {path.name}")
    if coords.size == 0:
        raise BoxError(f"No heavy atoms found in {path.name}")
    return coords


def radius_of_gyration(coords: np.ndarray) -> float:
    """Unweighted radius of gyration of a set of points (RG_DEFINITION)."""
    points = np.asarray(coords, dtype=float).reshape(-1, 3)
    if points.shape[0] == 0:
        raise BoxError("Radius of gyration needs at least one atom")
    centroid = points.mean(axis=0)
    return float(math.sqrt(np.mean(np.sum((points - centroid) ** 2, axis=1))))


def rg_scaled_edge(rg_angstrom: float) -> float:
    return float(RG_EDGE_FACTOR * rg_angstrom)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json_record(path: Path) -> Optional[Dict[str, object]]:
    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


def ligand_preparation_record(ligand_path: Path) -> Tuple[Optional[Dict[str, object]], str]:
    """Ligand preparation provenance of a prepared file: ``(record, reason)``.

    The record lives in ``<prepared ligands>/preparation_steps``. The file ``<stem>.json`` is tried first,
    then any record whose ``output_file`` has the same name. Zero or several matches give no record.
    """
    ligand_path = Path(ligand_path)
    directory = ligand_path.parent / PREPARATION_STEPS_DIRNAME
    if not directory.is_dir():
        return None, "no_preparation_steps_directory"
    direct = directory / f"{ligand_path.stem}.json"
    if direct.is_file():
        payload = _read_json_record(direct)
        if payload and Path(str(payload.get("output_file") or "")).name == ligand_path.name:
            return payload, ""
    matches = []
    for report in sorted(directory.glob("*.json")):
        payload = _read_json_record(report)
        if payload and Path(str(payload.get("output_file") or "")).name == ligand_path.name:
            matches.append(payload)
    if not matches:
        return None, "no_preparation_record_for_output"
    if len(matches) > 1:
        return None, "ambiguous_preparation_records"
    return matches[0], ""


def rg_from_preparation(ligand_path: Path) -> Tuple[Optional[float], Optional[int], str]:
    """Rg recorded by ligand preparation for this exact file: ``(rg, heavy_atom_count, reason_if_none)``."""
    record, reason = ligand_preparation_record(ligand_path)
    if record is None:
        return None, None, f"no ligand preparation provenance ({reason})"
    rg = record.get("radius_of_gyration_angstrom")
    if rg is None or record.get("rg_method") != RG_DEFINITION:
        return None, None, "ligand preparation provenance has no radius_of_gyration_angstrom with rg_method heavy_atom_unweighted_v1"
    recorded_sha = str(record.get("output_sha256") or "")
    if not recorded_sha:
        return None, None, "ligand preparation provenance has no output_sha256 to tie the Rg to the file"
    if recorded_sha != file_sha256(ligand_path):
        return None, None, "the file changed after ligand preparation (output_sha256 differs)"
    count = record.get("rg_heavy_atom_count")
    return float(rg), (int(count) if count is not None else None), ""


def compute_ligand_box(ligand_path: Path, fixed_edge: Optional[float] = None) -> PairBox:
    """Box for one ligand.

    The Rg comes from the ligand preparation provenance when it is recorded for this file
    (``rg_source = ligand_preparation``). Otherwise it is computed here
    (``rg_source = computed_at_box_step``) and a note is added to ``warnings``.

    ``fixed_edge`` given -> ``user_fixed`` (explicit user value; Rg still recorded
    when the conformer can be read). ``fixed_edge`` None -> ``rg_scaled_v1``,
    which requires a readable conformer. There is no silent fallback size.
    """
    ligand_path = Path(ligand_path)
    rg, count, reason = rg_from_preparation(ligand_path)
    notes: List[str] = []
    if rg is not None:
        source = RG_SOURCE_PREPARATION
    else:
        source = RG_SOURCE_COMPUTED
        notes.append(
            f"[box] {ligand_path.name}: Rg computed at the box step ({reason}); "
            f"{RG_DEFINITION} from the file, not from ligand preparation."
        )
    if fixed_edge is not None:
        edge = float(fixed_edge)
        if not math.isfinite(edge) or edge <= 0:
            raise BoxError(f"Explicit box size must be a positive number of Angstrom, got {fixed_edge!r}")
        if rg is None:
            try:
                coords = heavy_atom_coordinates(ligand_path)
                rg = radius_of_gyration(coords)
                count = int(coords.shape[0])
            except BoxError:
                rg, count = None, 0
        elif count is None:
            count = int(heavy_atom_coordinates(ligand_path).shape[0])
        return PairBox(
            box_method=BOX_METHOD_USER_FIXED,
            ligand_rg_angstrom=rg,
            edge_angstrom=edge,
            ligand_heavy_atom_count=int(count or 0),
            ligand_file=ligand_path.name,
            warnings=notes,
            rg_source=source,
        )
    if rg is None:
        coords = heavy_atom_coordinates(ligand_path)
        rg = radius_of_gyration(coords)
        count = int(coords.shape[0])
    elif count is None:
        count = int(heavy_atom_coordinates(ligand_path).shape[0])
    return PairBox(
        box_method=BOX_METHOD_RG,
        ligand_rg_angstrom=rg,
        edge_angstrom=rg_scaled_edge(rg),
        ligand_heavy_atom_count=int(count),
        ligand_file=ligand_path.name,
        warnings=notes,
        rg_source=source,
    )


def containment_status(
    center: Sequence[float],
    edge: float,
    reference_coords: Optional[np.ndarray],
) -> Dict[str, object]:
    """Check centre +/- edge/2 against reference heavy atoms.

    Returns status plus the indices of reference atoms outside the box.
    """
    if reference_coords is None:
        return {"status": CONTAINMENT_NOT_EVALUATED, "outside_count": 0, "reference_atom_count": 0}
    points = np.asarray(reference_coords, dtype=float).reshape(-1, 3)
    half = float(edge) / 2.0
    centre = np.asarray(center, dtype=float).reshape(3)
    inside = np.all(np.abs(points - centre) <= half + 1e-9, axis=1)
    outside = int((~inside).sum())
    return {
        "status": CONTAINMENT_CONTAINED if outside == 0 else CONTAINMENT_NOT_CONTAINED,
        "outside_count": outside,
        "reference_atom_count": int(points.shape[0]),
    }


def box_warnings(
    *,
    pair_label: str,
    edge_angstrom: float,
    containment: Dict[str, object],
    reference_label: str = "",
) -> List[str]:
    messages: List[str] = []
    if containment.get("status") == CONTAINMENT_NOT_CONTAINED:
        messages.append(
            f"[box] {pair_label}: box does not contain all reference heavy atoms"
            f" ({containment.get('outside_count')}/{containment.get('reference_atom_count')} outside)"
            f" of reference ligand {reference_label or '?'}."
        )
    if float(edge_angstrom) > EDGE_WARNING_ANGSTROM:
        messages.append(
            f"[box] {pair_label}: edge {float(edge_angstrom):.2f} A exceeds {EDGE_WARNING_ANGSTROM:.0f} A"
            " (AutoDock Vina FAQ: large search spaces reduce sampling efficiency)."
        )
    return messages


def load_reference_coordinates(path: Optional[Path]) -> Optional[np.ndarray]:
    if path is None:
        return None
    return heavy_atom_coordinates(Path(path))


def describe_method(method: str) -> str:
    return {
        BOX_METHOD_RG: f"cubic edge = {RG_EDGE_FACTOR} x Rg of ligand heavy atoms ({RG_DEFINITION})",
        BOX_METHOD_USER_FIXED: "explicit user-fixed edge",
    }.get(method, "unrecorded")


def join_messages(messages: Iterable[str]) -> str:
    return " | ".join(str(message) for message in messages if str(message))
