"""Strict, provenance-gated export of selected docking poses for MD consumers."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import dataclasses

import pandas as pd

from docking.project_layout import (
    ProtonationPolicyConflict,
    ensure_numbered_output_layout,
    load_manifest,
    load_protonation_policy,
    manifest_path,
)
from post_docking_analysis.atom_mapping import (
    AUTODOCK_ELEMENTS as _AUTODOCK_ELEMENTS,
    LINEAGE_METHOD,
    GraphPose,
    LineageError,
    compare_lineage_graph_poses,
    load_pdbqt_lineage_pose,
    load_sdf_graph_pose,
)
from post_docking_analysis.complex_validation import validate_md_receptor_pdb, validate_md_system_pdb
from post_docking_analysis.docking_parser import parse_vina_pdbqt
from post_docking_analysis.md_chemistry import (
    ChemistryBackend,
    PreparedLigand,
    docked_microspecies_from_pdbqt,
    docked_microspecies_from_sdf,
    get_openbabel_backend,
    openbabel_capability,
)


SCHEMA_VERSION = "md_inputs_manifest_v1"
PROVENANCE_SCHEMA_VERSION = "md_inputs_provenance_v1"
METHOD_VERSION = "md_inputs_charmm_gui_cgenff_v1_1.1"
CONSUMER_PROFILE = "charmm_gui_cgenff_v1"
# Spec 036 R1b: the default MD ligand is the docked microspecies. The legacy flag value
# "openbabel_predicted" and --reprotonate-at-ph are the explicit, recorded override path.
DOCKED_STATE_POLICY = "docked_state"
LEGACY_REPROTONATION_POLICY = "openbabel_predicted"
AUTO_MAPS_SCHEMA_VERSION = "md_inputs_auto_maps_v1"
AUTO_MAPS_DIRNAME = "md_inputs_maps"
_CHARGE_BLIND_MAPPING = "charge_blind_protonation_state_v1"
SUCCESS_STATUSES = {"completed"}
ROW_STATUSES = {
    "completed",
    "failed",
    "skipped_disabled",
    "skipped_missing_dependency",
    "skipped_missing_configuration",
    "not_comparable",
}
_SAFE_KEY = re.compile(r"[^A-Za-z0-9._-]+")
@dataclass(frozen=True)
class MDInputsRequest:
    project_dir: Path
    engine: str
    tags_file: Path
    receptor_map: Path
    pH: Optional[float] = None  # the explicit re-protonation pH; None for the docked-state default
    protonation_policy: str = DOCKED_STATE_POLICY
    topology_map: Optional[Path] = None
    charge_map: Optional[Path] = None
    ligand_formats: Tuple[str, ...] = ("mol2",)
    consumer_profile: str = CONSUMER_PROFILE
    protonate: bool = True
    force: bool = False
    override_source: str = ""  # "", "reprotonate_at_ph" or "legacy_openbabel_predicted"

    @property
    def reprotonation_override(self) -> bool:
        return self.pH is not None

    @classmethod
    def from_dict(cls, payload: Dict[str, object]) -> "MDInputsRequest":
        root = Path(str(payload.get("project_dir") or "")).expanduser().resolve()
        engine = str(payload.get("engine") or "").strip().lower()
        if engine not in {"gnina", "vina", "smina"}:
            raise ValueError("--engine must be one of: gnina, vina, smina")
        consumer = str(payload.get("consumer_profile") or CONSUMER_PROFILE).strip().lower()
        if consumer != CONSUMER_PROFILE:
            raise ValueError(f"unsupported_consumer_profile:{consumer}")
        raw_policy = str(payload.get("protonation_policy") or "").strip().lower()
        if raw_policy not in ("", DOCKED_STATE_POLICY, LEGACY_REPROTONATION_POLICY):
            raise ValueError(f"unsupported_protonation_policy:{raw_policy}")
        ph, override_source = _reprotonation_override(payload, raw_policy)
        protonation_policy = LEGACY_REPROTONATION_POLICY if ph is not None else DOCKED_STATE_POLICY
        formats = tuple(
            sorted(
                {str(item).strip().lower() for item in (payload.get("ligand_formats") or ["mol2"]) if str(item).strip()}
                | {"mol2"}
            )
        )
        if any(item not in {"mol2", "sdf", "pdb"} for item in formats):
            raise ValueError("ligand formats must be mol2, sdf, or pdb")

        def required_path(key: str) -> Path:
            raw = str(payload.get(key) or "").strip()
            if not raw:
                raise ValueError(f"--{key.replace('_', '-')} is required")
            return Path(raw).expanduser().resolve()

        def optional_path(key: str) -> Optional[Path]:
            raw = str(payload.get(key) or "").strip()
            return Path(raw).expanduser().resolve() if raw else None

        request = cls(
            project_dir=root,
            engine=engine,
            tags_file=required_path("tags_file"),
            receptor_map=required_path("receptor_map"),
            topology_map=optional_path("topology_map"),
            charge_map=optional_path("charge_map"),
            pH=ph,
            protonation_policy=protonation_policy,
            ligand_formats=formats,
            consumer_profile=consumer,
            protonate=bool(payload.get("protonate", True)),
            force=bool(payload.get("force", False)),
            override_source=override_source,
        )
        for label, path in (("project", request.project_dir), ("tags", request.tags_file), ("receptor_map", request.receptor_map)):
            if not path.exists():
                raise FileNotFoundError(f"{label}_path_missing:{path}")
        for label, path in (("topology_map", request.topology_map), ("charge_map", request.charge_map)):
            if path is not None and not path.is_file():
                raise FileNotFoundError(f"{label}_missing:{path}")
        return request

    def to_dict(self) -> Dict[str, object]:
        return {
            "project_dir": str(self.project_dir),
            "engine": self.engine,
            "tags_file": str(self.tags_file),
            "receptor_map": str(self.receptor_map),
            "topology_map": str(self.topology_map) if self.topology_map else "",
            "charge_map": str(self.charge_map) if self.charge_map else "",
            "pH": self.pH,
            "protonation_policy": self.protonation_policy,
            "override_source": self.override_source,
            "ligand_formats": list(self.ligand_formats),
            "consumer_profile": self.consumer_profile,
            "protonate": self.protonate,
            "force": self.force,
        }


def _validate_ph(value: object, label: str) -> float:
    try:
        ph = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid_{label}") from exc
    if not math.isfinite(ph) or ph < 0.0 or ph > 14.0:
        raise ValueError(f"{label} must be a finite value from 0 through 14")
    return ph


def _reprotonation_override(payload: Dict[str, object], raw_policy: str) -> Tuple[Optional[float], str]:
    """Spec 036 R1b: explicit re-protonation request. None means the docked-state default."""
    legacy = raw_policy == LEGACY_REPROTONATION_POLICY
    explicit_raw = payload.get("reprotonate_at_ph")
    ph_raw = payload.get("pH")
    has_explicit = explicit_raw not in (None, "")
    has_ph = ph_raw not in (None, "")
    if not has_explicit and has_ph and str(payload.get("override_source") or "") == "reprotonate_at_ph":
        # A request rebuilt from its own manifest (DAG path) keeps its explicit-override label.
        explicit_raw, has_explicit, has_ph = ph_raw, True, False
    if not legacy and not has_explicit:
        if has_ph:
            raise ValueError(
                "a bare pH is ambiguous: the default MD export uses the docked microspecies. "
                "Pass --reprotonate-at-ph <pH> for an explicit override, or --protonation-policy openbabel_predicted with --ph"
            )
        return None, ""
    if legacy and not has_ph and not has_explicit:
        raise ValueError("--protonation-policy openbabel_predicted requires an explicit --ph")
    values: List[Tuple[str, float]] = []
    if has_explicit:
        values.append(("reprotonate_at_ph", _validate_ph(explicit_raw, "reprotonate_at_ph")))
    if legacy and has_ph:
        values.append(("legacy_openbabel_predicted", _validate_ph(ph_raw, "pH")))
    if len({round(value, 9) for _, value in values}) > 1:
        raise ValueError(f"conflicting_reprotonation_pH:{values[0][1]}!={values[1][1]}")
    return values[0][1], ("reprotonate_at_ph" if has_explicit else "legacy_openbabel_predicted")


class MDExportError(RuntimeError):
    def __init__(self, status: str, gate: str, reason: str):
        if status not in ROW_STATUSES:
            status = "failed"
        self.status = status
        self.gate = gate
        self.reason = reason
        super().__init__(reason)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_hash(payload: object) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def _safe_key(value: object) -> str:
    normalized = _SAFE_KEY.sub("_", str(value or "").strip()).strip("._")
    if not normalized:
        raise MDExportError("skipped_missing_configuration", "G1", "empty_output_identity")
    return normalized[:120]


def load_exact_tags(path: Path) -> List[str]:
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"tags_file_missing:{source}")
    if source.suffix.lower() == ".csv":
        frame = pd.read_csv(source, dtype=str).fillna("")
        if "tag" not in frame.columns:
            raise ValueError("tags CSV must contain a 'tag' column")
        raw = frame["tag"].tolist()
    else:
        raw = source.read_text(encoding="utf-8", errors="strict").splitlines()
    tags: List[str] = []
    seen: Set[str] = set()
    for value in raw:
        tag = str(value).strip()
        if not tag or tag.startswith("#"):
            continue
        if tag in seen:
            raise ValueError(f"duplicate_requested_tag:{tag}")
        seen.add(tag)
        tags.append(tag)
    if not tags:
        raise ValueError("tags file contains no tags")
    return tags


def _load_map(path: Optional[Path]) -> List[Dict[str, str]]:
    if path is None:
        return []
    frame = pd.read_csv(path, dtype=str).fillna("")
    return [{str(key): str(value).strip() for key, value in row.items()} for row in frame.to_dict(orient="records")]


def _match_map(rows: Sequence[Dict[str, str]], selected: Dict[str, object], file_field: str) -> Optional[Dict[str, str]]:
    tag = str(selected.get("tag") or "")
    ligand = str(selected.get("ligand") or "")
    protein = str(selected.get("protein") or selected.get("receptor") or "")
    exact = [row for row in rows if row.get("tag") == tag and row.get(file_field)]
    if len(exact) == 1:
        return exact[0]
    if len(exact) > 1:
        raise MDExportError("skipped_missing_configuration", "G2", f"ambiguous_{file_field}_map_for_tag:{tag}")
    candidates = [
        row for row in rows
        if row.get(file_field)
        and (not row.get("ligand") or row.get("ligand") == ligand)
        and (not row.get("protein") and not row.get("receptor") or row.get("protein") == protein or row.get("receptor") == protein)
    ]
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        raise MDExportError("skipped_missing_configuration", "G2", f"ambiguous_{file_field}_map_for_tag:{tag}")
    return None


def _resolve_mapped_path(raw: str, map_file: Path, project_dir: Path) -> Path:
    candidate = Path(str(raw).strip()).expanduser()
    if not candidate.is_absolute():
        beside_map = (map_file.parent / candidate).resolve()
        candidate = beside_map if beside_map.exists() else (project_dir / candidate).resolve()
    else:
        candidate = candidate.resolve()
    return candidate


def _resolve_selection_path(raw: object, project_dir: Path) -> Path:
    candidate = Path(str(raw or "").strip()).expanduser()
    return candidate.resolve() if candidate.is_absolute() else (project_dir / candidate).resolve()


def _selection_file(project_dir: Path) -> Path:
    return ensure_numbered_output_layout(project_dir)["post_scores_unified"] / "best_pose_per_tag_by_engine.csv"


def load_selected_rows(request: MDInputsRequest) -> List[Dict[str, object]]:
    selection_file = _selection_file(request.project_dir)
    if not selection_file.is_file():
        raise FileNotFoundError(f"best_pose_selection_missing:{selection_file}")
    frame = pd.read_csv(selection_file)
    required = {"engine", "tag", "protein", "ligand", "pose", "pose_file", "affinity_kcal_mol"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"best_pose_selection_missing_columns:{','.join(missing)}")
    frame["engine"] = frame["engine"].astype(str).str.strip().str.lower()
    frame["tag"] = frame["tag"].astype(str).str.strip()
    frame = frame[frame["engine"] == request.engine]
    tags = load_exact_tags(request.tags_file)
    rows: List[Dict[str, object]] = []
    for tag in tags:
        matches = frame[frame["tag"] == tag]
        if len(matches) != 1:
            rows.append({"tag": tag, "engine": request.engine, "_selection_error": f"selected_row_count:{len(matches)}"})
        else:
            rows.append(matches.iloc[0].to_dict())
    return rows


def _sdf_records(path: Path) -> List[str]:
    return [record for record in path.read_text(encoding="utf-8", errors="strict").split("$$$$") if record.strip()]


def _sdf_properties(path: Path, pose: int) -> Dict[str, str]:
    records = _sdf_records(path)
    if pose < 1 or pose > len(records):
        raise MDExportError("failed", "G1", f"pose_index_out_of_range:{pose}>{len(records)}")
    lines = records[pose - 1].splitlines()
    properties: Dict[str, str] = {}
    index = 0
    while index < len(lines):
        match = re.match(r"^>\s*<([^>]+)>", lines[index].strip())
        if not match:
            index += 1
            continue
        values: List[str] = []
        index += 1
        while index < len(lines) and lines[index].strip():
            values.append(lines[index].strip())
            index += 1
        properties[match.group(1)] = "\n".join(values)
    return properties


def _native_score(pose_file: Path, engine: str, pose: int) -> Tuple[Optional[float], str]:
    if engine == "gnina":
        properties = _sdf_properties(pose_file, pose)
        for field in ("minimizedAffinity", "vina_affinity", "affinity"):
            if field in properties:
                try:
                    return float(properties[field].split()[0]), field
                except ValueError:
                    raise MDExportError("failed", "G1", f"invalid_native_score:{field}")
        return None, "missing"
    parsed = parse_vina_pdbqt(pose_file)
    if parsed.empty:
        raise MDExportError("failed", "G1", "native_pose_score_parse_empty")
    match = parsed[pd.to_numeric(parsed["pose"], errors="coerce") == int(pose)]
    if len(match) != 1:
        raise MDExportError("failed", "G1", f"native_pose_score_row_count:{len(match)}")
    return float(match.iloc[0]["vina_affinity"]), "vina_affinity"


def _pdbqt_pose_coordinates(path: Path, pose: int) -> Tuple[List[Tuple[float, float, float]], List[str]]:
    coordinates: List[Tuple[float, float, float]] = []
    elements: List[str] = []
    current = 0
    has_models = False
    capture = pose == 1
    with path.open("r", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                has_models = True
                try:
                    current = int(line.split()[1])
                except (IndexError, ValueError):
                    current += 1
                capture = current == pose
                continue
            if line.startswith("ENDMDL"):
                if capture:
                    break
                capture = False
                continue
            if has_models and not capture:
                continue
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                parts = line.split()
                try:
                    x, y, z = float(parts[6]), float(parts[7]), float(parts[8])
                except (IndexError, ValueError) as exc:
                    raise MDExportError("failed", "G3", "invalid_pdbqt_coordinates") from exc
            token = (line.split()[-1] if line.split() else "").upper()
            element = _AUTODOCK_ELEMENTS.get(token)
            if not element:
                raw_element = line[76:78].strip().upper() if len(line) >= 78 else ""
                element = _AUTODOCK_ELEMENTS.get(raw_element)
            if not element:
                raise MDExportError("not_comparable", "G3", f"unknown_autodock_atom_type:{token}")
            if element.upper() == "H":
                continue
            coordinates.append((x, y, z))
            elements.append(element)
    if not coordinates:
        raise MDExportError("failed", "G1", f"selected_pose_has_no_heavy_atoms:{pose}")
    return coordinates, elements


def _parse_atom_map(raw: str, atom_count: int) -> List[int]:
    value = str(raw or "").strip()
    if not value:
        raise MDExportError("not_comparable", "G3", "missing_pose_to_topology_atom_map")
    try:
        parsed = json.loads(value)
        if not isinstance(parsed, list):
            raise ValueError
        mapping = [int(item) for item in parsed]
    except Exception:
        try:
            mapping = [int(item.strip()) for item in value.split(",") if item.strip()]
        except ValueError as exc:
            raise MDExportError("not_comparable", "G3", "invalid_pose_to_topology_atom_map") from exc
    if mapping and min(mapping) == 1 and max(mapping) == atom_count:
        mapping = [item - 1 for item in mapping]
    return _validate_atom_permutation(mapping, atom_count)


def _validate_atom_permutation(mapping: Sequence[int], atom_count: int) -> List[int]:
    """Pose heavy-atom index -> topology heavy-atom index; must be a 0-based permutation."""
    mapping = [int(item) for item in mapping]
    if len(mapping) != atom_count or sorted(mapping) != list(range(atom_count)):
        raise MDExportError("not_comparable", "G3", "pose_to_topology_atom_map_not_permutation")
    return mapping


def _derive_lineage_atom_map(
    pose_file: Path,
    pose: int,
    topology_file: Path,
    topology_pose: int,
) -> Tuple[List[Tuple[float, float, float]], List[str], List[int], Dict[str, object]]:
    """Spec 034 R1c: pose-to-topology atom_map from the Meeko SMILES/IDX lineage.

    Uses the same ``compare_lineage_graph_poses`` / selected mapping as the Spec 031 module.
    The returned permutation has the same semantics as ``_parse_atom_map``:
    ``atom_map[pose_heavy_index] == topology_heavy_index``.
    """
    try:
        lineage_pose = load_pdbqt_lineage_pose(pose_file, pose)
    except LineageError as exc:
        if exc.reason == "no_lineage_remarks":
            raise MDExportError("not_comparable", "G3", "missing_pose_to_topology_atom_map:no_lineage_remarks") from exc
        raise MDExportError("not_comparable", "G3", exc.reason) from exc
    try:
        topology_graph = load_sdf_graph_pose(topology_file, topology_pose)
    except (OSError, RuntimeError, ValueError) as exc:
        raise MDExportError("not_comparable", "G3", f"topology_not_comparable:{exc}") from exc
    result = compare_lineage_graph_poses(lineage_pose, topology_graph)
    if not result.comparable or result.selected_mapping is None:
        raise MDExportError("not_comparable", "G3", f"lineage_mapping_not_comparable:{result.reason}")
    atom_count = len(lineage_pose.elements)
    atom_map = _validate_atom_permutation(list(result.selected_mapping), atom_count)
    coordinates = [tuple(float(value) for value in row) for row in lineage_pose.coords]
    evidence = {
        "lineage_source": result.lineage_source,
        "reference_lineage_source": result.reference_lineage_source,
        "mapping_reason": result.reason,
        "mapped_heavy_atoms": int(result.mapped_heavy_atoms),
        "mapping_coverage": float(result.mapping_coverage),
        "valid_mapping_count": int(result.valid_mapping_count),
    }
    return coordinates, list(lineage_pose.elements), atom_map, evidence


def _charge_for_row(rows: Sequence[Dict[str, str]], selected: Dict[str, object]) -> Optional[int]:
    if not rows:
        return None
    match = _match_map(rows, selected, "net_charge")
    if match is None:
        raise MDExportError("skipped_missing_configuration", "G6", f"charge_map_missing_row:{selected.get('tag')}")
    try:
        return int(match["net_charge"])
    except ValueError as exc:
        raise MDExportError("skipped_missing_configuration", "G6", "invalid_net_charge") from exc


def _input_descriptor(path: Path) -> Dict[str, object]:
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise MDExportError("skipped_missing_configuration", "G2", f"input_file_missing:{source}")
    return {"path": str(source), "sha256": sha256_file(source), "size_bytes": source.stat().st_size}


def _materialize_receptor_pdb(source: Path, project_dir: Path) -> Tuple[Path, Dict[str, object]]:
    receptor = Path(source).expanduser().resolve()
    if receptor.suffix.lower() == ".pdb":
        return receptor, {"source_format": "pdb", "conversion": "none"}
    if receptor.suffix.lower() not in {".cif", ".mmcif"}:
        raise MDExportError("not_comparable", "G5", f"unsupported_receptor_format:{receptor.suffix.lower()}")
    try:
        from Bio.PDB import MMCIFParser, PDBIO
        import Bio
    except Exception as exc:
        raise MDExportError("skipped_missing_dependency", "G5", "biopython_required_for_mmcif") from exc
    cache_root = ensure_numbered_output_layout(project_dir)["post_metadata"] / "md_inputs_receptor_cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    source_hash = sha256_file(receptor)
    destination = cache_root / f"{source_hash}.pdb"
    if not destination.is_file():
        try:
            structure = MMCIFParser(QUIET=True).get_structure(receptor.stem, str(receptor))
            for model in structure:
                for chain in model:
                    if len(str(chain.id).strip()) != 1:
                        raise ValueError(f"mmcif_chain_not_pdb_compatible:{chain.id}")
                    for residue in chain:
                        for atom in residue:
                            if str(atom.get_altloc() or " ") != " ":
                                raise ValueError(f"unresolved_altloc:{atom.get_altloc()}")
            temporary = destination.with_suffix(".pdb.tmp")
            writer = PDBIO()
            writer.set_structure(structure)
            writer.save(str(temporary), preserve_atom_numbering=True)
            os.replace(temporary, destination)
        except MDExportError:
            raise
        except Exception as exc:
            raise MDExportError("not_comparable", "G5", f"mmcif_to_pdb_failed:{exc}") from exc
    return destination, {
        "source_format": receptor.suffix.lower().lstrip("."),
        "conversion": "biopython_mmcif_to_pdb",
        "conversion_backend": "biopython",
        "conversion_backend_version": str(getattr(Bio, "__version__", "unknown")),
        "materialized_pdb": str(destination),
        "materialized_pdb_sha256": sha256_file(destination),
    }


def _gate(status: str, reason: str = "", **evidence: object) -> Dict[str, object]:
    return {"status": status, "reason": reason, **evidence}


def _choose_ligand_chain(receptor_pdb: Path) -> str:
    used = {
        line[21]
        for line in receptor_pdb.read_text(encoding="utf-8", errors="strict").splitlines()
        if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21 and line[21] != " "
    }
    for candidate in "ZYXWVUTSRQPONMLKJIHGFEDCBA9876543210":
        if candidate not in used:
            return candidate
    raise MDExportError("failed", "G8", "no_available_single_character_ligand_chain")


def _ligand_pdb_line(serial: int, atom, chain: str) -> str:
    chars = [" "] * 80
    chars[0:6] = list("HETATM")
    chars[6:11] = list(f"{serial:5d}")
    chars[12:16] = list(f"{atom.name:>4s}")
    chars[16] = " "
    chars[17:20] = list("LIG")
    chars[21] = chain
    chars[22:26] = list(f"{1:4d}")
    chars[30:38] = list(f"{atom.x:8.3f}")
    chars[38:46] = list(f"{atom.y:8.3f}")
    chars[46:54] = list(f"{atom.z:8.3f}")
    chars[54:60] = list(f"{1.0:6.2f}")
    chars[60:66] = list(f"{0.0:6.2f}")
    chars[76:78] = list(f"{atom.element:>2s}"[-2:])
    return "".join(chars)


def _write_system_pdb(receptor: Path, ligand: PreparedLigand, output: Path) -> Tuple[int, str]:
    receptor_lines: List[str] = []
    serial = 0
    for line in receptor.read_text(encoding="utf-8", errors="strict").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        serial += 1
        if serial > 99999:
            raise MDExportError("failed", "G4", "pdb_atom_serial_overflow")
        padded = list(line.ljust(80)[:80])
        padded[6:11] = list(f"{serial:5d}")
        receptor_lines.append("".join(padded))
    receptor_count = serial
    chain = _choose_ligand_chain(receptor)
    ligand_serials: Dict[int, int] = {}
    ligand_lines: List[str] = []
    for atom in ligand.atoms:
        serial += 1
        if serial > 99999:
            raise MDExportError("failed", "G4", "pdb_atom_serial_overflow")
        ligand_serials[int(atom.index)] = serial
        ligand_lines.append(_ligand_pdb_line(serial, atom, chain))
    adjacency: Dict[int, Set[int]] = {value: set() for value in ligand_serials.values()}
    for bond in ligand.bonds:
        left = ligand_serials[int(bond.begin)]
        right = ligand_serials[int(bond.end)]
        adjacency[left].add(right)
        adjacency[right].add(left)
    conect: List[str] = []
    for source in sorted(adjacency):
        targets = sorted(adjacency[source])
        for start in range(0, len(targets), 4):
            conect.append("CONECT" + f"{source:5d}" + "".join(f"{target:5d}" for target in targets[start : start + 4]))
    output.write_text("\n".join(receptor_lines + [f"TER   {receptor_count + 1:5d}"] + ligand_lines + conect + ["END"]) + "\n", encoding="utf-8")
    return receptor_count, chain


def _write_ligand_only_pdb(ligand: PreparedLigand, output: Path, chain: str) -> None:
    serials = {int(atom.index): int(atom.index) for atom in ligand.atoms}
    atom_lines = [_ligand_pdb_line(int(atom.index), atom, chain) for atom in ligand.atoms]
    adjacency: Dict[int, Set[int]] = {serial: set() for serial in serials.values()}
    for bond in ligand.bonds:
        left, right = serials[int(bond.begin)], serials[int(bond.end)]
        adjacency[left].add(right)
        adjacency[right].add(left)
    conect: List[str] = []
    for source in sorted(adjacency):
        targets = sorted(adjacency[source])
        for start in range(0, len(targets), 4):
            conect.append("CONECT" + f"{source:5d}" + "".join(f"{target:5d}" for target in targets[start : start + 4]))
    output.write_text("\n".join(atom_lines + [f"TER   {len(atom_lines) + 1:5d}"] + conect + ["END"]) + "\n", encoding="utf-8")


def _validate_mol2_contract(path: Path, ligand: PreparedLigand) -> Dict[str, object]:
    lines = Path(path).read_text(encoding="utf-8", errors="strict").splitlines()
    section = ""
    atom_rows: List[List[str]] = []
    bond_rows: List[List[str]] = []
    for line in lines:
        if line.startswith("@<TRIPOS>"):
            section = line.strip().upper()
            continue
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if section == "@<TRIPOS>ATOM":
            atom_rows.append(line.split())
        elif section == "@<TRIPOS>BOND":
            bond_rows.append(line.split())
    errors: List[str] = []
    names = [row[1] for row in atom_rows if len(row) >= 2]
    expected_names = [atom.name for atom in ligand.atoms]
    if names != expected_names:
        errors.append("mol2_atom_names_do_not_match_deterministic_names")
    if len(atom_rows) != len(ligand.atoms):
        errors.append(f"mol2_atom_count_mismatch:{len(atom_rows)}!={len(ligand.atoms)}")
    if len(bond_rows) != len(ligand.bonds):
        errors.append(f"mol2_bond_count_mismatch:{len(bond_rows)}!={len(ligand.bonds)}")
    return {
        "status": "completed" if not errors else "failed",
        "is_valid": not errors,
        "errors": errors,
        "atom_count": len(atom_rows),
        "bond_count": len(bond_rows),
        "unique_atom_name_count": len(set(names)),
    }


def _outputs_valid_for_cache(provenance: Dict[str, object]) -> bool:
    outputs = provenance.get("outputs") or {}
    if not isinstance(outputs, dict) or not outputs:
        return False
    for record in outputs.values():
        if not isinstance(record, dict):
            return False
        path = Path(str(record.get("path") or ""))
        if not path.is_file() or sha256_file(path) != str(record.get("sha256") or ""):
            return False
    return True


def _atomic_publish(staging: Path, destination: Path) -> None:
    root = destination.parents[2].resolve()
    if root.name != "md_inputs" or root not in destination.resolve().parents:
        raise RuntimeError(f"unsafe_md_inputs_destination:{destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    backup = destination.with_name(destination.name + ".previous")
    if backup.exists():
        shutil.rmtree(backup)
    if destination.exists():
        destination.replace(backup)
    try:
        staging.replace(destination)
    except Exception:
        if backup.exists() and not destination.exists():
            backup.replace(destination)
        raise
    if backup.exists():
        shutil.rmtree(backup)


def _row_identity(row: Dict[str, object], engine: str) -> Tuple[str, str, str]:
    protein = _safe_key(row.get("protein") or row.get("receptor"))
    ligand = _safe_key(row.get("ligand"))
    return protein, ligand, _safe_key(engine)


def build_md_inputs_config_payload(request: MDInputsRequest) -> Dict[str, object]:
    """Normalize a request and fingerprint every currently resolvable dependency."""
    dependency_paths: Set[Path] = {_selection_file(request.project_dir), request.tags_file, request.receptor_map}
    project_manifest_file = manifest_path(request.project_dir)
    if project_manifest_file.is_file():
        dependency_paths.add(project_manifest_file)
    if request.topology_map:
        dependency_paths.add(request.topology_map)
    if request.charge_map:
        dependency_paths.add(request.charge_map)
    receptor_rows = _load_map(request.receptor_map)
    topology_rows = _load_map(request.topology_map)
    try:
        selected_rows = load_selected_rows(request)
        resolution_error = ""
    except Exception as exc:
        selected_rows = []
        resolution_error = f"{exc.__class__.__name__}:{exc}"
    for row in selected_rows:
        if row.get("_selection_error"):
            continue
        pose_file = _resolve_selection_path(row.get("pose_file"), request.project_dir)
        dependency_paths.add(pose_file)
        for record_file, _ in _ligand_provenance_records(request.project_dir, row.get("ligand")):
            dependency_paths.add(record_file)
        receptor_row = _match_map(receptor_rows, row, "receptor_file")
        if receptor_row:
            dependency_paths.add(_resolve_mapped_path(receptor_row["receptor_file"], request.receptor_map, request.project_dir))
        topology_row = _match_map(topology_rows, row, "topology_file") if request.topology_map else None
        if topology_row and request.topology_map:
            dependency_paths.add(_resolve_mapped_path(topology_row["topology_file"], request.topology_map, request.project_dir))
    dependencies: List[Dict[str, object]] = []
    for path in sorted(dependency_paths, key=lambda item: str(item)):
        dependencies.append(
            {"path": str(path), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
            if path.is_file()
            else {"path": str(path), "sha256": "", "missing": True}
        )
    payload = {
        "schema_version": SCHEMA_VERSION,
        "method_version": METHOD_VERSION,
        "request": request.to_dict(),
        "dependencies": dependencies,
        "chemistry_capability": openbabel_capability(),
        "resolution_error": resolution_error,
    }
    payload["request_fingerprint"] = _json_hash(payload)
    return payload


def _docking_provenance(request: MDInputsRequest, selected: Dict[str, object]) -> Dict[str, object]:
    try:
        manifest = load_manifest(request.project_dir)
    except Exception:
        manifest = {}
    settings = dict((manifest.get("engine_settings") or {}).get(request.engine, {}) or {})

    def first(*values: object) -> object:
        for value in values:
            if value is None:
                continue
            try:
                if bool(pd.isna(value)):
                    continue
            except (TypeError, ValueError):
                pass
            if str(value).strip():
                return value
        return "unknown"

    return {
        "engine": request.engine,
        "engine_version": first(settings.get("version"), settings.get("engine_version")),
        "executable_path": first(settings.get("executable_path"), settings.get("binary")),
        "container_digest": first(settings.get("container_digest"), settings.get("image_digest")),
        "scoring_mode": first(selected.get("scoring_mode"), settings.get("scoring_mode"), settings.get("cnn_scoring")),
        "seed": first(selected.get("seed"), settings.get("seed")),
        "exhaustiveness": first(selected.get("exhaustiveness"), settings.get("exhaustiveness")),
        "box": {
            key: first(selected.get(key))
            for key in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")
        },
    }


# ---------------------------------------------------------------- Spec 036 R1b: docked state

_LIGAND_SUFFIXES = (".pdbqt", ".sdf", ".mol2", ".mol", ".pdb")
_RECEPTOR_SUFFIXES = (".pdbqt",)
_PREPARATION_POLICIES = ("ph_model", "explicit_state", "as_input")


class MDMapError(RuntimeError):
    """An md-input map cannot be generated from project evidence; ``reason`` is machine-readable (R6)."""

    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


def _name_key(value: object, suffixes: Tuple[str, ...]) -> str:
    text = Path(str(value or "").strip()).name
    for suffix in suffixes:
        if text.lower().endswith(suffix):
            return text[: -len(suffix)]
    return text


def _ligand_provenance_records(project_dir: Path, ligand: object) -> List[Tuple[Path, Dict[str, object]]]:
    """Ligand preparation step reports (Spec 034 R2b) whose input or output name matches ``ligand``."""
    key = _name_key(ligand, _LIGAND_SUFFIXES)
    directory = ensure_numbered_output_layout(project_dir)["prep_ligands_prepared"] / "preparation_steps"
    records: List[Tuple[Path, Dict[str, object]]] = []
    if not key or not directory.is_dir():
        return records
    for report in sorted(directory.glob("*.json")):
        try:
            payload = json.loads(report.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if not isinstance(payload, dict):
            continue
        names = {_name_key(payload.get("input_file"), _LIGAND_SUFFIXES), _name_key(payload.get("output_file"), _LIGAND_SUFFIXES)}
        if key in names:
            records.append((report, payload))
    return records


def _project_policy_identity(project_dir: Path) -> Optional[Dict[str, object]]:
    block = load_protonation_policy(project_dir)
    if not block:
        return None
    return {
        "ligand_policy": block.get("ligand_policy"),
        "ligand_ph": block.get("ligand_ph"),
        "receptor_ph": block.get("receptor_ph"),
        "receptor_force_field": block.get("receptor_force_field"),
        "source": block.get("source"),
        "set_at": block.get("set_at"),
        "sha256": _json_hash(block),
    }


def _docked_state_for_row(request: MDInputsRequest, selected: Dict[str, object], pose_file: Path, pose: int) -> Dict[str, object]:
    """Spec 036 R1b: the docked microspecies of the selected pose, cross-checked with preparation provenance.

    Vina/Smina: the Meeko ``REMARK SMILES`` lineage (formal charges and hydrogens). GNINA: the pose SDF record.
    Both must agree with the ligand preparation step report (``protonation.prepared_state`` and, for
    PDBQT, ``protonation.pdbqt_state``). Missing or disagreeing evidence fails closed.
    """
    ligand = selected.get("ligand")
    records = _ligand_provenance_records(request.project_dir, ligand)
    key = _name_key(ligand, _LIGAND_SUFFIXES)
    if not records:
        raise MDExportError("skipped_missing_configuration", "G6", f"ligand_provenance_missing:{key}")
    if len(records) > 1:
        raise MDExportError("skipped_missing_configuration", "G6", f"ligand_provenance_ambiguous:{key}:{len(records)}")
    provenance_file, provenance = records[0]
    try:
        if request.engine == "gnina":
            state = docked_microspecies_from_sdf(pose_file, pose)
        else:
            state = docked_microspecies_from_pdbqt(pose_file, pose)
    except RuntimeError as exc:
        raise MDExportError("skipped_missing_dependency", "G6", str(exc)) from exc
    except ValueError as exc:
        raise MDExportError("not_comparable", "G6", f"docked_state_unavailable:{exc}") from exc
    protocol = provenance.get("protonation") if isinstance(provenance.get("protonation"), dict) else {}
    policy = str(protocol.get("policy") or "")
    if policy not in _PREPARATION_POLICIES:
        raise MDExportError("skipped_missing_configuration", "G6", f"preparation_policy_missing:{provenance_file.name}")
    prepared = protocol.get("prepared_state")
    if not isinstance(prepared, dict) or prepared.get("net_formal_charge") is None:
        raise MDExportError("skipped_missing_configuration", "G6", f"preparation_prepared_state_missing:{provenance_file.name}")
    docked_net = int(state["net_formal_charge"])
    prepared_net = int(prepared["net_formal_charge"])
    if docked_net != prepared_net:
        raise MDExportError("failed", "G6", f"docked_state_differs_from_preparation_provenance:net_charge:{docked_net}!={prepared_net}")
    if request.engine != "gnina":
        pdbqt_state = protocol.get("pdbqt_state")
        if not isinstance(pdbqt_state, dict) or not pdbqt_state.get("microspecies_smiles"):
            raise MDExportError("skipped_missing_configuration", "G6", f"preparation_pdbqt_state_missing:{provenance_file.name}")
        if str(pdbqt_state["microspecies_smiles"]) != str(state["canonical_smiles"]):
            raise MDExportError("failed", "G6", "docked_state_differs_from_preparation_provenance:microspecies")
    else:
        prepared_smiles = str(prepared.get("microspecies_smiles") or "")
        if prepared_smiles and state["canonical_smiles"] and prepared_smiles != str(state["canonical_smiles"]):
            raise MDExportError("failed", "G6", "docked_state_differs_from_preparation_provenance:microspecies")
    state.update(
        {
            "ligand_policy": policy,
            "preparation_ph": protocol.get("ph") if protocol.get("ph_model_run") else None,
            "ligand_preparation_record": {"path": str(provenance_file), "sha256": sha256_file(provenance_file)},
            "charge_authority": "predicted" if policy == "ph_model" else "docked_state",
            "scientific_review": "human_review_required" if policy in {"ph_model", "as_input"} else "completed",
        }
    )
    return state


def _docked_state_best_effort(request: MDInputsRequest, selected: Dict[str, object], pose_file: Path, pose: int):
    """Override rows record the docked state when it can be read, without making it a precondition."""
    try:
        return _docked_state_for_row(request, selected, pose_file, pose), ""
    except MDExportError as exc:
        return None, f"{exc.gate}:{exc.reason}"


def _check_row_against_project_policy(project_policy: Optional[Dict[str, object]], docked_state: Dict[str, object]) -> None:
    if not project_policy:
        return
    stored = project_policy.get("ligand_policy")
    docked = docked_state.get("ligand_policy")
    if stored != docked:
        raise MDExportError(
            "failed",
            "G6",
            f"protonation_policy_conflict:ligand_policy:project={stored}!=docked_state={docked}",
        )


def _neutral_parent_pose(pose: GraphPose) -> GraphPose:
    """Charge-blind copy of a heavy-atom graph: each charged atom is taken to its neutral parent.

    The Spec 031 node rule compares formal charges, which differ between protonation states of
    the same molecule. Removing the charge and one proton per unit of charge gives the parent
    heavy-atom graph that the crystal topology describes. The Spec 031 rules are unchanged.
    """
    graph = pose.graph.copy()
    for node in graph.nodes:
        data = graph.nodes[node]
        charge = int(data.get("formal_charge", 0))
        if "total_h" in data:
            data["total_h"] = int(data["total_h"]) - charge
        data["formal_charge"] = 0
    return dataclasses.replace(pose, graph=graph)


def _derive_charge_blind_atom_map(
    pose_file: Path,
    pose: int,
    topology_file: Path,
    topology_pose: int,
) -> Tuple[List[Tuple[float, float, float]], List[str], List[int], Dict[str, object]]:
    """Spec 036 R1b G3 for a docked microspecies: the Spec 034 lineage map on charge-blind parent graphs."""
    try:
        lineage_pose = load_pdbqt_lineage_pose(pose_file, pose)
    except LineageError as exc:
        if exc.reason == "no_lineage_remarks":
            raise MDExportError("not_comparable", "G3", "missing_pose_to_topology_atom_map:no_lineage_remarks") from exc
        raise MDExportError("not_comparable", "G3", exc.reason) from exc
    try:
        topology_graph = load_sdf_graph_pose(topology_file, topology_pose)
    except (OSError, RuntimeError, ValueError) as exc:
        raise MDExportError("not_comparable", "G3", f"topology_not_comparable:{exc}") from exc
    result = compare_lineage_graph_poses(_neutral_parent_pose(lineage_pose), _neutral_parent_pose(topology_graph))
    if not result.comparable or result.selected_mapping is None:
        raise MDExportError("not_comparable", "G3", f"lineage_mapping_not_comparable:{result.reason}")
    atom_count = len(lineage_pose.elements)
    atom_map = _validate_atom_permutation(list(result.selected_mapping), atom_count)
    coordinates = [tuple(float(value) for value in row) for row in lineage_pose.coords]
    evidence = {
        "lineage_source": result.lineage_source,
        "reference_lineage_source": result.reference_lineage_source,
        "mapping_reason": result.reason,
        "mapping_coverage": float(result.mapping_coverage),
        "valid_mapping_count": int(result.valid_mapping_count),
        "mapping_charge_policy": _CHARGE_BLIND_MAPPING,
    }
    return coordinates, list(lineage_pose.elements), atom_map, evidence


def _prepare_docked_ligand(backend, request: MDInputsRequest, pose_file: Path, pose: int, docked_state: Dict[str, object],
                           pose_coords, pose_elements):
    try:
        if request.engine == "gnina":
            return backend.prepare_docked_state(
                pose_file,
                source_pose=pose,
                expected_net_charge=int(docked_state["net_formal_charge"]),
                add_hydrogens=request.protonate,
            )
        return backend.prepare_docked_state(
            pose_file,
            source_pose=pose,
            expected_net_charge=int(docked_state["net_formal_charge"]),
            add_hydrogens=request.protonate,
            pose_heavy_coordinates=pose_coords,
            pose_heavy_elements=pose_elements,
            pose_to_smiles=docked_state["pose_to_smiles"],
            microspecies_smiles=str(docked_state["smiles"]),
        )
    except MDExportError:
        raise
    except ValueError as exc:
        raise MDExportError("failed", "G6", str(exc)) from exc


def _docked_gate_evidence(docked_state: Optional[Dict[str, object]], reason: str = "") -> Optional[Dict[str, object]]:
    if docked_state is None:
        return {"status": "unavailable", "reason": reason} if reason else None
    return {
        "source": docked_state["source"],
        "net_formal_charge": int(docked_state["net_formal_charge"]),
        "microspecies_smiles": docked_state["canonical_smiles"],
        "ligand_policy": docked_state.get("ligand_policy"),
        "preparation_record": docked_state.get("ligand_preparation_record"),
        "mapping_charge_policy": _CHARGE_BLIND_MAPPING if docked_state["source"] == "meeko_remark_smiles" else "not_applicable",
    }


def _protonation_fingerprint(request: MDInputsRequest, docked_state: Optional[Dict[str, object]], project_policy) -> Dict[str, object]:
    return {
        "mode": "reprotonation_override" if request.reprotonation_override else "docked_state",
        "override_pH": request.pH,
        "override_source": request.override_source,
        "docked_net_charge": int(docked_state["net_formal_charge"]) if docked_state else None,
        "docked_smiles": str(docked_state["canonical_smiles"]) if docked_state else "",
        "docked_ligand_policy": str(docked_state.get("ligand_policy")) if docked_state else "",
        "ligand_preparation_sha256": str((docked_state or {}).get("ligand_preparation_record", {}).get("sha256", "")),
        "project_policy_sha256": str((project_policy or {}).get("sha256", "")),
    }


def check_request_protonation_against_project(project_dir: Path, request: MDInputsRequest) -> Optional[Dict[str, object]]:
    """Explicit re-protonation values must agree with a stored project policy (R1a)."""
    block = load_protonation_policy(project_dir)
    if not block or not request.reprotonation_override:
        return block
    stored_policy = str(block.get("ligand_policy") or "")
    requested = f"reprotonate at pH {request.pH:.2f} ({request.override_source})"
    if stored_policy != "ph_model":
        raise ProtonationPolicyConflict("ligand_policy", stored_policy, requested)
    stored_ph = block.get("ligand_ph")
    if stored_ph is None or abs(float(stored_ph) - float(request.pH)) > 1e-9:
        raise ProtonationPolicyConflict("ligand_ph", stored_ph, f"{request.pH:.2f} ({request.override_source})")
    return block


# ---------------------------------------------------------------- Spec 036 R6: auto md-input maps


def _best_pose_rows(project_dir: Path, engine: str) -> "pd.DataFrame":
    selection_file = _selection_file(project_dir)
    if not selection_file.is_file():
        raise MDMapError(f"best_pose_selection_missing:{selection_file}")
    frame = pd.read_csv(selection_file)
    required = {"engine", "tag", "protein", "ligand", "pose", "pose_file"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise MDMapError(f"best_pose_selection_missing_columns:{','.join(missing)}")
    frame["engine"] = frame["engine"].astype(str).str.strip().str.lower()
    frame["tag"] = frame["tag"].astype(str).str.strip()
    rows = frame[frame["engine"] == engine]
    if rows.empty:
        raise MDMapError(f"no_best_pose_rows_for_engine:{engine}")
    duplicated = rows["tag"][rows["tag"].duplicated()].tolist()
    if duplicated:
        raise MDMapError(f"duplicate_best_pose_tag:{duplicated[0]}")
    return rows


def _receptor_lineage_for(project_dir: Path, protein: object) -> Tuple[Path, Dict[str, object]]:
    """Strict receptor lineage (Spec 033): prepared PDBQT -> ``.preparation.json`` -> the chain PDB it was prepared from."""
    key = _name_key(protein, _RECEPTOR_SUFFIXES)
    prepared_dir = ensure_numbered_output_layout(project_dir)["prep_receptors_prepared"]
    candidates = sorted(item for item in prepared_dir.glob("*.pdbqt") if item.stem == key) if key else []
    if not candidates:
        raise MDMapError(f"receptor_lineage_missing:{key}")
    if len(candidates) > 1:
        raise MDMapError(f"receptor_lineage_ambiguous:{key}:{len(candidates)}")
    prepared = candidates[0]
    record = prepared.parent / (prepared.name + ".preparation.json")
    if not record.is_file():
        raise MDMapError(f"receptor_preparation_record_missing:{record.name}")
    try:
        payload = json.loads(record.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise MDMapError(f"receptor_preparation_record_unreadable:{record.name}") from exc
    source = Path(str(payload.get("input_file") or "")).expanduser()
    if not str(source) or not source.is_file():
        raise MDMapError(f"receptor_source_missing:{source}")
    if source.suffix.lower() not in {".pdb", ".cif", ".mmcif"}:
        raise MDMapError(f"receptor_source_format_unsupported:{source.suffix.lower()}")
    source = source.resolve()
    return source, {
        "protein": key,
        "prepared_pdbqt": str(prepared),
        "preparation_record": str(record),
        "preparation_record_sha256": sha256_file(record),
        "receptor_source": str(source),
        "receptor_source_sha256": sha256_file(source),
    }


def _ligand_topology_for(project_dir: Path, ligand: object) -> Tuple[Path, Dict[str, object]]:
    """Connectivity-bearing topology (Spec 032 G3): the original input SDF recorded in the ligand provenance."""
    key = _name_key(ligand, _LIGAND_SUFFIXES)
    records = _ligand_provenance_records(project_dir, ligand)
    if not records:
        raise MDMapError(f"ligand_provenance_missing:{key}")
    if len(records) > 1:
        raise MDMapError(f"ligand_provenance_ambiguous:{key}:{len(records)}")
    report, payload = records[0]
    source = Path(str(payload.get("input_file") or "")).expanduser()
    if not str(source) or not source.is_file():
        raise MDMapError(f"ligand_topology_source_missing:{source}")
    if source.suffix.lower() not in {".sdf", ".mol", ".mol2"}:
        raise MDMapError(f"ligand_topology_format_unsupported:{source.suffix.lower()}")
    source = source.resolve()
    return source, {
        "ligand": key,
        "preparation_record": str(report),
        "preparation_record_sha256": sha256_file(report),
        "topology_source": str(source),
        "topology_source_sha256": sha256_file(source),
    }


def resolve_md_input_maps(
    project_dir: Path,
    engine: str,
    *,
    tags_file: Optional[Path] = None,
    receptor_map: Optional[Path] = None,
    topology_map: Optional[Path] = None,
    auto_maps: bool = False,
) -> Dict[str, object]:
    """Spec 036 R6: generate the md-input maps that the caller did not supply.

    Generated maps are written to ``.meta/md_inputs_maps/`` with the source hashes in
    ``auto_maps_manifest.json``. A missing or ambiguous source fails with an explicit reason.
    ``auto_maps`` generates every map and refuses explicit map files as ambiguous.
    """
    root = Path(project_dir).expanduser().resolve()
    if auto_maps and any(value is not None for value in (tags_file, receptor_map, topology_map)):
        raise MDMapError("auto_maps_conflicts_with_explicit_map_file")
    need_tags = auto_maps or tags_file is None
    need_receptor = auto_maps or receptor_map is None
    need_topology = (auto_maps or topology_map is None) and engine != "gnina"
    sources: Dict[str, str] = {
        "tags": "explicit" if tags_file is not None else "auto_generated",
        "receptor_map": "explicit" if receptor_map is not None else "auto_generated",
        "topology_map": (
            "explicit" if topology_map is not None else ("auto_generated" if need_topology else "not_required_for_gnina")
        ),
    }
    if not (need_tags or need_receptor or need_topology):
        return {
            "tags_file": str(tags_file),
            "receptor_map": str(receptor_map),
            "topology_map": str(topology_map) if topology_map is not None else "",
            "sources": sources,
            "manifest_file": "",
        }
    maps_dir = ensure_numbered_output_layout(root)["meta"] / AUTO_MAPS_DIRNAME
    maps_dir.mkdir(parents=True, exist_ok=True)
    rows = _best_pose_rows(root, engine)
    if need_tags:
        tag_values = [str(tag) for tag in rows["tag"].tolist()]
    else:
        tag_values = load_exact_tags(tags_file)
        unknown = [tag for tag in tag_values if tag not in set(rows["tag"].tolist())]
        if unknown:
            raise MDMapError(f"tag_not_in_best_pose_table:{unknown[0]}")
    selected = rows[rows["tag"].isin(tag_values)].set_index("tag").loc[tag_values].reset_index()
    generated: Dict[str, Dict[str, object]] = {}
    evidence: Dict[str, object] = {"receptors": [], "ligands": []}
    if need_tags:
        path = maps_dir / "tags.csv"
        pd.DataFrame({"tag": tag_values}).to_csv(path, index=False)
        generated["tags"] = {"path": str(path), "rows": len(tag_values)}
    if need_receptor:
        receptor_rows = []
        for _, row in selected.iterrows():
            source, record = _receptor_lineage_for(root, row["protein"])
            receptor_rows.append({"tag": row["tag"], "protein": record["protein"], "receptor_file": str(source)})
            evidence["receptors"].append(record)
        path = maps_dir / "receptor_map.csv"
        pd.DataFrame(receptor_rows, columns=["tag", "protein", "receptor_file"]).to_csv(path, index=False)
        generated["receptor_map"] = {"path": str(path), "rows": len(receptor_rows)}
    if need_topology:
        topology_rows = []
        for _, row in selected.iterrows():
            source, record = _ligand_topology_for(root, row["ligand"])
            topology_rows.append(
                {
                    "tag": row["tag"],
                    "ligand": record["ligand"],
                    "topology_file": str(source),
                    "topology_pose": 1,
                    "atom_map": "",
                }
            )
            evidence["ligands"].append(record)
        path = maps_dir / "topology_map.csv"
        pd.DataFrame(
            topology_rows, columns=["tag", "ligand", "topology_file", "topology_pose", "atom_map"]
        ).to_csv(path, index=False)
        generated["topology_map"] = {"path": str(path), "rows": len(topology_rows)}
    for entry in generated.values():
        entry["sha256"] = sha256_file(Path(str(entry["path"])))
    manifest_file = maps_dir / "auto_maps_manifest.json"
    manifest = {
        "schema_version": AUTO_MAPS_SCHEMA_VERSION,
        "generated_at": _utc_now(),
        "engine": engine,
        "maps": generated,
        "sources": evidence,
        "best_pose_selection": {"path": str(_selection_file(root)), "sha256": sha256_file(_selection_file(root))},
    }
    manifest_file.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "tags_file": str(generated["tags"]["path"]) if need_tags else str(tags_file),
        "receptor_map": str(generated["receptor_map"]["path"]) if need_receptor else str(receptor_map),
        "topology_map": (
            str(generated["topology_map"]["path"]) if need_topology else (str(topology_map) if topology_map else "")
        ),
        "sources": sources,
        "manifest_file": str(manifest_file),
        "maps": generated,
    }


def export_md_inputs(
    request: MDInputsRequest,
    *,
    chemistry_backend: Optional[ChemistryBackend] = None,
) -> Dict[str, object]:
    numbered = ensure_numbered_output_layout(request.project_dir)
    output_root = numbered["analysis_root_numbered"] / "md_inputs"
    output_root.mkdir(parents=True, exist_ok=True)
    selection_file = _selection_file(request.project_dir)
    receptor_rows = _load_map(request.receptor_map)
    topology_rows = _load_map(request.topology_map)
    charge_rows = _load_map(request.charge_map)
    check_request_protonation_against_project(request.project_dir, request)
    project_policy = _project_policy_identity(request.project_dir)
    backend = chemistry_backend if chemistry_backend is not None else get_openbabel_backend()
    capability = (
        {"status": "completed", "backend": backend.name, "backend_version": backend.version, "reason": ""}
        if backend is not None
        else openbabel_capability()
    )
    config_payload = build_md_inputs_config_payload(request)
    config_file = numbered["meta"] / "md_inputs_config.json"
    config_file.parent.mkdir(parents=True, exist_ok=True)
    config_serialized = json.dumps(config_payload, indent=2, sort_keys=True)
    if not config_file.exists() or config_file.read_text(encoding="utf-8", errors="replace") != config_serialized:
        config_file.write_text(config_serialized, encoding="utf-8")
    results: List[Dict[str, object]] = []

    try:
        selected_rows = load_selected_rows(request)
    except Exception as exc:
        selected_rows = [{"tag": "", "engine": request.engine, "_selection_error": str(exc)}]

    for selected in selected_rows:
        tag = str(selected.get("tag") or "")
        gates: Dict[str, Dict[str, object]] = {}
        row_result: Dict[str, object] = {
            "tag": tag,
            "engine": request.engine,
            "status": "failed",
            "reason": "",
            "gates": gates,
            "cache_reused": False,
        }
        try:
            if selected.get("_selection_error"):
                raise MDExportError("skipped_missing_configuration", "G1", str(selected["_selection_error"]))
            pose = int(pd.to_numeric(selected.get("pose"), errors="raise"))
            if pose < 1:
                raise MDExportError("failed", "G1", f"invalid_pose_number:{pose}")
            pose_file = _resolve_selection_path(selected.get("pose_file"), request.project_dir)
            if not pose_file.is_file():
                raise MDExportError("skipped_missing_configuration", "G1", f"pose_file_missing:{pose_file}")
            expected_score = float(selected.get("affinity_kcal_mol"))
            native_score, native_field = _native_score(pose_file, request.engine, pose)
            if native_score is not None and abs(native_score - expected_score) >= 0.02:
                raise MDExportError("failed", "G1", f"native_score_mismatch:{native_score}:{expected_score}")
            gates["G1"] = _gate(
                "completed",
                pose=pose,
                expected_affinity_kcal_mol=expected_score,
                native_score=native_score,
                native_score_field=native_field,
                tolerance=0.02,
            )

            receptor_row = _match_map(receptor_rows, selected, "receptor_file")
            if receptor_row is None:
                raise MDExportError("skipped_missing_configuration", "G2", f"receptor_map_missing_row:{tag}")
            receptor_file = _resolve_mapped_path(receptor_row["receptor_file"], request.receptor_map, request.project_dir)
            pose_input = _input_descriptor(pose_file)
            receptor_input = _input_descriptor(receptor_file)

            topology_file = pose_file
            topology_pose = pose
            topology_row: Optional[Dict[str, str]] = None
            pose_coords: Optional[List[Tuple[float, float, float]]] = None
            pose_elements: Optional[List[str]] = None
            atom_map: Optional[List[int]] = None
            atom_map_source = "native_sdf_order"
            lineage_evidence: Optional[Dict[str, object]] = None
            if request.engine != "gnina":
                if request.topology_map is None:
                    raise MDExportError("not_comparable", "G3", "pdbqt_requires_explicit_topology_map")
                topology_row = _match_map(topology_rows, selected, "topology_file")
                if topology_row is None:
                    raise MDExportError("not_comparable", "G3", f"topology_map_missing_row:{tag}")
                topology_file = _resolve_mapped_path(topology_row["topology_file"], request.topology_map, request.project_dir)
                topology_pose = int(topology_row.get("topology_pose") or 1)
                explicit_atom_map = str(topology_row.get("atom_map") or "").strip()
                if explicit_atom_map:
                    # An explicit atom_map always wins over lineage derivation.
                    pose_coords, pose_elements = _pdbqt_pose_coordinates(pose_file, pose)
                    atom_map = _parse_atom_map(explicit_atom_map, len(pose_coords))
                    atom_map_source = "explicit_topology_map"
                elif request.reprotonation_override:
                    # Spec 034 R1c: derive the permutation from the Meeko SMILES/IDX lineage (predicted override path).
                    pose_coords, pose_elements, atom_map, lineage_evidence = _derive_lineage_atom_map(
                        pose_file, pose, topology_file, topology_pose
                    )
                    atom_map_source = LINEAGE_METHOD
                else:
                    # Spec 036 R1b: docked microspecies. The Spec 034 lineage map is taken on charge-blind parent graphs.
                    pose_coords, pose_elements, atom_map, lineage_evidence = _derive_charge_blind_atom_map(
                        pose_file, pose, topology_file, topology_pose
                    )
                    atom_map_source = LINEAGE_METHOD
            elif pose_file.suffix.lower() != ".sdf":
                raise MDExportError("not_comparable", "G3", "gnina_pose_source_must_be_sdf")
            topology_input = _input_descriptor(topology_file)
            gates["G2"] = _gate(
                "completed",
                pose_file=pose_input,
                topology_file=topology_input,
                receptor_file=receptor_input,
                selection_file=_input_descriptor(selection_file),
            )
            gates["G3"] = _gate(
                "completed",
                topology_source="native_pose_sdf" if request.engine == "gnina" else "explicit_topology_map",
                topology_pose=topology_pose,
                atom_mapping="native_sdf_order" if request.engine == "gnina" else "explicit_pose_to_topology_permutation",
                atom_map=atom_map,
                atom_map_source=atom_map_source,
                lineage=lineage_evidence,
            )

            receptor_pdb, receptor_materialization = _materialize_receptor_pdb(receptor_file, request.project_dir)
            receptor_validation = validate_md_receptor_pdb(receptor_pdb)
            if not receptor_validation["is_valid"]:
                raise MDExportError(str(receptor_validation["status"]), "G5", ";".join(receptor_validation["errors"]))
            gates["G5"] = _gate(
                "completed",
                **{key: value for key, value in receptor_validation.items() if key not in {"status", "is_valid"}},
                **receptor_materialization,
            )

            if request.reprotonation_override:
                docked_state, docked_reason = _docked_state_best_effort(request, selected, pose_file, pose)
            else:
                docked_state = _docked_state_for_row(request, selected, pose_file, pose)
                docked_reason = ""
                _check_row_against_project_policy(project_policy, docked_state)
            if backend is None:
                raise MDExportError("skipped_missing_dependency", "G6", str(capability.get("reason") or "openbabel_missing"))
            expected_charge = _charge_for_row(charge_rows, selected)
            if expected_charge is not None:
                charge_authority = "human_approved"
                review_status = "completed"
            elif request.reprotonation_override:
                charge_authority = "predicted"
                review_status = "human_review_required"
            else:
                charge_authority = str(docked_state["charge_authority"])
                review_status = str(docked_state["scientific_review"])
            docking_provenance = _docking_provenance(request, selected)
            project_manifest_file = manifest_path(request.project_dir)
            project_manifest_input = (
                _input_descriptor(project_manifest_file)
                if project_manifest_file.is_file()
                else {"path": str(project_manifest_file), "sha256": "", "status": "unknown_missing"}
            )
            protein_key, ligand_key, engine_key = _row_identity(selected, request.engine)
            destination = output_root / protein_key / ligand_key / engine_key
            fingerprint_payload = {
                "method_version": METHOD_VERSION,
                "request": {key: value for key, value in request.to_dict().items() if key != "force"},
                "tag": tag,
                "pose": pose,
                "expected_score": expected_score,
                "inputs": {
                    "pose": pose_input,
                    "topology": topology_input,
                    "receptor": receptor_input,
                    "project_manifest": project_manifest_input,
                },
                "backend": {"name": backend.name, "version": backend.version},
                "expected_charge": expected_charge,
                "protonation": _protonation_fingerprint(request, docked_state, project_policy),
                "docking_provenance": docking_provenance,
                "atom_map": atom_map,
                "atom_map_source": atom_map_source,
            }
            row_fingerprint = _json_hash(fingerprint_payload)
            existing_provenance = destination / "provenance.json"
            if not request.force and existing_provenance.is_file():
                try:
                    cached = json.loads(existing_provenance.read_text(encoding="utf-8"))
                except Exception:
                    cached = {}
                if cached.get("row_fingerprint") == row_fingerprint and _outputs_valid_for_cache(cached):
                    row_result.update(
                        {
                            "status": "completed",
                            "reason": "cache_reused_content_hash_match",
                            "cache_reused": True,
                            "output_dir": str(destination),
                            "provenance_file": str(existing_provenance),
                            "row_fingerprint": row_fingerprint,
                            "gates": cached.get("gates", gates),
                        }
                    )
                    results.append(row_result)
                    continue
            if request.reprotonation_override:
                ligand = backend.prepare(
                    topology_file,
                    source_pose=topology_pose,
                    pH=request.pH,
                    protonate=request.protonate,
                    pose_heavy_coordinates=pose_coords,
                    pose_heavy_elements=pose_elements,
                    pose_to_topology=atom_map,
                )
            else:
                ligand = _prepare_docked_ligand(backend, request, pose_file, pose, docked_state, pose_coords, pose_elements)
                docked_net = int(docked_state["net_formal_charge"])
                if int(ligand.net_charge) != docked_net:
                    raise MDExportError("failed", "G6", f"md_state_differs_from_docked_state:{int(ligand.net_charge)}!={docked_net}")
            if expected_charge is not None and int(ligand.net_charge) != int(expected_charge):
                raise MDExportError("failed", "G6", f"approved_charge_mismatch:{ligand.net_charge}!={expected_charge}")
            if not request.protonate and not any(atom.element.upper() == "H" for atom in ligand.atoms):
                raise MDExportError("failed", "G6", "no_protonate_requires_explicit_hydrogens")
            ligand_ph = request.pH if request.reprotonation_override else docked_state.get("preparation_ph")
            ph_source = (
                request.override_source
                if request.reprotonation_override
                else ("preparation_provenance" if ligand_ph is not None else "not_applicable_no_ph_model")
            )
            override_evidence: Optional[Dict[str, object]] = None
            if request.reprotonation_override:
                docked_net_known = int(docked_state["net_formal_charge"]) if docked_state else None
                override_evidence = {
                    "pH": request.pH,
                    "source": request.override_source,
                    "docked_net_charge": docked_net_known,
                    "docked_state_unavailable": docked_reason or None,
                    "md_net_charge": int(ligand.net_charge),
                    "charge_changed_from_docked": (
                        int(ligand.net_charge) != docked_net_known if docked_net_known is not None else None
                    ),
                }
            gates["G6"] = _gate(
                "completed",
                pH=ligand_ph,
                pH_source=ph_source,
                protonated=request.protonate,
                docked_state=_docked_gate_evidence(docked_state, docked_reason),
                reprotonation_override=override_evidence,
                net_charge=int(ligand.net_charge),
                charge_authority=charge_authority,
                scientific_review=review_status,
                implicit_hydrogen_count=int(ligand.implicit_hydrogen_count),
                heavy_coordinate_max_delta=float(ligand.heavy_coordinate_max_delta),
                chemistry_backend=ligand.backend,
                chemistry_backend_version=ligand.backend_version,
            )
            if len({atom.name for atom in ligand.atoms}) != len(ligand.atoms):
                raise MDExportError("failed", "G7", "duplicate_ligand_atom_names")
            gates["G7"] = _gate(
                "completed",
                ligand_atom_count=len(ligand.atoms),
                ligand_bond_count=len(ligand.bonds),
                unique_atom_name_count=len({atom.name for atom in ligand.atoms}),
            )

            destination.parent.mkdir(parents=True, exist_ok=True)
            staging_parent = Path(tempfile.mkdtemp(prefix=f".{engine_key}_", dir=str(destination.parent)))
            staging = staging_parent / engine_key
            staging.mkdir()
            try:
                shutil.copy2(receptor_pdb, staging / "receptor.pdb")
                written = backend.write(ligand, staging, request.ligand_formats)
                receptor_count, ligand_chain = _write_system_pdb(receptor_pdb, ligand, staging / "system.pdb")
                if "pdb" in written:
                    _write_ligand_only_pdb(ligand, written["pdb"], ligand_chain)
                mol2_validation = _validate_mol2_contract(written["mol2"], ligand)
                if not mol2_validation["is_valid"]:
                    raise MDExportError("failed", "G7", ";".join(mol2_validation["errors"]))
                system_validation = validate_md_system_pdb(
                    staging / "system.pdb",
                    expected_receptor_atoms=receptor_count,
                    expected_ligand_atoms=len(ligand.atoms),
                    ligand_chain=ligand_chain,
                    expected_ligand_bonds=[(bond.begin, bond.end) for bond in ligand.bonds],
                )
                if not system_validation["is_valid"]:
                    validation_errors = [str(error) for error in system_validation["errors"]]
                    g4_errors = [
                        error for error in validation_errors
                        if "atom_count_mismatch" in error or "duplicate_atom_serial" in error
                    ]
                    g7_errors = [error for error in validation_errors if error.startswith("conect_")]
                    if g4_errors:
                        raise MDExportError("failed", "G4", ";".join(g4_errors))
                    if g7_errors:
                        raise MDExportError("failed", "G7", ";".join(g7_errors))
                    raise MDExportError("failed", "G8", ";".join(validation_errors))
                gates["G4"] = _gate(
                    "completed",
                    receptor_atom_count=receptor_count,
                    ligand_atom_count=len(ligand.atoms),
                    system_atom_count=receptor_count + len(ligand.atoms),
                )
                gates["G7"].update(
                    {
                        "mol2_validation": mol2_validation,
                        "conect_edge_count": system_validation["conect_edge_count"],
                    }
                )
                gates["G8"] = _gate("completed", **{key: value for key, value in system_validation.items() if key not in {"status", "is_valid"}})
                output_paths = {"receptor_pdb": staging / "receptor.pdb", "system_pdb": staging / "system.pdb"}
                output_paths.update({f"ligand_{fmt}": path for fmt, path in written.items()})
                provenance = {
                    "schema_version": PROVENANCE_SCHEMA_VERSION,
                    "method_version": METHOD_VERSION,
                    "generated_at": _utc_now(),
                    "status": "completed",
                    "tag": tag,
                    "engine": request.engine,
                    "protein": str(selected.get("protein") or ""),
                    "ligand": str(selected.get("ligand") or ""),
                    "pose_no": pose,
                    "affinity_kcal_mol": expected_score,
                    "native_score": native_score,
                    "native_score_field": native_field,
                    "pH": ligand_ph,
                    "pH_source": ph_source,
                    "net_charge": int(ligand.net_charge),
                    "charge_authority": charge_authority,
                    "scientific_review": review_status,
                    "protonation_policy": project_policy,
                    "docked_state": _docked_gate_evidence(docked_state, docked_reason),
                    "reprotonation_override": override_evidence,
                    "consumer_profile": request.consumer_profile,
                    "force_field": "unknown_not_parameterized",
                    "receptor_protonation_performed": False,
                    "row_fingerprint": row_fingerprint,
                    "inputs": fingerprint_payload["inputs"],
                    "chemistry_backend": fingerprint_payload["backend"],
                    "docking_provenance": docking_provenance,
                    "atom_map_source": atom_map_source,
                    "gates": gates,
                    "outputs": {},
                }
                for key, path in output_paths.items():
                    provenance["outputs"][key] = {"path": str(destination / path.name), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
                (staging / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True), encoding="utf-8")
                _atomic_publish(staging, destination)
            finally:
                if staging_parent.exists():
                    shutil.rmtree(staging_parent, ignore_errors=True)
            row_result.update(
                {
                    "status": "completed",
                    "reason": "",
                    "output_dir": str(destination),
                    "provenance_file": str(destination / "provenance.json"),
                    "row_fingerprint": row_fingerprint,
                }
            )
        except MDExportError as exc:
            gates.setdefault(exc.gate, _gate(exc.status, exc.reason))
            row_result.update({"status": exc.status, "reason": exc.reason})
        except Exception as exc:
            gates.setdefault("internal", _gate("failed", f"{exc.__class__.__name__}:{exc}"))
            row_result.update({"status": "failed", "reason": f"{exc.__class__.__name__}:{exc}"})
        results.append(row_result)

    counts = {status: sum(1 for row in results if row["status"] == status) for status in sorted(ROW_STATUSES)}
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "method_version": METHOD_VERSION,
        "generated_at": _utc_now(),
        "status": "completed" if results and all(row["status"] == "completed" for row in results) else "failed",
        "request": request.to_dict(),
        "request_fingerprint": config_payload["request_fingerprint"],
        "chemistry_capability": capability,
        "protonation_policy": project_policy,
        "selection_file": {"path": str(selection_file), "sha256": sha256_file(selection_file) if selection_file.is_file() else ""},
        "counts": counts,
        "rows": results,
    }
    manifest_file = output_root / "md_inputs_manifest.json"
    temporary = manifest_file.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(temporary, manifest_file)
    run_tracking_file = output_root / "run_tracking" / "manifest.json"
    run_tracking_file.parent.mkdir(parents=True, exist_ok=True)
    run_tracking_file.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    manifest["manifest_file"] = str(manifest_file)
    manifest["run_tracking_manifest_file"] = str(run_tracking_file)
    manifest["config_file"] = str(config_file)
    manifest["output_root"] = str(output_root)
    return manifest
