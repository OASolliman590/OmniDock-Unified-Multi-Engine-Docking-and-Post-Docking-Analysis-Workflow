"""
Engine detection and auto-routing for DockForge post-docking analysis.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

from docking.project_layout import (
    ensure_engine_layout,
    ensure_numbered_output_layout,
    load_manifest,
    pairlist_path,
    save_manifest,
)
from post_docking_analysis.docking_parser import parse_autodock4_dlg, parse_vina_pdbqt
from post_docking_analysis.generate_scores_csv import parse_gnina_log, load_pairlist_mapping
from docking.runners.job_contract import parse_sdf_scores
from post_docking_analysis.pose_geometry import content_hash
from post_docking_analysis.score_import import read_explicit_score_import


SUPPORTED_ENGINES = ("gnina", "smina", "vina", "autodock4")
GNINA_LOG_MARKERS = ("cnn_affinity", "cnnscore", "cnn score", "no gpu detected")
SMINA_WEIGHT_PATTERN = re.compile(
    r"^\s*(gauss|repulsion|hydrophobic|non_dir_h_bond|num_tors_div)\s*[:=]?\s*([+-]?\d+(?:\.\d+)?)",
    re.IGNORECASE,
)


def _safe_read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return ""


def _candidate_engines_from_manifest(project_dir: Path) -> List[str]:
    try:
        manifest = load_manifest(project_dir)
    except Exception:
        return []
    raw = manifest.get("engines", []) or []
    seen: List[str] = []
    for engine in raw:
        token = str(engine or "").strip().lower()
        if token in SUPPORTED_ENGINES and token not in seen:
            seen.append(token)
    return seen


def _candidate_engines_from_filesystem(project_dir: Path) -> List[str]:
    candidates: List[str] = []
    dock_root = project_dir / "4-Docking" if (project_dir / "4-Docking").exists() else project_dir
    for engine in SUPPORTED_ENGINES:
        if (dock_root / f"{engine}_out").exists():
            candidates.append(engine)
    logs_dir = dock_root / "logs"
    if logs_dir.exists():
        for log_file in sorted(logs_dir.glob("*.log")):
            content = _safe_read_text(log_file).lower()
            if any(marker in content for marker in GNINA_LOG_MARKERS):
                if "gnina" not in candidates:
                    candidates.append("gnina")
    return candidates


def _pairlist_size(project_dir: Path) -> int:
    try:
        pairlist = pairlist_path(project_dir)
        if pairlist.exists():
            return int(len(pd.read_csv(pairlist)))
    except Exception:
        return 0
    return 0


def _gnina_details(project_dir: Path) -> Dict[str, object]:
    layout = ensure_engine_layout(project_dir, "gnina")
    mapping = load_pairlist_mapping(pairlist_path(project_dir)) if pairlist_path(project_dir).is_file() else None
    pose_files = [path for path in sorted(layout["poses"].glob("*.sdf")) if mapping is None or path.stem in mapping]
    valid_score_rows = 0
    matched_marker = ""
    score_columns = []
    for pose_file in pose_files:
        try:
            rows = parse_sdf_scores(pose_file)
        except (ValueError, KeyError, OSError):
            rows = []
        if rows:
            valid_score_rows += len(rows)
            score_columns = ["vina_affinity", "cnn_score", "cnn_affinity"]
            matched_marker = "SDF properties"
    if not valid_score_rows:
        for log_file in sorted(layout["logs"].glob("*.log")):
            if not (layout["poses"] / f"{log_file.stem}.sdf").is_file():
                continue
            rows, _ = parse_gnina_log(log_file, mapping)
            valid_score_rows += len(rows)
            if rows:
                score_columns = ["vina_affinity", "cnn_score", "cnn_affinity"]
                matched_marker = "GNINA score table"
    return {
        "pose_files_found": len(pose_files),
        "pose_format": "sdf",
        "valid_score_rows": int(valid_score_rows),
        "score_columns": score_columns,
        "log_marker_matched": matched_marker,
        "rmsd_available": False,
    }


def _extract_smina_weights(log_files: Iterable[Path]) -> Tuple[Dict[str, float], bool]:
    merged: Dict[str, float] = {}
    inconsistent = False
    for log_file in log_files:
        weights: Dict[str, float] = {}
        for line in _safe_read_text(log_file).splitlines():
            match = SMINA_WEIGHT_PATTERN.match(line)
            if not match:
                continue
            try:
                weights[match.group(1).lower()] = float(match.group(2))
            except ValueError:
                continue
        if not weights:
            continue
        if not merged:
            merged = dict(weights)
            continue
        if weights != merged:
            inconsistent = True
    return merged, inconsistent


def _vina_family_details(project_dir: Path, engine: str) -> Dict[str, object]:
    layout = ensure_engine_layout(project_dir, engine)
    pose_files = sorted(layout["poses"].glob("*.pdbqt"))
    score_file = layout["scores"] / "normalized_scores.csv"
    log_files = sorted(layout["logs"].glob("*.log"))
    valid_score_rows = 0
    score_columns = ["vina_affinity"]
    rmsd_available = engine == "vina"
    if score_file.exists():
        try:
            frame = pd.read_csv(score_file)
            valid_score_rows = int(len(frame))
            if engine == "smina" and {"rmsd_lb", "rmsd_ub"}.issubset(frame.columns):
                rmsd_series = pd.to_numeric(frame["rmsd_lb"], errors="coerce").dropna()
                if not rmsd_series.empty and (rmsd_series == -1.0).all():
                    rmsd_available = False
        except Exception:
            valid_score_rows = 0
    if valid_score_rows == 0:
        for pose_file in pose_files[:3]:
            parsed = parse_vina_pdbqt(pose_file)
            if not parsed.empty:
                valid_score_rows += int(len(parsed))
                if engine == "smina":
                    rmsd_series = pd.to_numeric(parsed.get("rmsd_lb"), errors="coerce").dropna()
                    if not rmsd_series.empty and (rmsd_series == -1.0).all():
                        rmsd_available = False
                break
    payload: Dict[str, object] = {
        "pose_files_found": len(pose_files),
        "pose_format": "pdbqt",
        "valid_score_rows": int(valid_score_rows),
        "score_columns": score_columns,
        "rmsd_available": bool(rmsd_available),
    }
    if engine == "smina":
        weights, inconsistent = _extract_smina_weights(log_files)
        payload["smina_scoring_weights"] = weights
        payload["inconsistent_scoring_weights"] = bool(inconsistent)
    return payload


def _autodock4_details(project_dir: Path) -> Dict[str, object]:
    layout = ensure_engine_layout(project_dir, "autodock4")
    pose_files = sorted(layout["poses"].glob("*.dlg"))
    score_file = layout["scores"] / "normalized_scores.csv"
    valid_score_rows = 0
    score_columns = ["autodock4_affinity"]
    if score_file.exists():
        try:
            frame = pd.read_csv(score_file)
            valid_score_rows = int(len(frame))
        except Exception:
            valid_score_rows = 0
    if valid_score_rows == 0:
        for pose_file in pose_files[:3]:
            parsed = parse_autodock4_dlg(pose_file)
            if not parsed.empty:
                valid_score_rows += int(len(parsed))
                break
    return {
        "pose_files_found": len(pose_files),
        "pose_format": "dlg",
        "valid_score_rows": int(valid_score_rows),
        "score_columns": score_columns,
        "rmsd_available": False,
    }


def _engine_status(details: Dict[str, object], *, total_pairs: int, min_coverage_pct: float) -> str:
    pose_files = int(details.get("pose_files_found", 0) or 0)
    valid_score_rows = int(details.get("valid_score_rows", 0) or 0)
    if pose_files <= 0 or valid_score_rows <= 0:
        return "absent"
    coverage = float(details.get("coverage_pct", 0.0) or 0.0)
    if total_pairs <= 0:
        return "valid"
    return "partial" if coverage < float(min_coverage_pct) else "valid"


def detect_engines(
    project_dir: str | Path,
    *,
    override_engine: Optional[str] = None,
    redetect: bool = False,
    min_coverage_pct: float = 30.0,
) -> Dict[str, object]:
    root = Path(project_dir).expanduser().resolve()
    manifest: Dict[str, object] = {}
    try:
        manifest = load_manifest(root)
    except Exception:
        manifest = {}

    import hashlib
    digest = hashlib.sha256()
    for engine in SUPPORTED_ENGINES:
        layout = ensure_engine_layout(root, engine)
        for directory in (layout["poses"], layout["logs"], layout["scores"]):
            for path in sorted(directory.glob("*")):
                if path.is_file():
                    digest.update(str(path).encode())
                    digest.update(content_hash(path).encode())
    if pairlist_path(root).is_file():
        digest.update(content_hash(pairlist_path(root)).encode())
    digest.update(content_hash(Path(__file__)).encode())
    source_fingerprint = digest.hexdigest()
    if not redetect and isinstance(manifest.get("detected_engines"), dict) and manifest["detected_engines"].get("source_fingerprint") == source_fingerprint:
        cached = dict(manifest["detected_engines"])
        cached["cache_hit"] = True
        override = str(override_engine or "").strip().lower()
        if override:
            cached["detection_override"] = True
            cached["override_reason"] = f"--engine {override}"
            cached["routed_to_engine"] = override
            cached["routing_decision"] = "single_engine"
        return cached

    candidates = _candidate_engines_from_manifest(root)
    detection_method = "manifest" if candidates else "filesystem"
    if not candidates:
        candidates = _candidate_engines_from_filesystem(root)

    total_pairs = _pairlist_size(root)
    engines_payload: Dict[str, Dict[str, object]] = {}
    valid_engines: List[str] = []

    for engine in candidates:
        if engine == "gnina":
            details = _gnina_details(root)
        elif engine == "autodock4":
            details = _autodock4_details(root)
        else:
            details = _vina_family_details(root, engine)
        layout = ensure_engine_layout(root, engine)
        if not details["pose_files_found"]:
            imported = read_explicit_score_import(layout["scores"])
            if imported is not None and not imported.empty:
                if not imported["engine"].eq(engine).all():
                    raise ValueError(f"Imported score engine must match folder engine {engine}")
                details["pose_files_found"] = int(imported["tag"].nunique())
                details["valid_score_rows"] = len(imported)
                details["pose_format"] = "explicit_score_import_geometry_unavailable"
        coverage_base = total_pairs if total_pairs > 0 else int(details.get("pose_files_found", 0) or 0)
        pose_files_found = int(details.get("pose_files_found", 0) or 0)
        coverage_pct = 0.0 if coverage_base <= 0 else round((pose_files_found / float(coverage_base)) * 100.0, 2)
        details["coverage_pct"] = coverage_pct
        details["status"] = _engine_status(details, total_pairs=total_pairs, min_coverage_pct=min_coverage_pct)
        engines_payload[engine] = details
        if details["status"] != "absent":
            valid_engines.append(engine)

    override = str(override_engine or "").strip().lower()
    detection_override = False
    routed_to_engine: Optional[str] = None
    if override:
        detection_override = True
        routed_to_engine = override
        routing_decision = "single_engine"
    elif len(valid_engines) == 1:
        routed_to_engine = valid_engines[0]
        routing_decision = "single_engine"
    elif len(valid_engines) >= 2:
        routing_decision = "multi_engine"
    else:
        routing_decision = "no_valid_engines"

    report = {
        "source_fingerprint": source_fingerprint,
        "detection_method": detection_method,
        "detection_override": detection_override,
        "override_reason": f"--engine {override}" if detection_override else None,
        "routing_decision": routing_decision,
        "routed_to_engine": routed_to_engine,
        "valid_engines": valid_engines,
        "engines": engines_payload,
        "min_coverage_pct": float(min_coverage_pct),
        "total_pairlist_pairs": int(total_pairs),
        "cache_hit": False,
    }

    numbered = ensure_numbered_output_layout(root)
    report_path = numbered["post_metadata"] / "engine_detection_report.json"
    report["report_file"] = str(report_path)
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    if manifest:
        manifest["detected_engines"] = report
        try:
            save_manifest(root, manifest)
        except Exception:
            pass
    return report
