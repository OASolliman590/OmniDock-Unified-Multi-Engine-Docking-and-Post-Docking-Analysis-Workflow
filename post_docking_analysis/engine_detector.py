"""
Engine detection and auto-routing for DockForge post-docking analysis.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

from docking.project_layout import (
    ensure_numbered_output_layout,
    load_manifest,
    pairlist_path,
    save_manifest,
)
from post_docking_analysis.docking_parser import parse_autodock4_dlg, parse_vina_pdbqt
from post_docking_analysis.engine_hpc_adapter import detect_engine_layout
from post_docking_analysis.generate_scores_csv import parse_gnina_log


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
    for engine in SUPPORTED_ENGINES:
        detected = detect_engine_layout(project_dir, engine)
        if detected.get("pose_folder") or detected.get("log_folder"):
            candidates.append(engine)
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
    detected = detect_engine_layout(project_dir, "gnina")
    pose_dir = Path(detected["pose_folder"]) if detected.get("pose_folder") else None
    log_dir = Path(detected["log_folder"]) if detected.get("log_folder") else None
    pose_files = sorted(pose_dir.glob("*.sdf")) if pose_dir else []
    log_files = (
        sorted(
            set(log_dir.glob("*.log")).union(log_dir.glob("*_log")),
            key=lambda path: str(path),
        )
        if log_dir
        else []
    )
    valid_score_rows = 0
    matched_marker = ""
    score_columns: List[str] = []
    for log_file in log_files:
        content = _safe_read_text(log_file)
        lower = content.lower()
        if not matched_marker:
            for marker in GNINA_LOG_MARKERS:
                if marker in lower:
                    matched_marker = marker
                    break
        rows, _ = parse_gnina_log(log_file, {})
        if rows:
            valid_score_rows += len(rows)
            score_columns = ["vina_affinity", "cnn_score", "cnn_affinity"]
            break
    return {
        "pose_files_found": len(pose_files),
        "pose_format": "sdf",
        "valid_score_rows": int(valid_score_rows),
        "score_columns": score_columns,
        "log_marker_matched": matched_marker,
        "rmsd_available": False,
        "input_layout": _layout_diagnostics(detected),
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
    detected = detect_engine_layout(project_dir, engine)
    pose_dir = Path(detected["pose_folder"]) if detected.get("pose_folder") else None
    log_dir = Path(detected["log_folder"]) if detected.get("log_folder") else None
    score_dir = Path(detected["score_folder"]) if detected.get("score_folder") else None
    pose_files = sorted(pose_dir.glob("*.pdbqt")) if pose_dir else []
    score_file = score_dir / "normalized_scores.csv" if score_dir else None
    log_files = (
        sorted(
            set(log_dir.glob("*.log")).union(log_dir.glob("*_log")),
            key=lambda path: str(path),
        )
        if log_dir
        else []
    )
    valid_score_rows = 0
    score_columns = ["vina_affinity"]
    rmsd_available = engine == "vina"
    if score_file is not None and score_file.exists():
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
        "input_layout": _layout_diagnostics(detected),
    }
    if engine == "smina":
        weights, inconsistent = _extract_smina_weights(log_files)
        payload["smina_scoring_weights"] = weights
        payload["inconsistent_scoring_weights"] = bool(inconsistent)
    return payload


def _autodock4_details(project_dir: Path) -> Dict[str, object]:
    detected = detect_engine_layout(project_dir, "autodock4")
    pose_dir = Path(detected["pose_folder"]) if detected.get("pose_folder") else None
    score_dir = Path(detected["score_folder"]) if detected.get("score_folder") else None
    pose_files = sorted(pose_dir.glob("*.dlg")) if pose_dir else []
    score_file = score_dir / "normalized_scores.csv" if score_dir else None
    valid_score_rows = 0
    score_columns = ["autodock4_affinity"]
    if score_file is not None and score_file.exists():
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
        "input_layout": _layout_diagnostics(detected),
    }


def _layout_diagnostics(detected: Dict[str, object]) -> Dict[str, object]:
    return {
        "pose_folder": str(detected.get("pose_folder") or ""),
        "log_folder": str(detected.get("log_folder") or ""),
        "score_folder": str(detected.get("score_folder") or ""),
        "layout": str(detected.get("layout") or ""),
        "ambiguous": bool(detected.get("ambiguous", False)),
        "warnings": list(detected.get("warnings") or []),
        "artifact_resolution": dict(detected.get("artifact_resolution") or {}),
    }


def _detection_source_signature(project_dir: Path, manifest: Dict[str, object]) -> str:
    rows = [
        "manifest_engines="
        + ",".join(
            sorted(
                str(engine).strip().lower()
                for engine in (manifest.get("engines", []) or [])
                if str(engine).strip()
            )
        )
    ]
    suffixes = {
        "poses": {"gnina": (".sdf",), "vina": (".pdbqt",), "smina": (".pdbqt",), "autodock4": (".dlg",)},
        "logs": {engine: (".log", "_log") for engine in SUPPORTED_ENGINES},
        "scores": {engine: (".csv",) for engine in SUPPORTED_ENGINES},
    }
    files = set()
    for engine in SUPPORTED_ENGINES:
        detected = detect_engine_layout(project_dir, engine)
        resolution = dict(detected.get("artifact_resolution") or {})
        for artifact in ("poses", "logs", "scores"):
            artifact_row = dict(resolution.get(artifact) or {})
            for raw_dir in artifact_row.get("populated_candidates", []) or []:
                directory = Path(str(raw_dir))
                if not directory.is_dir():
                    continue
                for suffix in suffixes[artifact][engine]:
                    files.update(path.resolve() for path in directory.glob(f"*{suffix}") if path.is_file())
    try:
        pairlist = pairlist_path(project_dir)
        if pairlist.is_file():
            files.add(pairlist.resolve())
    except Exception:
        pass
    for path in sorted(files, key=lambda item: str(item)):
        try:
            stat = path.stat()
            rows.append(f"{path}|{stat.st_mtime_ns}|{stat.st_size}")
        except OSError:
            rows.append(f"{path}|stat_error")
    return hashlib.sha256("\n".join(rows).encode("utf-8")).hexdigest()


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

    source_signature = _detection_source_signature(root, manifest)
    if not redetect and isinstance(manifest.get("detected_engines"), dict):
        cached = dict(manifest["detected_engines"])
        if str(cached.get("source_signature") or "") == source_signature:
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
        "detection_method": detection_method,
        "detection_override": detection_override,
        "override_reason": f"--engine {override}" if detection_override else None,
        "routing_decision": routing_decision,
        "routed_to_engine": routed_to_engine,
        "valid_engines": valid_engines,
        "engines": engines_payload,
        "min_coverage_pct": float(min_coverage_pct),
        "total_pairlist_pairs": int(total_pairs),
        "source_signature": source_signature,
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
