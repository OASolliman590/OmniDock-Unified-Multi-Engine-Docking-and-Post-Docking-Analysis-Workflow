"""
Engine HPC Adapter — auto-detect local/HPC directory layouts per engine.

This module extends the original GNINA-only adapter concept so Vina, Smina,
and AutoDock4 resolve the same core path contract:

- pose folder
- log folder
- receptors folder
- pairlist file

Callers can use ``detect_engine_layout(project_dir, engine)`` directly or the
engine-specific convenience helpers.
"""

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union


SUPPORTED_ENGINES = ("gnina", "vina", "smina", "autodock4")


def detect_engine_layout(project_dir: Union[str, Path], engine: str) -> Dict[str, object]:
    """
    Detect docking output layout under *project_dir* for a specific engine.

    Returns a dict with resolved paths and layout metadata. For compatibility
    with existing GNINA callers, ``sdf_folder`` is kept as an alias to the
    resolved pose folder.
    """
    token = str(engine or "").strip().lower()
    if token not in SUPPORTED_ENGINES:
        raise ValueError(f"Unsupported engine for layout detection: {engine}")

    root = Path(project_dir).resolve()
    dock_root = root / "4-Docking" if (root / "4-Docking").is_dir() else root

    pose_resolution = _resolve_artifact_directory(
        _pose_dir_candidates(root, dock_root, token),
        _pose_suffixes(token),
        artifact="poses",
    )
    log_resolution = _resolve_artifact_directory(
        _log_dir_candidates(root, dock_root, token),
        [".log", "_log"],
        artifact="logs",
    )
    score_resolution = _resolve_artifact_directory(
        _score_dir_candidates(root, dock_root, token),
        [".csv"],
        artifact="scores",
    )
    pose_dir = _optional_path(pose_resolution.get("selected"))
    logs_dir = _optional_path(log_resolution.get("selected"))
    receptors = _first_existing(root, ["receptors", "receptor"]) or _first_existing(dock_root, ["receptors", "receptor"])
    pairlist = _find_pairlist_file(dock_root)

    if logs_dir is not None:
        layout = "hpc"
        log_folder = logs_dir
    elif pose_dir is not None and _matching_files(pose_dir, [".log", "_log"]):
        layout = "local"
        log_folder = pose_dir
        log_resolution = _resolution_for_selected_directory(
            pose_dir,
            artifact="logs",
            suffixes=[".log", "_log"],
            reason="logs_colocated_with_poses",
        )
    else:
        layout = "local"
        log_folder = None

    warnings: List[str] = []
    for resolution in (pose_resolution, log_resolution, score_resolution):
        warning = str(resolution.get("warning") or "")
        if warning:
            warnings.append(warning)

    return {
        "pose_folder": pose_dir,
        "sdf_folder": pose_dir,  # compatibility alias with GNINA adapter contract
        "log_folder": log_folder,
        "score_folder": _optional_path(score_resolution.get("selected")),
        "receptors_folder": receptors,
        "pairlist_file": pairlist,
        "layout": layout,
        "engine": token,
        "artifact_resolution": {
            "poses": pose_resolution,
            "logs": log_resolution,
            "scores": score_resolution,
        },
        "ambiguous": bool(warnings),
        "warnings": warnings,
    }


def detect_gnina_layout(project_dir: Union[str, Path]) -> Dict[str, object]:
    return detect_engine_layout(project_dir, "gnina")


def detect_vina_layout(project_dir: Union[str, Path]) -> Dict[str, object]:
    return detect_engine_layout(project_dir, "vina")


def detect_smina_layout(project_dir: Union[str, Path]) -> Dict[str, object]:
    return detect_engine_layout(project_dir, "smina")


def detect_autodock4_layout(project_dir: Union[str, Path]) -> Dict[str, object]:
    return detect_engine_layout(project_dir, "autodock4")


def layout_to_cli_args(layout: Dict[str, object]) -> Dict[str, str]:
    pose = layout.get("pose_folder") or layout.get("sdf_folder")
    return {
        "pose_folder": str(pose) if pose else "",
        "sdf_folder": str(layout.get("sdf_folder")) if layout.get("sdf_folder") else "",
        "log_folder": str(layout.get("log_folder")) if layout.get("log_folder") else "",
        "receptors_folder": str(layout.get("receptors_folder")) if layout.get("receptors_folder") else "",
        "pairlist_file": str(layout.get("pairlist_file")) if layout.get("pairlist_file") else "",
    }


def _pose_suffixes(engine: str) -> List[str]:
    if engine == "gnina":
        return [".sdf"]
    if engine in {"vina", "smina"}:
        return [".pdbqt"]
    if engine == "autodock4":
        return [".dlg"]
    return []


def _pose_dir_candidates(root: Path, dock_root: Path, engine: str) -> List[Path]:
    if engine == "gnina":
        return [
            dock_root / "gnina_out" / "poses",
            dock_root / "gnina_out",
            dock_root / "docking_output",
            dock_root / "output",
            root / "3-Docking" / "gnina" / "poses",
            root / "engines" / "gnina" / "poses",
            root / "3-Docking" / "engines" / "gnina" / "poses",
        ]

    legacy_root = dock_root / f"{engine}_out"
    return [
        legacy_root / "poses",
        legacy_root,
        root / "3-Docking" / engine / "poses",
        root / "engines" / engine / "poses",
        root / "3-Docking" / "engines" / engine / "poses",
    ]


def _log_dir_candidates(root: Path, dock_root: Path, engine: str) -> List[Path]:
    if engine == "gnina":
        return [
            dock_root / "logs",
            dock_root / "log",
            root / "3-Docking" / "gnina" / "logs",
            root / "engines" / "gnina" / "logs",
        ]

    legacy_root = dock_root / f"{engine}_out"
    return [
        legacy_root / "logs",
        root / "3-Docking" / engine / "logs",
        root / "engines" / engine / "logs",
        root / "3-Docking" / "engines" / engine / "logs",
        dock_root / "logs",
    ]


def _score_dir_candidates(root: Path, dock_root: Path, engine: str) -> List[Path]:
    legacy_root = dock_root if engine == "gnina" else dock_root / f"{engine}_out"
    legacy_scores = legacy_root / ("results" if engine == "gnina" else "scores")
    return [
        legacy_scores,
        root / "3-Docking" / engine / "scores",
        root / "engines" / engine / "scores",
        root / "3-Docking" / "engines" / engine / "scores",
    ]


def _first_existing(root: Path, candidates: List[str]) -> Optional[Path]:
    for name in candidates:
        p = root / name
        if p.is_dir():
            return p
    return None


def _first_existing_dir(candidates: List[Path]) -> Optional[Path]:
    for path in candidates:
        if path.is_dir():
            return path
    return None


def _deduplicate_paths(candidates: Sequence[Path]) -> List[Path]:
    result: List[Path] = []
    seen = set()
    for candidate in candidates:
        path = Path(candidate).expanduser().resolve()
        key = str(path).casefold()
        if key in seen:
            continue
        seen.add(key)
        result.append(path)
    return result


def _matching_files(directory: Path, suffixes: Sequence[str]) -> List[Path]:
    if not directory.is_dir():
        return []
    matches: List[Path] = []
    for suffix in suffixes:
        pattern = f"*{suffix}"
        matches.extend(path for path in directory.glob(pattern) if path.is_file())
    return sorted(set(matches), key=lambda path: str(path))


def _resolve_artifact_directory(
    candidates: Sequence[Path],
    suffixes: Sequence[str],
    *,
    artifact: str,
) -> Dict[str, object]:
    """Resolve one input artifact without creating or modifying directories."""
    ordered = _deduplicate_paths(candidates)
    existing = [path for path in ordered if path.is_dir()]
    populated = [path for path in existing if _matching_files(path, suffixes)]
    selected = populated[0] if populated else None
    warning = ""
    status = "resolved" if selected is not None else "absent"
    if len(populated) > 1:
        status = "ambiguous"
        warning = (
            f"Multiple populated {artifact} directories detected; selected {selected} by deterministic "
            f"candidate precedence and retained diagnostics for: {', '.join(str(path) for path in populated)}"
        )
    return {
        "artifact": artifact,
        "status": status,
        "selected": str(selected) if selected is not None else "",
        "selection_reason": "first_populated_candidate" if selected is not None else "no_populated_candidate",
        "candidate_precedence": [str(path) for path in ordered],
        "existing_candidates": [str(path) for path in existing],
        "populated_candidates": [str(path) for path in populated],
        "warning": warning,
    }


def _resolution_for_selected_directory(
    selected: Path,
    *,
    artifact: str,
    suffixes: Sequence[str],
    reason: str,
) -> Dict[str, object]:
    matches = _matching_files(selected, suffixes)
    return {
        "artifact": artifact,
        "status": "resolved" if matches else "absent",
        "selected": str(selected) if matches else "",
        "selection_reason": reason if matches else "no_populated_candidate",
        "candidate_precedence": [str(selected)],
        "existing_candidates": [str(selected)] if selected.is_dir() else [],
        "populated_candidates": [str(selected)] if matches else [],
        "warning": "",
    }


def _optional_path(value: object) -> Optional[Path]:
    token = str(value or "").strip()
    return Path(token) if token else None


def _find_pairlist_file(root: Path) -> Optional[Path]:
    roots: List[Path] = []
    current = root
    for _ in range(5):
        roots.append(current)
        parent = current.parent
        if parent == current:
            break
        current = parent

    for candidate_root in roots:
        direct = candidate_root / "pairlist.csv"
        if direct.is_file():
            return direct

    for candidate_root in roots:
        for path in sorted(candidate_root.glob("*pairlist*.csv")):
            if path.is_file():
                return path

    return None


def _has_log_files(directory: Path) -> bool:
    return any(directory.glob("*.log"))


def _has_files_with_suffixes(directory: Path, suffixes: List[str]) -> bool:
    for suffix in suffixes:
        if any(directory.glob(f"*{suffix}")):
            return True
    return False
