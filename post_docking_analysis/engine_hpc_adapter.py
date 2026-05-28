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
from typing import Dict, List, Optional, Union


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

    pose_dir = _first_existing_dir(_pose_dir_candidates(root, dock_root, token))
    logs_dir = _first_existing_dir(_log_dir_candidates(root, dock_root, token))
    receptors = _first_existing(root, ["receptors", "receptor"]) or _first_existing(dock_root, ["receptors", "receptor"])
    pairlist = _find_pairlist_file(dock_root)

    pose_suffixes = _pose_suffixes(token)
    if logs_dir is not None and _has_log_files(logs_dir):
        layout = "hpc"
        log_folder = logs_dir
    elif pose_dir is not None and _has_log_files(pose_dir):
        layout = "local"
        log_folder = pose_dir
    else:
        layout = "hpc" if logs_dir is not None else "local"
        log_folder = logs_dir or pose_dir

    # If the detected directory is an engine root, prefer nested poses/ when valid.
    if pose_dir is not None and pose_dir.is_dir():
        nested_poses = pose_dir / "poses"
        if nested_poses.is_dir() and _has_files_with_suffixes(nested_poses, pose_suffixes):
            pose_dir = nested_poses

    if pose_dir is not None and not _has_files_with_suffixes(pose_dir, pose_suffixes):
        pose_dir = None

    return {
        "pose_folder": pose_dir,
        "sdf_folder": pose_dir,  # compatibility alias with GNINA adapter contract
        "log_folder": log_folder,
        "receptors_folder": receptors,
        "pairlist_file": pairlist,
        "layout": layout,
        "engine": token,
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
            dock_root / "gnina_out",
            dock_root / "docking_output",
            dock_root / "output",
            root / "3-Docking" / "gnina" / "poses",
            root / "engines" / "gnina" / "poses",
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
