"""Smina HPC adapter wrapper for engine-specific layout detection."""

from pathlib import Path
from typing import Dict, Optional, Union

from .engine_hpc_adapter import detect_smina_layout as _detect_smina_layout
from .engine_hpc_adapter import layout_to_cli_args as _layout_to_cli_args


def detect_smina_layout(project_dir: Union[str, Path]) -> Dict[str, Optional[Path]]:
    return _detect_smina_layout(project_dir)


def layout_to_cli_args(layout: Dict[str, Optional[Path]]) -> Dict[str, str]:
    raw = _layout_to_cli_args(layout)
    return {
        "pose_folder": raw.get("pose_folder", ""),
        "log_folder": raw.get("log_folder", ""),
        "receptors_folder": raw.get("receptors_folder", ""),
        "pairlist_file": raw.get("pairlist_file", ""),
    }
