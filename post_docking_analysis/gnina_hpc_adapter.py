"""
Backward-compatible GNINA adapter wrapper.

The canonical engine-aware adapter now lives in
``post_docking_analysis.engine_hpc_adapter``.
"""

from pathlib import Path
from typing import Dict, Optional, Union

from .engine_hpc_adapter import detect_gnina_layout as _detect_gnina_layout
from .engine_hpc_adapter import layout_to_cli_args as _layout_to_cli_args


def detect_gnina_layout(project_dir: Union[str, Path]) -> Dict[str, Optional[Path]]:
    return _detect_gnina_layout(project_dir)


def layout_to_cli_args(layout: Dict[str, Optional[Path]]) -> Dict[str, str]:
    raw = _layout_to_cli_args(layout)
    return {
        "sdf_folder": raw.get("sdf_folder", ""),
        "log_folder": raw.get("log_folder", ""),
        "receptors_folder": raw.get("receptors_folder", ""),
        "pairlist_file": raw.get("pairlist_file", ""),
    }
