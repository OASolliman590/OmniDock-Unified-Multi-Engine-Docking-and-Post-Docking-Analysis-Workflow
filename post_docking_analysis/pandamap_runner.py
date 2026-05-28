"""
Shared PandaMap CLI runner utilities.

This module centralizes PandaMap command execution and CLI mode detection so
all analyzers behave consistently.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional
import logging
import os
import re
import subprocess


class PandaMapRunner:
    """Execute PandaMap commands and detect supported CLI style."""

    def __init__(self, conda_env: str = "pandamap", logger: Optional[logging.Logger] = None):
        self.conda_env = conda_env
        self.logger = logger or logging.getLogger(__name__)
        self.conda_executable = self._resolve_conda_executable()
        self._cli_mode: Optional[str] = None
        self._supported_options: Optional[set[str]] = None
        self._help_text: Optional[str] = None
        self._version_text: Optional[str] = None

    @staticmethod
    def _resolve_conda_executable() -> str:
        """
        Resolve a working conda binary path.

        `conda` is not always on PATH in non-interactive shells used by automation.
        """
        conda_env = os.environ.get("CONDA_EXE")
        candidates = []
        if conda_env:
            candidates.append(conda_env)
        home = Path.home()
        candidates.extend(
            [
                "conda",
                str(home / "opt" / "miniconda3" / "bin" / "conda"),
                str(home / "miniconda3" / "bin" / "conda"),
                str(home / "anaconda3" / "bin" / "conda"),
            ]
        )
        for candidate in candidates:
            if candidate == "conda":
                return candidate
            if Path(candidate).exists():
                return candidate
        return "conda"

    def run(
        self,
        args: List[str],
        timeout: int = 300,
        cwd: Optional[Path] = None
    ) -> subprocess.CompletedProcess:
        """Run PandaMap command inside the configured conda environment."""
        cmd = [self.conda_executable, "run", "-n", self.conda_env, "pandamap", *args]
        return subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=str(cwd) if cwd else None,
            check=False,
        )

    def get_help_text(self) -> str:
        """Return CLI help text for the installed PandaMap command."""
        if self._help_text is not None:
            return self._help_text
        try:
            result = self.run(["-h"], timeout=60)
            self._help_text = f"{result.stdout}\n{result.stderr}".strip()
        except Exception:
            self._help_text = ""
        return self._help_text

    def get_version(self) -> str:
        """Return installed PandaMap version string when available."""
        if self._version_text is not None:
            return self._version_text
        for args in (["--version"], ["-v"]):
            try:
                result = self.run(args, timeout=30)
                text = f"{result.stdout}\n{result.stderr}".strip()
                if text:
                    match = re.search(r"(\d+\.\d+\.\d+)", text)
                    self._version_text = match.group(1) if match else text.splitlines()[0].strip()
                    return self._version_text
            except Exception:
                continue
        self._version_text = ""
        return self._version_text

    def detect_cli_mode(self) -> str:
        """
        Detect PandaMap CLI style.

        Returns
        -------
        str
            "single" for positional structure-file mode, otherwise
            "subcommand" for generate/visualize style.
        """
        if self._cli_mode:
            return self._cli_mode

        try:
            help_text = self.get_help_text().lower()
            if "structure_file" in help_text and "--3d-output" in help_text:
                self._cli_mode = "single"
            else:
                self._cli_mode = "subcommand"
        except Exception:
            self._cli_mode = "subcommand"

        self.logger.info(f"🐼 Detected PandaMap CLI mode: {self._cli_mode}")
        return self._cli_mode

    def get_supported_options(self) -> set[str]:
        """
        Return the set of long-form CLI options supported by installed pandamap.
        """
        if self._supported_options is not None:
            return self._supported_options

        options: set[str] = set()
        try:
            help_text = self.get_help_text()
            for token in re.findall(r"--[A-Za-z0-9][A-Za-z0-9-]*", help_text):
                options.add(token.strip())
        except Exception:
            options = set()

        self._supported_options = options
        return self._supported_options

    def describe_capabilities(self) -> Dict[str, object]:
        """
        Return a normalized snapshot of PandaMap runtime capabilities.
        """
        options = self.get_supported_options()
        option_matrix = {
            "ligand_selection": "--ligand" in options,
            "custom_2d_output": "--output" in options,
            "text_report": "--report" in options,
            "custom_report_path": "--report-file" in options,
            "interactive_3d": "--3d" in options,
            "custom_3d_output": "--3d-output" in options or "--output" in options,
            "delta_g_estimation": "--deltaG" in options,
            "custom_title": "--title" in options,
            "dpi_control": "--dpi" in options,
            "viewer_width": "--width" in options,
            "viewer_height": "--height" in options,
            "hide_surface": "--no-surface" in options,
            "disable_3d_cues": "--no-3d-cues" in options,
        }
        return {
            "conda_env": self.conda_env,
            "conda_executable": self.conda_executable,
            "version": self.get_version(),
            "cli_mode": self.detect_cli_mode(),
            "supported_options": sorted(options),
            "option_matrix": option_matrix,
            "documented_input_formats": ["pdb", "mmcif/cif", "pdbqt"],
        }
