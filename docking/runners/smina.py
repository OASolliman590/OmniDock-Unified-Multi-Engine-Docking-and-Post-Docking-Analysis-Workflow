from __future__ import annotations

from pathlib import Path
from typing import List

from ..models import PairlistRow
from .vina import VinaRunner


class SminaRunner(VinaRunner):
    name = "smina"
    binary_name = "smina"

    def build_command(self, row: PairlistRow, pose_file: Path, log_file: Path) -> List[str]:
        command = super().build_command(row, pose_file=pose_file, log_file=log_file)
        scoring = str(self.runtime.get("scoring") or "").strip()
        if scoring:
            command.extend(["--scoring", scoring])
        return command
