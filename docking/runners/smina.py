from __future__ import annotations

from .vina import VinaRunner


class SminaRunner(VinaRunner):
    name = "smina"
    binary_name = "smina"
