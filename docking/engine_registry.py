from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable

from .runners.autodock4 import AutoDock4Runner
from .runners.gnina import GninaRunner
from .runners.smina import SminaRunner
from .runners.vina import VinaRunner


ENGINE_REGISTRY = {
    "gnina": GninaRunner,
    "vina": VinaRunner,
    "smina": SminaRunner,
    "autodock4": AutoDock4Runner,
}


def normalize_engines(raw_value: Iterable[str]) -> list[str]:
    engines = []
    for item in raw_value:
        for token in str(item).split(","):
            candidate = token.strip().lower()
            if candidate:
                engines.append(candidate)
    invalid = sorted(set(engine for engine in engines if engine not in ENGINE_REGISTRY))
    if invalid:
        raise ValueError(f"Unsupported engines: {', '.join(invalid)}")
    return list(dict.fromkeys(engines))


def build_runner(engine: str, project_root: Path, runtime: Dict[str, object]):
    return ENGINE_REGISTRY[engine](project_root, runtime=runtime)
