from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
import math
from typing import Dict, List, Literal, Optional, Tuple

from .models import PairlistRow


ParameterMode = Literal["basic", "advanced"]
ParameterPreset = Literal["screening_fast", "balanced", "exhaustive"]

BASIC_PRESETS: Dict[str, Dict[str, int]] = {
    "screening_fast": {"exhaustiveness": 8, "num_modes": 10},
    "balanced": {"exhaustiveness": 16, "num_modes": 20},
    "exhaustive": {"exhaustiveness": 32, "num_modes": 40},
}


@dataclass
class CommonDockingParameters:
    exhaustiveness: int = 16
    num_modes: int = 20
    seed: Optional[int] = None
    box_scale: float = 1.0
    box_padding: float = 0.0

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class DockingParameterSchema:
    mode: str = "basic"
    preset: str = "balanced"
    common: CommonDockingParameters = field(default_factory=CommonDockingParameters)
    engine_extensions: Dict[str, Dict[str, object]] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, object]:
        payload = asdict(self)
        payload["common"] = self.common.to_dict()
        return payload


def resolve_parameter_schema(
    *,
    mode: str,
    preset: str,
    exhaustiveness: int,
    num_modes: int,
    seed: Optional[int],
    box_scale: float,
    box_padding: float,
    runtime_by_engine: Dict[str, Dict[str, object]],
) -> DockingParameterSchema:
    normalized_mode = str(mode or "basic").strip().lower()
    normalized_mode = normalized_mode if normalized_mode in {"basic", "advanced"} else "basic"
    normalized_preset = str(preset or "balanced").strip().lower()
    normalized_preset = normalized_preset if normalized_preset in BASIC_PRESETS else "balanced"

    common = CommonDockingParameters(
        exhaustiveness=int(exhaustiveness),
        num_modes=int(num_modes),
        seed=int(seed) if seed is not None else None,
        box_scale=float(box_scale),
        box_padding=float(box_padding),
    )
    if normalized_mode == "basic":
        preset_values = BASIC_PRESETS[normalized_preset]
        common.exhaustiveness = int(preset_values["exhaustiveness"])
        common.num_modes = int(preset_values["num_modes"])

    extensions: Dict[str, Dict[str, object]] = {}
    for engine, payload in runtime_by_engine.items():
        runtime = dict(payload or {})
        extensions[engine] = {
            "cpu": runtime.get("cpu"),
            "binary": runtime.get("binary"),
            "use_gpu": runtime.get("use_gpu"),
            "device": runtime.get("device"),
            "cnn_scoring": runtime.get("cnn_scoring"),
            "scoring": runtime.get("scoring"),
        }
    return DockingParameterSchema(
        mode=normalized_mode,
        preset=normalized_preset,
        common=common,
        engine_extensions=extensions,
    )


def validate_parameter_schema(schema: DockingParameterSchema, selected_engines: List[str]) -> Tuple[List[str], List[str]]:
    errors: List[str] = []
    warnings: List[str] = []

    common = schema.common
    if common.exhaustiveness < 1:
        errors.append("`exhaustiveness` must be >= 1.")
    if common.num_modes < 1:
        errors.append("`num_modes` must be >= 1.")
    if common.seed is not None and common.seed < 0:
        errors.append("`seed` must be >= 0.")
    if not math.isfinite(common.box_scale) or common.box_scale <= 0:
        errors.append("`box_scale` must be > 0.")
    if not math.isfinite(common.box_padding) or common.box_padding < 0:
        errors.append("`box_padding` must be >= 0.")

    if schema.mode == "basic" and schema.preset not in BASIC_PRESETS:
        errors.append(f"Unknown basic preset: {schema.preset}")
    if schema.mode == "advanced" and schema.preset in BASIC_PRESETS:
        warnings.append("Advanced mode selected; basic preset is kept only for traceability.")

    if "gnina" in selected_engines:
        gnina = schema.engine_extensions.get("gnina", {})
        cnn_scoring = str(gnina.get("cnn_scoring") or "").strip()
        if not cnn_scoring:
            warnings.append("GNINA `cnn_scoring` is empty; GNINA defaults will be used.")

    return errors, warnings


def apply_schema_to_runtime(
    schema: DockingParameterSchema,
    runtime_by_engine: Dict[str, Dict[str, object]],
) -> Dict[str, Dict[str, object]]:
    updated: Dict[str, Dict[str, object]] = {}
    common = schema.common
    for engine, payload in runtime_by_engine.items():
        runtime = dict(payload or {})
        runtime["exhaustiveness"] = int(common.exhaustiveness)
        runtime["num_modes"] = int(common.num_modes)
        if common.seed is not None:
            runtime["seed"] = int(common.seed)
        runtime["box_scale"] = float(common.box_scale)
        runtime["box_padding"] = float(common.box_padding)
        runtime["parameter_mode"] = schema.mode
        runtime["parameter_preset"] = schema.preset
        updated[engine] = runtime
    return updated


def transform_pairlist_rows(rows: List[PairlistRow], schema: DockingParameterSchema) -> List[PairlistRow]:
    scale = float(schema.common.box_scale)
    padding = float(schema.common.box_padding)
    transformed: List[PairlistRow] = []
    for row in rows:
        transformed.append(
            replace(
                row,
                center_x=float(row.center_x),
                center_y=float(row.center_y),
                center_z=float(row.center_z),
                size_x=float(row.size_x) * scale + (2.0 * padding),
                size_y=float(row.size_y) * scale + (2.0 * padding),
                size_z=float(row.size_z) * scale + (2.0 * padding),
            )
        )
    return transformed


def validate_pairlist_geometry(rows: List[PairlistRow]) -> List[str]:
    errors = []
    seen = set()
    for row in rows:
        values = [row.center_x, row.center_y, row.center_z, row.size_x, row.size_y, row.size_z]
        if not all(math.isfinite(float(value)) for value in values) or any(float(value) <= 0 for value in values[3:]):
            errors.append(f"{row.tag}: box centers must be finite and sizes must be finite and > 0.")
        if row.tag in seen:
            errors.append(f"Duplicate docking pair tag: {row.tag}")
        seen.add(row.tag)
    return errors
