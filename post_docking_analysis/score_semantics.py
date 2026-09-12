"""Score meaning shared by selection, normalization, and reporting."""
from dataclasses import dataclass
import math

@dataclass(frozen=True)
class ScoreSpec:
    name: str
    unit: str
    lower_is_better: bool
    purpose: str

SCORES = {
    "affinity_kcal_mol": ScoreSpec("affinity_kcal_mol", "kcal/mol", True, "empirical score"),
    "vina_affinity": ScoreSpec("vina_affinity", "kcal/mol", True, "empirical score"),
    "autodock4_affinity": ScoreSpec("autodock4_affinity", "kcal/mol", True, "empirical score"),
    "cnn_affinity": ScoreSpec("cnn_affinity", "pK", False, "predicted affinity"),
    "cnn_score": ScoreSpec("cnn_score", "dimensionless", False, "pose ranking"),
    "consensus_score": ScoreSpec("consensus_score", "dimensionless", False, "relative prioritization"),
}

def score_spec(name: str) -> ScoreSpec:
    try:
        return SCORES[str(name)]
    except KeyError as exc:
        raise ValueError(f"Unknown score semantics: {name}") from exc

def sort_value(value, metric: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.inf
    if not math.isfinite(number):
        return math.inf
    return number if score_spec(metric).lower_is_better else -number

def explicit_true(value) -> bool:
    """Missing/unknown flags are never positive evidence."""
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    try:
        return bool(math.isfinite(float(value)) and float(value) == 1.0)
    except (TypeError, ValueError):
        return False
