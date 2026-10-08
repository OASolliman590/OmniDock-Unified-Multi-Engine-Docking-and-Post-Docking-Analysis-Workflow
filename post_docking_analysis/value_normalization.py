"""Shared scalar normalization helpers for analysis metadata contracts."""

from __future__ import annotations

import math
import numbers
from typing import Any

import numpy as np
import pandas as pd


_TRUE_TOKENS = {"1", "true", "t", "yes", "y", "on"}
_FALSE_TOKENS = {"", "0", "false", "f", "no", "n", "off", "none", "null", "nan", "na", "n/a"}


def normalize_boolean(value: Any, *, default: bool = False) -> bool:
    """Normalize serialized boolean values without treating non-empty text as true."""
    if value is None:
        return bool(default)
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, str):
        token = value.strip().lower()
        if token in _TRUE_TOKENS:
            return True
        if token in _FALSE_TOKENS:
            return False
        return bool(default)
    if isinstance(value, numbers.Number):
        try:
            if math.isnan(float(value)):
                return bool(default)
        except (TypeError, ValueError):
            return bool(default)
        return bool(value)
    try:
        if pd.isna(value):
            return bool(default)
    except (TypeError, ValueError):
        return bool(default)
    return bool(default)


def normalize_boolean_series(series: pd.Series, *, default: bool = False) -> pd.Series:
    """Return a boolean Series using :func:`normalize_boolean` elementwise."""
    return series.map(lambda value: normalize_boolean(value, default=default)).astype(bool)
