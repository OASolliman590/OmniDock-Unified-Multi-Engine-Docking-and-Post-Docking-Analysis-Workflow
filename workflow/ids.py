from __future__ import annotations

import re
from datetime import datetime, timezone
from typing import Optional


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def iso_timestamp(value: Optional[datetime] = None) -> str:
    current = value or utc_now()
    if current.tzinfo is None:
        current = current.replace(tzinfo=timezone.utc)
    return current.astimezone(timezone.utc).isoformat()


def format_run_timestamp(value: Optional[datetime] = None) -> str:
    current = value or utc_now()
    if current.tzinfo is None:
        current = current.replace(tzinfo=timezone.utc)
    return current.astimezone(timezone.utc).strftime("%Y%m%d_%H%M%S")


def build_run_id(prefix: str, *, value: Optional[datetime] = None, suffix: str = "") -> str:
    cleaned_prefix = re.sub(r"[^a-z0-9]+", "_", str(prefix or "run").strip().lower()).strip("_") or "run"
    token = format_run_timestamp(value)
    cleaned_suffix = re.sub(r"[^a-zA-Z0-9]+", "_", str(suffix or "").strip()).strip("_")
    if cleaned_suffix:
        return f"{cleaned_prefix}_{token}_{cleaned_suffix}"
    return f"{cleaned_prefix}_{token}"

