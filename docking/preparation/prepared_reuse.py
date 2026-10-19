"""Prepare once: the reuse decision shared by receptor and ligand preparation (Spec 037 R4a).

An existing prepared output is reused (``reused_unchanged``) only when its provenance record matches
all of the following at the time of the new request:

- the resolved protonation policy hash (``policy_hash``; ``None`` when no project policy is stored);
- the request fingerprint (the effective values that shape the output and that the policy hash does
  not cover, for example the pH and force field of a run without a project policy);
- the input file SHA-256 (``input_sha256``);
- the output file SHA-256 (``output_sha256``), so a tampered or truncated output is re-prepared.

``force`` (``--force`` or ``preparation.force``) always re-prepares. A record written before R4a has no
``request_fingerprint`` and is therefore re-prepared once.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

REUSED_UNCHANGED = "reused_unchanged"
PREPARED = "prepared"
REQUEST_FINGERPRINT_VERSION = "prepare_once_request_v1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def request_fingerprint(request: Dict[str, Any]) -> str:
    """Canonical hash of the effective request values (sorted keys, no defaults filled in here)."""
    payload = json.dumps(
        {"version": REQUEST_FINGERPRINT_VERSION, "request": request},
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def load_provenance(path: Path) -> Tuple[Optional[Dict[str, Any]], str]:
    """Provenance record at ``path``: ``(record, reason)``. A missing or unreadable record gives no record."""
    if not Path(path).is_file():
        return None, "no_provenance"
    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None, "provenance_unreadable"
    if not isinstance(payload, dict):
        return None, "provenance_unreadable"
    return payload, ""


def reuse_block_reason(
    provenance: Optional[Dict[str, Any]],
    *,
    output: Path,
    input_file: Path,
    policy_hash: Optional[str],
    fingerprint: str,
    force: bool,
) -> str:
    """Why the existing output must be prepared again. An empty string means it may be reused."""
    if force:
        return "forced"
    if provenance is None:
        return "no_provenance"
    if not Path(output).is_file():
        return "output_missing"
    if "policy_hash" not in provenance or provenance.get("policy_hash") != policy_hash:
        return "policy_changed"
    if provenance.get("request_fingerprint") != fingerprint:
        return "request_changed" if provenance.get("request_fingerprint") else "provenance_predates_prepare_once"
    try:
        current_input = sha256_file(input_file)
    except OSError:
        return "input_missing"
    if provenance.get("input_sha256") != current_input:
        return "input_changed"
    if not provenance.get("output_sha256") or provenance.get("output_sha256") != sha256_file(output):
        return "output_changed"
    return ""
