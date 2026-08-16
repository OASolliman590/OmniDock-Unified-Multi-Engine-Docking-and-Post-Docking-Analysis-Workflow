"""HPC profile discovery and redaction.

Profiles carry site-specific secrets -- SSH targets, cluster accounts, home
directory paths. The security model in HPC_DEPLOYMENT_GUIDE.md keeps public
templates in ``examples/hpc_profiles/`` and real settings in a project's
``.workflow/hpc_profiles/``. The web layer must never render or log those
real values (FR-032, SC-009), so everything returned to the browser passes
through ``redact`` first.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from .config import REPO_ROOT

BUILTIN_PROFILE_DIR = REPO_ROOT / "examples" / "hpc_profiles"
PROJECT_PROFILE_SUBDIR = Path(".workflow") / "hpc_profiles"

REDACTED = "[redacted]"

# Keys whose values are always secret, matched case-insensitively anywhere
# in the nested structure.
_SECRET_KEYS = {
    "account", "slurm_account", "ssh_target", "user", "username",
    "password", "token", "key", "api_key", "secret",
    "project_root_base", "remote_project_dir", "home", "binary", "image",
}

# Values that look like a location or credential even under a benign key.
_SECRET_VALUE_PATTERNS = [
    re.compile(r"^[\w.\-]+@[\w.\-]+$"),        # user@host
    re.compile(r"^(/home/|/cluster/|/users/)"),  # absolute home-ish paths
    re.compile(r"^~[\\/]"),                     # ~/...
]


class HpcError(Exception):
    """Raised when a profile cannot be read."""


def _is_secret_key(key: str) -> bool:
    lowered = key.lower()
    return any(secret in lowered for secret in _SECRET_KEYS)


def _looks_secret(value: str) -> bool:
    return any(pattern.search(value) for pattern in _SECRET_VALUE_PATTERNS)


def redact(payload):
    """Recursively replace secret-shaped values with a placeholder."""
    if isinstance(payload, dict):
        return {
            key: (REDACTED if _is_secret_key(str(key)) else redact(value))
            for key, value in payload.items()
        }
    if isinstance(payload, list):
        return [redact(item) for item in payload]
    if isinstance(payload, str) and _looks_secret(payload):
        return REDACTED
    return payload


def _read_profile(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        raise HpcError(f"Could not read profile {path.name}: {exc}") from exc


def list_profiles(project_root: str | Path | None = None) -> list[dict]:
    """Profiles from the public templates and the project's own directory.

    Every returned profile is redacted -- callers cannot accidentally leak
    one by forgetting to filter.
    """
    profiles: list[dict] = []

    def collect(directory: Path, source: str) -> None:
        if not directory.is_dir():
            return
        for path in sorted(directory.glob("*.json")):
            try:
                raw = _read_profile(path)
            except HpcError:
                continue
            name = raw.get("name") or path.stem.replace(".template", "")
            profiles.append({
                "name": name,
                "source": source,
                "file": path.name,
                "description": raw.get("description", ""),
                "engines": sorted((raw.get("engines") or {}).keys()),
                "fields": redact(raw),
            })

    collect(BUILTIN_PROFILE_DIR, "builtin")
    if project_root:
        collect(Path(project_root) / PROJECT_PROFILE_SUBDIR, "project")

    return profiles


def profile_remote_target(project_root: str | Path | None, name: str) -> str | None:
    """The unredacted ssh target for a profile, for confirmation matching.

    Returned to server-side code only -- never serialized to the browser.
    """
    directories = [BUILTIN_PROFILE_DIR]
    if project_root:
        directories.insert(0, Path(project_root) / PROJECT_PROFILE_SUBDIR)

    for directory in directories:
        if not directory.is_dir():
            continue
        for path in directory.glob("*.json"):
            try:
                raw = _read_profile(path)
            except HpcError:
                continue
            candidate = raw.get("name") or path.stem.replace(".template", "")
            if candidate == name:
                target = (raw.get("remote") or {}).get("ssh_target")
                return str(target) if target else None
    return None
