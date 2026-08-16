"""Locations and defaults for the web UI's own data.

Everything the web layer owns lives outside project directories so that a
project tree is never polluted by UI bookkeeping, and so a project can be
registered even when it is read-only (FR-035).
"""

from __future__ import annotations

import os
from pathlib import Path

# Repo root: the directory containing main.py, used as cwd fallback and to
# locate the CLI entry point when building job argv.
REPO_ROOT = Path(__file__).resolve().parent.parent

_ENV_HOME = "OMNIDOCK_WEBUI_HOME"


def app_home() -> Path:
    """Root of the web UI's own data directory."""
    override = os.environ.get(_ENV_HOME)
    if override:
        return Path(override).expanduser().resolve()
    return Path.home() / ".omnidock_webui"


def projects_file() -> Path:
    return app_home() / "projects.json"


def jobs_dir() -> Path:
    return app_home() / "jobs"


def logs_dir() -> Path:
    return app_home() / "logs"


def ensure_app_dirs() -> Path:
    """Create the app data directories if they do not exist yet."""
    home = app_home()
    home.mkdir(parents=True, exist_ok=True)
    jobs_dir().mkdir(parents=True, exist_ok=True)
    logs_dir().mkdir(parents=True, exist_ok=True)
    return home


# Binding defaults. Off-localhost binding requires an explicit opt-in flag
# on the CLI entry point -- there is no authentication in this feature, so
# the default must never be reachable from the network (FR-034).
HOST = "127.0.0.1"
PORT = 8770
