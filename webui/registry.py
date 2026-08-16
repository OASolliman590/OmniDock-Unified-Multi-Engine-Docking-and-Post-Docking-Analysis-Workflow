"""Registry of OmniDock project directories known to the web UI.

The registry is the web layer's own data -- it lives in the app data
directory, never inside a project (FR-035). A project is identified by the
canonical form of its path, so the same directory cannot be registered
twice under different spellings.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

from .config import projects_file


class RegistryError(Exception):
    """Raised when a project cannot be registered or found."""


@dataclass
class Project:
    id: str
    path: str
    name: str
    registered_at: str

    @property
    def root(self) -> Path:
        return Path(self.path)

    @property
    def available(self) -> bool:
        """Whether the directory still exists and is readable.

        Never raises: a project whose directory was deleted or unmounted is
        reported as unavailable rather than breaking the page (FR-030).
        """
        try:
            return self.root.is_dir() and os.access(self.root, os.R_OK)
        except OSError:
            return False

    @property
    def has_workflow(self) -> bool:
        try:
            return (self.root / ".workflow" / "state.json").is_file()
        except OSError:
            return False

    def to_dict(self) -> dict:
        data = asdict(self)
        data["available"] = self.available
        data["has_workflow"] = self.has_workflow
        return data


def canonicalize(path: str | os.PathLike) -> Path:
    """Reduce a user-supplied path to one canonical form.

    Expands ``~``, makes it absolute, resolves symlinks, and normalizes case
    on Windows -- so "C:\\Proj", "c:\\proj\\", and a symlink to it all map to
    a single registry entry.
    """
    resolved = Path(path).expanduser()
    try:
        resolved = resolved.resolve(strict=False)
    except OSError as exc:
        raise RegistryError(f"Could not resolve path: {path} ({exc})") from exc

    if sys.platform == "win32":
        resolved = Path(os.path.normcase(str(resolved)))
    return resolved


def project_id_for(canonical: Path) -> str:
    """Stable id derived from the canonical path (see data-model.md)."""
    digest = hashlib.sha256(str(canonical).encode("utf-8")).hexdigest()
    return digest[:8]


def _read_all() -> list[Project]:
    path = projects_file()
    if not path.exists():
        return []
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return []
    projects = []
    for entry in raw:
        try:
            projects.append(Project(**entry))
        except TypeError:
            continue
    return projects


def _write_all(projects: list[Project]) -> None:
    path = projects_file()
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = [
        {"id": p.id, "path": p.path, "name": p.name, "registered_at": p.registered_at}
        for p in projects
    ]
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(path)


def register(path: str | os.PathLike, name: str | None = None) -> Project:
    """Register a project directory, or return the existing entry for it."""
    canonical = canonicalize(path)

    if not canonical.exists():
        raise RegistryError(f"Directory does not exist: {canonical}")
    if not canonical.is_dir():
        raise RegistryError(f"Not a directory: {canonical}")
    if not os.access(canonical, os.R_OK):
        raise RegistryError(f"Directory is not readable: {canonical}")

    pid = project_id_for(canonical)
    projects = _read_all()
    for existing in projects:
        if existing.id == pid:
            return existing

    project = Project(
        id=pid,
        path=str(canonical),
        name=(name or "").strip() or canonical.name or str(canonical),
        registered_at=datetime.now(timezone.utc).isoformat(),
    )
    projects.append(project)
    _write_all(projects)
    return project


def list_projects() -> list[Project]:
    return sorted(_read_all(), key=lambda p: p.name.lower())


def get_project(project_id: str) -> Project:
    for project in _read_all():
        if project.id == project_id:
            return project
    raise RegistryError(f"No registered project with id {project_id}")


def unregister(project_id: str) -> None:
    """Remove a project from the registry. Never touches the directory itself."""
    projects = _read_all()
    remaining = [p for p in projects if p.id != project_id]
    if len(remaining) == len(projects):
        raise RegistryError(f"No registered project with id {project_id}")
    _write_all(remaining)
