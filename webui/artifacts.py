"""Discovery and safe serving of a run's output artifacts.

Serving files from a user-supplied project root over HTTP introduces a
path-traversal risk the terminal app never had, so every request goes
through ``safe_resolve`` (FR-022).
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

# Canonical output topology (README "Output Topology"). Longest prefixes
# first so `interactions/prolif` wins over `interactions`.
_CATEGORY_PREFIXES: list[tuple[str, str]] = [
    ("interactions/pandamap", "interactions.pandamap"),
    ("interactions/prolif", "interactions.prolif"),
    ("interactions/ligplot", "interactions.ligplot"),
    ("interactions/poseview", "interactions.poseview"),
    ("interactions", "interactions"),
    ("rmsd_analysis", "rmsd_analysis"),
    ("3d_visualizations", "3d_visualizations"),
    ("visualizations", "visualizations"),
    ("best_poses", "best_poses"),
    ("complexes", "complexes"),
    ("raw_data", "raw_data"),
    ("analysis", "analysis"),
    ("reports", "reports"),
    ("run_tracking", "run_tracking"),
]

_MEDIA_KINDS: dict[str, str] = {
    ".csv": "table", ".tsv": "table",
    ".png": "image", ".svg": "image", ".jpg": "image", ".jpeg": "image", ".gif": "image",
    ".html": "html", ".htm": "html",
    ".pdb": "structure", ".sdf": "structure", ".pdbqt": "structure", ".mol2": "structure",
    ".md": "text", ".txt": "text", ".log": "text", ".json": "text", ".yaml": "text", ".yml": "text",
}

# Directories that are never interesting as results.
_SKIP_DIRS = {".git", "__pycache__", ".workflow", ".pytest_cache"}


class PathEscape(Exception):
    """Raised when a requested path resolves outside its project root."""


class ArtifactError(Exception):
    """Raised when artifacts or a run summary cannot be read."""


@dataclass
class Artifact:
    path: str  # project-relative, forward slashes
    category: str
    media_kind: str
    size: int
    modified_at: float

    def to_dict(self) -> dict:
        return {
            "path": self.path,
            "category": self.category,
            "media_kind": self.media_kind,
            "size": self.size,
            "modified_at": self.modified_at,
        }


def safe_resolve(project_root: str | Path, rel_path: str) -> Path:
    """Resolve ``rel_path`` inside ``project_root`` or raise PathEscape.

    Resolution happens *after* joining and follows symlinks, so a link
    inside the project pointing outside it is still caught.
    """
    root = Path(project_root).resolve()
    candidate = Path(str(rel_path).replace("\\", "/"))

    if candidate.is_absolute() or candidate.drive:
        raise PathEscape(f"Absolute paths are not allowed: {rel_path}")
    if any(part == ".." for part in candidate.parts):
        raise PathEscape(f"Path traversal is not allowed: {rel_path}")

    resolved = (root / candidate).resolve()
    try:
        resolved.relative_to(root)
    except ValueError:
        raise PathEscape(f"Path escapes the project root: {rel_path}") from None
    return resolved


def categorize(rel_path: str) -> str:
    lowered = rel_path.replace("\\", "/").lower()
    for prefix, category in _CATEGORY_PREFIXES:
        if f"/{prefix}/" in f"/{lowered}" or lowered.startswith(f"{prefix}/"):
            return category
    return "other"


def media_kind(rel_path: str) -> str:
    return _MEDIA_KINDS.get(Path(rel_path).suffix.lower(), "other")


def _as_artifact(root: Path, absolute: Path) -> Artifact | None:
    try:
        stat = absolute.stat()
    except OSError:
        return None
    rel = absolute.relative_to(root).as_posix()
    return Artifact(
        path=rel,
        category=categorize(rel),
        media_kind=media_kind(rel),
        size=stat.st_size,
        modified_at=stat.st_mtime,
    )


def find_runs(project_root: str | Path, max_depth: int = 4) -> list[dict]:
    """Locate directories containing a run_tracking/ folder."""
    root = Path(project_root)
    runs: list[dict] = []
    if not root.is_dir():
        return runs

    for tracking in root.rglob("run_tracking"):
        if not tracking.is_dir():
            continue
        run_dir = tracking.parent
        try:
            rel = run_dir.relative_to(root).as_posix() or "."
        except ValueError:
            continue
        if len(Path(rel).parts) > max_depth:
            continue
        try:
            modified = tracking.stat().st_mtime
        except OSError:
            modified = 0.0
        runs.append({
            "run_dir": rel,
            "has_manifest": (tracking / "run_manifest.json").is_file(),
            "modified_at": modified,
        })

    return sorted(runs, key=lambda r: r["modified_at"], reverse=True)


def load_outputs_index(project_root: str | Path, run_dir: str) -> list[Artifact]:
    """Artifacts for a run, from outputs_index.json when available.

    Falls back to walking the run directory, because a run that failed
    partway may never have written an index (FR-020).
    """
    root = Path(project_root).resolve()
    run_abs = safe_resolve(root, run_dir) if run_dir not in ("", ".") else root
    index_file = run_abs / "run_tracking" / "outputs_index.json"

    artifacts: list[Artifact] = []
    seen: set[str] = set()

    if index_file.is_file():
        try:
            rows = json.loads(index_file.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            rows = []
        for row in rows if isinstance(rows, list) else []:
            raw = None
            if isinstance(row, dict):
                for key in ("path", "file", "output", "output_path", "artifact"):
                    if row.get(key):
                        raw = str(row[key])
                        break
            elif isinstance(row, str):
                raw = row
            if not raw:
                continue

            candidate = Path(raw)
            absolute = candidate if candidate.is_absolute() else (run_abs / candidate)
            try:
                absolute = absolute.resolve()
                absolute.relative_to(root)  # index entries must stay in-project
            except (OSError, ValueError):
                continue
            if not absolute.is_file():
                continue
            artifact = _as_artifact(root, absolute)
            if artifact and artifact.path not in seen:
                seen.add(artifact.path)
                artifacts.append(artifact)

    if not artifacts:
        for absolute in run_abs.rglob("*"):
            if not absolute.is_file():
                continue
            if any(part in _SKIP_DIRS for part in absolute.parts):
                continue
            artifact = _as_artifact(root, absolute)
            if artifact and artifact.path not in seen:
                seen.add(artifact.path)
                artifacts.append(artifact)

    return sorted(artifacts, key=lambda a: (a.category, a.path))


def load_run_summary(project_root: str | Path, run_dir: str) -> dict:
    """Parse run_manifest.json into the browser-facing run summary."""
    run_abs = safe_resolve(project_root, run_dir) if run_dir not in ("", ".") else Path(project_root)
    manifest_file = run_abs / "run_tracking" / "run_manifest.json"
    if not manifest_file.is_file():
        return {"available": False, "engine": None, "stage_contract": [],
                "steps": [], "optional_features": []}

    try:
        manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        raise ArtifactError(f"Could not read run manifest: {exc}") from exc

    return {
        "available": True,
        "engine": manifest.get("engine") or manifest.get("engine_name"),
        "stage_contract": _stage_contract_rows(manifest),
        "steps": _step_rows(run_abs, manifest),
        "optional_features": _optional_feature_rows(manifest),
        "raw": manifest,
    }


def _stage_contract_rows(manifest: dict) -> list[dict]:
    """Normalize stage_contract into rows flagging enforced != requested.

    The manifest's shape has varied across specs, so both the mapping and
    the list-of-dicts forms are accepted.
    """
    contract = manifest.get("stage_contract") or {}
    rows: list[dict] = []

    if isinstance(contract, dict):
        for stage, value in contract.items():
            if isinstance(value, dict):
                requested = value.get("requested")
                enforced = value.get("enforced", value.get("effective"))
            else:
                requested = enforced = value
            rows.append({
                "stage": stage, "requested": requested, "enforced": enforced,
                "matches": requested == enforced,
            })
    elif isinstance(contract, list):
        for entry in contract:
            if not isinstance(entry, dict):
                continue
            requested = entry.get("requested")
            enforced = entry.get("enforced", entry.get("effective"))
            rows.append({
                "stage": entry.get("stage") or entry.get("name") or "?",
                "requested": requested, "enforced": enforced,
                "matches": requested == enforced,
            })

    return rows


def _step_rows(run_abs: Path, manifest: dict) -> list[dict]:
    steps = manifest.get("step_status") or manifest.get("steps")
    if steps is None:
        status_file = run_abs / "run_tracking" / "step_status.json"
        if status_file.is_file():
            try:
                steps = json.loads(status_file.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError):
                steps = []

    rows: list[dict] = []
    if isinstance(steps, list):
        for entry in steps:
            if isinstance(entry, dict):
                rows.append({
                    "name": entry.get("step") or entry.get("name") or "?",
                    "status": str(entry.get("status") or "unknown"),
                    "detail": str(entry.get("detail") or entry.get("note") or ""),
                })
    elif isinstance(steps, dict):
        for name, value in steps.items():
            status = value.get("status") if isinstance(value, dict) else value
            rows.append({"name": name, "status": str(status), "detail": ""})
    return rows


def _optional_feature_rows(manifest: dict) -> list[dict]:
    """Per-optional-feature classification (spec 029 FR-008)."""
    features = (
        manifest.get("optional_features")
        or manifest.get("optional_feature_status")
        or {}
    )
    rows: list[dict] = []
    if isinstance(features, dict):
        for name, value in features.items():
            if isinstance(value, dict):
                rows.append({
                    "name": name,
                    "classification": str(value.get("classification")
                                          or value.get("status") or "unknown"),
                    "detail": str(value.get("detail") or value.get("reason") or ""),
                })
            else:
                rows.append({"name": name, "classification": str(value), "detail": ""})
    elif isinstance(features, list):
        for entry in features:
            if isinstance(entry, dict):
                rows.append({
                    "name": entry.get("name") or "?",
                    "classification": str(entry.get("classification")
                                          or entry.get("status") or "unknown"),
                    "detail": str(entry.get("detail") or entry.get("reason") or ""),
                })
    return rows


def read_table(absolute: Path, page: int = 1, page_size: int = 100) -> dict:
    """Read a CSV/TSV into paginated rows without loading the whole file."""
    delimiter = "\t" if absolute.suffix.lower() == ".tsv" else ","
    page = max(1, int(page))
    page_size = max(1, min(int(page_size), 1000))
    start = (page - 1) * page_size

    columns: list[str] = []
    rows: list[list] = []
    total = 0

    try:
        with open(absolute, "r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            for index, record in enumerate(reader):
                if index == 0:
                    columns = record
                    continue
                total += 1
                if start <= total - 1 < start + page_size:
                    rows.append(record)
    except OSError as exc:
        raise ArtifactError(f"Could not read table: {exc}") from exc

    return {
        "columns": columns, "rows": rows,
        "page": page, "page_size": page_size, "total_rows": total,
    }
