"""Read-only bridge onto ``workflow.state``.

Two spec requirements meet here and pull in opposite directions:

* FR-002 -- read project state only through the ``workflow.state``
  accessors, never by parsing ``.workflow/state.json`` ourselves.
* FR-035 -- the web layer must never write into a project directory.

The accessors (`ensure_state`, `load_state`, `list_steps`,
`summarize_state`, `list_background_tasks`) all *create* the state file
when it is missing: `ensure_state` writes defaults, and `load_state`
delegates to it. Calling them against an uninitialized directory would
therefore write.

Resolution: guard on the state file's existence first. When it exists the
accessors only read, so FR-002 holds. When it does not, report
``has_workflow: False`` and leave initialization to the explicit
``workflow init`` action, which runs as a CLI subprocess (subprocesses may
write; the web process may not).
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

from workflow.interactive import TIMELINE_STEPS
from workflow.state import (
    list_background_tasks,
    list_steps,
    load_state,
    state_path,
    summarize_state,
)

# How long to wait before retrying a read that failed. The writer is atomic
# (`_atomic_dump_json`), but a reader can still catch a rename mid-flight on
# some filesystems.
_RETRY_DELAY_SECONDS = 0.15


class StateUnreadable(Exception):
    """Raised when project state exists but could not be read."""


@dataclass
class Timeline:
    steps: list[dict]


def has_workflow(root: Path) -> bool:
    """Whether this project has been initialized, without creating anything."""
    try:
        return state_path(Path(root)).is_file()
    except OSError:
        return False


def _read_with_retry(fn, root: Path, what: str):
    """Call a state accessor, retrying once on a transient read failure."""
    try:
        return fn(root)
    except Exception:
        time.sleep(_RETRY_DELAY_SECONDS)
        try:
            return fn(root)
        except Exception as exc:
            raise StateUnreadable(f"Could not read {what} for {root}: {exc}") from exc


def project_state(root: Path) -> dict:
    """Context, flags, artifacts, and a human summary for one project.

    Returns ``has_workflow: False`` and empty sections for an uninitialized
    directory rather than initializing it.
    """
    root = Path(root)
    if not has_workflow(root):
        return {
            "has_workflow": False,
            "current_context": {},
            "feature_flags": {},
            "artifacts": {},
            "summary": [],
            "background_tasks": [],
        }

    payload = _read_with_retry(load_state, root, "workflow state")
    summary = _read_with_retry(summarize_state, root, "state summary")
    tasks = _read_with_retry(
        lambda r: list_background_tasks(r, include_finished=True), root, "background tasks"
    )

    return {
        "has_workflow": True,
        "current_context": payload.get("current_context", {}) or {},
        "feature_flags": payload.get("feature_flags", {}) or {},
        "artifacts": payload.get("artifacts", {}) or {},
        "summary": list(summary or []),
        # Interactive-session tasks. Displayed for context only -- web jobs
        # live in webui.jobs and are never sourced from here (data-model.md).
        "background_tasks": list(tasks or []),
    }


def timeline(root: Path) -> list[dict]:
    """The workflow timeline in TIMELINE_STEPS order, joined with recorded status.

    Order and labels come from ``workflow.interactive.TIMELINE_STEPS`` so the
    web UI cannot drift from the terminal UI (FR-003, FR-007).
    """
    root = Path(root)
    recorded: dict = {}
    if has_workflow(root):
        recorded = _read_with_retry(list_steps, root, "workflow steps") or {}

    rows = []
    for key, label in TIMELINE_STEPS:
        step = recorded.get(key) or {}
        rows.append(
            {
                "key": key,
                "label": label,
                "status": str(step.get("status") or "pending"),
                "timestamp": step.get("finished_at") or step.get("started_at") or None,
                "note": "; ".join(str(n) for n in (step.get("notes") or [])),
            }
        )
    return rows
