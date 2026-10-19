"""Read-only protonation view for the web UI (Spec 038 R4).

The view is built by ``workflow.interactive.build_protonation_view``, the function the interactive
shell prints, and that reads through the Spec 037 resolver in ``docking.project_layout``. This module
adds no policy logic and writes nothing. Setting an override is a launched CLI job
(``workflow protonation-policy``, see ``targets.PROTONATION_TARGET``), like every other write in the UI.
"""

from __future__ import annotations

from pathlib import Path

from workflow.interactive import PROTONATION_CLEAR_NOTE, PROTONATION_STALE_WARNING, build_protonation_view

TARGET = "workflow.protonation-policy"
STALE_WARNING = PROTONATION_STALE_WARNING
CLEAR_NOTE = PROTONATION_CLEAR_NOTE


def read_view(project_root: Path) -> dict:
    """The resolved policy of a project as JSON-ready data. Opens nothing for writing."""
    try:
        return build_protonation_view(Path(project_root))
    except (OSError, ValueError) as exc:
        return {"has_manifest": False, "error": f"protonation_view_failed: {exc}", "stale": [], "receptors": [], "ligands": []}
