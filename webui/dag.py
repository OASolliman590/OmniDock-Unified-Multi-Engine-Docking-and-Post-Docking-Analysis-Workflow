"""View model for the artifact DAG (spec 024).

Reads the report written by ``ArtifactGraph._write_report`` plus the cache
file it references. The status vocabulary is imported from
``artifact_graph`` rather than restated, so a new status cannot silently
render as "unknown" (FR-025).

Note on edges: the execution report records per-node status but not the
nodes' declared inputs/outputs -- those live on the ArtifactNode objects,
which exist only inside the process that ran the graph. What the report
does carry is the tier grouping (nodes within a tier have no
inter-dependencies, and each tier depends on the ones before it) and
per-node ``blocked_by``. The view is therefore tier-layered with explicit
blocking edges, which is what the data actually supports.
"""

from __future__ import annotations

import json
from pathlib import Path

from post_docking_analysis.artifact_graph import NODE_STATUSES

# Default filenames from ArtifactGraph: the cache lives at
# 4-Working/dag_cache.json and the report sits beside it.
CACHE_FILENAME = "dag_cache.json"
REPORT_FILENAME = "dag_execution_report.json"

# Search order for a project's DAG report.
_REPORT_LOCATIONS = [
    Path("4-Working") / REPORT_FILENAME,
    Path(REPORT_FILENAME),
    Path("run_tracking") / REPORT_FILENAME,
]


class DagError(Exception):
    """Raised when a DAG report exists but cannot be read."""


def find_report(project_root: str | Path, run_dir: str | None = None) -> Path | None:
    """Locate a DAG execution report for a project or one of its runs."""
    root = Path(project_root)
    bases = [root]
    if run_dir and run_dir not in ("", "."):
        bases.insert(0, root / run_dir)

    for base in bases:
        for relative in _REPORT_LOCATIONS:
            candidate = base / relative
            if candidate.is_file():
                return candidate

    # Fall back to a bounded search rather than missing a non-standard layout.
    try:
        for candidate in root.rglob(REPORT_FILENAME):
            if candidate.is_file():
                return candidate
    except OSError:
        pass
    return None


def load_report(report_path: str | Path) -> dict:
    path = Path(report_path)
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        raise DagError(f"Could not read DAG report {path}: {exc}") from exc


def _load_cache(report: dict, report_path: Path) -> dict:
    """Cache payload referenced by the report, for cache-hit detail."""
    cache_file = report.get("cache_file")
    candidates = []
    if cache_file:
        candidates.append(Path(cache_file))
    candidates.append(report_path.with_name(CACHE_FILENAME))

    for candidate in candidates:
        try:
            if candidate.is_file():
                return json.loads(candidate.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            continue
    return {}


def to_view_model(report: dict, cache: dict | None = None) -> dict:
    """Normalize a raw report into what the template renders."""
    cache_nodes = (cache or {}).get("nodes", {}) or {}
    raw_nodes = report.get("nodes", {}) or {}

    nodes: dict[str, dict] = {}
    for name, entry in raw_nodes.items():
        if not isinstance(entry, dict):
            continue
        status = str(entry.get("status") or "pending")
        cached = bool(entry.get("cache_hit")) or status == "cache_hit"
        nodes[name] = {
            "name": name,
            "status": status,
            # An unrecognized status means artifact_graph gained one this
            # module has not been taught to style; surface it, don't hide it.
            "known_status": status in NODE_STATUSES,
            "cached": cached,
            "optional": bool(entry.get("optional")),
            "blocked_by": list(entry.get("blocked_by") or []),
            "details": str(entry.get("details") or ""),
            "duration_s": float(entry.get("duration_s") or 0.0),
            "cache_key": (cache_nodes.get(name) or {}).get("cache_key", ""),
        }

    tiers = []
    for tier in report.get("tiers", []) or []:
        if not isinstance(tier, dict):
            continue
        names = [n for n in (tier.get("nodes") or []) if n in nodes]
        tiers.append({
            "tier_index": tier.get("tier_index"),
            "nodes": [nodes[n] for n in names],
            "wall_clock_s": tier.get("wall_clock_s"),
            "parallel_time_saved_s": tier.get("parallel_time_saved_s"),
        })

    # Any node missing from the tier listing still has to be shown.
    listed = {n["name"] for tier in tiers for n in tier["nodes"]}
    orphans = [node for name, node in nodes.items() if name not in listed]
    if orphans:
        tiers.append({"tier_index": None, "nodes": orphans,
                      "wall_clock_s": None, "parallel_time_saved_s": None})

    edges = [
        {"from": blocker, "to": node["name"], "kind": "blocked_by"}
        for node in nodes.values()
        for blocker in node["blocked_by"]
        if blocker in nodes
    ]

    counts = report.get("status_counts") or {}
    if not counts:
        counts = {}
        for node in nodes.values():
            counts[node["status"]] = counts.get(node["status"], 0) + 1

    return {
        "available": True,
        "generated_at": report.get("generated_at", ""),
        "requested_artifact": report.get("requested_artifact", ""),
        "node_count": report.get("node_count", len(nodes)),
        "wall_clock_s": report.get("wall_clock_s"),
        "parallel_time_saved_s": report.get("parallel_time_saved_s"),
        "status_counts": counts,
        "known_statuses": sorted(NODE_STATUSES),
        "tiers": tiers,
        "edges": edges,
    }


def load_dag(project_root: str | Path, run_dir: str | None = None) -> dict:
    """Find, read, and normalize a project's DAG report."""
    report_path = find_report(project_root, run_dir)
    if report_path is None:
        return {"available": False, "tiers": [], "edges": [],
                "status_counts": {}, "known_statuses": sorted(NODE_STATUSES)}
    report = load_report(report_path)
    cache = _load_cache(report, report_path)
    model = to_view_model(report, cache)
    model["report_file"] = str(report_path)
    return model
