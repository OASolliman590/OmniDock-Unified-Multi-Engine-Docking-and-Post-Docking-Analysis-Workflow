"""
Lightweight artifact-DAG executor for post-docking analysis.

This module provides the execution spine for Spec 024. Nodes declare the
artifacts they consume and produce, and the graph executor derives tiered
execution, cache hits, and failure isolation from those contracts.
"""

from __future__ import annotations

import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Set


NODE_STATUSES = {
    "pending",
    "running",
    "completed",
    "completed_with_warnings",
    "failed",
    "skipped_optional",
    "blocked_by_failure",
    "cache_hit",
}


@dataclass(frozen=True)
class ArtifactNode:
    name: str
    inputs: List[str]
    outputs: List[str]
    compute: Callable[[], object]
    optional: bool = False
    cacheable: bool = True


class ArtifactGraph:
    def __init__(
        self,
        *,
        cache_file: str | Path,
        report_file: Optional[str | Path] = None,
        max_workers: int = 4,
    ) -> None:
        self.cache_file = Path(cache_file).expanduser().resolve()
        self.report_file = (
            Path(report_file).expanduser().resolve()
            if report_file is not None
            else self.cache_file.with_name("dag_execution_report.json")
        )
        self.max_workers = max(1, int(max_workers or 1))
        self.nodes: Dict[str, ArtifactNode] = {}
        self.output_to_node: Dict[str, str] = {}
        self.last_execution_report: Dict[str, object] = {}

    def register(self, node: ArtifactNode) -> None:
        if node.name in self.nodes:
            raise ValueError(f"Duplicate node registration: {node.name}")
        if not node.outputs:
            raise ValueError(f"Node {node.name} must declare at least one output")
        for output in node.outputs:
            if output in self.output_to_node:
                raise ValueError(
                    f"Artifact output already registered: {output} "
                    f"(owner={self.output_to_node[output]})"
                )
        self.nodes[node.name] = node
        for output in node.outputs:
            self.output_to_node[output] = node.name

    def request(self, artifact: str, force: bool = False) -> Dict[str, object]:
        producer = self.output_to_node.get(artifact)
        if not producer:
            raise KeyError(f"No node produces requested artifact: {artifact}")
        tiers = self._topo_sort(artifact)
        report = self._execute(tiers, artifact=artifact, force=force)
        self.last_execution_report = report
        self._write_report(report)
        return report

    def _collect_required_nodes(self, artifact: str, seen: Optional[Set[str]] = None) -> Set[str]:
        producer = self.output_to_node.get(artifact)
        if not producer:
            return set()
        if seen is None:
            seen = set()
        if producer in seen:
            return seen
        seen.add(producer)
        node = self.nodes[producer]
        for input_artifact in node.inputs:
            upstream = self.output_to_node.get(input_artifact)
            if upstream:
                self._collect_required_nodes(input_artifact, seen)
        return seen

    def _topo_sort(self, target: str) -> List[List[ArtifactNode]]:
        required_names = self._collect_required_nodes(target)
        if not required_names:
            raise KeyError(f"No dependency subgraph found for target: {target}")

        deps: Dict[str, Set[str]] = {}
        reverse: Dict[str, Set[str]] = {name: set() for name in required_names}
        for name in required_names:
            node = self.nodes[name]
            upstream_names = {
                self.output_to_node[input_artifact]
                for input_artifact in node.inputs
                if input_artifact in self.output_to_node and self.output_to_node[input_artifact] in required_names
            }
            deps[name] = set(upstream_names)
            for upstream_name in upstream_names:
                reverse.setdefault(upstream_name, set()).add(name)

        tiers: List[List[ArtifactNode]] = []
        ready = sorted(name for name, node_deps in deps.items() if not node_deps)
        visited: Set[str] = set()
        while ready:
            current_names = list(ready)
            ready = []
            tiers.append([self.nodes[name] for name in current_names])
            for name in current_names:
                visited.add(name)
                for dependent in sorted(reverse.get(name, set())):
                    if dependent in visited:
                        continue
                    deps[dependent].discard(name)
                    if not deps[dependent]:
                        ready.append(dependent)
            ready = sorted(set(ready))

        if visited != required_names:
            unresolved = ", ".join(sorted(required_names - visited))
            raise ValueError(f"Cycle or unresolved dependency detected in artifact graph: {unresolved}")
        return tiers

    def _load_cache(self) -> Dict[str, object]:
        if not self.cache_file.exists():
            return {"nodes": {}}
        try:
            payload = json.loads(self.cache_file.read_text(encoding="utf-8"))
        except Exception:
            return {"nodes": {}}
        if not isinstance(payload, dict):
            return {"nodes": {}}
        if not isinstance(payload.get("nodes"), dict):
            payload["nodes"] = {}
        return payload

    def _save_cache(self, payload: Dict[str, object]) -> None:
        self.cache_file.parent.mkdir(parents=True, exist_ok=True)
        self.cache_file.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    @staticmethod
    def _utc_now_iso() -> str:
        return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")

    @staticmethod
    def _artifact_path_rows(artifacts: Iterable[str]) -> List[Path]:
        rows: List[Path] = []
        for artifact in artifacts:
            rows.append(Path(str(artifact)).expanduser().resolve())
        return rows

    def _compute_cache_key(self, node: ArtifactNode) -> str:
        rows: List[str] = []
        for path in sorted(self._artifact_path_rows(node.inputs), key=lambda item: str(item)):
            if path.exists():
                try:
                    stat = path.stat()
                    rows.append(f"{path}|exists|{stat.st_mtime_ns}|{stat.st_size}")
                except OSError:
                    rows.append(f"{path}|exists|stat_error")
            else:
                rows.append(f"{path}|missing")
        return hashlib.sha256("\n".join(rows).encode("utf-8")).hexdigest()

    def _outputs_exist(self, node: ArtifactNode) -> bool:
        return all(path.exists() for path in self._artifact_path_rows(node.outputs))

    def _is_cached(self, node: ArtifactNode, cache_payload: Dict[str, object]) -> bool:
        if not node.cacheable:
            return False
        node_cache = ((cache_payload.get("nodes") or {}).get(node.name) or {}) if isinstance(cache_payload, dict) else {}
        expected_key = str(node_cache.get("cache_key") or "")
        if not expected_key:
            return False
        if not self._outputs_exist(node):
            return False
        return expected_key == self._compute_cache_key(node)

    def _write_report(self, report: Dict[str, object]) -> None:
        self.report_file.parent.mkdir(parents=True, exist_ok=True)
        self.report_file.write_text(json.dumps(report, indent=2), encoding="utf-8")

    def _execute(self, tiers: List[List[ArtifactNode]], *, artifact: str, force: bool) -> Dict[str, object]:
        cache_payload = self._load_cache()
        cache_nodes = cache_payload.setdefault("nodes", {})
        cache_payload["generated_at"] = self._utc_now_iso()
        cache_payload["cache_file"] = str(self.cache_file)
        node_records: Dict[str, Dict[str, object]] = {}
        required_names = {node.name for tier in tiers for node in tier}
        dependency_failures: Dict[str, Set[str]] = {name: set() for name in required_names}
        wall_started = time.perf_counter()
        tier_reports: List[Dict[str, object]] = []

        def _record(name: str) -> Dict[str, object]:
            existing = node_records.get(name)
            if existing:
                return existing
            entry = {
                "name": name,
                "status": "pending",
                "started_at": "",
                "ended_at": "",
                "duration_s": 0.0,
                "cache_hit": False,
                "optional": bool(self.nodes[name].optional),
                "blocked_by": [],
                "details": "",
            }
            node_records[name] = entry
            return entry

        for tier_index, tier in enumerate(tiers):
            tier_started = time.perf_counter()
            tier_node_names = [node.name for node in tier]
            future_map = {}
            tier_report = {"tier_index": tier_index, "nodes": tier_node_names}

            with ThreadPoolExecutor(max_workers=min(self.max_workers, max(1, len(tier)))) as pool:
                for node in tier:
                    entry = _record(node.name)
                    blocking = sorted(dependency_failures.get(node.name, set()))
                    if blocking:
                        entry["status"] = "blocked_by_failure"
                        entry["blocked_by"] = blocking
                        entry["details"] = "blocked by dependency failure"
                        entry["started_at"] = self._utc_now_iso()
                        entry["ended_at"] = entry["started_at"]
                        continue
                    if not force and self._is_cached(node, cache_payload):
                        entry["status"] = "cache_hit"
                        entry["cache_hit"] = True
                        entry["started_at"] = self._utc_now_iso()
                        entry["ended_at"] = entry["started_at"]
                        entry["duration_s"] = 0.0
                        continue

                    entry["status"] = "running"
                    entry["started_at"] = self._utc_now_iso()
                    future_map[pool.submit(node.compute)] = (node, time.perf_counter())

                for future in as_completed(future_map):
                    node, started_at = future_map[future]
                    entry = _record(node.name)
                    duration = time.perf_counter() - started_at
                    entry["duration_s"] = round(duration, 6)
                    entry["ended_at"] = self._utc_now_iso()
                    try:
                        result = future.result()
                    except Exception as exc:
                        if node.optional:
                            entry["status"] = "skipped_optional"
                            entry["details"] = str(exc)
                        else:
                            entry["status"] = "failed"
                            entry["details"] = str(exc)
                    else:
                        warnings_present = False
                        details = ""
                        if isinstance(result, dict):
                            warnings_present = bool(result.get("warnings"))
                            details = str(result.get("details") or "")
                        elif result not in (None, True, False):
                            details = str(result)
                        entry["status"] = "completed_with_warnings" if warnings_present else "completed"
                        entry["details"] = details
                        if node.cacheable:
                            cache_nodes[node.name] = {
                                "cache_key": self._compute_cache_key(node),
                                "outputs": list(node.outputs),
                                "updated_at": self._utc_now_iso(),
                            }

            for node in tier:
                entry = _record(node.name)
                if entry["status"] == "failed":
                    for dependent_name in required_names:
                        dependent = self.nodes[dependent_name]
                        if node.name in {
                            self.output_to_node[input_artifact]
                            for input_artifact in dependent.inputs
                            if input_artifact in self.output_to_node
                        }:
                            dependency_failures.setdefault(dependent_name, set()).add(node.name)

            tier_wall = time.perf_counter() - tier_started
            tier_sum = sum(float(_record(node.name).get("duration_s", 0.0) or 0.0) for node in tier)
            tier_report["wall_clock_s"] = round(tier_wall, 6)
            tier_report["sum_node_duration_s"] = round(tier_sum, 6)
            tier_report["parallel_time_saved_s"] = round(max(0.0, tier_sum - tier_wall), 6)
            tier_reports.append(tier_report)

        wall_clock = time.perf_counter() - wall_started
        parallel_saved = round(sum(float(item.get("parallel_time_saved_s", 0.0) or 0.0) for item in tier_reports), 6)
        status_counts: Dict[str, int] = {}
        for entry in node_records.values():
            status = str(entry.get("status") or "pending")
            status_counts[status] = status_counts.get(status, 0) + 1

        report = {
            "generated_at": self._utc_now_iso(),
            "requested_artifact": artifact,
            "report_file": str(self.report_file),
            "cache_file": str(self.cache_file),
            "force": bool(force),
            "node_count": len(required_names),
            "wall_clock_s": round(wall_clock, 6),
            "parallel_time_saved_s": parallel_saved,
            "status_counts": status_counts,
            "tiers": tier_reports,
            "nodes": node_records,
        }
        self._save_cache(cache_payload)
        return report


__all__ = ["ArtifactGraph", "ArtifactNode", "NODE_STATUSES"]
