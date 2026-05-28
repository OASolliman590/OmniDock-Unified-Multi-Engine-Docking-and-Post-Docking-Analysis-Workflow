"""
Unified post-docking analysis entrypoint.

This module defines the canonical pipeline class used by workflow and CLI.
Legacy split implementations remain importable for compatibility, but user-facing
entrypoints should route through this module.
"""

from __future__ import annotations

from typing import Dict, List, Optional

from post_docking_analysis.multi_engine_pipeline import (
    UNIFIED_COMPAT_COLUMNS,
    MultiEngineAnalysisPipeline,
)
from post_docking_analysis.engine_detector import detect_engines


def _normalize_engine_tokens(raw_value: object) -> List[str]:
    if raw_value is None:
        return []
    if isinstance(raw_value, str):
        raw_items = [token.strip().lower() for token in raw_value.split(",")]
    else:
        raw_items = [str(token).strip().lower() for token in raw_value]
    normalized: List[str] = []
    for token in raw_items:
        if token and token not in normalized:
            normalized.append(token)
    return normalized


def _resolve_engine_scope(
    *,
    requested_engines: object,
    preset_name: Optional[str],
    manifest: Dict[str, object],
    detection: Dict[str, object],
    override_engine: Optional[str],
) -> Dict[str, object]:
    valid_engines = [
        str(item).strip().lower()
        for item in (detection.get("valid_engines") or [])
        if str(item).strip()
    ]
    preset_map = manifest.get("engine_presets") if isinstance(manifest.get("engine_presets"), dict) else {}
    preset_value = None
    scope_source = "detected_all"

    if preset_name:
        preset_key = str(preset_name).strip()
        if preset_key not in preset_map:
            raise ValueError(f"Unknown engine preset: {preset_key}")
        preset_value = preset_map.get(preset_key)
        scope_source = f"preset:{preset_key}"

    requested = _normalize_engine_tokens(requested_engines if requested_engines not in (None, "") else preset_value)
    if override_engine:
        requested = [override_engine]
        scope_source = "engine_override"

    if not requested:
        requested = list(valid_engines)

    absent = [engine for engine in requested if engine not in valid_engines]
    if absent:
        raise ValueError(
            "Requested engines are not available for analysis in this project: "
            + ", ".join(absent)
        )

    excluded = [engine for engine in valid_engines if engine not in requested]
    exclusion_reason = (
        "cli_override" if scope_source == "engine_override" else "user_excluded" if excluded else ""
    )
    excluded_payload = [
        {"engine": engine, "exclusion_reason": exclusion_reason or "user_excluded"}
        for engine in excluded
    ]

    scoped_detection = dict(detection)
    scoped_engines_payload = {}
    for engine, payload in (detection.get("engines") or {}).items():
        row = dict(payload or {})
        if engine in excluded and row.get("status") != "absent":
            row["excluded_this_session"] = True
        scoped_engines_payload[engine] = row
    scoped_detection["engines"] = scoped_engines_payload

    return {
        "engines_in_scope": requested,
        "excluded_engines": excluded_payload,
        "scope_source": scope_source,
        "preset_name": str(preset_name or "").strip(),
        "scoped_detection": scoped_detection,
        "scoped_engine_count": len(requested),
        "detected_engine_count": len(valid_engines),
    }


class UnifiedPostDockingPipeline(MultiEngineAnalysisPipeline):
    """Canonical unified post-docking pipeline."""

    def __init__(self, *args, **kwargs):
        project_dir = kwargs.get("project_dir")
        dag_scope = str(kwargs.pop("dag_scope", "") or "").strip().lower()
        dag_force = bool(kwargs.pop("dag_force", False))
        requested_mode = str(kwargs.get("analysis_mode") or "auto").strip().lower()
        override_engine = str(kwargs.get("engine") or "").strip().lower() or None
        favorite_engine = str(kwargs.get("favorite_engine") or "").strip().lower() or None
        requested_engines = kwargs.pop("engines", None)
        engine_preset = str(kwargs.pop("engine_preset", "") or "").strip()
        redetect = bool(kwargs.pop("redetect", False))
        interactive_fallback = bool(kwargs.pop("interactive_fallback", False))
        manifest = {}
        try:
            from docking.project_layout import load_manifest
            manifest = load_manifest(project_dir) if project_dir else {}
        except Exception:
            manifest = {}

        detection = detect_engines(
            project_dir,
            override_engine=override_engine,
            redetect=redetect,
        ) if project_dir else {
            "routing_decision": "single_engine" if override_engine else "multi_engine",
            "routed_to_engine": override_engine,
            "valid_engines": [override_engine] if override_engine else [],
            "detection_override": bool(override_engine),
        }

        routing_decision = str(detection.get("routing_decision") or "")
        valid_engines = [str(item).strip().lower() for item in (detection.get("valid_engines") or []) if str(item).strip()]
        routed_to_engine = str(detection.get("routed_to_engine") or "").strip().lower() or None

        explicit_engine_hint = override_engine or favorite_engine or str(kwargs.get("engine") or "").strip().lower() or None

        if routing_decision == "no_valid_engines" and explicit_engine_hint:
            detection["routing_decision"] = "single_engine"
            detection["routed_to_engine"] = explicit_engine_hint
            detection["detection_override"] = bool(override_engine)
            routing_decision = "single_engine"
            routed_to_engine = explicit_engine_hint

        if routing_decision == "no_valid_engines" and not interactive_fallback:
            raise ValueError(
                "No valid docking-engine outputs were detected for this project. "
                "Check 4-Working/metadata/engine_detection_report.json for coverage and parsing details, "
                "or rerun with a valid --engine override once outputs exist."
            )

        scope_resolution = _resolve_engine_scope(
            requested_engines=requested_engines,
            preset_name=engine_preset or None,
            manifest=manifest,
            detection=detection,
            override_engine=override_engine,
        )
        scoped_engines = list(scope_resolution["engines_in_scope"])
        kwargs["engines_in_scope"] = scoped_engines
        kwargs["excluded_engines"] = list(scope_resolution["excluded_engines"])
        kwargs["scope_source"] = str(scope_resolution["scope_source"])
        kwargs["engine_preset_name"] = str(scope_resolution["preset_name"])
        kwargs["detected_engine_count"] = int(scope_resolution["detected_engine_count"])

        effective_mode = requested_mode
        if override_engine:
            effective_mode = "single_engine"
            kwargs["engine"] = override_engine
            kwargs["favorite_engine"] = override_engine
        elif len(scoped_engines) == 1:
            scoped_engine = scoped_engines[0]
            effective_mode = "single_engine"
            kwargs["engine"] = scoped_engine
            kwargs["favorite_engine"] = scoped_engine
        elif requested_mode in {"", "auto"}:
            effective_mode = routing_decision if routing_decision in {"single_engine", "multi_engine"} else "comparative_all_engines"
        elif requested_mode == "comparative_all_engines" and routing_decision == "single_engine" and routed_to_engine:
            effective_mode = "single_engine"
            kwargs["engine"] = routed_to_engine
            kwargs["favorite_engine"] = routed_to_engine
        elif requested_mode == "favorite_engine_continue" and not favorite_engine and routed_to_engine:
            kwargs["favorite_engine"] = routed_to_engine
        elif requested_mode == "single_engine" and not override_engine and not kwargs.get("engine") and routed_to_engine:
            kwargs["engine"] = routed_to_engine

        if effective_mode == "multi_engine":
            effective_mode = "comparative_all_engines"
        kwargs["analysis_mode"] = effective_mode
        self.engine_detection_report = dict(scope_resolution["scoped_detection"])
        self.analysis_mode_requested = requested_mode
        super().__init__(*args, **kwargs)
        self.engine_detection_report = dict(scope_resolution["scoped_detection"])
        self.manifest.setdefault("detected_engines", detection)
        report_file = self.engine_detection_report.get("report_file")
        if report_file:
            try:
                from pathlib import Path
                Path(str(report_file)).write_text(__import__("json").dumps(self.engine_detection_report, indent=2), encoding="utf-8")
            except Exception:
                pass
        self.analysis_mode_requested = requested_mode
        self.engines_in_scope = scoped_engines
        self.excluded_engines = list(scope_resolution["excluded_engines"])
        self.scope_source = str(scope_resolution["scope_source"])
        self.engine_preset_name = str(scope_resolution["preset_name"])
        self.detected_engine_count = int(scope_resolution["detected_engine_count"])
        self.dag_scope = dag_scope
        self.dag_force = dag_force
        self.dag_execution_report = {}
        self.dag_execution_report_file = ""
        if effective_mode == "single_engine":
            detected_engine = str(kwargs.get("engine") or kwargs.get("favorite_engine") or routed_to_engine or "").strip().lower()
            if detected_engine:
                self.engine = detected_engine
                self.favorite_engine = detected_engine
        self.single_engine_mode = self.analysis_mode == "single_engine"

    def _effective_dag_scope(self) -> str:
        requested = str(self.dag_scope or "").strip().lower()
        if requested:
            return requested
        if self.analysis_mode == "single_engine":
            return "comparison_only"
        fallback_scope = str(getattr(self, "analysis_scope", "") or "full").strip().lower()
        if fallback_scope in {"rescoring_only", "qc_only"}:
            return "report_only"
        return fallback_scope or "full"

    def run(self) -> bool:
        scope = self._effective_dag_scope()
        run_started_at = self._utc_now_iso()
        report = self.request_artifact_scope(scope, force=self.dag_force)
        self.dag_execution_report = dict(report)
        self.dag_execution_report_file = str(self._dag_artifact_paths()["dag_execution_report"])
        failed_statuses = {"failed", "blocked_by_failure"}
        node_rows = report.get("nodes") or {}
        succeeded = not any(
            str((payload or {}).get("status") or "") in failed_statuses
            for payload in node_rows.values()
        )
        if succeeded and self.analysis_mode == "single_engine":
            target_engine = str(self.engine or self.favorite_engine or "").strip().lower()
            if target_engine:
                scores = self._load_or_build_scores()
                self._write_single_engine_reports(scores, target_engine, self.output_dir / target_engine)
        try:
            step_states = [
                {
                    "step": str(name),
                    "status": str((payload or {}).get("status") or ""),
                    "started_at": str((payload or {}).get("started_at") or ""),
                    "ended_at": str((payload or {}).get("ended_at") or ""),
                    "details": str((payload or {}).get("details") or ""),
                }
                for name, payload in node_rows.items()
            ]
            failed_nodes = [
                str(name)
                for name, payload in node_rows.items()
                if str((payload or {}).get("status") or "") in failed_statuses
            ]
            self._write_run_tracking_artifacts(
                run_started_at=run_started_at,
                run_status="validated" if succeeded else "failed",
                error_message="; ".join(failed_nodes),
                step_states=step_states,
            )
        except Exception:
            pass
        return succeeded


__all__ = ["UnifiedPostDockingPipeline", "UNIFIED_COMPAT_COLUMNS"]
