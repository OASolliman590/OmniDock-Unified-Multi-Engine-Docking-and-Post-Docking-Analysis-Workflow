from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
import tempfile
import threading
from datetime import datetime, timezone
from functools import wraps
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .ids import build_run_id, iso_timestamp
from .models import CheckpointMetadataRecord, DockForgeFeatureFlags
from .locking import file_lock

WORKFLOW_DIRNAME = ".workflow"
STATE_FILENAME = "state.json"
STATE_BACKUP_SUFFIX = ".bak"
_STATE_IO_LOCK = threading.RLock()
META_DIRNAME = ".meta"


def _state_transaction(fn):
    """Serialize the entire read/update/write transaction, including processes."""
    @wraps(fn)
    def locked(root, *args, **kwargs):
        with file_lock(Path(root).expanduser().resolve() / WORKFLOW_DIRNAME / "state.lock"):
            return fn(root, *args, **kwargs)
    return locked


def _now_iso() -> str:
    return iso_timestamp()


def _default_feature_flags() -> Dict[str, bool]:
    return DockForgeFeatureFlags().to_dict()


def workflow_dir(root: Path) -> Path:
    return Path(root) / WORKFLOW_DIRNAME


def state_path(root: Path) -> Path:
    return workflow_dir(root) / STATE_FILENAME


def _backup_state_path(path: Path) -> Path:
    return path.with_name(f"{path.name}{STATE_BACKUP_SUFFIX}")


def _atomic_dump_json(path: Path, payload: Dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_path = tempfile.mkstemp(prefix=f"{path.name}.", suffix=".tmp", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_path, path)
    finally:
        if os.path.exists(temp_path):
            os.remove(temp_path)


def meta_dir(root: Path) -> Path:
    return Path(root).expanduser().resolve() / META_DIRNAME


def _sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    hasher = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            hasher.update(chunk)
    return hasher.hexdigest()


def _run_version_command(command: List[str], timeout_seconds: int = 8) -> str:
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout_seconds,
        )
    except Exception:
        return ""
    output = (completed.stdout or completed.stderr or "").strip()
    if not output:
        return ""
    return output.splitlines()[0].strip()


def _dump_yaml_or_json(path: Path, payload: Dict[str, object]) -> None:
    try:
        import yaml  # type: ignore

        rendered = yaml.safe_dump(payload, sort_keys=False)
    except Exception:
        rendered = json.dumps(payload, indent=2)
    path.write_text(rendered, encoding="utf-8")


def _default_state(root: Path) -> Dict[str, object]:
    timestamp = _now_iso()
    return {
        "version": 2,
        "project_root": str(root),
        "created_at": timestamp,
        "updated_at": timestamp,
        "current_context": {
            "project_root": str(root),
            "analysis_output_dir": "",
            "favorite_engine": "",
            "top_pose_selection_policy": "",
            "top_pose_global_aggregation": "",
            "selected_engines": [],
            "enabled_panels": [],
            "pair_mode": "",
            "active_target": "",
            "layout_profile": "",
            "raw_proteins_dir": "",
            "raw_ligands_dir": "",
            "raw_ligands_sdf_dir": "",
            "prepared_proteins_dir": "",
            "prepared_ligands_dir": "",
            "docking_root": "",
            "post_docking_root": "",
        },
        "artifacts": {},
        "steps": {},
        "background_tasks": {},
        "feature_flags": _default_feature_flags(),
        "checkpoint_metadata": {
            "latest_checkpoint_id": "",
            "records": [],
        },
    }


@_state_transaction
def ensure_state(root: Path) -> Dict[str, object]:
    root = Path(root).expanduser().resolve()
    with _STATE_IO_LOCK:
        path = state_path(root)
        if path.exists():
            return load_state(root)
        workflow_dir(root).mkdir(parents=True, exist_ok=True)
        payload = _default_state(root)
        save_state(root, payload)
        return payload


@_state_transaction
def load_state(root: Path) -> Dict[str, object]:
    root = Path(root).expanduser().resolve()
    with _STATE_IO_LOCK:
        path = state_path(root)
        if not path.exists():
            return ensure_state(root)
        try:
            with open(path, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
        except json.JSONDecodeError:
            backup_path = _backup_state_path(path)
            if backup_path.exists():
                with open(backup_path, "r", encoding="utf-8") as handle:
                    payload = json.load(handle)
                _atomic_dump_json(path, payload)
            else:
                payload = _default_state(root)
                _atomic_dump_json(path, payload)
    # Backward-compatible migration for older state files.
    payload.setdefault("background_tasks", {})
    payload.setdefault("steps", {})
    payload.setdefault("artifacts", {})
    payload.setdefault("current_context", {})
    payload.setdefault("feature_flags", _default_feature_flags())
    checkpoint_metadata = payload.setdefault("checkpoint_metadata", {})
    if not isinstance(checkpoint_metadata, dict):
        checkpoint_metadata = {}
    checkpoint_metadata.setdefault("latest_checkpoint_id", "")
    checkpoint_metadata.setdefault("records", [])
    payload["checkpoint_metadata"] = checkpoint_metadata
    return payload


@_state_transaction
def save_state(root: Path, payload: Dict[str, object]) -> Path:
    root = Path(root).expanduser().resolve()
    with _STATE_IO_LOCK:
        workflow_dir(root).mkdir(parents=True, exist_ok=True)
        payload["updated_at"] = _now_iso()
        payload.setdefault("background_tasks", {})
        payload.setdefault("steps", {})
        payload.setdefault("feature_flags", _default_feature_flags())
        checkpoint_metadata = payload.setdefault("checkpoint_metadata", {})
        if isinstance(checkpoint_metadata, dict):
            checkpoint_metadata.setdefault("latest_checkpoint_id", "")
            checkpoint_metadata.setdefault("records", [])
        path = state_path(root)
        backup_path = _backup_state_path(path)
        if path.exists():
            try:
                shutil.copy2(path, backup_path)
            except OSError:
                pass
        _atomic_dump_json(path, payload)
        return path


def get_feature_flags(root: Path) -> Dict[str, bool]:
    payload = ensure_state(root)
    flags = payload.get("feature_flags", {})
    if not isinstance(flags, dict):
        flags = {}
    merged: Dict[str, bool] = _default_feature_flags()
    for key, default_value in merged.items():
        value = flags.get(key, default_value)
        merged[key] = bool(value)
    return merged


@_state_transaction
def update_feature_flags(root: Path, **kwargs: object) -> Dict[str, bool]:
    payload = ensure_state(root)
    existing = payload.setdefault("feature_flags", _default_feature_flags())
    if not isinstance(existing, dict):
        existing = _default_feature_flags()
    supported = set(_default_feature_flags().keys())
    for key, value in kwargs.items():
        if key in supported and value is not None:
            existing[key] = bool(value)
    payload["feature_flags"] = existing
    save_state(root, payload)
    return get_feature_flags(root)


@_state_transaction
def update_context(root: Path, **kwargs: object) -> Dict[str, object]:
    payload = ensure_state(root)
    payload.setdefault("current_context", {}).update(
        {key: value for key, value in kwargs.items() if value not in (None, "")}
    )
    save_state(root, payload)
    return payload


@_state_transaction
def update_artifacts(root: Path, **kwargs: object) -> Dict[str, object]:
    payload = ensure_state(root)
    payload.setdefault("artifacts", {}).update(
        {key: value for key, value in kwargs.items() if value not in (None, "")}
    )
    save_state(root, payload)
    return payload


@_state_transaction
def record_step(
    root: Path,
    target: str,
    status: str,
    inputs: Optional[Dict[str, object]] = None,
    outputs: Optional[Dict[str, object]] = None,
    notes: Optional[List[str]] = None,
) -> Dict[str, object]:
    payload = ensure_state(root)
    payload.setdefault("current_context", {})["active_target"] = target
    step_record = payload.setdefault("steps", {}).get(target, {})
    started_at = step_record.get("started_at")
    if not started_at:
        started_at = _now_iso()
    payload["steps"][target] = {
        "status": status,
        "started_at": started_at,
        "finished_at": _now_iso(),
        "inputs": inputs or {},
        "outputs": outputs or {},
        "notes": notes or [],
    }
    save_state(root, payload)
    return payload


def list_steps(root: Path) -> Dict[str, Dict[str, object]]:
    payload = ensure_state(root)
    steps = payload.get("steps", {})
    if isinstance(steps, dict):
        return {str(key): value for key, value in steps.items() if isinstance(value, dict)}
    return {}


@_state_transaction
def update_background_task(
    root: Path,
    task_id: str,
    *,
    label: str,
    target: str,
    status: str,
    progress: int = 0,
    note: str = "",
    result_status: str = "",
    error: str = "",
) -> Dict[str, object]:
    payload = ensure_state(root)
    tasks = payload.setdefault("background_tasks", {})
    existing = tasks.get(task_id, {}) if isinstance(tasks, dict) else {}
    timestamp = _now_iso()
    if not isinstance(existing, dict):
        existing = {}
    tasks[task_id] = {
        "task_id": task_id,
        "label": label,
        "target": target,
        "status": status,
        "progress": max(0, min(int(progress), 100)),
        "note": note,
        "result_status": result_status,
        "error": error,
        "created_at": existing.get("created_at", timestamp),
        "updated_at": timestamp,
        "root": str(Path(root).expanduser().resolve()),
    }
    save_state(root, payload)
    return tasks[task_id]


def list_background_tasks(root: Path, *, include_finished: bool = True) -> List[Dict[str, object]]:
    payload = ensure_state(root)
    tasks = payload.get("background_tasks", {})
    if not isinstance(tasks, dict):
        return []
    rows = [value for value in tasks.values() if isinstance(value, dict)]
    if not include_finished:
        rows = [row for row in rows if str(row.get("status", "")).strip() not in {"completed", "failed", "cancelled"}]
    rows.sort(key=lambda row: str(row.get("updated_at", "")), reverse=True)
    return rows


@_state_transaction
def record_checkpoint_metadata(
    root: Path,
    *,
    source_project_dir: str,
    target_project_dir: str,
    layout_profile: str,
    marker_file: str = "",
    note: str = "",
    metadata: Optional[Dict[str, Any]] = None,
    checkpoint_id: Optional[str] = None,
) -> Dict[str, object]:
    """Persist checkpoint lineage metadata for Checkpoint & Revise operations."""

    payload = ensure_state(root)
    checkpoint_store = payload.setdefault("checkpoint_metadata", {})
    if not isinstance(checkpoint_store, dict):
        checkpoint_store = {}
    records = checkpoint_store.setdefault("records", [])
    if not isinstance(records, list):
        records = []

    checkpoint_token = str(checkpoint_id or build_run_id("checkpoint"))
    record = CheckpointMetadataRecord(
        checkpoint_id=checkpoint_token,
        source_project_dir=str(source_project_dir),
        target_project_dir=str(target_project_dir),
        created_at=_now_iso(),
        layout_profile=str(layout_profile or ""),
        marker_file=str(marker_file or ""),
        note=str(note or ""),
        metadata=dict(metadata or {}),
    ).to_dict()
    records.append(record)
    checkpoint_store["records"] = records
    checkpoint_store["latest_checkpoint_id"] = checkpoint_token
    payload["checkpoint_metadata"] = checkpoint_store
    save_state(root, payload)
    return record


def list_checkpoint_metadata(root: Path) -> List[Dict[str, object]]:
    payload = ensure_state(root)
    checkpoint_store = payload.get("checkpoint_metadata", {})
    if not isinstance(checkpoint_store, dict):
        return []
    records = checkpoint_store.get("records", [])
    if not isinstance(records, list):
        return []
    # Stable sort then uses reverse insertion order as the tie-breaker, so
    # multiple checkpoints recorded within one timestamp tick stay newest-first.
    cleaned = [row for row in reversed(records) if isinstance(row, dict)]
    cleaned.sort(key=lambda row: str(row.get("created_at", "")), reverse=True)
    return cleaned


def summarize_state(root: Path) -> List[str]:
    payload = ensure_state(root)
    flags = get_feature_flags(root)
    enabled_flags = sorted([name for name, enabled in flags.items() if enabled])
    checkpoint_records = list_checkpoint_metadata(root)
    lines = [
        f"Workflow root: {payload.get('project_root', root)}",
        f"Updated: {payload.get('updated_at', '')}",
    ]
    if enabled_flags:
        lines.append(f"Feature flags enabled: {', '.join(enabled_flags)}")
    if checkpoint_records:
        latest_checkpoint = checkpoint_records[0]
        lines.append(
            "Latest checkpoint: "
            f"{latest_checkpoint.get('checkpoint_id', '')} -> {latest_checkpoint.get('target_project_dir', '')}"
        )
    current = payload.get("current_context", {})
    if current.get("favorite_engine"):
        lines.append(f"Favorite engine: {current['favorite_engine']}")
    if current.get("top_pose_selection_policy"):
        lines.append(f"Top-pose policy: {current['top_pose_selection_policy']}")
    if current.get("top_pose_global_aggregation"):
        lines.append(f"Top-pose aggregation: {current['top_pose_global_aggregation']}")
    if current.get("selected_engines"):
        lines.append(f"Engines: {', '.join(current['selected_engines'])}")
    if current.get("analysis_output_dir"):
        lines.append(f"Analysis output: {current['analysis_output_dir']}")
    active_tasks = list_background_tasks(root, include_finished=False)
    if active_tasks:
        lines.append(f"Active background tasks: {len(active_tasks)}")
    steps = payload.get("steps", {})
    if not steps:
        lines.append("No workflow steps recorded yet.")
        return lines
    lines.append("Recorded steps:")
    for target, record in sorted(steps.items()):
        lines.append(f"  - {target}: {record.get('status', 'unknown')}")
    return lines


def write_meta_config_freeze(
    root: Path,
    *,
    config: Dict[str, object],
    source: str = "",
    note: str = "",
) -> Path:
    """
    Write a frozen config snapshot under `.meta/config.yaml`.
    """
    root = Path(root).expanduser().resolve()
    target_dir = meta_dir(root)
    target_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "generated_at": _now_iso(),
        "project_root": str(root),
        "source": str(source or ""),
        "note": str(note or ""),
        "config": dict(config or {}),
    }
    output_file = target_dir / "config.yaml"
    _dump_yaml_or_json(output_file, payload)
    return output_file


def write_meta_run_manifest(
    root: Path,
    *,
    input_paths: Optional[List[Path]] = None,
    run_config: Optional[Dict[str, object]] = None,
    tool_versions: Optional[Dict[str, str]] = None,
) -> Path:
    """
    Write `.meta/run_manifest.json` with checksums, versions, and timestamps.
    """
    root = Path(root).expanduser().resolve()
    target_dir = meta_dir(root)
    target_dir.mkdir(parents=True, exist_ok=True)

    checksums: List[Dict[str, object]] = []
    missing_inputs: List[str] = []
    normalized_inputs = sorted(
        {
            str(Path(candidate).expanduser().resolve())
            for candidate in (input_paths or [])
        }
    )
    for raw_path in normalized_inputs:
        path = Path(raw_path)
        if not path.exists() or not path.is_file():
            missing_inputs.append(str(path))
            continue
        try:
            relative = path.relative_to(root)
            relative_path = str(relative)
        except ValueError:
            relative_path = str(path)
        checksums.append(
            {
                "path": relative_path,
                "sha256": _sha256_file(path),
                "size_bytes": int(path.stat().st_size),
                "modified_at": datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat(),
            }
        )
    checksums.sort(key=lambda row: str(row.get("path", "")))

    detected_versions = dict(tool_versions or {})
    if not detected_versions:
        detected_versions.update(
            {
                "python": _run_version_command(["python", "--version"]) or platform.python_version(),
                "obabel": _run_version_command(["obabel", "-V"]),
                "vina": _run_version_command(["vina", "--version"]),
                "smina": _run_version_command(["smina", "--version"]),
                "gnina": _run_version_command(["gnina", "--version"]),
                "autodock4": _run_version_command(["autodock4", "--version"]),
                "autogrid4": _run_version_command(["autogrid4", "--version"]),
            }
        )
    detected_versions = {key: value for key, value in detected_versions.items() if str(value).strip()}

    payload = {
        "generated_at": _now_iso(),
        "project_root": str(root),
        "state_file": str(state_path(root)),
        "input_paths": [row.get("path", "") for row in checksums],
        "missing_input_paths": missing_inputs,
        "input_count": int(len(checksums)),
        "missing_input_count": int(len(missing_inputs)),
        "input_checksum_digest": hashlib.sha256(
            json.dumps(checksums, sort_keys=True).encode("utf-8")
        ).hexdigest()
        if checksums
        else "",
        "input_checksums": checksums,
        "tool_versions": detected_versions,
        "run_config": dict(run_config or {}),
    }
    output_file = target_dir / "run_manifest.json"
    output_file.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return output_file


def write_meta_env_lock(
    root: Path,
    *,
    preferred: Optional[str] = None,
) -> Path:
    """
    Write `.meta/env.lock` snapshot with graceful fallback.
    Preference order: explicit `preferred`, then conda export, then pip freeze.
    """
    root = Path(root).expanduser().resolve()
    target_dir = meta_dir(root)
    target_dir.mkdir(parents=True, exist_ok=True)

    preferred_token = str(preferred or "").strip().lower()
    commands: List[Tuple[str, List[str]]] = []
    if preferred_token == "conda":
        commands.append(("conda", ["conda", "env", "export", "--no-builds"]))
    elif preferred_token == "pip":
        commands.append(("pip", ["python", "-m", "pip", "freeze"]))
    commands.extend(
        [
            ("conda", ["conda", "env", "export", "--no-builds"]),
            ("pip", ["python", "-m", "pip", "freeze"]),
        ]
    )

    snapshot_source = "unavailable"
    snapshot_text = ""
    errors: List[str] = []
    seen = set()
    for source, command in commands:
        key = tuple(command)
        if key in seen:
            continue
        seen.add(key)
        try:
            completed = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=False,
                timeout=20,
            )
        except Exception as exc:
            errors.append(f"{source}: {exc}")
            continue
        if completed.returncode != 0:
            stderr = (completed.stderr or "").strip()
            errors.append(f"{source}: {stderr or f'returncode={completed.returncode}'}")
            continue
        output = (completed.stdout or "").strip()
        if output:
            snapshot_source = source
            snapshot_text = output
            break
        errors.append(f"{source}: empty output")

    if not snapshot_text:
        snapshot_text = "\n".join(
            [
                "# Environment snapshot unavailable",
                "# Neither conda env export nor pip freeze returned usable output.",
            ]
            + [f"# {message}" for message in errors]
        )

    output_file = target_dir / "env.lock"
    header = "\n".join(
        [
            f"# generated_at: {_now_iso()}",
            f"# project_root: {root}",
            f"# source: {snapshot_source}",
            "",
        ]
    )
    output_file.write_text(header + snapshot_text + ("\n" if not snapshot_text.endswith("\n") else ""), encoding="utf-8")
    return output_file
