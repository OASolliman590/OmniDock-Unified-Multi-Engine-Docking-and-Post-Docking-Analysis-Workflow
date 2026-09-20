from __future__ import annotations

import json
import re
import shlex
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Callable, Dict, List, Optional, Tuple

from .engine_registry import build_runner
from .hpc_profiles import merge_condor_with_profile, merge_slurm_with_profile
from .models import PairlistRow
from .project_layout import deployment_root, detect_layout_profile, ensure_project_layout
from .runners import job_contract


DEFAULT_SLURM = {
    "time": "08:00:00",
    "mem": "16G",
    "cpus_per_task": 4,
    "gpus": 0,
    "partition": "",
    "account": "",
    "job_name_prefix": "dock",
    "submit_mode": "slurm_array",
    "array_parallelism": 0,
    "extra_args": [],
}


def _sanitize_token(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "-", str(value).strip())
    return cleaned.strip("-") or "job"


def _normalize_extra_args(raw_value: object) -> List[str]:
    if raw_value in (None, ""):
        return []
    if isinstance(raw_value, list):
        return [str(item).strip() for item in raw_value if str(item).strip()]
    return [item.strip() for item in str(raw_value).split(";") if item.strip()]


def _is_relative_to(path: Path, base: Path) -> bool:
    try:
        path.relative_to(base)
        return True
    except ValueError:
        return False


def _map_project_path(raw_value: object, local_root: Path, execution_root: Path) -> str:
    text = str(raw_value)
    local_prefix = str(local_root)
    execution_prefix = str(execution_root)
    if local_prefix and local_prefix in text:
        mapped = text.replace(local_prefix, execution_prefix)
        if isinstance(execution_root, PurePosixPath) and "\\" in local_prefix:
            mapped = mapped.replace("\\", "/")
        return mapped
    candidate = Path(text)
    if not candidate.is_absolute():
        return text
    if not _is_relative_to(candidate, local_root):
        return text
    relative = candidate.relative_to(local_root)
    return str(execution_root / relative)


def _map_command(command: List[str], local_root: Path, execution_root: Path) -> List[str]:
    return [_map_project_path(token, local_root, execution_root) for token in command]


def _execution_root(local_root, remote_root):
    if not remote_root:
        return local_root
    remote = PurePosixPath(str(remote_root).replace("\\", "/"))
    if not remote.is_absolute():
        raise ValueError("Remote project root must be an absolute Linux/POSIX path")
    return remote


def _mapped_path(value, local_root, execution_root):
    mapped = _map_project_path(value, local_root, execution_root)
    return PurePosixPath(mapped) if type(execution_root) is PurePosixPath else Path(mapped)


def _map_job(job, local_root, execution_root):
    payload = job.to_dict()
    for key in ("pose_file", "log_file", "completion_file"):
        payload[key] = _map_project_path(payload[key], local_root, execution_root)
    payload["command"] = _map_command(payload["command"], local_root, execution_root)
    payload["input_files"] = [{**item, "path": _map_project_path(item["path"], local_root, execution_root)} for item in payload["input_files"]]
    payload["executables"] = [_map_project_path(value, local_root, execution_root) for value in payload["executables"]]
    return payload


def _write_job_executor(scripts_dir):
    target = scripts_dir / "execute_job.py"
    target.write_text(Path(job_contract.__file__).read_text(encoding="utf-8"), encoding="utf-8")


def _apply_allocated_resources(engine, runtime, settings, cpu_key):
    allocated = int(settings[cpu_key])
    requested = int(runtime.get("cpu") or allocated)
    if allocated < 1 or requested < 1 or requested > allocated:
        raise ValueError(f"{engine}: requested CPUs ({requested}) must be between 1 and scheduler allocation ({allocated})")
    runtime["cpu"] = requested
    if engine == "gnina":
        runtime["use_gpu"] = int(settings.get("gpus") or 0) > 0


def _slurm_settings_for_engine(
    engine: str,
    runtime: Dict[str, object],
    slurm_options: Dict[str, object],
    profile: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    merged = dict(DEFAULT_SLURM)
    merged.update(merge_slurm_with_profile(engine, runtime, slurm_options, profile or {}))
    merged["extra_args"] = _normalize_extra_args(merged.get("extra_args"))
    if not merged.get("cpus_per_task"):
        merged["cpus_per_task"] = int(runtime.get("cpu") or 4)
    if engine == "gnina" and runtime.get("use_gpu") is not False:
        merged["gpus"] = int(slurm_options["gpus"] if slurm_options.get("gpus") is not None else (merged.get("gpus") or 1))
    else:
        merged["gpus"] = int(slurm_options.get("gpus") or 0)
    return merged


def _validate_slurm_settings(
    engine: str,
    settings: Dict[str, object],
    profile: Optional[Dict[str, object]] = None,
) -> None:
    partition = str(settings.get("partition") or "").strip()
    if partition:
        return
    profile_name = str((profile or {}).get("name", "") or "").strip() or "custom"
    raise ValueError(
        "No Slurm partition configured for "
        f"{engine}. Configure profile slurm.partition or engines.{engine}.slurm.partition, "
        "or pass --slurm-partition."
        f" (active profile: {profile_name})"
    )


def _render_slurm_header(
    job_name: str,
    stdout_path: Path,
    stderr_path: Path,
    settings: Dict[str, object],
    array_spec: str = "",
) -> List[str]:
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --output={stdout_path}",
        f"#SBATCH --error={stderr_path}",
        f"#SBATCH --time={settings['time']}",
        f"#SBATCH --mem={settings['mem']}",
        f"#SBATCH --cpus-per-task={int(settings['cpus_per_task'])}",
    ]
    if settings.get("gpus"):
        lines.append(f"#SBATCH --gpus={int(settings['gpus'])}")
    if settings.get("partition"):
        lines.append(f"#SBATCH --partition={settings['partition']}")
    if settings.get("account"):
        lines.append(f"#SBATCH --account={settings['account']}")
    if array_spec:
        lines.append(f"#SBATCH --array={array_spec}")
    for extra in settings.get("extra_args", []):
        lines.append(f"#SBATCH {extra}")
    return lines


def _render_array_slurm_script(
    execution_project_root: Path,
    engine: str,
    job_name: str,
    execution_manifest_path: Path,
    log_dir: Path,
    settings: Dict[str, object],
    job_count: int,
    preamble_lines: Optional[List[str]] = None,
) -> str:
    native_pair_logs = True  # The portable executor captures per-engine logs.
    stdout_path = log_dir / f"{job_name}-%A_%a.out"
    stderr_path = log_dir / f"{job_name}-%A_%a.err"
    array_parallelism = int(settings.get("array_parallelism") or 0)
    array_spec = f"1-{job_count}"
    if array_parallelism > 0:
        array_spec = f"{array_spec}%{array_parallelism}"
    lines = _render_slurm_header(job_name, stdout_path, stderr_path, settings, array_spec=array_spec)
    lines.extend(
        [
            "",
            "set -eo pipefail",
            f'manifest_file={shlex.quote(str(execution_manifest_path))}',
            'line="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$manifest_file")"',
            'if [ -z "$line" ]; then',
            '  echo "No manifest entry found for task ${SLURM_ARRAY_TASK_ID}" >&2',
            "  exit 1",
            "fi",
            "IFS=$'\\t' read -r pair_job_name pair_tag pair_command pair_log_file <<< \"$line\"",
            f"cd {shlex.quote(str(execution_project_root))}",
            'if [ -n "$pair_log_file" ]; then',
            '  mkdir -p "$(dirname "$pair_log_file")"',
            "fi",
            "",
        ]
    )
    for line in preamble_lines or []:
        stripped = str(line).rstrip()
        if stripped:
            if stripped.startswith("if [ -f ~/.bashrc ]") or "source ~/.bashrc" in stripped:
                lines.append("set +u")
            lines.append(stripped)
            if stripped.startswith("if [ -f ~/.bashrc ]") or "source ~/.bashrc" in stripped:
                lines.append("set -u")
    if preamble_lines:
        lines.append("")
    lines.extend(
        [
            f"echo \"[$(date -Is)] Starting {engine} array task ${'{'}SLURM_ARRAY_TASK_ID{'}'}: ${'{'}pair_tag{'}'}\"",
            'eval "$pair_command"' if native_pair_logs else 'eval "$pair_command" >"$pair_log_file" 2>&1',
            f"echo \"[$(date -Is)] Finished {engine} array task ${'{'}SLURM_ARRAY_TASK_ID{'}'}: ${'{'}pair_tag{'}'}\"",
            "",
        ]
    )
    return "\n".join(lines)


def _render_single_job_slurm_script(
    execution_project_root: Path,
    engine: str,
    job_name: str,
    execution_manifest_path: Path,
    log_dir: Path,
    settings: Dict[str, object],
    preamble_lines: Optional[List[str]] = None,
) -> str:
    native_pair_logs = True  # The portable executor captures per-engine logs.
    stdout_path = log_dir / f"{job_name}-%j.out"
    stderr_path = log_dir / f"{job_name}-%j.err"
    lines = _render_slurm_header(job_name, stdout_path, stderr_path, settings)
    lines.extend(
        [
            "",
            "set -uo pipefail",
            f'manifest_file={shlex.quote(str(execution_manifest_path))}',
            f"cd {shlex.quote(str(execution_project_root))}",
            "total=0",
            "failed=0",
            "",
        ]
    )
    for line in preamble_lines or []:
        stripped = str(line).rstrip()
        if stripped:
            if stripped.startswith("if [ -f ~/.bashrc ]") or "source ~/.bashrc" in stripped:
                lines.append("set +u")
            lines.append(stripped)
            if stripped.startswith("if [ -f ~/.bashrc ]") or "source ~/.bashrc" in stripped:
                lines.append("set -u")
    if preamble_lines:
        lines.append("")
    lines.extend(
        [
            "while IFS=$'\\t' read -r pair_job_name pair_tag pair_command pair_log_file; do",
            "  [ -n \"$pair_tag\" ] || continue",
            "  total=$((total + 1))",
            f"  echo \"[$(date -Is)] Starting {engine} job ${'{'}total{'}'}: ${'{'}pair_tag{'}'}\"",
            "  mkdir -p \"$(dirname \"$pair_log_file\")\"",
            "  if eval \"$pair_command\"; then" if native_pair_logs else "  if eval \"$pair_command\" >\"$pair_log_file\" 2>&1; then",
            f"    echo \"[$(date -Is)] Finished {engine} job ${'{'}total{'}'}: ${'{'}pair_tag{'}'}\"",
            "  else",
            "    failed=$((failed + 1))",
            f"    echo \"[$(date -Is)] Failed {engine} job ${'{'}total{'}'}: ${'{'}pair_tag{'}'}\" >&2",
            "  fi",
            "done < \"$manifest_file\"",
            "",
            f"echo \"[$(date -Is)] {engine} manifest complete: total=${'{'}total{'}'} failed=${'{'}failed{'}'}\"",
            "if [ \"$failed\" -gt 0 ]; then",
            "  exit 1",
            "fi",
            "",
        ]
    )
    return "\n".join(lines)


DEFAULT_CONDOR = {
    "cpus": 8,
    "memory": "32GB",
    "disk": "20GB",
    "gpus": 0,
    "requirements": "",
    "job_name_prefix": "dock",
    "submit_mode": "condor_array",
    "parallelism": 0,
}


def _condor_settings_for_engine(
    engine: str,
    runtime: Dict[str, object],
    condor_options: Dict[str, object],
    profile: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    merged = dict(DEFAULT_CONDOR)
    merged.update(merge_condor_with_profile(engine, runtime, condor_options, profile or {}))
    if not merged.get("cpus"):
        merged["cpus"] = int(runtime.get("cpu") or 8)
    if engine == "gnina" and runtime.get("use_gpu") is not False:
        merged["gpus"] = int(condor_options["gpus"] if condor_options.get("gpus") is not None else (merged.get("gpus") or 1))
    else:
        merged["gpus"] = int(condor_options.get("gpus") or 0)
    return merged


def _validate_condor_settings(
    engine: str,
    settings: Dict[str, object],
    profile: Optional[Dict[str, object]] = None,
) -> None:
    submit_mode = str(settings.get("submit_mode") or "").strip()
    if submit_mode not in {"condor_array", "single_job"}:
        profile_name = str((profile or {}).get("name", "") or "").strip() or "custom"
        raise ValueError(
            f"Unsupported condor submit_mode for {engine}: {submit_mode!r}. "
            "Use 'condor_array' or 'single_job'."
            f" (active profile: {profile_name})"
        )


def _render_condor_submit_file(
    executable_path: Path,
    log_dir: Path,
    settings: Dict[str, object],
    job_count: int = 1,
    array_mode: bool = True,
) -> str:
    parallelism = int(settings.get("parallelism") or 0)
    if array_mode:
        arguments_line = "arguments               = $(Process)"
        output_line = f"output                  = {log_dir}/job-$(Process).out"
        error_line = f"error                   = {log_dir}/job-$(Process).err"
        queue_lines = []
        if parallelism > 0:
            queue_lines.append(f"max_materialize         = {parallelism}")
        queue_lines.append(f"queue {job_count}")
    else:
        arguments_line = "arguments               ="
        output_line = f"output                  = {log_dir}/job.out"
        error_line = f"error                   = {log_dir}/job.err"
        queue_lines = ["queue 1"]

    lines = [
        f"executable              = {executable_path}",
        arguments_line,
        f"log                     = {log_dir}/job.log",
        output_line,
        error_line,
        "",
        f"request_cpus            = {int(settings['cpus'])}",
        f"request_memory          = {settings['memory']}",
        f"request_disk            = {settings['disk']}",
    ]
    if settings.get("gpus") and int(settings["gpus"]) > 0:
        lines.append(f"request_gpus            = {int(settings['gpus'])}")
    if settings.get("requirements"):
        lines.append(f"requirements            = {settings['requirements']}")
    lines.extend([
        "",
        "getenv                  = True",
        "should_transfer_files   = NO",
        "transfer_executable     = FALSE",
        "",
        *queue_lines,
        "",
    ])
    return "\n".join(lines)


def _render_array_condor_script(
    execution_project_root: Path,
    engine: str,
    execution_manifest_path: Path,
    preamble_lines: Optional[List[str]] = None,
) -> str:
    native_pair_logs = True  # The portable executor captures per-engine logs.
    lines = [
        "#!/bin/bash",
        "set -eo pipefail",
        "task_id=$(($1 + 1))",
        f'manifest_file={shlex.quote(str(execution_manifest_path))}',
        'line="$(sed -n "${task_id}p" "$manifest_file")"',
        'if [ -z "$line" ]; then',
        '  echo "No manifest entry for process $1 (task ${task_id})" >&2',
        "  exit 1",
        "fi",
        "IFS=$'\\t' read -r pair_job_name pair_tag pair_command pair_log_file <<< \"$line\"",
        f"cd {shlex.quote(str(execution_project_root))}",
        'if [ -n "$pair_log_file" ]; then',
        '  mkdir -p "$(dirname "$pair_log_file")"',
        "fi",
        "",
    ]
    for line in preamble_lines or []:
        stripped = str(line).rstrip()
        if stripped:
            if stripped.startswith("if [ -f ~/.bashrc ]") or "source ~/.bashrc" in stripped:
                lines.append("set +u")
            lines.append(stripped)
            if stripped.startswith("if [ -f ~/.bashrc ]") or "source ~/.bashrc" in stripped:
                lines.append("set -u")
    if preamble_lines:
        lines.append("")
    lines.extend([
        f"echo \"[$(date -Is)] Starting {engine} process $1: ${'{'}pair_tag{'}'}\"",
        'eval "$pair_command"' if native_pair_logs else 'eval "$pair_command" >"$pair_log_file" 2>&1',
        f"echo \"[$(date -Is)] Finished {engine} process $1: ${'{'}pair_tag{'}'}\"",
        "",
    ])
    return "\n".join(lines)


def _render_single_job_condor_script(
    execution_project_root: Path,
    engine: str,
    execution_manifest_path: Path,
    preamble_lines: Optional[List[str]] = None,
) -> str:
    native_pair_logs = True  # The portable executor captures per-engine logs.
    lines = [
        "#!/bin/bash",
        "set -uo pipefail",
        f'manifest_file={shlex.quote(str(execution_manifest_path))}',
        f"cd {shlex.quote(str(execution_project_root))}",
        "total=0",
        "failed=0",
        "",
    ]
    for line in preamble_lines or []:
        stripped = str(line).rstrip()
        if stripped:
            if stripped.startswith("if [ -f ~/.bashrc ]") or "source ~/.bashrc" in stripped:
                lines.append("set +u")
            lines.append(stripped)
            if stripped.startswith("if [ -f ~/.bashrc ]") or "source ~/.bashrc" in stripped:
                lines.append("set -u")
    if preamble_lines:
        lines.append("")
    lines.extend([
        "while IFS=$'\\t' read -r pair_job_name pair_tag pair_command pair_log_file; do",
        "  [ -n \"$pair_tag\" ] || continue",
        "  total=$((total + 1))",
        f"  echo \"[$(date -Is)] Starting {engine} job ${'{'}total{'}'}: ${'{'}pair_tag{'}'}\"",
        "  mkdir -p \"$(dirname \"$pair_log_file\")\"",
        "  if eval \"$pair_command\"; then" if native_pair_logs else "  if eval \"$pair_command\" >\"$pair_log_file\" 2>&1; then",
        f"    echo \"[$(date -Is)] Finished {engine} job ${'{'}total{'}'}: ${'{'}pair_tag{'}'}\"",
        "  else",
        "    failed=$((failed + 1))",
        f"    echo \"[$(date -Is)] Failed {engine} job ${'{'}total{'}'}: ${'{'}pair_tag{'}'}\" >&2",
        "  fi",
        "done < \"$manifest_file\"",
        "",
        f"echo \"[$(date -Is)] {engine} manifest complete: total=${'{'}total{'}'} failed=${'{'}failed{'}'}\"",
        "if [ \"$failed\" -gt 0 ]; then",
        "  exit 1",
        "fi",
        "",
    ])
    return "\n".join(lines)


@dataclass(frozen=True)
class _DeploymentContext:
    root: Path
    execution_root: Path | PurePosixPath
    layout_profile: str
    stage_root: Path
    execution_stage_root: Path | PurePosixPath
    round_id: str
    docking_mode: str
    hpc_profile: Dict[str, str]
    remote_project_root: str
    skip_completed: bool
    pair_source_file: str
    rerun_manifest_file: str
    pairlist_rows: List[PairlistRow]


@dataclass(frozen=True)
class _EngineWorkspace:
    engine: str
    runtime: Dict[str, object]
    mapped_runtime: Dict[str, object]
    jobs: List[object]
    engine_root: Path
    execution_engine_root: Path | PurePosixPath
    scripts_dir: Path
    manifest_dir: Path
    logs_dir: Path
    execution_logs_dir: Path | PurePosixPath


def _open_deployment_context(
    project_root: Path,
    remote_project_root: str,
    round_id: str,
    docking_mode: str,
    hpc_profile: Optional[Dict[str, object]],
    skip_completed: bool,
    pair_source_file: str,
    rerun_manifest_file: str,
    pairlist_rows: List[PairlistRow],
) -> _DeploymentContext:
    root = Path(project_root).expanduser().resolve()
    execution_root = _execution_root(root, remote_project_root)
    layout_profile = detect_layout_profile(root)
    ensure_project_layout(root, layout_profile)
    stage_root = deployment_root(root, layout_profile) / _sanitize_token(round_id) / _sanitize_token(docking_mode)
    execution_stage_root = _mapped_path(stage_root, root, execution_root)
    stage_root.mkdir(parents=True, exist_ok=True)
    return _DeploymentContext(
        root=root,
        execution_root=execution_root,
        layout_profile=layout_profile,
        stage_root=stage_root,
        execution_stage_root=execution_stage_root,
        round_id=round_id,
        docking_mode=docking_mode,
        hpc_profile={
            "name": str((hpc_profile or {}).get("name", "")),
            "source": str((hpc_profile or {}).get("source", "")),
        },
        remote_project_root=remote_project_root or "",
        skip_completed=skip_completed,
        pair_source_file=pair_source_file,
        rerun_manifest_file=rerun_manifest_file,
        pairlist_rows=pairlist_rows,
    )


def _base_deployment_manifest(
    ctx: _DeploymentContext,
    scheduler: str,
    engines: List[str],
) -> Dict[str, object]:
    return {
        "generation_project_root": str(ctx.root),
        "project_root": str(ctx.execution_root),
        "execution_project_root": str(ctx.execution_root),
        "layout_profile": ctx.layout_profile,
        "round_id": ctx.round_id,
        "docking_mode": ctx.docking_mode,
        "engines": engines,
        "pair_count": len(ctx.pairlist_rows),
        "pair_source_file": _map_project_path(ctx.pair_source_file, ctx.root, ctx.execution_root) if ctx.pair_source_file else "",
        "rerun_manifest_file": _map_project_path(ctx.rerun_manifest_file, ctx.root, ctx.execution_root) if ctx.rerun_manifest_file else "",
        "generation_stage_root": str(ctx.stage_root),
        "stage_root": str(ctx.execution_stage_root),
        "scheduler": scheduler,
        "hpc_profile": ctx.hpc_profile,
        "engine_jobs": {},
    }


def _prepare_engine_workspace(
    ctx: _DeploymentContext,
    engine: str,
    runtime_by_engine: Dict[str, Dict[str, object]],
    settings_for_engine: Callable[[str, Dict[str, object]], Dict[str, object]],
    cpu_key: str,
    log_dir_name: str,
) -> Tuple[_EngineWorkspace, Dict[str, object]]:
    runtime = dict(runtime_by_engine.get(engine, {}))
    runtime["round_id"] = ctx.round_id
    runtime["docking_mode"] = ctx.docking_mode
    if ctx.rerun_manifest_file:
        runtime["rerun_manifest_file"] = _map_project_path(ctx.rerun_manifest_file, ctx.root, ctx.execution_root)
    engine_settings = settings_for_engine(engine, runtime)
    _apply_allocated_resources(engine, runtime, engine_settings, cpu_key)
    runner = build_runner(engine, ctx.root, runtime)
    jobs = runner.plan_jobs(ctx.pairlist_rows, skip_completed=ctx.skip_completed and not ctx.remote_project_root)
    for job in jobs:
        job.skip_completed = ctx.skip_completed
    engine_root = ctx.stage_root / engine
    execution_engine_root = ctx.execution_stage_root / engine
    scripts_dir = engine_root / "job_scripts"
    manifest_dir = engine_root / "manifest"
    logs_dir = engine_root / log_dir_name
    execution_logs_dir = execution_engine_root / log_dir_name
    scripts_dir.mkdir(parents=True, exist_ok=True)
    _write_job_executor(scripts_dir)
    manifest_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    mapped_runtime = {
        key: _map_project_path(value, ctx.root, ctx.execution_root) if isinstance(value, str) else value
        for key, value in runtime.items()
    }
    return (
        _EngineWorkspace(
            engine=engine,
            runtime=runtime,
            mapped_runtime=mapped_runtime,
            jobs=jobs,
            engine_root=engine_root,
            execution_engine_root=execution_engine_root,
            scripts_dir=scripts_dir,
            manifest_dir=manifest_dir,
            logs_dir=logs_dir,
            execution_logs_dir=execution_logs_dir,
        ),
        engine_settings,
    )


def _engine_payload(
    ctx: _DeploymentContext,
    workspace: _EngineWorkspace,
    settings_key: str,
    engine_settings: Dict[str, object],
) -> Dict[str, object]:
    return {
        "engine": workspace.engine,
        "generation_project_root": str(ctx.root),
        "project_root": str(ctx.execution_root),
        "execution_project_root": str(ctx.execution_root),
        "round_id": ctx.round_id,
        "docking_mode": ctx.docking_mode,
        "runtime": workspace.mapped_runtime,
        settings_key: engine_settings,
        "hpc_profile": ctx.hpc_profile,
        "jobs": [],
    }


def _planned_job_payloads(
    ctx: _DeploymentContext,
    workspace: _EngineWorkspace,
    engine_settings: Dict[str, object],
    script_field: str,
    execution_script_path: Path | PurePosixPath,
) -> Tuple[List[str], int, List[Dict[str, object]]]:
    planned_rows: List[str] = []
    planned_count = 0
    job_payloads: List[Dict[str, object]] = []
    for index, job in enumerate(workspace.jobs, start=1):
        pair_job_name = _sanitize_token(f"{engine_settings['job_name_prefix']}-{workspace.engine}-{ctx.round_id}-{index:03d}")
        mapped_command = _map_command(job.command, ctx.root, ctx.execution_root)
        wrapper_command = [
            str(workspace.runtime.get("job_python") or "python3"),
            str(workspace.execution_engine_root / "job_scripts" / "execute_job.py"),
            str(workspace.execution_engine_root / "manifest" / "deployment_manifest.json"),
            str(index - 1),
        ]
        mapped_pose_file = _map_project_path(job.pose_file, ctx.root, ctx.execution_root)
        mapped_log_file = _map_project_path(job.log_file, ctx.root, ctx.execution_root)
        array_index = None
        if job.status == "planned":
            planned_count += 1
            array_index = planned_count
            planned_rows.append(
                "\t".join(
                    [
                        pair_job_name,
                        job.tag,
                        shlex.join(wrapper_command),
                        mapped_log_file,
                    ]
                )
            )
        job_payloads.append(
            {
                **{
                    **_map_job(job, ctx.root, ctx.execution_root),
                    "command": mapped_command,
                    "pose_file": mapped_pose_file,
                    "log_file": mapped_log_file,
                },
                "job_name": pair_job_name,
                "array_index": array_index,
                script_field: str(execution_script_path),
            }
        )
    return planned_rows, planned_count, job_payloads


def _write_project_manifest(ctx: _DeploymentContext, deployment_manifest: Dict[str, object]) -> Dict[str, object]:
    project_manifest_path = ctx.stage_root / "deployment_manifest.json"
    deployment_manifest["deployment_manifest"] = str(ctx.execution_stage_root / "deployment_manifest.json")
    with open(project_manifest_path, "w", encoding="utf-8") as handle:
        json.dump(deployment_manifest, handle, indent=2)
    return deployment_manifest


def generate_condor_deployment(
    project_root: Path,
    engines: List[str],
    pairlist_rows: List[PairlistRow],
    runtime_by_engine: Dict[str, Dict[str, object]],
    round_id: str,
    docking_mode: str,
    condor_options: Optional[Dict[str, object]] = None,
    hpc_profile: Optional[Dict[str, object]] = None,
    remote_project_root: str = "",
    skip_completed: bool = False,
    pair_source_file: str = "",
    rerun_manifest_file: str = "",
) -> Dict[str, object]:
    ctx = _open_deployment_context(
        project_root,
        remote_project_root,
        round_id,
        docking_mode,
        hpc_profile,
        skip_completed,
        pair_source_file,
        rerun_manifest_file,
        pairlist_rows,
    )
    deployment_manifest = _base_deployment_manifest(ctx, "condor", engines)

    def settings_for_engine(engine: str, runtime: Dict[str, object]) -> Dict[str, object]:
        settings = _condor_settings_for_engine(engine, runtime, condor_options or {}, profile=hpc_profile)
        _validate_condor_settings(engine, settings, profile=hpc_profile)
        return settings

    for engine in engines:
        workspace, engine_settings = _prepare_engine_workspace(
            ctx, engine, runtime_by_engine, settings_for_engine, "cpus", "condor_logs"
        )
        engine_payload = _engine_payload(ctx, workspace, "condor", engine_settings)
        submit_lines = [
            "#!/bin/bash",
            "set -euo pipefail",
            'script_dir="$(cd "$(dirname "$0")" && pwd)"',
            "",
        ]
        submit_mode = str(engine_settings.get("submit_mode") or "condor_array").strip() or "condor_array"
        script_suffix = "array" if submit_mode == "condor_array" else "batch"
        condor_script_name = f"{engine}_{script_suffix}.sh"
        condor_script_path = workspace.scripts_dir / condor_script_name
        execution_condor_script_path = workspace.execution_engine_root / "job_scripts" / condor_script_name
        condor_submit_path = workspace.scripts_dir / "job.sub"
        execution_condor_submit_path = workspace.execution_engine_root / "job_scripts" / "job.sub"
        manifest_file_name = "job_array.tsv" if submit_mode == "condor_array" else "job_manifest.tsv"
        manifest_path = workspace.manifest_dir / manifest_file_name
        execution_manifest_path = workspace.execution_engine_root / "manifest" / manifest_path.name
        planned_rows, planned_count, engine_payload["jobs"] = _planned_job_payloads(
            ctx, workspace, engine_settings, "condor_script", execution_condor_script_path
        )

        if planned_rows:
            manifest_path.write_text("\n".join(planned_rows) + "\n", encoding="utf-8")
            if submit_mode == "condor_array":
                script_text = _render_array_condor_script(
                    ctx.execution_root,
                    engine,
                    execution_manifest_path,
                    preamble_lines=list(workspace.runtime.get("deployment_preamble") or []),
                )
                submit_file_text = _render_condor_submit_file(
                    execution_condor_script_path,
                    workspace.execution_logs_dir,
                    engine_settings,
                    job_count=planned_count,
                    array_mode=True,
                )
                submit_lines.extend([
                    f"echo \"Submitting {engine} condor array job with {planned_count} task(s)\"",
                    f"condor_submit \"$script_dir/job_scripts/job.sub\"",
                ])
            else:
                script_text = _render_single_job_condor_script(
                    ctx.execution_root,
                    engine,
                    execution_manifest_path,
                    preamble_lines=list(workspace.runtime.get("deployment_preamble") or []),
                )
                submit_file_text = _render_condor_submit_file(
                    execution_condor_script_path,
                    workspace.execution_logs_dir,
                    engine_settings,
                    job_count=1,
                    array_mode=False,
                )
                submit_lines.extend([
                    f"echo \"Submitting {engine} condor single-job with {planned_count} docking command(s)\"",
                    f"condor_submit \"$script_dir/job_scripts/job.sub\"",
                ])
            condor_script_path.write_text(script_text, encoding="utf-8")
            condor_script_path.chmod(0o755)
            condor_submit_path.write_text(submit_file_text, encoding="utf-8")
        else:
            submit_lines.append(f"echo \"No planned {engine} jobs to submit\"")

        submit_path = workspace.engine_root / "submit_all.sh"
        submit_path.write_text("\n".join(submit_lines) + "\n", encoding="utf-8")
        submit_path.chmod(0o755)
        engine_payload["submit_script"] = str(workspace.execution_engine_root / "submit_all.sh")
        engine_payload["submit_mode"] = submit_mode
        engine_payload["job_count"] = planned_count
        engine_payload["job_manifest"] = str(execution_manifest_path)
        engine_payload["condor_script"] = str(execution_condor_script_path)
        engine_payload["condor_submit"] = str(execution_condor_submit_path)
        if submit_mode == "condor_array":
            engine_payload["array_job_count"] = planned_count
            engine_payload["array_parallelism"] = int(engine_settings.get("parallelism") or 0)
            engine_payload["array_manifest"] = str(execution_manifest_path)
            engine_payload["array_script"] = str(execution_condor_script_path)

        engine_manifest_path = workspace.manifest_dir / "deployment_manifest.json"
        with open(engine_manifest_path, "w", encoding="utf-8") as handle:
            json.dump(engine_payload, handle, indent=2)
        deployment_manifest["engine_jobs"][engine] = {
            "planned": planned_count,
            "skipped": sum(1 for job in workspace.jobs if job.status == "skipped"),
            "manifest": str(workspace.execution_engine_root / "manifest" / "deployment_manifest.json"),
            "submit_script": str(workspace.execution_engine_root / "submit_all.sh"),
            "submit_mode": submit_mode,
            "condor_script": str(execution_condor_script_path),
            "condor_submit": str(execution_condor_submit_path),
            "job_manifest": str(execution_manifest_path),
        }
        if submit_mode == "condor_array":
            deployment_manifest["engine_jobs"][engine]["array_script"] = str(execution_condor_script_path)
            deployment_manifest["engine_jobs"][engine]["array_manifest"] = str(execution_manifest_path)
            deployment_manifest["engine_jobs"][engine]["array_parallelism"] = int(engine_settings.get("parallelism") or 0)

    return _write_project_manifest(ctx, deployment_manifest)


def generate_slurm_deployment(
    project_root: Path,
    engines: List[str],
    pairlist_rows: List[PairlistRow],
    runtime_by_engine: Dict[str, Dict[str, object]],
    round_id: str,
    docking_mode: str,
    slurm_options: Optional[Dict[str, object]] = None,
    hpc_profile: Optional[Dict[str, object]] = None,
    remote_project_root: str = "",
    skip_completed: bool = False,
    pair_source_file: str = "",
    rerun_manifest_file: str = "",
) -> Dict[str, object]:
    ctx = _open_deployment_context(
        project_root,
        remote_project_root,
        round_id,
        docking_mode,
        hpc_profile,
        skip_completed,
        pair_source_file,
        rerun_manifest_file,
        pairlist_rows,
    )
    deployment_manifest = _base_deployment_manifest(ctx, "slurm", engines)

    def settings_for_engine(engine: str, runtime: Dict[str, object]) -> Dict[str, object]:
        settings = _slurm_settings_for_engine(engine, runtime, slurm_options or {}, profile=hpc_profile)
        _validate_slurm_settings(engine, settings, profile=hpc_profile)
        return settings

    for engine in engines:
        workspace, engine_settings = _prepare_engine_workspace(
            ctx, engine, runtime_by_engine, settings_for_engine, "cpus_per_task", "slurm_logs"
        )
        engine_payload = _engine_payload(ctx, workspace, "slurm", engine_settings)
        submit_lines = [
            "#!/bin/bash",
            "set -euo pipefail",
            'script_dir="$(cd "$(dirname "$0")" && pwd)"',
            "",
        ]
        submit_mode = str(engine_settings.get("submit_mode") or "slurm_array").strip() or "slurm_array"
        if submit_mode not in {"slurm_array", "single_job"}:
            raise ValueError(f"Unsupported slurm submit_mode for {engine}: {submit_mode}")
        job_name = _sanitize_token(f"{engine_settings['job_name_prefix']}-{engine}-{round_id}")
        script_suffix = "array" if submit_mode == "slurm_array" else "batch"
        slurm_script_path = workspace.scripts_dir / f"{engine}_{script_suffix}.slurm"
        execution_slurm_script_path = workspace.execution_engine_root / "job_scripts" / slurm_script_path.name
        manifest_file_name = "job_array.tsv" if submit_mode == "slurm_array" else "job_manifest.tsv"
        manifest_path = workspace.manifest_dir / manifest_file_name
        execution_manifest_path = workspace.execution_engine_root / "manifest" / manifest_path.name
        planned_rows, planned_count, engine_payload["jobs"] = _planned_job_payloads(
            ctx, workspace, engine_settings, "slurm_script", execution_slurm_script_path
        )
        if planned_rows:
            manifest_path.write_text("\n".join(planned_rows) + "\n", encoding="utf-8")
            if submit_mode == "slurm_array":
                slurm_script_text = _render_array_slurm_script(
                    ctx.execution_root,
                    engine,
                    job_name,
                    execution_manifest_path,
                    workspace.execution_logs_dir,
                    engine_settings,
                    job_count=planned_count,
                    preamble_lines=list(workspace.runtime.get("deployment_preamble") or []),
                )
                submit_lines.extend(
                    [
                        f"echo \"Submitting {engine} array job with {planned_count} task(s)\"",
                        f"sbatch \"$script_dir/job_scripts/{slurm_script_path.name}\"",
                    ]
                )
            else:
                slurm_script_text = _render_single_job_slurm_script(
                    ctx.execution_root,
                    engine,
                    job_name,
                    execution_manifest_path,
                    workspace.execution_logs_dir,
                    engine_settings,
                    preamble_lines=list(workspace.runtime.get("deployment_preamble") or []),
                )
                submit_lines.extend(
                    [
                        f"echo \"Submitting {engine} single-job batch with {planned_count} docking command(s)\"",
                        f"sbatch \"$script_dir/job_scripts/{slurm_script_path.name}\"",
                    ]
                )
            slurm_script_path.write_text(slurm_script_text, encoding="utf-8")
            slurm_script_path.chmod(0o755)
        else:
            submit_lines.append(f"echo \"No planned {engine} jobs to submit\"")
        submit_path = workspace.engine_root / "submit_all.sh"
        submit_path.write_text("\n".join(submit_lines) + "\n", encoding="utf-8")
        submit_path.chmod(0o755)
        engine_payload["submit_script"] = str(workspace.execution_engine_root / "submit_all.sh")
        engine_payload["submit_mode"] = submit_mode
        engine_payload["job_count"] = planned_count
        engine_payload["job_manifest"] = str(execution_manifest_path)
        engine_payload["slurm_script"] = str(execution_slurm_script_path)
        if submit_mode == "slurm_array":
            engine_payload["array_job_count"] = planned_count
            engine_payload["array_parallelism"] = int(engine_settings.get("array_parallelism") or 0)
            engine_payload["array_manifest"] = str(execution_manifest_path)
            engine_payload["array_script"] = str(execution_slurm_script_path)

        engine_manifest_path = workspace.manifest_dir / "deployment_manifest.json"
        with open(engine_manifest_path, "w", encoding="utf-8") as handle:
            json.dump(engine_payload, handle, indent=2)
        deployment_manifest["engine_jobs"][engine] = {
            "planned": planned_count,
            "skipped": sum(1 for job in workspace.jobs if job.status == "skipped"),
            "manifest": str(workspace.execution_engine_root / "manifest" / "deployment_manifest.json"),
            "submit_script": str(workspace.execution_engine_root / "submit_all.sh"),
            "submit_mode": submit_mode,
            "slurm_script": str(execution_slurm_script_path),
            "job_manifest": str(execution_manifest_path),
        }
        if submit_mode == "slurm_array":
            deployment_manifest["engine_jobs"][engine]["array_script"] = str(execution_slurm_script_path)
            deployment_manifest["engine_jobs"][engine]["array_manifest"] = str(execution_manifest_path)
            deployment_manifest["engine_jobs"][engine]["array_parallelism"] = int(engine_settings.get("array_parallelism") or 0)

    return _write_project_manifest(ctx, deployment_manifest)


def auto_execute_rerun_manifest(
    project_root: Path,
    *,
    rerun_manifest_file: Path,
    execution_mode: str = "plan_only",
    submit: bool = False,
) -> Dict[str, object]:
    """
    Future-facing stub for rerun-manifest execution orchestration.

    This intentionally returns a structured not-yet-implemented payload so that
    interactive/CLI surfaces can expose the capability without hard-wiring
    behavior until execution adapters are finalized.
    """
    root = Path(project_root).expanduser().resolve()
    rerun_file = Path(rerun_manifest_file).expanduser().resolve()
    return {
        "status": "not_implemented",
        "project_root": str(root),
        "rerun_manifest_file": str(rerun_file),
        "execution_mode": str(execution_mode or "plan_only"),
        "submit_requested": bool(submit),
        "message": (
            "Rerun-manifest auto execution is scaffolded but not active yet. "
            "Use standard deployment generation + submit commands for now."
        ),
    }
