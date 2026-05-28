from __future__ import annotations

import re
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .hpc_profiles import resolve_remote_project_dir, resolve_ssh_target


DEFAULT_SYNC_EXCLUDES = [
    ".DS_Store",
    "__pycache__/",
    ".workflow/state.json",
]
SBATCH_JOB_RE = re.compile(r"Submitted batch job (\d+)")
CONDOR_CLUSTER_RE = re.compile(r"submitted to cluster (\d+)")


def _merge_sync_excludes(profile: Dict[str, object], extra_excludes: Optional[List[str]] = None) -> List[str]:
    remote = profile.get("remote", {}) if isinstance(profile, dict) else {}
    configured = remote.get("sync_excludes", []) if isinstance(remote, dict) else []
    seen = set()
    merged: List[str] = []
    for candidate in [*DEFAULT_SYNC_EXCLUDES, *configured, *(extra_excludes or [])]:
        token = str(candidate).strip()
        if token and token not in seen:
            seen.add(token)
            merged.append(token)
    return merged


def resolve_remote_target(
    local_project_root: Path,
    profile: Dict[str, object],
    ssh_target: str = "",
    remote_project_dir: str = "",
) -> Tuple[str, str]:
    resolved_ssh_target = resolve_ssh_target(profile, ssh_target)
    resolved_remote_dir = resolve_remote_project_dir(local_project_root, profile, remote_project_dir)
    if not resolved_ssh_target:
        raise ValueError("No SSH target configured. Provide --ssh-target or set remote.ssh_target in the local HPC profile.")
    if not resolved_remote_dir:
        raise ValueError(
            "No remote project directory configured. Provide --remote-project-dir or set remote.project_root/project_root_base in the local HPC profile."
        )
    return resolved_ssh_target, resolved_remote_dir


def build_sync_commands(
    local_project_root: Path,
    ssh_target: str,
    remote_project_dir: str,
    profile: Dict[str, object],
    delete: bool = False,
    extra_excludes: Optional[List[str]] = None,
) -> Dict[str, List[str]]:
    local_root = Path(local_project_root).expanduser().resolve()
    if shutil.which("rsync") is None:
        raise RuntimeError("rsync is required for dock sync but was not found on PATH.")

    mkdir_cmd = ["ssh", ssh_target, "mkdir", "-p", remote_project_dir]
    rsync_cmd = ["rsync", "-az"]
    if delete:
        rsync_cmd.append("--delete")
    for pattern in _merge_sync_excludes(profile, extra_excludes):
        rsync_cmd.extend(["--exclude", pattern])
    rsync_cmd.extend([f"{str(local_root)}/", f"{ssh_target}:{remote_project_dir}/"])
    return {"mkdir": mkdir_cmd, "rsync": rsync_cmd}


def sync_project_to_remote(
    local_project_root: Path,
    ssh_target: str,
    remote_project_dir: str,
    profile: Dict[str, object],
    delete: bool = False,
    dry_run: bool = False,
    extra_excludes: Optional[List[str]] = None,
) -> Dict[str, object]:
    commands = build_sync_commands(
        local_project_root=local_project_root,
        ssh_target=ssh_target,
        remote_project_dir=remote_project_dir,
        profile=profile,
        delete=delete,
        extra_excludes=extra_excludes,
    )
    payload = {
        "ssh_target": ssh_target,
        "remote_project_dir": remote_project_dir,
        "commands": {name: shlex.join(cmd) for name, cmd in commands.items()},
        "dry_run": dry_run,
    }
    if dry_run:
        payload["status"] = "dry_run"
        return payload

    mkdir_result = subprocess.run(commands["mkdir"], capture_output=True, text=True)
    if mkdir_result.returncode != 0:
        payload["status"] = "failed"
        payload["failed_command"] = "mkdir"
        payload["stderr"] = mkdir_result.stderr.strip()
        return payload

    rsync_result = subprocess.run(commands["rsync"], capture_output=True, text=True)
    payload["status"] = "completed" if rsync_result.returncode == 0 else "failed"
    payload["stdout"] = rsync_result.stdout
    payload["stderr"] = rsync_result.stderr.strip()
    payload["returncode"] = rsync_result.returncode
    return payload


def infer_round_mode_from_deployment_root(deployment_root: str) -> Tuple[str, str]:
    path = Path(deployment_root)
    if not deployment_root:
        return "", ""
    mode = path.name if path.name in {"screen", "exhaustive"} else ""
    round_id = path.parent.name if mode and path.parent.name else ""
    return round_id, mode


def build_submit_commands(
    ssh_target: str,
    remote_project_dir: str,
    round_id: str,
    mode: str,
    engines: List[str],
) -> List[List[str]]:
    commands: List[List[str]] = []
    for engine in engines:
        remote_script = (
            Path(remote_project_dir)
            / "4-Docking"
            / "deployments"
            / round_id
            / mode
            / engine
            / "submit_all.sh"
        )
        remote_command = f"test -f {shlex.quote(str(remote_script))} && bash {shlex.quote(str(remote_script))}"
        commands.append(["ssh", ssh_target, f"bash -lc {shlex.quote(remote_command)}"])
    return commands


def submit_remote_deployment(
    ssh_target: str,
    remote_project_dir: str,
    round_id: str,
    mode: str,
    engines: List[str],
    dry_run: bool = False,
) -> Dict[str, object]:
    commands = build_submit_commands(
        ssh_target=ssh_target,
        remote_project_dir=remote_project_dir,
        round_id=round_id,
        mode=mode,
        engines=engines,
    )
    payload = {
        "ssh_target": ssh_target,
        "remote_project_dir": remote_project_dir,
        "round_id": round_id,
        "mode": mode,
        "engines": engines,
        "commands": [shlex.join(cmd) for cmd in commands],
        "dry_run": dry_run,
        "results": [],
    }
    if dry_run:
        payload["status"] = "dry_run"
        return payload

    overall = "completed"
    for engine, command in zip(engines, commands):
        result = subprocess.run(command, capture_output=True, text=True)
        stdout = result.stdout or ""
        job_ids = SBATCH_JOB_RE.findall(stdout) or CONDOR_CLUSTER_RE.findall(stdout)
        payload["results"].append(
            {
                "engine": engine,
                "command": shlex.join(command),
                "status": "completed" if result.returncode == 0 else "failed",
                "returncode": result.returncode,
                "stdout": stdout,
                "stderr": result.stderr.strip(),
                "job_ids": job_ids,
                "job_id": job_ids[-1] if job_ids else "",
            }
        )
        if result.returncode != 0:
            overall = "failed"
    payload["status"] = overall
    return payload
