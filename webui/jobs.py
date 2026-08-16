"""Durable job records for work launched from the web UI.

Why this exists rather than reusing ``BackgroundTaskManager``: that manager
wraps an in-process ``ThreadPoolExecutor``, so its work dies with the host
process. Docking runs last hours and the web server must be restartable, so
web-launched work runs as a detached subprocess whose record lives on disk
and is reconciled against OS process state at startup (FR-012, FR-013).

Jobs always invoke the existing CLI. No pipeline stage is executed in this
process (FR-010).
"""

from __future__ import annotations

import json
import os
import shlex
import signal
import subprocess
import sys
import threading
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

from .config import REPO_ROOT, jobs_dir, logs_dir

try:  # optional -- only used to verify process identity
    import psutil  # type: ignore
except ImportError:  # pragma: no cover - environment dependent
    psutil = None  # type: ignore

TERMINAL_STATUSES = {"completed", "failed", "cancelled"}
VALID_STATUSES = {"queued", "running"} | TERMINAL_STATUSES

# Identity check outcomes. "unverifiable" is deliberately distinct from
# False: without a probe we must not declare a live process dead (FR-014).
IDENTITY_MATCH = "match"
IDENTITY_MISMATCH = "mismatch"
IDENTITY_UNVERIFIABLE = "unverifiable"

_TERMINATE_GRACE_SECONDS = 5.0


class JobError(Exception):
    """Raised when a job cannot be created, found, or controlled."""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass
class Job:
    job_id: str
    project_id: str
    target: str
    argv: list[str]
    command: str
    cwd: str
    log_path: str
    pid: int | None = None
    proc_start_time: float | None = None
    status: str = "queued"
    exit_code: int | None = None
    created_at: str = field(default_factory=_now)
    started_at: str | None = None
    finished_at: str | None = None
    note: str = ""

    @property
    def is_terminal(self) -> bool:
        return self.status in TERMINAL_STATUSES

    def to_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------- storage


def _job_file(job_id: str) -> Path:
    return jobs_dir() / f"{job_id}.json"


def _exit_marker(job_id: str) -> Path:
    return jobs_dir() / f"{job_id}.exit"


def save_job(job: Job) -> None:
    if job.status not in VALID_STATUSES:
        raise JobError(f"Invalid job status: {job.status}")
    path = _job_file(job.job_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(job.to_dict(), indent=2), encoding="utf-8")
    tmp.replace(path)


def load_job(job_id: str) -> Job:
    path = _job_file(job_id)
    if not path.exists():
        raise JobError(f"No job with id {job_id}")
    try:
        return Job(**json.loads(path.read_text(encoding="utf-8")))
    except (json.JSONDecodeError, TypeError) as exc:
        raise JobError(f"Job record for {job_id} is unreadable: {exc}") from exc


def list_jobs(project_id: str | None = None) -> list[Job]:
    directory = jobs_dir()
    if not directory.exists():
        return []
    jobs: list[Job] = []
    for path in directory.glob("*.json"):
        if path.name.endswith(".json.tmp"):
            continue
        try:
            job = Job(**json.loads(path.read_text(encoding="utf-8")))
        except (json.JSONDecodeError, TypeError):
            continue
        if project_id is None or job.project_id == project_id:
            jobs.append(job)
    return sorted(jobs, key=lambda j: j.created_at, reverse=True)


def has_running_job(project_id: str) -> bool:
    """Whether this project already has work in flight (FR-018)."""
    return any(not j.is_terminal for j in list_jobs(project_id))


# ---------------------------------------------------------------- process


def _process_alive(pid: int) -> bool:
    if pid is None or pid <= 0:
        return False
    if psutil is not None:
        try:
            return psutil.pid_exists(pid) and psutil.Process(pid).status() != psutil.STATUS_ZOMBIE
        except Exception:
            return False
    if sys.platform == "win32":  # pragma: no cover - platform specific
        try:
            out = subprocess.run(
                ["tasklist", "/FI", f"PID eq {pid}", "/NH"],
                capture_output=True, text=True, timeout=10,
            )
            return str(pid) in out.stdout
        except Exception:
            return False
    try:  # POSIX
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def _process_start_time(pid: int) -> float | None:
    if psutil is None:
        return None
    try:
        return float(psutil.Process(pid).create_time())
    except Exception:
        return None


def _identity_matches(job: Job) -> str:
    """Whether the live PID is still *this* job's process.

    PID existence alone is not enough: the OS recycles PIDs, so a dead job
    can otherwise read as running forever (FR-014).
    """
    if job.pid is None:
        return IDENTITY_MISMATCH
    observed = _process_start_time(job.pid)
    if observed is None or job.proc_start_time is None:
        return IDENTITY_UNVERIFIABLE
    # Filesystem/platform clocks vary slightly; a 2s window is generous for
    # identity while still far below any realistic PID-recycle interval.
    return IDENTITY_MATCH if abs(observed - job.proc_start_time) < 2.0 else IDENTITY_MISMATCH


# ---------------------------------------------------------------- lifecycle


def create_job(project_id: str, target: str, argv: list[str], cwd: str) -> Job:
    """Persist a job record before anything is launched (FR-012)."""
    job_id = f"job_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')}_{uuid4().hex[:8]}"
    job = Job(
        job_id=job_id,
        project_id=project_id,
        target=target,
        argv=list(argv),
        command=" ".join(shlex.quote(a) for a in argv),
        cwd=str(cwd),
        log_path=str(logs_dir() / f"{job_id}.log"),
    )
    save_job(job)
    return job


def _job_env() -> dict:
    """Environment for a launched job.

    PYTHONIOENCODING is set because main.py prints a non-ASCII banner; on a
    Windows console defaulting to cp1252 that raises UnicodeEncodeError and
    every job would die before doing any work.

    OMNIDOCK_NON_INTERACTIVE plus stdin=DEVNULL means a stray prompt fails
    fast instead of hanging forever holding the job open (FR-019).
    """
    env = os.environ.copy()
    env["OMNIDOCK_NON_INTERACTIVE"] = "1"
    env["PYTHONUNBUFFERED"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"
    return env


def launch(job: Job) -> Job:
    """Start the job's subprocess, detached from the request cycle."""
    log_file = Path(job.log_path)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    try:
        handle = open(log_file, "ab", buffering=0)
    except OSError as exc:
        job.status = "failed"
        job.note = f"could not open log file: {exc}"
        job.finished_at = _now()
        save_job(job)
        raise JobError(job.note) from exc

    try:
        proc = subprocess.Popen(
            job.argv,
            cwd=job.cwd,
            stdout=handle,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            env=_job_env(),
            close_fds=True,
        )
    except OSError as exc:
        handle.close()
        job.status = "failed"
        job.note = f"launch failed: {exc}"
        job.finished_at = _now()
        save_job(job)
        raise JobError(job.note) from exc

    job.pid = proc.pid
    job.proc_start_time = _process_start_time(proc.pid)
    job.status = "running"
    job.started_at = _now()
    save_job(job)

    _watch(job.job_id, proc, handle)
    return job


def _watch(job_id: str, proc: subprocess.Popen, handle) -> None:
    """Record the exit code when the process ends.

    The marker file is what lets a *restarted* server learn the outcome of a
    process it never owned a handle to (T032).
    """

    def _run() -> None:
        try:
            code = proc.wait()
        finally:
            try:
                handle.close()
            except Exception:
                pass
        try:
            _exit_marker(job_id).write_text(str(code), encoding="utf-8")
        except OSError:
            pass
        try:
            job = load_job(job_id)
        except JobError:
            return
        if job.status == "cancelled":
            return
        job.exit_code = code
        job.status = "completed" if code == 0 else "failed"
        if code != 0 and not job.note:
            job.note = f"exited with code {code}"
        job.finished_at = _now()
        save_job(job)

    threading.Thread(target=_run, name=f"job-watch-{job_id}", daemon=True).start()


def cancel(job_id: str) -> Job:
    """Terminate a running job, leaving any artifacts it produced in place (FR-016)."""
    job = load_job(job_id)
    if job.is_terminal:
        return job

    if job.pid and _process_alive(job.pid):
        try:
            if sys.platform == "win32":  # pragma: no cover - platform specific
                subprocess.run(["taskkill", "/PID", str(job.pid), "/T", "/F"],
                               capture_output=True, timeout=15)
            else:
                os.kill(job.pid, signal.SIGTERM)
                deadline = time.time() + _TERMINATE_GRACE_SECONDS
                while time.time() < deadline and _process_alive(job.pid):
                    time.sleep(0.2)
                if _process_alive(job.pid):
                    os.kill(job.pid, signal.SIGKILL)
        except (OSError, subprocess.SubprocessError) as exc:
            job.note = f"cancel signal failed: {exc}"

    job.status = "cancelled"
    job.finished_at = _now()
    if not job.note:
        job.note = "cancelled by user"
    save_job(job)
    return job


def reconcile_all() -> list[Job]:
    """Bring non-terminal job records back in line with reality at startup.

    Implements the four-case rule in plan.md / data-model.md.
    """
    reconciled: list[Job] = []
    for job in list_jobs():
        if job.is_terminal:
            continue

        marker = _exit_marker(job.job_id)

        if job.pid is None:
            job.status = "failed"
            job.note = "never started"
            job.finished_at = _now()

        elif not _process_alive(job.pid):
            if marker.exists():
                try:
                    code = int(marker.read_text(encoding="utf-8").strip())
                except (OSError, ValueError):
                    code = None
                if code is None:
                    job.status = "failed"
                    job.note = "process ended; exit code unreadable"
                else:
                    job.exit_code = code
                    job.status = "completed" if code == 0 else "failed"
                    if code != 0:
                        job.note = f"exited with code {code}"
            else:
                job.status = "failed"
                job.note = "process vanished"
            job.finished_at = _now()

        else:
            identity = _identity_matches(job)
            if identity == IDENTITY_MISMATCH:
                job.status = "failed"
                job.note = "pid recycled by an unrelated process"
                job.finished_at = _now()
            elif identity == IDENTITY_UNVERIFIABLE:
                # Alive but unprovable: keep it running rather than killing a
                # real run on a guess. Record that verification was impossible.
                job.note = "still running (identity unverifiable: psutil unavailable)"
            else:
                continue  # genuinely still running -- nothing to change

        save_job(job)
        reconciled.append(job)
    return reconciled
