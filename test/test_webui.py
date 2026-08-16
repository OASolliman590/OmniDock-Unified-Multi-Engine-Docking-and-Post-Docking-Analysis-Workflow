"""Tests for the OmniDock web UI layer.

Run with:  pytest test/test_webui.py -v
"""

from __future__ import annotations

import importlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

FIXTURE_PROJECT = REPO_ROOT / "test" / "fixtures" / "webui_project"


@pytest.fixture()
def app_home(monkeypatch, tmp_path):
    """Isolate the web UI's data directory for each test."""
    home = tmp_path / "webui_home"
    monkeypatch.setenv("OMNIDOCK_WEBUI_HOME", str(home))
    # config reads the env var lazily, but modules cache nothing -- reimport
    # to be safe if a previous test already touched them.
    import webui.config
    importlib.reload(webui.config)
    import webui.registry
    import webui.jobs
    importlib.reload(webui.registry)
    importlib.reload(webui.jobs)
    webui.config.ensure_app_dirs()
    return home


@pytest.fixture()
def client(app_home):
    import webui.app
    importlib.reload(webui.app)
    app = webui.app.create_app()
    app.config.update(TESTING=True)
    return app.test_client()


# ------------------------------------------------------------------ registry


def test_health(client):
    res = client.get("/api/health")
    assert res.status_code == 200
    assert res.get_json()["ok"] is True


def test_register_project(client):
    res = client.post("/api/projects", json={"path": str(FIXTURE_PROJECT)})
    assert res.status_code == 201
    body = res.get_json()
    assert body["has_workflow"] is True
    assert body["available"] is True


def test_duplicate_path_spellings_register_once(client):
    """Trailing separators and case variants must not create extra entries."""
    base = str(FIXTURE_PROJECT)
    client.post("/api/projects", json={"path": base})
    client.post("/api/projects", json={"path": base + os.sep})
    client.post("/api/projects", json={"path": base.upper()})
    assert len(client.get("/api/projects").get_json()) == 1


def test_register_nonexistent_path_rejected(client):
    res = client.post("/api/projects", json={"path": str(FIXTURE_PROJECT / "nope")})
    assert res.status_code == 400
    assert res.get_json()["code"] == "INVALID_INPUT"


def test_unavailable_project_does_not_error(client, tmp_path):
    target = tmp_path / "vanishing"
    target.mkdir()
    pid = client.post("/api/projects", json={"path": str(target)}).get_json()["id"]
    shutil.rmtree(target)
    assert client.get("/api/projects").get_json()[0]["available"] is False
    assert client.get(f"/project/{pid}").status_code == 200


# ------------------------------------------------------------------ state


def test_timeline_matches_timeline_steps_order(client):
    from workflow.interactive import TIMELINE_STEPS

    pid = client.post("/api/projects", json={"path": str(FIXTURE_PROJECT)}).get_json()["id"]
    rows = client.get(f"/api/projects/{pid}/timeline").get_json()
    assert [r["key"] for r in rows] == [k for k, _ in TIMELINE_STEPS]


def test_timeline_reflects_recorded_status(client):
    pid = client.post("/api/projects", json={"path": str(FIXTURE_PROJECT)}).get_json()["id"]
    rows = {r["key"]: r for r in client.get(f"/api/projects/{pid}/timeline").get_json()}
    assert rows["pdb.collect"]["status"] == "completed"
    assert rows["dock.run"]["status"] == "failed"
    assert rows["pdb.prepare_both"]["status"] == "pending"


def test_readonly_routes_never_write_into_project(client, tmp_path):
    """FR-035: the web process must not create .workflow/ in a project.

    The workflow.state accessors create state.json when it is absent, so
    this guards the state_adapter's existence check.
    """
    project = tmp_path / "untouched"
    project.mkdir()
    pid = client.post("/api/projects", json={"path": str(project)}).get_json()["id"]

    client.get(f"/api/projects/{pid}/state")
    client.get(f"/api/projects/{pid}/timeline")
    client.get(f"/project/{pid}")

    assert list(project.iterdir()) == [], "web layer wrote into the project directory"


def test_uninitialized_project_reports_no_workflow(client, tmp_path):
    project = tmp_path / "fresh"
    project.mkdir()
    pid = client.post("/api/projects", json={"path": str(project)}).get_json()["id"]
    assert client.get(f"/api/projects/{pid}/state").get_json()["has_workflow"] is False


# ------------------------------------------------------------------ jobs


def _sleep_argv(seconds: float) -> list[str]:
    return [sys.executable, "-c", f"import time; time.sleep({seconds})"]


def _print_argv(text: str) -> list[str]:
    return [sys.executable, "-c", f"print({text!r})"]


def test_job_record_persisted_before_launch(app_home):
    import webui.jobs as J

    job = J.create_job("proj1", "test.target", _print_argv("hi"), str(REPO_ROOT))
    assert job.status == "queued"
    assert J.load_job(job.job_id).status == "queued", "record must exist before launch"


def test_job_completes_with_exit_zero(app_home):
    import webui.jobs as J

    job = J.launch(J.create_job("proj1", "t", _print_argv("hello"), str(REPO_ROOT)))
    for _ in range(100):
        job = J.load_job(job.job_id)
        if job.is_terminal:
            break
        time.sleep(0.1)
    assert job.status == "completed"
    assert job.exit_code == 0
    assert "hello" in Path(job.log_path).read_text(encoding="utf-8")


def test_job_failure_records_exit_code(app_home):
    import webui.jobs as J

    argv = [sys.executable, "-c", "import sys; sys.exit(3)"]
    job = J.launch(J.create_job("proj1", "t", argv, str(REPO_ROOT)))
    for _ in range(100):
        job = J.load_job(job.job_id)
        if job.is_terminal:
            break
        time.sleep(0.1)
    assert job.status == "failed"
    assert job.exit_code == 3


def test_restart_durability(app_home):
    """SC-003: a job must survive a web-server restart with correct status."""
    import webui.jobs as J

    job = J.launch(J.create_job("proj1", "t", _sleep_argv(3), str(REPO_ROOT)))
    assert job.status == "running"

    # Simulate a server restart: drop all in-process state (the watcher
    # thread included) by reimporting the module, then reconcile from disk.
    importlib.reload(J)
    reconciled = J.reconcile_all()
    after = J.load_job(job.job_id)
    assert after.status == "running", (
        f"job should still be running after restart, got {after.status} ({after.note})"
    )
    assert not any(r.job_id == job.job_id and r.status == "failed" for r in reconciled)

    # Let it finish; a fresh reconcile must pick up the real outcome.
    deadline = time.time() + 20
    while time.time() < deadline:
        if not J._process_alive(after.pid):
            break
        time.sleep(0.2)
    time.sleep(0.5)
    J.reconcile_all()
    final = J.load_job(job.job_id)
    assert final.status == "completed", f"expected completed, got {final.status} ({final.note})"


def test_reconcile_marks_recycled_pid_as_failed(app_home):
    """FR-014: PID liveness alone must not be trusted."""
    import webui.jobs as J

    job = J.create_job("proj1", "t", _sleep_argv(0.1), str(REPO_ROOT))
    # A live PID (this test process) recorded with a start time that cannot
    # belong to it -- exactly the recycled-PID situation.
    job.pid = os.getpid()
    job.proc_start_time = 1.0
    job.status = "running"
    job.started_at = "2020-01-01T00:00:00+00:00"
    J.save_job(job)

    J.reconcile_all()
    after = J.load_job(job.job_id)
    assert after.status == "failed"
    assert "recycled" in after.note


def test_reconcile_missing_pid_marks_failed(app_home):
    import webui.jobs as J

    job = J.create_job("proj1", "t", _print_argv("x"), str(REPO_ROOT))
    job.status = "running"  # never actually launched
    J.save_job(job)
    J.reconcile_all()
    after = J.load_job(job.job_id)
    assert after.status == "failed"
    assert after.note == "never started"


def test_cancel_running_job(app_home):
    import webui.jobs as J

    job = J.launch(J.create_job("proj1", "t", _sleep_argv(30), str(REPO_ROOT)))
    cancelled = J.cancel(job.job_id)
    assert cancelled.status == "cancelled"
    time.sleep(0.5)
    assert J.load_job(job.job_id).status == "cancelled", "watcher must not overwrite a cancel"


def test_has_running_job(app_home):
    import webui.jobs as J

    assert J.has_running_job("projX") is False
    J.launch(J.create_job("projX", "t", _sleep_argv(5), str(REPO_ROOT)))
    assert J.has_running_job("projX") is True


def test_jobs_scoped_by_project(app_home):
    """SC-007: one project's job list must never include another's."""
    import webui.jobs as J

    J.launch(J.create_job("projA", "t", _print_argv("a"), str(REPO_ROOT)))
    J.launch(J.create_job("projB", "t", _print_argv("b"), str(REPO_ROOT)))
    assert {j.project_id for j in J.list_jobs("projA")} == {"projA"}
    assert {j.project_id for j in J.list_jobs("projB")} == {"projB"}


def test_job_env_is_non_interactive_and_utf8(app_home):
    """FR-019 plus the Windows banner-encoding guard."""
    import webui.jobs as J

    env = J._job_env()
    assert env["OMNIDOCK_NON_INTERACTIVE"] == "1"
    assert env["PYTHONUNBUFFERED"] == "1"
    # main.py prints a non-ASCII banner; without this every job dies on a
    # cp1252 Windows console before doing any work.
    assert env["PYTHONIOENCODING"] == "utf-8"


def test_job_stdin_is_closed(app_home):
    """A stray prompt must fail fast rather than hang (FR-019)."""
    import webui.jobs as J

    argv = [sys.executable, "-c", "input('prompt: ')"]
    job = J.launch(J.create_job("proj1", "t", argv, str(REPO_ROOT)))
    deadline = time.time() + 20
    while time.time() < deadline:
        job = J.load_job(job.job_id)
        if job.is_terminal:
            break
        time.sleep(0.1)
    assert job.is_terminal, "job reading stdin should not hang"
    assert job.status == "failed"


def test_main_py_runs_under_job_env(app_home):
    """The real CLI must survive the job environment on this platform."""
    import webui.jobs as J

    argv = [sys.executable, "main.py", "--help"]
    job = J.launch(J.create_job("proj1", "cli.help", argv, str(REPO_ROOT)))
    deadline = time.time() + 60
    while time.time() < deadline:
        job = J.load_job(job.job_id)
        if job.is_terminal:
            break
        time.sleep(0.2)
    log = Path(job.log_path).read_text(encoding="utf-8", errors="replace")
    assert "UnicodeEncodeError" not in log, f"job env failed to handle CLI output:\n{log[:500]}"


# ------------------------------------------------------------------ api


def test_job_api_roundtrip(client):
    import webui.jobs as J

    job = J.launch(J.create_job("proj1", "t", _print_argv("api"), str(REPO_ROOT)))
    res = client.get(f"/api/jobs/{job.job_id}")
    assert res.status_code == 200
    assert res.get_json()["job_id"] == job.job_id
    assert client.get("/api/jobs?project_id=proj1").status_code == 200


def test_unknown_job_returns_404(client):
    res = client.get("/api/jobs/job_does_not_exist")
    assert res.status_code == 404
    assert res.get_json()["code"] == "NOT_FOUND"
