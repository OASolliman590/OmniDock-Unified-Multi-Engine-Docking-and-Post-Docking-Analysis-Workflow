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
    """Isolate the web UI's data directory for each test.

    Any job still running at teardown is killed. Several tests launch
    deliberately long-lived processes; leaving them behind leaks processes
    and, on Windows, keeps the log file open so pytest cannot remove
    tmp_path (a PermissionError that surfaces as an unrelated test failing).
    """
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

    yield home

    import webui.jobs as J
    for job in J.list_jobs():
        if job.is_terminal:
            continue
        try:
            J.cancel(job.job_id)
        except Exception:
            pass
    # Give the OS a moment to release the log file handles.
    for _ in range(20):
        if not any(j.pid and J._process_alive(j.pid) for j in J.list_jobs()):
            break
        time.sleep(0.1)


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


def test_exit_marker_written_when_launcher_is_gone(app_home, tmp_path):
    """A job must resolve correctly even if no watcher thread survives.

    Regression: the marker used to be written by a watcher thread inside the
    launching process. That thread dies with the server -- exactly the case
    the marker exists for -- so a finished, successful job reconciled as
    "process vanished" / failed.
    """
    import webui.jobs as J

    # Launch from a *separate* interpreter that exits immediately, leaving
    # no watcher behind, then reconcile as a fresh server would.
    script = f"""
import sys
sys.path.insert(0, {str(REPO_ROOT)!r})
import os
os.environ["OMNIDOCK_WEBUI_HOME"] = {str(app_home)!r}
import webui.config, importlib
importlib.reload(webui.config)
import webui.jobs as J
importlib.reload(J)
job = J.create_job("projZ", "t", [sys.executable, "-c", "import time; time.sleep(2)"], {str(REPO_ROOT)!r})
J.launch(job)
print(job.job_id)
"""
    out = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, timeout=60)
    assert out.returncode == 0, out.stderr
    job_id = out.stdout.strip().splitlines()[-1]

    # Wait for the supervised process to finish on its own.
    deadline = time.time() + 30
    while time.time() < deadline:
        if J._exit_marker(job_id).exists():
            break
        time.sleep(0.3)

    assert J._exit_marker(job_id).exists(), "supervisor must write the exit marker itself"

    J.reconcile_all()
    final = J.load_job(job_id)
    assert final.status == "completed", f"expected completed, got {final.status} ({final.note})"
    assert final.exit_code == 0


def test_refresh_resolves_a_job_with_no_surviving_watcher(app_home):
    """A job must resolve on read, not only at server startup.

    A job outlives the server that launched it -- the point of the design --
    so after a restart there is no watcher thread for it. Without on-read
    reconciliation the record sits at "running" forever once it finishes.
    """
    import webui.jobs as J

    job = J.launch(J.create_job("projR", "t", _sleep_argv(1), str(REPO_ROOT)))

    deadline = time.time() + 30
    while time.time() < deadline:
        if J._exit_marker(job.job_id).exists():
            break
        time.sleep(0.2)
    assert J._exit_marker(job.job_id).exists()

    # Reproduce what a restarted server sees: the record still says running
    # (the old watcher died with the previous process) but the supervisor
    # left a marker behind.
    stale = J.load_job(job.job_id)
    stale.status = "running"
    stale.exit_code = None
    stale.finished_at = None
    J.save_job(stale)

    refreshed = J.refresh(job.job_id)
    assert refreshed.status == "completed"
    assert refreshed.exit_code == 0


def test_stale_running_record_does_not_block_new_launches(app_home):
    """has_running_job must reconcile, or one stale record blocks forever."""
    import webui.jobs as J

    job = J.launch(J.create_job("projS", "t", _sleep_argv(1), str(REPO_ROOT)))

    deadline = time.time() + 30
    while time.time() < deadline:
        if J._exit_marker(job.job_id).exists():
            break
        time.sleep(0.2)

    stale = J.load_job(job.job_id)
    stale.status = "running"
    J.save_job(stale)

    assert J.has_running_job("projS") is False


def test_exit_marker_wins_over_a_recycled_pid(app_home):
    """A recorded outcome must not be shadowed by PID reuse.

    Regression: liveness was checked before the marker, so once the
    supervisor's PID was recycled by an unrelated process, a finished job
    took the "pid recycled" branch and was reported as failed.
    """
    import webui.jobs as J

    job = J.launch(J.create_job("projM", "t", _print_argv("done"), str(REPO_ROOT)))
    deadline = time.time() + 30
    while time.time() < deadline:
        if J._exit_marker(job.job_id).exists():
            break
        time.sleep(0.2)
    assert J._exit_marker(job.job_id).exists()

    # Put the record back to running and point it at a live PID whose start
    # time cannot match -- i.e. the PID was reused.
    stale = J.load_job(job.job_id)
    stale.status = "running"
    stale.exit_code = None
    stale.finished_at = None
    stale.pid = os.getpid()
    stale.proc_start_time = 1.0
    J.save_job(stale)

    resolved = J.refresh(job.job_id)
    assert resolved.status == "completed", (
        f"marker should win over liveness, got {resolved.status} ({resolved.note})")
    assert resolved.exit_code == 0


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


# ------------------------------------------------------------------ artifacts

FIXTURE_RUN = "3-Results/comparative_demo"


@pytest.fixture()
def project_id(client):
    return client.post("/api/projects", json={"path": str(FIXTURE_PROJECT)}).get_json()["id"]


@pytest.mark.parametrize("attack", [
    "../../../etc/passwd",
    "..\\..\\..\\Windows\\win.ini",
    "../",
    "subdir/../../outside.txt",
    "%2e%2e%2fetc%2fpasswd",
    "....//....//etc/passwd",
    "/etc/passwd",
    "C:\\Windows\\win.ini",
    "\\\\server\\share\\file.txt",
])
def test_path_traversal_is_rejected(client, project_id, attack):
    """SC-005: every escape attempt must be refused, none served."""
    res = client.get(f"/api/projects/{project_id}/file?path={attack}")
    assert res.status_code in (403, 404), f"{attack} returned {res.status_code}"
    if res.status_code == 403:
        assert res.get_json()["code"] == "PATH_ESCAPE"


def test_symlink_escape_is_rejected(client, project_id, tmp_path):
    """A link inside the project pointing outside it must not be followed."""
    import webui.artifacts as A

    secret = tmp_path / "secret.txt"
    secret.write_text("classified", encoding="utf-8")
    link = FIXTURE_PROJECT / "escape_link"
    try:
        link.symlink_to(secret)
    except (OSError, NotImplementedError):
        pytest.skip("symlink creation requires elevation on this platform")

    try:
        with pytest.raises(A.PathEscape):
            A.safe_resolve(FIXTURE_PROJECT, "escape_link")
    finally:
        link.unlink(missing_ok=True)


def test_safe_resolve_allows_paths_inside_the_project():
    import webui.artifacts as A

    resolved = A.safe_resolve(FIXTURE_PROJECT, f"{FIXTURE_RUN}/reports/consensus_ranked.csv")
    assert resolved.is_file()


def test_find_runs_locates_run_tracking(client, project_id):
    runs = client.get(f"/api/projects/{project_id}/runs").get_json()
    assert any(r["run_dir"] == FIXTURE_RUN for r in runs)
    assert all(r["has_manifest"] for r in runs if r["run_dir"] == FIXTURE_RUN)


def test_artifacts_include_every_outputs_index_entry(client, project_id):
    """T080: nothing declared in outputs_index.json may be dropped."""
    import json as _json

    index = _json.loads(
        (FIXTURE_PROJECT / FIXTURE_RUN / "run_tracking" / "outputs_index.json")
        .read_text(encoding="utf-8"))
    expected = {f"{FIXTURE_RUN}/{row['path']}" for row in index}

    listed = {a["path"] for a in client.get(
        f"/api/projects/{project_id}/artifacts?run={FIXTURE_RUN}").get_json()}
    assert expected.issubset(listed)


def test_artifacts_are_categorized_by_topology(client, project_id):
    items = {a["path"]: a for a in client.get(
        f"/api/projects/{project_id}/artifacts?run={FIXTURE_RUN}").get_json()}
    assert items[f"{FIXTURE_RUN}/reports/consensus_ranked.csv"]["category"] == "reports"
    assert items[f"{FIXTURE_RUN}/analysis/affinity_summary.csv"]["category"] == "analysis"
    assert items[f"{FIXTURE_RUN}/visualizations/plot.svg"]["category"] == "visualizations"
    assert items[f"{FIXTURE_RUN}/interactions/prolif/barcode_qc_rep1.csv"]["category"] \
        == "interactions.prolif"


def test_media_kinds_are_detected(client, project_id):
    items = {a["path"]: a for a in client.get(
        f"/api/projects/{project_id}/artifacts?run={FIXTURE_RUN}").get_json()}
    assert items[f"{FIXTURE_RUN}/reports/consensus_ranked.csv"]["media_kind"] == "table"
    assert items[f"{FIXTURE_RUN}/reports/summary.html"]["media_kind"] == "html"
    assert items[f"{FIXTURE_RUN}/visualizations/plot.svg"]["media_kind"] == "image"


def test_artifacts_fall_back_to_directory_walk(client, project_id, tmp_path):
    """FR-020: a run that never wrote an index must still list files."""
    import webui.artifacts as A

    run = FIXTURE_PROJECT / "3-Results" / "no_index_run"
    (run / "run_tracking").mkdir(parents=True, exist_ok=True)
    (run / "reports").mkdir(parents=True, exist_ok=True)
    (run / "reports" / "orphan.csv").write_text("a,b\n1,2\n", encoding="utf-8")
    try:
        found = A.load_outputs_index(FIXTURE_PROJECT, "3-Results/no_index_run")
        assert any(a.path.endswith("orphan.csv") for a in found)
    finally:
        import shutil
        shutil.rmtree(run, ignore_errors=True)


def test_run_summary_flags_enforced_not_matching_requested(client, project_id):
    """FR-023: an overridden stage must be visibly called out."""
    summary = client.get(
        f"/api/projects/{project_id}/runs/{FIXTURE_RUN}/summary").get_json()
    assert summary["available"] is True
    assert summary["engine"] == "gnina"

    rows = {r["stage"]: r for r in summary["stage_contract"]}
    assert rows["run_rmsd"]["requested"] is False
    assert rows["run_rmsd"]["enforced"] is True
    assert rows["run_rmsd"]["matches"] is False
    assert rows["run_visualizations"]["matches"] is True


def test_run_summary_classifies_optional_features(client, project_id):
    """FR-024: skipped-disabled, skipped-missing-dependency and failed differ."""
    summary = client.get(
        f"/api/projects/{project_id}/runs/{FIXTURE_RUN}/summary").get_json()
    features = {f["name"]: f["classification"] for f in summary["optional_features"]}
    assert features["poseview"] == "skipped_missing_dependency"
    assert features["pandamap"] == "skipped_disabled"
    assert features["py3dmol"] == "completed"
    assert features["pymol"] == "failed_error"


def test_table_endpoint_paginates(client, project_id):
    body = client.get(
        f"/api/projects/{project_id}/table"
        f"?path={FIXTURE_RUN}/reports/consensus_ranked.csv&page=1&page_size=2").get_json()
    assert body["columns"] == ["ligand", "receptor", "score", "rank"]
    assert body["total_rows"] == 3
    assert len(body["rows"]) == 2

    page2 = client.get(
        f"/api/projects/{project_id}/table"
        f"?path={FIXTURE_RUN}/reports/consensus_ranked.csv&page=2&page_size=2").get_json()
    assert len(page2["rows"]) == 1


def test_file_endpoint_serves_an_artifact(client, project_id):
    res = client.get(
        f"/api/projects/{project_id}/file?path={FIXTURE_RUN}/reports/consensus_ranked.csv")
    assert res.status_code == 200
    assert b"ligA" in res.data


def test_file_endpoint_404_for_missing_file(client, project_id):
    res = client.get(f"/api/projects/{project_id}/file?path={FIXTURE_RUN}/nope.csv")
    assert res.status_code == 404


def test_results_page_renders(client, project_id):
    assert client.get(f"/project/{project_id}/results").status_code == 200


# ------------------------------------------------------------------ multi-project


@pytest.fixture()
def two_projects(client, tmp_path):
    """Two registered projects, each a real directory."""
    ids = []
    for name in ("alpha", "beta"):
        directory = tmp_path / name
        (directory / ".workflow").mkdir(parents=True)
        (directory / ".workflow" / "state.json").write_text(
            json.dumps({
                "version": 2, "project_root": str(directory),
                "current_context": {"favorite_engine": name[:1] and "gnina"},
                "artifacts": {}, "steps": {}, "background_tasks": {},
                "feature_flags": {}, "checkpoint_metadata": {},
            }), encoding="utf-8")
        ids.append(client.post("/api/projects",
                               json={"path": str(directory), "name": name}).get_json()["id"])
    return ids


def test_concurrent_jobs_across_projects_stay_isolated(client, app_home, two_projects):
    """SC-007: two projects run independently with no cross-contamination."""
    import webui.jobs as J

    pid_a, pid_b = two_projects
    job_a = J.launch(J.create_job(pid_a, "t.a", _sleep_argv(5), str(REPO_ROOT)))
    job_b = J.launch(J.create_job(pid_b, "t.b", _sleep_argv(5), str(REPO_ROOT)))

    # Both progress; neither blocks the other.
    assert J.load_job(job_a.job_id).status == "running"
    assert J.load_job(job_b.job_id).status == "running"

    listed_a = {j["job_id"] for j in client.get(f"/api/jobs?project_id={pid_a}").get_json()}
    listed_b = {j["job_id"] for j in client.get(f"/api/jobs?project_id={pid_b}").get_json()}

    assert job_a.job_id in listed_a and job_a.job_id not in listed_b
    assert job_b.job_id in listed_b and job_b.job_id not in listed_a
    assert not (listed_a & listed_b), "project job lists must not intersect"

    # Separate log files.
    assert job_a.log_path != job_b.log_path


def test_concurrent_warning_is_scoped_per_project(client, app_home, two_projects):
    """A running job in one project must not gate launches in another."""
    import webui.jobs as J

    pid_a, pid_b = two_projects
    J.launch(J.create_job(pid_a, "t.a", _sleep_argv(20), str(REPO_ROOT)))

    assert J.has_running_job(pid_a) is True
    assert J.has_running_job(pid_b) is False


def test_dashboard_shows_only_its_own_projects_jobs(client, app_home, two_projects):
    import webui.jobs as J

    pid_a, pid_b = two_projects
    J.launch(J.create_job(pid_a, "only.in.alpha", _print_argv("a"), str(REPO_ROOT)))

    page_a = client.get(f"/project/{pid_a}").get_data(as_text=True)
    page_b = client.get(f"/project/{pid_b}").get_data(as_text=True)

    assert "only.in.alpha" in page_a
    assert "only.in.alpha" not in page_b


def test_project_switcher_lists_all_projects(client, two_projects):
    pid_a, _pid_b = two_projects
    page = client.get(f"/project/{pid_a}").get_data(as_text=True)
    assert "alpha" in page and "beta" in page


# ------------------------------------------------------------------ dag


def test_dag_reports_all_node_statuses(client, project_id):
    """SC-008: every NODE_STATUSES member must survive into the view model."""
    from post_docking_analysis.artifact_graph import NODE_STATUSES

    model = client.get(f"/api/projects/{project_id}/dag").get_json()
    assert model["available"] is True

    seen = {n["status"] for tier in model["tiers"] for n in tier["nodes"]}
    assert seen == set(NODE_STATUSES), f"missing statuses: {set(NODE_STATUSES) - seen}"
    assert all(n["known_status"] for tier in model["tiers"] for n in tier["nodes"])


def test_every_node_status_has_a_visually_distinct_style():
    """SC-008: distinguishable means distinct, not merely present.

    Parses app.css and asserts no two NODE_STATUSES share an identical
    (background, color, border) triple -- cache_hit vs completed and
    blocked_by_failure vs failed were originally rendered identically.
    """
    import re
    from post_docking_analysis.artifact_graph import NODE_STATUSES

    css = (REPO_ROOT / "webui" / "static" / "app.css").read_text(encoding="utf-8")

    styles: dict[str, str] = {}
    for status in NODE_STATUSES:
        # Match the rule block whose selector list includes .badge-<status>
        pattern = rf"(^|\n)([^{{}}]*\.badge-{re.escape(status)}\b[^{{}}]*)\{{([^}}]*)\}}"
        match = re.search(pattern, css)
        assert match, f"no CSS rule for .badge-{status}"
        body = " ".join(match.group(3).split())
        styles[status] = body

    duplicates: dict[str, list[str]] = {}
    for status, body in styles.items():
        duplicates.setdefault(body, []).append(status)

    clashes = {body: names for body, names in duplicates.items() if len(names) > 1}
    assert not clashes, f"statuses share an identical style: {clashes}"


def test_dag_status_vocabulary_comes_from_artifact_graph():
    """FR-025: the status list is imported, not restated."""
    import webui.dag as D
    from post_docking_analysis.artifact_graph import NODE_STATUSES

    model = D.to_view_model({"nodes": {}, "tiers": []})
    assert set(model["known_statuses"]) == set(NODE_STATUSES)


def test_dag_marks_cache_hits_distinctly(client, project_id):
    """FR-027: a cached node must be distinguishable from a recomputed one."""
    model = client.get(f"/api/projects/{project_id}/dag").get_json()
    nodes = {n["name"]: n for tier in model["tiers"] for n in tier["nodes"]}

    assert nodes["node_cache_hit"]["cached"] is True
    assert nodes["node_completed"]["cached"] is False
    # The cache key is read from dag_cache.json, proving the two files join.
    assert nodes["node_cache_hit"]["cache_key"] == "sha256:deadbeef"


def test_dag_exposes_blocking_edges(client, project_id):
    model = client.get(f"/api/projects/{project_id}/dag").get_json()
    nodes = {n["name"]: n for tier in model["tiers"] for n in tier["nodes"]}

    assert nodes["node_blocked_by_failure"]["blocked_by"] == ["node_failed"]
    assert {"from": "node_failed", "to": "node_blocked_by_failure",
            "kind": "blocked_by"} in model["edges"]


def test_dag_groups_nodes_into_tiers(client, project_id):
    model = client.get(f"/api/projects/{project_id}/dag").get_json()
    assert len(model["tiers"]) >= 2
    assert all("tier_index" in tier for tier in model["tiers"])


def test_dag_flags_an_unrecognized_status():
    """A status artifact_graph gains later must be surfaced, not hidden."""
    import webui.dag as D

    model = D.to_view_model({
        "nodes": {"weird": {"name": "weird", "status": "quantum_superposition"}},
        "tiers": [{"tier_index": 0, "nodes": ["weird"]}],
    })
    node = model["tiers"][0]["nodes"][0]
    assert node["status"] == "quantum_superposition"
    assert node["known_status"] is False


def test_dag_absent_report_is_not_an_error(client, tmp_path):
    project = tmp_path / "no_dag"
    project.mkdir()
    pid = client.post("/api/projects", json={"path": str(project)}).get_json()["id"]
    model = client.get(f"/api/projects/{pid}/dag").get_json()
    assert model["available"] is False
    assert client.get(f"/project/{pid}/dag").status_code == 200


def test_dag_page_renders(client, project_id):
    assert client.get(f"/project/{project_id}/dag").status_code == 200


# ------------------------------------------------------------------ hpc

SENTINELS = {
    "ssh_target": "secretuser@secret-cluster.example.edu",
    "account": "SENTINEL-SLURM-ACCOUNT",
    "project_root_base": "/cluster/users/secretuser/docking",
    "binary": "/home/secretuser/bin/gnina",
}


@pytest.fixture()
def project_with_secret_profile(client, tmp_path):
    """A project carrying an HPC profile full of known sentinel secrets."""
    project = tmp_path / "hpc_project"
    profiles = project / ".workflow" / "hpc_profiles"
    profiles.mkdir(parents=True)
    (project / ".workflow" / "state.json").write_text(json.dumps({
        "version": 2, "project_root": str(project), "current_context": {},
        "artifacts": {}, "steps": {}, "background_tasks": {},
        "feature_flags": {}, "checkpoint_metadata": {},
    }), encoding="utf-8")

    (profiles / "secret-site.json").write_text(json.dumps({
        "name": "secret-site",
        "description": "fixture profile",
        "remote": {
            "ssh_target": SENTINELS["ssh_target"],
            "project_root_base": SENTINELS["project_root_base"],
        },
        "engines": {
            "gnina": {
                "runtime": {"binary": SENTINELS["binary"]},
                "slurm_gpu": {"account": SENTINELS["account"], "partition": "gpu"},
            }
        },
    }, indent=2), encoding="utf-8")

    pid = client.post("/api/projects", json={"path": str(project)}).get_json()["id"]
    return pid, project


def test_hpc_profiles_are_redacted(client, project_with_secret_profile):
    """SC-009: no sentinel secret may appear in the API response."""
    pid, _project = project_with_secret_profile
    res = client.get(f"/api/projects/{pid}/hpc/profiles")
    assert res.status_code == 200
    body = res.get_data(as_text=True)

    for label, secret in SENTINELS.items():
        assert secret not in body, f"{label} leaked into the profiles response"
    assert "secret-site" in body  # the profile is still listed, just redacted


def test_hpc_page_does_not_render_secrets(client, project_with_secret_profile):
    pid, _project = project_with_secret_profile
    page = client.get(f"/project/{pid}/hpc").get_data(as_text=True)
    for label, secret in SENTINELS.items():
        assert secret not in page, f"{label} leaked into the rendered page"


def test_sync_requires_confirmation(client, project_with_secret_profile):
    """FR-033: a remote-touching action needs an explicit target confirmation."""
    pid, _project = project_with_secret_profile
    res = client.post(f"/api/projects/{pid}/hpc/sync", json={"profile": "secret-site"})
    assert res.status_code == 400
    assert res.get_json()["code"] == "CONFIRM_REQUIRED"


def test_submit_requires_confirmation(client, project_with_secret_profile):
    pid, _project = project_with_secret_profile
    res = client.post(f"/api/projects/{pid}/hpc/submit", json={"profile": "secret-site"})
    assert res.get_json()["code"] == "CONFIRM_REQUIRED"


def test_sync_rejects_a_mismatched_confirmation(client, project_with_secret_profile):
    pid, _project = project_with_secret_profile
    res = client.post(f"/api/projects/{pid}/hpc/sync", json={
        "profile": "secret-site", "confirm_target": "wrong@host",
    })
    assert res.status_code == 400
    assert res.get_json()["code"] == "CONFIRM_REQUIRED"


def test_deploy_does_not_require_confirmation(client, app_home, project_with_secret_profile):
    """deploy only writes local assets, so it needs no remote confirmation."""
    pid, _project = project_with_secret_profile
    res = client.post(f"/api/projects/{pid}/hpc/deploy", json={"profile": "secret-site"})
    assert res.status_code == 202, res.get_data(as_text=True)


def test_hpc_unknown_action_rejected(client, project_with_secret_profile):
    pid, _project = project_with_secret_profile
    res = client.post(f"/api/projects/{pid}/hpc/rm-rf", json={"profile": "secret-site"})
    assert res.status_code == 400


def test_redact_handles_nested_and_listed_secrets():
    import webui.hpc as H

    payload = {
        "safe": "value",
        "remote": {"ssh_target": "me@host", "port": 22},
        "list": [{"account": "acct-1"}, "plain", "/home/me/secret"],
    }
    out = H.redact(payload)
    assert out["safe"] == "value"
    assert out["remote"]["ssh_target"] == H.REDACTED
    assert out["remote"]["port"] == 22
    assert out["list"][0]["account"] == H.REDACTED
    assert out["list"][1] == "plain"
    assert out["list"][2] == H.REDACTED


def test_builtin_templates_are_listed(client, project_id):
    """The public-safe templates must be discoverable."""
    profiles = client.get(f"/api/projects/{project_id}/hpc/profiles").get_json()
    names = {p["name"] for p in profiles}
    assert any("bibalex" in n for n in names)
    assert all(p["source"] in ("builtin", "project") for p in profiles)


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


# ------------------------------------------------------------------ hardening


def test_webui_does_not_import_pipeline_internals():
    """SC-004: no pipeline stage may be reimplemented in the web layer.

    Walks each module's AST. post_docking_analysis is allowed only for
    artifact_graph (read-only status constants); workflow.execution -- where
    the stages actually run -- is off limits entirely.
    """
    import ast

    allowed_post_docking = {"post_docking_analysis.artifact_graph"}
    violations: list[str] = []

    for path in sorted((REPO_ROOT / "webui").glob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            names: list[str] = []
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                names = [node.module]

            for name in names:
                if name.startswith("workflow.execution"):
                    violations.append(f"{path.name}: imports {name}")
                if name.startswith("post_docking_analysis") and name not in allowed_post_docking:
                    violations.append(f"{path.name}: imports {name}")

    assert not violations, "web layer must not import pipeline internals: " + "; ".join(violations)


def test_every_job_runs_through_the_cli():
    """SC-004 the other way: argv is always a main.py invocation (FR-010)."""
    import webui.targets as T

    for key in T.all_targets():
        argv = T.build_argv(key, {}, "/tmp/project")
        assert argv[1] == "main.py", f"{key} does not invoke the CLI: {argv[:3]}"


def test_cli_survives_a_legacy_console_encoding():
    """`python main.py ...` must not die on a cp1252 console.

    Regression: main.py prints a non-ASCII banner, so on a Windows console
    every command -- including --help -- raised UnicodeEncodeError before
    doing any work. Runs a real subprocess with a legacy encoding forced
    and no PYTHONIOENCODING escape hatch.
    """
    env = dict(os.environ)
    env.pop("PYTHONIOENCODING", None)
    env["PYTHONLEGACYWINDOWSSTDIO"] = "1"

    result = subprocess.run(
        [sys.executable, "main.py", "webui", "--help"],
        cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=120, env=env,
    )
    combined = (result.stdout or "") + (result.stderr or "")
    assert "UnicodeEncodeError" not in combined, combined[:600]
    assert result.returncode == 0, combined[:600]


def test_webui_launcher_script_exists():
    """A launcher so the UI can be started without knowing the venv path."""
    launcher = REPO_ROOT / "start_webui.bat"
    assert launcher.is_file()
    body = launcher.read_text(encoding="utf-8")
    assert "main.py webui" in body
    assert ".venv" in body


def test_default_bind_is_localhost():
    """FR-034: the UI has no auth, so it must not default to a public bind."""
    from webui.config import HOST

    assert HOST == "127.0.0.1"


def test_remote_bind_requires_explicit_flag():
    """Binding off localhost must be refused without --allow-remote."""
    import webui_cli

    assert webui_cli.main(["--host", "0.0.0.0"]) == 2
    assert webui_cli.main(["--host", "192.168.1.10"]) == 2


def test_readonly_browsing_never_writes_into_a_project(client, tmp_path):
    """FR-035: exercise every read-only route, assert the tree is unchanged."""
    project = tmp_path / "snapshot_project"
    (project / ".workflow").mkdir(parents=True)
    (project / ".workflow" / "state.json").write_text(json.dumps({
        "version": 2, "project_root": str(project), "current_context": {},
        "artifacts": {}, "steps": {}, "background_tasks": {},
        "feature_flags": {}, "checkpoint_metadata": {},
    }), encoding="utf-8")
    (project / "reports").mkdir()
    (project / "reports" / "a.csv").write_text("x,y\n1,2\n", encoding="utf-8")

    def snapshot():
        return {
            p.relative_to(project).as_posix(): p.stat().st_mtime_ns
            for p in sorted(project.rglob("*")) if p.is_file()
        }

    before = snapshot()
    pid = client.post("/api/projects", json={"path": str(project)}).get_json()["id"]

    for route in [
        f"/api/projects/{pid}/state",
        f"/api/projects/{pid}/timeline",
        f"/api/projects/{pid}/runs",
        f"/api/projects/{pid}/artifacts",
        f"/api/projects/{pid}/dag",
        f"/api/projects/{pid}/hpc/profiles",
        f"/api/projects/{pid}/file?path=reports/a.csv",
        f"/api/projects/{pid}/table?path=reports/a.csv",
        f"/project/{pid}",
        f"/project/{pid}/launch",
        f"/project/{pid}/results",
        f"/project/{pid}/dag",
        f"/project/{pid}/hpc",
    ]:
        client.get(route)

    assert snapshot() == before, "a read-only route modified the project tree"


def test_only_app_module_imports_flask():
    """Adapters stay framework-free so they can be reused or ported."""
    import ast

    offenders = []
    for path in sorted((REPO_ROOT / "webui").glob("*.py")):
        if path.name == "app.py":
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                names = [node.module]
            if any(n.split(".")[0] == "flask" for n in names):
                offenders.append(path.name)

    assert not offenders, f"only app.py may import Flask, found in: {offenders}"


# ------------------------------------------------------------------ targets


def test_every_analysis_label_becomes_a_target(client):
    """SC-002: the catalog is derived from ANALYSIS_LABELS, not hardcoded."""
    from workflow.interactive import ANALYSIS_LABELS

    keys = {t["key"] for t in client.get("/api/targets").get_json()}
    assert set(ANALYSIS_LABELS).issubset(keys)


def test_new_analysis_label_surfaces_without_code_change(client, monkeypatch):
    """SC-002 the hard way: add a label at runtime, expect it in the API."""
    from workflow import interactive

    patched = dict(interactive.ANALYSIS_LABELS)
    patched["analyze.stage.fixture_probe"] = "Fixture probe stage"
    monkeypatch.setattr(interactive, "ANALYSIS_LABELS", patched)

    keys = {t["key"] for t in client.get("/api/targets").get_json()}
    assert "analyze.stage.fixture_probe" in keys


def test_every_target_has_a_form_schema(client):
    import webui.targets as T

    for key in T.all_targets():
        res = client.get(f"/api/targets/{key}/form")
        assert res.status_code == 200, key
        body = res.get_json()
        assert body["key"] == key
        assert isinstance(body["fields"], list)


def test_generated_argv_parses_against_real_cli():
    """Every target's argv must be accepted by workflow/cli.py's parser.

    This is what catches template drift -- e.g. the stage subcommand being
    `structure-quality` while the label key is `structure_quality`.
    """
    import webui.targets as T

    for key, target in T.all_targets().items():
        # Supply a plausible value for each required field so the parser
        # sees a complete command; we are testing the template, not
        # validation (which has its own tests).
        values = {}
        for f in target.fields:
            if not f.required:
                continue
            if f.type == "multiselect":
                values[f.name] = [f.choices[0]["value"]] if f.choices else ["gnina"]
            elif f.type == "select" and f.choices:
                values[f.name] = f.choices[0]["value"]
            elif f.type in ("int", "float"):
                values[f.name] = 1
            else:
                values[f.name] = str(FIXTURE_PROJECT)

        argv = T.build_argv(key, values, str(FIXTURE_PROJECT))
        T.assert_argv_parses(argv)  # raises TargetError on mismatch


def test_generated_argv_parses_with_all_optional_fields_supplied():
    """Optional flags must also match the CLI, not just the bare command."""
    import webui.targets as T

    for key, target in T.all_targets().items():
        values = {}
        for f in target.fields:
            if f.type == "bool":
                values[f.name] = True
            elif f.type == "multiselect":
                values[f.name] = [f.choices[0]["value"]] if f.choices else ["gnina"]
            elif f.type == "select" and f.choices:
                values[f.name] = f.choices[0]["value"]
            elif f.type in ("int", "float"):
                values[f.name] = 2
            else:
                values[f.name] = str(FIXTURE_PROJECT)

        argv = T.build_argv(key, values, str(FIXTURE_PROJECT))
        T.assert_argv_parses(argv)


def test_structure_quality_maps_to_hyphenated_subcommand():
    import webui.targets as T

    argv = T.build_argv("analyze.stage.structure_quality", {}, "/tmp/p")
    assert argv[2:5] == ["analyze", "stage", "structure-quality"]


def test_favorite_engine_maps_to_hyphenated_subcommand():
    import webui.targets as T

    argv = T.build_argv("analyze.favorite_engine", {"favorite_engine": "gnina"}, "/tmp/p")
    assert argv[2:4] == ["analyze", "favorite-engine"]


def test_interaction_aliases_declare_clean_routing(client):
    """FR-011: the UI states the reroute rather than silently substituting."""
    for key in ["analyze.interactions.prolif", "analyze.interactions.ligplot",
                "analyze.interactions.pandamap"]:
        body = client.get(f"/api/targets/{key}/form").get_json()
        assert body["routes_to"] == "analyze.interactions.clean", key


def test_form_defaults_come_from_project_context(client):
    pid = client.post("/api/projects", json={"path": str(FIXTURE_PROJECT)}).get_json()["id"]
    body = client.get(
        f"/api/targets/analyze.favorite_engine/form?project_id={pid}"
    ).get_json()
    fav = next(f for f in body["fields"] if f["name"] == "favorite_engine")
    assert fav["default"] == "gnina"  # from the fixture's current_context


def test_engine_choices_come_from_interactive_constant():
    import webui.targets as T
    from workflow.interactive import ENGINE_CHOICES

    values = {c["value"] for c in T._engine_choices()}
    expected = {v for v, _ in ENGINE_CHOICES if v != "all"}
    assert values == expected


# ------------------------------------------------------------------ launch


def test_launch_rejects_unknown_target(client):
    pid = client.post("/api/projects", json={"path": str(FIXTURE_PROJECT)}).get_json()["id"]
    res = client.post(f"/api/projects/{pid}/jobs", json={"target": "nope.nope", "values": {}})
    assert res.status_code == 400
    assert res.get_json()["code"] == "INVALID_INPUT"


def test_launch_rejects_missing_required_field(client):
    pid = client.post("/api/projects", json={"path": str(FIXTURE_PROJECT)}).get_json()["id"]
    res = client.post(f"/api/projects/{pid}/jobs",
                      json={"target": "analyze.favorite_engine", "values": {}})
    assert res.status_code == 400
    assert "favorite_engine" in res.get_json()["fields"]


def test_launch_rejects_nonexistent_path(client):
    pid = client.post("/api/projects", json={"path": str(FIXTURE_PROJECT)}).get_json()["id"]
    res = client.post(f"/api/projects/{pid}/jobs", json={
        "target": "analyze.comparative",
        "values": {"config_file": "/definitely/not/here.yaml"},
    })
    assert res.status_code == 400
    assert "config_file" in res.get_json()["fields"]


def test_launch_rejects_unknown_engine(client):
    pid = client.post("/api/projects", json={"path": str(FIXTURE_PROJECT)}).get_json()["id"]
    res = client.post(f"/api/projects/{pid}/jobs", json={
        "target": "dock.run", "values": {"engines": ["gnina", "notanengine"]},
    })
    assert res.status_code == 400
    assert "engines" in res.get_json()["fields"]


def test_launch_creates_job_with_expected_command(client):
    pid = client.post("/api/projects", json={"path": str(FIXTURE_PROJECT)}).get_json()["id"]
    res = client.post(f"/api/projects/{pid}/jobs", json={
        "target": "analyze.comparative", "values": {},
    })
    assert res.status_code == 202
    job = res.get_json()
    assert job["argv"][1:5] == ["main.py", "analyze", "comparative", "--project-dir"]
    assert job["target"] == "analyze.comparative"
    # The command shown to the user is the command that ran (FR-008)
    assert "analyze" in job["command"] and "comparative" in job["command"]


def test_second_concurrent_job_requires_confirmation(client, app_home):
    """FR-018: warn before a second job against the same project."""
    import webui.jobs as J

    pid = client.post("/api/projects", json={"path": str(FIXTURE_PROJECT)}).get_json()["id"]
    J.launch(J.create_job(pid, "t", _sleep_argv(20), str(REPO_ROOT)))

    res = client.post(f"/api/projects/{pid}/jobs",
                      json={"target": "analyze.comparative", "values": {}})
    assert res.status_code == 409
    assert res.get_json()["code"] == "JOB_RUNNING"

    ok = client.post(f"/api/projects/{pid}/jobs", json={
        "target": "analyze.comparative", "values": {}, "confirm_concurrent": True,
    })
    assert ok.status_code == 202


def test_launch_page_renders(client):
    pid = client.post("/api/projects", json={"path": str(FIXTURE_PROJECT)}).get_json()["id"]
    assert client.get(f"/project/{pid}/launch").status_code == 200


def test_preview_matches_launched_command(client):
    """FR-008: the preview must be the command that actually runs."""
    pid = client.post("/api/projects", json={"path": str(FIXTURE_PROJECT)}).get_json()["id"]
    payload = {"target": "analyze.comparative", "values": {"rescoring_top_n": 5}}

    preview = client.post(f"/api/projects/{pid}/preview", json=payload).get_json()
    launched = client.post(f"/api/projects/{pid}/jobs", json=payload).get_json()

    assert preview["argv"] == launched["argv"]
    assert preview["command"] == launched["command"]


# ------------------------------------------------------------------ logs


def test_read_from_is_incremental(tmp_path):
    """FR-015: a second read returns only the new bytes."""
    import webui.logs as L

    # Binary mode, matching how jobs.py opens the log -- text mode would
    # translate newlines on Windows and the byte offsets would not line up.
    log = tmp_path / "job.log"
    log.write_bytes(b"first\n")
    text1, offset1, _ = L.read_from(log, 0)
    assert text1 == "first\n"

    with open(log, "ab") as fh:
        fh.write(b"second\n")

    text2, offset2, size = L.read_from(log, offset1)
    assert text2 == "second\n", "second read must not repeat earlier bytes"
    assert offset2 == size


def test_read_from_past_eof_returns_nothing(tmp_path):
    import webui.logs as L

    log = tmp_path / "job.log"
    log.write_text("abc", encoding="utf-8")
    text, offset, size = L.read_from(log, 999)
    assert text == ""
    assert offset == size == 3


def test_read_from_missing_file_is_safe(tmp_path):
    import webui.logs as L

    text, offset, size = L.read_from(tmp_path / "nope.log", 0)
    assert (text, offset, size) == ("", 0, 0)


def test_read_from_respects_max_bytes(tmp_path):
    import webui.logs as L

    log = tmp_path / "big.log"
    log.write_bytes(b"x" * 10_000)
    text, offset, size = L.read_from(log, 0, max_bytes=100)
    assert len(text) == 100
    assert offset == 100
    assert size == 10_000


def test_tail_reads_only_the_end_of_a_large_log(tmp_path):
    """A multi-hour log must not be loaded whole."""
    import webui.logs as L

    log = tmp_path / "big.log"
    log.write_bytes(b"y" * (L.DEFAULT_TAIL_BYTES * 3))
    text, offset, size = L.tail(log)
    assert size == L.DEFAULT_TAIL_BYTES * 3
    assert len(text) < size
    assert offset == size


def test_log_endpoint_returns_offset_and_status(client, app_home):
    import webui.jobs as J

    job = J.launch(J.create_job("proj1", "t", _print_argv("log line"), str(REPO_ROOT)))
    for _ in range(100):
        if J.load_job(job.job_id).is_terminal:
            break
        time.sleep(0.1)

    body = client.get(f"/api/jobs/{job.job_id}/log").get_json()
    assert "log line" in body["text"]
    assert body["status"] == "completed"
    assert body["offset"] == body["size"]

    # Reading from the end returns nothing new.
    again = client.get(f"/api/jobs/{job.job_id}/log?offset={body['offset']}").get_json()
    assert again["text"] == ""


def test_log_endpoint_rejects_bad_offset(client, app_home):
    import webui.jobs as J

    job = J.create_job("proj1", "t", _print_argv("x"), str(REPO_ROOT))
    res = client.get(f"/api/jobs/{job.job_id}/log?offset=abc")
    assert res.status_code == 400


def test_stream_emits_status_event_for_finished_job(client, app_home):
    import webui.jobs as J

    job = J.launch(J.create_job("proj1", "t", _print_argv("streamed"), str(REPO_ROOT)))
    for _ in range(100):
        if J.load_job(job.job_id).is_terminal:
            break
        time.sleep(0.1)

    res = client.get(f"/api/jobs/{job.job_id}/stream")
    assert res.status_code == 200
    assert res.mimetype == "text/event-stream"
    payload = res.get_data(as_text=True)
    assert "event: status" in payload
    assert "completed" in payload


def test_job_page_renders(client, app_home):
    import webui.jobs as J

    job = J.create_job("proj1", "t", _print_argv("x"), str(REPO_ROOT))
    assert client.get(f"/job/{job.job_id}").status_code == 200


def test_job_page_404_for_unknown_job(client):
    assert client.get("/job/job_nope").status_code == 404


def test_preview_preserves_windows_style_paths(client):
    """A path with backslashes must survive into the preview intact.

    Regression: embedding the path in a JS template literal let \\t and \\f
    act as escape sequences and mangled the displayed command.
    """
    pid = client.post("/api/projects", json={"path": str(FIXTURE_PROJECT)}).get_json()["id"]
    body = client.post(f"/api/projects/{pid}/preview",
                       json={"target": "analyze.comparative", "values": {}}).get_json()
    project_path = client.get("/api/projects").get_json()[0]["path"]
    assert project_path in " ".join(body["argv"])
    assert "\t" not in body["command"] and "\f" not in body["command"]
