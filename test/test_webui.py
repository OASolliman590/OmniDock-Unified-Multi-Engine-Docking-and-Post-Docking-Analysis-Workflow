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
