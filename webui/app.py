"""Flask application factory and route registration.

Routes only -- all behavior lives in the adapter modules. This is the one
module in the package permitted to import Flask.

State is read fresh on every request. The web UI must never serve a cached
copy, because the CLI and the interactive shell write the same project
state concurrently from other processes (FR-005).
"""

from __future__ import annotations

import json
import shlex
import sys
import time

from flask import Flask, Response, abort, jsonify, render_template, request, send_file

from . import __version__
from . import artifacts as artifacts_mod
from . import jobs as jobs_mod
from . import logs as logs_mod
from . import registry
from . import state_adapter
from . import targets as targets_mod
from .config import REPO_ROOT, ensure_app_dirs


def _project_or_404(pid: str) -> registry.Project:
    try:
        return registry.get_project(pid)
    except registry.RegistryError:
        abort(404, description=f"No registered project with id {pid}")


def _err(message: str, code: str, status: int):
    return jsonify({"ok": False, "error": message, "code": code}), status


def create_app() -> Flask:
    app = Flask(__name__)
    ensure_app_dirs()

    # Bring job records back in line with OS reality. Jobs outlive this
    # process, so records left mid-flight by a previous run must be
    # resolved before anything is served (FR-013).
    jobs_mod.reconcile_all()

    # ---------- health ----------

    @app.get("/api/health")
    def api_health():
        return jsonify({"ok": True, "version": __version__})

    # ---------- projects ----------

    @app.get("/api/projects")
    def api_projects_list():
        return jsonify([p.to_dict() for p in registry.list_projects()])

    @app.post("/api/projects")
    def api_projects_create():
        body = request.get_json(silent=True) or {}
        path = (body.get("path") or "").strip()
        if not path:
            return _err("A project path is required.", "INVALID_INPUT", 400)
        try:
            project = registry.register(path, body.get("name"))
        except registry.RegistryError as exc:
            return _err(str(exc), "INVALID_INPUT", 400)
        return jsonify(project.to_dict()), 201

    @app.delete("/api/projects/<pid>")
    def api_projects_delete(pid: str):
        try:
            registry.unregister(pid)
        except registry.RegistryError as exc:
            return _err(str(exc), "NOT_FOUND", 404)
        return jsonify({"ok": True})

    @app.get("/api/projects/<pid>/state")
    def api_project_state(pid: str):
        project = _project_or_404(pid)
        if not project.available:
            return _err("Project directory is unavailable.", "PROJECT_UNAVAILABLE", 409)
        try:
            return jsonify(state_adapter.project_state(project.root))
        except state_adapter.StateUnreadable as exc:
            return _err(str(exc), "STATE_UNREADABLE", 503)

    @app.get("/api/projects/<pid>/timeline")
    def api_project_timeline(pid: str):
        project = _project_or_404(pid)
        if not project.available:
            return _err("Project directory is unavailable.", "PROJECT_UNAVAILABLE", 409)
        try:
            return jsonify(state_adapter.timeline(project.root))
        except state_adapter.StateUnreadable as exc:
            return _err(str(exc), "STATE_UNREADABLE", 503)

    @app.post("/api/projects/<pid>/init")
    def api_project_init(pid: str):
        """Initialize workflow state via the CLI (FR-004).

        Runs as a subprocess: the web process itself must not write into a
        project directory (FR-035).
        """
        project = _project_or_404(pid)
        if not project.available:
            return _err("Project directory is unavailable.", "PROJECT_UNAVAILABLE", 409)

        argv = [
            sys.executable, "main.py", "workflow", "init",
            "--project-dir", project.path,
        ]
        job = jobs_mod.create_job(project.id, "workflow.init", argv, str(REPO_ROOT))
        try:
            job = jobs_mod.launch(job)
        except jobs_mod.JobError as exc:
            return _err(str(exc), "INVALID_INPUT", 500)
        return jsonify(job.to_dict()), 202

    # ---------- targets ----------

    @app.get("/api/targets")
    def api_targets():
        return jsonify([
            {"key": t.key, "label": t.label, "group": t.group, "routes_to": t.routes_to}
            for t in targets_mod.all_targets().values()
        ])

    @app.get("/api/targets/<path:key>/form")
    def api_target_form(key: str):
        context: dict = {}
        project_id = request.args.get("project_id")
        if project_id:
            try:
                project = registry.get_project(project_id)
                if project.available:
                    context = state_adapter.project_state(project.root).get("current_context", {})
            except (registry.RegistryError, state_adapter.StateUnreadable):
                context = {}
        try:
            return jsonify(targets_mod.form_schema(key, context))
        except targets_mod.TargetError as exc:
            return _err(str(exc), "NOT_FOUND", 404)

    @app.post("/api/projects/<pid>/preview")
    def api_launch_preview(pid: str):
        """Render the exact command a launch would run.

        Built by the same build_argv the launcher uses, so the preview can
        never drift from what actually executes (FR-008). Reconstructing it
        in the browser would be a second implementation.
        """
        project = _project_or_404(pid)
        body = request.get_json(silent=True) or {}
        key = (body.get("target") or "").strip()
        try:
            argv = targets_mod.build_argv(key, body.get("values") or {}, project.path)
        except targets_mod.TargetError as exc:
            return _err(str(exc), "INVALID_INPUT", 400)
        return jsonify({"argv": argv, "command": " ".join(shlex.quote(a) for a in argv)})

    # ---------- jobs ----------

    @app.post("/api/projects/<pid>/jobs")
    def api_project_launch(pid: str):
        project = _project_or_404(pid)
        if not project.available:
            return _err("Project directory is unavailable.", "PROJECT_UNAVAILABLE", 409)

        body = request.get_json(silent=True) or {}
        key = (body.get("target") or "").strip()
        values = body.get("values") or {}

        try:
            target = targets_mod.get_target(key)
        except targets_mod.TargetError as exc:
            return _err(str(exc), "INVALID_INPUT", 400)

        errors = targets_mod.validate(key, values)
        if errors:
            return jsonify({
                "ok": False, "code": "INVALID_INPUT",
                "error": "Some fields need attention.", "fields": errors,
            }), 400

        # Pipeline stages are not guaranteed concurrency-safe within one
        # project, so a second job needs explicit confirmation (FR-018).
        if jobs_mod.has_running_job(pid) and not body.get("confirm_concurrent"):
            return jsonify({
                "ok": False, "code": "JOB_RUNNING",
                "error": "This project already has a job running. "
                         "Resubmit with confirm_concurrent to run both.",
            }), 409

        argv = targets_mod.build_argv(key, values, project.path)
        try:
            targets_mod.assert_argv_parses(argv)
        except targets_mod.TargetError as exc:
            return _err(str(exc), "INVALID_INPUT", 400)

        job = jobs_mod.create_job(pid, target.key, argv, str(REPO_ROOT))
        try:
            job = jobs_mod.launch(job)
        except jobs_mod.JobError as exc:
            return _err(str(exc), "INVALID_INPUT", 500)
        return jsonify(job.to_dict()), 202

    @app.get("/api/jobs")
    def api_jobs_list():
        project_id = request.args.get("project_id")
        jobs = [jobs_mod.refresh(j.job_id) for j in jobs_mod.list_jobs(project_id)]
        return jsonify([j.to_dict() for j in jobs])

    @app.get("/api/jobs/<job_id>")
    def api_job_detail(job_id: str):
        try:
            return jsonify(jobs_mod.refresh(job_id).to_dict())
        except jobs_mod.JobError as exc:
            return _err(str(exc), "NOT_FOUND", 404)

    @app.post("/api/jobs/<job_id>/cancel")
    def api_job_cancel(job_id: str):
        try:
            return jsonify(jobs_mod.cancel(job_id).to_dict())
        except jobs_mod.JobError as exc:
            return _err(str(exc), "NOT_FOUND", 404)

    @app.get("/api/jobs/<job_id>/log")
    def api_job_log(job_id: str):
        try:
            job = jobs_mod.refresh(job_id)
        except jobs_mod.JobError as exc:
            return _err(str(exc), "NOT_FOUND", 404)

        raw_offset = request.args.get("offset")
        if raw_offset is None:
            text, offset, size = logs_mod.tail(job.log_path)
        else:
            try:
                start = int(raw_offset)
            except ValueError:
                return _err("offset must be an integer", "INVALID_INPUT", 400)
            text, offset, size = logs_mod.read_from(job.log_path, start)

        return jsonify({"text": text, "offset": offset, "size": size, "status": job.status})

    @app.get("/api/jobs/<job_id>/stream")
    def api_job_stream(job_id: str):
        """Server-sent events: new log chunks, then a terminal status event."""
        try:
            job = jobs_mod.load_job(job_id)
        except jobs_mod.JobError as exc:
            return _err(str(exc), "NOT_FOUND", 404)

        try:
            start = int(request.args.get("offset", "0"))
        except ValueError:
            start = 0

        def generate():
            offset = start
            idle = 0.0
            while True:
                text, offset, _size = logs_mod.read_from(job.log_path, offset)
                if text:
                    idle = 0.0
                    yield f"event: log\ndata: {json.dumps({'text': text, 'offset': offset})}\n\n"
                else:
                    idle += 0.5

                current = jobs_mod.refresh(job_id)
                if current.is_terminal:
                    # Drain anything written between the last read and exit.
                    tail_text, offset, _ = logs_mod.read_from(job.log_path, offset)
                    if tail_text:
                        yield f"event: log\ndata: {json.dumps({'text': tail_text, 'offset': offset})}\n\n"
                    payload = {"status": current.status, "exit_code": current.exit_code,
                               "note": current.note}
                    yield f"event: status\ndata: {json.dumps(payload)}\n\n"
                    return

                # Keep proxies and browsers from closing an idle stream.
                if idle >= 15.0:
                    idle = 0.0
                    yield ": keepalive\n\n"
                time.sleep(0.5)

        return Response(
            generate(),
            mimetype="text/event-stream",
            headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
        )

    # ---------- results ----------

    @app.get("/api/projects/<pid>/runs")
    def api_project_runs(pid: str):
        project = _project_or_404(pid)
        if not project.available:
            return _err("Project directory is unavailable.", "PROJECT_UNAVAILABLE", 409)
        return jsonify(artifacts_mod.find_runs(project.root))

    @app.get("/api/projects/<pid>/artifacts")
    def api_project_artifacts(pid: str):
        project = _project_or_404(pid)
        if not project.available:
            return _err("Project directory is unavailable.", "PROJECT_UNAVAILABLE", 409)
        run = request.args.get("run", ".")
        try:
            found = artifacts_mod.load_outputs_index(project.root, run)
        except artifacts_mod.PathEscape as exc:
            return _err(str(exc), "PATH_ESCAPE", 403)
        return jsonify([a.to_dict() for a in found])

    @app.get("/api/projects/<pid>/runs/<path:run>/summary")
    def api_run_summary(pid: str, run: str):
        project = _project_or_404(pid)
        if not project.available:
            return _err("Project directory is unavailable.", "PROJECT_UNAVAILABLE", 409)
        try:
            return jsonify(artifacts_mod.load_run_summary(project.root, run))
        except artifacts_mod.PathEscape as exc:
            return _err(str(exc), "PATH_ESCAPE", 403)
        except artifacts_mod.ArtifactError as exc:
            return _err(str(exc), "INVALID_INPUT", 500)

    @app.get("/api/projects/<pid>/file")
    def api_project_file(pid: str):
        project = _project_or_404(pid)
        rel = request.args.get("path", "")
        if not rel:
            return _err("A path is required.", "INVALID_INPUT", 400)
        try:
            absolute = artifacts_mod.safe_resolve(project.root, rel)
        except artifacts_mod.PathEscape as exc:
            return _err(str(exc), "PATH_ESCAPE", 403)
        if not absolute.is_file():
            return _err(f"No such file: {rel}", "NOT_FOUND", 404)
        # as_attachment=False so images and reports render in place.
        return send_file(absolute, as_attachment=False,
                         download_name=absolute.name, max_age=0)

    @app.get("/api/projects/<pid>/table")
    def api_project_table(pid: str):
        project = _project_or_404(pid)
        rel = request.args.get("path", "")
        if not rel:
            return _err("A path is required.", "INVALID_INPUT", 400)
        try:
            absolute = artifacts_mod.safe_resolve(project.root, rel)
        except artifacts_mod.PathEscape as exc:
            return _err(str(exc), "PATH_ESCAPE", 403)
        if not absolute.is_file():
            return _err(f"No such file: {rel}", "NOT_FOUND", 404)
        try:
            page = int(request.args.get("page", "1"))
            page_size = int(request.args.get("page_size", "100"))
        except ValueError:
            return _err("page and page_size must be integers", "INVALID_INPUT", 400)
        try:
            return jsonify(artifacts_mod.read_table(absolute, page, page_size))
        except artifacts_mod.ArtifactError as exc:
            return _err(str(exc), "INVALID_INPUT", 500)

    # ---------- pages ----------

    @app.get("/")
    def page_projects():
        return render_template("projects.html", projects=registry.list_projects())

    @app.get("/project/<pid>")
    def page_dashboard(pid: str):
        project = _project_or_404(pid)
        error = None
        state: dict = {}
        steps: list = []
        if project.available:
            try:
                state = state_adapter.project_state(project.root)
                steps = state_adapter.timeline(project.root)
            except state_adapter.StateUnreadable as exc:
                error = str(exc)
        else:
            error = "This project directory is no longer available."

        return render_template(
            "dashboard.html",
            project=project,
            all_projects=registry.list_projects(),
            state=state,
            steps=steps,
            error=error,
            jobs=jobs_mod.list_jobs(pid)[:15],
        )

    # Page endpoints for later phases are registered now so base.html's nav
    # always resolves. Each is filled in by its own phase.

    @app.get("/project/<pid>/launch")
    def page_launch(pid: str):
        project = _project_or_404(pid)
        catalog = targets_mod.all_targets().values()
        groups: dict[str, list] = {}
        for target in catalog:
            groups.setdefault(target.group, []).append(target)
        return render_template(
            "launch.html",
            project=project,
            all_projects=registry.list_projects(),
            groups=groups,
            has_running=jobs_mod.has_running_job(pid),
        )

    @app.get("/job/<job_id>")
    def page_job(job_id: str):
        try:
            job = jobs_mod.load_job(job_id)
        except jobs_mod.JobError:
            abort(404, description=f"No job with id {job_id}")

        project = None
        try:
            project = registry.get_project(job.project_id)
        except registry.RegistryError:
            pass  # job for an unregistered/removed project still viewable

        return render_template(
            "job.html",
            job=job,
            project=project,
            all_projects=registry.list_projects(),
        )

    @app.get("/project/<pid>/results")
    def page_results(pid: str):
        project = _project_or_404(pid)
        runs = artifacts_mod.find_runs(project.root) if project.available else []
        return render_template(
            "results.html",
            project=project,
            all_projects=registry.list_projects(),
            runs=runs,
        )

    @app.get("/project/<pid>/dag")
    def page_dag(pid: str):
        project = _project_or_404(pid)
        return render_template(
            "placeholder.html",
            project=project,
            all_projects=registry.list_projects(),
            section="DAG",
            phase="Phase 7",
        )

    return app
