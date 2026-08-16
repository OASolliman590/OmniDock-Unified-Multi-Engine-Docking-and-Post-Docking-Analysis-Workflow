"""Flask application factory and route registration.

Routes only -- all behavior lives in the adapter modules. This is the one
module in the package permitted to import Flask.

State is read fresh on every request. The web UI must never serve a cached
copy, because the CLI and the interactive shell write the same project
state concurrently from other processes (FR-005).
"""

from __future__ import annotations

import sys

from flask import Flask, abort, jsonify, render_template, request

from . import __version__
from . import jobs as jobs_mod
from . import registry
from . import state_adapter
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

    # ---------- jobs ----------

    @app.get("/api/jobs")
    def api_jobs_list():
        project_id = request.args.get("project_id")
        return jsonify([j.to_dict() for j in jobs_mod.list_jobs(project_id)])

    @app.get("/api/jobs/<job_id>")
    def api_job_detail(job_id: str):
        try:
            return jsonify(jobs_mod.load_job(job_id).to_dict())
        except jobs_mod.JobError as exc:
            return _err(str(exc), "NOT_FOUND", 404)

    @app.post("/api/jobs/<job_id>/cancel")
    def api_job_cancel(job_id: str):
        try:
            return jsonify(jobs_mod.cancel(job_id).to_dict())
        except jobs_mod.JobError as exc:
            return _err(str(exc), "NOT_FOUND", 404)

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
        )

    # Page endpoints for later phases are registered now so base.html's nav
    # always resolves. Each is filled in by its own phase.

    @app.get("/project/<pid>/launch")
    def page_launch(pid: str):
        project = _project_or_404(pid)
        return render_template(
            "placeholder.html",
            project=project,
            all_projects=registry.list_projects(),
            section="Launch",
            phase="Phase 4",
        )

    @app.get("/project/<pid>/results")
    def page_results(pid: str):
        project = _project_or_404(pid)
        return render_template(
            "placeholder.html",
            project=project,
            all_projects=registry.list_projects(),
            section="Results",
            phase="Phase 6",
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
