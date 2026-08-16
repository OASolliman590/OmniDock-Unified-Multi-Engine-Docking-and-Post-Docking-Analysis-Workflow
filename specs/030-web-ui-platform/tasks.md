# Tasks: OmniDock Web UI Platform

**Input**: Design documents from `/specs/030-web-ui-platform/`
**Prerequisites**: `spec.md`, `plan.md`, `contracts/api.md`, `data-model.md`

## Format: `[ID] [P?] Description`
- **[P]**: parallelizable — touches files no other in-flight task touches
- Every task names the file(s) it changes. Do not touch files outside the named set.
- After each phase, run its verification block before starting the next phase.

## Ground Rules (read before starting any task)

1. **Never reimplement a pipeline stage.** All scientific work runs via `subprocess` calling `python main.py …`.
2. **Import vocabulary, don't copy it.** `ENGINE_CHOICES`, `ANALYSIS_LABELS`, `TIMELINE_STEPS` are imported from `workflow/interactive.py` at runtime. Never hardcode those lists.
3. **Read state through accessors.** Use `workflow/state.py` functions; never `json.load(".workflow/state.json")` directly.
4. **The web layer never writes into a project directory.** Only launched subprocesses do.
5. **Adapters stay framework-free.** Only `webui/app.py` imports Flask.

---

## Phase 1: Scaffold

- [x] T001 Add `flask>=3.0` to `requirements.txt`.
- [x] T002 Create `webui/__init__.py` with `__version__ = "0.1.0"`.
- [x] T003 Create `webui/config.py` exposing `APP_HOME` (default `~/.omnidock_webui`, overridable via `OMNIDOCK_WEBUI_HOME`), `PROJECTS_FILE`, `JOBS_DIR`, `LOGS_DIR`, `HOST = "127.0.0.1"`, `PORT = 8770`; add `ensure_app_dirs()` that creates them.
- [x] T004 Create `webui/app.py` with a `create_app()` factory returning a Flask app, calling `ensure_app_dirs()` at startup.
- [x] T005 Add `GET /api/health` returning `{"ok": true, "version": <webui version>}`.
- [x] T006 Create `webui_cli.py` at repo root that runs `create_app().run(host=HOST, port=PORT)`; support `--host`/`--port` overrides and a `--allow-remote` flag that is the *only* way to bind off-localhost (FR-034).
- [x] T007 Register a `webui` subcommand in `workflow/cli.py` that invokes the same entry point, so `python main.py webui` works.
- [x] T008 Create `webui/templates/base.html`: page shell with a header, a project switcher placeholder, a nav bar, and a `{% block content %}`.
- [x] T009 Create `webui/static/app.css` with a minimal readable stylesheet (no framework, no CDN).

**Verify Phase 1**: `python webui_cli.py` starts; `curl 127.0.0.1:8770/api/health` returns ok; `python main.py webui --help` works.

---

## Phase 2: Project Registry & Dashboard (US1, P1)

- [x] T010 Create `webui/registry.py` with a `Project` dataclass (`id`, `path`, `name`, `registered_at`).
- [x] T011 Implement `canonicalize(path)` in `registry.py`: expand `~`, resolve to absolute, resolve symlinks, and normalize case on Windows — so one directory cannot register twice (Edge Case: duplicate spellings).
- [x] T012 Implement `register(path, name=None)`: validate the directory exists and is readable, generate a stable id (hash of canonical path), persist to `projects.json`; return existing entry if already registered.
- [x] T013 Implement `list_projects()`, `get_project(id)`, `unregister(id)` with atomic writes to `projects.json`.
- [x] T014 Add `available` computation to project listing: `False` when the directory no longer exists (FR-030) — never raise.
- [x] T015 Create `webui/state_adapter.py` importing `ensure_state`, `load_state`, `list_steps`, `summarize_state`, `list_background_tasks` from `workflow.state`.
- [x] T016 Implement `project_state(project)` in `state_adapter.py` returning a plain dict: `current_context`, `feature_flags`, `artifacts`, plus `has_workflow: bool` (False when `.workflow/` is absent).
- [x] T017 Implement `timeline(project)` in `state_adapter.py`: import `TIMELINE_STEPS` from `workflow.interactive`, join with `list_steps()`, return ordered `[{key, label, status, timestamp, note}]` (FR-003).
- [x] T018 Make all `state_adapter` reads tolerant of a transient partial `state.json`: retry once after a short delay, then surface a typed error (Edge Case: mid-write read).
- [x] T019 Add `GET /api/projects` and `POST /api/projects` (body `{path, name?}`) to `app.py`, delegating to `registry`.
- [x] T020 Add `GET /api/projects/<pid>/state` and `GET /api/projects/<pid>/timeline`.
- [x] T021 Add `POST /api/projects/<pid>/init` that shells out to `python main.py workflow init --project-dir <path>` (FR-004).
- [x] T022 Create `webui/templates/dashboard.html` rendering the timeline as an ordered list with per-step status badges, plus a context panel (favorite engine, selected engines, docking root, analysis output dir).
- [x] T023 Create `webui/templates/projects.html` with a registration form and the registered-project list; mark unavailable ones visibly.
- [x] T024 Wire routes `GET /` (projects page) and `GET /project/<pid>` (dashboard).
- [x] T025 Ensure no route caches state between requests (FR-005) — read fresh on every call; add a comment stating why.

**Verify Phase 2**: register a fixture project; dashboard timeline matches `.workflow/state.json`; register the same path with a trailing slash and assert one entry; delete the directory and assert the UI shows it unavailable.

---

## Phase 3: Durable Job Runner (US3, P1) — build before launch UI

- [x] T026 Create `webui/jobs.py` with a `Job` dataclass: `job_id`, `project_id`, `target`, `argv` (list), `cwd`, `log_path`, `pid`, `proc_start_time`, `status`, `exit_code`, `created_at`, `started_at`, `finished_at`, `note`.
- [x] T027 Define `TERMINAL_STATUSES = {"completed", "failed", "cancelled"}` and `VALID_STATUSES` including `queued`, `running` (FR-017).
- [x] T028 Implement `save_job(job)` / `load_job(job_id)` / `list_jobs(project_id=None)` reading and writing `jobs/<job_id>.json` atomically.
- [x] T029 Implement `create_job(project, target, argv)`: persist the record with status `queued` **before** any process starts (FR-012).
- [x] T030 Implement `launch(job)`: open `logs/<job_id>.log` for append, `subprocess.Popen(argv, cwd=project.path, stdout=log, stderr=STDOUT, stdin=DEVNULL)`, record `pid` and `proc_start_time`, set status `running`, persist.
- [x] T031 In `launch()`, set a non-interactive environment: copy `os.environ`, set `OMNIDOCK_NON_INTERACTIVE=1` and `PYTHONUNBUFFERED=1`, and pass `stdin=DEVNULL` so a stray prompt fails instead of hanging (FR-019).
- [x] T032 Write an exit-marker file `jobs/<job_id>.exit` containing the return code when the process finishes, via a small watcher thread — so a restarted server can recover the code without the original `Popen` object.
- [x] T033 Implement `_process_alive(pid)` using `psutil` when importable, else an OS-level fallback (`os.kill(pid, 0)` on POSIX, `OpenProcess` query on Windows).
- [x] T034 Implement `_identity_matches(job)`: compare recorded `proc_start_time` against the live process's start time; return `unverifiable` (distinct from True/False) when no probe is available (FR-014).
- [x] T035 Implement `reconcile_all()` applying the four-case rule from `plan.md` (no pid → failed; not alive → exit marker or failed; alive+identity match → running; alive+mismatch → failed "pid recycled").
- [x] T036 Call `reconcile_all()` once during `create_app()` startup (FR-013).
- [x] T037 Implement `cancel(job_id)`: terminate the process (SIGTERM, then SIGKILL after a grace period), set status `cancelled`, never delete produced artifacts (FR-016).
- [x] T038 Implement `has_running_job(project_id)` used later to warn on a second concurrent job (FR-018).
- [x] T039 Add `GET /api/jobs?project_id=`, `GET /api/jobs/<job_id>`, `POST /api/jobs/<job_id>/cancel`.
- [x] T040 Add `test/test_webui.py` with a restart-durability test: launch a long fixture job, simulate a server restart by calling `reconcile_all()` in a fresh module state, assert the status stays `running`, then let it exit and assert `completed` (SC-003).
- [x] T041 [P] Add a reconciliation test for the recycled-PID case: fabricate a record with a live PID but a mismatched `proc_start_time`; assert `failed` (FR-014).

**Verify Phase 3**: all Phase-3 tests pass; a job launched from a Python shell survives the shell exiting.

---

## Phase 4: Launch Targets & Forms (US2, P1)

- [ ] T042 Create `webui/targets.py` importing `ANALYSIS_LABELS` and `ENGINE_CHOICES` from `workflow.interactive` (FR-007).
- [ ] T043 Define a `FormField` dataclass (`name`, `label`, `type` in {`text`,`path`,`select`,`multiselect`,`bool`,`int`}, `required`, `default`, `choices`) and a `LaunchTarget` dataclass (`key`, `label`, `group`, `fields`, `argv_template`).
- [ ] T044 Build the analysis-target catalog by iterating `ANALYSIS_LABELS` — do not enumerate keys by hand (SC-002).
- [ ] T045 Add pipeline targets not in `ANALYSIS_LABELS`: `pdb.collect`, `pdb.prepare_both`, `prep.pairlist`, `prep.project`, `dock.run`, `dock.dry-run`, mapping each to its `workflow/cli.py` subcommand.
- [ ] T046 Implement `build_argv(target, form_values, project)` producing the full argv list (`["python", "main.py", …]`) with `--project-dir` always set to the project path.
- [ ] T047 Implement server-side validation `validate(target, form_values)`: required fields present, paths exist, engines within `ENGINE_CHOICES`, unknown target rejected (FR-009).
- [ ] T048 Mark the three interaction aliases (`analyze.interactions.prolif`, `.ligplot`, `.pandamap`) with a `routes_to: "analyze.interactions.clean"` note so the UI can state the routing (FR-011, spec 029 FR-002).
- [ ] T049 Pre-fill form defaults from the project's `current_context` (favorite engine, selected engines, output dir) in the form-schema endpoint.
- [ ] T050 Add `GET /api/targets` (catalog) and `GET /api/targets/<key>/form?project_id=` (schema with defaults applied) (FR-006).
- [ ] T051 Add `POST /api/projects/<pid>/jobs` — validate, build argv, `create_job`, `launch`, return the job record; include the resolved command string in the response (FR-008).
- [ ] T052 Return HTTP 409 with a warning payload when `has_running_job(project_id)` is true, unless the request sets `confirm_concurrent: true` (FR-018).
- [ ] T053 Create `webui/templates/launch.html`: target picker grouped by `group`, dynamic form rendered from the schema, and a live **command preview** box showing the exact argv before launch.
- [ ] T054 Show the `routes_to` note in the form when present (FR-011).
- [ ] T055 [P] Add tests: every `ANALYSIS_LABELS` key yields a form schema; a monkeypatched extra label appears in `/api/targets` without web-side changes (SC-002).
- [ ] T056 [P] Add a validation test matrix: missing required field, non-existent path, unknown engine, unknown target — each rejected with a field-level message and no job created.

**Verify Phase 4**: launching `analyze.comparative` from the browser against a fixture project creates a job whose argv matches the documented CLI invocation.

---

## Phase 5: Live Log Streaming & Job UI (US3, P1)

- [ ] T057 Create `webui/logs.py` with `read_from(log_path, offset, max_bytes)` returning `(text, new_offset, size)` — incremental byte-offset reads, never loading the whole file (FR-015).
- [ ] T058 Implement `tail(log_path, max_bytes)` returning only the last N bytes for initial page load (Edge Case: very large logs).
- [ ] T059 Add `GET /api/jobs/<job_id>/log?offset=` returning `{text, offset, size, status}`.
- [ ] T060 Add `GET /api/jobs/<job_id>/stream` as an SSE endpoint emitting new log chunks and terminal-status events; close the stream when the job reaches a terminal status.
- [ ] T061 Create `webui/templates/job.html`: status header (target, status badge, exit code, timestamps), the resolved command, a log pane, and a Cancel button.
- [ ] T062 Add `webui/static/logstream.js`: connect to the SSE endpoint, append chunks, auto-scroll unless the user has scrolled up, and fall back to polling `/log?offset=` if `EventSource` errors.
- [ ] T063 Cap the in-DOM log buffer (for example 5000 lines, dropping from the top) so long runs cannot freeze the tab.
- [ ] T064 Add a jobs list view to the dashboard showing this project's jobs with status and start time.
- [ ] T065 [P] Add a test asserting incremental reads: write to a log file in two steps and assert the second read returns only the new bytes.

**Verify Phase 5**: launch a job that prints steadily; log lines appear in the browser within ~2s (SC-006); Cancel transitions it to `cancelled`.

---

## Phase 6: Results Browser (US4, P1)

- [ ] T066 Create `webui/artifacts.py` with `safe_resolve(project_root, rel_path)`: resolve the joined path, resolve symlinks, and assert the result is inside `project_root`; raise `PathEscape` otherwise (FR-022).
- [ ] T067 Implement `load_outputs_index(project, run_dir)` reading `run_tracking/outputs_index.json`, falling back to a directory walk when absent (FR-020).
- [ ] T068 Implement `categorize(rel_path)` mapping into the canonical topology buckets: `analysis`, `complexes`, `best_poses`, `reports`, `rmsd_analysis`, `interactions/{prolif,ligplot,pandamap,poseview}`, `3d_visualizations`, `visualizations`, `raw_data`, `other`.
- [ ] T069 Implement `media_kind(path)` → `table` (.csv/.tsv), `image` (.png/.svg/.jpg), `html` (.html), `structure` (.pdb/.sdf/.pdbqt), `other`.
- [ ] T070 Implement `load_run_summary(run_dir)` parsing `run_manifest.json`: engine, `stage_contract` requested-vs-enforced, per-step status, optional-feature classifications (FR-023, FR-024).
- [ ] T071 Add `GET /api/projects/<pid>/runs` listing discovered run directories containing `run_tracking/`.
- [ ] T072 Add `GET /api/projects/<pid>/artifacts?run=` returning categorized artifact entries.
- [ ] T073 Add `GET /api/projects/<pid>/file?path=` serving a file through `safe_resolve`, returning 403 on escape and 404 on missing.
- [ ] T074 Add `GET /api/projects/<pid>/table?path=` returning parsed CSV rows with pagination (`page`, `page_size`) using pandas.
- [ ] T075 Create `webui/templates/results.html`: run selector, category accordion, and a viewer pane.
- [ ] T076 Add `webui/static/results.js`: render tables (sortable, paginated), images inline, HTML in an iframe sandboxed to the served origin (FR-021).
- [ ] T077 Render the run summary panel with `enforced != requested` rows visually flagged (FR-023).
- [ ] T078 Render optional-feature classifications with distinct treatment for `skipped_disabled` vs `skipped_missing_dependency` vs `failed_error` (FR-024).
- [ ] T079 [P] Add the path-traversal test matrix: `../`, `..\\`, URL-encoded separators, an absolute outside path, and a symlink escaping the root — all must return 403 (SC-005).
- [ ] T080 [P] Add a test that every entry in a fixture `outputs_index.json` appears in the artifacts response.

**Verify Phase 6**: a completed fixture run renders tables, images, and its report; all traversal attempts return 403.

---

## Phase 7: Artifact DAG View (US5, P2)

- [ ] T081 Create `webui/dag.py` importing `NODE_STATUSES` from `post_docking_analysis.artifact_graph` (FR-025) — never hardcode the status strings.
- [ ] T082 Implement `load_dag_report(project, run_dir)` reading the report written by `ArtifactGraph._write_report`, plus `4-Working/dag_cache.json` for cache indicators.
- [ ] T083 Implement `to_view_model(report)` producing `{nodes: [{name, status, inputs, outputs, cached}], edges: [{from, to}]}` with edges derived from input/output relationships.
- [ ] T084 Add `GET /api/projects/<pid>/dag?run=`.
- [ ] T085 Create `webui/templates/dag.html` rendering nodes in dependency tiers with one distinct visual treatment per status (SC-008).
- [ ] T086 Show cache hits distinctly from recomputed nodes (FR-027).
- [ ] T087 Link each node's outputs into the results browser (FR-026).
- [ ] T088 [P] Add a fixture DAG report exercising all eight statuses and a test asserting each renders distinctly.

**Verify Phase 7**: the fixture DAG shows all eight statuses; a failed required node's dependents display as `blocked_by_failure`.

---

## Phase 8: Multi-Project & Concurrency (US6, P2)

- [ ] T089 Add a project switcher to `base.html` populated from `/api/projects`, preserving the current page type when switching.
- [ ] T090 Scope every project-specific route and template strictly by `project_id`; add a test that two projects' job lists never intersect (SC-007).
- [ ] T091 Surface the concurrent-job warning from T052 in the launch UI with an explicit confirm control.
- [ ] T092 [P] Add a concurrency test: launch a job in each of two projects and assert both progress independently.

**Verify Phase 8**: two projects run jobs simultaneously with independent logs and no state bleed.

---

## Phase 9: HPC Panel (US7, P3)

- [ ] T093 Create `webui/hpc.py` with `list_profiles(project)` reading `examples/hpc_profiles/` and `<project>/.workflow/hpc_profiles/` (FR-031).
- [ ] T094 Implement `redact(profile)` stripping accounts, home paths, and any credential-shaped values before the profile is returned to the browser (FR-032).
- [ ] T095 Add `GET /api/projects/<pid>/hpc/profiles` returning redacted profiles only.
- [ ] T096 Add launch targets `dock.deploy`, `dock.sync`, `dock.submit` routed through the same `jobs` runner.
- [ ] T097 Require a confirmation payload naming the remote target for `sync` and `submit`; reject without it (FR-033).
- [ ] T098 Create `webui/templates/hpc.html` with profile selection, a deploy/sync/submit panel, and the confirmation step.
- [ ] T099 [P] Add a secret-hygiene test: a fixture profile with sentinel secrets must not appear in rendered HTML or job logs (SC-009).

**Verify Phase 9**: deploy with the public-safe template produces the same assets as the CLI; sentinels never leak.

---

## Phase 10: Hardening, Tests, Docs

- [ ] T100 Add a static import-discipline test walking `webui/*.py` ASTs and failing on any import from `post_docking_analysis` other than `artifact_graph`, or any import of `workflow.execution` (SC-004).
- [ ] T101 Add a test asserting `create_app().run` defaults bind to `127.0.0.1` and that off-localhost requires `--allow-remote` (FR-034).
- [ ] T102 Add a test that the web layer writes nothing into a project directory: snapshot a fixture project tree, exercise every read-only route, assert the tree is unchanged (FR-035).
- [ ] T103 Add a non-interactive fail-fast test: a target that would prompt exits with a clear error rather than hanging (FR-019).
- [ ] T104 Create `test/fixtures/webui_project/` with `.workflow/state.json`, a `run_tracking/` set (manifest + step_status + outputs_index), and a DAG report covering all statuses.
- [ ] T105 Write `docs/WEB_UI_GUIDE.md`: install, launch, register a project, run an analysis, read results, and the explicit security posture (localhost-only, no auth).
- [ ] T106 Update `README.md` with a Web UI section pointing to the guide.
- [ ] T107 Update `CHANGELOG.md` with the feature entry.
- [ ] T108 Add the web test module to the existing smoke harness registration in `test/test_dockforge_smoke.py`.
- [ ] T109 Final pass: run `python -m py_compile` across `webui/`, run the full `pytest`, and walk every command in `docs/WEB_UI_GUIDE.md` against a fresh checkout.

**Verify Phase 10**: full suite green; docs commands all execute; SC-001 through SC-009 each have a passing verification.

---

## Suggested Commit Boundaries

One commit per phase, each independently reviewable:

1. `feat(webui): scaffold flask app and entry point`
2. `feat(webui): project registry and state dashboard`
3. `feat(webui): durable job runner with restart reconciliation`
4. `feat(webui): launch targets and forms`
5. `feat(webui): live log streaming and job views`
6. `feat(webui): results browser with safe artifact serving`
7. `feat(webui): artifact DAG visualization`
8. `feat(webui): multi-project switching and concurrency guards`
9. `feat(webui): HPC deployment panel`
10. `test(webui): hardening, fixtures, and documentation`
