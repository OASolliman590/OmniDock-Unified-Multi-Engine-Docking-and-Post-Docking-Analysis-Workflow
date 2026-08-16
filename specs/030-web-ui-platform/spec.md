# Feature Specification: OmniDock Web UI Platform

**Feature Branch**: `030-web-ui-platform`
**Created**: 2026-08-16
**Status**: Draft
**Input**: User description: "thorough analysis of omnidock repo and plan how it can be turned into this interactive version or web app" — turn the questionary-driven terminal workflow into a browser-based interactive platform.

## Context

OmniDock today is a terminal application. `main.py` normalizes legacy argv and delegates to `workflow/cli.py:main()`, which exposes five command groups (`workflow`, `pdb`, `prep`, `dock`, `analyze`). The guided experience lives in `workflow/interactive.py` (~132 KB), which drives `questionary` prompts and calls into `workflow/execution.py` (~123 KB) to run targets.

Three existing subsystems already provide almost everything a web UI needs, and this feature is primarily about **exposing them over HTTP rather than rebuilding them**:

1. **Durable project state** — `workflow/state.py` maintains `.workflow/state.json` (schema `version: 2`) with `current_context`, `artifacts`, `steps`, `background_tasks`, `feature_flags`, and `checkpoint_metadata`. Writes are atomic and lock-guarded (`_atomic_dump_json`, `_STATE_IO_LOCK`), so an external reader is safe.
2. **Machine-readable run tracking** — `post_docking_analysis/simplified_pipeline_impl.py` writes `run_tracking/run_manifest.json`, `step_status.json`, and `outputs_index.{csv,json}`. The manifest already carries `stage_contract` (requested vs enforced) and per-optional-feature classification (`completed` / `skipped_disabled` / `skipped_missing_dependency` / `failed_error`, spec 029 FR-008).
3. **Artifact DAG with status and caching** — `post_docking_analysis/artifact_graph.py` executes tier-grouped nodes through a `ThreadPoolExecutor` and records per-node status from the set `pending | running | completed | completed_with_warnings | failed | skipped_optional | blocked_by_failure | cache_hit`, with a hash-based cache in `4-Working/dag_cache.json` (spec 024).

The gaps that make a browser UI non-trivial:

- **`BackgroundTaskManager` (`workflow/execution.py:132`) is in-process only.** It wraps a `ThreadPoolExecutor(max_workers=2)`; if the host process exits, running work is lost and `background_tasks` entries in `state.json` are stranded in `running`. Docking runs last hours, so a web server process cannot own them.
- **The interactive surface is not reusable.** `workflow/interactive.py` interleaves prompting, validation, and dispatch; its prompt logic cannot be called headlessly. The web UI must target `workflow/execution.py` / the CLI, not `interactive.py`.
- **No HTTP surface exists.** No web framework is present in `requirements.txt`. Spec 018 FR-070 reserved the idea ("Architecture MUST expose extension points for dashboard-driven project inspection") but nothing implements it. This feature fulfills FR-070.
- **Artifact paths are user-supplied.** Serving files from a project root over HTTP introduces path-traversal risk that the terminal app never had.

**Intended outcome:** a locally-run web application that lets a researcher register a project directory, inspect its workflow state and timeline, launch preparation/docking/analysis runs through forms instead of questionary prompts, watch those runs stream progress live, and browse the resulting artifacts, reports, and consensus tables in the browser — with several projects and several runs in flight concurrently, and with no duplication of pipeline logic.

## Non-Goals

- **Not a multi-tenant or public web service.** Binds to `127.0.0.1` only; there is no authentication, no user accounts, and no network exposure in this feature.
- **Not a replacement for the CLI or the interactive shell.** Both remain fully supported; the web UI is an additional surface over the same execution layer.
- **Not a reimplementation of any pipeline stage.** Every scientific computation continues to run through existing `workflow/execution.py` targets and CLI commands.
- **Not a molecular structure editor.** 3D pose display is read-only rendering of existing artifacts.

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Register a project and see its true state (Priority: P1)

A researcher starts the web app, points it at an existing OmniDock project directory, and immediately sees which timeline steps have completed, which engines are configured, what the favorite engine is, and where the analysis outputs live — without reading `state.json` by hand or re-running `workflow status`.

**Why this priority**: Everything else needs a project context. This story alone replaces `python main.py workflow status` with a readable view and proves the state-adapter layer works.

**Independent Test**: Register a fixture project directory through the UI; assert the rendered timeline matches `TIMELINE_STEPS` (`workflow/interactive.py:86`) with per-step status drawn from `list_steps()`, and that `current_context` fields (favorite engine, selected engines, docking root) match `.workflow/state.json` byte-for-byte.

**Acceptance Scenarios**:
1. **Given** a directory that already contains `.workflow/state.json`, **When** the user registers it, **Then** the dashboard lists every timeline step with its recorded status and timestamps.
2. **Given** a directory with no `.workflow/`, **When** the user registers it, **Then** the UI offers to initialize it (equivalent to `workflow init`) and does not crash or write partial state.
3. **Given** a path that does not exist or is not readable, **When** the user submits it, **Then** the UI shows a specific error and registers nothing.
4. **Given** a registered project whose `state.json` is modified externally (by a CLI run in another terminal), **When** the user refreshes the dashboard, **Then** the new state is shown — the UI holds no stale cached copy.

---

### User Story 2 — Launch a run from a form instead of questionary prompts (Priority: P1)

A researcher selects an analysis target (for example `analyze.comparative`), fills in a form that presents the same choices the terminal wizard would have asked for — engines, favorite engine, output directory, feature toggles — and starts the run. The equivalent CLI command is shown before launch so the choice is auditable and reproducible.

**Why this priority**: This is the actual "interactive version" the feature exists to deliver. Without it the app is a read-only viewer.

**Independent Test**: For each target in `ANALYSIS_LABELS` (`workflow/interactive.py:68`), request its form schema and assert every required field is declared; submit a valid payload against a fixture project and assert the launched command matches the documented CLI invocation for that target.

**Acceptance Scenarios**:
1. **Given** a registered project, **When** the user opens the launch page, **Then** every target from `ANALYSIS_LABELS` plus the pipeline commands (`pdb.collect`, `prep.pairlist`, `prep.project`, `dock.run`) are selectable with human-readable labels.
2. **Given** a selected target, **When** the form renders, **Then** engine choices come from `ENGINE_CHOICES` (`workflow/interactive.py:60`) and defaults are pre-filled from the project's `current_context`.
3. **Given** a completed form, **When** the user clicks Launch, **Then** the exact CLI command is displayed and a job record is created before the subprocess starts.
4. **Given** an invalid form (missing required field, non-existent path), **When** the user submits, **Then** validation fails in the UI with a field-level message and no job is created.
5. **Given** a target routed to the clean interactions contract (`analyze.interactions.prolif` / `.ligplot` / `.pandamap`, spec 029 FR-002), **When** the user selects it, **Then** the UI states that it routes to `analyze.interactions.clean` rather than silently substituting.

---

### User Story 3 — Watch a long run without holding the browser or the server hostage (Priority: P1)

A researcher launches a multi-hour docking run, watches live log output and step progress in the browser, closes the laptop lid, and returns later — from a fresh browser session, or after restarting the web server — to find the run still progressing and its history intact.

**Why this priority**: Docking runs outlive both HTTP requests and server processes. If job durability is not designed in from the start, every later story inherits a broken execution model. This is the single highest-risk piece of the feature.

**Independent Test**: Launch a job; kill the web server process; restart it; assert the job list still shows the job, that its status is reconciled from disk (still `running` if the OS process is alive, `failed` with a reason if it died), and that the log view continues from the persisted log file.

**Acceptance Scenarios**:
1. **Given** a launched job, **When** the user watches its detail page, **Then** new log lines appear without a manual refresh.
2. **Given** a running job, **When** the web server is restarted, **Then** the job record survives and its live status is recovered from the persisted record plus an OS liveness check.
3. **Given** a job whose process exited non-zero, **When** the user views it, **Then** the status is `failed`, the exit code is shown, and the log tail containing the error is displayed.
4. **Given** a running job, **When** the user requests cancellation, **Then** the process is terminated, the status becomes `cancelled`, and partial artifacts are left in place rather than deleted.
5. **Given** two jobs launched against different projects, **When** both run, **Then** neither blocks the other and each writes to its own log file.
6. **Given** a job for a project, **When** it completes, **Then** the project's timeline reflects the new step status on the next dashboard load.

---

### User Story 4 — Browse results, tables, and reports in the browser (Priority: P1)

After an analysis completes, the researcher explores its outputs in the UI: consensus and ranking tables rendered as sortable tables, generated plots and interaction diagrams shown inline, and HTML reports opened in place — without hunting through the numbered output topology in a file manager.

**Why this priority**: Producing results the user still has to find on disk only half-solves the problem. Reading `outputs_index.json` makes this generic rather than per-stage bespoke.

**Independent Test**: Against a fixture project with a completed analysis, assert the results view lists every entry from `run_tracking/outputs_index.json`, that each CSV renders as a table, and that requesting a path outside the project root returns HTTP 403.

**Acceptance Scenarios**:
1. **Given** a completed run, **When** the user opens results, **Then** artifacts are grouped by the canonical output topology (`analysis/`, `complexes/`, `best_poses/`, `reports/`, `rmsd_analysis/`, `interactions/{prolif,ligplot,pandamap,poseview}/`, `3d_visualizations/`, `visualizations/`, `raw_data/`).
2. **Given** a CSV artifact, **When** the user opens it, **Then** it renders as a paginated, sortable table rather than raw text.
3. **Given** a PNG or SVG artifact, **When** the user opens it, **Then** it displays inline at readable size.
4. **Given** an HTML report, **When** the user opens it, **Then** it renders in the results pane.
5. **Given** a crafted request containing `../` or an absolute path outside the registered project root, **When** the server handles it, **Then** it returns 403 and serves nothing.
6. **Given** a run whose `stage_contract` shows `enforced != requested`, **When** the user views the run summary, **Then** those differences are called out explicitly (spec 029 FR-007 parity).

---

### User Story 5 — See the artifact DAG and understand what ran, cached, or failed (Priority: P2)

A researcher views the artifact graph for a run: which nodes completed, which were cache hits, which were skipped as optional, and which failed and blocked their dependents — turning the DAG that spec 024 already executes into something legible.

**Why this priority**: High explanatory value and moderate cost, since `artifact_graph.py` already records exactly this. It is P2 because the run is still usable through the results browser without it.

**Independent Test**: Run the DAG against a fixture project with one deliberately failing optional node; assert the rendered graph shows that node as `skipped_optional` and its dependents as not `blocked_by_failure`.

**Acceptance Scenarios**:
1. **Given** a completed DAG run, **When** the user opens the graph view, **Then** every node shows one of the eight statuses from `NODE_STATUSES` (`artifact_graph.py:21`) with a distinct visual treatment.
2. **Given** a node that was a cache hit, **When** the user inspects it, **Then** the UI indicates it was served from `4-Working/dag_cache.json` and did not recompute.
3. **Given** a failed required node, **When** the user inspects downstream nodes, **Then** they are shown as `blocked_by_failure` with a pointer to the failing node.
4. **Given** a node with outputs, **When** the user clicks it, **Then** its output artifacts link into the results browser.

---

### User Story 6 — Work on several projects and runs at once (Priority: P2)

A researcher keeps multiple projects registered and several runs in flight, switching between them via tabs without losing scroll position, form state, or live log connections.

**Why this priority**: This is how the work is actually done (the MRSA, K. pneumoniae, and Acinetobacter projects proceed in parallel), but a single-project UI is still useful, so it follows the core loop.

**Independent Test**: Register two projects, launch one job in each, and assert both job detail views stream independently and that the project switcher never mixes state between them.

**Acceptance Scenarios**:
1. **Given** several registered projects, **When** the user switches projects, **Then** the dashboard, jobs, and results all scope to the selected project.
2. **Given** runs in flight in two projects, **When** the user switches between them, **Then** both continue running and each log view resumes correctly.
3. **Given** a registered project directory that has been deleted from disk, **When** the UI loads, **Then** it is marked unavailable rather than causing an error page.

---

### User Story 7 — Prepare and submit HPC work from the browser (Priority: P3)

A researcher generates Slurm deployment assets, syncs a prepared project to the cluster, and submits it, using the same `dock deploy` / `dock sync` / `dock submit` commands the CLI exposes, with profile selection surfaced as a form.

**Why this priority**: Valuable but the highest-risk surface (it touches remote systems and credentials), and the CLI path already works. It should land only after the local loop is proven.

**Independent Test**: With a public-safe profile (`examples/hpc_profiles/bibalex-apptainer-conda.template.json`), run deploy in a dry mode and assert the generated asset list matches the CLI's output for the same inputs.

**Acceptance Scenarios**:
1. **Given** a prepared project, **When** the user opens the HPC panel, **Then** available profiles are listed from `examples/hpc_profiles/` and `.workflow/hpc_profiles/`.
2. **Given** a selected profile, **When** the user generates deployment assets, **Then** the same files are produced as the equivalent `dock deploy` invocation.
3. **Given** any HPC action, **When** the UI displays a profile, **Then** no secret values (accounts, home paths, credentials) are rendered in the browser or written to the job log.
4. **Given** a sync or submit action, **When** the user triggers it, **Then** a confirmation step naming the remote target is required before execution.

---

### Edge Cases

- **Project registered twice under different path spellings** (trailing slash, different case on Windows, symlink): must resolve to one canonical entry, not two.
- **`state.json` mid-write when the UI reads it**: reader must tolerate a transient partial file and retry rather than surfacing a JSON parse error (writer is atomic, but the backup path `_backup_state_path` exists for a reason).
- **Job log grows very large** (multi-hour docking): the log view must tail rather than load the whole file into memory or the DOM.
- **Two jobs writing to the same project simultaneously**: the UI must warn when a second job targets a project that already has a running job, since pipeline stages are not guaranteed concurrency-safe within one project.
- **Orphaned `running` entries in `background_tasks`** left by a previously killed interactive session: reconciliation must not treat these as live web jobs.
- **Server restarted while a job's PID has been reused** by an unrelated OS process: liveness check must verify identity (recorded start time or command line), not PID existence alone.
- **Artifact file replaced while being served**: partial reads must not be cached and presented as complete.
- **A target that requires a TTY** (any residual `questionary` prompt reached through a non-interactive path): must fail fast with a clear message rather than hanging forever waiting on stdin.
- **Non-ASCII or space-containing project paths** on Windows: subprocess launch must quote correctly.
- **Disk full during a run**: job must record a `failed` status with the underlying error rather than appearing to hang.

## Requirements *(mandatory)*

### Functional Requirements

**Project registry and state (P1):**
- **FR-001**: System MUST let a user register an OmniDock project by absolute filesystem path, canonicalize it, and persist the registry across server restarts.
- **FR-002**: System MUST read project state exclusively through `workflow/state.py` accessors (`ensure_state`, `load_state`, `list_steps`, `summarize_state`, `list_background_tasks`) rather than parsing `.workflow/state.json` directly.
- **FR-003**: System MUST render the timeline using `TIMELINE_STEPS` and per-step status, showing status, timestamp, and any recorded note.
- **FR-004**: System MUST offer project initialization (equivalent to `workflow init`) for a directory that has no `.workflow/` yet.
- **FR-005**: System MUST re-read state on each request and MUST NOT serve a cached copy that could diverge from concurrent CLI usage.

**Run launching (P1):**
- **FR-006**: System MUST expose a machine-readable form schema per launchable target, declaring each field's name, type, required-ness, default, and allowed values.
- **FR-007**: System MUST source engine choices from `ENGINE_CHOICES` and target labels from `ANALYSIS_LABELS` so the web UI cannot drift from the terminal UI.
- **FR-008**: System MUST display the exact CLI command that will be executed before launching it, and MUST store that command on the job record.
- **FR-009**: System MUST validate all form input server-side before creating a job, rejecting non-existent paths and unknown targets.
- **FR-010**: System MUST execute runs by invoking OmniDock's existing CLI (`python main.py …`) as a subprocess, and MUST NOT import pipeline internals to re-run stages in the web process.
- **FR-011**: System MUST state, in the UI, when a selected target routes to a different canonical target (interactions aliases → `analyze.interactions.clean`).

**Job durability and monitoring (P1):**
- **FR-012**: System MUST persist every job record to disk at creation, including job id, project id, command, working directory, log path, PID, process start time, and status.
- **FR-013**: System MUST survive a web-server restart with all job records intact and MUST reconcile each non-terminal job's status against OS process liveness on startup.
- **FR-014**: System MUST verify process identity (recorded start time and/or command line), not PID existence alone, when reconciling.
- **FR-015**: System MUST stream new log output to the browser without full-page refresh, reading incrementally from the log file by byte offset.
- **FR-016**: System MUST support cancelling a running job by terminating its process and recording status `cancelled`, leaving produced artifacts in place.
- **FR-017**: System MUST record terminal status as one of `completed` (exit 0), `failed` (non-zero exit, with code), or `cancelled`.
- **FR-018**: System MUST warn before launching a second job against a project that already has a running job.
- **FR-019**: System MUST run jobs with a non-interactive environment so that any residual TTY prompt fails fast rather than blocking.

**Results browsing (P1):**
- **FR-020**: System MUST list run artifacts from `run_tracking/outputs_index.json` when present, falling back to a directory walk of the canonical output topology.
- **FR-021**: System MUST render CSV/TSV artifacts as paginated sortable tables, images (PNG/SVG) inline, and HTML reports in place.
- **FR-022**: System MUST serve artifact files only from within a registered project root, rejecting traversal and absolute paths outside it with HTTP 403.
- **FR-023**: System MUST surface the run's `stage_contract` and explicitly highlight entries where `enforced != requested`.
- **FR-024**: System MUST show each optional feature's classification (`completed` / `skipped_disabled` / `skipped_missing_dependency` / `failed_error`) from the run manifest.

**DAG visualization (P2):**
- **FR-025**: System MUST render the artifact DAG with one visually distinct treatment per status in `NODE_STATUSES`.
- **FR-026**: System MUST link each DAG node's declared outputs into the results browser.
- **FR-027**: System MUST indicate cache hits distinctly from recomputed nodes.

**Multi-project and concurrency (P2):**
- **FR-028**: System MUST support multiple registered projects with independent dashboards, job lists, and results views.
- **FR-029**: System MUST support concurrent jobs across different projects without serialization.
- **FR-030**: System MUST mark a registered project whose directory is missing as unavailable without erroring the page.

**HPC (P3):**
- **FR-031**: System MUST list HPC profiles from `examples/hpc_profiles/` and `.workflow/hpc_profiles/` and drive `dock deploy` / `dock sync` / `dock submit` through the same CLI path.
- **FR-032**: System MUST NOT render or log secret profile values (accounts, home-directory paths, credentials).
- **FR-033**: System MUST require an explicit confirmation naming the remote target before any sync or submit action.

**Safety and scope (all phases):**
- **FR-034**: System MUST bind to `127.0.0.1` by default and MUST NOT enable a network-exposed bind without an explicit, documented opt-in flag.
- **FR-035**: System MUST NOT modify project files except through the CLI subprocesses it launches; the web layer itself is read-plus-launch only.

### Key Entities

- **Project Registration**: a canonicalized absolute path to an OmniDock project root plus a stable generated id and display name; persisted in the app's own registry file, independent of any project's `.workflow/state.json`.
- **Launch Target**: a runnable unit exposed to the UI — either an analysis target token (`analyze.comparative`, `analyze.interactions.clean`, …) or a pipeline command (`pdb.collect`, `prep.pairlist`, `prep.project`, `dock.run`, `dock.deploy`) — carrying a label, a form schema, and a CLI command template.
- **Job Record**: one launched subprocess: id, project id, target, resolved argv, cwd, log file path, PID, process start time, status (`queued` | `running` | `completed` | `failed` | `cancelled`), exit code, created/started/finished timestamps. Persisted at creation and updated on transition; this is the durable replacement for in-memory `BackgroundTaskManager` futures for web-launched work.
- **Artifact Entry**: one output file surfaced to the user: project-relative path, category (from output topology), media kind (`table` | `image` | `html` | `structure` | `other`), size, and modified time; sourced from `outputs_index.json` where available.
- **Run Summary**: the browser-facing view of `run_tracking/run_manifest.json` — engine identity, `stage_contract` requested-vs-enforced, per-step status, and optional-feature classifications.
- **DAG View Model**: nodes with name, status, inputs, outputs, and cache indicator, derived from the artifact graph report; edges from declared input/output relationships.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: A researcher can go from launching the app to a running analysis on a registered project in under 2 minutes, without typing a CLI command or editing a config file.
- **SC-002**: 100% of launchable targets shown in the UI are drawn from `ANALYSIS_LABELS` / `ENGINE_CHOICES` at runtime, so adding a target in the terminal UI surfaces it in the web UI with no web-side code change (verified by adding a fixture target and asserting it appears).
- **SC-003**: 100% of jobs survive a web-server restart with correct reconciled status (verified by a restart test that launches a sleep-like job, restarts the server, and asserts status continuity).
- **SC-004**: 0 pipeline stages are reimplemented in the web layer — every scientific computation runs via a CLI subprocess (verified by asserting the web package imports nothing from `post_docking_analysis` except read-only manifest/DAG parsing helpers).
- **SC-005**: 100% of path-traversal attempts against the artifact endpoint return 403 (verified by a test matrix of `../`, absolute paths, symlink escapes, and encoded separators).
- **SC-006**: Live log latency from process write to browser display is under 2 seconds for a run producing steady output.
- **SC-007**: Two projects each running a job show independent, correct progress with no cross-contamination of state (verified by a concurrency test).
- **SC-008**: Every one of the eight `NODE_STATUSES` values is visually distinguishable in the DAG view (verified by a fixture DAG report exercising all eight).
- **SC-009**: No secret values from any HPC profile appear in rendered HTML or job logs (verified by grepping a captured session against a known-secret fixture profile).

## Critical Files To Read (not modify)

The web layer reads these; it must not fork or duplicate their logic:

- `workflow/state.py` — state accessors (`ensure_state:138`, `load_state:150`, `list_steps:276`, `summarize_state:384`, `list_background_tasks:319`), `.workflow/state.json` schema (`_default_state:101`)
- `workflow/interactive.py` — `ENGINE_CHOICES:60`, `ANALYSIS_LABELS:68`, `TIMELINE_STEPS:86` (the single source of truth for labels and ordering)
- `workflow/cli.py` — command-group and flag definitions, the contract the web layer shells out to
- `post_docking_analysis/artifact_graph.py` — `NODE_STATUSES:21`, DAG report shape
- `post_docking_analysis/simplified_pipeline_impl.py` — `run_tracking/` writers (`run_manifest.json:242`, `step_status.json:448`, `outputs_index:473`)
- `README.md` — canonical output topology (lines ~652-667)

## New Files To Create

```
webui/
├── __init__.py
├── app.py                  # Flask app factory + route registration
├── config.py               # host/port/registry location settings
├── registry.py             # project registration, canonicalization, persistence
├── state_adapter.py        # read-only bridge onto workflow/state.py
├── targets.py              # launch-target catalog + form schemas (from interactive.py constants)
├── jobs.py                 # durable job records, subprocess launch, reconciliation, cancel
├── logs.py                 # offset-based incremental log reads
├── artifacts.py            # outputs_index parsing, topology grouping, safe path resolution
├── dag.py                  # artifact-graph report → view model
├── hpc.py                  # profile listing + deploy/sync/submit launch (P3)
├── templates/              # base, dashboard, launch, job detail, results, dag
└── static/                 # css + minimal js (log streaming, table sorting)
webui_cli.py                # `python webui_cli.py` / `python main.py webui` entry point
test/test_webui.py          # web-layer test suite
docs/WEB_UI_GUIDE.md        # setup and usage
```

## Existing Utilities To Reuse

- `workflow.state.ensure_state / load_state / list_steps / summarize_state` — FR-002/FR-003 read through these; do not re-parse the JSON.
- `workflow.interactive.ENGINE_CHOICES / ANALYSIS_LABELS / TIMELINE_STEPS` — import these constants directly (they are module-level and side-effect free) so FR-007 and SC-002 hold by construction.
- `workflow.cli` argument definitions — the source of truth for building CLI command templates in `targets.py`; keep template flags aligned with the parsers at `workflow/cli.py:414-660`.
- `post_docking_analysis.artifact_graph.NODE_STATUSES` — import for FR-025 rather than hardcoding status strings.
- `run_tracking/outputs_index.json` and `run_manifest.json` — already structured for machine reading; FR-020/FR-023/FR-024 parse rather than recompute.
- `4-Working/dag_cache.json` and the DAG report written by `ArtifactGraph._write_report:198` — the data source for FR-025/FR-027.
- Spec-kit templates and the existing smoke harness (`test/test_dockforge_smoke.py`) — register new web tests against the same runner.

## Verification

1. **State fidelity**: register a fixture project; diff the UI-rendered context against `.workflow/state.json`; run a CLI command in a separate shell; refresh; assert the UI reflects the change (FR-005).
2. **Target parity**: assert the set of targets offered by the UI equals `set(ANALYSIS_LABELS)` plus the declared pipeline commands; add a fixture entry to `ANALYSIS_LABELS` and assert it appears without web-side changes (SC-002).
3. **Job durability**: launch a long-running fixture job; `kill` the server; restart; assert status is still `running` and the log resumes; let it finish; assert `completed` (SC-003).
4. **Identity reconciliation**: fabricate a job record whose PID has been recycled by another process; assert reconciliation marks it `failed`, not `running` (FR-014).
5. **Path safety**: request artifacts using `../`, `..\\`, URL-encoded separators, an absolute path outside the root, and a symlink pointing outside; assert 403 in every case (SC-005).
6. **Concurrency**: launch jobs in two projects; assert both progress and that each project's dashboard shows only its own job (SC-007).
7. **DAG rendering**: feed a fixture report containing all eight `NODE_STATUSES`; assert each renders distinctly and cache hits are labelled (SC-008).
8. **No logic duplication**: static check that `webui/` imports nothing from `post_docking_analysis` beyond `artifact_graph` constants and manifest reading (SC-004).
9. **Secret hygiene**: load a fixture HPC profile containing sentinel secrets; exercise the HPC panel; grep rendered HTML and job logs for the sentinels; assert zero hits (SC-009).
10. **Non-interactive safety**: launch a target that would prompt on a TTY; assert it exits with a clear error rather than hanging (FR-019).
