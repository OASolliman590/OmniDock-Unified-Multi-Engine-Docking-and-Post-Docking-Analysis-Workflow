# Implementation Plan: OmniDock Web UI Platform

**Branch**: `030-web-ui-platform` | **Date**: 2026-08-16 | **Spec**: `specs/030-web-ui-platform/spec.md`
**Input**: Feature specification from `/specs/030-web-ui-platform/spec.md`

## Summary

Add a locally-hosted web application (`webui/`) that exposes OmniDock's existing workflow state, execution targets, and output artifacts through a browser. The web layer is a **thin adapter, not a second implementation**: it reads state through `workflow/state.py`, sources its vocabulary from `workflow/interactive.py` constants, executes work by shelling out to the existing CLI, and reads results from `run_tracking/` and the artifact-DAG report.

The one genuinely new subsystem is a **durable job runner** (`webui/jobs.py`). It exists because `BackgroundTaskManager` (`workflow/execution.py:132`) is an in-process `ThreadPoolExecutor` whose work dies with the host process, which is unacceptable for hour-long docking runs behind a restartable web server.

## Technical Context

**Language/Version**: Python 3.11 (matches existing project)
**Web framework**: Flask >= 3.0
**Frontend**: Server-rendered Jinja2 templates + vanilla JS. No build step, no npm, no SPA framework.
**Live updates**: Server-Sent Events (`text/event-stream`) for log tailing; short polling for job/DAG status.
**Storage**: JSON files on disk — project registry and job records under an app data directory; no database.
**Process model**: `subprocess.Popen` per job, detached from the request cycle, stdout+stderr redirected to a per-job log file.
**Testing**: `pytest` against a Flask test client, plus fixture project directories.
**Target Platform**: Local workstation (Windows, macOS, Linux) — bound to `127.0.0.1`.
**Project Type**: Additive package inside the existing single-project repo.

### Framework decision (Flask over FastAPI)

Chosen for consistency and implementation cost, not performance:
- The workload is a handful of local requests plus file tailing; async throughput is irrelevant here.
- Server-rendered templates avoid an entire frontend toolchain, which matters when the implementation is intended for a smaller model working task-by-task.
- Flask handles SSE adequately via a generator response, which covers the only streaming requirement.
- A sibling internal tool (`claritydynamics_ui`) already uses Flask + Jinja2 + vanilla JS, so patterns transfer.

If a future feature needs many concurrent streaming clients or a typed API surface, FastAPI is the migration target; the adapter modules (`registry`, `state_adapter`, `targets`, `jobs`, `artifacts`, `dag`) are deliberately framework-free so only `app.py` would change.

### Dependency additions

Add to `requirements.txt`:
```
flask>=3.0
```
Everything else (`pandas` for CSV rendering, `pyyaml`) is already present. No new scientific dependencies.

## Architecture

```
Browser (Jinja2 pages + fetch/SSE)
        |
webui/app.py  — routes only, no logic
        |
   +----+----+----------+-----------+---------+-------+
   |         |          |           |         |       |
registry  state_    targets       jobs    artifacts  dag
  .py     adapter     .py         .py       .py      .py
   |         |          |           |         |       |
   |   workflow/    workflow/   subprocess  run_    artifact
   |    state.py   interactive   Popen →   tracking/ graph
   |               .py + cli.py  main.py    files    report
   |
app data dir (projects.json, jobs/*.json, logs/*.log)
```

**Hard rule enforced by SC-004**: `webui/` may import `workflow.state`, `workflow.interactive` (constants only), and `post_docking_analysis.artifact_graph.NODE_STATUSES`. It may not import pipeline execution modules. All computation happens in subprocesses.

### Why shelling out rather than importing `workflow.execution`

1. **Isolation** — a crash or memory blowup in a docking stage cannot take down the web server.
2. **Durability** — a subprocess outlives the Flask process, satisfying FR-013.
3. **Auditability** — every action maps to a copy-pasteable CLI command (FR-008), so the UI can never do something the CLI cannot.
4. **Zero drift** — CLI argument parsing stays the single contract; the web layer builds argv and nothing more.

The cost is per-run interpreter startup, which is negligible against multi-minute-to-multi-hour stages.

### Job lifecycle

```
create record (status=queued, persisted)  ──►  Popen  ──►  status=running (+pid, +start_time)
                                                              │
                        ┌─────────────────────────────────────┼──────────────────────┐
                        ▼                                     ▼                      ▼
              exit 0 → completed                  exit != 0 → failed(code)     cancel → cancelled
```

Reconciliation on server startup, for every record not in a terminal state:
1. No PID recorded → `failed` ("never started").
2. PID not alive → read exit marker file if present, else `failed` ("process vanished").
3. PID alive **and** identity matches recorded process start time → keep `running`.
4. PID alive but identity mismatch (recycled PID) → `failed` ("pid recycled"). This is FR-014 and must not be skipped.

Identity check uses `psutil` if available; otherwise fall back to a platform-specific process start-time probe, and if neither is available, treat "alive but unverifiable" as `running` while recording that verification was unavailable. Do not add `psutil` as a hard dependency for this alone — make it optional.

## Data layout

App data directory (default `~/.omnidock_webui/`, overridable by `OMNIDOCK_WEBUI_HOME`):

```
projects.json          # [{id, path, name, registered_at}]
jobs/<job_id>.json     # one durable job record
logs/<job_id>.log      # combined stdout+stderr for that job
```

Kept outside project directories so the UI never pollutes scientific output trees (FR-035), and so a project can be registered read-only.

## Phasing

Phases map to spec priorities. Each phase is independently demoable and testable.

| Phase | Stories | Delivers |
|---|---|---|
| 1 | infrastructure | Flask scaffold, config, app data dir, health endpoint |
| 2 | US1 (P1) | Project registry + dashboard + timeline |
| 3 | US3 (P1) | Durable job runner + reconciliation (built *before* launch UI) |
| 4 | US2 (P1) | Target catalog, launch forms, command preview |
| 5 | US3 (P1) | Live log streaming + cancel |
| 6 | US4 (P1) | Results browser + safe artifact serving + run summary |
| 7 | US5 (P2) | DAG view |
| 8 | US6 (P2) | Multi-project switching + concurrency guards |
| 9 | US7 (P3) | HPC panel |
| 10 | — | Tests, docs, hardening |

**Ordering note**: Phase 3 (job durability) precedes Phase 4 (launch UI) deliberately. Building the launch forms first would invite an in-process execution shortcut that the durability requirement later has to tear out.

## Validation Strategy

1. **Per-phase**: `python -m py_compile` on touched modules; `pytest test/test_webui.py -k <phase>`.
2. **Fixture projects**: create `test/fixtures/webui_project/` containing a minimal `.workflow/state.json`, a `run_tracking/` set, and a DAG report exercising all eight `NODE_STATUSES`.
3. **Restart test**: the durability guarantee (SC-003) is verified by an actual process kill/restart cycle, not a mock.
4. **Security tests**: path-traversal matrix (SC-005) and secret-leak grep (SC-009) run in CI-equivalent form before the feature is considered closed.
5. **Parity check**: SC-002 asserted by monkeypatching an extra entry into `ANALYSIS_LABELS` and asserting it surfaces.
6. **Import discipline**: SC-004 asserted by a static test that walks `webui/` imports and fails on disallowed modules.

## Risks

| Risk | Impact | Mitigation |
|---|---|---|
| PID reuse causes a dead job to read as running | Wrong status, user waits forever | FR-014 identity verification; explicit test |
| A target still expects a TTY prompt | Job hangs forever holding a slot | FR-019 non-interactive env; fail-fast test |
| Two concurrent jobs corrupt one project's outputs | Scientific data loss | FR-018 warn-on-second-job; document as a guard, not a lock |
| Path traversal exposes files outside a project | Local file disclosure | FR-022 resolve-and-contain check with a dedicated test matrix |
| Large log files exhaust memory | Server unresponsive | FR-015 offset reads with a bounded tail window |
| Web UI drifts from terminal UI vocabulary | Two divergent products | FR-007 import constants at runtime; SC-002 test |
| Scope creep into a hosted multi-user service | Security surface the design never intended | Non-Goals section; FR-034 localhost bind |

## Out of Scope for This Plan

Authentication, HTTPS, multi-user access control, remote hosting, database-backed storage, real-time collaborative editing, and any change to scientific pipeline behavior. If a hosted deployment is ever wanted, it needs its own spec covering authn/authz, transport security, and resource isolation — none of which this feature attempts.
