# Data Model: OmniDock Web UI

Two categories of data:

- **Owned** — created and written by the web layer, stored in the app data directory.
- **Borrowed** — produced by OmniDock itself, read-only to the web layer.

The web layer never writes borrowed data (FR-035).

---

## Owned data

Location: `$OMNIDOCK_WEBUI_HOME` (default `~/.omnidock_webui/`), deliberately outside project trees.

```
~/.omnidock_webui/
├── projects.json           # registry
├── jobs/
│   ├── <job_id>.json       # one durable job record
│   └── <job_id>.exit       # exit-code marker written when the process ends
└── logs/
    └── <job_id>.log        # combined stdout+stderr
```

### Project

`projects.json` holds an array of:

| Field | Type | Notes |
|---|---|---|
| `id` | string | stable; first 8 hex chars of sha256 of the canonical path |
| `path` | string | canonical absolute path (T011: `~` expanded, symlinks resolved, Windows case-normalized) |
| `name` | string | display name; defaults to the directory basename |
| `registered_at` | string | ISO 8601 UTC |

Computed at read time, never stored:

| Field | Type | Notes |
|---|---|---|
| `available` | bool | directory currently exists and is readable (FR-030) |
| `has_workflow` | bool | `.workflow/state.json` present |

**Invariant**: two registrations whose canonical paths are equal must produce one entry. `id` derives from the canonical path, so this holds by construction.

### Job

One file per job at `jobs/<job_id>.json`.

| Field | Type | Notes |
|---|---|---|
| `job_id` | string | `job_<UTC timestamp>_<8 hex>` — sortable and unique |
| `project_id` | string | FK to a Project `id` |
| `target` | string | launch-target key, e.g. `analyze.comparative` |
| `argv` | string[] | the exact argument vector passed to `Popen` |
| `command` | string | display-only rendering of `argv` (FR-008) |
| `cwd` | string | always the project path |
| `log_path` | string | absolute path to the job log |
| `pid` | int \| null | null until launched |
| `proc_start_time` | float \| null | process creation time; the identity token for FR-014 |
| `status` | enum | `queued` \| `running` \| `completed` \| `failed` \| `cancelled` |
| `exit_code` | int \| null | set on terminal transition |
| `created_at` | string | ISO 8601 UTC, set before launch (FR-012) |
| `started_at` | string \| null | when `Popen` returned |
| `finished_at` | string \| null | terminal transition time |
| `note` | string | reason text for `failed` / `cancelled` (e.g. `"pid recycled"`) |

**State machine** — the only legal transitions:

```
queued ──► running ──► completed        (exit 0)
   │          ├──────► failed           (exit != 0, vanished, pid recycled)
   │          └──────► cancelled        (user request)
   └─────────────────► failed           (launch raised)
```

`TERMINAL_STATUSES = {completed, failed, cancelled}`. A terminal job is never re-reconciled.

**Reconciliation rules** (run once at server startup, T035):

| Recorded state | Observation | Result |
|---|---|---|
| non-terminal | `pid` is null | `failed`, note `"never started"` |
| non-terminal | process not alive, exit marker present | `completed` / `failed` from the marker code |
| non-terminal | process not alive, no marker | `failed`, note `"process vanished"` |
| non-terminal | alive, identity matches | stays `running` |
| non-terminal | alive, identity mismatch | `failed`, note `"pid recycled"` |
| non-terminal | alive, identity unverifiable | stays `running`, note records that verification was unavailable |

**Invariants**:
1. A job record exists on disk before its process starts — a crash between the two leaves a `queued` record, never an untracked process.
2. `argv[0:2] == ["python", "main.py"]` — every job is a CLI invocation (FR-010, SC-004).
3. `cwd` always equals the owning project's path.

---

## Borrowed data (read-only)

### Workflow state — `<project>/.workflow/state.json`

Schema `version: 2`, written by `workflow/state.py`. Read only through its accessors (FR-002). Fields consumed by the UI:

| Path | Used for |
|---|---|
| `current_context.favorite_engine` | form defaults, dashboard |
| `current_context.selected_engines` | form defaults |
| `current_context.analysis_output_dir` | results discovery |
| `current_context.docking_root`, `post_docking_root` | dashboard, artifact roots |
| `current_context.pair_mode`, `layout_profile` | dashboard |
| `steps` | timeline status (via `list_steps()`) |
| `feature_flags` | dashboard |
| `background_tasks` | display only — **not** the web job source of truth |
| `checkpoint_metadata` | dashboard (future) |

**Important**: `background_tasks` records interactive-session tasks from `BackgroundTaskManager`. Web jobs live in the owned store. Reconciliation must never treat a stale `background_tasks` entry as a web job (Edge Case: orphaned running entries).

### Run tracking — `<run_dir>/run_tracking/`

| File | Used for |
|---|---|
| `run_manifest.json` | run summary: engine, `stage_contract`, optional-feature classification (FR-023/FR-024) |
| `step_status.json` | per-step status list |
| `outputs_index.json` | authoritative artifact list (FR-020) |
| `outputs_index.csv` | ignored — JSON preferred |

Optional-feature classifications (spec 029 FR-008): `completed`, `skipped_disabled`, `skipped_missing_dependency`, `failed_error`.

### Artifact DAG

| Source | Used for |
|---|---|
| DAG report from `ArtifactGraph._write_report` | node names, statuses, inputs/outputs |
| `4-Working/dag_cache.json` | cache-hit indicators (FR-027) |

Node statuses are exactly `post_docking_analysis.artifact_graph.NODE_STATUSES`: `pending`, `running`, `completed`, `completed_with_warnings`, `failed`, `skipped_optional`, `blocked_by_failure`, `cache_hit`. Import the constant; never restate the list (T081).

### Output topology

Artifact categories map to the canonical layout:

| Directory | Category |
|---|---|
| `analysis/` | `analysis` |
| `complexes/` | `complexes` |
| `best_poses/` | `best_poses` |
| `reports/` | `reports` |
| `rmsd_analysis/` | `rmsd_analysis` |
| `interactions/{prolif,ligplot,pandamap,poseview}/` | `interactions.<tool>` |
| `3d_visualizations/` | `3d_visualizations` |
| `visualizations/` | `visualizations` |
| `raw_data/` | `raw_data` |
| anything else | `other` |

Media kinds: `.csv`/`.tsv` → `table`; `.png`/`.svg`/`.jpg` → `image`; `.html` → `html`; `.pdb`/`.sdf`/`.pdbqt` → `structure`; otherwise `other`.

---

## Vocabulary imported at runtime

Never duplicated in the web layer (FR-007, SC-002):

| Constant | Module | Used for |
|---|---|---|
| `ENGINE_CHOICES` | `workflow/interactive.py:60` | engine select options |
| `ANALYSIS_LABELS` | `workflow/interactive.py:68` | analysis target catalog and labels |
| `TIMELINE_STEPS` | `workflow/interactive.py:86` | dashboard timeline order and labels |
| `NODE_STATUSES` | `post_docking_analysis/artifact_graph.py:21` | DAG status vocabulary |

---

## Path safety model

Every artifact request passes through `safe_resolve(project_root, rel_path)` (T066):

1. Reject absolute inputs and any component equal to `..` before joining.
2. Join to `project_root`, then fully resolve (following symlinks).
3. Assert the resolved path is a descendant of the resolved `project_root`.
4. Raise `PathEscape` → HTTP 403 otherwise.

Resolution happens **after** joining, so a symlink inside the project pointing outside it is still caught (FR-022, SC-005).
