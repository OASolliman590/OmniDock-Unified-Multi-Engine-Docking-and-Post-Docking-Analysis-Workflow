# API Contract: OmniDock Web UI

All endpoints are served from `127.0.0.1` only (FR-034). All responses are JSON unless stated otherwise. All errors use the shape:

```json
{ "ok": false, "error": "human-readable message", "code": "MACHINE_CODE" }
```

Error codes: `NOT_FOUND`, `INVALID_INPUT`, `PATH_ESCAPE`, `PROJECT_UNAVAILABLE`, `JOB_RUNNING`, `STATE_UNREADABLE`, `CONFIRM_REQUIRED`.

---

## Health

### `GET /api/health`
```json
{ "ok": true, "version": "0.1.0" }
```

---

## Projects

### `GET /api/projects`
```json
[
  { "id": "a1b2c3d4", "path": "/abs/path/proj", "name": "MRSA X19",
    "registered_at": "2026-08-16T10:00:00Z", "available": true, "has_workflow": true }
]
```

### `POST /api/projects`
Request: `{ "path": "/abs/path/proj", "name": "optional" }`
Response `201`: the project object. Re-registering an equivalent path returns `200` with the existing entry (T011/T012).
Errors: `INVALID_INPUT` (missing/nonexistent/unreadable path).

### `GET /api/projects/<pid>/state`
```json
{
  "has_workflow": true,
  "current_context": { "favorite_engine": "gnina", "selected_engines": ["gnina","vina"],
                       "docking_root": "...", "analysis_output_dir": "...", "pair_mode": "..." },
  "feature_flags": { },
  "artifacts": { },
  "summary": ["line 1", "line 2"]
}
```
`summary` comes from `workflow.state.summarize_state`. Errors: `STATE_UNREADABLE` after one retry.

### `GET /api/projects/<pid>/timeline`
Order follows `TIMELINE_STEPS` exactly (FR-003).
```json
[
  { "key": "pdb.collect", "label": "Fetch proteins from PDB",
    "status": "completed", "timestamp": "2026-08-01T09:00:00Z", "note": "" },
  { "key": "dock.run", "label": "Run docking", "status": "pending", "timestamp": null, "note": "" }
]
```
`status` is whatever `list_steps()` recorded, or `"pending"` when absent.

### `POST /api/projects/<pid>/init`
Launches `workflow init` as a job. Response `202`: the job object.

---

## Targets

### `GET /api/targets`
```json
[
  { "key": "analyze.comparative", "label": "Comparative analysis", "group": "analysis" },
  { "key": "analyze.interactions.prolif", "label": "ProLIF target (routes to clean interactions pipeline)",
    "group": "analysis", "routes_to": "analyze.interactions.clean" },
  { "key": "dock.run", "label": "Run docking", "group": "docking" }
]
```
Built by iterating `ANALYSIS_LABELS` plus declared pipeline targets — never a hardcoded list (SC-002).

### `GET /api/targets/<key>/form?project_id=<pid>`
```json
{
  "key": "analyze.comparative",
  "label": "Comparative analysis",
  "routes_to": null,
  "fields": [
    { "name": "engines", "label": "Engines", "type": "multiselect", "required": true,
      "default": ["gnina","vina"],
      "choices": [ { "value": "gnina", "label": "GNINA" }, { "value": "vina", "label": "Vina" } ] },
    { "name": "output_dir", "label": "Output directory", "type": "path",
      "required": false, "default": "/abs/path/proj/3-Results" }
  ]
}
```
`choices` for engine fields come from `ENGINE_CHOICES`. Defaults are pre-filled from the project's `current_context` (T049).

---

## Jobs

### `POST /api/projects/<pid>/jobs`
Request:
```json
{ "target": "analyze.comparative",
  "values": { "engines": ["gnina","vina"], "output_dir": "/abs/path/proj/3-Results" },
  "confirm_concurrent": false }
```
Response `202`:
```json
{ "job_id": "job_20260816T101500_ab12cd34", "project_id": "a1b2c3d4",
  "target": "analyze.comparative",
  "argv": ["python","main.py","analyze","comparative","--project-dir","/abs/path/proj","--engines","gnina,vina"],
  "command": "python main.py analyze comparative --project-dir /abs/path/proj --engines gnina,vina",
  "status": "running", "pid": 48213, "log_path": "/home/u/.omnidock_webui/logs/job_….log",
  "created_at": "2026-08-16T10:15:00Z", "started_at": "2026-08-16T10:15:00Z" }
```
Errors:
- `400 INVALID_INPUT` — validation failed; body includes `"fields": {"engines": "unknown engine 'foo'"}`.
- `409 JOB_RUNNING` — a job is already running for this project; retry with `confirm_concurrent: true` (FR-018).

### `GET /api/jobs?project_id=<pid>`
Array of job objects, newest first.

### `GET /api/jobs/<job_id>`
One job object, including `exit_code`, `finished_at`, and `note` when terminal.

### `POST /api/jobs/<job_id>/cancel`
Response `200`: the updated job with `status: "cancelled"`. Terminating an already-terminal job returns `200` unchanged.

### `GET /api/jobs/<job_id>/log?offset=<int>`
```json
{ "text": "…new bytes…", "offset": 40960, "size": 40960, "status": "running" }
```
Incremental only — never returns the whole file when `offset` is supplied (FR-015).

### `GET /api/jobs/<job_id>/stream`
`Content-Type: text/event-stream`. Events:
```
event: log
data: {"text":"…","offset":41200}

event: status
data: {"status":"completed","exit_code":0}
```
Server closes the stream once a terminal status is emitted.

---

## Runs, Artifacts, Files

### `GET /api/projects/<pid>/runs`
```json
[ { "run_dir": "3-Results/comparative_20260816", "has_manifest": true,
    "modified_at": "2026-08-16T12:00:00Z" } ]
```

### `GET /api/projects/<pid>/runs/<run>/summary`
```json
{
  "engine": "gnina",
  "stage_contract": [
    { "stage": "run_rmsd", "requested": true, "enforced": true, "matches": true },
    { "stage": "run_visualizations", "requested": false, "enforced": true, "matches": false }
  ],
  "steps": [ { "name": "affinity", "status": "completed" } ],
  "optional_features": [
    { "name": "poseview", "classification": "skipped_missing_dependency", "detail": "binary absent" }
  ]
}
```
Rows with `matches: false` must be visually flagged (FR-023).

### `GET /api/projects/<pid>/artifacts?run=<run>`
```json
[ { "path": "3-Results/…/reports/consensus_ranked.csv", "category": "reports",
    "media_kind": "table", "size": 20481, "modified_at": "2026-08-16T12:00:00Z" } ]
```
Sourced from `run_tracking/outputs_index.json`, falling back to a directory walk (FR-020).

### `GET /api/projects/<pid>/file?path=<project-relative>`
Returns raw file bytes with a guessed content type.
Errors: `403 PATH_ESCAPE` for any path resolving outside the project root — including `../`, encoded separators, absolute paths, and symlink escapes (FR-022, SC-005). `404 NOT_FOUND` when missing.

### `GET /api/projects/<pid>/table?path=<p>&page=<n>&page_size=<n>`
```json
{ "columns": ["ligand","receptor","score"],
  "rows": [["lig1","rec1",-9.2]],
  "page": 1, "page_size": 100, "total_rows": 4821 }
```

---

## DAG

### `GET /api/projects/<pid>/dag?run=<run>`
```json
{
  "nodes": [ { "name": "consensus_ranked", "status": "completed", "cached": false,
               "inputs": ["normalized_scores.csv"], "outputs": ["consensus_ranked.csv"] } ],
  "edges": [ { "from": "normalized_scores", "to": "consensus_ranked" } ],
  "status_counts": { "completed": 12, "cache_hit": 3, "skipped_optional": 2 }
}
```
`status` values are exactly the members of `NODE_STATUSES` (FR-025).

---

## HPC (P3)

### `GET /api/projects/<pid>/hpc/profiles`
```json
[ { "name": "bibalex-apptainer-conda", "source": "builtin",
    "fields": { "partition": "gpu", "account": "[redacted]" } } ]
```
Secret-shaped values are always redacted before returning (FR-032, SC-009).

### `POST /api/projects/<pid>/hpc/<action>`
`action` ∈ `deploy` | `sync` | `submit`. Request:
```json
{ "profile": "bibalex-apptainer-conda", "engines": ["gnina"], "confirm_target": "cluster.example.edu" }
```
`sync` and `submit` require `confirm_target` matching the profile's remote host, else `400 CONFIRM_REQUIRED` (FR-033).
Response `202`: a job object.

---

## Page Routes (HTML)

| Route | Template | Purpose |
|---|---|---|
| `GET /` | `projects.html` | Register and list projects |
| `GET /project/<pid>` | `dashboard.html` | Timeline, context, recent jobs |
| `GET /project/<pid>/launch` | `launch.html` | Target picker + form + command preview |
| `GET /project/<pid>/results` | `results.html` | Run selector, artifacts, run summary |
| `GET /project/<pid>/dag` | `dag.html` | Artifact DAG |
| `GET /project/<pid>/hpc` | `hpc.html` | HPC panel (P3) |
| `GET /job/<job_id>` | `job.html` | Status, command, live log, cancel |
