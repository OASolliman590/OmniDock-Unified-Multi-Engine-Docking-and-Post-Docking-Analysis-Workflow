# OmniDock Web UI

A browser interface over the existing OmniDock workflow: register a project,
inspect its state, launch preparation/docking/analysis runs from forms instead
of terminal prompts, watch them stream live, and browse the results.

The web layer is a thin adapter. It reads project state through
`workflow/state.py`, takes its vocabulary from `workflow/interactive.py`, and
executes everything by launching the existing CLI as a subprocess. No pipeline
stage is reimplemented, so anything you can do here you could also do from the
command line — and the UI shows you that command before it runs.

## Security posture — read this first

- **Localhost only.** The server binds `127.0.0.1` by default. There is **no
  authentication of any kind**. Anyone who can reach the port can launch jobs
  and read files inside registered projects.
- Binding elsewhere requires the explicit `--allow-remote` flag, and you should
  not use it on an untrusted network.
- The web process never writes into a project directory. Only the CLI
  subprocesses it launches do.
- HPC profile values (accounts, SSH targets, home paths) are redacted before
  reaching the browser and are never written to a job log.

## Install

From the repository root:

```bash
pip install -r requirements.txt
```

`flask` is the only addition for the web UI. `psutil` is optional but
recommended — without it the server cannot verify that a recorded process ID
still belongs to the job it launched, and will say so rather than guess.

## Run

```bash
python main.py webui
```

or equivalently:

```bash
python webui_cli.py
```

Then open <http://127.0.0.1:8770>.

Options: `--host`, `--port`, `--debug`, and `--allow-remote` (see above).

### Windows note

`main.py` prints a non-ASCII banner. On a console defaulting to cp1252 that
raises `UnicodeEncodeError` before any work happens, so set:

```bash
set PYTHONIOENCODING=utf-8
```

before running the CLI directly. Jobs launched *by the web UI* already set this
in their environment, so they are unaffected.

## Register a project

On the landing page, enter the absolute path to an OmniDock project root.

Paths are canonicalized, so the same directory cannot be registered twice under
different spellings (trailing separator, different case on Windows, or via a
symlink).

If the directory has no `.workflow/` yet, the dashboard offers to initialize it,
which runs `main.py workflow init` as a job.

> On Windows, `workflow init` creates a backward-compatibility symlink and
> therefore needs Developer Mode or an elevated shell. Without one it fails with
> `WinError 1314`. Existing projects are unaffected.

## Dashboard

Shows the workflow timeline (the same steps and ordering as the terminal UI),
the project's current context — favorite engine, selected engines, docking and
analysis roots — and this project's recent jobs.

State is re-read on every request, so a run started from a terminal in another
window shows up here on refresh.

## Launch

Pick a target, fill in the form, and check the command preview before starting.
The preview is rendered by the same code that builds the real invocation, so it
is exactly what will run.

Targets come from `ANALYSIS_LABELS` and `ENGINE_CHOICES` at runtime — add one to
the terminal UI and it appears here with no change to the web code.

Notes:

- Fields are pre-filled from the project's context where possible.
- Some interaction targets (`prolif`, `ligplot`, `pandamap`) are routed by the
  pipeline through `analyze.interactions.clean`. The form says so rather than
  silently substituting.
- Launching a second job while one is already running for that project asks for
  confirmation first: pipeline stages are not guaranteed safe to run
  concurrently within one project.

## Jobs

Every launch becomes a durable job record, written to disk **before** the
process starts.

- Live log output streams to the browser; the view tails by byte offset, so a
  multi-hour log will not exhaust memory.
- Jobs survive a server restart. They run detached under a small supervisor
  that records the exit code, so a restarted server can recover the outcome of
  a process it never held a handle to.
- Reconciliation verifies process *identity*, not just that some process holds
  the recorded PID — the OS reuses PIDs, and a dead job would otherwise read as
  running forever.
- Cancelling terminates the process group and leaves any files already produced
  in place.

Statuses: `queued`, `running`, `completed`, `failed` (with exit code),
`cancelled`.

## Results

Select a run (any directory containing `run_tracking/`) to see:

- **Run summary** — engine, the stage contract with any `enforced != requested`
  override called out, per-step status, and each optional feature classified as
  completed, skipped (disabled), skipped (missing dependency), or failed.
- **Artifacts** grouped by the canonical output topology, with CSV/TSV rendered
  as sortable paginated tables, images inline, and HTML reports in place.

Artifacts come from `run_tracking/outputs_index.json` when present, falling back
to a directory walk so a run that failed partway still shows its files.

Every file request is confined to the project root; traversal attempts are
refused.

## DAG

If a run used the artifact graph, its execution report renders as dependency
tiers: everything in a tier can run in parallel, and each tier waits for the
ones before it. Each node shows its status, whether it was a cache hit, whether
it is optional, and what blocked it.

The report records status and blocking relationships but not each node's
declared inputs and outputs, so per-artifact edges are not drawn.

## HPC

Lists profiles from `examples/hpc_profiles/` (public templates) and the
project's `.workflow/hpc_profiles/` (real settings), always redacted.

- `deploy` generates Slurm assets locally and needs no confirmation.
- `sync` and `submit` contact a remote system and require you to type the
  profile's SSH target to confirm.

## Where the UI keeps its own data

Outside your projects, under `~/.omnidock_webui/` (override with
`OMNIDOCK_WEBUI_HOME`):

```
projects.json        registry of registered project paths
jobs/<id>.json       one durable job record
jobs/<id>.exit       exit-code marker written by the job's supervisor
logs/<id>.log        combined stdout and stderr for that job
```

Deleting this directory forgets registrations and job history; it never touches
project data.

## Tests

```bash
pytest test/test_webui.py -v
```

## Specification

The governing spec, plan, task list, API contract and data model live in
[`specs/030-web-ui-platform/`](../specs/030-web-ui-platform/).
