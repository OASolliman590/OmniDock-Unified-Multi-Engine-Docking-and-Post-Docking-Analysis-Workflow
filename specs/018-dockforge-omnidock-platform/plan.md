# Implementation Plan: Omni-DockForge: End-to-End Docking, Consensus by Design.

**Branch**: `018-dockforge-omnidock-platform` | **Date**: 2026-03-30 | **Spec**: `/specs/018-dockforge-omnidock-platform/spec.md`  
**Input**: Feature specification from `/specs/018-dockforge-omnidock-platform/spec.md`

## Summary

This plan upgrades the existing docking workflow into DockForge, a modular omni-docking platform with:
- non-blocking interactive workflow orchestration,
- engine-aware ligand preparation branching,
- environment-aware docking execution,
- post-docking analysis as first-class project entry,
- cross-engine allscore consensus with explainability,
- reproducibility via checkpointing/provenance manifests,
- future-ready extension points (ensemble, dashboard, rerun automation).

The implementation is staged to keep current workflows operational while introducing new contracts incrementally.

## Technical Context

**Language/Version**: Python 3.9+  
**Primary Dependencies**: pandas, numpy, matplotlib, seaborn, RDKit, Open Babel (optional), Meeko (optional), AutoDockTools prepare scripts (optional), existing docking engine binaries  
**Storage**: File-first (CSV + JSON + project folders), optional SQLite backend in later milestone  
**Testing**: Existing script-based smoke tests (`test_pipeline.py`) + new targeted integration/smoke scripts  
**Target Platform**: macOS/Linux local environments, remote/HPC adapter-compatible execution  
**Project Type**: Single Python project with CLI + interactive orchestration + packaged post-docking module  
**Performance Goals**:
- interactive progress refresh <= 2s under normal local loads,
- no UI blocking for long-running steps,
- preflight validation before expensive docking tasks
**Constraints**:
- preserve backwards compatibility for existing projects where feasible,
- graceful degradation for missing optional tools,
- deterministic output structure for auditability
**Scale/Scope**:
- multi-protein multi-ligand projects,
- multi-engine comparative workflows,
- long-running docking/analysis sessions with resume support

## Constitution Check

*GATE: Must pass before implementation starts. Re-check after design updates.*

1. **Pipeline-First Structure**: PASS
- Plan keeps entry scripts thin and pushes new logic into reusable modules (`workflow/`, `docking/`, `post_docking_analysis/`).

2. **Reproducible Tooling and Configuration**: PASS
- Introduces `.meta` run manifests, config freezes, environment locks, and documented defaults.

3. **Scripted Verification**: PASS
- Adds milestone-specific smoke/integration commands and parity checks for key flows.

4. **Non-Destructive Outputs and Naming**: PASS
- Preserves deterministic, non-destructive project topology and sessionized outputs.

5. **Graceful Optional-Tool Behavior**: PASS
- All optional tools stay capability-gated with explicit fallback/error messaging.

## Project Structure

### Documentation (this feature)

```text
specs/018-dockforge-omnidock-platform/
├── spec.md
├── plan.md
└── tasks.md
```

### Source Code (repository root)

```text
main.py
interactive_pipeline.py
core_pipeline.py
cli_pipeline.py

workflow/
├── cli.py
├── interactive.py
├── execution.py
├── models.py
└── state.py

docking/
├── cli.py
├── models.py
├── engine_registry.py
├── deployment.py
├── remote_ops.py
├── project_layout.py
├── preparation/
│   ├── ligand_preparation.py
│   ├── ligand_quality.py
│   └── pairlist_builder.py
└── runners/
    ├── base.py
    ├── vina.py
    ├── smina.py
    ├── gnina.py
    └── autodock4.py

autodock_preparation.py
prep_autodock_enhanced.sh

post_docking_analysis/
├── cli.py
├── pipeline.py
├── multi_engine_pipeline.py
├── consensus.py
├── structure_quality.py
├── report_generator.py
├── visualizer.py
├── py3dmol_visualizer.py
├── pymol_visualizer.py
└── plugins/
```

**Structure Decision**: Keep the current single-project Python structure and add new capabilities by expanding existing workflow/docking/post_docking modules instead of introducing a new service split.

## Architecture and Component Boundaries

1. **Interactive Orchestration Layer**
- Owns menus, panel navigation, timeline rendering, and user prompts.
- Delegates execution to background task manager.

2. **Workflow State Layer**
- Owns durable stage-step statuses, validation outcomes, and session persistence.
- Provides read APIs for timeline and progress views.

3. **Task Execution Layer**
- Owns asynchronous job lifecycle, progress callbacks, and cancellation hooks.
- Abstracts local/env-specific execution runtime.

4. **Preparation Orchestration Layer**
- Owns ligand prep graph selection and engine compatibility checks.
- Coordinates Open Babel, Meeko, and AutoDockTools adapters.

5. **Docking Engine Layer**
- Keeps runner interfaces uniform while preserving per-engine parameter extensions.
- Supports environment/profile-specific command builders.

6. **Post-Docking Analysis Layer**
- Owns score ingestion, normalization, consensus, QC modules, correlation, and reporting.
- Runs independent from docking stage when requested.

7. **Reproducibility and Data Layer**
- Owns `.meta` manifests, run hashes, checkpoint lineage, optional SQLite mirror.

## Milestone Implementation Phases

### Milestone 1: Core Workflow and UX Scaffolding

Deliverables:
- timeline/progress model integrated in interactive workflow,
- universal `Back` action in all panels,
- non-blocking background task manager,
- `Checkpoint & Revise` terminology and state model,
- execution abstraction skeleton.

Dependencies:
- none beyond existing workflow framework.

### Milestone 2: Ligand Preparation Refactor

Deliverables:
- Open Babel enrichment pipeline (H/protonation/3D/minimize/charges),
- engine-aware preparation graph and compatibility validator,
- redesigned preparation menu/modes,
- PLIP removed from preparation stage,
- output validation reports per ligand.

Dependencies:
- Milestone 1 state/task scaffolding.

### Milestone 3: Docking Orchestration and Runtime Profiles

Deliverables:
- environment-aware run docking flow,
- preflight checks for required files and engine assets,
- common + engine-specific parameter schemas,
- GNINA GPU/CNN parameter controls,
- local execution path plus backend adapter hooks.

Dependencies:
- Milestone 2 preparation contracts.

### Milestone 4: Post-Docking Promotion

Deliverables:
- project-centric analysis entry and modular analysis scopes,
- raw score ingestion contracts,
- normalization + allscore + consensus outputs,
- cross-engine correlation outputs,
- validated complex structure exports.

Dependencies:
- Milestone 3 docking output contracts.

### Milestone 5: QC, Biology Integration, and Rescoring

Deliverables:
- upstream ADMET and ligand QC filters,
- receptor QC gate before docking,
- decoy/reference benchmarking (EF/ROC where available),
- external biology file mapping + correlation reports,
- OnionNet2 plugin integration point,
- configurable target-aware strong/moderate/weak framework.

Dependencies:
- Milestone 4 allscore foundation.

### Milestone 6: Reproducibility and Reporting

Deliverables:
- `.meta` manifest system and environment snapshots,
- durable session/checkpoint lineage,
- SQLite backend evaluation/adoption with CSV parity,
- consolidated reporting outputs,
- extension points for ensemble/dash/rerun automation.

Dependencies:
- Milestones 1-5 stabilized output contracts.

## Dependency Ordering and Critical Path

Critical path:
1. workflow state + background execution primitives,
2. preparation compatibility and output contracts,
3. runtime/environment abstraction,
4. allscore/consensus core,
5. QC/classification integration,
6. reproducibility/storage hardening.

Parallelizable lanes:
- UX timeline rendering and session persistence (after state model exists),
- engine parameter schema design and execution adapter skeleton,
- post-docking module refactor and correlation module development,
- reporting templates and `.meta` manifest generation.

## Migration Strategy (Current Pipeline -> DockForge)

1. **Alias then rename**
- Introduce DockForge naming while maintaining temporary compatibility labels where needed.

2. **Additive migration of folder topology**
- New projects use canonical structure immediately.
- Existing projects gain migration utility that maps old outputs into canonical folders non-destructively.

3. **Dual-mode compatibility window**
- Keep legacy prep/run flows available behind compatibility toggle while new engine-aware graph stabilizes.

4. **Session state backfill**
- Existing sessions can be imported with inferred statuses (`needs_review` where uncertainty exists).

5. **Consensus/data parity checks**
- Validate new allscore outputs against prior comparative outputs before default switch.

## Data Storage Plan

Phase strategy:
- Phase A (Milestones 1-4): file-first storage (`CSV/JSON`) with stricter schema contracts and canonical paths.
- Phase B (Milestone 6): optional SQLite backend (`results.db`) as query accelerator and cross-run index.

SQLite adoption criteria:
- zero regression in ranking outputs vs CSV parity tests,
- acceptable performance gain for multi-protein joins/correlations,
- clear export path to existing CSV reports.

## Job Execution Model

- Introduce `TaskManager` abstraction with task lifecycle:
  - queued -> in_progress -> completed/failed/cancelled.
- Heavy jobs (docking/post-docking modules) run in background workers.
- Interactive layer polls/receives progress snapshots.
- Task metadata persists in `sessions/<interactive_session>/task_state.json`.
- Execution backends use adapter interface:
  - `LocalCPUAdapter`, `LocalGPUAdapter`, `CondaAdapter`, `ContainerAdapter`, `RemoteAdapter` (scaffold first, expand iteratively).

## Testing and Validation Strategy

1. **Workflow UX smoke tests**
- Non-blocking navigation during running tasks.
- Universal back navigation across all panels.
- Timeline statuses + validated step rendering.

2. **Preparation matrix smoke tests**
- Every preparation mode with representative ligands.
- Invalid combination rejection tests.
- Tool-missing graceful-failure tests.

3. **Docking execution tests**
- Local run preflight validation,
- parameter schema validation (basic vs advanced),
- GNINA GPU/CNN option validation.

4. **Post-docking and consensus tests**
- Analyze existing project without docking rerun,
- allscore normalization outputs,
- per-protein/global consensus + correlation outputs,
- complex structure generation validation.

5. **QC/biology/rescoring tests**
- ADMET and receptor QC gate behavior,
- external biology mapping and unresolved entity reports,
- OnionNet2 plugin capability-gated execution.

6. **Reproducibility tests**
- `.meta` manifest completeness,
- checkpoint lineage integrity,
- CSV vs SQLite parity checks (if SQLite enabled).

## Rollout and Change Management

- Rollout is milestone-gated with feature flags for high-risk changes:
  - `dockforge_task_manager`,
  - `dockforge_prep_graph`,
  - `dockforge_allscore_v2`,
  - `dockforge_sqlite_backend`.
- Default switches only after smoke/parity checks pass.
- Update user-facing guides (`README.md`, `POST_DOCKING_ANALYSIS_GUIDE.md`, `AUTODOCK_PREPARATION_GUIDE.md`, `HPC_DEPLOYMENT_GUIDE.md`) per milestone.

## Complexity Tracking

No constitution violations currently required. Complexity remains within existing single-repo architecture by enforcing modular boundaries and staged rollout.
