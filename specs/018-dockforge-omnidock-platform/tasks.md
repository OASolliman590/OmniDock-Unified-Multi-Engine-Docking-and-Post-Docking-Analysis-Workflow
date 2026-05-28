# Tasks: Omni-DockForge: End-to-End Docking, Consensus by Design.

**Input**: Design documents from `/specs/018-dockforge-omnidock-platform/`  
**Prerequisites**: `spec.md`, `plan.md`  
**Organization**: Tasks are grouped by milestone and mapped to user stories and platform foundations.

## Status Snapshot (2026-03-31)

- Reconciled against current repository implementation and smoke tests.
- Completed: **106 / 106**
- Open: **0 / 106**

## Format: `[ID] [P?] [Area] Description`

- **[P]**: Parallelizable (different files/no hard dependency)
- **[Area]**: `FOUNDATION`, `M1`, `M2`, `M3`, `M4`, `M5`, `M6`

## Phase 1: Foundation and Spec Alignment (Blocking)

**Purpose**: Establish shared contracts before milestone work.

- [x] T001 [FOUNDATION] Create DockForge feature flag/config section in `workflow/models.py` and `workflow/state.py` for staged rollout toggles.
- [x] T002 [P] [FOUNDATION] Add canonical status enum (`not_started`, `in_progress`, `completed`, `validated`, `failed`, `skipped`, `needs_review`) in `workflow/models.py`.
- [x] T003 [P] [FOUNDATION] Add schema constants for canonical folder topology in `docking/project_layout.py`.
- [x] T004 [FOUNDATION] Create migration notes doc `docs/dockforge_migration_notes.md` (or root guide equivalent) describing old->new mapping.
- [x] T005 [FOUNDATION] Add shared utility for deterministic timestamp/run-id formatting in `workflow/state.py` or dedicated helper module.
- [x] T006 [FOUNDATION] Add scripted smoke entry scaffold `test/test_dockforge_smoke.py` for upcoming milestone checks.

**Blocker Rule**: No milestone coding starts until T001-T006 are complete.

---

## Milestone 1: Core Workflow and UX Scaffolding

**Goal**: non-blocking UX, timeline visibility, universal back navigation, Checkpoint & Revise naming.

### M1.1 Interactive Timeline + State Rendering

- [x] T007 [M1] Implement workflow timeline renderer in `workflow/interactive.py`.
- [x] T008 [P] [M1] Add timeline state query API in `workflow/state.py`.
- [x] T009 [P] [M1] Add validated-step visual rendering helper (green state) in `workflow/interactive.py`.
- [x] T010 [M1] Thread timeline display into all interactive stage menus in `interactive_pipeline.py`.

### M1.2 Universal Back Navigation

- [x] T011 [M1] Introduce shared back-navigation contract in `workflow/cli.py` and `workflow/interactive.py`.
- [x] T012 [M1] Update each interactive panel handler in `interactive_pipeline.py` to register `Back` action.
- [x] T013 [P] [M1] Add regression guard checks for panel backflow in `test/test_dockforge_smoke.py`.

### M1.3 Background Task Manager

- [x] T014 [M1] Add `BackgroundTask` model in `workflow/models.py`.
- [x] T015 [M1] Implement `TaskManager` lifecycle methods in `workflow/execution.py`.
- [x] T016 [P] [M1] Add session task-state persistence in `workflow/state.py`.
- [x] T017 [M1] Integrate background task launch hooks into long-running actions in `interactive_pipeline.py`.
- [x] T018 [P] [M1] Add task monitor panel in `workflow/interactive.py`.

### M1.4 Checkpoint & Revise Rename

- [x] T019 [M1] Replace `maturation` labels with `Checkpoint & Revise` in `interactive_pipeline.py`.
- [x] T020 [P] [M1] Update related CLI text/help in `main.py`, `cli_pipeline.py`, and `README.md`.
- [x] T021 [M1] Add checkpoint metadata model and serialization scaffold in `workflow/state.py`.

### M1 Verification

- [x] T022 [M1] Add smoke test: start long task, navigate away, verify task continues (`test/test_dockforge_smoke.py`).
- [x] T023 [M1] Add smoke test: universal back from every major panel (`test/test_dockforge_smoke.py`).

**Milestone 1 Blockers**: T007, T011, T014, T019 must complete before Milestone 2.

---

## Milestone 2: Ligand Preparation Refactor

**Goal**: engine-aware prep graph with Open Babel enrichment and PLIP removal from prep.

### M2.1 Remove PLIP from Preparation

- [x] T024 [M2] Remove PLIP prep prompt/options from `core_pipeline.py` and `interactive_pipeline.py`.
- [x] T025 [P] [M2] Ensure PLIP remains available only in post-docking modules (`post_docking_analysis/plip_integration.py`, related CLI prompts).
- [x] T026 [M2] Update preparation docs to remove PLIP mention in `AUTODOCK_PREPARATION_GUIDE.md` and `README.md`.

### M2.2 Open Babel Enrichment Pipeline

- [x] T027 [M2] Implement Open Babel enrichment step sequence in `docking/preparation/ligand_preparation.py`:
  - add hydrogens,
  - protonation at pH,
  - 3D conformer generation,
  - minimization,
  - partial charges.
- [x] T028 [P] [M2] Add pH default (7.4) and validation rules in `docking/models.py` and `docking/preparation/ligand_preparation.py`.
- [x] T029 [M2] Add per-ligand step report output writer in `docking/preparation/ligand_preparation.py`.

### M2.3 Engine-Aware Preparation Graph

- [x] T030 [M2] Define preparation mode enum and graph contracts in `docking/models.py`.
- [x] T031 [M2] Implement mode/engine compatibility validator in `docking/preparation/ligand_preparation.py`.
- [x] T032 [P] [M2] Add preparation menu redesign in `interactive_pipeline.py` and `workflow/interactive.py`.
- [x] T033 [M2] Implement `Engine-aware full preparation` orchestrator in `docking/preparation/ligand_preparation.py`.
- [x] T034 [M2] Add invalid-combination failure messaging in `docking/cli.py` and `interactive_pipeline.py`.

### M2.4 Preparation Output Contracts

- [x] T035 [M2] Define expected output artifacts per engine family in `docking/project_layout.py`.
- [x] T036 [P] [M2] Add preparation output validator for each mode in `docking/preparation/ligand_preparation.py`.
- [x] T037 [M2] Add preflight summary report generation (`CSV/JSON`) for preparation results.

### M2 Verification

- [x] T038 [M2] Add smoke tests for all preparation modes in `test/test_dockforge_smoke.py`.
- [x] T039 [M2] Add tool-missing graceful degradation tests for Open Babel/Meeko/ADT paths in `test/test_dockforge_smoke.py`.

**Milestone 2 Blockers**: T027, T030, T031, T033 complete before Milestone 3.

---

## Milestone 3: Docking Orchestration and Runtime Profiles

**Goal**: environment-aware run docking, parameter schemas, GPU-aware GNINA controls.

### M3.1 Execution Environment Abstraction

- [x] T040 [M3] Add execution environment models in `docking/models.py`.
- [x] T041 [M3] Implement adapter interface in `workflow/execution.py` or `docking/deployment.py`.
- [x] T042 [P] [M3] Implement local CPU adapter in `docking/deployment.py`.
- [x] T043 [P] [M3] Implement local GPU adapter skeleton in `docking/deployment.py`.
- [x] T044 [P] [M3] Implement conda/container adapter skeletons in `docking/deployment.py`.
- [x] T045 [P] [M3] Define remote/HPC adapter contract in `docking/remote_ops.py`.

### M3.2 Run Docking UX and Preflight

- [x] T046 [M3] Redesign `Run Docking` interactive flow for environment/folder selection in `interactive_pipeline.py`.
- [x] T047 [M3] Add prerequisites folder/workdir prompts in `workflow/interactive.py` and `interactive_pipeline.py`.
- [x] T048 [M3] Implement preflight required-input validation per engine in `docking/cli.py` and `docking/engine_registry.py`.
- [x] T049 [P] [M3] Add preflight report outputs under canonical structure in `docking/project_layout.py`.

### M3.3 Parameter Schema and Presets

- [x] T050 [M3] Define common docking parameter schema in `docking/models.py`.
- [x] T051 [M3] Add engine-specific extension schemas for Vina/Smina/GNINA/AD4 in `docking/models.py` and runner files.
- [x] T052 [P] [M3] Add `basic` vs `advanced` parameter mode UX in `interactive_pipeline.py`.
- [x] T053 [M3] Add screening/exhaustive default presets in `docking/engine_registry.py`.
- [x] T054 [M3] Implement GNINA CNN/GPU validation in `docking/runners/gnina.py`.
- [x] T055 [P] [M3] Add parameter serialization into run metadata files in `workflow/state.py`.

### M3 Verification

- [x] T056 [M3] Add smoke tests for environment selection and preflight behavior in `test/test_dockforge_smoke.py`.
- [x] T057 [M3] Add schema validation tests (basic/advanced + engine extension constraints).

**Milestone 3 Blockers**: T040, T046, T050, T051 before Milestone 4.

---

## Milestone 4: Post-Docking Promotion and AllScore Core

**Goal**: project-centric analysis entry, normalization/consensus, complex export reliability.

### M4.1 Project-Centric Analysis Entry

- [x] T058 [M4] Add analysis-only project entry path in `interactive_pipeline.py`.
- [x] T059 [P] [M4] Add modular analysis scope selector (comparison/rescoring/QC/report/full) in `post_docking_analysis/cli.py`.
- [x] T060 [M4] Persist analysis session state and reopen support in `workflow/state.py` and `post_docking_analysis/pipeline.py`.

### M4.2 AllScore Ingestion and Normalization

- [x] T061 [M4] Create unified score ingestion layer in `post_docking_analysis/multi_engine_pipeline.py`.
- [x] T062 [M4] Add normalization strategy module (configurable) in `post_docking_analysis/consensus.py`.
- [x] T063 [P] [M4] Emit raw/normalized/consensus outputs under `4-PostDocking/scores/{raw,unified,consensus}`.
- [x] T064 [M4] Add consensus explainability payload expansion in `post_docking_analysis/consensus.py`.

### M4.3 Correlation and Complex Export

- [x] T065 [M4] Add per-protein/global engine rank correlation outputs in `post_docking_analysis/correlation_analyzer.py`.
- [x] T066 [M4] Harden protein+pose complex builder validation in `post_docking_analysis/structure_processor.py`.
- [x] T067 [P] [M4] Add complex export indexing for downstream rescoring/visualization in `post_docking_analysis/pipeline.py`.

### M4 Verification

- [x] T068 [M4] Add smoke test: analyze existing project without rerunning docking.
- [x] T069 [M4] Add smoke test: allscore output contract generation and explainability artifacts.
- [x] T070 [M4] Add smoke test: complex export integrity checks (chain/residue/pose validation).

**Milestone 4 Blockers**: T061, T062, T065 before Milestone 5.

---

## Milestone 5: QC, Biology Integration, and Rescoring Plugins

**Goal**: stronger scientific gates and external-correlation capabilities.

### M5.1 Upstream QC Gates

- [x] T071 [M5] Add receptor QC gate before docking in `core_pipeline.py` and `post_docking_analysis/structure_quality.py` (or preparation QC module).
- [x] T072 [M5] Add upstream ligand QC/ADMET filter integration in `docking/preparation/ligand_quality.py`.
- [x] T073 [P] [M5] Add configurable QC thresholds to config schema in `post_docking_analysis/config/schema.yaml` and relevant workflow config loader.

### M5.2 Benchmarking and Classification

- [x] T074 [M5] Implement optional decoy/reference benchmarking (EF1%, EF5%, ROC-AUC) in `post_docking_analysis/plugins/enrichment_plugin.py`.
- [x] T075 [M5] Implement target-aware strong/moderate/weak classifier in `post_docking_analysis/consensus.py` or dedicated module.
- [x] T076 [P] [M5] Add classifier configuration + report outputs in `post_docking_analysis/report_generator.py`.

### M5.3 External Biology Correlation

- [x] T077 [M5] Add biology file ingestion (CSV/TSV/JSON) in new module `post_docking_analysis/biology_integration.py`.
- [x] T078 [M5] Implement entity mapping and unresolved mapping report generation.
- [x] T079 [P] [M5] Add docking-biology correlation outputs to `5-Analysis/polypharmacology/` and report summaries.

### M5.4 OnionNet2 Plugin Point

- [x] T080 [M5] Define rescoring plugin interface update in `post_docking_analysis/plugin_manager.py`.
- [x] T081 [M5] Add OnionNet2 plugin scaffold in `post_docking_analysis/plugins/onionnet2_plugin.py`.
- [x] T082 [P] [M5] Add capability-gated execution + graceful fallback messaging for missing OnionNet2 binary/env.

### M5 Verification

- [x] T083 [M5] Add smoke tests for QC gate pass/fail behavior.
- [x] T084 [M5] Add smoke tests for classification policy outputs (target-aware, configurable).
- [x] T085 [M5] Add smoke tests for biology mapping unresolved-entity reporting.

**Milestone 5 Blockers**: Cleared (T071, T072, T075, T080 completed).

---

## Milestone 6: Reproducibility, Storage, and Reporting Finalization

**Goal**: auditable manifests/session durability, optional SQLite backend, extensibility hooks.

### M6.1 `.meta` Provenance and Checkpoint Lineage

- [x] T086 [M6] Implement `.meta/config.yaml` freeze writer in `workflow/state.py`.
- [x] T087 [M6] Implement `.meta/run_manifest.json` with checksums/tool versions/timestamps in `workflow/state.py`.
- [x] T088 [P] [M6] Implement `.meta/env.lock` snapshot generation (graceful fallback if unavailable).
- [x] T089 [M6] Implement checkpoint lineage records for `Checkpoint & Revise` in `workflow/state.py`.

### M6.2 SQLite Backend (Optional, Parity-Gated)

- [x] T090 [M6] Add SQLite schema module `post_docking_analysis/storage_sqlite.py`.
- [x] T091 [M6] Implement dual-write pipeline (CSV + SQLite) behind feature flag.
- [x] T092 [P] [M6] Add CSV-vs-SQLite parity validator command/script in `test/test_dockforge_smoke.py` or dedicated tool.

### M6.3 Reporting and Canonical Structure Consolidation

- [x] T093 [M6] Enforce canonical numbered output folders in `docking/project_layout.py`.
- [x] T094 [M6] Add final report index (`7-Reports/START_HERE.md`) generation in `post_docking_analysis/report_generator.py`.
- [x] T095 [P] [M6] Add consolidated run summary across engines/QC/classification in report outputs.

### M6.4 Future Extension Hooks

- [x] T096 [M6] Add ensemble receptor extension interface stubs in `docking/models.py` and `docking/runners/base.py`.
- [x] T097 [P] [M6] Add rerun-manifest auto-execution interface stub in `docking/deployment.py`.
- [x] T098 [P] [M6] Add dashboard-ready export contract doc and JSON index emitters in `post_docking_analysis/report_generator.py`.

### M6 Verification

- [x] T099 [M6] Smoke test: checkpoint -> revise -> rerun lineage integrity.
- [x] T100 [M6] Smoke test: manifest completeness and deterministic output paths.
- [x] T101 [M6] Smoke test: SQLite parity checks (if enabled) or explicit skip with rationale.

---

## Cross-Cutting Documentation and Release Tasks

- [x] T102 [P] [M6] Update `README.md` with DockForge naming and workflow map.
- [x] T103 [P] [M6] Update `POST_DOCKING_ANALYSIS_GUIDE.md` for project-centric analysis entry.
- [x] T104 [P] [M6] Update `AUTODOCK_PREPARATION_GUIDE.md` with engine-aware preparation graph and Open Babel enrichment.
- [x] T105 [P] [M6] Update `HPC_DEPLOYMENT_GUIDE.md` to align with environment abstraction and preflight behavior.
- [x] T106 [M6] Add release checklist doc for feature-flag rollout and fallback strategy.

---

## Dependencies & Execution Order

### Milestone Dependencies

- Foundation -> M1 -> M2 -> M3 -> M4 -> M5 -> M6
- M4 can begin partially once M3 output contracts stabilize.
- M5 requires allscore normalization outputs from M4.
- M6 requires stable outputs from M1-M5.

### Parallelization Opportunities

- Within each milestone, tasks marked `[P]` are parallel-safe.
- Suggested parallel teams:
  - Team A: interactive/task manager/state,
  - Team B: preparation/runtime adapters,
  - Team C: post-docking consensus/QC/reporting.

### Explicit Blockers

- No M2 without M1 task manager/state timeline primitives.
- No M3 without M2 preparation compatibility validator.
- No M4 without M3 parameter and output contracts.
- No M5 without M4 allscore ingestion/normalization core.
- No M6 default switch without parity + smoke validation.

---

## Implementation Strategy

### MVP-First Path

1. Complete Foundation + Milestone 1.
2. Complete Milestone 2 + Milestone 3 to ensure reliable preparation/execution.
3. Deliver Milestone 4 for first-class post-docking and consensus baseline.
4. Add Milestones 5 and 6 for scientific rigor and reproducibility hardening.

### Validation Gates

- Gate A (after M1): Non-blocking UX + back navigation proven.
- Gate B (after M3): Engine-aware prep and docking preflight proven.
- Gate C (after M4): Analysis-only entry + allscore outputs proven.
- Gate D (after M6): Reproducibility manifests + docs + parity checks proven.
