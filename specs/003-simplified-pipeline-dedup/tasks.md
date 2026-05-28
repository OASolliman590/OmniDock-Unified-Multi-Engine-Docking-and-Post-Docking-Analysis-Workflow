# Tasks: Simplified Pipeline Deduplication

**Input**: Design documents from `/specs/003-simplified-pipeline-dedup/`
**Prerequisites**: `plan.md`, `spec.md`
**Status (2026-04-01)**: Closed as superseded by `018-dockforge-omnidock-platform` and `020-post-docking-scientific-foundation`. Remaining historical checklist items were reconciled into newer specs.

**Tests**: Validation uses `python -m py_compile` and a simplified CLI smoke
run in the `pdb-prepare-wizard` environment.

**Organization**: These tasks target the parts of the deduplication spec that
still appear incomplete in the current code.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel
- **[Story]**: Workstream label (`US1` or `US2`)

## Phase 1: Setup

**Purpose**: Confirm the remaining deduplication targets

- [x] T001 Capture all target stages listed in `specs/003-simplified-pipeline-dedup/spec.md`
- [x] T002 Note the verification commands and expected smoke path in `specs/003-simplified-pipeline-dedup/plan.md`

---

## Phase 2: Foundational Review

**Purpose**: Establish the current helper contracts before refactoring

- [x] T003 Audit `_best_poses_from_scores()` in `post_docking_analysis/simplified_pipeline.py`
- [x] T004 [P] Audit `_get_complex_pdb_files()` usage in `post_docking_analysis/simplified_pipeline.py`

---

## Phase 3: User Story 1 - Centralize Best-Pose Selection

**Goal**: Use one helper for all best-pose selection logic

**Independent Test**: Confirm the RMSD, report, py3Dmol, ProLIF, LigPlot, and
PandaMap stages all resolve best poses through the helper.

### Implementation for User Story 1

- [x] T005 [US1] Finalize the shared contract for `_best_poses_from_scores()` in `post_docking_analysis/simplified_pipeline.py`
- [x] T006 [P] [US1] Replace remaining best-pose selection in RMSD and report stages within `post_docking_analysis/simplified_pipeline.py`
- [x] T007 [P] [US1] Replace remaining best-pose selection in py3Dmol and ProLIF stages within `post_docking_analysis/simplified_pipeline.py`
- [x] T008 [P] [US1] Replace remaining best-pose selection in LigPlot and PandaMap stages within `post_docking_analysis/simplified_pipeline.py`

---

## Phase 4: User Story 2 - Centralize Complex PDB Discovery

**Goal**: Use one helper for generated complex file discovery

**Independent Test**: Confirm all downstream stages resolve complex PDB files
through `_get_complex_pdb_files()`.

### Implementation for User Story 2

- [x] T009 [US2] Finalize the shared contract for `_get_complex_pdb_files()` in `post_docking_analysis/simplified_pipeline.py`
- [x] T010 [P] [US2] Replace remaining inline complex-file discovery in extraction and report paths within `post_docking_analysis/simplified_pipeline.py`
- [x] T011 [P] [US2] Replace remaining inline complex-file discovery in py3Dmol and ProLIF paths within `post_docking_analysis/simplified_pipeline.py`
- [x] T012 [P] [US2] Replace remaining inline complex-file discovery in LigPlot and PandaMap paths within `post_docking_analysis/simplified_pipeline.py`

---

## Phase 5: Validation & Closeout

**Purpose**: Prove behavior stayed unchanged

- [x] T013 [P] Run `python -m py_compile post_docking_analysis/simplified_pipeline.py`
- [x] T014 [P] Run the simplified CLI smoke command from `specs/003-simplified-pipeline-dedup/plan.md`
- [x] T015 [P] Record unchanged successful-path and missing-input-path behavior in `specs/003-simplified-pipeline-dedup/plan.md`
