# Tasks: Pipeline Unification and Deduplication

**Input**: Design documents from `/specs/002-pipeline-unification/`
**Prerequisites**: `plan.md`, `spec.md`
**Status (2026-04-01)**: Closed as superseded by `018-dockforge-omnidock-platform` and `020-post-docking-scientific-foundation`. Remaining historical checklist items were reconciled into newer specs.

**Tests**: Use the repo's command-based validation path: `py_compile` plus one
targeted simplified pipeline smoke run.

**Organization**: These tasks reflect the remaining audit and closeout work for
this feature based on the current repository state.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel
- **[Story]**: User story label (`US1` to `US3`)

## Phase 1: Setup

**Purpose**: Lock down the verification scope for the refactor

- [x] T001 Record the validation commands and expected outcomes in `specs/002-pipeline-unification/plan.md`
- [x] T002 Capture the touched-file audit scope in `specs/002-pipeline-unification/spec.md`

---

## Phase 2: Foundational Review

**Purpose**: Establish a clear baseline before sign-off

- [x] T003 Audit current PandaMap runner usage across `post_docking_analysis/pandamap_runner.py`, `post_docking_analysis/pandamap_integration.py`, and `post_docking_analysis/publication_pandamap.py`
- [x] T004 [P] Audit the current canonical pipeline execution path in `post_docking_analysis/pipeline.py`
- [x] T005 [P] Audit remaining duplication or import ambiguity in `post_docking_analysis/simplified_pipeline.py`

---

## Phase 3: User Story 1 - Consistent PandaMap Runtime (Priority: P1)

**Goal**: Confirm both PandaMap analyzers use one shared runtime path

**Independent Test**: Verify both analyzers resolve CLI mode and command
execution through the shared runner.

### Implementation for User Story 1

- [x] T006 [US1] Consolidate any remaining PandaMap CLI detection logic in `post_docking_analysis/pandamap_runner.py`
- [x] T007 [P] [US1] Remove any remaining runner drift from `post_docking_analysis/pandamap_integration.py`
- [x] T008 [P] [US1] Remove any remaining runner drift from `post_docking_analysis/publication_pandamap.py`

---

## Phase 4: User Story 2 - Remove Dead Pipeline Flow (Priority: P1)

**Goal**: Confirm `pipeline.py` has one authoritative execution path

**Independent Test**: Inspect `_generate_all_scores_csv()` and run a GNINA
fast-path smoke command without behavior regression.

### Implementation for User Story 2

- [x] T009 [US2] Remove any residual unreachable legacy block from `post_docking_analysis/pipeline.py`
- [x] T010 [US2] Verify the GNINA fast-path behavior remains unchanged in `post_docking_analysis/pipeline.py`

---

## Phase 5: User Story 3 - Eliminate Local Duplication and Ambiguity (Priority: P2)

**Goal**: Remove remaining local duplication in the simplified pipeline

**Independent Test**: Import the module cleanly and confirm naming and output
organization still use one shared helper path.

### Implementation for User Story 3

- [x] T011 [US3] Remove any duplicate or overwritten imports in `post_docking_analysis/simplified_pipeline.py`
- [x] T012 [US3] Consolidate any remaining small helper duplication in `post_docking_analysis/simplified_pipeline.py`

---

## Phase 6: Validation & Closeout

**Purpose**: Prove the refactor is complete

- [x] T013 [P] Run `python -m py_compile` for `post_docking_analysis/pandamap_runner.py`, `post_docking_analysis/pandamap_integration.py`, `post_docking_analysis/publication_pandamap.py`, `post_docking_analysis/pipeline.py`, and `post_docking_analysis/simplified_pipeline.py`
- [x] T014 [P] Run one simplified pipeline smoke command and record the outcome in `specs/002-pipeline-unification/plan.md`
- [x] T015 [P] Update any refactor notes in `post_docking_analysis/README.md`
