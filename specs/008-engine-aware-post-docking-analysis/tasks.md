# Tasks: Engine-Aware Post-Docking Analysis

**Input**: Design documents from `/specs/008-engine-aware-post-docking-analysis/`
**Prerequisites**: `plan.md`, `spec.md`

**Tests**: Use `py_compile`, CLI help checks, and canonical project smoke runs.

## Phase 1: Setup

- [x] T001 Create a normalized multi-engine score schema in `post_docking_analysis/multi_engine_pipeline.py`
- [x] T002 Add canonical project score discovery for GNINA, Vina, and Smina

## Phase 2: User Story 1 - Comparative Reports

- [x] T003 [US1] Implement comparative report generation in `post_docking_analysis/multi_engine_pipeline.py`
- [x] T004 [US1] Emit per-engine and cross-engine best-pose summaries
- [x] T005 [US1] Emit unified score tables with compatibility columns for all engines

## Phase 3: User Story 2 - Single/Favorite Engine Filtering

- [x] T006 [US2] Implement `single_engine` filtering in `post_docking_analysis/multi_engine_pipeline.py`
- [x] T007 [US2] Implement `favorite_engine_continue` filtering and outputs
- [x] T007a [US2] Bridge GNINA favorites into the simplified structural pipeline
- [x] T007b [US2] Bridge Vina/Smina favorites into best-pose PDB extraction plus continued downstream report/quality outputs
- [x] T007c [US2] Reuse the safe subset of the simplified downstream stack for Vina/Smina favorites
- [x] T007d [US2] Add a non-GNINA PDB-based RMSD bridge for best-pose and per-complex scopes where comparable structures exist
- [x] T007e [US2] Emit explicit bridge stage-status notes for optional and partially available stages

## Phase 4: User Story 3 - CLI Integration

- [x] T008 [US3] Extend `post_docking_analysis/cli.py` with engine-aware arguments
- [x] T009 [US3] Preserve legacy non-project analysis behavior

## Phase 5: Validation

- [x] T010 Run `python -m py_compile` for touched post-docking modules
- [x] T011 Run `python -m post_docking_analysis --help`
- [x] T012 Run engine-aware analysis on a canonical project
- [x] T013 Run favorite-engine continuation smoke coverage for at least one PDBQT-based engine
