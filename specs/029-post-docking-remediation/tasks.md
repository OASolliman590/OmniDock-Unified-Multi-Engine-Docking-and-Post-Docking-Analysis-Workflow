# Tasks: Post-Docking Analysis Remediation & Contract Uniformity

**Input**: Design documents from `/specs/029-post-docking-remediation/`
**Prerequisites**: `spec.md`, `plan.md`

## Phase 1: Spec Setup

- [x] T001 Materialize remediation spec under `specs/029-post-docking-remediation/`
- [x] T002 Create implementation plan and task tracker for this cycle

## Phase 2: P1 Contract Enforcement

- [x] T003 FR-001 Enforce `run_rmsd=True` / `run_visualizations=True` for non-GNINA simplified pipeline construction in `workflow/execution.py`
- [x] T004 FR-002 Route standalone interaction aliases (`prolif` / `ligplot` / `pandamap`) through `analyze.interactions.clean` across interactive + non-interactive dispatch
- [x] T005 FR-003 Block legacy `PostDockingAnalysisPipeline` stage routes (`analyze.stage.structure_quality`, `analyze.visuals.pymol`) with migration guidance and smoke coverage

## Phase 3: P1 Engine Parity

- [x] T006 FR-004 Extend `post_docking_analysis/engine_detector.py` with AutoDock4 detection/details parity and score-schema reporting
- [x] T007 FR-004 Extend multi-engine normalization bridge in `post_docking_analysis/_legacy_multi_engine_pipeline_impl.py` to parse AutoDock4 `.dlg` outputs
- [x] T008 FR-004 Add smoke coverage in `test/test_dockforge_smoke.py` asserting four-engine detection parity (`gnina/vina/smina/autodock4`)

## Phase 4: Next P1 Slice

- [x] T009 FR-005 Design and implement per-engine HPC adapter parity (`vina`, `smina`, `autodock4`) analogous to `gnina_hpc_adapter.py`
- [x] T010 FR-006 Assert run-tracking schema consistency across all engines

## Phase 5: P2 Observability Follow-Through

- [x] T011 FR-007 Surface run-tracking manifest path + one-line stage-contract diff in interactive analysis completion output with smoke coverage
- [x] T012 FR-008 Normalize optional-feature manifest taxonomy (`completed` / `skipped_disabled` / `skipped_missing_dependency` / `failed_error`) with smoke coverage

## Phase 6: P2 Favorite-Flow Alignment

- [x] T013 FR-009 Surface `analyze.interactions.clean` inside favorite-engine continuation flow and validate with interactive smoke coverage

## Phase 7: P2 Legacy-Naming Cleanup

- [x] T014 FR-010 Rename `_legacy_*_impl.py` files to canonical `*_pipeline_impl.py`, collapse wrapper wording, and rewire smoke imports

## Phase 8: P2 Dead-Code Removal

- [x] T015 FR-011 Remove `rmsd_analyzer.py`, `affinity_analyzer.py`, `plugin_manager.py`, and `plugins/`, with migration rewires + retirement smoke coverage
