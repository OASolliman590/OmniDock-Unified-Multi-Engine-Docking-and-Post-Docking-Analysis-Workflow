# Tasks: Unified Visualization Overhaul

**Input**: Design documents from `/specs/001-unify-and-optimize/`
**Prerequisites**: `plan.md`, `spec.md`
**Status (2026-04-01)**: Closed as superseded by `018-dockforge-omnidock-platform` and `020-post-docking-scientific-foundation`. Remaining historical checklist items were reconciled into newer specs.

**Tests**: Validation is command-based for this repository. Use syntax checks and
smoke runs captured in `specs/001-unify-and-optimize/plan.md`.

**Organization**: Tasks are grouped by user story so each outcome can be
verified independently.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no blocking dependency)
- **[Story]**: User story label (`US1` to `US4`)

## Phase 1: Setup (Shared Context)

**Purpose**: Lock the feature scope and verification path before implementation

- [x] T001 Confirm touched modules and validation commands in `specs/001-unify-and-optimize/plan.md`
- [x] T002 Capture dataset assumptions and expected output folders in `specs/001-unify-and-optimize/spec.md`

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Shared utilities and safety rails required by all user stories

- [x] T003 Create shared receptor and protein naming helpers in `post_docking_analysis/protein_naming.py`
- [x] T004 [P] Add output organization and manifest helpers in `post_docking_analysis/simplified_pipeline.py`
- [x] T005 [P] Add reusable protein-label formatting hooks in `post_docking_analysis/hierarchical_analyzer.py`
- [x] T006 [P] Add optional-tool guardrails for PandaMap execution in `post_docking_analysis/publication_pandamap.py`
- [x] T007 [P] Add optional-tool guardrails for interaction-map generation in `post_docking_analysis/prolif_interaction_maps.py`
- [x] T008 [P] Add optional-tool guardrails for LigPlot execution in `post_docking_analysis/ligplot_integration.py`

**Checkpoint**: Shared naming, output routing, and optional-tool behavior are
ready for feature work

---

## Phase 3: User Story 1 - Protein Naming In Visuals (Priority: P1) 🎯 MVP

**Goal**: Use human-readable protein labels across generated visual outputs

**Independent Test**: Run the simplified pipeline on receptors named by PDB code
and verify figure labels and manifests use display names.

### Implementation for User Story 1

- [x] T009 [US1] Persist receptor-to-display-name mappings in `post_docking_analysis/protein_naming.py`
- [x] T010 [P] [US1] Apply protein display labels to simplified visual outputs in `post_docking_analysis/simplified_pipeline.py`
- [x] T011 [P] [US1] Apply protein display labels to hierarchical outputs in `post_docking_analysis/hierarchical_analyzer.py`
- [x] T012 [US1] Document naming behavior and fallback rules in `POST_DOCKING_ANALYSIS_GUIDE.md`

**Checkpoint**: Display names are visible in figures and exported metadata

---

## Phase 4: User Story 2 - Better Interaction Visual Quality (Priority: P1)

**Goal**: Improve PandaMap, ProLIF, and LigPlot defaults for publication-ready
interaction outputs

**Independent Test**: Generate PandaMap, ProLIF, and LigPlot artifacts and
verify configured defaults produce non-empty outputs without regressing batch
behavior.

### Implementation for User Story 2

- [x] T013 [US2] Tune PandaMap shape and surface defaults in `post_docking_analysis/publication_pandamap.py`
- [x] T014 [P] [US2] Tune ProLIF batch defaults and output behavior in `post_docking_analysis/prolif_interaction_maps.py`
- [x] T015 [P] [US2] Tune LigPlot defaults and robust execution behavior in `post_docking_analysis/ligplot_integration.py`
- [x] T016 [US2] Wire interaction-quality settings through the main pipeline in `post_docking_analysis/simplified_pipeline.py`

**Checkpoint**: Interaction outputs reflect the new defaults and still degrade
gracefully when tools are unavailable

---

## Phase 5: User Story 3 - Correct and Readable RMSD Figures (Priority: P1)

**Goal**: Make RMSD outputs readable and robust for valid and edge-case inputs

**Independent Test**: Run the RMSD stage on both populated and sparse datasets
and verify plots render without label collisions or crash paths.

### Implementation for User Story 3

- [x] T017 [US3] Refine RMSD matrix preparation and scope handling in `post_docking_analysis/enhanced_rmsd_analyzer.py`
- [x] T018 [P] [US3] Improve RMSD figure labels, ticks, and layout in `post_docking_analysis/enhanced_rmsd_analyzer.py`
- [x] T019 [US3] Add sparse-data and degenerate-matrix fallbacks in `post_docking_analysis/simplified_pipeline.py`

**Checkpoint**: RMSD stage completes cleanly with readable output on target
datasets

---

## Phase 6: User Story 4 - Unified Output Topology (Priority: P1)

**Goal**: Consolidate visuals and raw artifacts into deterministic folders with
a linkage manifest

**Independent Test**: Run the pipeline and verify `visualizations/`,
`raw_data/`, and the manifest are populated consistently.

### Implementation for User Story 4

- [x] T020 [US4] Consolidate rendered figures under `visualizations/` in `post_docking_analysis/simplified_pipeline.py`
- [x] T021 [P] [US4] Consolidate raw CSV, JSON, TXT, and log outputs under `raw_data/` in `post_docking_analysis/simplified_pipeline.py`
- [x] T022 [US4] Emit the visualization linkage manifest in `post_docking_analysis/simplified_pipeline.py`
- [x] T023 [US4] Document output topology and manifest usage in `POST_DOCKING_ANALYSIS_GUIDE.md`

**Checkpoint**: Users can navigate one visualization tree and one raw-data tree
with an auditable manifest

---

## Phase 7: Polish & Validation

**Purpose**: Close the feature with repo-standard validation and documentation

- [x] T024 [P] Run syntax validation commands listed in `specs/001-unify-and-optimize/plan.md`
- [x] T025 [P] Run smoke validation commands and record outcome notes in `specs/001-unify-and-optimize/plan.md`
- [x] T026 [P] Update user-facing guide references for the final output layout in `post_docking_analysis/README.md`

---

## Dependencies & Execution Order

### Phase Dependencies

- Phase 1 has no dependencies.
- Phase 2 depends on Phase 1 and blocks all user stories.
- Phases 3 to 6 depend on Phase 2 and can proceed in priority order.
- Phase 7 depends on the desired user stories being complete.

### User Story Dependencies

- **US1** starts immediately after Phase 2 and provides the MVP.
- **US2** depends on the shared optional-tool guardrails from Phase 2.
- **US3** depends on the shared naming and pipeline helpers from Phase 2.
- **US4** depends on the output organization helpers from Phase 2 and should be
  finalized after the main visualization changes land.

### Parallel Opportunities

- T004 to T008 can run in parallel after T003.
- T010 and T011 can run in parallel.
- T014 and T015 can run in parallel.
- T018 can run in parallel with T019 once T017 defines the RMSD scope.
- T021 and T022 can run in parallel once T020 establishes the output layout.

---

## Implementation Strategy

### MVP First

1. Complete Phases 1 and 2.
2. Complete Phase 3.
3. Validate protein display naming end to end before moving on.

### Incremental Delivery

1. Land naming changes first.
2. Add interaction quality improvements.
3. Finish RMSD readability work.
4. Finalize output consolidation and documentation.
5. Close with syntax and smoke validation.

---

## Notes

- This repository uses smoke-style validation rather than a full `pytest` suite.
- Keep optional external tools non-fatal and document degraded behavior.
- Do not commit generated docking, visualization, or analysis output files.
