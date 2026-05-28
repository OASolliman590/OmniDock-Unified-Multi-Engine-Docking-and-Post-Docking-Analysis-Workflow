# Tasks: Unified Output Topology

**Input**: `/specs/021-unified-output-topology/spec.md`
**Scope**: Collapse fragmented post-docking output directories into one canonical `5-Analysis/` root.

## Format: `[ID] [P?] Description`
- **[P]**: parallelizable

## Phase 1: Output Root Switch

- [x] T001 Audit all write paths in `post_docking_analysis/` that reference `5-Post_Docking_Analysis` and collect into a migration list.
- [x] T002 Redirect session-level write paths to `5-Analysis/sessions/<session_id>/` in unified pipeline.
- [x] T003 Add compatibility shim: create symlink `5-Post_Docking_Analysis → 5-Analysis` when old dir does not exist; emit deprecation log notice.
- [x] T004 When `5-Post_Docking_Analysis/` already exists from prior run, log notice and leave it untouched.

## Phase 2: Hierarchy Enforcement

- [x] T005 [P] Enforce `scores/`, `best_poses/`, `complexes/`, `rmsd_analysis/` write into `sessions/<id>/`.
- [x] T006 [P] Enforce `interactions/` and `visualizations/` write into `sessions/<id>/`.
- [x] T007 Enforce `comparative/`, `top_pose_ligand_performance/`, `polypharmacology/`, `structure_quality/` write at `5-Analysis/` root (aggregate, not per-session).
- [x] T008 Update `LATEST_SESSION` symlink to point to most recently completed session after each run.

## Phase 3: Working Data Rename

- [x] T009 Rename `4-PostDocking/` → `4-Working/` across all pipeline write paths.
- [x] T010 Update any references to `4-PostDocking` in CLI help text, log messages, and docs.

## Phase 4: Navigation Document

- [x] T011 Generate `5-Analysis/START_HERE.md` in consolidation step with: session ID, key artifact paths, complexes path for MD users, top-pose path, link to `7-Reports/`.
- [x] T012 Include warning in `START_HERE.md` if `5-Post_Docking_Analysis/` exists with prior data.

## Phase 5: Verification

- [x] T013 Smoke test: after a run, assert zero files written to `5-Post_Docking_Analysis/`.
- [x] T014 Smoke test: assert `5-Analysis/START_HERE.md` exists and contains `complexes/`, `top_pose_ligand_performance/`, `7-Reports/` paths.
- [x] T015 Smoke test: assert `LATEST_SESSION` symlink resolves correctly.
- [x] T016 Smoke test: assert `4-Working/` created, not `4-PostDocking/`, on new runs.
- [x] T017 Grep assertion: `grep -r "5-Post_Docking_Analysis" post_docking_analysis/` returns only shim/deprecation notice lines.
