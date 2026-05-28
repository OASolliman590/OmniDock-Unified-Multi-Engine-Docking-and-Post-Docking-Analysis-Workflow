# Tasks: Visualization + RMSD Overhaul Closeout

**Input**: Design documents from `/specs/004-visualization-rmsd-overhaul/`
**Prerequisites**: `plan.md`, `spec.md`
**Status (2026-04-01)**: Closed as superseded by `018-dockforge-omnidock-platform` and `020-post-docking-scientific-foundation`.

**Tests**: Use `py_compile`, simplified CLI help checks, and one simplified
pipeline smoke run against the documented output tree.

## Phase 1: Audit Complete

- [x] T001 Confirm the current label-formatting contract in simplified and hierarchical outputs
- [x] T002 Confirm the current RMSD output tree and scope split
- [x] T003 Confirm canonical interaction output routing under `interactions/`
- [x] T004 Align user-facing docs with the current CLI and visualization behavior

## Phase 2: Verification Completed

- [x] T005 Run `python -m py_compile` for the touched visualization/RMSD modules
- [x] T006 Recheck the simplified CLI help surface for `--prompt-protein-names` and `--enable-poseview`

## Phase 3: Remaining Proof

- [x] T007 Run one simplified pipeline smoke command and inspect the output tree against the spec
- [x] T008 Record final optional-tool caveats and close the spec
