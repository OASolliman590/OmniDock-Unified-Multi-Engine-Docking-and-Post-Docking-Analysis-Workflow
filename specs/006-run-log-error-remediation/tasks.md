# Tasks: Simplified Pipeline Run-Log Error Remediation Closeout

**Input**: Design documents from `/specs/006-run-log-error-remediation/`
**Prerequisites**: `plan.md`, `spec.md`
**Status (2026-04-01)**: Closed as superseded by `018-dockforge-omnidock-platform` and `020-post-docking-scientific-foundation`.

**Tests**: Use `py_compile`, targeted regression proofs, and one real rerun of
 the simplified pipeline command.

## Phase 1: Audit Complete

- [x] T001 Confirm PoseView status normalization and terminal-state handling
- [x] T002 Confirm LigPlot validation and stage-level failure accounting
- [x] T003 Confirm ProLIF frame-network guardrails
- [x] T004 Confirm pairlist auto-discovery and logging paths

## Phase 2: Verification Completed

- [x] T005 Run `python -m py_compile` for the touched remediation modules
- [x] T006 Run a targeted PoseView completion-state regression proof
- [x] T007 Run a targeted deep-layout pairlist regression proof

## Phase 3: Remaining Proof

- [x] T008 Re-run the target simplified pipeline command on a representative GNINA dataset
- [x] T009 Record final warning/error distributions and external-tool caveats in `plan.md`
