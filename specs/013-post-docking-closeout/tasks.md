# Tasks: Post-Docking Closeout

**Input**: Design documents from `/specs/013-post-docking-closeout/`
**Prerequisites**: `plan.md`, `spec.md`

**Tests**: Use `py_compile`, CLI help checks, non-TTY interactive-route checks,
 and a targeted PoseView completion-state regression proof.

## Phase 1: Setup

- [x] T001 Create the numbered branch `013-post-docking-closeout`
- [x] T002 Scaffold `spec.md`, `plan.md`, and `tasks.md`

## Phase 2: Entry-Point Usability

- [x] T003 Make `python main.py` no-argument behavior safe for non-TTY sessions
- [x] T004 Lazy-load heavy workflow and simplified post-docking imports for help surfaces
- [x] T005 Add clear non-TTY guards for interactive-only workflow and PDB routes

## Phase 3: Run-Time Closeout

- [x] T006 Fix PoseView terminal success handling to trust normalized completion state
- [x] T007 Decide and implement the final pairlist auto-discovery contract
- [x] T008 Align simplified post-docking docs with the current CLI and visualization behavior
- [x] T009 Rewrite stale closeout tracking in specs `004`, `005`, and `006`

## Phase 4: Validation

- [x] T010 Run `python -m py_compile main.py workflow/cli.py post_docking_analysis/simplified_cli.py post_docking_analysis/poseview_integration.py`
- [x] T011 Run help and non-TTY guard checks for the updated entry points
- [x] T012 Run a targeted PoseView completion-state regression proof
- [x] T012a Run a targeted deep-layout pairlist regression proof
- [x] T013 Re-run any additional validations required by T007 to T009
