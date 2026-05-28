# Tasks: Interactive Simplified CLI Prompts Closeout

**Input**: Design documents from `/specs/005-interactive-cli-prompts/`
**Prerequisites**: `plan.md`, `spec.md`

**Tests**: Use `py_compile`, simplified CLI help checks, and non-TTY route
checks for interactive commands.

## Phase 1: Audit Complete

- [x] T001 Record the interactive and non-interactive validation matrix in `plan.md`
- [x] T002 Confirm the guided prompt flow and fallback behavior in the current CLI/pipeline code

## Phase 2: Verification Completed

- [x] T003 Confirm guided prompting in `post_docking_analysis/simplified_cli.py`
- [x] T004 Confirm per-target protein naming prompts in `post_docking_analysis/simplified_pipeline.py`
- [x] T005 Confirm non-TTY runs fail cleanly instead of prompting
- [x] T006 Confirm explicit help and fully specified invocation surfaces remain intact

## Phase 3: Closeout

- [x] T007 Record final verification outcomes in `specs/005-interactive-cli-prompts/plan.md`
