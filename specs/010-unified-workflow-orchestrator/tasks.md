# Tasks: Unified Workflow Orchestrator

**Input**: Design documents from `/specs/010-unified-workflow-orchestrator/`
**Prerequisites**: `plan.md`, `spec.md`

**Tests**: Use `py_compile`, grouped CLI help checks, and legacy alias help
checks.

## Phase 1: Foundations

- [x] T001 Create the `workflow/` package with models and state persistence
- [x] T002 Add grouped CLI parsing for `workflow`, `pdb`, `prep`, `dock`, and `analyze`
- [x] T003 Replace `main.py` with the grouped router and legacy alias normalization

## Phase 2: Workflow Shell

- [x] T004 Implement the guided interactive workflow shell
- [x] T005 Add workflow state summaries and resume support
- [x] T006 Add jump-to-function support in the interactive shell

## Phase 3: Execution Adapters

- [x] T007 Wrap existing PDB preparation surfaces for grouped CLI use
- [x] T008 Wrap docking project preparation and docking execution for grouped CLI use
- [x] T009 Add high-level comparative and favorite-engine analysis wrappers
- [x] T010 Add stage-level post-docking wrappers for RMSD, reports, visualizations, interactions, structure quality, and PyMOL

## Phase 4: Documentation

- [x] T011 Update the root README with grouped workflow examples
- [x] T012 Add Spec Kit coverage for the unified orchestrator
- [x] T017 Add `workflow init` and GNINA HPC-compatibility coverage to docs/specs

## Phase 5: Validation

- [x] T013 Run `python -m py_compile main.py workflow/cli.py workflow/interactive.py workflow/execution.py workflow/state.py workflow/models.py`
- [x] T014 Run grouped CLI help checks
- [x] T015 Run legacy alias help checks
- [x] T016 Add and validate `python main.py workflow jump --help`
- [x] T018 Validate `workflow init` creates a GNINA-compatible scaffold and workflow state
