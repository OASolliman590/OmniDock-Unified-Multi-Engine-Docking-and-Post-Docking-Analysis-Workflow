# Tasks: Docking Preparation Automation Closeout

**Input**: Design documents from `/specs/009-docking-preparation-automation/`
**Prerequisites**: `plan.md`, `spec.md`

**Tests**: Use `py_compile`, grouped help checks, and a synthetic staged
project smoke flow. Keep live `pdb collect` verification as the remaining proof
gate.

## Phase 1: Audit Current Implementation

- [x] T001 Confirm staged layout and manifest behavior in `docking/project_layout.py` and `workflow/execution.py`
- [x] T002 Confirm Excel summary parsing, pairlist building, and project materialization in `docking/preparation/`

## Phase 2: Verification Completed

- [x] T003 Run `python -m py_compile main.py docking/*.py docking/preparation/*.py docking/runners/*.py workflow/*.py post_docking_analysis/multi_engine_pipeline.py`
- [x] T004 Run `python main.py workflow init --help`, `python main.py prep pairlist --help`, and `python main.py prep project --help`
- [x] T005 Run a synthetic staged-project smoke flow through `workflow init`, `prep pairlist`, and `prep project`
- [x] T006 Update the spec and plan wording to match the grouped prep CLI and legacy `prepare-docking` alias behavior

## Phase 3: Remaining Proof

- [x] T007 Revalidate `python main.py pdb collect --project-dir <dir> --pdbs ...` on a live structure fetch/workbook update path
- [x] T008 Record final `pdb collect` evidence, Open Babel status, and workbook-fix note in `specs/009-docking-preparation-automation/plan.md`
