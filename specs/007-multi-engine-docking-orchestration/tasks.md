# Tasks: Multi-Engine Docking Orchestration Closeout

**Input**: Design documents from `/specs/007-multi-engine-docking-orchestration/`
**Prerequisites**: `plan.md`, `spec.md`

**Tests**: Use `py_compile`, grouped/legacy help checks, and a synthetic staged
project `dock run --dry-run` smoke flow. Keep score-normalization reruns as the
remaining proof gate.

## Phase 1: Audit Current Implementation

- [x] T001 Confirm the grouped CLI and legacy alias contract in `main.py`, `workflow/cli.py`, and `docking/cli.py`
- [x] T002 Confirm the current project layout and engine runner locations in `docking/project_layout.py` and `docking/runners/`

## Phase 2: Verification Completed

- [x] T003 Run `python -m py_compile main.py docking/*.py docking/preparation/*.py docking/runners/*.py workflow/*.py post_docking_analysis/multi_engine_pipeline.py`
- [x] T004 Run `python main.py dock --help`
- [x] T005 Run a synthetic staged-project smoke flow through `workflow init`, `prep pairlist`, `prep project`, and `dock run --dry-run`
- [x] T006 Update the spec and plan wording to match the grouped `dock run` surface and staged compatibility layout

## Phase 3: Remaining Proof

- [x] T007 Revalidate GNINA log-backed normalization into `normalized_scores.csv`
- [x] T008 Revalidate Vina/Smina pose-backed normalization into `normalized_scores.csv`
- [x] T009 Record the final normalization evidence and parser-fix note in `specs/007-multi-engine-docking-orchestration/plan.md`
