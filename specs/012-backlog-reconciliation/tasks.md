# Tasks: Backlog Reconciliation After Workflow Foundation

**Input**: Design documents from `/specs/012-backlog-reconciliation/`
**Prerequisites**: `plan.md`, `spec.md`

**Tests**: Use `git branch --show-current`, `python -m py_compile`, and the
already-recorded verification commands from specs `007` and `009`.

## Phase 1: Setup

- [x] T001 Create the numbered branch `012-backlog-reconciliation`
- [x] T002 Scaffold `spec.md`, `plan.md`, and `tasks.md` under `specs/012-backlog-reconciliation/`
- [x] T003 Record the validated `007`/`009` closeout state in `specs/012-backlog-reconciliation/plan.md`

## Phase 2: Legacy Spec Audit

- [x] T004 Audit `specs/001-unify-and-optimize/` against the current post-docking code and docs
- [x] T005 Audit `specs/002-pipeline-unification/` against the current PandaMap/pipeline code paths
- [x] T006 Audit `specs/003-simplified-pipeline-dedup/` against the current simplified pipeline helpers
- [x] T007 Audit `specs/004-visualization-rmsd-overhaul/` against the current visualization/RMSD outputs and docs
- [x] T008 Audit `specs/005-interactive-cli-prompts/` against the current grouped CLI and workflow shell
- [x] T009 Audit `specs/006-run-log-error-remediation/` against the current optional-tool and run-log behavior

## Phase 3: Reconciliation Output

- [x] T010 Build a status matrix for specs `001` through `011` in `specs/012-backlog-reconciliation/plan.md`
- [x] T011 Classify specs `001` through `006` as `completed`, `superseded`, or `open` in `specs/012-backlog-reconciliation/spec.md`
- [x] T012 Choose one next implementation-ready feature after reconciliation and record it in `specs/012-backlog-reconciliation/spec.md`

## Phase 4: Handoff

- [x] T013 Define the initial validation path for the selected next feature in `specs/012-backlog-reconciliation/plan.md`
- [x] T014 Decide whether the selected next feature should continue on `012-backlog-reconciliation` or be promoted into the next numbered feature branch

## Outcome

- `012-backlog-reconciliation` is now the authoritative audit/decision record.
- The next implementation work should move to `013-post-docking-closeout`.
