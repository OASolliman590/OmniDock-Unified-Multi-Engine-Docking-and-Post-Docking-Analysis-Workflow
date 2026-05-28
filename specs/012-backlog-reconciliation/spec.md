# Feature Specification: Backlog Reconciliation After Workflow Foundation

**Feature Branch**: `012-backlog-reconciliation`  
**Created**: 2026-03-08  
**Status**: Completed  
**Input**: Current repository state after validating specs `007`, `009`, `010`, and `011`

## Context

The repository’s implementation state advanced past the active branch name.
Core workflow, docking-preparation, docking-execution, and engine-aware
analysis foundations are now in place, but the early specs (`001` to `006`)
still need an authoritative status pass before the next feature branch is
chosen.

## Goals

1. Reconcile legacy specs `001` to `006` against the current codebase and docs.
2. Preserve the validated status of specs `007`, `009`, `010`, and `011`.
3. Define exactly one next implementation-ready feature after reconciliation.
4. Restore numbered branch/spec alignment so Spec Kit tooling resolves the
   active feature correctly.

## Functional Requirements

- FR-001: The active Git branch MUST be `012-backlog-reconciliation`.
- FR-002: `specs/012-backlog-reconciliation/` MUST contain `spec.md`,
  `plan.md`, and `tasks.md`.
- FR-003: The reconciliation output MUST classify each legacy spec from
  `001` to `006` as `completed`, `superseded`, or `open`, with a short
  rationale grounded in current code/docs.
- FR-004: The reconciliation output MUST record the currently validated
  foundation:
  - `007-multi-engine-docking-orchestration`
  - `009-docking-preparation-automation`
  - `010-unified-workflow-orchestrator`
  - `011-docking-first-project-workflow`
- FR-005: The reconciliation output MUST identify one next
  implementation-ready feature to continue after this audit.
- FR-006: The reconciliation output MUST define the initial validation path for
  that next feature.

## Acceptance Criteria

1. `git branch --show-current` returns `012-backlog-reconciliation`.
2. The `specs/012-backlog-reconciliation/` directory contains `spec.md`,
   `plan.md`, and `tasks.md`.
3. The plan documents the verified `007`/`009` evidence and the legacy-spec
   audit scope.
4. The task list closes reconciliation work before any new implementation
   branch is proposed.

## Reconciliation Outcome

| Spec | Status | Rationale |
| --- | --- | --- |
| `001-unify-and-optimize` | `superseded` | This was an umbrella visualization overhaul. Its implemented parts were carried forward into later narrower specs, especially `004`, while its own task list never became the authoritative closeout record. |
| `002-pipeline-unification` | `completed` | Shared PandaMap execution is centralized, the duplicate pipeline flow is removed, and the affected modules compile cleanly. |
| `003-simplified-pipeline-dedup` | `superseded` | The core helper consolidation landed, but the original spec was overtaken by later output-contract and resilience work in `004` and `006`. |
| `004-visualization-rmsd-overhaul` | `open` | Core implementation is present, but the closeout record is stale and user-facing docs still drift from the current behavior. |
| `005-interactive-cli-prompts` | `completed` | The simplified CLI prompts missing inputs interactively, fails clearly in non-TTY mode, and supports per-target protein naming prompts. |
| `006-run-log-error-remediation` | `open` | Most remediation is implemented, but rerun-backed closeout evidence is missing and PoseView completion handling still has a raw status-code edge. |

## Foundation Confirmed

- `007-multi-engine-docking-orchestration`: completed and revalidated
- `008-engine-aware-post-docking-analysis`: completed
- `009-docking-preparation-automation`: completed and revalidated
- `010-unified-workflow-orchestrator`: completed
- `011-docking-first-project-workflow`: completed

## Next Feature Decision

- Selected next feature: `013-post-docking-closeout`
- Reason: the only materially open backlog after reconciliation is the residual
  post-docking closeout work split across `004` and `006`
- Initial scope:
  - align simplified CLI and visualization docs with current behavior
  - close the remaining PoseView completion-state edge
  - harden or explicitly bound pairlist auto-discovery behavior
  - record rerun-backed closeout evidence for the post-docking path
- Delivery rule: keep `012` as the audit branch and promote implementation into
  the next numbered branch instead of continuing feature work on `012`
