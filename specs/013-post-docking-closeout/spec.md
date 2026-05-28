# Feature Specification: Post-Docking Closeout

**Feature Branch**: `013-post-docking-closeout`  
**Created**: 2026-03-08  
**Status**: Completed  
**Input**: Reconciliation outcome from `012-backlog-reconciliation` plus observed pipeline usability failures in the current CLI surfaces

## Context

`012-backlog-reconciliation` identified the remaining real backlog as closeout work
spread across `004-visualization-rmsd-overhaul` and
`006-run-log-error-remediation`. While continuing that closeout, the current
entry points also showed concrete usability failures:

- `python main.py` with no arguments failed in non-interactive shells instead of
  falling back to help.
- interactive-only routes were not clearly guarded for non-TTY sessions.
- `--help` on key CLI surfaces imported the heavy execution stack and emitted
  avoidable stderr noise.
- PoseView completion handling still trusted a raw `status_code == 200` gate
  even though normalized completion state was already available.

## Goals

1. Make the main entry points predictable in both TTY and non-TTY contexts.
2. Keep CLI help surfaces lightweight and readable.
3. Finish the remaining PoseView completion-state fix from spec `006`.
4. Close the remaining documentation and closeout drift from specs `004` to `006`.
5. Either harden or explicitly bound the pairlist auto-discovery contract.

## Functional Requirements

- FR-001: `python main.py` with no arguments MUST open the interactive workflow
  only in a real TTY session.
- FR-002: In non-TTY contexts, interactive-only commands MUST fail with a clear,
  actionable message instead of surfacing prompt-toolkit or stdin errors.
- FR-003: `python main.py --help`, `python main.py analyze --help`, and
  `python -m post_docking_analysis.simplified_cli --help` MUST avoid importing
  the heavy execution path when help text is sufficient.
- FR-004: PoseView terminal success handling MUST rely on normalized completion
  state rather than a raw `status_code == 200` check.
- FR-005: The remaining doc drift from specs `004`, `005`, and `006` MUST be
  recorded and closed explicitly.
- FR-006: The pairlist auto-discovery depth MUST support the expected staged
  project ancestry and be documented explicitly in user-facing docs.

## Acceptance Criteria

1. `python main.py` with no args in a non-TTY session prints top-level help and
   exits successfully.
2. `python main.py workflow interactive`, `python main.py workflow resume`,
   `python main.py pdb legacy-interactive`, and
   `python main.py pdb collect --selection-mode interactive ...` fail cleanly in
   non-TTY mode.
3. The key help surfaces above print without the previous import-time execution
   noise.
4. A bounded PoseView regression proof demonstrates that a normalized-completed
   result with a non-200 status code is treated as success.
5. `specs/004-visualization-rmsd-overhaul/`,
   `specs/005-interactive-cli-prompts/`, and
   `specs/006-run-log-error-remediation/` are left with accurate remaining-work
   tracking.
