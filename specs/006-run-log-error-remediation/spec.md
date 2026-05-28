# Feature Specification: Simplified Pipeline Run-Log Error Remediation

## Context

A full run recorded in
`<project-root>/5-Post_Docking_Analysis/simplified_pipeline.log`
completes successfully but still emits high error/warning volume and misleading
stage success summaries.

## Run Appraisal (Observed in This Log)

The run has 23,490 log lines with 672 warnings and 72 errors, while still ending
with `Simplified pipeline completed successfully`.

The dominant hard error is repeated 72 times during PoseView processing:
`invalid literal for int() with base 10: 'accepted'`.
This is consistent with API payloads that return textual status values.

LigPlot stage shows high retry noise:
178 `HBADD failed` warnings and 166 `LIGPLOT returned non-zero` warnings,
yet the stage summary reports `LigPlot successful runs: 40` and
`LigPlot exported diagrams: 39`.
This indicates current success accounting can report success despite repeated
non-zero tool exits.

ProLIF stage reports 9 warnings for frame LigNetwork generation:
`None of ['ligand', 'protein', 'interaction', 'atoms'] are in the columns`.
The stage still generates output files, so this is a robustness/verbosity issue
rather than a pipeline-stopper.

The run also logs `No pairlist.csv - falling back to basic analysis`, which
downgrades comparative/hierarchical analysis when pairlist discovery is expected
from common GNINA layouts.

## Root Causes (Code-Level)

1. PoseView status parsing assumes `status_code` is numeric and directly casts
   strings via `int(...)` in multiple locations in
   `post_docking_analysis/poseview_integration.py`.
2. LigPlot success semantics in `post_docking_analysis/ligplot_integration.py`
   accept stage success if any `ligplot.*` file exists, even after repeated
   non-zero exits and fallback attempts.
3. ProLIF frame LigNetwork call in
   `post_docking_analysis/prolif_interaction_maps.py` executes without checking
   whether frame-level interaction columns exist.
4. Pairlist availability in explicit-folder mode is not resilient when
   `--pairlist` is omitted, causing avoidable fallback to basic analysis.

## Goals

1. Eliminate PoseView crash-path errors caused by textual API statuses.
2. Make LigPlot success/failure accounting truthful and auditable.
3. Reduce avoidable warning storms from repeated retries with identical failures.
4. Make ProLIF frame-network generation graceful when interaction dataframe
   schema is incomplete.
5. Improve pairlist resolution so hierarchical analysis is used whenever a valid
   pairlist is discoverable from standard project ancestry.

## Out of Scope

- Redesigning external tool algorithms (PoseView API, HBPLUS/HBADD internals,
  ProLIF internals).
- Changing docking scoring logic or affinity thresholds.
- Replacing LigPlot/PoseView/ProLIF with different interaction engines.

## Functional Requirements

- FR-001: PoseView client must normalize status values from both numeric and
  textual payload forms (`accepted`, `queued`, `processing`, `completed`, etc.)
  without raising conversion exceptions.
- FR-002: PoseView stage must report per-entry terminal outcomes with explicit
  failure reasons and total counts that match emitted outputs.
- FR-003: LigPlot runner must mark a complex as successful only when required
  output artifacts pass validation after a concrete attempt.
- FR-004: LigPlot stage summary must include processed, successful, failed, and
  exported counts with deterministic definitions.
- FR-005: LigPlot retry logging must be condensed to avoid repeated identical
  warnings for the same complex and failure signature.
- FR-006: ProLIF frame LigNetwork generation must be skipped gracefully when
  frame-level interaction schema is missing; this should not emit hard warnings
  per complex.
- FR-007: If `pairlist_file` is not explicitly provided, pipeline startup must
  attempt deterministic pairlist auto-discovery using SDF/log/receptor parent
  ancestry before falling back to basic analysis.

## Acceptance Criteria

1. Re-running the same dataset no longer produces PoseView `invalid literal for
   int() with base 10: 'accepted'` errors in `simplified_pipeline.log`.
2. PoseView summary counts are internally consistent with generated files and
   include non-empty error_reason distribution for failures.
3. LigPlot stage reports failures when valid outputs are not produced, and
   `successful` does not exceed count of validated result directories.
4. Warning volume from repeated LigPlot retries is substantially reduced by
   de-duplication/aggregation per complex.
5. ProLIF frame LigNetwork missing-column cases are handled via graceful skip
   and do not interrupt map generation.
6. When a pairlist exists in standard ancestry and is not explicitly passed,
   pipeline uses it and does not emit `No pairlist.csv` fallback warning.
7. `python -m py_compile` passes for all touched modules.

## Risks and Mitigations

API payload drift on proteins.plus can still occur; mitigation is to centralize
status normalization and treat unknown statuses conservatively (pending vs
failed) with explicit logs.

LigPlot output validity checks may be too strict or too loose; mitigation is to
define minimal required artifacts (`ligplot.drw` and non-empty `ligplot.ps`) and
record per-attempt diagnostics for tuning.
