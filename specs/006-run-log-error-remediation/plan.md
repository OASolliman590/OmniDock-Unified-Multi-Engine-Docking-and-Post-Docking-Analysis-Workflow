# Implementation Plan: Simplified Pipeline Run-Log Error Remediation Closeout

## Current Status

Most of this remediation work is now implemented in the current codebase:

- PoseView status normalization handles mixed textual and numeric responses
- PoseView analyzer success now trusts normalized terminal state rather than a
  raw `status_code == 200` gate
- LigPlot success accounting is tied to validated outputs and stage-level
  `processed/successful/failed/exported/failure_reasons` summaries
- ProLIF frame LigNetwork missing-schema cases are skipped gracefully
- pairlist auto-discovery is implemented and now searches deeper ancestry from
  the resolved project/input roots

## Verification Completed

1. `python -m py_compile post_docking_analysis/poseview_integration.py post_docking_analysis/ligplot_integration.py post_docking_analysis/prolif_interaction_maps.py post_docking_analysis/simplified_pipeline.py post_docking_analysis/simplified_cli.py post_docking_analysis/simplified_input_handler.py post_docking_analysis/gnina_hpc_adapter.py`
2. Targeted PoseView regression proof using a fake completed response with
   `status_code=201` and `status_state=completed`
3. Targeted deep-layout pairlist regression proof covering both
   `auto_detect_pairlist_file()` and `_find_pairlist_file()`
4. Code-path audit confirmed:
   - condensed LigPlot failure accounting
   - ProLIF frame guardrails
   - pairlist auto-detection logging in the simplified pipeline

## Remaining Gap

This spec remains open because the original rerun-backed sign-off has not been
re-recorded against a real simplified-pipeline dataset after the closeout fixes.

## Remaining Validation Path

1. Re-run the target simplified pipeline command on a representative GNINA run
2. Inspect the resulting log for:
   - PoseView terminal outcomes and failure reasons
   - LigPlot summary consistency
   - reduced warning duplication
   - auto-resolved pairlist behavior when discoverable
3. Record final external-tool caveats after that rerun
