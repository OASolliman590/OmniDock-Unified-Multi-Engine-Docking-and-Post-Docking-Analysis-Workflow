# Implementation Plan: Post-Docking Closeout

## Current Checkpoint

- Active branch: `013-post-docking-closeout`
- Reconciliation source: `012-backlog-reconciliation`
- First closeout fixes already landed in this branch:
  - `main.py`: no-argument entry now opens interactive workflow only in TTY
    sessions and falls back to top-level help otherwise.
  - `workflow/cli.py`: heavy execution imports are now lazy, non-TTY guards were
    added for interactive-only workflow/PDB routes, and analysis help surfaces
    no longer pay the execution import cost.
  - `post_docking_analysis/simplified_cli.py`: the simplified pipeline import is
    now lazy, so help text no longer imports the analysis stack.
  - `post_docking_analysis/poseview_integration.py`: analyzer success now trusts
    normalized completion state instead of raw `status_code == 200`.
  - `post_docking_analysis/simplified_input_handler.py` and
    `post_docking_analysis/gnina_hpc_adapter.py`: pairlist auto-detection now
    searches up to four parent levels above the starting directory/project root.

## Closeout Result

This branch is complete:

- user-facing post-docking docs now match the current CLI, visualization, and
  pairlist behavior
- the stale `004`, `005`, and `006` plan/task files are rewritten as accurate
  closeout records
- the only remaining backlog is the real smoke/rerun proof still called out in
  `004` and `006`

## Validation Completed So Far

1. `python -m py_compile main.py workflow/cli.py post_docking_analysis/simplified_cli.py post_docking_analysis/poseview_integration.py`
2. `python main.py --help`
3. `python main.py`
4. `python main.py workflow interactive`
5. `python main.py pdb legacy-interactive`
6. `python main.py pdb collect --project-dir /tmp/pdbw_dummy -p 1STP --selection-mode interactive`
7. `python main.py analyze --help`
8. `python -m post_docking_analysis.simplified_cli --help`
9. Targeted PoseView regression proof using a fake completed response with
   `status_code=201` and `status_state=completed`
10. Targeted deep-layout pairlist regression proof covering both
    `auto_detect_pairlist_file()` and `_find_pairlist_file()`

## Final Validation Path

1. `python -m py_compile` for the touched CLI/runtime modules
2. Re-run help and non-TTY guard checks after the entry-point fixes
3. Run targeted PoseView and pairlist regression proofs
4. Verify `specs/013-post-docking-closeout/` resolves as the active feature dir

## Expected Deliverable

1. A usable CLI entry surface for both scripted and interactive contexts
2. A closed PoseView completion-state edge
3. Updated docs and legacy spec records that reflect the actual remaining work
