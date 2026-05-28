# Implementation Plan: Interactive Simplified CLI Prompts Closeout

## Current Status

This feature is functionally complete in the current codebase:

- `post_docking_analysis/simplified_cli.py` supports guided prompting for
  missing required inputs during TTY sessions
- non-interactive runs fail with clear, deterministic errors instead of trying
  to prompt
- per-target protein naming prompts are implemented in the simplified pipeline
- `workflow interactive` is now the main interactive entry point for the wider
  application, while the simplified CLI remains available directly
- help surfaces are now lightweight because the simplified pipeline import is
  lazy on `--help`

## Verification Completed

1. `python -m py_compile main.py workflow/cli.py post_docking_analysis/simplified_cli.py post_docking_analysis/simplified_pipeline.py`
2. `python -m post_docking_analysis.simplified_cli --help`
3. `python main.py`
4. `python main.py workflow interactive`
5. Code-path audit of per-target protein naming prompts in
   `post_docking_analysis/simplified_pipeline.py`

## Remaining Gap

No implementation gap remains for this feature. Any future work here is normal
UX/documentation maintenance rather than unfinished spec scope.
