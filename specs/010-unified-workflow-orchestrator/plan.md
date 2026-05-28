# Implementation Plan: Unified Workflow Orchestrator

1. Add a new `workflow/` package for:
   - grouped CLI parsing
   - interactive workflow prompts
   - workflow-state persistence
   - execution adapters over existing pipeline modules
2. Replace `main.py` with a grouped command router and legacy-alias
   normalization.
3. Wrap existing PDB, docking-prep, docking, and post-docking surfaces instead
   of reimplementing their scientific logic.
4. Add stage-level analysis adapters that expose GNINA simplified-pipeline
   stages and canonical-engine bridge stages through stable workflow targets.
5. Persist command outcomes to `.workflow/state.json` and expose `workflow
   status`, `workflow resume`, and `workflow jump`.
6. Add `workflow init` so project roots can be bootstrapped before interactive
   or CLI execution. GNINA-enabled projects must also expose root-level
   GNINA HPC-compatible paths while preserving the canonical manifest-driven
   multi-engine layout.
7. Update user-facing docs and add Spec Kit coverage for the orchestration
   feature.
8. Validate with:
   - `python -m py_compile main.py workflow/*.py`
   - grouped CLI help checks
   - legacy alias help checks
