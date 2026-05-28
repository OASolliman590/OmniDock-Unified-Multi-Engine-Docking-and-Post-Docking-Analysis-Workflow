# Implementation Plan: Docking Preparation Automation

## Current Checkpoint

- The staged docking-preparation flow is already landed across
  `docking/preparation/`, `docking/project_layout.py`, `workflow/execution.py`,
  and `workflow/interactive.py`.
- The current public contract is the grouped CLI:
  - `python main.py workflow init ...`
  - `python main.py prep pairlist ...`
  - `python main.py prep project ...`
- Legacy `python main.py prepare-docking ...` continues to normalize to
  `prep project`.

## Verified In This Audit

1. `python -m py_compile main.py docking/*.py docking/preparation/*.py docking/runners/*.py workflow/*.py post_docking_analysis/multi_engine_pipeline.py`
2. Help checks:
   - `python main.py workflow init --help`
   - `python main.py prep pairlist --help`
   - `python main.py prep project --help`
3. Synthetic staged-project smoke flow:
   - `python main.py workflow init --project-dir /tmp/pdbw_spec_resume_fixture/project`
   - `python main.py prep pairlist --project-dir /tmp/pdbw_spec_resume_fixture/project --mode cocrystal_only`
   - `python main.py prep project --project-dir /tmp/pdbw_spec_resume_fixture/project`
4. Live `pdb collect` rerun on `/tmp/pdbw_collect_fixture`:
   - `python main.py workflow init --project-dir /tmp/pdbw_collect_fixture`
   - `python main.py pdb collect --project-dir /tmp/pdbw_collect_fixture -p 1STP --selection-mode heuristic`
   - Result: downloaded `1STP`, extracted `BTN_A_300`, created
     `1-Raw_Ligand/1STP_ligand_BTN_A_300.pdb`,
     `1-Raw_Ligand/1STP_ligand_BTN_A_300.sdf`, and populated
     `2-Raw_Protien/multi_pdb_analysis.xlsx` with the `Summary` sheet.
5. Validation uncovered and fixed a workbook-save failure path in
   `workflow/execution.py` so all-failed collect runs now preserve a visible
   `Summary` sheet instead of crashing with `At least one sheet must be visible`.

## Remaining Validation Gap

- None in this closeout pass.
