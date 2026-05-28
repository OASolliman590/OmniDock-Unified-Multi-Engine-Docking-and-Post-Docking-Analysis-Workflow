# Implementation Plan: Multi-Engine Docking Orchestration

## Current Checkpoint

- The core implementation is already landed in `docking/`, `main.py`, and
  `workflow/cli.py`.
- The current public contract is the grouped CLI surface
  `python main.py dock run ...`, with legacy `python main.py dock ...`
  invocations normalized to `dock run`.
- The runner stack already provides:
  - project layout helpers and manifest models
  - engine registry plus GNINA, Vina, and Smina runners
  - GNINA log normalization through `generate_all_scores_csv`
  - Vina/Smina PDBQT normalization through `post_docking_analysis/docking_parser.py`

## Verified In This Audit

1. `python -m py_compile main.py docking/*.py docking/preparation/*.py docking/runners/*.py workflow/*.py post_docking_analysis/multi_engine_pipeline.py`
2. `python main.py dock --help`
3. Synthetic staged-project smoke flow:
   - `python main.py workflow init --project-dir /tmp/pdbw_spec_resume_fixture/project`
   - `python main.py prep pairlist --project-dir /tmp/pdbw_spec_resume_fixture/project --mode cocrystal_only`
   - `python main.py prep project --project-dir /tmp/pdbw_spec_resume_fixture/project`
   - `python main.py dock run --project-dir /tmp/pdbw_spec_resume_fixture/project --dry-run`
4. Synthetic normalization rerun under `/tmp/pdbw_norm_fixture`:
   - GNINA log-backed normalization produced `4-Docking/results/normalized_scores.csv` with 2 rows.
   - Vina pose-backed normalization produced `4-Docking/vina_out/scores/normalized_scores.csv` with 2 rows.
   - Smina pose-backed normalization produced `4-Docking/smina_out/scores/normalized_scores.csv` with 2 rows.
5. Validation uncovered and fixed a GNINA parsing bug in `post_docking_analysis/generate_scores_csv.py` where score extraction stopped after the first pose row.

## Remaining Validation Gap

- None in this closeout pass.
