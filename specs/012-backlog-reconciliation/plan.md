# Implementation Plan: Backlog Reconciliation After Workflow Foundation

## Current Checkpoint

- The active branch is now `012-backlog-reconciliation`.
- `git branch --show-current` confirms the numbered branch alignment.
- `python -m py_compile` passes for the audited post-docking modules:
  - `post_docking_analysis/protein_naming.py`
  - `post_docking_analysis/pandamap_runner.py`
  - `post_docking_analysis/pandamap_integration.py`
  - `post_docking_analysis/publication_pandamap.py`
  - `post_docking_analysis/hierarchical_analyzer.py`
  - `post_docking_analysis/simplified_input_handler.py`
  - `post_docking_analysis/simplified_cli.py`
  - `post_docking_analysis/simplified_pipeline.py`
  - `post_docking_analysis/poseview_integration.py`
  - `post_docking_analysis/ligplot_integration.py`
  - `post_docking_analysis/prolif_interaction_maps.py`
  - `post_docking_analysis/pipeline.py`
- `python -m post_docking_analysis.simplified_cli --help` confirms the current
  simplified CLI surface, including `--prompt-protein-names`,
  `--enable-poseview`, and `--interactive`.
- Specs `007`, `009`, `010`, and `011` represent the validated workflow
  foundation for the next planning cycle.
- This closeout pass also fixed two defects discovered during validation:
  - `post_docking_analysis/generate_scores_csv.py`: GNINA score parsing no
    longer stops after the first pose row.
  - `workflow/execution.py`: all-failed `pdb collect` runs now preserve a
    visible `Summary` sheet instead of crashing during workbook save.

## Verified Foundation To Carry Forward

1. `specs/007-multi-engine-docking-orchestration`
   - `python main.py dock --help`
   - synthetic `dock run --dry-run`
   - synthetic GNINA/Vina/Smina normalization into `normalized_scores.csv`
2. `specs/009-docking-preparation-automation`
   - `python main.py workflow init --help`
   - `python main.py prep pairlist --help`
   - `python main.py prep project --help`
   - live `python main.py pdb collect --project-dir /tmp/pdbw_collect_fixture -p 1STP --selection-mode heuristic`
3. `specs/010-unified-workflow-orchestrator`
   - grouped CLI and workflow shell are already checked off in the existing task list
4. `specs/011-docking-first-project-workflow`
   - staged docking project flow is already checked off in the existing task list

## Legacy Spec Matrix

| Spec | Status | Basis |
| --- | --- | --- |
| `001-unify-and-optimize` | `superseded` | Its broad visualization/output scope was replaced by later narrower specs. Protein naming, RMSD topology, and output consolidation are already present in the current code. |
| `002-pipeline-unification` | `completed` | Shared PandaMap execution lives in `post_docking_analysis/pandamap_runner.py`, the analyzer callers use it, and `pipeline.py` no longer carries the dead duplicate execution path. |
| `003-simplified-pipeline-dedup` | `superseded` | Best-pose and complex-file helper consolidation exists, but the original tracking is overtaken by later output-contract work in `004` and resilience work in `006`. |
| `004-visualization-rmsd-overhaul` | `open` | Labels, affinity visuals, RMSD scope, and interaction routing are implemented, but spec closeout remains stale and docs still describe older behavior in places. |
| `005-interactive-cli-prompts` | `completed` | Guided prompts, non-TTY failure behavior, and per-target protein naming prompts are implemented; the grouped workflow shell now covers the primary interactive entry path. |
| `006-run-log-error-remediation` | `open` | PoseView normalization, LigPlot accounting, ProLIF frame guards, and pairlist auto-discovery are implemented, but rerun-backed closeout proof is missing and one PoseView completion edge remains. |
| `007-multi-engine-docking-orchestration` | `completed` | Revalidated with grouped CLI help, staged dry-run flow, and GNINA/Vina/Smina normalization outputs. |
| `008-engine-aware-post-docking-analysis` | `completed` | Task list is fully checked and the engine-aware continuation/reporting path is already captured as done. |
| `009-docking-preparation-automation` | `completed` | Revalidated with grouped prep help, staged materialization smoke flow, and live `pdb collect`. |
| `010-unified-workflow-orchestrator` | `completed` | Grouped CLI, workflow shell, wrappers, and validation tasks are already checked off. |
| `011-docking-first-project-workflow` | `completed` | Docking-project-first layout, `pdb collect`, `prep pairlist`, `prep project`, and interactive shell updates are already checked off. |

## Validation Path

1. `git branch --show-current`
2. `python -m py_compile post_docking_analysis/generate_scores_csv.py workflow/execution.py`
3. Reuse the validated commands already recorded in:
   - `specs/007-multi-engine-docking-orchestration/plan.md`
   - `specs/009-docking-preparation-automation/plan.md`

## Reconciliation Findings

1. The early backlog is not uniformly stale:
   - `002` and `005` are functionally complete
   - `001` and `003` are superseded tracking shells
   - `004` and `006` still contain real closeout work
2. The validated workflow foundation now starts at `007` and continues through
   `011` without an unresolved branch/spec mismatch.
3. The next feature should not be a net-new capability. The real next step is
   to close the remaining post-docking gaps left in `004` and `006`.

## Selected Next Feature

- Feature branch to create: `013-post-docking-closeout`
- Why this branch exists:
  - `004` still has documentation drift around visualization outputs and
    simplified CLI flags
  - `006` still lacks rerun-backed closeout evidence and has a remaining
    PoseView completion-state edge
- Proposed scope:
  - align `POST_DOCKING_ANALYSIS_GUIDE.md`,
    `post_docking_analysis/README.md`, and
    `post_docking_analysis/SIMPLIFIED_IMPLEMENTATION_STATUS.md` with the
    current simplified pipeline behavior
  - update `post_docking_analysis/poseview_integration.py` so terminal success
    relies on normalized completion state rather than a raw `status_code == 200`
    gate
  - either widen `auto_detect_pairlist_file()` for deeper common layouts or
    explicitly document the supported ancestry contract
  - record the actual closeout evidence in the legacy `004` and `006` plan/task
    files

## Initial Validation Path For `013-post-docking-closeout`

1. `python -m py_compile post_docking_analysis/simplified_cli.py post_docking_analysis/simplified_pipeline.py post_docking_analysis/poseview_integration.py post_docking_analysis/simplified_input_handler.py`
2. `python -m post_docking_analysis.simplified_cli --help`
3. Run one non-interactive missing-argument check and confirm a clear failure
   message instead of prompts
4. Run one canonical simplified-pipeline smoke flow without an explicit
   `--pairlist` and confirm either:
   - automatic pairlist resolution is logged, or
   - the ancestry limit is documented as the intended contract
5. If network access is available, run one bounded PoseView smoke or otherwise
   record a fixture-backed completion-state proof in the plan

## Promotion Decision

- `012-backlog-reconciliation` stops at the audit and planning boundary.
- The selected next work should be promoted into a new numbered branch rather
  than implemented on `012`.
