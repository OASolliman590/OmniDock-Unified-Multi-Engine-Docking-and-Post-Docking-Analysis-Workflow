# Spec 038 R3 (T004): post-docking inventory and STOP report

**Status (2026-10-19): STOPPED for Lead decision.** The golden baseline is captured and reproducible. The two
post-docking paths disagree on outputs that change scientific or status content, so R3 code changes are not made.
No file in `post_docking_analysis/` was changed. Section 5 lists the decisions needed before R3 can continue.

Scratch artefacts (not in the repo): golden captures and harness under the session scratchpad,
`spec038_golden/baseline` (frozen), `spec038_golden/baseline_repeat` and `spec038_golden/recheck2` (determinism checks),
`spec038_harness/capture_golden.py` and `compare_golden.py`.

## 1. Golden baseline

- **Cases (16 runs).** Each run uses a fresh project copy and is run through two entry paths.
  - The DAG path is `UnifiedPostDockingPipeline.run()`.
  - The legacy path is `MultiEngineAnalysisPipeline.run()`.
  - `smoke_seeded`: the smoke test's own seeded 3-engine project (gnina, vina, smina; one receptor, one ligand). Runs: comparative full, single vina, single gnina.
  - `synthetic`: 2 receptors x 4 ligands x 3 engines x 3 poses. Ligand L4 has positive affinities. Runs: comparative full, comparative full with the filter on, single smina, single gnina, favorite_engine_continue (vina).
- **Captured.** Every file under the project and the session. Text files (CSV, JSON, TXT, MD, HTML) are stored normalised. Masked: the work root, timestamps, run ids, durations, and cache keys, source signatures and fingerprints. Structural outputs are recorded but not compared: `.artifact_history/` (random ids) and `run_tracking/outputs_index.*` (sizes and mtimes). PNG files are compared by decoded pixels.
- **Determinism.** `PYTHONHASHSEED=0`. Two captures of the unchanged code: 2573 compared files identical. The only differences are 36 PNG pixel differences in two seaborn figures, `cross_engine_engine_distribution.png` and `affinity_by_protein.png`. The cause is unseeded jitter in `sns.stripplot(jitter=0.18)` (`multi_engine_pipeline_impl.py`, around line 4158). Their CSV inputs are compared and identical.
- **Recheck on the current tree** (after the other agents' edits and the brief `git stash`): identical to the baseline on all 2573 compared files. `post_docking_analysis/` has no changes from other agents.
- **Baseline runs on the untouched code.**
  - Post-docking subset (`test_postdocking_correctness.py`, `test_spec031_*`, `035_*`, `036_*`, `037_*`): 319 passed, 1 skipped, 39 xfailed, 5 xpassed.
  - Contract smoke (`--skip-all-engines --skip-prep-matrix`): exit 0.
  - Full smoke (no flags, includes the all-engines dry run): exit 0.
  - Full pytest: 724 passed, 3 skipped, 49 xfailed, 5 xpassed, exit 0.

## 2. Entry points and callers

| Module | Lines | Role | Callers | Proposed R3 action |
|---|---|---|---|---|
| `multi_engine_pipeline_impl.py` | 6597 | `MultiEngineAnalysisPipeline` (`run()` 222 lines), every DAG node, the comparative writer, the single-engine writers | `unified_pipeline.py`, `workflow/execution.py` (md_inputs), tests | Keep as the implementation. `run()` becomes a DAG wrapper once D-A is decided. |
| `multi_engine_pipeline.py` | 9 | re-export shim | `unified_pipeline.py` | Keep |
| `unified_pipeline.py` | 279 | `UnifiedPostDockingPipeline`; its `run()` already executes the DAG | `cli.py`, `workflow/execution.py` (lines ~1626, ~1978, ~2862), tests | Keep as the single entry class |
| `simplified_pipeline_impl.py` | 5183 | `SimplifiedPostDockingPipeline`. Its own RMSD, PLIP, PoseView, PyMOL, PandaMap and visualisation stages. Also used as a bridge by the multi-engine module (lines ~963, ~6353, ~6563) | `simplified_cli.py`, `workflow/execution.py` (~2672, ~2899) | Compatibility-layer target, blocked by the bridge (section 4, E-2) |
| `simplified_pipeline.py` | 9 | re-export shim | `workflow/execution.py` | Keep |
| `simplified_cli.py` | 536 | CLI `python -m post_docking_analysis.simplified_cli`. Delegates to the unified path for manifest projects. | Documented in `README.md`, `POST_DOCKING_ANALYSIS_GUIDE.md`, `VISUALIZATION_GUIDE.md`; tests | Keep as a thin shim over the unified path |
| `cli.py` | 351 | console script `post-docking-analysis` (`setup.py`) and `python -m post_docking_analysis`. Uses `UnifiedPostDockingPipeline`. | `setup.py`, docs, `skills/post-docking-contract-review` | Keep |
| `pipeline.py` | 1484 | legacy standalone `PostDockingAnalysisPipeline` (own `sys.path` hack to a `docking_analysis` directory) | `examples/run_pipeline.py`, `example_usage.py`, `example_pairlist.py` only (no workflow, webui or test callers) | Retirement candidate. Examples must be updated first. |

Shared modules named in Spec 038 R3: `pose_selection.py` (109 lines), `consensus.py`, `replicates.py` (209), `docking_parser.py`, `docked_microspecies.py`. `multi_engine_pipeline_impl.py` imports all of them. `simplified_pipeline_impl.py` imports none of them.

## 3. Duplicate functions

Found by AST comparison of function bodies (names and docstrings ignored) and by name across the two implementations.

| Function | simplified_pipeline_impl.py | multi_engine_pipeline_impl.py | Identical? | Action |
|---|---|---|---|---|
| `_parse_atom_serial` | 1347 | 4912 | yes (AST) | Keep one shared copy |
| `_normalize_rmsd_scopes` | 105 | 276 | same role, bodies differ | Consolidate |
| `_optional_feature_record` | 126 | 302 | differs | Consolidate |
| `_infer_element_from_pdbqt` | 1268 (60 lines) | 4807 (54 lines) | differs | Consolidate into one PDBQT helper |
| `_normalize_ligand_resname` | 1384 | 4982 | differs | Consolidate |
| `_derive_ligand_identity` | 1427 | 5016 | differs | Consolidate |
| `_build_complex_metadata` | 3157 | 1043 | differs | Consolidate |
| `_classify_output_path` | 393 | 1807 | differs | Consolidate |
| `_infer_optional_feature_statuses` | 536 (144 lines) | 2033 (26 lines) | differs | Consolidate |
| `_prune_empty_directories` | 5003 | 1849 | differs | Consolidate |
| `_scope_signature`, `_persist_scope_state` | 2469, 2510 | 5897, 5918 | differs | Consolidate |
| best-pose selection | `_best_pose_per_protein_ligand` (3206, idxmin on `vina_affinity`), `_best_poses_from_scores` (3047) | `_best_rows_by_group` (5222), `_best_by_tag` (5236), `_select_best_pose_rows` (5276), which wrap shared `pose_selection.select_best_pose_rows` | different rules | Do not collapse until the single-engine GNINA rule is confirmed (section 5, D-C) |
| score normalisers | none | `_build_normalized_frame` (2699) | n/a | Source of truth is `consensus.py` (`normalize_engine_scores`, `_safe_rank_score`, `_safe_minmax_scale`, `_safe_zscore`). `correlation_analyzer._rank_normalize` (41) is a second rank normaliser for a diagnostic. Check semantics before merging. |
| score parsing | `_generate_scores_csv` (1074) uses `generate_scores_csv` only | `_build_normalized_frame` uses `docking_parser` (Vina, AD4) and `generate_scores_csv` | n/a | Keep `docking_parser` as the parser for Vina and AD4 |
| summary writers | `_generate_reports` (3022) | `_write_comparative_reports` (2916; `summary.txt` at 3543), `_write_single_engine_reports` (4641; `summary.txt` at 4755), `report_generator.generate_summary_report` (415, shared) | differs | One writer (section 4, D-A1) |
| `main` (CLI) | `simplified_cli.py:261` (272 lines) | `cli.py:53` (296 lines); `pipeline.py:1471` (11 lines) | differs | One argument parser |
| `run` | `SimplifiedPostDockingPipeline.run` 681 (105 lines) | `MultiEngineAnalysisPipeline.run` 2188 (222 lines); `UnifiedPostDockingPipeline.run` 233 (44 lines, DAG) | differs | `run()` = DAG wrapper (D-A) |

## 4. DAG and legacy `run()` on the same inputs (golden evidence)

Comparison of session outputs, after the masking in section 1. Scope columns (`scope_source`, `excluded_engines` and related) are excluded where stated.

**Identical between the two paths:**
- `reports/consensus_ranked_hits.csv`, `reports/engine_summary.csv`, `reports/best_engine_per_complex.csv`, `reports/combined_engine_scores.csv` (8, 3, 8 and 72 rows): byte-identical after normalisation.
- Top-pose files: consensus rank, winner engine and consensus score are identical for every ligand.
- Single-engine vina and gnina: `best_poses.csv`, `scores.csv`, `protein_summary.csv`, `best_pose_per_tag_by_engine.csv` identical, apart from scope metadata. The scope difference is an artefact of the harness: the legacy call did not pass the scope resolution that the unified wrapper computes.

**Divergences (STOP):**

| ID | Output | Legacy `run()` | DAG (`UnifiedPostDockingPipeline.run()`) | Class |
|---|---|---|---|---|
| D-A1 | `session/reports/summary.txt` (comparative) | full multi-engine summary | the placeholder "Report generation completed through comparative writer, but summary artifact was absent." in all three DAG captures, including the repeat and recheck | report content. No rank changes. |
| D-A2 | `top_pose_ligand_performance/top_pose_per_ligand_global.csv` and `..._per_ligand_per_protein.csv` | `normalized_affinity_score` = `per_engine_rank` values (1.0, 0.7727, 0.5, 0.1818 on the synthetic fixture) and `normalization_method` | `normalized_affinity_score` NaN on every row. Also missing: `redocking_validation_status`, `redocking_validation_reason`, `redocking_reference`, `pose_reproducibility_*`, `docked_microspecies_smiles`, `protonation_policy`, `docked_state_source`, `replicates_pooled`, `normalization_method`, engine-scope columns | score value and status columns. Ranks are unchanged. |
| D-A3 | `session/favorite_engine/` (15 files: `best_poses.csv`, `scores.csv`, `summary.txt`, `downstream_bridge.txt`, `deep_analysis/`) | written (favorite_engine_continue stage) | not written. The DAG has no favorite stage. | output set |
| D-B | `excluded_problematic_ligands.csv` and the consensus tables when the filter is enabled | filter on: excludes `R1_site_1_L4` and `R2_site_1_L4` (`best_affinity>0.0`); 91 text files change, including ranks | filter ignored. Filter-on output is identical to filter-off output. | opt-in (default off), scientific when on |
| E-1 | legacy single-engine smina | `AttributeError: engine_detection_report` (`_write_single_engine_reports`, line ~4677; the attribute is set only by the `UnifiedPostDockingPipeline` subclass) | completes | pre-existing crash |
| E-2 | `SimplifiedPostDockingPipeline` bridge | used by the multi-engine module (lines ~963, ~6353, ~6563) | same | blocks removing the simplified implementation |

**DAG-only outputs and copies:** `visualizations/`, `complexes/` and `interactions/` exist only on the DAG path. The comparative node also copies 16 files from `reports/` into `reports/comparative/` and `reports/biology_correlation/`. Those copies are a duplicate-writer symptom.

## 5. Single-writer map for shared `reports/` files (current)

| File or directory | Writers | Problem |
|---|---|---|
| `reports/summary.txt` | comparative writer (line ~3543). Called by the comparative node (line 1288), the polypharmacology fallback (1318) and the biology fallback (1328). The DAG `reports` node has the same output path (`_dag_artifact_paths()["reports"]`, line 584). | Two owners of one file. The graph retires the writer's `summary.txt` into `.artifact_history/` before the `reports` node runs (`.artifact_history/<id>/0/summary.txt` is present in every DAG capture). This is the D-A1 cause. |
| `reports/polypharmacology/`, `reports/biology_correlation/` | comparative writer. The polypharmacology and biology nodes have the same output directories. | Outputs are retired into `.artifact_history` and rewritten. Retirement history varies between identical runs. |
| `reports/*.csv`, `reports/*.json`, `VALIDATION_SUMMARY.md` | comparative writer, copied to `reports/comparative/` by the comparative node | Copies are duplicates of files in `reports/` |

**Target (after decisions D-A1 and D-A3):** the comparative node is the only caller of the comparative writer. The polypharmacology, biology and reports nodes only read. Their fallback calls are removed. The reports node owns `summary.txt`. The DAG graph must not retire files that another node wrote. Expected effect on golden outputs: no change to scientific files, and no `.artifact_history` churn for `reports/`. The summary text changes only if D-A1 is approved.

## 6. Decisions for the Lead (STOP)

1. **D-A (route `run()` through the DAG).** Routing `run()` through the DAG changes the legacy outputs (D-A1, D-A2, D-A3, D-B). The Lead needs to choose one of these:
   - (a) Keep `run()` as it is and mark it as a compatibility path. Recommended until (b) is ready.
   - (b) Fix the DAG first, then route `run()`. The fixes are: the summary is restored (D-A1); the atlas carries `normalized_affinity_score` and the status columns (D-A2); the favorite stage becomes a DAG node or is retired (D-A3). Golden comparison is the acceptance test.
   - (c) Accept the DAG output as it is. Not recommended: it replaces the full summary with a placeholder and drops status columns.
2. **D-B (problematic-ligand filter).** Proposed: a DAG node `problematic_ligand_filter`, after `raw_scores`, that writes `excluded_problematic_ligands.csv` and records `exclude_problematic_ligands`, `positive_affinity_threshold` and `minimum_pose_count` in `analysis_parameters.json`. Default stays off. Enabling it changes ranks on both paths, so this needs sign-off.
3. **D-C (best-pose rule for the single-engine GNINA path).** `_best_pose_per_protein_ligand` uses idxmin on `vina_affinity`, which is a different rule from `pose_selection`. The Lead should confirm which rule applies before the two are merged.
4. **E-1 (legacy smina crash).** Proposed one-line fix: initialise `self.engine_detection_report = {}` in `MultiEngineAnalysisPipeline.__init__`. Not made, because it changes a crash into outputs. If approved, add a golden comparison for smina.

Until these are decided, the duplicate removal in section 3 is not started, because removing those duplicates would have to preserve the outputs above.

## 7. Test file

`test/test_spec038_postdocking.py`: 16 tests.
- Passing (6): the compatibility entry points import; the filter default is off; consensus outputs are identical between the two paths; single-engine values are identical; the legacy favorite directory is written; the legacy filter is applied when enabled and not applied when disabled.
- Strict xfail (10): one per divergence above (D-A1, D-A2 twice, D-A3, D-B, E-1), the routing check, the single-writer checks (node source and retired history), and the whole-session equality check.
- Each xfail is strict. When a fix lands, its test fails, and it must then be flipped to a plain test.

The new module was run on the current tree with the other agents' changes present. Result: 6 passed, 10 xfailed.

## 8. Other notes

- Non-scientific nondeterminism: the seaborn jitter (section 1). Not fixed, because seeding would change the figures.
- A stray empty directory `/spec038_harness` at the filesystem root was created by mistake during the golden setup. A safety check blocked its removal. It should be removed with `rmdir /spec038_harness`. The harness itself lives in the scratchpad.
