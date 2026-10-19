# Spec 035 inventory: `main` (7fc130a, 8e36027) vs `chore/install-agentic-bioinformatics-stack`

**Status:** read-only analysis, 2026-10-19. No code, merge, rebase, cherry-pick or commit was made.
**Branch HEAD:** 64aebd6. **Merge base:** cfe3305. **main:** `origin/main` = 8d4434d, which contains 8e36027 and 7fc130a.
**Trial merge:** `git merge-tree --write-tree HEAD origin/main` (no write).

## 0. Findings that change the spec

1. **33 files conflict, not 15.** 32 content conflicts and 1 add/add (`docking/preparation/receptor_preparation.py`). The list is in section 2 (rows marked "conflict").
2. **Silent takeovers.** A plain `git merge` applies these with no conflict, and several change scientific behaviour:
   - Files the branch never touched: `docking/preparation/ligand_quality.py` (S14), `post_docking_analysis/prolif_interaction_maps.py` (O20), `biology_integration.py` (S25), `docking_parser.py` (S24), `docking/preflight.py`, `docking/remote_ops.py`.
   - Files the branch touched that auto-merge: `enhanced_rmsd_analyzer.py` (S5, Kabsch removed), `workflow/execution.py` (S15, raw PDB to SDF conversion removed), `simplified_pipeline_impl.py` (S26), `generate_scores_csv.py` (O18), `report_generator.py` (O21), `docking/deployment.py` (O9), `workflow/interactive.py` (O7).
   - Conclusion: the final merge must be done unit by unit, as Spec 035 says, not with a single `git merge`.
3. **Already present.** `docking/preparation/structure_contract.py` is byte-identical on both sides. `receptor_preparation.py` differs and the branch's Spec 033 version is authoritative.
4. **Parallel implementations.** `main` adds `pose_geometry.py`, `asset_identity.py` and `score_semantics.py` next to branch code that already implements the same accepted methods (`atom_mapping.py`, `reference_policy.py`, consensus v2). Those `main` modules are alternatives, not extensions.

## 1. Summary (counts by action)

| Action | Units | Rows |
|---|---|---|
| port | 27 | O1–O22, O24–O28 (operational, tests; start with O2–O5) |
| already_present | 1 | S17 |
| keep_branch | 11 | S1, S3–S6, S9, S11, S16, S21, O23, O30 (accepted method wins) |
| compare_then_decide | 18 | S2, S7, S8, S10, S12–S15, S18–S20, S22–S28 (needs the Lead, section 6) |
| review_public | 1 | O29 (`docs/correctness-update.md`); the `audit/` findings are in section 4 |

Total units: 58 (O1–O30, S1–S28). Counts were taken from the tables in section 2.

## 2. Unit inventory

Conflict = listed in the trial merge (section 0). "Spec/decision" is the accepted method the unit overlaps. Tests are from `main`'s `test/test_*_correctness.py` unless noted.

### 2a. Operational, tests and docs

| ID | Class | main files / functions | Branch equivalent and overlap | Merge | Action | Verifying test / contradicts accepted method? |
|---|---|---|---|---|---|---|
| O1 | operational | `.github/workflows/correctness.yml`, `pyproject.toml`, `pytest.ini`, `MANIFEST.in`, `scripts/check_wheel.py`, `setup.py` (`BuildWithPreparationScript`, `find_packages` excludes `audit`, `python_requires` 3.10) | none; branch `setup.py` changed 2 lines | new; setup.py auto | port | `python -m build && python scripts/check_wheel.py`. Confirm every `py_modules` entry exists before porting `setup.py`. No contradiction. |
| O2 | operational | 8e36027: `requests>=2.31.0` in both requirements files; `scikit-learn` in root | branch pins are stricter: `rdkit>=2024.3.1`, `openbabel==3.2.1`, pdb2pqr, meeko, `networkx<4`. Keep those. `requests` is imported by `post_docking_analysis/poseview_integration.py` on both sides | conflict | port `requests` only | `python -c "import requests"`. No contradiction. |
| O3 | operational | `.gitignore` (`.test-tmp*/`, `pytest-cache-files-*/`) | none | auto | port | none. |
| O4 | operational | `workflow/locking.py` (`file_lock`, `_thread_lock`) | none | new | port | `test_workflow_correctness.py::test_lock_timeout_also_bounds_competing_threads`. |
| O5 | operational | `workflow/state.py` (`_state_transaction` on 9 functions; checkpoints newest-first via `reversed`) | branch state updates are unlocked | conflict | port | `test_state_concurrent_processes_preserve_all_updates`, `test_checkpoint_order_is_newest_first`. |
| O6 | operational | `post_docking_analysis/artifact_graph.py`: lock around execute, `_atomic_json`, `_paths_fingerprint`, content-cache-v2 key (code and dependency versions), `_retire_outputs` (moves old outputs to `.artifact_history`), output fingerprint | branch artifact DAG (Spec 024/026); branch changed 10 lines | conflict | port | `test_cache_detects_*`, `test_cache_rejects_modified_output`. Behaviour change: any edit to a package `.py` invalidates all caches. |
| O7 | operational | `workflow/interactive.py`: symlink replace and `target_is_directory=True` | branch 218 lines changed | auto | port | smoke test. |
| O8 | operational | `docking/runners/job_contract.py` (new, stdlib only): `execute_job`, `completion_valid`, `validate_pose_file`, `parse_sdf_scores`, `execution_identity`, `file_hash`, `write_json` | none. Overlaps branch output acceptance (Spec 031 ingestion) | new | port | `test_invalid_output_and_startup_failure_never_complete`, `test_multimodel_sdf_and_ad4_completion_validation`, `test_resume_requires_exact_inputs_protocol_and_validated_output`. |
| O9 | operational | `docking/deployment.py`: `_execution_root`, `_mapped_path`, `_map_job`, `_write_job_executor`, `_apply_allocated_resources`. `docking/models.py`: `engine`, `fingerprint`, `completion_file`, `input_files`, `executables`, `skip_completed` | branch HPC deployment | deployment auto; models conflict | port after O8 | `test_portable_deployment_uses_shared_executor_and_posix_paths`. |
| O10 | operational | `docking/runners/base.py`: `accepted_pose_files` helper and completion checks. `runners/vina.py`: poses via `accepted_pose_files` | branch `base.py` has the Spec 036 replicate/ensemble code (260 lines changed) | conflict | port, manual merge around replicate code | `test_changed_input_after_plan_does_not_execute`, `test_empty_legacy_pose_does_not_resume`, `test_failed_rerun_retires_old_pose_and_clears_scores`. |
| O11 | operational | `runners/gnina.py`: `--no_gpu` when `use_gpu` is false; `--cpu` independent of GPU device. Fixes audit D6 | branch never passes `--no_gpu` (grep: absent) | conflict | port | `test_gnina_resource_flags_independent`. |
| O12 | operational | `runners/autodock4.py`: staged receptor, ligand and parameter copies with fixed names; `cd pair_dir`; `spacing:.12g`; sha256 of inputs. Grid rule is S28 | branch uses absolute paths | conflict | port (staging only) | `test_ad4_grid_and_parameter_files_are_portable_and_consistent`: its `npts` assertion conflicts with branch (S28). |
| O13 | operational | `docking/parameter_schema.py`: `validate_pairlist_geometry`, finite `box_scale` and `box_padding`. `docking/cli.py`: return 1 when any job failed | branch `cli.py` 478 lines changed | conflict | port | `test_nonfinite_box_schema_rejected`, `test_cli_returns_failure_for_failed_engine`. |
| O14 | operational | `docking/project_layout.py`: symlink `OSError` fallback (README plus `directory_fallback`) | branch Spec 033 layout (829 lines changed) | conflict | port (one hunk) | `test_symlink_privilege_failure_has_usable_fallback`. |
| O15 | operational | `docking/preparation/excel_sites.py`: reject duplicate PDB/property rows; ±inf to NaN; size columns must be positive | branch `excel_sites.py` and `binding_site_center.py` (Spec 031 grid center) | conflict | port guards only | `test_spreadsheet_rejects_duplicate_sites_and_nonfinite_centers`. |
| O16 | operational | `docking/preparation/asset_identity.py` `index_assets`: refuse ambiguous aliases. `active_assets`: prefer `.pdbqt` | branch `pairlist_builder.py` 293 lines changed | new; conflict | port the ambiguity refusal. The frame check is S18 | `test_builder_selects_cleaned_receptor_and_rejects_ambiguous_conformations`. |
| O17 | operational | `visualization_suite.py`: `applymap` replaced by `apply` (pandas deprecation) | branch changed 96 lines | conflict | port | smoke test. |
| O18 | operational | `generate_scores_csv.py`: skip `*.runner.log`; logging instead of print | auto | auto | port | `test_runner_and_unlisted_logs_do_not_enter_scores`. |
| O19 | operational | `engine_detector.py`: pose files limited to pairlist mapping (diff read) | branch `engine_detector.py` 268 lines changed | conflict | port after a check against Spec 022 engine detection | smoke test. |
| O20 | operational | `prolif_interaction_maps.py`: pose iteration changed (skimmed, not verified line by line) | branch untouched, so auto | auto | port after line review | smoke test. |
| O21 | docs | `report_generator.py`: "MD-ready" becomes "Receptor-pose complexes (no simulation topology)" | Spec 032 wording. Consistent with Spec 032 `reported_unverified` | auto | port | check `test_spec032_md_inputs.py` for the old wording. |
| O22 | scientific (definitions) | `score_semantics.py` (new): `SCORES` table, `sort_value`, `explicit_true` | branch has the same directions inline (`pose_selection.py`, consensus `_safe_rank_score`). Values match Spec 031 | new | port | `test_all_normalizers_preserve_scientific_direction`. |
| O23 | tests | `test/test_dockforge_smoke.py` (`_mark_score_import` added) | branch file is much larger (519 lines changed) | conflict | keep_branch; hand-port added contracts | smoke test. |
| O24 | tests | `test/test_docking_correctness.py` (24 tests) | none | new | port after edits | Contradicts branch: `npts 54 54 54` (lines 205, 220) versus branch `round()`; see S28. |
| O25 | tests | `test/test_workflow_correctness.py` (15 tests) | none | new | port | itself. |
| O26 | tests | `test/test_postdocking_correctness.py` (34 tests) | none | new | port after edits | Contradicts Lead decision 3: legacy modes at lines 206–207, 215–217, 223, 287. Its atlas test (about line 287) depends on S8. Also uses removed `compute_geometric_consensus`. See section 3. |
| O27 | tests | `test/test_preparation_correctness.py` (20 tests) | none | new | port after decisions | Contradicts branch: pocket metrics must be `None` (lines 112–119; see S20). Expects `RuntimeError` on pH failure with no fallback (line 181; see S13). |
| O28 | tests | `test/test_reference_chemistry.py` (4 tests) | none | new | port local tests; mark the network test skipped | itself. Depends on S19. |
| O29 | docs | `docs/correctness-update.md` (181 lines) | none | new | review_public | Line 42 says "AutoDock4 grids round outward", which contradicts branch `round()` unless D2 changes it. |
| O30 | docs | `README.md` (main adds favorite-engine text) | branch README has its own updates | auto | keep_branch; port hunks after decisions | none. |

### 2b. Scientific units

| ID | Class | main files / functions | Branch equivalent and overlap (spec/decision) | Merge | Action | Test / contradiction |
|---|---|---|---|---|---|---|
| S1 | scientific | `consensus.py`: `SUPPORTED_CONSENSUS_MODES` = `dockbox_geometric` (default), `weighted_hybrid`, `strict_consensus`, `favorite_guardrails` | removed on branch (Spec 036 R3a; Lead decision 3) | conflict | keep_branch | Contradicts Lead 3 in `test_postdocking_correctness.py`. Record, do not "fix" the branch. |
| S2 | scientific | `consensus.py` normalizers: missing becomes NaN, not 0; direction via `score_spec`; grouping by engine, protein, site and scoring function; winner by normalized percentile; raw mean and spread blanked across engines | branch `normalize_engine_scores` still fills 0.0 (`consensus.py:256–263`). Spec 031 forbids neutral imputation | conflict | compare_then_decide (D16) | `test_engine_score_offset_does_not_change_consensus_winner_or_rank`. |
| S3 | scientific | `geometric_consensus.py`: `compute_geometric_consensus` uses fixed-frame RMSD and a receptor-frame check. `_kabsch_rmsd` is a misnamed, non-fitting RMSD. `_pairwise_rmsd` is a stub | branch `compute_geometric_consensus_v2` (graph mapping, Kabsch; Spec 031 item 4; Spec 034 "Geometric consensus keeps Kabsch") | conflict | keep_branch | Tests call the removed API. The identical-pose tests pass under both frames. |
| S4 | scientific | `pose_geometry.py` (new): `pose_rmsd`, `fixed_frame_rmsd` (RDKit `CalcRMS` over all substructure matches), `load_pose_molecule` (Meeko `SMILES IDX`), `selected_record`, `content_hash` | branch `atom_mapping.py`: `compare_graph_poses_in_place` (`in_place_v1`), Kabsch secondary (Spec 031 mapping; Spec 034 T002a) | new | keep_branch. Port its edge-case tests against `atom_mapping` | `test_postdocking_correctness.py` 40–71 and the PDBQT/Meeko tests. |
| S5 | scientific | `enhanced_rmsd_analyzer.py`: Kabsch helper becomes a plain RMSD ("preserves the supplied frame") | branch Spec 036 R3d: labelled `diagnostic_not_spec031`, keeps Kabsch | auto (silent) | keep_branch | none. |
| S6 | scientific | `redocking_validation.py`: `_compute_pose_rmsd` is in-place only; `_experimental_provenance`; `_is_reference_candidate_row` is an explicit flag or a `pair_source` token; no fuzzy file match; receptor-frame check; missing-engine rows added | branch: in-place primary plus Kabsch secondary (Spec 034 T002a); `classify_reference_row` (`reference_policy.py:76`, Spec 031 `explicit_reference_metadata_v2`); fuzzy lookup only in legacy mode | conflict | keep_branch for RMSD frame and reference classification | `test_redocking_uses_selected_pose_and_requires_provenance`, `test_reference_discovery_never_uses_raw_conformer`. |
| S7 | scientific | redocking gate: `SUFFICIENT_PROTOCOL_CHECKS` when coverage is 1.0, at least 3 unique experimental complexes, and pass rate at least 0.80 | branch: HIGH/MEDIUM/LOW by pass rate; `validated` only when HIGH | conflict | compare_then_decide (D8) | `test_sparse_validation_does_not_claim_high_confidence` reads `gate["coverage"]`, which the branch does not emit. |
| S8 | scientific | `top_pose_selector.py`: multi-engine primary is the within-engine normalized score; winner engine breaks ties; raw affinity never breaks cross-engine ties | branch default `best_affinity` (min raw affinity across engines) unchanged; listed as open in Spec 036 | conflict | compare_then_decide (D9) | `test_atlas_cross_engine_identity_and_rank_ignore_raw_score_offset` uses removed `weighted_hybrid`. |
| S9 | scientific | `multi_engine_pipeline_impl.py`: DAG consensus inputs, `_run_dag_consensus_strategy`, per-engine winner, reference baselines per engine and scoring function, GNINA `cnn_affinity` sorted descending, `score_spec` idxmin/idxmax | branch Spec 036 R3b and R4 (`pose_selection.select_best_pose_rows`), v2 summaries | conflict | keep_branch | `test_gnina_selectors_prefer_higher_predicted_pk` (lines 199–204) selects on `cnn_affinity`. Spec 031 chooses the GNINA representative by `cnn_score`. Record as a conflict. |
| S10 | scientific | `multi_engine_pipeline_impl.py` `_build_pair_competition`: raises when `min_affinity_advantage` is above 0 across engines; `affinity_advantage` set to None | branch promotes on raw advantage (Spec 036 open list) | conflict | compare_then_decide (D10) | none in `main` tests. |
| S11 | scientific | `pose_extractor.py`: `sort_value` direction-aware selection; `pose_number` guard | branch `pose_selection.py` (Spec 031 GNINA rule, v2 ranking direction) | conflict | keep_branch | none. |
| S12 | scientific | `score_import.py` (new): `read_explicit_score_import`, hash-checked, marked unverified | none | new | compare_then_decide (D13) | `test_explicit_import_is_hash_checked_and_marked_unverified`, `test_explicit_import_rejects_invalid_identity_score_and_pose`. |
| S13 | scientific | `ligand_preparation.py`: one-molecule SDF supplier; refuses PDB; keeps 3D input (no gen3d, no MMFF/UFF); no RDKit AddHs fallback; OpenBabel pH failure raises | branch `normalize_ligand_to_sdf` falls back to RDKit when `require_ph_model` is false (`ligand_preparation.py:898–911`), recorded as `rdkit_fallback_no_ph_model`. Spec 034 protonation | conflict | compare_then_decide (D3, D4) | `test_ph_failure_does_not_fall_back_to_add_hydrogens`, `test_pdb_cannot_be_treated_as_authoritative_ligand_chemistry`, `test_missing_openbabel_does_not_substitute_untreated_rdkit_state`. |
| S14 | scientific | `ligand_quality.py`: `repair_project_ligands` no longer guesses atom types (`_infer_fixed_atom_type`, `_replace_atom_type` removed) | branch unchanged, still guesses (`ligand_quality.py:174–178`) | auto (silent) | compare_then_decide (D6) | `test_unknown_type_repair_cannot_delete_chemistry`. |
| S15 | scientific | `workflow/execution.py`: no OpenBabel PDB to SDF conversion (an inferred bond graph is not validated); an old SDF is moved to history | branch still converts (`execution.py:692`). Spec 031: PDB coordinates do not establish bond order | auto (silent) | compare_then_decide (D5) | Reference-chemistry tests. |
| S16 | scientific | `receptor_preparation.py`: PDB2PQR without `--noopt --nodebump`; `receptor_use_pdb2pqr` defaults to false; strict conservation | branch Spec 033 R2c accepted rule (`PDB2PQR_OPTIONS`, `receptor_preparation.py:57`) | add/add | keep_branch (Spec 033, option 1) | `test_receptor_backend_cannot_silently_drop_heavy_atoms`: check it against the branch's terminal-addition rule. |
| S17 | scientific | `docking/preparation/structure_contract.py` | identical | identical | already_present | none. |
| S18 | scientific | `asset_identity.py` `pair_reference_metadata`: copies reference fields from sidecars and clears them when `reference_frame_id` differs from `receptor_frame_id`. `project_builder.py` reference handling | branch `reference_policy.py` and `pair_curation.py` (Spec 031 explicit reference) | conflict (`pairlist_builder`, `project_builder`) | compare_then_decide (D17) | `test_project_build_preserves_verified_native_reference_and_frame`, `test_reference_metadata_does_not_mislabel_cross_docking`. |
| S19 | scientific | `workflow/reference_chemistry.py` (new): `retrieve_reference_sdf` fetches a ligand model SDF from models.rcsb.org over HTTPS | none. Spec 031 `skipped_network`; AGENTS section 5 (retrieval separate from execution) | new | compare_then_decide (D14) | `test/test_reference_chemistry.py`. |
| S20 | scientific | `core_pipeline.py` `analyze_pocket_properties`: pocket volume and druggability are `None` / "Not evaluated" | branch returns labelled uncalibrated heuristics (`core_pipeline.py:764–871`). AGENTS allows labelled heuristics | conflict | compare_then_decide (D7) | `test_unsupported_pocket_metrics_are_not_fabricated` asserts `None`. |
| S21 | scientific | `core_pipeline.py` `clean_pdb` and `extract_residue_level_coordinates`: mostly restructuring (not diffed line by line) | branch version | conflict | keep_branch | none. |
| S22 | scientific (QC) | `multi_engine_pipeline_impl.py` `_normalize_pdb_record(line, line[:6], line[21:22])` keeps record name and chain | branch hard-codes `"ATOM  "` and chain `"A"` (`multi_engine_pipeline_impl.py:4868, 4870`). Spec 037 multi-entity | conflict | compare_then_decide (D15) | none in `main` tests. |
| S23 | scientific (QC) | Vina `rmsd_lb` collapse warning removed (mode 1 is 0 by definition) | branch warns (`multi_engine_pipeline_impl.py:4554–4558`, `visualization_suite.py:328`) | conflict | compare_then_decide (D15) | none. |
| S24 | scientific (QC) | `docking_parser.py`: missing RMSD becomes NaN, not 0.0 | branch uses 0.0 | auto (silent) | compare_then_decide (D15) | `test_missing_smina_rmsd_is_not_zero`. |
| S25 | scientific (QC) | `biology_integration.py`: keys must be non-missing and unique; `many_to_one` validated (diff read) | branch allows duplicate keys | auto (silent) | compare_then_decide (D15) | `test_biology_replicates_cannot_duplicate_candidates`. |
| S26 | scientific | `simplified_pipeline_impl.py`: cross-engine polypharmacology writes `not_evaluable` status | branch: Spec 036 open list | auto | compare_then_decide (D10) | none. |
| S27 | scientific | `docking/cli.py` `_restore_saved_protocol`: saved scientific settings override defaults; an explicit `--exhaustiveness` switches to `advanced` | branch Lead decision 5: basic mode refuses explicit `--exhaustiveness` | conflict | compare_then_decide (D18) | none in `main` tests. |
| S28 | scientific | `runners/autodock4.py` `_grid_points`: `2*ceil(size/(2*spacing))` (even count). Box in Å unchanged | branch: `int(round(size/spacing))` (`autodock4.py:72–73`). Spec 031 says box dimensions are unchanged | conflict | compare_then_decide (D2) | `test_ad4_grid_and_parameter_files_are_portable_and_consistent` expects `npts 54 54 54`. For a 20 Å box at 0.375 Å, branch gives 53. |

## 3. Side-by-side for compare_then_decide units

- **S2 normalizers.** Branch: `normalize_engine_scores` fills missing values with 0.0 (`consensus.py:256–263`). Main: missing values stay NaN, and the table is grouped by engine, protein, site and scoring function. The function feeds `normalized_scores` (`multi_engine_pipeline_impl.py:2874`, `:5236`). v2 composites use rank percentiles, so they are not changed directly. Check each consumer of the zero-filled column before porting.
- **S7 redocking gate.** Branch: pass rate at least 0.80 gives `HIGH`, which is `validated`. Main: also requires coverage 1.0 and at least 3 unique experimental complexes, and emits `SUFFICIENT_PROTOCOL_CHECKS`. The thresholds 1.0 and 3 are new numbers.
- **S8 top pose.** Branch: `best_affinity` ranks by minimum raw kcal/mol across engines. Main: ranks by within-engine normalized score, and the v2 winner breaks ties.
- **S10 and S26 cross-engine raw comparisons.** Branch: promotes on raw advantage and writes polypharmacology. Main: raises for a positive advantage threshold across engines and writes `not_evaluable`.
- **S13 ligand preparation.** Branch: an RDKit fallback runs when the pH model is not required. Main: any pH failure raises, and PDB ligands are refused.
- **S14 ligand quality.** Branch: guesses unknown atom types and rewrites them. Main: refuses to rewrite; unknown types fail.
- **S15 raw ligand conversion.** Branch: OpenBabel converts raw PDB ligands to SDF with inferred bonds. Main: refuses, because bonds are not validated.
- **S18 reference frame.** Branch: reference status comes from explicit fields. Main: also clears a reference whose frame differs from the receptor.
- **S19 RCSB retrieval.** Main: fetches a reference SDF over the network. Branch: no network retrieval.
- **S20 pocket.** Branch: `pocket_volume_A3 = 4/3 pi 5^3`, labelled heuristic. Main: `None`, "Not evaluated".
- **S22 to S25 (QC).** Main: Vina missing RMSD is NaN; no collapse warning at mode 1; chains kept; biology keys unique. Branch: 0.0; warning; chain "A" hard-coded; duplicates allowed.
- **S27 protocol restore.** Main: saved settings override defaults, and explicit exhaustiveness switches to advanced. Branch: basic mode refuses explicit exhaustiveness (Lead 5).
- **S28 AD4 grid.** Branch: 53 points for 20 Å at 0.375 Å. Main: 54 points. Both keep the box at 20 Å.

## 4. Public-repository review (review_public)

Import no file in `audit/` as-is. Findings (file:line; counts are full-file totals from `git show 7fc130a`):

- `audit/build-correctness.log`: 6 lines with `<local Windows user path>` (2, 26, 198, 222, and two more). Local Windows path containing the operator's username.
- `audit/correctness-tests.log`: 147 lines with the same path, including line 2 (interpreter), line 962 (XML output path), and `AppData\Local\Temp` paths. Lines 53–54 name NMRbox HTCondor tests (names only).
- `audit/correctness-tests.xml`: line 1 (the whole file is one line). `a machine hostname` and a timestamp with timezone.
- `audit/docking-dry-run.log`: line 10, `<local Windows user path>`.
- `audit/docking-probe-results.json`: 32 lines with `<local Windows user path>` (e.g. 21, 22, 30). The probe-data directory is ignored, so these are local run paths.
- `audit/integrated-technical-results.json`: lines 15, 16, 25 with `<local Windows user path>`.
- `audit/preparation-audit.md:127`: "atom-type/root/TORSDOF" matches a `/root/` search. False positive, not a path.
- Not private, but environment-specific: `test/test_dockforge_smoke.py` (main) at 2607, 2611, 2626, 2646 hard-codes `/home/user/docking` and `/home/user/gnina`. Replace with `tmp_path` before porting.
- `docking/hpc_profiles.py` and `examples/hpc_profiles/*.template.*` use `<user>@sulfur.nmrbox.org` as a public placeholder. Already in the base commit. Acceptable.
- No unpublished scientific results were found: no 1IEP, STI or NMRbox project data, and no Notion or Omar content. Audit files reuse synthetic PDB IDs such as `1ABC` (`audit/preparation-audit.md`, `audit/preparation-probes.py`).

Required action: import the findings and conclusions from `AUDIT.md`, `preparation-audit.md`, `postdocking-audit.md`, `docking-audit.md` and `technical-audit.md` as sanitized text, with no logs, XML or JSON that carry the Windows paths. Regenerate any log that is kept, from a clean run on a neutral path. `audit/.gitignore` allows `correctness-tests.log`, `build-correctness.log`, `wheel-verification.log` and `docking-dry-run.log` through, so those are the ones at risk.

## 5. Ordered port plan

Each step is a small commit. Run it and confirm green before the next. Use the branch's own commands from `CLAUDE.md`, plus the `main` tests named above.

0. **Decisions** (section 6), at least D1, D12 and D19.
1. **Dependencies (O2, O3).** Verify: `python -c "import requests"`; `python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix`.
2. **Locking and state (O4, O5).** Verify: `pytest test/test_workflow_correctness.py -k "state or lock or checkpoint"`, then `pytest test/test_webui.py`.
3. **Artifact cache (O6, O7).** Verify: `pytest test/test_workflow_correctness.py -k cache`; `python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix`.
4. **Job contract (O8, O9, O10).** Verify: `pytest test/test_docking_correctness.py -k "deploy or resume or completion or multimodel"`; `pytest test/test_spec036_*.py test/test_spec031_*.py`.
5. **Runner fixes (O11, O13, O14, O18, O22).** Verify: `pytest test/test_docking_correctness.py -k "gnina or nonfinite or cli or symlink"`; `pytest test/test_spec031_*.py`.
6. **Preparation guards (O15, O16, O17, O19, O20, O21).** Verify: `pytest test/test_preparation_correctness.py -k "spreadsheet or ambiguous"`; `pytest test/test_spec033_*.py`.
7. **Packaging and CI (O1).** Verify: `python -m build && python scripts/check_wheel.py && python -m pytest test -q`.
8. **AutoDock4 staging (O12)**, after D2 is decided. Verify: `pytest test/test_docking_correctness.py -k ad4`.
9. **Scientific ports** (only the items the Lead approved in section 6). One commit per decision, each with its own `main` test case adapted to the accepted method. Verify each with the affected `test/test_spec0*` suite and `python test/test_dockforge_smoke.py`.
10. **Test files** (O24–O28). Port the files after steps 1–9. Remove the legacy-mode parametrisations, re-target the removed APIs, and mark the network test skipped. Record any remaining contradiction, do not change the branch to pass it.
11. **Merge.** Run `git merge origin/main` only when steps 0–10 are done. Re-run `git merge-tree --write-tree HEAD origin/main` and confirm that every conflict is resolved toward the decisions. Verify: `pytest test/test_webui.py`, `python test/test_dockforge_smoke.py`.

## 6. Decisions for the Lead

Each question is yes/no or A/B. The recommendation is in brackets.

- **D1.** Approve the phased port order in section 5, operational first? [yes]
- **D2.** AutoDock4 grid points: A = even ceiling `2*ceil(size/(2*spacing))` (main; 54 for 20 Å at 0.375 Å), or B = `round(size/spacing)` (branch; 53)? Box in Å is unchanged either way. [A, once the AutoGrid even-point convention is checked against its documentation.]
- **D3.** Ligand chemistry when the pH step fails: A = refuse (main), or B = continue with a labelled RDKit fallback (branch)? [A]
- **D4.** Refuse PDB-only ligand input, and stop re-minimising 3D input? [yes to refusing PDB; keep the branch's 3D handling until checked]
- **D5.** Remove the raw PDB to SDF conversion with inferred bonds in `execution.py`? [yes; Spec 031 says PDB coordinates do not establish bond order]
- **D6.** Remove the atom-type guessing in `ligand_quality.py`? [yes]
- **D7.** Pocket volume and druggability: A = labelled uncalibrated heuristic (branch), or B = "Not evaluated" (main)? [B]
- **D8.** Redocking gate: A = branch `HIGH` at pass rate 0.80, or B = main coverage 1.0, at least 3 unique experimental complexes, pass rate 0.80? [B, with the numbers 1.0 and 3 confirmed by the Lead]
- **D9.** Top-pose default: A = keep `best_affinity` (minimum raw affinity across engines), or B = within-engine normalized primary with the v2 winner? [B; Lead decision 6 and the score-incommensurability rule]
- **D10.** Cross-engine raw comparisons (advantage promotion, polypharmacology): A = keep, or B = refuse with `not_evaluable`? [B]
- **D11.** Confirm inter-engine geometric agreement stays Kabsch, per Spec 034? [yes]
- **D12.** Confirm Lead decision 3 reading: the legacy consensus modes stay removed, and the Spec 031 "preserve legacy names" note is superseded? [yes; this was requested in `consensus_inventory.md`]
- **D13.** Port the explicit score import (`score_import.py`) now? [no; defer to its own spec]
- **D14.** Port the RCSB reference retrieval (`reference_chemistry.py`) now? [no; needs an approved network policy]
- **D15.** QC group, four yes/no items: Vina missing RMSD becomes NaN; no collapse warning at mode 1; keep record name and chain in complex export; biology keys must be unique. [yes to all four]
- **D16.** Normalizers: missing values become NaN rather than 0.0 in the legacy normalizers? [yes; Spec 031 forbids neutral imputation]
- **D17.** Clear a reference whose frame differs from the receptor's frame? [yes]
- **D18.** Saved protocol vs explicit flags (`cli.py`): A = saved protocol overrides defaults (main), or B = Lead 5, basic mode refuses explicit `--exhaustiveness` (branch)? [B]
- **D19.** Import `audit/` into the public repository? [no as files; import only the sanitised findings in section 4]

## 7. Corrections and open items

- Section 0 replaces the spec's "15 files" with 33.
- Units not diffed line by line: O20, S21, and the `main` changes to `cli_pipeline.py`, `main.py`, `batch_pdb_preparation.py`, `autodock_preparation.py`, `test_pipeline.py` and `prep_autodock_enhanced.sh` (keep_branch by default). Review them before the final merge.
