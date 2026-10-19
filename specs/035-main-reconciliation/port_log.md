# Spec 035 port log (M1 operational, M2 decided scientific)

**Source:** `main@7fc130a` with `8e36027` (fetched `origin/main` = `ec327f4`, which removes the private audit logs). No `audit/` file was imported (D19; the directory is absent on this branch).
**Branch:** `chore/install-agentic-bioinformatics-stack`, from HEAD `7c3f34c`.
**Method:** hand port by reading `git show 7fc130a -- <file>` and `8e36027`. No merge, rebase or cherry-pick. Nothing pushed.
**Verification:** focused tests before each commit; full suite after every 3 to 4 commits; the contract smoke test (`--skip-all-engines --skip-prep-matrix`) before each code commit.

## 1. Commits (oldest first)

| Unit | Commit | Subject |
|---|---|---|
| O2 | a3272f3 | declare requests dependency |
| O3 | 4629765 | ignore pytest temporary directories |
| O4, O5 | 5b8da11 | cross-process workflow state locking |
| O6 | f3644d0 | artifact cache: locking, atomic writes, content keys |
| O8, O10 | 7366d1f | job contract: hashed inputs, completion-gated resume |
| O9 | 8ecd6f7 | portable HPC deployment with the shared job executor |
| O4 | 81e5480 | ignore workflow lock sidecars in test fixtures |
| O11 | b9db547 | gnina passes --no_gpu and --cpu independently of the device |
| O13 | c2c7922 | reject invalid box geometry; failed jobs exit non-zero |
| O14 | 193950e | explain the directory fallback when symlinks are unavailable |
| O18 | 13bf28a | exclude scheduler runner logs from gnina score tables |
| O15 | 9b39097 | refuse duplicate site rows and non-finite site values |
| O16 | 89688fe | refuse ambiguous prepared-asset aliases |
| O21 | 27a3fca | report complexes as receptor-pose complexes, not MD-ready |
| O7 | 604d7c1 | latest-session shortcuts: never delete a directory, text pointer fallback |
| O1 | 0444c3a | packaging, wheel check and CI |
| D2 | 64d1208 | AutoDock4 grid points use the even ceiling |
| O12 | 3f9a45d | AutoDock4 stages portable inputs in the pair directory |
| D3, D4 | 2234f50 | refuse failed pH and PDB ligands; keep 3D input as supplied |
| D5 | 48be806 | no SDF from PDB coordinates; stale SDFs go to history |
| D6 | 8fceddf | no atom-type guessing; a repair never converts raw PDB |
| D7 | 8dc9fd1 | pocket volume and druggability are "Not evaluated" |
| D8 | 9547c31 | redocking gate: full coverage, 3 complexes, pass rate 0.80 |
| D15 | bfb6125 | QC group: NaN RMSD, no Vina collapse warning, record/chain, unique biology keys |
| D16 | 0bf63fb | legacy normalisers keep a missing score as NaN |
| D17 | 44d24ef | a reference from another frame is cleared |
| O18 | 74cdf82 | gnina logs the pairlist does not list never enter scores (reverted next) |
| O18 | 5b12017 | revert the unlisted-log rule; keep the Spec 031 directional matching |
| tests | a66ca7c | main's correctness test files, adapted to accepted behaviour |

Commits 74cdf82 and 5b12017 are kept on purpose: the first was a port that broke an accepted smoke contract, and the revert records why.

## 2. Per unit

Status key: **ported**, **ported (partial)**, **present** (the branch already has it), **kept** (branch version stays), **deferred**, **not ported**.

### M1 operational

- **O1 packaging and CI: ported.** `setup.py` (py_modules, package_data, the shell-script build step, `python_requires >= 3.10`, extras; the branch's console entry points are kept), `pyproject.toml`, `MANIFEST.in`, `pytest.ini`, `scripts/check_wheel.py`, `.github/workflows/correctness.yml`. The `pytest.ini` `_smoke_*` pattern collects the 98 smoke contracts; they pass under pytest. Verified: `python -m build` and `python scripts/check_wheel.py` (exit 0). Build outputs removed from the tree.
- **O2 dependencies: ported.** `requests>=2.31.0` in both requirements files. scikit-learn was already declared. The branch's stricter pins are kept.
- **O3 `.gitignore`: ported**, plus `test/fixtures/**/.workflow/*.lock` (commit 81e5480). Without it, the lock sidecar created by fixture writes appears as untracked, because the fixture `.workflow` is re-included by a negation.
- **O4, O5 locking and state: ported (partial).** `workflow/locking.py` verbatim. **Adaptation:** only the writers are locked (`save_state`, `update_*`, `record_step`, `update_background_task`, `record_checkpoint_metadata`). Reads are unlocked because the state file is replaced atomically. A locked read would create `.workflow/state.lock` on read-only browsing, which the accepted FR-035 test forbids (`test_webui.py::test_readonly_browsing_never_writes_into_a_project`). `list_checkpoint_metadata` is kept: the branch already sorts newest-first with an index tie-breaker.
- **O6 artifact cache: ported (partial).** Ported: the file lock around `request`, atomic JSON writes, the content-cache-v2 key (node version, package code fingerprint, Python and dependency versions, content-hashed inputs), `optional_inputs`, missing required input fails, `False` or `status: failed` is a failure, outputs must exist and be non-empty, and output retirement to `.artifact_history`. `content_hash_inputs` is kept as a no-op for call-site compatibility.
  - **Deferred, not ported:** the output-content fingerprint, and the "inputs changed during computation" failure. Both are blocked by one pre-existing branch defect: `_write_comparative_reports` (called by the comparative, polypharmacology and biology-correlation nodes, which run concurrently in one tier) rewrites the top-pose atlas files that `top_pose_atlas` owns. The legacy `run()` path does the same. Main's fix (`preserve_dag_atlas`) changes which atlas is final and sits in S9, which is keep_branch. **Decision needed from the Lead.**
  - **Wiring:** `interactions` and `reports` declare their optional producers as `optional_inputs`, so the Spec 026 decoupling holds.
  - **Test adaptations:** `test/test_dockforge_smoke.py` report node declares `optional_inputs` (main's change). The parameter-invalidation check expects `validation_gate` to be a content-addressed cache hit, because `normalized_scores.csv` is byte-identical across the two methods in that fixture (verified by SHA-256).
  - **Test status:** main's `test_cache_rejects_modified_output` and `test_input_mutation_during_computation_is_not_cached_as_success` are strict xfail (deferred, reason above).
- **O7 latest-session shortcuts: ported.** Symlinks as directory links, no deletion of a real directory at the shortcut path, and a `LATEST_*.txt` pointer when symlinks fail. The branch's `rmtree` is removed.
- **O8 job contract: ported.** `docking/runners/job_contract.py` verbatim (stdlib only). Planning records a protocol fingerprint, completion file, hashed inputs and executables. `run()` executes through `execute_job`. `--skip-completed` reuses a job only with a valid completion record, not merely an existing pose.
- **O9 deployment: ported.** Executor bundling, POSIX path mapping, `_apply_allocated_resources`, job wrapper commands, and `skip_completed` only for local deployments. Spec 036 and 037 parameter records are kept.
- **O10 runner base and Vina collector: ported (partial).** Spec 036 replicate tags, seeds and R5d warnings are kept. `accepted_pose_files` is replicate-aware. Main's `PairlistRow` reference fields (S18) are not ported.
- **O11 gnina flags: ported.** `--no_gpu` and independent `--cpu`. **Not ported:** main's pose-based collector (it replaces the accepted Spec 031 log-based `all_scores` ingestion; S-keep).
- **O12 AutoDock4 staging: ported.** Fixed staged names in the pair directory, relative GPF/DPF names, `.12g` spacing, `parameter_file` staging, accepted-pose collection, and input hashes that include the staged copies.
- **O13 geometry and exit codes: ported (partial).** `validate_pairlist_geometry`, finite box checks, refusal before planning, `dock` returns 1 on any failed job, and the preflight geometry and fail-closed QC errors. **Not ported:** main's `_restore_saved_protocol` (S27; D18 keeps the branch's basic-mode refusal).
- **O14 symlink fallback: ported (partial).** The README is added. The branch's accepted status name `fallback_dir_created` is kept.
- **O15 spreadsheet guards: ported (partial).** Duplicate PDB/property rows are refused; infinite centres become missing; site sizes, when present, are checked positive. The size columns are not added to the branch's catalogue (`center_method` is accepted, Spec 036 R6).
- **O16 asset aliases: ported (partial).** The `_index_assets` and `_resolve_asset` ambiguity refusal, and the preference for `X_cleaned.pdbqt` over raw `X.pdbqt`. Main's `asset_identity.py` module is not ported (parallel implementation; S18).
- **O17 `applymap`: present.** The branch's `_map_frame_elements` already uses `DataFrame.map` with a fallback.
- **O18 score tables: ported (partial).** Runner-wrapper `*.runner.log` files are excluded. **Reverted (5b12017):** the unlisted-log rule, which contradicts the accepted Spec 031 directional-matching contract (an unmatched extended-name log keeps its own tag). **Not ported:** the pose-completion checks (they depend on S4) and the logging change.
- **O19 engine detection: deferred.** Main's gnina detection depends on `pose_geometry` (S4) and `score_import` (D13). Its pairlist filter by stem would drop every replicate pose, because branch names carry Spec 036 `__repNN` suffixes. A replicate-aware version needs a Spec 022 review.
- **O20 ProLIF maps: kept (flagged).** Main's diff changes interaction chemistry: implicit hydrogens, receptor standardisation, a required ligand SDF sidecar, a status file, and new interpretation labels. No D-decision covers it, so the branch version stays. A candidate for a follow-up decision.
- **O21 report wording: ported.**
- **O22 `score_semantics.py`: not ported.** After D16, no consumer on the branch needs it, and the branch's inline direction rules are accepted (Spec 031 values match). Porting it would add an unused module.
- **O23 smoke additions: partial.** Only the contracts needed by ported units are applied (the report node, the parameter invalidation, the Vina collapse). Main's `_mark_score_import` and the score-import fixtures are not ported (D13).
- **O30 README: kept.**
- **O29 `docs/correctness-update.md`: not ported.** It describes main's state, including deferred and rejected items, so it would claim behaviour this branch does not have. This port log is the branch's record.

### M2 scientific (Lead decisions D2 to D17)

- **D2 AutoDock4 even npts: ported.** `_grid_points = 2*ceil(size/(2*spacing))`; 20 Å at 0.375 Å gives 54. Verified by `test/test_spec035_autodock4_grid.py`: even count, box covered, at most one spacing of slack, and non-finite or non-positive inputs refused. **Not verified against the binary:** `autogrid4` is not installed here, so the even-npts convention is taken from the documented AutoGrid requirement.
- **D3 failed pH step refuses: ported.** The labelled RDKit fallback and its function are removed. A pH failure raises `ProtonationStateError` (`ph_model_unavailable`). The output contract treats the `rdkit` backend label as an error.
  - Tests: `test_spec034_protonation`, `test_spec036_protonation_policy`, `test_spec034_lineage` pass unchanged (none exercised the fallback). Main's `test_ph_failure_does_not_fall_back_to_add_hydrogens` is adapted to match `ph_model_unavailable` (the branch's accepted Spec 034 reason code). Main's `test_missing_openbabel_does_not_substitute_untreated_rdkit_state` is adapted to `RuntimeError` with "Open Babel" (the accepted refusal type).
- **D4 refuse PDB; keep 3D input: ported.** PDB is refused in `normalize_ligand_to_sdf` and `prepare_ligand_for_vina_family`. Input with 3D coordinates is not regenerated (no `--gen3d`) and not re-minimised (no MMFF/UFF). The payload records `had_3d_input` and `generated_3d`.
  - **Not ported (no decision covers them):** the one-molecule SDF rule, the sidecar reference-field copy, and the source-SMILES and requested-pH summary fields.
- **D5 no PDB to SDF with inferred bonds: ported.** `run_pdb_collect` no longer converts. An earlier SDF moves to `metadata/previous_ligand_references`, and the user is told to supply a verified SDF. **Not ported:** main's RCSB retrieval in the same block (D14).
- **D6 no atom-type guessing: ported.** `sanitize_prepared_ligand_pdbqt` validates and refuses; the guessing helpers are removed. `repair_project_ligands` uses the raw SDF whenever one exists, and refuses a PDB-only ligand (D4). The branch's REMARK-SMILES repair is kept (no decision). Main's full repair rewrite (sidecars, asset identity) is not ported.
- **D7 pocket "Not evaluated": ported.** Pocket volume and druggability report `None` and "Not evaluated" with method labels. The electrostatic and hydrophobic residue counts keep their labelled status.
  - **Test adapted (accepted Spec 031):** `test/test_spec031_remediation.py::test_pocket_outputs_disclose_heuristic_contract` now asserts the D7 contract.
- **D8 redocking gate: ported.** HIGH, the only route to `validated`, now also needs coverage 1.0, at least 3 unique experimental reference complexes (keyed by PDB id, frame and reference-pose SHA-256), and pass rate at least 0.80. Unmet criteria are named in `validation_gate_reason`. The labels HIGH/MEDIUM/LOW are kept. The gate marks and does not block (Spec 036 R2 is unchanged).
- **D9 top-pose default: already present.** Spec 037 R4b on this branch (7c3f34c) moves cross-engine raw-affinity uses to the v2 consensus rank. The default `best_affinity` is checked to be single-engine.
- **D10 cross-engine raw comparisons: partly present, partly open.** The positive advantage gate already raises (Spec 037 R4b). The simplified pipeline's polypharmacology `not_evaluable` status (S26) is **not** present. Not in the M2 list. **STOP: needs a decision before porting.**
- **D11, D12, D18 confirmations: branch kept.** Kabsch inter-engine agreement (D11). Legacy consensus modes removed (D12). Basic mode refuses explicit `--exhaustiveness` (D18). Main's atlas tests that use legacy modes are xfail (D12).
- **D13 score import: deferred.** Main's `score_import.py` is not ported. Its tests are xfail or absent.
- **D14 RCSB retrieval: not ported.** `test/test_reference_chemistry.py` is skipped at module level, because all four tests need `workflow.reference_chemistry`.
- **D15 QC group: ported (four items).** Missing AutoDock RMSD is NaN (`docking_parser`). No Vina collapse warning (smina's is kept). Receptor lines keep their record name and chain. Biology keys must be non-missing and unique, and the merge is many-to-one. **Test adapted (smoke):** the Vina solo summary assertion now requires the warning's absence.
- **D16 legacy normalisers: ported (partial).** A missing score stays NaN. Degenerate groups keep the branch's 0.0 for finite rows. Not ported: main's 0.5 midpoint and main's grouping and winner rules (no decision covers them).
- **D17 reference frame guard: ported, inert.** `classify_reference_row` clears a reference whose `reference_frame_id` and `receptor_frame_id` are both present and differ. The branch's pair rows do not yet carry those fields (S18), so the guard does nothing until they are produced. Checked by unit calls, not by a fixture run.
- **D19 audit: no import.** The `audit/` directory is absent; main's findings are not imported.

### Silent-merge scientific files (inventory section 0)

| File | Disposition |
|---|---|
| `docking/preparation/ligand_quality.py` (S14) | D6 ported |
| `post_docking_analysis/prolif_interaction_maps.py` (O20) | kept; no decision covers the chemistry (flagged) |
| `post_docking_analysis/biology_integration.py` (S25) | D15 ported |
| `post_docking_analysis/docking_parser.py` (S24) | D15 ported |
| `post_docking_analysis/enhanced_rmsd_analyzer.py` (S5) | kept (accepted Spec 036 R3d; Kabsch diagnostic) |
| `workflow/execution.py` (S15) | D5 ported |
| `post_docking_analysis/simplified_pipeline_impl.py` (S26) | **not ported: D10 open (see STOP)** |

## 3. Test adaptations (every change to an existing accepted test)

| File | Test | Change | Why |
|---|---|---|---|
| `test_spec036_docking_params.py` | `test_poses_returned_and_warning_when_fewer_than_num_modes` | fake poses gain an ATOM record; receptor and ligand staged | a pose without atoms is invalid under the job contract (O8); inputs must be hashable |
| `test_spec036_integration_fixes.py` | `test_replicate_jobs_record_wall_time_and_duration` | same | same |
| `test_spec031_remediation.py` | `test_pocket_outputs_disclose_heuristic_contract` | asserts D7 "Not evaluated" | D7 (Lead) |
| `test_dockforge_smoke.py` | artifact-graph report node | `optional_inputs` | main's contract; O6 |
| `test_dockforge_smoke.py` | parameter invalidation | `validation_gate` is a content-addressed hit when inputs are byte-identical | O6 (verified by SHA-256) |
| `test_dockforge_smoke.py` | Vina solo summary | must not warn on collapse | D15 (Lead) |
| `test_docking_correctness.py` | `test_deploy_preserves_reviewed_local_protocol` | advanced mode; deploy receives explicit values | D18 (basic refuses explicit exhaustiveness); S27 kept |
| `test_docking_correctness.py` | `test_nonfinite_box_schema_rejected` | adds a seed | Spec 036 R5b (seed required) |
| `test_docking_correctness.py` | `test_symlink_privilege_failure_has_usable_fallback` | status `fallback_dir_created`; README asserted | Spec 033 status name kept; O14 README |
| `test_preparation_correctness.py` | pH-failure and missing-OpenBabel tests | expect the accepted `ProtonationStateError` (`ph_model_unavailable`, `RuntimeError`) | D3; Spec 034 |
| `test_preparation_correctness.py` | two receptor tests | pH, `ph_source`, force field and its source supplied | Spec 033/034 accepted inputs |
| `test_preparation_correctness.py` | `asset_identity` import | removed (S18) | parallel module not ported |
| `test_postdocking_correctness.py` | imports of `pose_geometry`, `compute_geometric_consensus`, `_is_better_pose` | guarded | S4, S3, S11 kept |

### Strict or non-strict xfail (contradictions kept on purpose, with reason)

- **Workflow:** O6 integrity (2 tests).
- **Docking:** S18 `PairlistRow` reference fields (1); S-keep gnina log ingestion (1).
- **Preparation:** S21 core_pipeline (4); deferred altloc and engine-aware validation (4); deferred stale-copy refusal (1); S18 pair reference (2); S16 receptor backend (2).
- **Post-docking, strict:** S4 pose geometry (8); S3 geometric consensus (3); D13 score import (2); D12 legacy modes (1 test, 4 parameters); S6 redocking helpers (1 test with 2 parameters plus 1); S9 DAG and per-engine winner (7); O19 source fingerprint (1); O18 unlisted-log rule (1).
- **Post-docking, non-strict** (mixed outcomes; some parameters pass on the branch): S2 normaliser direction (1 test, 3 parameters) and S6 invalid pose index (1 test, 2 parameters).
- **Skipped:** `test_reference_chemistry.py` (D14).

## 4. Full suite and smokes

Baseline at 7c3f34c: 475 passed, 2 skipped (before `pytest.ini` collected the smoke contracts).
Final at a66ca7c and after: see section 5.

## 5. Open points for the Lead (STOP)

1. **D10 polypharmacology (S26).** Not in the M2 list, and the branch lacks the `not_evaluable` cross-engine status. Decide before porting.
2. **DAG write ownership (O6).** Three report nodes rewrite `top_pose_atlas` outputs. Blocks the output fingerprint, the changed-input failure, and main's `preserve_dag_atlas`. Needs an S9 decision.
3. **S18 pair reference fields.** Blocks D17 from taking effect and four deferred tests. Needs the PairlistRow reference-field decision.
4. **O20 ProLIF chemistry and O19 detection.** Candidate follow-ups; neither is covered by a D-decision.
