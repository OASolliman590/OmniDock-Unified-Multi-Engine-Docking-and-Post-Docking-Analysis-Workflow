# Tasks: Post-Docking Scientific Foundation

**Input**: `/specs/020-post-docking-scientific-foundation/spec.md`, `/specs/020-post-docking-scientific-foundation/plan.md`  
**Scope**: scientific correctness, statistical rigor, and reproducible consensus behavior

## Format: `[ID] [P?] Description`

- **[P]**: parallelizable

## Phase 1: Pipeline Unification

- [x] T001 Create `post_docking_analysis/unified_pipeline.py` entrypoint.
- [x] T002 Route workflow execution through unified post-docking orchestration.
- [x] T002a Route `post_docking_analysis.simplified_cli` to unified execution for manifest-backed DockForge projects (legacy fallback kept only for unsupported flag combinations).
- [x] T003 Retire legacy logic from `simplified_pipeline.py` and `multi_engine_pipeline.py` (keep wrappers only).
- [x] T004 Add explicit smoke assertion that all user-facing post-docking commands resolve to unified path.

## Phase 2: Critical Correctness Fixes

- [x] T005 Fix normalization direction in `consensus.py` (`per_engine_minmax`, `per_engine_zscore`).
- [x] T006 Keep rank-based normalization contract with strongest binder mapped highest.
- [x] T007 Replace ambiguous pairlist matching with strict deterministic mapping in `generate_scores_csv.py`.
- [x] T008 Emit unmatched and empty/malformed log reports (`unmatched_log_files.txt`, `empty_log_files.txt`).
- [x] T009 Enforce deterministic pose tie-break logic in `pose_extractor.py`.
- [x] T010 Keep dashboard/report newline contract fixed in `report_generator.py`.

## Phase 3: Validation Gate

- [x] T011 Emit `validation_gate_status.json` from redocking validation stage.
- [x] T012 Enforce effective policy fallback when `reference_anchor` is not allowed by validation gate.
- [x] T013 Mirror gate status into canonical consensus outputs.
- [x] T014 Include requested/effective hit-class policy and gate state in consolidated run context.
- [x] T015 Add smoke assertions for gate artifact presence and fallback behavior.

## Phase 4: Scoring and Classification Correctness

- [x] T016 Prevent single-engine low-agreement warning from triggering in single-engine runs.
- [x] T017 Keep atomic QC downgrade behavior (single final class decision).
- [x] T018 Enforce N<3 classification suppression (`classification_basis=insufficient_n`).
- [x] T019 Set RMSD default scopes to `per_complex` in CLI/interactive/pipeline defaults.
- [x] T020 Remove dead RMSD scope branches (`per_protein`, `global`) from user-facing execution contracts.

## Phase 5: Statistical Rigor

- [x] T021 Add minimum sample-size guard (`MIN_CORRELATION_N`) for correlation outputs.
- [x] T022 Add BH/FDR correction and q-values to correlation outputs.
- [x] T023 Use rank-normalized values for cross-engine correlation analysis.
- [x] T024 Add contract smoke for per-protein/global correlation outputs under sparse and mixed sample sizes.

## Phase 6: Structure Integrity and Complex Outputs

- [x] T025 Retire `structure_processor.py` stub and replace with active structure validation module.
- [x] T026 Preserve ligand residue identity in complex generation path (avoid forced `UNK`).
- [x] T027 Keep PDBQT fallback parser chain/residue contract verified by smoke checks.
- [x] T028 Add regression fixture with multi-ligand complex to ensure residue identity remains stable end-to-end.

## Phase 7: DockBox Consensus and Explainability

- [x] T029 Keep `dockbox_geometric` as default consensus mode with explainability output.
- [x] T030 Preserve weighted/strict/favorite modes as explicit alternatives.
- [x] T031 Add synthetic regression fixture proving geometric consensus and rank aggregation remain mode-isolated.

## Documentation and Release Readiness

- [x] T032 Add `specs/020-post-docking-scientific-foundation/plan.md`.
- [x] T033 Add `specs/020-post-docking-scientific-foundation/tasks.md`.
- [x] T034 Update user-facing post-docking guide with validation-gate behavior, effective policy fallback, and new artifacts.
- [x] T035 Publish migration note for legacy-to-unified post-docking path behavior.

## Verification Tasks

- [x] V001 `python -m compileall` on touched workflow/post-docking modules.
- [x] V002 `python test_pipeline.py`.
- [x] V003 `python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix`.
- [x] V004 Full `python test/test_dockforge_smoke.py` when optional tool dependencies are available.

## Exit Criteria

- [x] E001 All P0/P1 scientific defects in spec are closed with passing smoke checks.
- [x] E002 Unified pipeline is the only user-facing post-docking execution path.
- [x] E003 Validation gate + effective policy fallback are visible in reports and canonical mirrors.
- [x] E004 Correlation outputs are minimum-N safe and multiple-testing corrected.
- [x] E005 Documentation reflects final behavior and migration path.
