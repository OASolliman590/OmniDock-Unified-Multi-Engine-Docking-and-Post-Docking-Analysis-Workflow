# Spec 031 — Appraisal Remediation Task Ledger

**Status:** Gate 1 and Gate 2 implementation scopes completed and accepted by the Scientific Lead on 2026-08-21  
**Source of truth:** `spec.md` and `plan.md` in this directory  
**Task count:** 32  
**Execution rule:** Do not start a task whose state is `blocked_gate_2` or whose dependencies are incomplete.

## Implementation Snapshot — 2026-08-21

- **Completed Gate 1 scope:** T001–T011, T020–T025, and the applicable engineering portions of T027–T031.
- **Evidence still unavailable:** T026 retains `evidence_pending` for the sanitized NMRbox project, commands, engine identities, score inputs, pose checksums, and truth labels.
- **Gate 2 implemented:** T012–T019 are complete under the dated records in `decisions/`; NetworkX `>=3.2,<4` is the approved bounded atom-mapping backend.
- **Human review:** T032 is complete under `decisions/gate-2-implementation-acceptance.md`.

## State Definitions

- `ready_after_gate_1`: engineering work may begin after the revised specification is approved.
- `ready_after_gate_2`: scientific-method implementation may begin because its dated decision record is approved.
- `blocked_gate_2`: implementation depends on a Scientific Lead policy decision recorded at Gate 2.
- `evidence_pending`: the reported observation must be reproduced or imported before a fix is selected.
- `closeout`: verification or human review performed after its implementation dependencies complete.
- `done`: task and stated verification have completed with retained evidence.

## Phase 0 — Evidence Reconciliation and Reproduction

### T001 — Freeze the evidence and fixture inventory

- **State:** `ready_after_gate_1`
- **Covers:** Evidence classification for all 18 findings: F-01 through F-07, NMR-01 through NMR-09, and M-01/M-02.
- **Work:** Add a sanitized evidence manifest for every finding, recording `verified_code`, `observed_artifact`, `reported_unverified`, or `superseded_claim`. Record the assessed commit, current commit, input/output paths or absence reason, SHA-256 values, command, tool versions, and status vocabulary.
- **Files:** `specs/031-appraisal-remediation-p0/evidence/README.md`, `specs/031-appraisal-remediation-p0/evidence/manifest.json`.
- **Dependencies:** Gate 1 approval.
- **Verify:** JSON schema/field assertions in `test/test_spec031_remediation.py`; no private or credential-bearing paths in committed evidence.

### T002 — Create the discoverable remediation regression suite

- **State:** `ready_after_gate_1`
- **Covers:** FR-003 and the regression harness used by all approved requirements.
- **Work:** Create `test/test_spec031_remediation.py` with fixture helpers and initially failing or skipped-by-explicit-reason tests mapped to requirements. Tests must use the repository status vocabulary; a missing artifact must never pass.
- **Files:** `test/test_spec031_remediation.py`.
- **Dependencies:** T001.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py --collect-only`.

### T003 — Reproduce and classify the Vina zero-row failure

- **State:** `evidence_pending`
- **Covers:** FR-001; F-01/F-02 and NMR-04.
- **Work:** Run the same Vina artifact through `parse_vina_pdbqt()` directly and through the canonical/unified loader. Test canonical, flat, legacy, and mixed per-artifact layouts plus pairlist-tag/pose-stem mismatches. Preserve diagnostics and classify the failing layer before editing production code. Do not encode the superseded GNINA-parser explanation.
- **Files inspected:** `docking/runners/vina.py`, `post_docking_analysis/unified_analysis.py`, `post_docking_analysis/engine_hpc_adapter.py`, `docking/project_layout.py`.
- **Dependencies:** T001, T002.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "vina and reproduction"` plus the evidence-manifest entry for NMR-02.

### T004 — Repair the smoke harness call contract

- **State:** `ready_after_gate_1`
- **Covers:** FR-003; F-03.
- **Work:** Make `_prepare_non_gnina_context()` and all call sites agree on `analysis_config`, preserving backward compatibility where appropriate. Add a regression covering the failing call path.
- **Files:** `test/test_dockforge_smoke.py` and the production helper only if the reproduced contract requires it.
- **Dependencies:** T002.
- **Verify:** `python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix`.

## Phase A — Evidence-Backed P0/P1 Engineering Corrections

### T005 — Fix the proven ingestion or artifact-resolution defect

- **State:** `evidence_pending`
- **Covers:** FR-001, OR-003; F-01/F-02 and NMR-04.
- **Work:** Apply the smallest correction at the layer demonstrated by T003. If the defect is layout resolution, reuse or extend the read-only engine layout detector and resolve each artifact independently; do not call a mutating migration helper while reading. If T003 identifies another layer, update this task record and evidence before implementation.
- **Files:** Limited to the production files identified by T003 and focused tests.
- **Dependencies:** T003.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "vina and (ingestion or layout)"`.

### T006 — Fail explicitly on empty or partially missing ingestion

- **State:** `ready_after_gate_1`
- **Covers:** FR-002; F-02.
- **Work:** Replace silent zero-row success with an explicit status and actionable diagnostics containing engine, expected artifact type, resolved locations, ambiguity, parser outcome, and missing-dependency/configuration distinction. Preserve successful partial-engine disclosure where policy allows it.
- **Files:** Unified analysis/orchestration status path selected in T003; focused tests.
- **Dependencies:** T003; coordinate with T005.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "empty_ingestion or ingestion_status"`.

### T007 — Correct pandas compatibility and exercise real plotting

- **State:** `ready_after_gate_1`
- **Covers:** FR-004; F-06.
- **Work:** Replace the incompatible DataFrame operation with the supported equivalent and add a test that constructs the affected frame and executes the plotting path. A mere import test is insufficient.
- **Files:** Plotting/analysis module named by the NMR-01 evidence entry; focused tests.
- **Dependencies:** T001, T002.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "pandas or plot"`.

### T008 — Normalize manifest paths to POSIX form

- **State:** `ready_after_gate_1`
- **Covers:** FR-005; F-07.
- **Work:** Serialize manifest-relative paths with `/` on every platform while preserving correct native filesystem resolution at runtime. Cover Windows-shaped and POSIX-shaped inputs without changing scientific identifiers.
- **Files:** Manifest serialization/loading modules identified by evidence; focused tests.
- **Dependencies:** T001, T002.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "manifest and path"`.

### T009 — Share model-aware pose extraction and honor selected pose

- **State:** `ready_after_gate_1`
- **Covers:** FR-006; NMR-06.
- **Work:** Introduce or consolidate a model-aware PDBQT pose extractor used by both unified analysis and redocking. Reject an out-of-range selected pose explicitly; never silently fall back to model 1. Preserve pose/model provenance.
- **Files:** `post_docking_analysis/unified_analysis.py`, redocking module identified by evidence, focused tests.
- **Dependencies:** T002.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "selected_pose or model_boundary or redock"`.

### T010 — Make co-crystal metadata propagation atomic

- **State:** `ready_after_gate_1`
- **Covers:** FR-007, OR-005; NMR-05 and M-01.
- **Work:** Normalize booleans so `False`, `0`, empty, and missing remain false, preserve required metadata from `pair_intent.csv` into the relevant canonical contract, and avoid partial updates. Do not invent source lineage from file naming. Preserve current legacy classification behavior until SR-004 is approved; this task supplies the metadata needed to change that behavior safely.
- **Files:** Pairlist builder/validator, project materialization path, schemas, focused tests.
- **Dependencies:** T001, T002.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "reference_ligand or cocrystal or boolean"`.

### T011 — Run the offline 1IEP Phase A regression

- **State:** `closeout`
- **Covers:** FR-001 through FR-007 and the offline portion of OR-004.
- **Work:** Exercise the sanitized 1IEP fixture without network access or docking-engine execution. Record receptor `1IEP` chain `A`, ligand `STI` residue `A:201`, Vina `1.2.7` Python mode with `sf=vina`, legacy center `[16.757738, 56.582794, 13.778150]` Å, box `20 × 20 × 20` Å, exhaustiveness `8`, seed `12345`, CPU `2`, and requested poses `5`. Mark pH, tautomer, charge, hydrogen, and exact Open Babel version as unknown unless supplied by evidence.
- **Files:** Sanitized fixture and `test/test_spec031_remediation.py`.
- **Dependencies:** T005 through T010.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "1iep and phase_a"`.

## Phase B — Scientific Method Corrections

### T012 — Approve the atom-mapping and RMSD comparability policy

- **State:** `done`
- **Covers:** SR-001; NMR-09.
- **Decision:** Approved `spec031-atom-mapping-v1`: heavy atoms only, symmetry-aware bond-graph mapping, complete coverage, minimum Kabsch RMSD, and `not_comparable` for every unproven mapping. Centroid-based correspondence is prohibited.
- **Artifact:** `specs/031-appraisal-remediation-p0/decisions/atom-mapping.md`.
- **Dependencies:** Gate 2 approval recorded 2026-08-21; T009 supplies pose semantics.
- **Verify:** Decision record includes RMSD definition, atom inclusion, symmetry, alternate location, protonation/topology mismatch, and failure policy.

### T013 — Implement approved symmetry-aware RMSD comparability

- **State:** `done`
- **Covers:** SR-001; NMR-09.
- **Work:** Implement only the policy approved in T012. Emit mapping method, mapped/total atom counts, ambiguity, and `not_comparable` reason. Never report a numeric RMSD for an invalid mapping.
- **Files:** RMSD/geometry module selected after T012; schema and focused tests.
- **Dependencies:** T012.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "rmsd or atom_mapping or symmetry"`.

### T014 — Approve the consensus-ranking policy and version boundary

- **State:** `done`
- **Covers:** SR-002; NMR-09 and M-02.
- **Decision:** Approved `consensus_rank_geometry_qc_v2`: equal-weight complete-case within-engine/within-target rank percentiles; geometry is comparable-only QC/tie-breaker; missing components are never imputed or renormalized.
- **Artifact:** `specs/031-appraisal-remediation-p0/decisions/consensus-v2.md`.
- **Dependencies:** Gate 2 approval recorded 2026-08-21; T012 because geometry comparability affects eligibility.
- **Verify:** Decision record explicitly forbids raw cross-engine score equivalence and cross-target acceptance tests.

### T015 — Implement versioned consensus semantics

- **State:** `done`
- **Covers:** SR-002; NMR-09 and M-02.
- **Work:** Implement the approved policy as an explicitly versioned method, fix score direction consistently in top-pose selection, disclose included/failed/missing engines, and preserve or migrate legacy output only as approved. Do not use the unrelated `7MYL` and `4ZBE` targets as an equality or exact-rank acceptance criterion.
- **Files:** Consensus/top-pose modules, output schema, report templates, focused tests.
- **Dependencies:** T013, T014.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "consensus or top_pose or score_direction"`.

### T016 — Approve grid-center policy and provenance

- **State:** `done`
- **Covers:** SR-003; NMR-07.
- **Decision:** Approved `binding_site_center_v1`: selected co-crystal ligand heavy-atom centroid; apo/no-reference cases require explicit user coordinates or return `skipped_missing_configuration`; no cavity or PLIP contact-atom fallback.
- **Artifact:** `specs/031-appraisal-remediation-p0/decisions/grid-center.md`.
- **Dependencies:** Gate 2 approval recorded 2026-08-21; T010 for authoritative reference metadata.
- **Verify:** Decision records coordinate system, units, source selection, atom selection, fallback/skip behavior, and legacy alias policy.

### T017 — Implement approved grid-center provenance

- **State:** `done`
- **Covers:** SR-003; NMR-07.
- **Work:** Implement the approved strategy and persist center coordinates, units, source method, source ligand/residue, selected atom count, input checksum, tool/version where applicable, and legacy aliases. Migrate callers without silently changing existing projects.
- **Files:** Grid-generation/materialization modules, schemas, reports, focused tests.
- **Dependencies:** T010, T016.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "grid_center or center_provenance"`.

### T018 — Approve the reference-classification migration policy

- **State:** `done`
- **Covers:** SR-004; NMR-05 and M-01.
- **Decision:** Approved `explicit_reference_metadata_v2`: explicit metadata is authoritative; name heuristics are disabled by default; opt-in legacy mode creates warned migration candidates but cannot establish a current validation anchor before confirmation.
- **Artifact:** `specs/031-appraisal-remediation-p0/decisions/reference-classification.md`.
- **Dependencies:** Gate 2 approval recorded 2026-08-21; T010 established the metadata propagation baseline.
- **Verify:** Decision record defines authority, legacy-mode scope, missing-lineage behavior, warnings, method/version, and historical reproducibility.

### T019 — Implement the approved reference-classification policy

- **State:** `done`
- **Covers:** SR-004; NMR-05 and M-01.
- **Work:** Consolidate reference/source-layer inference behind the approved policy. Explicit metadata is authoritative; missing lineage is `unclassified`/non-reference; any approved legacy migration is anchored, warned, and versioned. Verify that `novel_ligand_1`, `"False"`, empty, and NaN never become references accidentally.
- **Files:** Reference inference sites in `redocking_validation.py`, `simplified_pipeline_impl.py`, and the shared policy module selected after T018; schemas/reports and focused tests.
- **Dependencies:** T010, T018.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "reference_classification or legacy_inference"`.

## Phase C — Operational Hardening and Communication

### T020 — Add exact pairlist deployment

- **State:** `ready_after_gate_1`
- **Covers:** OR-001; NMR-01.
- **Work:** Add `--pairlist-file` to the deploy interface at `docking/cli.py::deploy_main`, `workflow/cli.py::_add_deploy_args`, and `workflow/cli.py::_build_deploy_argv`. Make it mutually exclusive with rerun-manifest selection, preserve canonical pairlist/state, record the selected file and SHA-256, and generate only the exact requested jobs. Do not submit jobs.
- **Files:** `docking/cli.py`, `workflow/cli.py`, deploy tests and schema/provenance code.
- **Dependencies:** T001, T002.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "exact_pairlist or deploy"`.

### T021 — Resolve mixed engine layouts read-only and per artifact

- **State:** `ready_after_gate_1`
- **Covers:** OR-003; NMR-04.
- **Work:** Extend `post_docking_analysis/engine_hpc_adapter.py::detect_engine_layout()` or a shared read-only resolver so scores, poses, logs, and metadata can each resolve independently across canonical/flat/legacy layouts. Report ambiguity; do not mutate or migrate layouts during analysis reads.
- **Files:** `post_docking_analysis/engine_hpc_adapter.py`, unified analysis loader, focused tests.
- **Dependencies:** T003; coordinate with T005 to avoid duplicate resolution logic.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "mixed_layout or ambiguity or artifact_resolution"`.

### T022 — Invalidate every cache layer and propagate force

- **State:** `ready_after_gate_1`
- **Covers:** OR-002; NMR-03.
- **Work:** Define cache keys and invalidation for raw GNINA `all_scores.csv`, normalized analysis cache, DAG raw-input dependencies, and layout-detection cache. Propagate force/recompute consistently and record why a cache entry was accepted or invalidated.
- **Files:** GNINA runner/aggregation, unified analysis/cache, DAG dependency declarations, layout detector, focused tests.
- **Dependencies:** T021.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "cache or force or resume"`.

### T023 — Align environment dependency manifests

- **State:** `ready_after_gate_1`
- **Covers:** OR-004; F-04.
- **Work:** Reconcile runtime and development dependency declarations with the imports and test/plot paths actually exercised by Spec 031. Add no dependency without a demonstrated requirement; pin or bound versions according to repository policy.
- **Files:** Existing environment/requirements manifests only; dependency contract test if present.
- **Dependencies:** T007 and import evidence from T002.
- **Verify:** Repository-supported environment validation plus `python -m pytest -q test/test_spec031_remediation.py -k "dependency_manifest"`.

### T024 — Repair the installed entry-point contract

- **State:** `ready_after_gate_1`
- **Covers:** OR-004; F-05.
- **Work:** Point the package entry point at an importable callable with preserved CLI behavior and add a test that resolves and invokes its help/version path without executing docking.
- **Files:** Package metadata and the owning CLI module; focused tests.
- **Dependencies:** T002.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py -k "entry_point"`.

### T025 — Correct scientific and operational communication

- **State:** `ready_after_gate_1`
- **Covers:** Section 10; NMR-02, NMR-08, M-02, and linked scientific-communication findings.
- **Work:** Label fixed-sphere volumes, empirical residue counts, and charge ratios as heuristic descriptors; document ML classification basis; state engine-score incommensurability and consensus limitations; align `submit_mode` and CLI documentation with executable behavior. Do not claim affinity, selectivity, mechanism, or clinical efficacy.
- **Files:** Report templates, CLI help, schema descriptions, and relevant user documentation.
- **Dependencies:** Gate 1 for independent wording; T014/T015 only before finalizing consensus wording.
- **Verify:** Documentation assertions/snapshots in `test/test_spec031_remediation.py` and `python -m pytest -q test/test_webui.py`.

### T026 — Retain or explicitly classify the NMR evidence set

- **State:** `evidence_pending`
- **Covers:** FR-001 evidence integrity; NMR-01 through NMR-09.
- **Work:** Import sanitized NMR reproduction artifacts into the Spec 031 evidence directory. If an artifact cannot be retained, record `reported_unverified` with owner, reason, and the exact reproduction needed; never convert absence into success.
- **Files:** `specs/031-appraisal-remediation-p0/evidence/`.
- **Dependencies:** T003 through T010 as applicable.
- **Verify:** Evidence-manifest completeness assertions in `test/test_spec031_remediation.py`.

## Closeout — Spec Kit Gates 3 through 11

### T027 — Run syntax, formatting, and diff checks

- **State:** `done`
- **Covers:** Gate 5.
- **Work:** Run the repository-supported formatter/linter/type checks if configured. Treat `compileall` as syntax/import-bytecode verification only, not as lint or type checking.
- **Dependencies:** All implementation tasks in the approved branch scope.
- **Verify:** `python -m compileall -q .`; `git diff --check`; configured lint/type commands recorded separately.

### T028 — Run focused, web UI, smoke, and full tests

- **State:** `done`
- **Covers:** Gates 6 and 7.
- **Dependencies:** T027.
- **Verify:** `python -m pytest -q test/test_spec031_remediation.py`; `python -m pytest -q test/test_webui.py`; `python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix`; then `python -m pytest -q`.

### T029 — Validate schemas and canonical output topology

- **State:** `done`
- **Covers:** Gate 8.
- **Work:** Validate changed inputs/outputs and confirm preservation of `0-Input/` through `7-Reports/`, `.meta/`, and `sessions/`, plus `.workflow/state.json` transitions.
- **Dependencies:** T028.
- **Verify:** Project/schema validation commands and fixture assertions retained in the evidence manifest.

### T030 — Validate workflow lint, resume, cache, and force behavior

- **State:** `done`
- **Covers:** Gate 9.
- **Work:** Prove unchanged inputs reuse valid cache entries, changed inputs invalidate all dependent layers, `force` recomputes the requested scope, and mixed-layout detection does not mutate data.
- **Dependencies:** T021, T022, T029.
- **Verify:** Focused cache/resume tests plus the repository workflow-lint command recorded in evidence.

### T031 — Validate provenance and run-manifest integrity

- **State:** `done`
- **Covers:** Gate 10.
- **Work:** Confirm SHA-256 inputs, tool/engine versions, executable paths, container digest when applicable, seed, box definition, parameters, selected pairlist, method versions, status, and generated artifacts in `run_tracking/manifest.json` and `.meta/`.
- **Dependencies:** T011, T015, T017, T019, T020, T030 as included in the approved scope.
- **Verify:** Manifest/schema assertions and checksum recomputation retained with the fixture result.

### T032 — Obtain Human Scientific Lead review

- **State:** `done`
- **Covers:** Gate 11 and final acceptance.
- **Work:** Present completed evidence, deviations, unresolved unknowns, scientific-method versions, failed/missing engines, and limitations. The Scientific Lead accepts, requests revision, or rejects. No external submission, upload, HPC launch, cloud action, or scientific interpretation is authorized by this task.
- **Dependencies:** T031 and all approved Gate 2 decision records.
- **Verify:** Dated approval/revision record in `specs/031-appraisal-remediation-p0/decisions/`.
- **Acceptance:** Approved 2026-08-21; see `decisions/gate-2-implementation-acceptance.md`.

## Dependency Summary

```text
Gate 1 → T001 → T002 ─┬→ T003 → T005/T006/T021 → T022
                      ├→ T004
                      ├→ T007/T008/T009/T010 → T011
                      ├→ T020/T023/T024
                      └→ evidence completion T026

Gate 2 → T012 → T013 ─┐
         T014 ─────────┴→ T015
         T016 + T010 ───→ T017
         T018 + T010 ───→ T019

Approved implementation scope → T027 → T028 → T029 → T030 → T031 → T032
```

## Scientific Assumptions and Non-Goals

- The retained 1IEP parameters describe a regression fixture, not a newly selected biological protocol.
- Protonation pH, tautomer states, charge model, hydrogen policy, and exact preparation-tool versions remain unknown until the Scientific Lead or preserved provenance supplies them.
- No task compares raw scores across engines or targets, infers experimental affinity, or treats docking as proof of mechanism or efficacy.
- No task launches docking engines, submits HPC jobs, uploads data, creates cloud resources, or changes remote state.
