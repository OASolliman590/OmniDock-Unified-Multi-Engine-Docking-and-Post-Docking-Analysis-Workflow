# Feature Specification: Appraisal Remediation and Multi-Engine Pipeline Hardening

**Feature ID**: `031-appraisal-remediation-p0`  
**Created**: 2026-08-19  
**Revised**: 2026-08-21 after Gate 2 scientific approval  
**Status**: Gate 1 engineering and Gate 2 scientific-method scopes implemented and accepted by the Scientific Lead on 2026-08-21  
**Assessed baseline**: `cfe3305` (`main`)  
**Current audit baseline**: `9d7400c` on `chore/install-agentic-bioinformatics-stack`

## 1. Decision and Gate State

Gate 1 was approved and its engineering-safe scope was implemented and verified on 2026-08-19. On 2026-08-21 the Scientific Lead approved the recommended Gate 2 baseline, including the bounded chemistry-dependency policy. The normative method records are `decisions/atom-mapping.md`, `decisions/consensus-v2.md`, `decisions/grid-center.md`, and `decisions/reference-classification.md`. Phase B was then implemented and locally verified under `phase_b_plan.md`, and the Scientific Lead accepted the implementation in `decisions/gate-2-implementation-acceptance.md`.

No HPC submission, cloud resource, remote upload, live engine run, or other remote mutation is authorized by this specification. Verification uses retained or synthetic local fixtures unless the scientific lead separately approves external execution.

## 2. Evidence Integrity

Evidence is classified as follows:

- `verified_code`: directly confirmed in current source and, where relevant, commit `cfe3305`.
- `observed_artifact`: supported by a retained log, manifest, or fixture, but the root cause may remain unresolved.
- `reported_unverified`: reported from the NMRbox session but not reproducible from artifacts currently retained in this repository.
- `superseded_claim`: contradicted by the assessed source and prohibited as an implementation premise.

### 2.1 Critical correction to the original plan

The retained appraisal proves that the canonical Vina analysis produced zero score rows and later failed in `classified_hits`. It does **not** prove that Vina was routed through the GNINA parser:

- `docking/runners/vina.py::VinaRunner.collect_normalized_scores()` calls `parse_vina_pdbqt()` in both current source and `cfe3305`.
- `MultiEngineAnalysisPipeline._build_normalized_frame()` also calls `parse_vina_pdbqt()` for Vina/Smina PDBQT files.
- `SimplifiedPostDockingPipeline` is GNINA-specific and has no engine field.
- `parse_smina_pdbqt()` and `parse_ad4_dlg()` do not exist; the supported functions are `parse_vina_pdbqt()` and `parse_autodock4_dlg()`.

Therefore, “Vina routes through the GNINA parser” is a `superseded_claim`. The valid finding is: **canonical ingestion produced zero rows; exact cause pending fixture-backed reproduction**. Layout resolution and pose-stem/pairlist-tag mismatch are plausible hypotheses, not accepted facts.

### 2.2 Evidence still missing

The repository does not retain the NMRbox 21-complex project tree, exact commands, engine versions, profile, normalized inputs, pose checksums, 4ZBE/6LL5 fixtures, or expected outputs. NMRbox empirical impact statements remain `reported_unverified` until a sanitized evidence package is imported. Missing evidence cannot count as a passing test.

## 3. Scope Inventory

This cycle tracks **18 uniquely identified items**: 9 P0, 6 P1, and 3 P2.

| ID | Priority | Evidence | Corrected scope |
|---|---|---|---|
| F-01 | P0 | observed_artifact | Reproduce and classify canonical Vina zero-row ingestion |
| F-02 | P0 | observed_artifact / verified_code | Fail explicitly at ingestion when requested outputs yield zero valid rows |
| F-03 | P0 | verified_code | Repair smoke API drift and add pytest-discoverable remediation tests |
| F-04 | P1 | verified_code | Complete or clearly separate runtime/test dependency declarations |
| F-05 | P1 | verified_code | Repair or remove the broken `pdb-prepare-wizard` entry point |
| F-06 | P0 | verified_code | Replace pandas `DataFrame.applymap` with a supported compatibility path and execute the plot in tests |
| F-07 | P0 | verified_code | Normalize contract-defined relative manifest paths to POSIX form; source is the Windows follow-up audit, not the original appraisal |
| NMR-01 | P1 | verified_code / reported_unverified impact | Add an exact deploy input at the deploy layer without mutating canonical curation inputs |
| NMR-02 | P2 | reported_unverified | Document and emit deterministic submission-mode selection |
| NMR-03 | P1 | verified_code | Invalidate GNINA intermediate, normalized-score, and DAG raw-source caches correctly |
| NMR-04 | P1 | verified_code / reported_unverified impact | Resolve poses, logs, and scores per artifact with ambiguity diagnostics |
| NMR-05 | P1 | verified_code corrected | Preserve and propagate existing co-crystal metadata; add explicit lineage where available |
| NMR-06 | P0 | verified_code | Select the requested PDBQT model and pass the selected pose through redocking RMSD |
| NMR-07 | P0 | verified_code | Replace mislabeled residue-centroid box provenance after scientific approval |
| NMR-08 | P2 | verified_code | Reconcile CLI aliases and improve discoverability without breaking existing flags |
| NMR-09 | P0 | verified_code | Replace unproven atom correspondence and version the consensus policy after scientific approval |
| M-01 | P0 | verified_code | Remove broad `_ligand_` reference inference atomically with metadata migration |
| M-02 | P2 | verified_code corrected | Decide whether GNINA CNN values participate in multi-engine consensus; single-GNINA paths already use CNN affinity |

The original appraisal’s heuristic-descriptor labeling, provenance strengthening, classification disclosure, and target-specific benchmark recommendations remain linked requirements in Sections 8 and 9; they are not silently discarded.

## 4. Goals

1. Make every reported correctness failure reproducible or explicitly evidence-pending.
2. Repair implementation-safe correctness defects without changing scientific policy.
3. Prevent parsing or mapping failures from masquerading as scientific disagreement.
4. Version all changes that alter ranking semantics.
5. Preserve canonical output topology, workflow state, schemas, cache behavior, and provenance.
6. Produce pytest-discoverable regression coverage plus a retained minimal bioinformatics fixture.
7. Keep scientific decisions explicit and reviewable.

## 5. Non-Goals

- Do not infer affinity, potency, selectivity, mechanism, or clinical efficacy from docking output.
- Do not compare raw scores across engines or targets as if they were a common physical scale.
- Do not silently choose consensus weights, RMSD mapping policy, pH, tautomer, charge model, or docking-box center.
- Do not make the GNINA-specific simplified pipeline generically engine-aware unless reproduction proves that this is the failing layer.
- Do not integrate fpocket/P2Rank, add a new workflow engine, or perform a broad pipeline refactor in this cycle.
- Do not treat one desired rank inversion as validation of a consensus method.
- Do not submit or sync NMRbox/HPC jobs during local verification.

## 6. Scientific Fixture and Assumption Record

### 6.1 Retained 1IEP fixture

| Parameter | Recorded value |
|---|---|
| Receptor | RCSB PDB `1IEP`, C-ABL kinase, chain A |
| Ligand | STI-571; residue `STI A 201` |
| Coordinate units | Cartesian Ångström coordinates from the retained PDB/PDBQT files |
| Engine | AutoDock Vina Python interface 1.2.7, scoring function `vina` |
| Old box center | `[16.757738, 56.582794, 13.778150]` Å; label as the historical PLIP full-residue-derived center |
| Box dimensions | `20 × 20 × 20` Å |
| Exhaustiveness / seed / CPU | `8` / `12345` / `2` |
| Requested poses | `5`; the retained direct parser evidence currently records two parsed score rows and must not claim five parsed rows |
| Protonation pH | Unknown — must not be inferred |
| Tautomer state | Unknown — must not be inferred |
| Charge and hydrogen policy | Unknown — must not be inferred |
| Conversion tool/version/command | Open Babel was used; exact version and complete command are not retained |

The fixture manifest must record full SHA-256 values for the PDB, receptor PDBQT, ligand PDBQT, and Vina pose PDBQT. Unknown preparation parameters remain explicit `unknown` values.

### 6.2 NMRbox fixture

Before NMRbox results can become pass/fail acceptance evidence, import a sanitized manifest containing the exact 21-row normalized inputs, pose references/checksums, engine versions, profile, commands, box definitions, seeds, scoring modes, formula configuration, and expected outputs. Until then, those cases may generate a comparison report only.

## 7. P0 Engineering Requirements

### FR-001 — Evidence-backed Vina ingestion diagnosis

Reproduce the zero-row symptom from a clean assessed baseline using an exact project-tree fixture. Capture filenames, pairlist tags, layout profile, detected artifact paths, direct parser row count, unified-loader row count, command, commit, and environment.

Acceptance:

1. A correctly tagged Vina multi-model PDBQT produces the expected normalized rows through the canonical loader.
2. A tag mismatch and each supported layout variant produce either valid rows or an explicit, actionable diagnostic.
3. The implementation patch targets only the layer demonstrated to fail.

### FR-002 — Explicit empty-ingestion failure

If requested engine outputs exist but yield zero valid rows, the raw-score stage must be `failed` with engine, parser, resolved paths, files examined, and reasons. It must not write a headerless CSV, mark the stage completed, or defer failure to pandas.

### FR-003 — Discoverable test baseline

Repair the `_prepare_non_gnina_context()` smoke signature drift and add `test/test_spec031_remediation.py`. The focused suite must be discoverable by `pytest` and must not rely only on the monolithic smoke script.

### FR-004 — pandas compatibility

Execute `_plot_hit_class_matrix()` on supported pandas versions, assert that the figure exists and is non-empty, and ensure a failed required figure cannot be reported as a successful scientific artifact.

### FR-005 — POSIX relative-path contract

Contract-defined relative paths in run-tracking indexes and `.meta` manifests must use `/` on every OS. Absolute external paths require an explicitly documented serialization policy and are not made portable by relabeling them as relative.

### FR-006 — Shared deterministic pose extraction

Use one PDB/PDBQT/SDF pose-extraction contract with 1-based `pose_index`. Redocking RMSD must pass the selected row’s pose instead of hardcoding pose 1.

Required behavior:

- No-MODEL files support pose 1 only.
- An EOF may close a final model missing `ENDMDL`.
- Out-of-range indices return `pose_index_out_of_range`; they never silently select model 1.
- Malformed nested/incomplete model structure returns an explicit parse failure.
- Outputs record requested pose, selected model/ordinal, format, heavy-atom count, parser version, and error.

### FR-007 — Atomic reference-metadata propagation

Preserve and propagate existing co-crystal metadata before changing any scientific reference-classification policy.

Requirements:

- Normalize booleans so `"False"`, empty strings, and NaN are false.
- Preserve existing co-crystal auto-tagging.
- Propagate `pair_source`, `is_cocrystal_benchmark`, co-crystal identity, and available source lineage through pair intent, pairlist, curation, materialization, filtering, and analysis joins.
- Append metadata compatibly and update strict validators.
- Never fabricate receptor/ligand source lineage from filename resemblance.
- Do not remove or reinterpret legacy classification heuristics under this engineering requirement; SR-004 governs that semantic change.

## 8. Gate 2 Scientific Work Packages

### SR-001 — Atom correspondence and RMSD

**Approved decision:** shared canonical topology plus symmetry-aware graph mapping; every unproven or incomplete mapping is `not_comparable`. See `decisions/atom-mapping.md`.

Required contract:

- Heavy-atom policy, hydrogen policy, graph/bond representation, symmetry treatment, alignment, cutoff, mapping version, mapped/total atoms, and coverage are recorded.
- Exact atom order is a fast path only when identity is proven.
- Element multisets or centroid-radius sorting do not establish correspondence.
- Same-formula constitutional isomers are `not_comparable`.
- Atom-count mismatch is `not_comparable` unless an explicit approved common-atom mapping exists.

This requirement supersedes the centroid-sorted agreement behavior accepted by Spec 027.

### SR-002 — Versioned consensus policy

**Approved decision:** equal-weight within-engine/within-target complete-case rank-percentile consensus; geometry is comparable-only QC/tie-breaker; no missing-component renormalization. See `decisions/consensus-v2.md`.

Required contract:

- Preserve legacy `dockbox_geometric` for reproducibility and add a distinct formula version/mode.
- Use per-engine/per-target normalized signals; no raw cross-target affinity comparisons.
- Keep proven disagreement distinct from `not_comparable` geometry.
- Do not impute an invented neutral geometry value.
- Emit formula version, score direction, components/weights, normalization, missing-data policy, tie rule, engines requested/included/missing, pose-selection rule, RMSD cutoff, and mapping version.
- Every selector and report must use the declared score direction. The current direction inconsistency between ranking and top-pose selection must be resolved.

The 7MYL/T2Z14 versus 4ZBE/OX-11 cross-target inversion is not a valid pass/fail criterion. A sanitized 21-complex dataset may be used only for a before/after comparison report unless truth labels and an approved validation design are supplied.

### SR-003 — Binding-site-center provenance

**Approved decision:** selected holo-ligand heavy-atom centroid; apo/no-reference cases require explicit user coordinates or return `skipped_missing_configuration`. See `decisions/grid-center.md`.

Required output fields:

- `ligand_centroid`
- `contact_residue_centroid`
- `binding_site_center`
- `center_method`
- contact definition and selection identity

Both active-site callers must consume `binding_site_center`. `overall_center` and `ligand_center` remain compatibility aliases for one migration window. Contact-atom centering is not in scope until PLIP atom/ring identifiers are parsed with a versioned, header-aware contract.

### SR-004 — Reference-ligand classification policy

**Approved decision:** explicit metadata is authoritative; scientific filename/text inference is removed by default; an opt-in warned legacy candidate/report mode cannot establish a current validation anchor before confirmation. See `decisions/reference-classification.md`.

Required contract:

- Explicit pair/extraction provenance is authoritative for scientific reference and self-redocking classification.
- Broad substring matching such as `_ligand_` cannot silently establish reference status.
- Imported/manual pairlists require explicit annotation to establish reference status.
- Missing lineage remains `unclassified`/non-reference with a warning.
- If an anchored migration mode is approved, emit its method/version and never represent inferred lineage as observed provenance.
- Consolidate reference/source-layer classification behind one policy implementation after metadata propagation is complete.
- Preserve legacy result reproducibility through an explicit mode or migration report; do not silently rewrite historical meaning.

## 9. P1 Operational and Packaging Requirements

### OR-001 — Exact deploy input

Add `--pairlist-file PATH` (or an approved equivalent) at the deploy layer. It deploys exactly those rows, is mutually exclusive with `--from-rerun-manifest`, records path and SHA-256, and does not change canonical `pairlist.csv` or pair-curation state. Deployment metadata may update normally. No remote submission is part of the test.

### OR-002 — Cache correctness

Cover four distinct invalidation surfaces:

1. GNINA intermediate `all_scores.csv` freshness.
2. Per-engine `normalized_scores.csv` freshness.
3. DAG raw-score input declarations and fingerprints.
4. Engine-layout/detection cache fingerprints.

Changing a pose, log, pairlist, pair-intent file, parser-relevant setting, or requested force flag must invalidate affected downstream artifacts. Provenance uses SHA-256; cache keys may use the repository’s documented deterministic fingerprint.

### OR-003 — Read-only per-artifact layout resolver

Resolve poses, logs, and existing scores independently. Detection must not create directories. Return selected candidate, populated alternatives, precedence reason, and ambiguity warnings. Conflicting populated candidates must not be chosen silently. Keep `ensure_engine_layout()` as the write-target API.

### OR-004 — Packaging baseline

Reconcile root and post-docking dependencies or define tested extras. Repair/remove the stale entry point and test every declared console script from a clean environment.

### OR-005 — Metadata lineage

Add source-structure lineage only when it is known from extraction/preparation provenance. Filename resemblance is not sufficient. Preserve optional column compatibility by appending metadata fields and updating schema validators.

## 10. P2 and Linked Scientific Communication

- Document and emit deterministic HPC `submit_mode`; do not submit during verification.
- Add backward-compatible CLI aliases and help text.
- Treat optional multi-engine CNN integration as a separate approved methodology version.
- Label fixed-sphere, residue-count, charge-ratio, and druggability outputs as heuristic descriptors.
- Make `classification_basis`, validation state, reference identity, and percentile fallback warnings prominent.
- Track a future known-active/decoy benchmark interface; docking remains hypothesis generation.
- Document score incommensurability across engines.

## 11. Status and Artifact Contracts

Every stage must distinguish at least:

- `completed`
- `failed`
- `skipped_disabled`
- `skipped_missing_dependency`
- `skipped_network`
- `skipped_missing_configuration`
- `not_comparable`

Missing required output is never success. Optional stages must record why they were skipped. Run manifests must retain input SHA-256 values, tool versions and executable/container identity, seed, box, scoring mode, preparation parameters, formula/mapping versions, and engine inclusion/failure state.

The canonical numbered topology and `.workflow/state.json` transition contract remain unchanged.

## 12. Verification Gates

| Gate | Required evidence |
|---|---|
| 1 Specification approval | Scientific lead, date, approved scope |
| 2 Architecture/scientific review | Decisions SR-001 through SR-004, including exact policies |
| 3 Tasks generated | `tasks.md` dependencies and commands reviewed |
| 4 Implementation | Surgical changes on an approved remediation branch |
| 5 Syntax/format/lint/type | `compileall` and `git diff --check`; configured lint/type tools must be added before claiming those checks |
| 6 Unit/contract tests | Focused pytest suite, web UI regression, smoke harness |
| 7 Minimal fixture | Offline 1IEP/Vina post-docking fixture |
| 8 Schema/topology | Pairlist, score, run-tracking, manifest, and numbered topology assertions |
| 9 Workflow/cache/resume | Unchanged rerun cache hit; mutated source invalidation; force behavior |
| 10 Provenance | Manifest and SHA-256 validation, including scientific parameters and unknowns |
| 11 Human scientific review | Before commit/push/merge or any external execution |

## 13. Traceability Summary

| Findings | Requirements | Planned tasks/tests |
|---|---|---|
| F-01/F-02, NMR-04 | FR-001, FR-002, OR-003 | Reproduction, resolver, explicit zero-row diagnostics |
| F-03 | FR-003 | Pytest remediation suite and smoke signature repair |
| F-06 | FR-004 | Executed plot regression |
| F-07 | FR-005 | Windows/POSIX schema tests |
| NMR-06 | FR-006 | Shared pose extraction and selected-pose redocking test |
| NMR-05, M-01 | FR-007, OR-005, SR-004 | Atomic metadata propagation followed by approved classification/migration tests |
| NMR-09 | SR-001, SR-002 | Mapping decision, symmetry tests, versioned consensus tests |
| NMR-07 | SR-003 | Center decision, provenance fields, caller migration tests |
| NMR-01 | OR-001 | Exact three-row local deployment test |
| NMR-03 | OR-002 | Four-surface cache invalidation tests |
| F-04/F-05 | OR-004 | Clean-install and console-entry tests |
| NMR-02/NMR-08/M-02 | Section 10 | Deferred documentation/methodology tasks |
