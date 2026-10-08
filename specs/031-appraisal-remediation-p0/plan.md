# Implementation Plan: Appraisal Remediation and Multi-Engine Pipeline Hardening

## 1. Plan Status

**Status**: Gate 1 and Gate 2 implementation completed and accepted by the Scientific Lead on 2026-08-21.  
**Scientific blockers**: Resolved by the four approved decision records in `decisions/`.  
**External execution**: out of scope unless separately approved.

The current checkout is `chore/install-agentic-bioinformatics-stack`, not a remediation branch. Existing user changes in `.gitignore`, `docking/project_layout.py`, and `test/test_dockforge_smoke.py` must be preserved. Spec 031 files are currently untracked; they become an implementation baseline only after human approval.

The Antigravity `implementation_plan.md` is a user-facing summary, not the source of truth. This file, `spec.md`, and `tasks.md` are the authoritative revised plan.

## 2. Branch and Change Strategy

The original branch strategy is retained as historical planning context. The current working tree contains the verified Gate 1 and Phase B implementations; closeout followed `phase_b_plan.md`:

| Phase | Proposed branch | Scope |
|---|---|---|
| A | `codex/031a-evidence-correctness` | Reproduction, discoverable tests, implementation-safe P0 fixes |
| B | `codex/031b-scientific-methods` | Only scientifically approved mapping, consensus, and center changes |
| C | `codex/031c-operational-hardening` | Exact deploy, cache/layout, packaging, documentation |

Phase B implemented only T013, T015, T017, and T019 after backend/schema contract tests. Phase C remains present in the current working tree and was regression-tested.

## 3. Architecture Principles

1. **Evidence before patching:** the Vina zero-row symptom is valid, but the original parser-dispatch cause is disproven. Reproduce first.
2. **Read resolution is not write layout:** use a read-only per-artifact resolver for existing outputs and retain `ensure_engine_layout()` for write targets.
3. **One canonical score contract:** engine adapters normalize into the existing canonical score columns. The GNINA legacy `all_scores.csv` schema is not silently redefined.
4. **One pose-extraction contract:** PDB/PDBQT/SDF pose selection shares 1-based ordinal behavior and explicit errors.
5. **No unproven atom mapping:** inability to prove correspondence becomes `not_comparable`.
6. **Version scientific semantics:** preserve legacy consensus output and emit formula version plus score direction.
7. **Provenance beats names:** explicit extraction/pair metadata is authoritative; legacy names are migration hints only.
8. **No missing-output success:** required empty artifacts fail at their producing stage.

## 4. Phase 0 — Evidence and Baseline Reconciliation

### A0.1 Retain a reproducible Spec 031 fixture

Create `test/fixtures/spec031/vina_1iep/fixture_manifest.json` and reference or copy the retained offline artifacts needed for post-docking regression. The manifest records accession, ligand residue/chain, old center provenance, box, engine/scoring mode, seed, exhaustiveness, pose count, hashes, and explicit unknown preparation parameters.

Do not run Vina to regenerate the fixture during tests. The fixture validates parsing and post-docking behavior only.

### A0.2 Reproduce and classify the Vina zero-row symptom

Add focused tests covering:

- canonical Vina pose path with exact pose-stem/pairlist-tag match;
- pose-stem/pairlist-tag mismatch;
- legacy nested, legacy flat, numbered/canonical, and mixed layouts;
- direct parser rows versus unified-loader rows;
- existing outputs that parse to zero rows.

Capture the exact failing layer. Only then choose among a tag diagnostic, layout resolver, cache correction, or other demonstrated fix. Do not make `SimplifiedPostDockingPipeline` generically engine-aware without evidence.

### A0.3 Restore the test baseline

Create `test/test_spec031_remediation.py` and fix the current `_prepare_non_gnina_context()` test call to pass `analysis_config`. Keep the monolithic smoke harness, but do not rely on it as the only regression suite.

### A0.4 Evidence disposition

Record every NMRbox item as `reported_unverified` until its sanitized evidence is retained. The 21-complex panel and 4ZBE/6LL5 comparisons are review reports, not automated pass/fail tests until their inputs and provenance exist.

## 5. Phase A — Implementation-Safe Correctness

### A1. Explicit score-ingestion failure

At the demonstrated canonical loader/DAG layer:

- distinguish missing outputs from present-but-unparseable outputs;
- report engine, parser, resolved paths, files examined, tag-filter decisions, and zero-row reasons;
- fail `raw_scores` immediately;
- write schema-bearing diagnostic artifacts, never a headerless empty CSV;
- block dependent required stages with recorded cause.

### A2. pandas compatibility

Replace direct `applymap` usage with a compatibility helper that supports the declared pandas range. Execute `_plot_hit_class_matrix()` in pytest and assert a non-empty image plus the correct manifest status.

### A3. POSIX relative paths

Change the three verified relative-path serialization sites to `.as_posix()`. Add platform-independent unit tests using Windows-style path semantics and generated run-tracking artifacts. Keep absolute external-path serialization explicitly separate.

### A4. Shared multi-model pose extraction

Refactor or expose a shared extractor used by geometric consensus and redocking validation:

- 1-based pose ordinal;
- no silent fallback to model 1;
- final missing `ENDMDL` may close at EOF;
- malformed/nested models fail explicitly;
- selected best-row pose is passed end to end;
- parser provenance is recorded.

### A5. Atomic reference-metadata propagation

Complete the implementation-safe part of NMR-05/M-01 first:

- preserve existing cocrystal auto-tagging;
- append optional provenance fields to canonical pairlist output and validators;
- normalize boolean values explicitly;
- propagate metadata through pair intent, curation, materialization, filters, and analysis joins;
- add regression fixtures for false-like values and metadata survival.

Do not claim source-PDB lineage when it is unavailable, and do not change legacy scientific classification semantics until B4 is approved.

### A6. Phase A fixture regression

The offline 1IEP project must reach the approved analysis scope with valid score rows, stable schemas, canonical topology, explicit optional-stage statuses, and a complete manifest. This test proves post-docking ingestion, not a new docking calculation or biological validity.

## 6. Phase B — Scientific Methods (Implemented and Accepted)

The authoritative execution sequence, test matrix, stop conditions, and file boundaries are in `phase_b_plan.md`. The four normative policies are in `decisions/`.

### B1. Atom correspondence decision and implementation

Approved policy:

- shared prepared-ligand topology/canonical atom map;
- graph-isomorphism and symmetry-aware minimum heavy-atom RMSD when graph parsing is available;
- `not_comparable` when correspondence cannot be proven;
- centroid sorting retained only as a labeled diagnostic that cannot create scientific agreement.

Required tests cover permuted atoms, symmetry-equivalent atoms, same-formula constitutional isomers, atom-count mismatch, and mapping coverage. Amend/supersede Spec 027’s 44-vs-43 soft-agreement contract.

### B2. Versioned consensus policy

Approved policy is equal-weight, within-engine/within-target complete-case rank-percentile consensus as the ranking signal, with geometry as comparable-only QC/tie-breaker. Weighted geometry, missing-component renormalization, and neutral imputation are prohibited in v2.

Implementation rules:

- preserve legacy `dockbox_geometric` for reproduction;
- add a distinct v2 name or formula version;
- no neutral numeric imputation for missing geometry;
- distinguish disagreement from `not_comparable`;
- emit score direction and use it in every selector/report;
- record engines included/failed, normalization, ties, pose selection, weights/components, mapping version, and cutoff.

Tests use multiple ligands within the same target. Cross-target raw affinity inversions are not acceptance criteria.

### B3. Binding-site center provenance

Approved holo policy: selected ligand heavy-atom centroid. Required migration:

- add `ligand_centroid`, `contact_residue_centroid`, `binding_site_center`, and `center_method`;
- update both callers to use `binding_site_center`;
- keep `overall_center` and `ligand_center` aliases for one compatibility window;
- make ligand/model/altloc/hydrogen policies explicit;
- select an explicit no-ligand behavior: user box, approved cavity result, or `skipped_missing_configuration`.

PLIP contact-atom centering is deferred until the parser retains versioned atom/ring identifiers for each interaction class.

### B4. Reference-classification migration

Approved policy: explicit extraction/pair provenance is authoritative; broad filename/text heuristics never establish scientific reference status. Legacy projects may use only an opt-in, warned candidate/report mode that cannot establish current validation anchors before explicit user confirmation.

During implementation, consolidate all reference/source-layer inference behind one policy, emit classification method/version and warnings, keep missing lineage `unclassified`/non-reference, and verify that `novel_ligand_1` is never a reference. A5 has already preserved the metadata needed to avoid accidental benchmark loss.

## 7. Phase C — Operational and Packaging Hardening

### C1. Exact scoped deployment

Implement at:

- `docking/cli.py::deploy_main`
- `workflow/cli.py::_add_deploy_args`
- `workflow/cli.py::_build_deploy_argv`

Use `--pairlist-file PATH` or an approved equivalent. It is mutually exclusive with `--from-rerun-manifest`, schema-validates rows, deploys exactly those rows, records path/SHA-256, and does not mutate canonical pairlist or curation state. Test local asset generation only; do not sync or submit.

### C2. Cache correctness across four invalidation surfaces

Correct separately:

1. GNINA `all_scores.csv` intermediate freshness.
2. Per-engine `normalized_scores.csv` freshness.
3. DAG `raw_scores` input declarations/fingerprints.
4. Engine-layout/detection cache fingerprints.

Enumerate resolved pose/log/score and pair-metadata sources as DAG inputs. Unchanged reruns are cache hits; changed raw data or relevant configuration invalidates downstream nodes. `--force-rescore` must propagate through every layer if retained as an API.

### C3. Read-only artifact resolver

Extend `post_docking_analysis/engine_hpc_adapter.py::detect_engine_layout()` or a shared resolver to independently detect pose, log, and existing-score sources. Return all populated candidates, selected candidate, precedence reason, and ambiguity. Use it in engine detection and score ingestion. It must perform no filesystem writes.

### C4. Packaging and entry points

Reconcile dependency manifests or define installable/tested extras. Repair or remove `pdb-prepare-wizard=main:main`. Test every declared console entry point in a clean environment.

### C5. Deferred communication and methodology

Document submission-mode selection, CLI aliases, heuristic descriptor limitations, classification basis/fallback, and engine-score incommensurability. Treat optional multi-engine CNN integration and known-active/decoy benchmarking as separate approved methodology work.

## 8. Primary File Map

| Area | Files/symbols |
|---|---|
| Score diagnosis/ingestion | `docking/runners/vina.py`, `post_docking_analysis/unified_analysis.py`, `multi_engine_pipeline_impl.py`, `engine_hpc_adapter.py`, `engine_detector.py`, `docking_parser.py` |
| Empty-stage handling | `multi_engine_pipeline_impl.py` DAG raw/normalized/consensus nodes |
| Test baseline | `test/test_spec031_remediation.py`, `test/test_dockforge_smoke.py` |
| pandas | `post_docking_analysis/visualization_suite.py::_plot_hit_class_matrix` |
| POSIX paths | `multi_engine_pipeline_impl.py`, `simplified_pipeline_impl.py`, `workflow/state.py` |
| Pose extraction/RMSD | `unified_analysis.py`, `geometric_consensus.py`, `redocking_validation.py` |
| Reference provenance | `pairlist_builder.py`, `pair_curation.py`, `simplified_input_handler.py`, `simplified_pipeline_impl.py`, `redocking_validation.py` |
| Consensus | `consensus.py`, `top_pose_selector.py`, downstream selectors/reports |
| Center provenance | `core_pipeline.py`, `batch_pdb_preparation.py`, `cli_pipeline.py` |
| Exact deploy | `docking/cli.py`, `workflow/cli.py` |
| Packaging | `requirements.txt`, `post_docking_analysis/requirements.txt`, `environment.yml`, `setup.py` |

## 9. Compatibility and Migration

- Preserve the legacy consensus mode and add a versioned replacement.
- Preserve legacy center aliases for one compatibility window.
- Append optional pairlist metadata fields after required fields and update strict validators.
- Keep write layout canonical; read legacy/mixed artifacts without silently moving them.
- Do not choose silently between conflicting populated layout candidates.
- Record legacy reference-name inference as a warning and migration source.
- Preserve existing comments/docstrings and avoid unrelated refactors.

## 10. Risks and Mitigations

| Risk | Mitigation |
|---|---|
| Fix targets the wrong Vina layer | Reproduction gate with direct parser versus canonical loader evidence |
| Invalid atom mapping produces plausible RMSD | Graph/symmetry-aware proof or `not_comparable`; no centroid-derived agreement |
| Ranking semantics silently change | Versioned formula, explicit direction, full component provenance |
| Missing geometry becomes fake evidence | No neutral numeric imputation; separate missing/disagreement states |
| Grid boxes move | Gate 2 approval, old/new comparison report, center provenance and aliases |
| Metadata removal hides legacy redocking | Atomic propagation/migration plus anchored warning fallback |
| Layout resolver chooses stale/conflicting artifacts | Read-only detection, deterministic precedence, ambiguity failure/warning |
| Cache remains stale at another layer | Test GNINA intermediate, normalized, DAG, and detection layers separately |
| Strict pairlist consumers break | Append optional columns and update schema compatibility tests |
| Tests overstate scientific validity | State explicitly that docking is hypothesis generation and fixtures prove software contracts |

## 11. Verification Commands and Gates

Manual Gates 1, 2, and 11 require approver, date, and referenced evidence.

```powershell
# Gate 5: configured repository checks
python -m compileall -q .
git diff --check

# Gate 6: focused and existing unit/contract suites
python -m pytest -q test/test_spec031_remediation.py
python -m pytest -q test/test_webui.py
python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix

# Gate 7: offline minimal bioinformatics fixture
python -m pytest -q test/test_spec031_remediation.py -k "vina_1iep"

# Gate 8: schema and canonical topology
python -m pytest -q test/test_spec031_remediation.py -k "schema or topology or posix"

# Gate 9: cache, resume, and force behavior
python -m pytest -q test/test_spec031_remediation.py -k "cache or resume or force_rescore"

# Gate 10: manifest, checksum, and scientific provenance
python -m pytest -q test/test_spec031_remediation.py -k "manifest or checksum or provenance"

# Final discoverability/regression check
python -m pytest -q
```

`compileall` is syntax checking only. The repository currently has no configured formatter, linter, or type checker. Gate 5 must not be reported as full format/lint/type coverage unless pinned tooling/configuration is approved and added.

## 12. Completion Criteria

The remediation is complete only when:

1. Every implemented finding maps to a requirement, task, test, and retained evidence artifact.
2. No required output is missing or headerless.
3. All relevant stage statuses and skip reasons are explicit.
4. The 1IEP fixture passes offline with complete provenance.
5. Cache invalidation and unchanged-input cache hits are both demonstrated.
6. Schema, topology, resume, and manifest gates pass.
7. Scientific-method changes carry approved decision records and versioned semantics.
8. The scientific lead completes Gate 11 before commit/push/merge or any external execution.
