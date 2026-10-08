# Spec 031 — Phase B Scientific Methods Implementation Plan

**Status:** Implemented, locally verified, and accepted by the Scientific Lead on 2026-08-21.  
**Approved baseline:** Decision records dated 2026-08-21 in `decisions/`.  
**External execution:** Not authorized. Verification is offline and local.

## Scope

Implement T013, T015, T017, and T019 only. Preserve the canonical numbered topology, current engine execution parameters, legacy scientific modes, and all Gate 1 fixes. Do not add cavity detection, common-substructure RMSD, calibrated affinity claims, remote execution, or broad refactors.

## Implementation Outcome

- **B0/B1 completed:** NetworkX `>=3.2,<4` is the bounded graph backend. `spec031-atom-mapping-v1` is shared by redocking and v2 inter-engine geometry, enforces complete bond-labelled heavy-atom isomorphism, symmetry enumeration, Kabsch alignment, and non-numeric `not_comparable` failures.
- **B2 completed:** `consensus_rank_geometry_qc_v2` is an opt-in CLI/DAG mode with engine-specific pose selection, within-target complete cases, equal-weight percentiles, missing-engine disclosure, no imputation, and comparable-only geometry tie handling. Legacy modes are unchanged.
- **B3 completed:** `binding_site_center_v1` resolves explicit user coordinates or the explicitly selected holo-ligand heavy-atom centroid, records checksum/model/altloc provenance, separates contact centroids, and retains compatibility aliases.
- **B4 completed:** `explicit_reference_metadata_v2` is the shared scientific policy. Names cannot establish reference status; opt-in legacy inference produces non-anchoring migration candidates.
- **B5 completed locally:** 54 Spec 031 tests and the complete smoke suite pass; the full repository suite is 161 passed and 1 skipped. External execution remains `skipped_disabled`.

## Implementation Sequence

### B0 — Backend and Schema Spike

1. Add failing contract tests and minimal fixtures for molecular graph mapping, consensus v2, center provenance, and reference classification.
2. Audit the declared Open Babel environment against the approved graph labels, automorphism enumeration, and deterministic mapping requirements.
3. If insufficient, evaluate one narrowly scoped chemistry-aware dependency, record why it is required, version-bound it in every applicable environment manifest, and add an import/capability test.
4. Define additive schemas/method IDs before production logic. Missing backend must resolve to `skipped_missing_dependency`/`not_comparable`.

**Gate:** No production RMSD implementation until the backend proves the permuted-atom and symmetry fixtures.

### B1 — Shared Atom Mapping and RMSD

1. Introduce one shared mapping/RMSD module; do not duplicate mapping in redocking and geometric consensus.
2. Implement canonical-lineage fast path, graph-isomorphism fallback, symmetry enumeration, Kabsch alignment, and complete-coverage enforcement.
3. Replace centroid-sorted soft agreement in current scientific decisions with `not_comparable`; retain legacy diagnostic fields only under their legacy method.
4. Route redocking and inter-engine geometry through the shared result schema and propagate method/status/provenance.

**Primary files:** new shared module under `post_docking_analysis/`, `geometric_consensus.py`, `redocking_validation.py`, applicable schemas/manifests.  
**Tests:** `test/test_spec031_atom_mapping.py` plus focused integration cases.

### B2 — Consensus v2

1. Add `consensus_rank_geometry_qc_v2` without changing legacy mode formulas.
2. Centralize engine metric/direction and deterministic pose-selection policy.
3. Build the per-target complete-case ligand universe, average-rank percentiles, equal-weight composite, strict incomplete status, and declared higher-is-better direction.
4. Apply geometry only as the approved comparable-only tie-break/QC signal.
5. Update selectors, schemas, manifests, reports, and user-facing limitations so every consumer honors the same direction and method version.

**Primary files:** `consensus.py`, `top_pose_selector.py`, engine score contract helpers, pipeline/report integration.  
**Tests:** `test/test_spec031_consensus_v2.py`; synthetic multi-engine, ties, missing engine, single-engine, and row-order cases.

### B3 — Binding-Site Center Provenance

1. Add a pure heavy-atom centroid/provenance helper with explicit model/chain/residue and deterministic altloc selection.
2. Emit the new authoritative fields and make both active-site callers consume `binding_site_center`.
3. Retain `overall_center`/`ligand_center` aliases for one migration window.
4. Implement explicit-user center handling and `skipped_missing_configuration`; do not add a cavity fallback.
5. Prove box dimensions and other docking parameters remain unchanged.

**Primary files:** `core_pipeline.py`, `batch_pdb_preparation.py`, `cli_pipeline.py`, grid schemas/manifests.  
**Tests:** `test/test_spec031_grid_center.py` with holo, multi-ligand, altloc, hydrogen, and apo fixtures.

### B4 — Reference Classification and Migration

1. Add one shared reference-policy API and additive classification provenance fields.
2. Route redocking eligibility, reference anchoring, simplified analysis, and scientific visualization masks through the explicit metadata policy.
3. Separate display-name inference from scientific classification.
4. Add opt-in legacy candidate/report mode; prohibit it from creating current validation anchors.
5. Produce an old/new migration report without rewriting canonical pairlist or historical artifacts.

**Primary files:** new shared policy module, `redocking_validation.py`, `simplified_pipeline_impl.py`, `visualization_suite.py`, classification/report integration.  
**Tests:** `test/test_spec031_reference_policy.py` with explicit, incomplete, false-like, misleading-name, and legacy-mode cases.

### B5 — Integration and Closeout

1. Run the four focused Phase B suites and the existing 31 Spec 031 Gate 1 tests.
2. Run web UI, complete DockForge smoke, and full pytest regression.
3. Validate additive CSV/JSON schemas, canonical topology, cache/resume/force behavior, method/version fields, SHA-256 provenance, and legacy-mode isolation.
4. Generate an offline before/after report using retained/synthetic fixtures. Do not interpret ranks as experimental affinity or use cross-target inversions as acceptance criteria.
5. Present results and unresolved `not_comparable`, missing dependency, and unclassified cases for Human Scientific Lead review.

## Acceptance Matrix

| Package | Required proof |
|---|---|
| Atom mapping | Permutation and symmetry pass; isomer/count/topology/backend failures emit no numeric RMSD. |
| Consensus v2 | Correct engine directions; equal-weight complete-case formula; no imputation; deterministic ties; legacy unchanged. |
| Grid center | Correct heavy-atom centroid/provenance; deterministic altloc; apo skip; unchanged box parameters. |
| Reference policy | Metadata authority; misleading names remain non-reference; legacy candidates cannot anchor validation. |
| Integration | Additive schemas, canonical topology, manifests, cache/resume, reports, smoke, and full suite pass. |

## Stop Conditions

- The chosen chemistry backend cannot prove the approved mapping contract.
- Required topology lineage is unavailable for retained fixtures and no valid graph can be constructed.
- A schema change would require silently reinterpreting historical outputs.
- Any implementation would alter box dimensions, scoring modes, engine parameters, or external execution state.
- Existing user changes overlap a Phase B edit and cannot be preserved surgically.

At a stop condition, record `not_comparable`, `skipped_missing_dependency`, or the applicable explicit status and return to the Scientific Lead; do not weaken the policy.
