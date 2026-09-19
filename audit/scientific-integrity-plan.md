# Scientific integrity correction plan

## Scope and baseline

User-authorized bounded scientific corrections across ligand preparation,
protein preparation, and post-docking analysis. Starting remote main:
`8d4434d3a33eb83b1e12cad82944b02c83270e47`, inspected 2026-09-14.
Branch: `codex/scientific-integrity`. Remote inspection found only merged PR #1
and its source branch; no newer remote changes were available.

The design is to strengthen existing preparation/analysis contracts with small
guards and reproducible positive/negative fixtures, preserving native frames,
authoritative graphs, content provenance, fail-closed validation, and honest
unavailable statuses. It does not add a universal target preparation policy or
unvalidated physical/affinity assessment. The inherited user brief explicitly
authorizes implementation in bounded batches after independent assessment.

## Ownership and verification

- Luna max preparation researcher: current execution and primary-source matrix.
- Luna max post-docking researcher: current execution and primary-source matrix.
- Separate Luna max tester: current baseline, original failures, corrected tests.
- Luna max workers: bounded implementation after root selects supported fixes.
- Independent Astra low reviewer: scientific/code risks, no implementation edits.
- Root: scope decisions, integration, diff inspection, final verification/commits.

Cleanup task `01a09f06-8fde-7c40-a904-5eb391851a27` owns general documentation,
repository hygiene, packaging and workflow structural cleanup. It reserves
preparation and post-docking scientific execution files pending exact boundaries.
This task adds its own scientific methods and evidence documents. The parent
checkout and benchmark artifacts are read-only; the benchmark task stays paused.

## Acceptance sequence

- [x] Clone current main, pin SHA, inspect branches/PRs, coordinate cleanup.
- [x] Trace current paths and cite primary scientific/backend sources.
- [x] Reproduce high-priority failures and select bounded corrections.
- [x] Implement with targeted positive and negative chemistry fixtures.
- [x] Separate tester verifies available baseline failures and corrected behavior.
- [x] Independent Astra review; resolve material findings.
- [x] Execute local full-regression attempt, filtered regression, build, wheel,
  shell syntax and diff gates; record unavailable chemistry separately.
- [x] Commit reviewed implementation changes.
- [x] Publish the reviewed branch at `f5760bf45c69d42a4420a866c77aa8f6e74875cc`.
- [x] Execute the mandatory chemistry CI matrix on Ubuntu and Windows, for
  Python 3.10 and 3.12.
- [x] Open [draft PR #2](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/pull/2).

Tests blocked by host policy will be recorded separately from code assertions.
No tests will be weakened to manufacture pass rates. No remote benchmark or
NMRBox execution is authorized by this task.

## Selected bounded batches

1. **Receptor and structure contracts.** `structure_contract.py` must preserve
   selected TER boundaries and handle explicit LINK/SSBOND/CONECT without
   silently severing meaningful connectivity. `receptor_preparation.py` must
   request PROPKA explicitly when the user selects PDB2PQR, retain chains, and
   distinguish upstream titration from unverified final state. Conservation
   must check atom/residue identities as well as coordinates and emit an atom
   map. `test/test_scientific_receptor_contract.py` covers intact selections,
   severed links, serial renumbering, swapped identities and displaced atoms.
2. **Pose and consensus contracts.** `pose_geometry.py` must preserve SDF record
   ordinals and constrain PDBQT chemistry/mapping metadata to one valid scope.
   The root probe `work/pose_selection_probe.py` reproduced selecting the second
   record as pose one after an empty record. `consensus.py` must reject duplicate
   engine/tag inputs, inconsistent candidate identity, undeclared engine scope,
   and mixed scoring functions within one calibration group instead of silently
   aggregating. New record/consensus tests preserve valid metadata layouts,
   score directions, missing-engine coverage and unique-row behavior.
3. **Ligand contract.** Assessment covers source components, direct-call
   configuration defaults and post-transformation chemistry. The ligand worker
   owns `ligand_preparation.py`, an optional `ligand_identity.py` helper, finite
   pH/raw-profile validation at the public seam in `docking/models.py`, and new
   ligand/configuration tests. Reject disconnected ordinary inputs; validate
   heavy-atom graph/stereochemistry and original 3D coordinates after
   normalization, allowing only recorded H/formal-charge protonation changes.
   Persist the normalized graph artifact and input/output hashes. Verify mapped
   Meeko output against the actual prepared state; report unsupported mapping
   scope honestly for other backends. Do not add tautomer enumeration or an
   unreviewed alternate RMSD metric.

Each worker must demonstrate original failure before claiming a correction.
Chemistry-dependent verification blocked locally is deferred to the existing CI
environment, not reclassified as passing. The separate tester and independent
reviewer remain mandatory after implementation.

## Independent replay checkpoint

The Luna tester replayed the same five non-chemistry test modules against an
isolated archive of the pinned baseline: 15 passed and 37 failed behaviorally,
with no missing-symbol collection failures. The corrected checkout passed all
52 tests. Exact commands and per-module results are in
[the verification record](scientific-verification.md). This checkpoint does not
certify the remaining chemistry fixtures, full regression or packaging gates.

## Final local checkpoint and publication boundary

At source commit `48c0771`, the final filtered regression passed 199 tests and
failed only the two known RDKit-limited redocking smoke cases. Build, extracted
wheel verification and diff checks passed. Strict full collection remains
blocked in five chemistry modules by Windows Application Control. The two
inherited legacy shell parse errors remain; their retirement is coordinated in
the separate cleanup task, which reports reviewed commit `3c1cb00`.

The initial full filtered run revealed three empty-list consensus-scope
compatibility failures. The worker restored the established inference default,
the reviewer inspected that correction, and the final full filtered rerun
confirmed all three were resolved. See the complete verification record.

Earlier publication attempts were rejected by automatic approval review because
retrieved task history was not accepted as direct user authorization. That
historical rejection is resolved: the user subsequently authorized publication,
and the reviewed branch was pushed at
`f5760bf45c69d42a4420a866c77aa8f6e74875cc`.

The required GitHub Actions run completed successfully on Ubuntu and Windows
with Python 3.10 and 3.12. Every matrix job completed its chemistry extra
installation, full pytest (`305 passed` per job), build and wheel checks
successfully. See
[run 34876658350](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/actions/runs/34876658350).
This CI result resolves the local RDKit limitation for the CI environment; it
does not add a real-backend or benchmark validation claim. Publication is
represented by [draft PR #2](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/pull/2),
with remote `main` still at the baseline commit
`8d4434d3a33eb83b1e12cad82944b02c83270e47`. No merge was attempted.
