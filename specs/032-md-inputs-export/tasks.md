# Tasks: Strict MD-Consumer Input Export

## Gate 1–3 — Approval and Contracts

- [x] **T001** Scientific Lead approves the corrected specification, named consumer profile, and non-goals.
- [x] **T002** Scientific Lead records pH, charge/protonation authority, receptor policy, export selection scope, and receptor-protonation boundary.
- [x] **T003** Freeze versioned request, row-provenance, aggregate-manifest, and validation-evidence schemas.

## Gate 4 — Implementation

- [x] **T004** Implement exact best-by-engine row selection, tag/engine allowlisting, pose bounds, native-score verification, and pose SHA-256.
- [x] **T005** Implement explicit ligand topology and receptor lineage resolvers; reject PDBQT-only chemistry and receptor sources.
- [x] **T006** Implement topology-to-pose mapping using the approved Spec 031 graph mapping contract; preserve heavy-atom coordinates and fail ambiguous mappings.
- [x] **T007** Implement the chemistry-backend capability contract and approved protonation/charge path with versions and charge authority.
- [x] **T008** Implement deterministic unique atom naming and bond-order-preserving MOL2/SDF fan-out.
- [x] **T009** Implement the strict PDB writer for receptor, ligand, and canonical `system.pdb`, including `TER`, `CONECT`, and `END`.
- [x] **T010** Extend validation with G1–G8 evidence without changing legacy visualization validation semantics.
- [x] **T011** Implement temporary-workspace generation, fail-closed atomic publication, per-row provenance, and aggregate manifest.
- [x] **T012** Register `md_inputs` paths, config artifact, optional DAG node, explicit scope, and content-provenance cache guard.
- [x] **T013** Add `analyze md-inputs` CLI parsing and dispatch with nonzero required-failure behavior.

## Gates 5–10 — Verification

- [x] **T014** Add backend-independent unit tests for selection, score tolerance, hashes, mapping, naming, PDB columns, atom counts, and statuses.
- [x] **T015** Add fixtures for non-first GNINA pose, stale source, PDBQT-only rejection, topology-backed Vina/Smina, DNA-only receptor, missing backend, and shifted PDB columns.
- [x] **T016** Add dependency-backed protonation/MOL2 tests, including explicit-H completeness, formal charge, heavy-coordinate preservation, and unusual-group review flags.
- [x] **T017** Validate JSON schemas and recompute every input/output SHA-256 from the produced manifests.
- [x] **T018** Verify DAG direct request, optional failure isolation, cache hit, pose/config/topology/receptor invalidation, and `--force`.
- [x] **T019** Run syntax, focused regression, Spec 031 regression, smoke, full relevant pytest, and `git diff --check`; record exact pass/skip/failure counts.
- [x] **T020** Retain a sanitized named-consumer fixture or mark consumer validation `reported_unverified`; never substitute an Open Babel round-trip.

Verification evidence is recorded in `implementation_report.md`.

## Gate 11 — Scientific Review

- [ ] **T021** Present method versions, selected poses/scores, pH, predicted or approved charge states, receptor policy, every G1–G8 result, missing/failed rows, provenance hashes, and limitations for Scientific Lead acceptance.
