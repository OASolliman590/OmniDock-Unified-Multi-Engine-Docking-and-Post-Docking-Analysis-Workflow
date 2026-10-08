# Tasks: Spec 033 — Preparation and Layout Integrity

## Gate 1–3 — Approval and contracts

- [x] **T001** Scientific Lead directs root-cause fixes for P1–P3 and option A (PDB2PQR at explicit pH) for receptor protonation (2026-10-08).
- [ ] **T002** Scientific Lead decides reconciliation of this branch with `main@7fc130a` (see `decisions/main-divergence.md`). Not required for T003–T012.

## Gate 4 — Implementation

- [x] **T003** (R1) Single layout-profile resolver. Replace hard-coded `"docking_legacy"` in workflow steps (`workflow/execution.py`, `workflow/interactive.py`) with the manifest-recorded profile. A missing manifest keeps today's legacy initialisation.
- [x] **T004** (R2a) Port `receptor_preparation.py`, `structure_contract.py` and the needed `receptor_quality.py` changes from `main@7fc130a`, without touching Spec 031/032 modules.
- [x] **T005** (R2b) Route shell receptor preparation through the strict module. Remove the silent `obabel -xr` fallback and report missing backends explicitly.
- [x] **T006** (R2c) `receptor_use_pdb2pqr` defaults to `true`, with explicit pH and force field recorded in provenance. Option 1 of `decisions/pdb2pqr-heavy-atom-policy.md` (`--noopt --nodebump`, terminal-carboxylate-O-only additions) is implemented.
- [x] **T007** (R2d) Polar-hydrogen gate for Vina-family receptors.
- [x] **T008** (R2e) Declare `pdb2pqr` and `meeko` in the dependency manifests and add a capability/version check.
- [x] **T009** (R3) Write preparation config and logs into the project, never the working directory.

## Gates 5–10 — Verification

- [x] **T010** Focused tests: canonical `prep pairlist` → `dock run --dry-run`; legacy project unchanged; receptor without H rejected; missing PDB2PQR/Meeko rejected; heavy-atom conservation; no files written to the working directory.
- [x] **T011** Full pytest, contract smoke, `compileall`, `git diff --check`; record exact counts in `evidence/`.
- [x] **T012** Re-run the bounded Spec 032 Phase 2 check (1IEP/STI, Vina) on the fixed pipeline and record software and scientific evidence separately. Done: `evidence/phase2_1iep_run_20261008.md` (docking completed; Spec 031 mapping `not_comparable`; MD-export row `not_comparable` at G3).

## Gate 11 — Human review

- [ ] **T013** Scientific Lead reviews Spec 033 behavior and the Phase 2 evidence. The agent does not mark this done.
