# Tasks: Spec 034 — Ligand Lineage and Protonation Integrity

- [x] **T001** Scientific Lead directs lineage mapping, a protonation policy that fits any case, and explicit receptor pH entry (2026-10-10).
- [ ] **T002** Scientific Lead decides the redocking RMSD frame (`decisions/redocking-rmsd-frame.md`).
- [ ] **T003** (R1a–b) Meeko SMILES/IDX lineage reader and `meeko_smiles_idx_lineage_v1` mapping source.
- [ ] **T004** (R1c–d) Spec 032 atom_map derivation from lineage and redocking reuse.
- [ ] **T005** (R2a–b) Fix the protonation no-op; measured charge, charged atoms and microspecies in provenance.
- [ ] **T006** (R2c–d) `ph_model` / `explicit_state` / `as_input` policies, state map file, interactive confirmation, and microspecies carried to best-pose provenance.
- [ ] **T007** (R3) Explicit receptor pH and force field entry (CLI required, interactive prompt); `ph_source` in provenance.
- [ ] **T008** Focused tests for acceptance items 1–4; full suite and smoke test; evidence with exact counts.
- [ ] **T009** Re-run post-docking and MD export (no new docking) on the retained Phase 2 poses; record software and scientific evidence separately.
- [ ] **T010** Scientific Lead review (agent does not mark done).
