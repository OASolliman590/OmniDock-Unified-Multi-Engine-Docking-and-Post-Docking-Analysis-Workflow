# Tasks: Spec 036 — Scientific Consistency and Validity

- [x] **T001** Scientific Lead decisions D1–D5 recorded (2026-10-10).
- [x] **T002** (R1) Project-wide protonation policy; MD export reuses the docked microspecies; microspecies in the best-pose table.
- [x] **T003** (R2) Redocking validation status marked on downstream rows and reports.
- [x] **T004** (R3) Consensus mode inventory and removal; v2 default; explicit single-engine ranking; retire legacy RMSD paths.
- [x] **T005** (R4) One shared pose selector and the GNINA regression test.
- [x] **T006** (R5) Rg-scaled box, required seed with replicates and reproducibility metric, unified exhaustiveness, explicit `energy_range`.
- [x] **T007** (R6) Software fixes: single-engine best-pose table, tag stems, exit status, centre provenance check, legacy script retirement, single documented path, generated md-input maps.
- [ ] **T008** Partial: full suite 343 passed / 2 skipped and contract smoke pass at the implementation commit. Open follow-ups: protonation_policy dropped by bootstrap_project_layout; full dry-run smoke fails at `prep project` (pairlist DictWriter fieldnames); pose_extractor GNINA cnn_affinity export bug; best_engine_per_complex compares raw cross-engine affinities. Full suite, smoke test and evidence.
- [ ] **T009** Bounded 1IEP/STI re-dock on the new defaults; software and scientific evidence kept separate.
- [ ] **T010** Scientific Lead review (agent does not mark done).
