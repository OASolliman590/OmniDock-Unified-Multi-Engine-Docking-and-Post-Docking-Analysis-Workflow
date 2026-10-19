# Tasks: Spec 037

- [x] **T001** Lead direction recorded (2026-10-10).
- [x] **T002** (R1) Layered protonation policy, pair coherence flags and re-prep guard done. The prepare-once skip is done under T007 (R4a). **Open:** overrides are CLI-only (not in the interactive mode or web UI).
- [x] **T003** (R2) Rg computed at ligand preparation and read by the box step.
- [x] **T004** (R3) Basic-mode exhaustiveness refusal, v2 rank for top-overall summaries, GNINA summary label.
- [x] **T005** Tests: full suite 437 passed / 2 skipped; contract and full dry-run smoke exit 0.
- [x] **T006** Scientific Lead review. Lead decision 2026-10-19: accept after R4 (T007) is complete. R4 completed and verified 2026-10-19 (475 passed / 2 skipped); acceptance applies.
- [x] **T007** (R4) Prepare-once skip (R4a), remaining raw cross-engine uses moved to v2 (R4b), deploy parameter schema (R4c). Implemented and committed 2026-10-19. Tests: `test/test_spec037_r4a_prepare_once.py` (13), `test/test_spec037_r4b_consensus.py` (13), `test/test_spec037_r4c_deploy.py` (12). Full suite and both smokes pass. Choices for the Scientific Lead are listed in `specs/036-scientific-consistency/consensus_inventory.md` (Spec 037 R4b) and in the R4 report. The Lead review remains open under T006.
