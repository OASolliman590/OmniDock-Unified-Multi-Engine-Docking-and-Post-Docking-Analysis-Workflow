# Historical snapshot — 2026-03-31

This record preserves a release-candidate assessment made on the date shown.
Its claims were not rerun for the current tree and are not current readiness,
test, benchmark, or scientific-validation evidence. Use [docs/README.md](README.md)
and live CLI help for current guidance.

# Omni-DockForge Release Candidate Summary

Date: 2026-03-31  
Version target: RC (to be finalized)  
Project: Omni-DockForge (formerly PDB Prepare Wizard)

## 1) Release Scope

This release candidate includes the Wave 4 platform-finalization track:

- Reproducibility manifests under `.meta/`:
  - `config.yaml`
  - `run_manifest.json`
  - `env.lock`
- Optional SQLite dual-write backend for post-docking comparative outputs (feature-flag gated).
- CSV-vs-SQLite parity validation path and explicit skip rationale when SQLite is disabled.
- Canonical numbered output layout enforcement and report indexing:
  - `0-Input` ... `7-Reports`
  - `7-Reports/START_HERE.md`
- Dashboard-ready export contract and JSON artifact index.
- Consolidated run summary outputs.
- Extension scaffolds for:
  - ensemble receptor support contracts
  - rerun-manifest auto-execution interface stub

## 2) Validation Evidence

Engineering verification gates executed successfully:

1. `python test/test_dockforge_smoke.py --skip-prep-matrix --skip-all-engines` ✅
2. `python test/test_dockforge_smoke.py --skip-all-engines` ✅
3. `python test_pipeline.py` ✅

Related run record:
- `docs/dockforge_release_checklist_run_20260331.md`

Notes:
- The receptor-QC failure trace shown in smoke output is an expected negative-gate contract test and does not indicate suite failure.

## 3) Feature-Flag Rollout Recommendation

Recommended defaults for RC:

- `enable_timeline=true`
- `enable_background_tasks=true`
- `enable_checkpoint_revise=true`
- `enable_execution_environment_adapters=true`
- `enable_post_docking_first_class=true`
- `enable_sqlite_dual_write=false` (recommended default until project-level parity is explicitly required)

Rationale:
- CSV path remains the stable baseline.
- SQLite path is implemented and validated, but default-off reduces operational surprise during RC rollout.

## 4) Risks and Mitigations

Primary operational risks:

- Environment variability for optional tools (OpenBabel, Meeko, PandaMap, PLIP, etc.).
- HPC profile differences (partition/account/submit mode).
- User confusion when enabling SQLite dual-write without expecting additional storage artifacts.

Mitigations in this RC:

- Optional-stage graceful degradation (skip with explicit notes).
- Slurm/profile validation contracts.
- Clear fallback strategy documented in:
  - `docs/dockforge_release_checklist.md`
  - `docs/dockforge_release_checklist_run_20260331.md`

## 5) Required Manual Decisions Before Tag

Pending non-engineering items:

1. Release type decision: `patch` / `minor` / `major`
2. Audience scope: internal pilot vs broad rollout
3. Final sign-offs:
   - Product/Workflow
   - Scientific/QC
   - Engineering
   - Documentation

## 6) Suggested Tagging and Rollout Sequence

1. Freeze branch for RC.
2. Set feature flags to recommended defaults.
3. Run one canonical project pilot and one legacy-layout pilot.
4. Collect sign-offs.
5. Tag RC.
6. Promote to general release after pilot acceptance.

## 7) Sign-Off Block

- Product/Workflow: ____________________  Date: __________
- Scientific/QC: _______________________  Date: __________
- Engineering: _________________________  Date: __________
- Documentation: _______________________  Date: __________

Decision:
- [ ] Approve RC and tag
- [ ] Approve RC with conditions
- [ ] Hold release pending action items
