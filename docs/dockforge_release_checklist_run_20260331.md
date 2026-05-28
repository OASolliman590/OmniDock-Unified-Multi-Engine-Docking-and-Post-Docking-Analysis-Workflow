# DockForge Release Checklist Run - 2026-03-31

Source template: `docs/dockforge_release_checklist.md`

## 1) Scope and Release Type

- [ ] Confirm release type: `patch` / `minor` / `major`.
- [ ] Confirm target audience: internal only / pilot users / broad rollout.
- [ ] Confirm migration impact for legacy layout projects.

## 2) Feature-Flag Plan

- [x] Decide default-on/default-off state for each flag.
- [x] Document any temporary overrides in release notes.
- [x] Confirm `enable_sqlite_dual_write` policy (default `false` unless parity is required).

Notes:
- Active recommendation for release candidate: keep `enable_sqlite_dual_write=false` by default.
- Feature-flag fallback map is documented in the base checklist.

## 3) Verification Gates

- [x] `python test/test_dockforge_smoke.py --skip-prep-matrix --skip-all-engines`
- [x] `python test/test_dockforge_smoke.py --skip-all-engines`
- [x] `python test_pipeline.py`

Command outcomes:
- Exit code 0 for all three commands.
- Expected negative-gate trace observed in smoke output (`Receptor QC gate failed`) as part of gate-testing; suite still passes.

Milestone artifact checks:
- [x] `.meta/config.yaml`, `.meta/run_manifest.json`, `.meta/env.lock` contract validated by smoke suite.
- [x] `7-Reports/START_HERE.md` generation validated by smoke suite.
- [x] SQLite dual-write parity path validated (enabled and disabled/skip rationale paths).

## 4) Data and Backward Compatibility

- [x] Validate at least one canonical-layout project end-to-end.
- [x] Validate at least one legacy-layout project end-to-end.
- [x] Confirm optional dependencies degrade gracefully (stages skip, pipeline continues).
- [x] Confirm post-docking scoped runs produce expected outputs.

Evidence:
- Canonical and legacy workflows are covered by the DockForge smoke contracts and unified all-engines smoke (`test_pipeline.py`).

## 5) Deployment and HPC Safety

- [x] Validate deployment generation for CPU/GPU-aware profile paths in smoke/deploy contracts.
- [x] Validate profile partition/account requirement enforcement.
- [x] Validate submit mode behavior (`slurm_array` and `single_job` contract path).

Notes:
- Local contract tests passed. Cluster-side live submission validation remains environment-dependent and should be repeated per HPC target.

## 6) Rollback / Fallback Strategy

- [x] Fallback strategy documented and verified for feature-flag control.
- [x] CSV-first fallback path available when SQLite dual-write is disabled.

## 7) Release Notes Requirements

- [x] Summarize changed modules and user-visible workflow updates.
- [x] Include new optional dependencies and behavior.
- [x] Include migration notes for folder/output differences.
- [x] Include enabled/disabled feature-flag expectations.

## 8) Final Sign-Off

- [ ] Product/Workflow sign-off
- [ ] Scientific/QC sign-off
- [ ] Engineering sign-off
- [ ] Documentation sign-off

## Decision

Current status: **Engineering verification complete; ready for release-candidate sign-off**.

Blocking manual items:
- release type selection
- audience scope
- final sign-offs in section 8
