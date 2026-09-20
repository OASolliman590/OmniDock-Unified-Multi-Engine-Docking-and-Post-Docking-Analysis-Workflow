# Historical release-checklist snapshot

This checklist records an earlier release process. Unchecked or checked items do
not describe the current tree, and historical completion claims have not been
rerun. Use [docs/testing.md](testing.md) for current gates.

# DockForge Release Checklist (Feature-Flag Rollout + Fallback)

Use this checklist before cutting a release branch/tag for Omni-DockForge.

## 1) Scope and Release Type

- [ ] Confirm release type: `patch` / `minor` / `major`.
- [ ] Confirm target audience: internal only / pilot users / broad rollout.
- [ ] Confirm migration impact for legacy layout projects.

## 2) Feature-Flag Plan

Current DockForge rollout flags (from `workflow/models.py`):

- `enable_timeline`
- `enable_background_tasks`
- `enable_checkpoint_revise`
- `enable_execution_environment_adapters`
- `enable_post_docking_first_class`
- `enable_sqlite_dual_write`

Release decision matrix:

- [ ] Decide default-on/default-off state for each flag.
- [ ] Document any temporary overrides in release notes.
- [ ] Confirm `enable_sqlite_dual_write` policy (recommended: default `false` until parity is validated on target environment).

## 3) Verification Gates

Run required smoke checks:

- [ ] `python test/test_dockforge_smoke.py --skip-prep-matrix --skip-all-engines`
- [ ] `python test/test_dockforge_smoke.py --skip-all-engines`
- [ ] `python test_pipeline.py`

Confirm milestone-specific outputs:

- [ ] `.meta/config.yaml`, `.meta/run_manifest.json`, `.meta/env.lock` produced for representative run.
- [ ] `5-Analysis/START_HERE.md` generated and points to `LATEST_SESSION`, `complexes/`, and `7-Reports/`.
- [ ] `7-Reports/START_HERE.md` generated.
- [ ] If SQLite dual-write enabled: `4-Working/storage/results.db` and `csv_sqlite_parity.csv` generated with no blocking mismatches.

## 4) Data and Backward Compatibility

- [ ] Validate at least one canonical-layout project end-to-end.
- [ ] Validate at least one legacy-layout project end-to-end.
- [ ] Confirm optional dependencies degrade gracefully (stages skip, pipeline continues).
- [ ] Confirm post-docking scoped runs (`comparison_only`, `qc_only`, `report_only`, `rescoring_only`) produce expected outputs.

## 5) Deployment and HPC Safety

- [ ] Validate deployment generation for at least one CPU engine and one GPU-aware engine profile.
- [ ] Validate profile includes partition/account or explicit override.
- [ ] Validate submit mode behavior (`slurm_array` and/or `single_job`) for target profile.

## 6) Rollback / Fallback Strategy

If release regression is detected:

1. Disable risky rollout flags in workflow state (especially `enable_sqlite_dual_write`).
2. Continue CSV-first report path (always supported).
3. Keep interaction and visualization extras optional; skip unavailable stages.
4. Re-run smoke suite with safe flags and re-tag as hotfix.

Operational fallback map:

- Timeline/UI regression -> set `enable_timeline=false`
- Background task regression -> set `enable_background_tasks=false`
- Checkpoint flow regression -> set `enable_checkpoint_revise=false`
- Environment adapter regression -> set `enable_execution_environment_adapters=false`
- Analysis entry regression -> set `enable_post_docking_first_class=false`
- SQLite/parity regression -> set `enable_sqlite_dual_write=false`

## 7) Release Notes Requirements

- [ ] Summarize changed modules and user-visible workflow updates.
- [ ] Include any new required/optional dependencies.
- [ ] Include explicit migration notes for folder/output differences.
- [ ] Include known limitations and enabled/disabled feature flags.

## 8) Final Sign-Off

- [ ] Product/Workflow sign-off
- [ ] Scientific/QC sign-off
- [ ] Engineering sign-off
- [ ] Documentation sign-off

When all sections are checked, release can proceed.
