# Omni-DockForge Remaining Roadmap (Post-Reconciliation)

Date: 2026-03-31  
Source of truth: `specs/018-dockforge-omnidock-platform/tasks.md`

## Current Position

- Completed tasks: 79 / 106
- Open tasks: 27 / 106

Milestone status:
- Phase 1 Foundation: 6/6 complete (0 open)
- Milestone 1 UX scaffolding: 17/17 complete (0 open)
- Milestone 2 preparation refactor: 16/16 complete (0 open)
- Milestone 3 docking orchestration: 18/18 complete (0 open)
- Milestone 4 post-docking promotion: 13/13 complete (0 open)
- Milestone 5 QC + biology + rescoring: 6/15 complete (9 open)
- Milestone 6 reproducibility + storage + reporting: 0/16 complete (16 open)
- Cross-cutting docs/release: 3/5 complete (2 open)

---

## Priority Roadmap

### Wave 1: Stabilization and UX Closeout (Completed)

Status: Completed on 2026-03-31.

Delivered:
- feature flags scaffolding in workflow state/model contracts
- migration notes document for old-to-new workflow mapping
- deterministic run-id/timestamp utility and adoption in workflow runtime
- checkpoint metadata lineage scaffold and checkpoint id tracking
- interactive/CLI naming cleanup for Omni-DockForge
- smoke coverage for non-blocking background execution
- smoke coverage for universal back navigation prompt contracts

Next critical wave: Wave 3.

---

### Wave 2: Preparation Hardening and Complex Integrity (Completed)

Status: Completed on 2026-03-31.

Goal: remove remaining preparation ambiguity and finalize structure-valid downstream artifacts.

Delivered:
- PLIP preparation-surface removal finalized across workflow CLI/interactive/execution contracts
- preparation docs updated to keep interaction-analysis tooling out of preparation scope
- pH validation contract enforced for preparation entry points
- per-ligand preparation step reports auto-generated and persisted
- mode-specific ligand output contract validation implemented and aggregated
- preparation preflight summary artifacts emitted in CSV/JSON
- complex export validation hardening integrated into extraction/index manifests

Exit criteria:
- Preparation outputs are deterministic and validated per profile.
- Complex exports are chain/residue/pose-safe for downstream rescoring/visualization.

---

### Wave 3: Scientific Maturation (QC/Benchmarking/Biology/Rescoring)

Goal: complete the scientific-strength layer for defendable hit triage.

Scope:
- T072 upstream ligand QC + ADMET integration
- T073 configurable QC thresholds in schema/config loaders
- T074 decoy/reference benchmarking (EF1%, EF5%, ROC-AUC)
- T076 classifier outputs integrated into reporting
- T079 biology correlation mirrored into analysis/report surfaces
- T080 plugin interface upgrade for rescoring modules
- T081 OnionNet2 plugin scaffold
- T082 capability-gated OnionNet2 execution/fallback
- T085 unresolved biology-mapping smoke coverage

Exit criteria:
- Score/QC/ADMET/biology channels remain separate and auditable.
- Optional rescoring is pluggable and failure-tolerant.

---

### Wave 4: Reproducibility and Platform Finalization

Goal: make runs fully auditable and future-ready (storage, parity, extension hooks).

Scope:
- T086 `.meta/config.yaml`
- T087 `.meta/run_manifest.json` with checksums + versions
- T088 `.meta/env.lock`
- T089 checkpoint lineage tracking
- T090 SQLite schema module
- T091 dual-write CSV + SQLite
- T092 parity validator (CSV vs SQLite)
- T093 canonical numbered output consolidation
- T094 `7-Reports/START_HERE.md` generation
- T095 consolidated run summary output
- T096 ensemble receptor extension stubs
- T097 rerun-manifest auto-execution stub
- T098 dashboard-ready export contracts
- T099 lineage integrity smoke
- T100 manifest completeness smoke
- T101 SQLite parity smoke/skip rationale
- T102 README DockForge workflow map update
- T106 release checklist for staged rollout

Exit criteria:
- Every run is reproducible, queryable, and explainable.
- Future extension points exist without architectural rewrites.

---

## Recommended Execution Order (Critical Path)

1. Wave 3 (blocks scientific defensibility)
2. Wave 4 (blocks reproducibility publication/readiness)

Parallelizable clusters:
- Cluster C (science): T072, T073, T074, T076, T079, T085
- Cluster D (plugins): T080, T081, T082
- Cluster E (storage/reproducibility): T086-T092, T099-T101

---

## Near-Term Sprint Proposal

Sprint 1:
- T072, T073, T074, T076, T079, T085

Sprint 2:
- T080, T081, T082, T086, T087, T088, T089

Sprint 3:
- T090-T092, T093-T101, T102, T106

---

## Validation Gate Commands

- `python test/test_dockforge_smoke.py --skip-prep-matrix --skip-all-engines`
- `python test/test_dockforge_smoke.py --skip-all-engines`
- `python test_pipeline.py`

Add new smoke slices in each wave before moving to the next.
