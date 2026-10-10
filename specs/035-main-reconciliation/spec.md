# Spec 035 — Reconciling this branch with `main@7fc130a` (design)

**Status:** design only (Scientific Lead direction 2026-10-10, "design it to fit any case"). No merge, rebase or code port has been performed. Implementation needs a separate go-ahead.

## Situation
This branch (Spec 031–034) and `main` (`7fc130a`, `8e36027`) diverge at `cfe3305`. A trial merge conflicts in 15 files. See `specs/033-prep-layout-integrity/decisions/main-divergence.md`.

## Design: per-change reconciliation, not a blanket merge
1. **Inventory.** Split `7fc130a` into its units of change (module or function level) and classify each one:
   - `operational`: packaging, CI, locking, layout, CLI plumbing or logging;
   - `scientific`: it changes a computed value, selection, threshold, protonation or chemistry;
   - `evidence`: the `audit/` reports and probes.
2. **Authority rule.** Where a unit overlaps a method the Scientific Lead accepted (Spec 031 mapping/RMSD, consensus v2, binding-site centre, reference policy; Spec 032 export gates; Spec 033 receptor preparation), the accepted method is authoritative. The `main` version is ported only if a side-by-side comparison shows it fixes a defect, and that port is a scientific decision recorded per unit.
3. **Operational units** are ported directly in small commits, each with its tests from `main` (`test/test_*_correctness.py`). If a `main` test asserts behaviour that contradicts an accepted method, it is recorded as a conflict, not "fixed" to pass.
4. **Scientific units without overlap** (for example the job contract, reference chemistry and asset identity) get a one-page comparison and a Scientific Lead decision before porting.
5. **Evidence units** (`audit/`) are imported only after a public-repository review: no local paths, private runs or unpublished results.
6. **Finish.** Once every unit is ported, rejected or deferred, `main` is merged into this branch with all conflicts resolved toward the decisions. Each case falls into one of these, so the same procedure works for any future divergence.

## Deliverables (when approved)
- `inventory.md`: each unit with its class, overlap, proposed action and test.
- One decision record per scientific unit.
- Small port commits, each followed by a green full suite and smoke run.

## Known overlap to examine first
`post_docking_analysis/consensus.py`, `geometric_consensus.py`, `redocking_validation.py`, `top_pose_selector.py`, `multi_engine_pipeline_impl.py`, `artifact_graph.py`, `core_pipeline.py`, `docking/project_layout.py`, `test/test_dockforge_smoke.py`.
