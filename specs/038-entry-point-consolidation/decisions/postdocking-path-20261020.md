# Spec 038 R3 — Scientific Lead decision (2026-10-20): fix all as recommended

- **D-A1:** the DAG writes the full `summary.txt`, with one owner, so the comparative writer's file is no longer retired by the reports node.
- **D-A2:** DAG top-pose outputs carry normalised scores and every accepted status column: redocking validation, pose reproducibility, docked microspecies, pair protonation and the replicate source.
- **D-A3:** favorite-engine outputs are produced by a DAG node.
- **D-B:** the DAG honours the problematic-ligand filter. It stays off by default, and its use is recorded.
- **D-C:** GNINA single-engine best-pose selection uses the shared v2 selector (Lead decision 4).
- **E-1:** fix the legacy single-engine Smina crash.

After these fixes, the legacy `run()` becomes a thin wrapper over the DAG and duplicate post-docking functions are collapsed. Every remaining golden-output difference is listed as intended.
