# Plan 028: Visualization Corrections and Scientific Guardrails

1. Add shared label/layout helpers in `visualization_suite.py` for protein labels, reference ligand tags, grouped ligand ordering, and dense-axis sizing.
2. Update ligand/protein-centric plots to use those helpers and reduce annotation overlap.
3. Repair scientifically weak plots by tightening input requirements, redefining ligand-level agreement, and replacing collapsed diagnostics with explicit placeholders.
4. Patch the remaining legacy cross-engine visuals so reference ligands are grouped and positive affinities are explained.
5. Add smoke checks for agreement summaries, reference labels, pose-diversity diagnostics, and validation affinity fallback joins.
