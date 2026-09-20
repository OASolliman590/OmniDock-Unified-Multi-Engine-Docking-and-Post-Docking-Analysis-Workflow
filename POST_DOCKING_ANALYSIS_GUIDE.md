# Post-docking analysis

The canonical entrypoints are under `python main.py analyze`. The lower-level
`python -m post_docking_analysis` interface remains available for direct package
use, legacy layouts, preprocessing, and explicit engine/scope controls.

## Canonical modes

Comparative multi-engine analysis:

```bash
python main.py analyze comparative \
  --project-dir docking_project \
  --analysis-scope full \
  --normalization-method per_engine_rank
```

Favorite-engine continuation:

```bash
python main.py analyze favorite-engine \
  --project-dir docking_project \
  --favorite-engine gnina \
  --analysis-scope full
```

Stage-specific surfaces are `hierarchical`, `polypharmacology`, `rmsd`,
`reports`, `visualizations`, and `structure-quality`. Interaction tools are
`pandamap`, `prolif`, `ligplot`, `poseview`, and `clean`; visualization-only
tools are `py3dmol` and `pymol`. Check nested help for required arguments:

```bash
python main.py analyze stage rmsd --help
python main.py analyze interactions clean --help
```

## Inputs, outputs, and provenance

Prefer a canonical project with `project_manifest.json`, `pairlist.csv`, and
validated engine outputs. Analysis writes numbered working data under
`4-Working/`, human-facing sessions and tables under `5-Analysis/`, visualization
artifacts under `6-Visualizations/`, and reports under `7-Reports/`. Explicit
`--output` values may route a session elsewhere. Treat generated files as a
bundle: score tables, geometry status/reasons, coverage, mappings, manifests,
configuration, and logs belong together.

A deliberate score-only import requires `scores/normalized_scores.csv` plus a
matching `normalized_scores.import.json` digest marker. It is labeled
`explicit_import_unverified_execution`; it does not prove engine execution, and
missing geometry remains unevaluable.

## Score semantics

- Vina/Smina/AutoDock energy-like scores and GNINA affinity sort lower first.
- GNINA CNNscore and CNNaffinity sort higher first.
- OmniDock consensus score is higher-is-better.

Comparative analysis normalizes within engine/scoring function and target/site.
Do not rank raw unlike scoring functions together or interpret their numeric
offsets as affinity differences. Reference anchoring requires the same engine
and scoring function. Docking scores are predictions, not measurements.

## RMSD and evaluability

Redocking RMSD measures displacement in the common native receptor coordinate
frame. It does not superpose the predicted ligand onto the reference, because
that would hide a misplaced pose. RDKit graph identity/symmetry plus an
authoritative mapping determines atom correspondence. The requested pose/model
is parsed exactly.

Missing/malformed poses, mismatched graphs, absent mappings, coordinate-only
PDB/PDBQT chemistry, unsupported legacy AD4 DLG geometry, or missing reference
frame provenance produce an explicit `not_evaluable` reason. They are not
threshold passes. Coverage must include expected, evaluated, unevaluable, and
missing comparisons; absence of evidence cannot increase confidence.

## Interactions and visualizations

The clean interaction pipeline can perform complex fixes and optional layered
PLIP/ProLIF extraction. Use a dry run before external scripts:

```bash
python main.py analyze clean --project-dir docking_project --dry-run
```

PLIP, ProLIF, LigPlot+, PoseView, PandaMap, PyMOL, and py3Dmol are optional
stages with separate dependencies. An interaction omitted by one tool is not
proof that it is absent. Visualizations are inspection aids; verify atom naming,
bonding, protonation, receptor frame, and selected pose. See
[interaction types](PLIP_INTERACTION_TYPES_GUIDE.md) and
[visualization notes](post_docking_analysis/VISUALIZATION_GUIDE.md).

## Limitations

Software gates do not establish prospective enrichment, affinity accuracy,
clinical usefulness, pose physical validity, or robustness across targets.
`physical_validity_status` remains `not_evaluated` unless a named validated
estimator supplies evidence. Report engine/tool versions, exclusions, failures,
coverage, and `not_evaluable` reasons. For implemented corrections and pending
scientific validation, see [docs/correctness-update.md](docs/correctness-update.md).
