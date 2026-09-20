# Protein-ligand interaction types

OmniDock can route complex structures through optional PLIP and ProLIF stages.
Their annotations are geometric classifications from the supplied structure;
they do not prove binding, efficacy, or experimental significance.

## Common classes

Depending on the backend and input chemistry, reports may include hydrogen
bonds, hydrophobic contacts, salt bridges, pi-stacking, pi-cation contacts,
halogen bonds, water bridges, and metal coordination. Each class depends on
backend-specific distance/angle rules and atom perception. Compare like
versions/configurations only.

## Run the supported path

Plan the clean pipeline first:

```bash
python main.py analyze clean --project-dir docking_project --dry-run
```

The clean path can include complex/format fixes plus layered PLIP and ProLIF.
Use `--no-plip`, `--no-prolif`, `--skip-chem-fixes`, or
`--skip-layered-plip` only when the omission is intentional and recorded.
Specialized routes are listed by:

```bash
python main.py analyze interactions --help
```

## Interpretation checklist

- confirm the intended pose, receptor chain/model, ligand identity, and common
  coordinate frame;
- verify bond orders, hydrogens, protonation, formal charges, metals, cofactors,
  waters, and alternate locations;
- retain the backend name/version, command/configuration, source structure, and
  logs with the interaction table;
- distinguish a tool failure or skipped stage from a valid zero-interaction
  result;
- compare reports only after matching atom perception and configuration;
- inspect representative structures rather than relying on counts alone.

PLIP and ProLIF may disagree because their chemistry perception and geometric
rules differ. LigPlot+, PoseView, PandaMap, PyMOL, and py3Dmol are optional
presentation or inspection paths, not adjudicators of scientific truth.

See [post-docking analysis](POST_DOCKING_ANALYSIS_GUIDE.md) and
[package visualization notes](post_docking_analysis/VISUALIZATION_GUIDE.md).
