# Visualization stages

The canonical analysis guide is
[POST_DOCKING_ANALYSIS_GUIDE.md](../POST_DOCKING_ANALYSIS_GUIDE.md). This page
only covers package visualization routing.

List visualization surfaces:

```bash
python main.py analyze visuals --help
python main.py analyze stage visualizations --help
```

`analyze visuals` exposes PyMOL and py3Dmol routes. The broader stage can
generate plots from an analysis session; optional interaction backends may also
produce 2D maps or HTML views. Outputs normally live below the project
`6-Visualizations/` topology or an explicit session output.

Before interpreting an image, verify the source pose/model, ligand identity,
bond orders, protonation, receptor frame, and failed/omitted stages. A rendered
contact or appealing pose is not experimental evidence. Do not hide redocking
displacement by ligand-on-reference fitting; RMSD evaluation uses the native
receptor frame.

Optional programs and Python libraries must be installed separately. Retain
their version/configuration and the source structure with generated artifacts.
