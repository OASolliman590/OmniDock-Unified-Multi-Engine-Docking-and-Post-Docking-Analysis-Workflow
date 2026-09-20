# Post-docking analysis package

This directory implements engine detection, score parsing and normalization,
consensus/rerun selection, RMSD and structure-quality evaluation, interactions,
visualizations, and reporting. User-facing guidance lives in the top-level
[post-docking analysis guide](../POST_DOCKING_ANALYSIS_GUIDE.md).

Preferred canonical-project commands:

```bash
python main.py analyze comparative --project-dir docking_project
python main.py analyze favorite-engine --project-dir docking_project --favorite-engine gnina
```

The direct package interface remains supported for explicit input/output,
preprocessing, legacy layouts, and detailed mode/scope selection:

```bash
python -m post_docking_analysis --help
```

Do not compare unlike raw scoring functions. Preserve direction semantics,
provenance, coverage, and `not_evaluable` reasons. Redocking RMSD requires an
authoritative graph/mapping and a common receptor frame; ligand-on-reference
fitting is not used.

Specialized pages:

- [Visualization stages](VISUALIZATION_GUIDE.md)
- [Examples](examples/README.md)
- [Archived earlier-state documents](docs/_archive/README.md)
- [Top-level interaction guide](../PLIP_INTERACTION_TYPES_GUIDE.md)

For current implemented scientific behavior and limits, see
[the correctness update](../docs/correctness-update.md).
