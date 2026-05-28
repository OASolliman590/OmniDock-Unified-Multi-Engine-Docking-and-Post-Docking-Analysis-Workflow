# Visualization Guide

This guide covers the visualization toolkit used by `post_docking_analysis.simplified_cli`.

## Run Command

```bash
python -m post_docking_analysis.simplified_cli \
  --project-dir /path/to/GNINA_project \
  --output /path/to/output
```

Skip all visualization stages:

```bash
python -m post_docking_analysis.simplified_cli \
  --project-dir /path/to/GNINA_project \
  --output /path/to/output \
  --no-visualizations
```

## Visualization Kits

1. Affinity overview plots
- Output examples: `affinity_distribution.png`, `top_performers.png`, `affinity_heatmap.png`, `affinity_by_protein.png`
- Source stage: simplified pipeline overview plotting

2. Hierarchical analysis visualizations
- Output path: `analysis/visualizations/`
- Source stage: `HierarchicalDockingAnalyzer.create_visualizations()`

3. RMSD visualizations
- Output paths:
  - `rmsd_analysis/per_complex_all_poses/`
  - `rmsd_analysis/per_protein_best_poses/`
  - `rmsd_analysis/global_best_poses/`
- Per-complex scope stores all-pose matrices and heatmaps for each protein-ligand combination
- Per-protein and global scopes store best-pose matrices, summaries, and enhanced RMSD visuals

4. PandaMap publication kit
- Output path: `interactions/pandamap/`
- Subfolders: `maps_2d/`, `maps_3d/`
- Generates publication-friendly formats (`pdf`, `svg`, `png`) plus 3D views
- New metadata artifacts:
  - `pandamap_publication_analysis_summary.json`
  - `pandamap_analysis_results.csv` (includes ligand-resolution provenance + per-figure title)
  - `pandamap_capabilities.json` (installed CLI mode/options + documented feature baseline)

5. py3Dmol 3D kit
- Output path: `3d_visualizations/`
- Subfolders: `individual_complexes/`, `aggregated_by_protein/`
- Output format: interactive `.html`

6. ProLIF interaction maps
- Output path: `interactions/prolif/`
- Output format: `.png` maps per complex

7. LigPlot+ diagrams
- Output path: `interactions/ligplot/`
- One subfolder per complex with native LigPlot artifacts and exported figures
- Optional `--ligplus-root` can be used to force a specific LigPlus installation

8. PoseView diagrams
- Output path: `interactions/poseview/`
- Enabled with `--enable-poseview`
- Output formats: `.png`, `.svg`, `.pdf`

## Current Default Parameters

Simplified pipeline defaults are tuned for publication-quality output:

- PandaMap:
  - `dpi=350`
  - `formats=[pdf, svg, png]`
  - `figure_width=12`, `figure_height=9`
  - `show_surface=true`, `show_3d_cues=true`
  - `generate_text_reports=true` (when `--report` is available in installed PandaMap)
  - `estimate_delta_g=true` (when `--deltaG` is available in installed PandaMap)
- ProLIF:
  - `dpi=350`
  - `figsize=(14, 10)`
  - `max_complexes=60`
- LigPlot+:
  - `contact_type=2`
  - `hydrogenate=auto`
  - `strip_metals=auto`
  - `max_complexes=40`
  - `export_formats=[png, pdf]`
  - `export_dpi=350`

## PandaMap Feature Coverage

The wrapper now tracks two levels of capability:

1. **Documented baseline (library feature set)**
- Interaction classes (16): hydrogen bonds, π-π, cation-π, π-cation, carbon-π, donor-π, amide-π, alkyl-π, hydrophobic, ionic, salt bridges, halogen bonds, metal coordination, covalent contacts, attractive charge, repulsive charge.
- Analysis modes: single-structure interaction maps, trajectory occupancy summaries, empirical ΔG estimation.
- Outputs: 2D maps (`png`), interactive 3D (`html`), text reports.
- Inputs: `pdb`, `cif/mmcif`, `pdbqt`.

2. **Runtime capability (installed CLI version)**
- Captured in `interactions/pandamap/pandamap_capabilities.json`.
- Includes detected CLI mode (`single` vs `subcommand`), supported long options, and option-level feature matrix.

## Publication-Ready Chemistry Checklist

To ensure chemistry and residues are clearly represented in figures:

1. Use cleaned complex PDBs (fixed bond orders/atom naming) before PandaMap.
2. Confirm `pandamap_analysis_results.csv` shows `ligand_resolution_source` as `detected_from_pdb_hetatm` whenever possible.
3. Prefer `pdf`/`svg` for manuscript panels; keep `png` for slides.
4. Keep 3D cues on (`show_surface=true`, `show_3d_cues=true`) unless journal style requires minimal rendering.
5. Use the generated text sidecar reports (`*.txt`) to cite interacting residues and interaction chemistry in captions/tables.

## Protein Naming in Visuals

When receptor/pose names contain PDB-like identifiers, the pipeline resolves a protein display name and prefixes visualization filenames with that name where possible. This improves readability while preserving traceability to original filenames.

Related files:
- `protein_name_mapping.csv`
- `protein_name_overrides.csv`

## Consolidated Output Model

After generating toolkit-specific outputs, the pipeline consolidates files into:

- `visualizations/`: all rendered visual assets grouped by type
- `interactions/`: canonical PandaMap, ProLIF, LigPlot, and PoseView outputs
- `raw_data/`: CSV/JSON/log/text artifacts and manifests

The linkage is tracked in:
- `raw_data/visualization_manifest.csv`
- `raw_data/visualization_manifest.json`

Manifest fields include visualization type, source file, inferred PDB code, protein display name, and linked raw artifacts. Interaction diagrams remain under `interactions/` and are not recopied into `visualizations/`.

## Troubleshooting

1. PandaMap missing
- Ensure PandaMap exists in the configured conda env (`pandamap` by default).

2. py3Dmol missing
- Install `py3Dmol` in the active environment.

3. ProLIF missing
- Install ProLIF and RDKit in the active environment.

4. LigPlot+ skipped
- Configure `LIGPLUS_ROOT` / `LIGPLUS_HOME` or pass `--ligplus-root`.

5. Matplotlib cache warning
- Set a writable `MPLCONFIGDIR` to avoid startup warnings.
