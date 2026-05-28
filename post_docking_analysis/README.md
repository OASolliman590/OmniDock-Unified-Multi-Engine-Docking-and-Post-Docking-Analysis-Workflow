# Post-Docking Analysis

Toolkit for analyzing GNINA, Vina, and Smina docking outputs, with both a recommended GNINA-focused simplified pipeline and an engine-aware canonical-project workflow.

## Recommended Entry Point

Use the simplified CLI for GNINA projects:

```bash
python -m post_docking_analysis.simplified_cli \
  --project-dir /path/to/GNINA_project \
  --output /path/to/output
```

Explicit folders mode:

```bash
python -m post_docking_analysis.simplified_cli \
  --sdf-folder /path/to/gnina_out \
  --log-folder /path/to/logs \
  --receptors-folder /path/to/receptors \
  --output /path/to/output \
  --pairlist /path/to/pairlist.csv
```

Useful flags:
- `--no-rmsd`: skip RMSD stage
- `--no-visualizations`: skip visualization stages
- `--ligplus-root`: run LigPlot+ from a specific LigPlus installation
- `--interactive`: prompt for missing paths/settings in a TTY session
- `--prompt-protein-names`: prompt once per detected receptor/PDB target
- `--enable-poseview`: enable PoseView REST API interaction diagrams

When `--pairlist` is omitted, the simplified CLI auto-searches for
`pairlist.csv` from the resolved project or input roots up to four parent
levels above those paths.

Legacy/config pipeline remains available:

```bash
python -m post_docking_analysis -i /path/to/results -o /path/to/output
```

## Canonical Multi-Engine Entry Point

Use the engine-aware CLI mode when the project was built with `python main.py prepare-docking` and run with `python main.py dock`:

```bash
python -m post_docking_analysis \
  --project-dir /path/to/docking_project \
  --analysis-mode comparative_all_engines \
  --output /path/to/output

python -m post_docking_analysis \
  --project-dir /path/to/docking_project \
  --analysis-mode favorite_engine_continue \
  --favorite-engine smina \
  --output /path/to/output
```

Analysis modes:
- `comparative_all_engines`: writes combined reports and unified score tables for all detected engines
- `single_engine`: filters reports to one engine
- `favorite_engine_continue`: filters to one engine and continues into deeper downstream outputs

Favorite-engine continuation behavior:
- `gnina`: uses the existing simplified structural pipeline
- `vina` / `smina`: reuse the simplified downstream stack where safe. The bridge writes filtered reports, extracts best-pose receptor-ligand PDB complexes, runs hierarchical affinity analysis, polypharmacology, general visualizations, PandaMap, structural-quality summaries, artifact consolidation, mirrors unified raw score tables under `favorite_engine/deep_analysis/raw_data/`, and runs PDB-based RMSD analysis under `favorite_engine/deep_analysis/rmsd_analysis/`

Bridge stage notes use explicit statuses so optional outputs are easier to interpret:
- `completed`
- `disabled`
- `missing_dependency`
- `missing_configuration`
- `partial_scripts_only_missing_binary`

## Simplified Pipeline Stages

1. Detect input files and load optional `pairlist.csv`
2. Build `all_scores.csv` from GNINA logs
3. Match poses to receptors
4. Build receptor-ligand complex PDBs
5. Run hierarchical binding affinity analysis
6. Run enhanced RMSD analysis (optional)
7. Generate reports
8. Generate visualization kits (optional)
9. Generate optional PoseView diagrams when enabled
10. Consolidate outputs into `visualizations/` and `raw_data/`

## Visualization Kits

1. Affinity summary plots (`matplotlib` / `seaborn`)
2. Hierarchical analysis visualizations
3. Enhanced RMSD visualizations
4. PandaMap publication-quality 2D + 3D
5. py3Dmol 3D HTML visualizations
6. ProLIF interaction maps
7. LigPlot+ interaction diagrams

Dependency contract:
- required baseline: comparative reports, score normalization, standard plots,
  complex generation, report generation
- optional extras: `py3Dmol`, `ProLIF`, `LigPlot+`, `PoseView`

Optional extras are intentionally non-blocking. The pipeline should finish and
record a skipped state when those tools are not installed or configured.

Protein naming is automatically inferred from receptor identifiers/PDB-like codes, and visualization filenames are aliased with detected display names when possible.

`top_performers.png` is a per-protein best-performer chart, not a global top-N
bar chart.

## Output Organization

Current output is split into:
- `visualizations/`: all rendered visual assets grouped by toolkit/type
- `interactions/`: canonical PandaMap, ProLIF, LigPlot, and PoseView outputs
- `raw_data/`: CSV/JSON/log/text artifacts and derived metadata

For canonical multi-engine projects, the engine-aware layer also writes:
- `reports/combined_engine_scores.csv`
- `reports/best_pose_per_tag_by_engine.csv`
- `reports/best_engine_per_complex.csv`
- `raw_data/unified_all_scores.csv`
- `raw_data/unified_best_poses.csv`

When `favorite_engine_continue` is used:
- `favorite_engine/all_scores.csv` and `favorite_engine/best_poses_unified.csv` hold the filtered engine view
- `favorite_engine/deep_analysis/` holds the continued downstream outputs for the selected engine

A linkage manifest is written to:
- `raw_data/visualization_manifest.csv`
- `raw_data/visualization_manifest.json`

Each manifest row records visualization type, original source path, inferred PDB code/display name, and linked raw files. Interaction diagrams stay under `interactions/` and are not duplicated into `visualizations/`.

## Main Files

- `simplified_cli.py`: GNINA-focused CLI
- `simplified_pipeline.py`: orchestrates simplified workflow
- `gnina_hpc_adapter.py`: auto-detects local vs HPC GNINA layout
- `publication_pandamap.py`: publication-style PandaMap generation
- `prolif_interaction_maps.py`: ProLIF map generation
- `ligplot_integration.py`: LigPlot+ integration
- `enhanced_rmsd_analyzer.py`: RMSD/clustering/diversity analysis
