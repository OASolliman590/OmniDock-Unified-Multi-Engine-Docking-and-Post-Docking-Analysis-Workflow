# Post-Docking Analysis Guide

## Overview

This repository provides two user-facing post-docking entry surfaces:

1. `post_docking_analysis` engine-aware canonical-project workflow (primary path)
2. `post_docking_analysis.simplified_cli` legacy wrapper (auto-delegates to unified path for manifest-backed projects)

The unified workflow includes:
- automatic local/HPC layout detection
- score aggregation from logs
- receptor-ligand complex generation
- hierarchical affinity analysis
- per-complex RMSD analysis
- multi-tool visualization generation
- consolidated output packaging (`visualizations/` and `raw_data/`)
- validation-gated hit classification with effective-policy fallback metadata

For canonical projects containing GNINA, Vina, and/or Smina results under one manifest, the unified path adds:
- engine-normalized score loading
- comparative reporting across engines
- unified score exports for downstream compatibility
- favorite-engine continuation into deeper downstream analysis

## Recommended Workflow (Simplified GNINA Pipeline)

### Quick Start

```bash
python -m post_docking_analysis.simplified_cli \
  --project-dir /path/to/GNINA_project \
  --output /path/to/post_docking_output
```

### Explicit Input Mode

```bash
python -m post_docking_analysis.simplified_cli \
  --sdf-folder /path/to/gnina_out \
  --log-folder /path/to/logs \
  --receptors-folder /path/to/receptors \
  --output /path/to/post_docking_output \
  --pairlist /path/to/pairlist.csv
```

### Common Flags

- `--no-rmsd`: skip RMSD stage
- `--no-visualizations`: skip visualization stages
- `--ligplus-root /path/to/LigPlus`: run LigPlot+ with explicit LigPlus root
- `--interactive`: prompt for missing paths/settings in a TTY session
- `--prompt-protein-names`: prompt once per detected receptor/PDB target
- `--enable-poseview`: enable PoseView REST API diagrams

Note: user-facing RMSD scope is now `per_complex` only. Legacy global/per-protein RMSD scope toggles were removed from CLI contracts.

When `--pairlist` is omitted, the simplified CLI auto-searches for
`pairlist.csv` from the resolved project or input roots up to four parent
levels above those paths.

## Canonical Multi-Engine Workflow

### Quick Start

```bash
python -m post_docking_analysis \
  --project-dir /path/to/docking_project \
  --analysis-mode comparative_all_engines \
  --output /path/to/analysis_output
```

### Continue From A Favorite Engine

```bash
python -m post_docking_analysis \
  --project-dir /path/to/docking_project \
  --analysis-mode favorite_engine_continue \
  --favorite-engine vina \
  --output /path/to/analysis_output
```

### Analysis Modes

- `comparative_all_engines`: compare all detected engines and emit unified score tables
- `single_engine`: analyze only one engine
- `favorite_engine_continue`: compare/filter first, then continue downstream from one selected engine

### Validation Gate And Effective Policy Fallback

When hit classification policy is `reference_anchor`, DockForge evaluates the redocking validation gate first.

- Validation gate artifact: `reports/validation_gate_status.json`
- If gate status is not validated, classification policy automatically falls back to `target_percentile`
- Requested vs effective policy is mirrored into:
  - `reports/consolidated_run_summary.json`
  - canonical score outputs under `4-Working/scores/consensus/`

This prevents unvalidated reference ligands from being used as class anchors.

### Top-Pose Atlas Controls

Engine-aware CLI now supports explicit top-pose controls:
- `--top-pose-policy {hybrid,best_affinity,best_consensus}`
- `--top-pose-aggregation {best_target}`

Example:

```bash
python -m post_docking_analysis \
  --project-dir /path/to/docking_project \
  --analysis-mode comparative_all_engines \
  --analysis-scope top_pose_only \
  --top-pose-policy best_affinity \
  --top-pose-aggregation best_target
```

### Favorite-Engine Continuation Behavior

- `gnina`: uses the existing simplified GNINA structural workflow
- `vina` / `smina`: use a PDBQT-based structural bridge that reuses the simplified downstream stack where safe. It writes filtered reports, extracts best-pose complex PDBs, runs hierarchical affinity analysis, polypharmacology, general visualizations, PandaMap, structural-quality summaries, artifact consolidation, mirrors unified score tables under `deep_analysis/raw_data/`, and runs PDB-based RMSD analysis under `deep_analysis/rmsd_analysis/`

Bridge stage notes are explicit about optional stage outcomes:
- `completed`: stage ran and produced usable outputs
- `disabled`: stage was intentionally turned off by configuration
- `missing_dependency`: required Python or external dependency is unavailable
- `missing_configuration`: required tool configuration such as `LIGPLUS_ROOT` is unavailable
- `partial_scripts_only_missing_binary`: scripts/manifests were written but the external renderer binary was not available

## Input Layouts

### Local Layout

```text
project/
├── gnina_out/      # SDF + .log in same folder
├── receptors/
└── pairlist.csv    # optional but recommended
```

### HPC Layout

```text
project/
├── gnina_out/      # SDF files
├── logs/           # .log files
├── receptors/
└── pairlist.csv    # optional but recommended
```

Layout detection is handled by `post_docking_analysis/gnina_hpc_adapter.py`.

## Simplified Pipeline Stages

1. Input discovery and validation
2. `all_scores.csv` generation from GNINA logs
3. Pose-to-receptor matching
4. Complex PDB creation
5. Hierarchical affinity analysis
6. RMSD clustering/diversity analysis (optional)
7. Report generation
8. Visualization generation (optional)
9. Artifact consolidation + manifest generation

## Visualization Toolkit

The simplified pipeline can generate:

1. Overview affinity plots (`matplotlib`/`seaborn`)
2. Hierarchical analysis visualizations
3. RMSD visualizations
4. PandaMap 2D/3D outputs (publication-oriented)
5. py3Dmol interactive 3D HTML
6. ProLIF interaction maps
7. LigPlot+ interaction diagrams

Dependency contract:
- required: core scoring/reporting, comparative analysis, RMSD, static plots
- optional: `py3Dmol`, `ProLIF`, `LigPlot+`, `PoseView`

If an optional dependency is missing, the pipeline records a skipped stage and
continues. Those tools are not required for a validated baseline run.

See `post_docking_analysis/VISUALIZATION_GUIDE.md` for details.

## Protein Naming and Labeling

The pipeline builds a protein naming map from receptor names and optional `pairlist.csv`.

Generated files:
- `protein_name_mapping.csv`: detected name mapping
- `protein_name_overrides.csv`: editable template for manual overrides

Visualization filenames may be aliased with resolved protein display names when identifiers contain PDB-like codes.

`top_performers.png` now reflects one best ligand per protein rather than a
global top-N ranking.

## Output Structure (Current)

```text
post_docking_output/
├── all_scores.csv
├── analysis/
│   ├── best_poses.csv
│   ├── best_per_protein.csv
│   ├── cross_protein_affinity_matrix.csv
│   └── visualizations/
├── complexes/
├── best_poses/
│   ├── strong_binders/
│   ├── moderate_binders/
│   └── weak_binders/
├── reports/
├── rmsd_analysis/
│   └── per_complex_all_poses/
│       ├── per_complex_rmsd_summary.csv
│       └── <complex_slug>/
├── interactions/
│   ├── pandamap/
│   │   ├── 2d_interaction_maps/
│   │   └── 3d_visualizations/
│   ├── prolif/
│   ├── ligplot/
│   └── poseview/
├── 3d_visualizations/
├── visualizations/                 # consolidated visual assets
└── raw_data/                       # consolidated raw artifacts
    ├── visualization_manifest.csv
    └── visualization_manifest.json
```

## Consolidated Outputs

The pipeline now performs a final packaging step:

- `visualizations/` contains copied overview, hierarchical, RMSD, and 3D visual assets
- `interactions/` remains the canonical home for PandaMap, ProLIF, LigPlot, and PoseView outputs
- `raw_data/` contains tabular/log/text artifacts
- `raw_data/visualization_manifest.*` links visualization files to related raw files

This makes downstream report assembly and auditing easier.

## Canonical Multi-Engine Output Additions

Engine-aware analysis writes:
- `reports/combined_engine_scores.csv`
- `reports/best_pose_per_tag_by_engine.csv`
- `reports/engine_summary.csv`
- `reports/best_engine_per_complex.csv`
- `reports/validation_gate_status.json`
- `reports/consolidated_run_summary.json`
- `reports/consensus_ranked_hits_with_classes.csv`
- `raw_data/unified_all_scores.csv`
- `raw_data/unified_best_poses.csv`

When `favorite_engine_continue` is used, outputs are written under:
- `favorite_engine/`: filtered engine-specific score tables
- `favorite_engine/deep_analysis/`: continued downstream artifacts for the chosen engine

### Top-Pose Ligand Performance Atlas

Comparative/favorite engine-aware runs now emit a dedicated ligand-centric atlas:

- Session-local:
  - `top_pose_ligand_performance/top_pose_per_ligand_per_protein.csv`
  - `top_pose_ligand_performance/top_pose_per_ligand_global.csv`
  - `top_pose_ligand_performance/ligand_performance_summary.csv`
  - `top_pose_ligand_performance/top_pose_selection_manifest.json`
- Canonical mirror:
  - Legacy projects: `PROJECT_ROOT/5-Analysis/top_pose_ligand_performance/`
  - Canonical layout: `PROJECT_ROOT/analysis/5-Analysis/top_pose_ligand_performance/`

The atlas selection is deterministic and policy-aware (`hybrid` default), and rows include:
- pose identity (`engine`, `tag`, `protein`, `ligand`, `site_id`, `pose`, `pose_file`)
- consensus context (`consensus_score`, engine agreement/support fields)
- interpretability context (`docking_quality_class`, `qc_status`, `admet_status`, optional `bio_*`)

Use analysis scope `top_pose_only` when you want fast regeneration of this atlas without running heavy visualization/interaction stages.

## Legacy Pipeline

Legacy/config-driven command:

```bash
python -m post_docking_analysis -i /path/to/docking/results -o /path/to/output
python -m post_docking_analysis --config my_config.yaml -i /path/to/docking/results
```

This path is still available for broader configuration scenarios, but GNINA projects should prefer the simplified CLI.
For manifest-backed canonical projects, `simplified_cli` delegates to the unified pipeline path unless legacy-only flags (`--no-rmsd` / `--no-visualizations`) are requested.

## Programmatic Usage

### Simplified Pipeline

```python
from post_docking_analysis.simplified_pipeline import SimplifiedPostDockingPipeline

pipeline = SimplifiedPostDockingPipeline(
    sdf_folder="/path/to/gnina_out",
    log_folder="/path/to/logs",
    receptors_folder="/path/to/receptors",
    output_dir="/path/to/output",
    pairlist_file="/path/to/pairlist.csv",
    run_rmsd=True,
    run_visualizations=True,
)
success = pipeline.run()
```

### Legacy Pipeline

```python
from post_docking_analysis.pipeline import PostDockingAnalysisPipeline

pipeline = PostDockingAnalysisPipeline(
    input_dir="/path/to/docking/results",
    output_dir="/path/to/output",
    config_file="/path/to/config.yaml",
)
success = pipeline.run_pipeline()
```

## Troubleshooting

1. Empty SDF files
- The simplified pipeline skips empty SDFs automatically and records them in `skipped_empty_sdf_files.txt`.

2. PandaMap unavailable
- Ensure PandaMap is installed in the configured conda env (default env name: `pandamap`).

3. LigPlot skipped
- Set `LIGPLUS_ROOT`/`LIGPLUS_HOME` or pass `--ligplus-root`.

4. py3Dmol / ProLIF missing
- Install in active environment; pipeline will skip those stages gracefully if unavailable.

5. Matplotlib cache warning
- Set writable `MPLCONFIGDIR` in your environment.
