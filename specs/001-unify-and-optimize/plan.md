# Implementation Plan: Unified Visualization Overhaul

**Branch**: `001-unify-and-optimize` | **Date**: 2026-02-11 | **Spec**: `specs/001-unify-and-optimize/spec.md`
**Input**: Feature specification from `specs/001-unify-and-optimize/spec.md`

## Summary

Implement a visualization-focused refactor for post-docking analysis with five outcomes: protein display-name mapping, PandaMap shape-oriented tuning, ProLIF/LigPlot tuning, corrected RMSD figures, and output consolidation into `visualizations/` + `raw_data/` with manifest linkage.

## Technical Context

**Language/Version**: Python 3.9+  
**Primary Dependencies**: pandas, numpy, matplotlib, seaborn, RDKit, ProLIF, PandaMap, LigPlot+  
**Storage**: File-based pipeline outputs (CSV/JSON/TXT/PNG/SVG/PDF/HTML/PSE)  
**Testing**: Script-level smoke runs + syntax checks  
**Target Platform**: macOS + Conda environments (`pdb-prepare-wizard`, `pandamap`)  
**Project Type**: Single Python package  
**Performance Goals**: No regression in end-to-end run stability; faster downstream consumption via organized outputs  
**Constraints**: Optional external tools may be missing; pipeline should degrade gracefully  
**Scale/Scope**: GNINA post-docking result sets with tens to hundreds of complexes

## Constitution Check

- Keep backwards compatibility for existing output paths where possible.
- Avoid destructive output migration; consolidate by copying/organizing with deterministic structure.
- Preserve optional-tool graceful behavior.

## Project Structure

### Documentation (this feature)

```text
specs/001-unify-and-optimize/
├── spec.md
└── plan.md
```

### Source Code (repository root)

```text
post_docking_analysis/
├── simplified_pipeline.py
├── hierarchical_analyzer.py
├── enhanced_rmsd_analyzer.py
├── publication_pandamap.py
├── prolif_interaction_maps.py
├── ligplot_integration.py
└── [new] protein_naming.py
```

**Structure Decision**: Implement core behavior in `simplified_pipeline.py` and reusable mapping helpers in a new `protein_naming.py`; patch visualization generators to consume display names and tuned parameters.

## Phases

1. Add reusable protein naming utility + mapping persistence.
2. Wire display names into hierarchical/simplified visualization stages.
3. Tune PandaMap, ProLIF, LigPlot defaults and expose control knobs.
4. Correct RMSD plotting ergonomics and edge-case handling.
5. Add output consolidation (`visualizations/`, `raw_data/`) and linkage manifest.
6. Run syntax and smoke validation on modified modules.
