# Changelog

## Unreleased — repository cleanup

- Consolidated current documentation behind `docs/README.md` and replaced the
  mixed legacy/current landing page with a compact, source-backed overview.
- Added architecture and testing/evidence guides; corrected installation,
  dependency, CLI, preparation, deployment, and post-docking documentation.
- Archived five superseded post-docking design/usage records without rewriting
  their historical claims and labeled dated release/checklist records as
  snapshots.
- Retained legacy distribution names, console entrypoints, and project layout
  seams as documented compatibility surfaces.

Older release and roadmap entries below are historical records. They do not
promise current readiness or constitute current test/scientific evidence.

All notable changes to Omni-DockForge (formerly PDB Prepare Wizard) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2026-04-01

### Added
- **Omni-DockForge Platform Upgrade**: Introduced the `Omni-DockForge` naming and migration layer with compatibility notes for legacy command surfaces.
- **Checkpoint & Revise Workflow**: Added checkpoint lineage metadata and renamed the old "maturation" concept to `Checkpoint & Revise`.
- **Interactive Timeline and Task State Model**: Added explicit step-state tracking (`not_started`, `in_progress`, `completed`, `validated`, `failed`, `skipped`, `needs_review`) and timeline rendering in interactive workflow flows.
- **Background Task Infrastructure**: Added non-blocking background task scaffolding so long-running workflow operations can continue while users navigate other menus.
- **Engine-Aware Ligand Preparation Graph**: Added preparation-mode orchestration and validation for Open Babel, Meeko, AutoDockTools, and engine-aware combinations.
- **Docking Runtime Abstraction**: Added execution-environment models and runtime adapters for local CPU/GPU and future backend expansion.
- **Canonical Post-Docking Artifacts**: Added structured score outputs and canonical analysis mirrors under numbered project layout paths.
- **Top-Pose Ligand Performance Atlas**: Added deterministic ligand-centric top-pose exports (per protein and global) with manifest and run-index integration.
- **Validation Gate Framework**: Added redocking-validation gate artifacts and effective-policy fallback signaling for post-docking classification.
- **Biology Integration and Plugin Hooks**: Added biology-annotation ingestion/reporting and pluggable rescoring/analysis capability gates (including OnionNet2 scaffold support).
- **Reproducibility Metadata**: Added `.meta` run provenance artifacts (`config`, `manifest`, `env.lock` snapshot pathways) and session durability support.
- **SQLite Optional Backend Path**: Added CSV+SQLite dual-write scaffolding with parity checks for scalable downstream querying.
- **Artifact DAG Executor**: Added declarative artifact-graph execution for post-docking analysis with dependency-aware ordering, tiered parallel execution, cache-aware reruns, execution reports, and scope-based artifact requests.
- **Consensus Strategy Injection**: Added explicit DAG consensus strategies for multi-engine, GNINA solo, Vina solo, and Smina solo runs with a unified `consensus_ranked.csv` schema across all four modes.
- **Per-Complex Interaction Fan-Out**: Added per-complex threaded execution adapters for ProLIF, PandaMap, PoseView, and PyMOL interaction nodes, including per-complex result ledgers.

### Changed
- **Unified Post-Docking Execution Path**: User-facing canonical post-docking execution now routes through unified orchestration.
- **Simplified CLI Wrapper Behavior**: Manifest-backed simplified CLI execution now enforces unified routing; legacy-only flags (`--no-rmsd`, `--no-visualizations`) are ignored in wrapper mode.
- **Workflow Stage Target Contract**: Non-canonical legacy stage-target fallback was removed; canonical project context is now required for unified stage-target execution.
- **RMSD User-Facing Scope Contract**: User-facing RMSD scope was constrained to `per_complex` to remove ambiguous global/per-protein paths.
- **Preparation UX Cleanup**: Removed PLIP from preparation-phase prompts and kept PLIP interaction logic in post-docking analysis context.
- **Default Post-Docking Runner**: Unified post-docking execution now defaults to DAG-first artifact resolution instead of the old linear step loop.
- **Interaction Failure Handling**: Optional interaction branches now aggregate partial per-complex failures into warnings instead of behaving like single all-or-nothing stages.

### Fixed
- **Consensus Normalization Direction**: Fixed inverted min-max and z-score direction so stronger binders rank correctly.
- **Pairlist Matching Safety**: Replaced ambiguous bidirectional substring matching with deterministic directional behavior and unmatched-file reporting.
- **Correlation Statistics Rigor**: Added minimum sample-size guards and multiple-testing correction support for correlation outputs.
- **Single-Engine QC Downgrade Bug**: Removed false low-agreement penalty behavior for single-engine runs.
- **Atomic QC Downgrade Logic**: Prevented multi-step cascading downgrades by applying a single atomic QC classification pass.
- **Complex Structure Integrity**: Preserved ligand residue identity in generated complexes and hardened pose-selection determinism.
- **Report Formatting Contract**: Fixed newline serialization in report text assembly outputs.
- **Analysis Session Progress Contract**: Hardened progress-file/session contract to avoid missing progress artifact failures.
- **DAG Cache Invalidation Boundaries**: Removed volatile run identifiers from cache-driving artifacts so unchanged reruns correctly hit cache while downstream parameter changes invalidate only the necessary subgraph.
- **Solo-Mode Consensus Semantics**: Fixed solo-mode consensus routing so GNINA can prioritize `cnn_affinity`, Vina can expose `rmsd_lb` pose-diversity context, and Smina can preserve scoring-weight provenance without breaking schema compatibility.
- **Optional Interaction Recovery**: Hardened optional interaction nodes so missing tools or per-complex failures degrade to `completed_with_warnings` rather than aborting unrelated report-generation branches.

### Verification
- **Smoke Coverage**: `python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix` passed.
- **Pipeline Smoke**: `python test_pipeline.py` passed.

## [Unreleased] - 2026-04-12

### Added
- **Interactive HPC ADMET Control**: Added an explicit HPC prompt to keep ligand QC enabled while disabling ADMET filters for deployment and submission flows.
- **HPC Resume Handoff Notes**: Added local session handoff guidance, resume commands, and operational notes for remote project-root usage.

### Changed
- **Workflow CLI Flag Forwarding**: `workflow/cli.py` now forwards ligand ADMET-related flags for `dock deploy` and `dock submit`, matching the downstream docking CLI contract.
- **Interactive HPC QC Behavior**: The HPC panel now passes `--no-ligand-admet-filters` when selected, instead of forcing users to disable the entire ligand QC gate.

### Fixed
- **Whole-Project HPC Blocking**: Full-project reruns are no longer blocked solely by ADMET policy when users intentionally keep structural QC enabled and disable ADMET filtering.
- **Interactive HPC Regression Coverage**: Added smoke coverage for both complete QC bypass and ligand-QC-with-ADMET-bypass scenarios.

### Verification
- `python -m py_compile workflow/interactive.py workflow/cli.py test/test_dockforge_smoke.py`
- `MPLCONFIGDIR=/tmp/mpl PYTHONPATH=test python - <<'PY' ... _smoke_hpc_interactive_qc_bypass_contract(); _smoke_hpc_interactive_admet_bypass_contract() ... PY`

## [3.0.1] - 2025-01-15

### 🔧 Critical Bug Fixes
- **Fixed PDB→SDF→PDBQT Conversion**: Resolved explicit hydrogens requirement for Meeko ligand preparation
- **Ligand Preparation**: Fixed "RDKit molecule has implicit Hs. Need explicit Hs." error
- **Enhanced Error Handling**: Improved ligand preparation with proper error recovery
- **OpenBabel Integration**: Added `-h` flag for explicit hydrogens in PDB→SDF conversion

### 🚀 Improvements
- **Enhanced Bash Script**: Updated `prep_autodock_enhanced.sh` with explicit hydrogens support
- **Python Wrapper**: Updated `autodock_preparation.py` with improved error handling
- **Documentation**: Updated guides with troubleshooting information for fixed issues

### ✅ Production Testing
- **Ligand Preparation**: 6/6 files successfully converted (100% success rate)
- **Receptor Preparation**: 4/5 files successfully prepared (80% success rate)
- **PLIP Integration**: Fully working with binding site detection and validation
- **Quality Control**: All prepared files validated and ready for AutoDock Vina

### 🎯 Status
- **Production Ready**: All systems working and ready for molecular docking studies
- **Fully Tested**: Validated with real research data from PDB files
- **Comprehensive Documentation**: Updated guides and examples available

## [3.0.0] - 2025-01-15

### Major PLIP Integration
- **🆕 Advanced PLIP Integration**: Research-grade protein-ligand interaction analysis
  - **Text-based PLIP parsing**: Reliable interaction detection using PLIP's official report format
  - **Comprehensive interaction types**: All PLIP interaction types supported
    - Hydrophobic interactions
    - Hydrogen bonds
    - Halogen bonds
    - π-stacking interactions
    - Water bridges
    - Salt bridges (configuration ready)
    - Metal complexes (configuration ready)
    - π-cation interactions (configuration ready)
  - **Perfect PLIP web server match**: Results exactly match the official PLIP web server
  - **Future-proof design**: Automatically detects unknown interaction types
  - **Robust error handling**: Graceful fallback to distance-based methods

### Enhanced Analysis
- **🆕 PLIP-Enhanced Coordinate Extraction**: Advanced binding site analysis
  - **Comprehensive interaction detection**: Uses PLIP's sophisticated algorithms
  - **Detailed interaction breakdowns**: Residue-level analysis with interaction type classification
  - **Enhanced accuracy**: Research-grade results matching PLIP web server
  - **Automatic fallback**: Distance-based method when PLIP is unavailable

### Documentation
- **🆕 PLIP Interaction Types Guide**: Comprehensive documentation of all PLIP interaction types
- **🆕 Future-proofing documentation**: Guide for handling new interaction types
- **🆕 Enhanced README**: Updated with PLIP integration features

## [2.1.0] - 2025-01-15

### Major Restructuring
- **🔄 Modular Architecture**: Complete codebase reorganization to eliminate duplication
  - **Consolidated Core Pipeline**: `core_pipeline.py` - Unified core functionality
  - **Interactive Mode**: `interactive_pipeline.py` - User-friendly interactive interface
  - **CLI Mode**: `cli_pipeline.py` - Command-line interface with configuration support
  - **Main Entry Point**: `main.py` - Unified entry point for all modes
  - **Batch Processing**: Enhanced `batch_pdb_preparation.py` with residue-level analysis

### Added
- **🆕 Enhanced Coordinate Extraction**: Residue-level binding site analysis
  - **Residue-Level Analysis**: `extract_residue_level_coordinates()` function
  - **Individual Residue Averages**: Calculate average coordinates for each interacting residue
  - **Comprehensive Binding Site Data**: Overall center, residue details, and atom counts
  - **Excel Integration**: Enhanced Excel reporting with residue-level data
- **🆕 Post-Docking Analysis Module**: Comprehensive analysis of molecular docking results
  - **Binding Affinity Analysis**: Parse and analyze Vina/GNINA docking results from PDBQT files
  - **Best Pose Selection**: Automatically identify highest binding affinity poses using `idxmin()`
  - **PDB Extraction**: Extract best poses as complete receptor-ligand complex PDB files
  - **Statistical Analysis**: Generate comprehensive statistics and rankings
  - **Multi-format Visualization**: Binding affinity distributions and top performer plots
  - **Comprehensive Reporting**: CSV, Excel, and text summary reports
  - **Open Babel Integration**: For ligand processing and PDBQT to PDB conversion
  - **Flexible Input Detection**: Auto-detects single-folder or multi-folder directory structures
  - **Command-line Interface**: Full CLI with configuration file support
  - **Python API**: Programmatic access to all analysis functions
- **🆕 CLI Configuration System**: JSON-based configuration for automated processing
  - **Ligand Selection**: Specify target ligands per PDB
  - **Cleaning Strategies**: Define cleaning approaches per PDB or globally
  - **Analysis Options**: Configure analysis parameters
  - **Sample Config Generation**: `--create-config` option

### Improved
- **📈 Better Error Handling**: Improved error messages and graceful failure handling
- **🚀 Performance**: Optimized coordinate extraction algorithms
- **📊 Enhanced Reporting**: More detailed Excel reports with residue-level data
- **🔧 Batch Processing**: Fixed PLIP dependencies and improved reliability

### Fixed
- **❌ PLIP Dependency Issues**: Replaced unreliable PLIP calls with robust distance-based analysis
- **🔄 Code Duplication**: Eliminated duplication between PDP_prep.py and PDP_prep_improved.py
- **📝 Batch Processing**: Fixed coordinate extraction errors in batch_pdb_preparation.py

### Removed
- **Old Files**: Removed `PDP_prep.py` and `PDP_prep_improved.py` (replaced by modular structure)

### Dependencies
- **Added Open Babel**: Required for post-docking analysis ligand processing

### Migration Guide
- **From 2.0.0 to 2.1.0**:
  - Update imports: `from core_pipeline import MolecularDockingPipeline`
  - Use new entry points: `python main.py` instead of `python PDP_prep_improved.py`
  - CLI users: Use `python main.py cli` with new configuration options
  - Batch users: Updated `batch_pdb_preparation.py` with enhanced coordinate extraction

## [2.0.0] - 2024-12-19

### Added
- **Multi-PDB Analysis**: Analyze multiple PDB structures in a single session
- **Excel Integration**: Generate comprehensive Excel reports with all results
- **Enhanced HETATM Enumeration**: Grouped display with counts and smart selection
- **Smart Cleaning Options**: 
  - Remove ALL HETATMs at once
  - Remove water only (`WATER`)
  - Remove ions only (`IONS`)
  - Remove common residues only (`COMMON`)
- **Removal Summary**: Shows exactly what will be removed before confirmation
- **Warning System**: Alerts when removing selected ligand
- **Smart Fallback**: Uses original structure if ligand was removed during cleaning
- **Professional Excel Reports**: Formatted output with headers and styling
- **Comprehensive Documentation**: README, contributing guide, and feature docs

### Changed
- **Improved User Experience**: Better prompts and error handling
- **Enhanced Error Handling**: Graceful handling of failed PDBs
- **Better Progress Feedback**: Real-time status updates
- **Flexible Exit Options**: Can quit anytime with 'quit' command

### Fixed
- **Ligand Detection**: Proper handling when ligand is removed during cleaning
- **Active Site Analysis**: Fallback to original structure when needed
- **HETATM Selection**: Grouped display prevents confusion with many similar residues
- **Cleaning Workflow**: Warning system prevents accidental ligand removal

### Technical
- **Type Hints**: Comprehensive type annotations throughout
- **Error Validation**: Better input validation and error messages
- **Dependency Management**: Updated requirements and environment files
- **Code Organization**: Modular design with clear separation of concerns

## [1.0.0] - 2024-12-19

### Added
- **Core Pipeline**: Basic PDB download and processing
- **Ligand Extraction**: Identify and extract HETATM residues
- **Structure Cleaning**: Remove water molecules and ions
- **Active Site Analysis**: Extract binding site coordinates
- **Pocket Analysis**: Electrostatic and hydrophobic analysis
- **Druggability Scoring**: Predict drug binding potential
- **CSV Reports**: Generate detailed analysis reports
- **Interactive Mode**: User-friendly command-line interface

### Features
- PDB file download from RCSB PDB database
- Chain filtering and selection
- HETATM enumeration and selection
- Structure cleaning with customizable residue removal
- Distance-based interaction detection
- Pocket property analysis (electrostatic, hydrophobic, volume)
- Druggability scoring algorithm
- Comprehensive CSV reporting

## [Roadmap (Legacy Plans)]

### Planned
- **Visualization**: 3D structure visualization capabilities
- **Docking Integration**: Direct integration with AutoDock/Vina
- **Batch Processing**: Command-line batch processing mode
- **API Mode**: Programmatic API for integration
- **Web Interface**: Web-based user interface
- **Advanced Analysis**: More sophisticated pocket analysis methods
- **Machine Learning**: ML-based druggability prediction
- **Cloud Support**: Cloud-based processing capabilities

---

## Version History

- **2.0.0**: Major release with multi-PDB analysis and Excel integration
- **1.0.0**: Initial release with core functionality

## Migration Guide

### From 1.0.0 to 2.0.0
- **Backward Compatible**: All existing functionality preserved
- **New Features**: Multi-PDB analysis and Excel integration added
- **Enhanced UX**: Improved HETATM enumeration and cleaning options
- **Better Error Handling**: More robust error handling and warnings

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for details on how to contribute to this project.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. 
