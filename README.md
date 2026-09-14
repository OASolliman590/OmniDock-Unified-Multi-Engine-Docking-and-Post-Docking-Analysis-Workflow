# Omni-DockForge: End-to-End Docking, Consensus by Design.

Formerly **PDB Prepare Wizard**.

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Conda](https://img.shields.io/badge/conda-forge-blue.svg)](https://conda-forge.org/)
[![Biopython](https://img.shields.io/badge/Biopython-1.79+-green.svg)](https://biopython.org/)

A comprehensive platform for multi-engine docking project preparation, execution, and post-docking analysis. DockForge provides automated protein/ligand workflows, comparative docking analytics, and optional interaction-analysis stages.

The scientific correctness update changes preparation, resume, score ranking, and
reference validation. Read [the migration and validation notes](docs/correctness-update.md)
before reusing old results. Consensus ranks are relative prioritization signals;
they are not experimental binding affinities or evidence of biological efficacy.

## DockForge Workflow Map

DockForge is organized as a reusable project lifecycle rather than a one-shot script:

1. Project setup and `workflow init`
2. Receptor/ligand preparation (engine-aware)
3. Pairlist and docking-folder materialization
4. Docking execution (local or deployment-target aware)
5. Post-docking analysis (full/comparison/QC/rescoring/report scopes)
6. Checkpoint & Revise (freeze current progress, then continue from a controlled state)

Interactive and CLI modes share the same state backbone (`.workflow/state.json`), so you can move between guided and scripted operation without losing workflow context.

### Canonical Output Topology (Numbered)

DockForge supports legacy compatibility, but the canonical auditable layout is:

- `0-Input/`
- `1-Preparation/`
- `2-GridBoxes/`
- `3-Docking/`
- `4-Working/`
- `5-Analysis/`
- `6-Visualizations/`
- `7-Reports/`
- `.meta/`
- `sessions/`

`5-Analysis/START_HERE.md` is generated as the main post-docking navigation index, and `7-Reports/START_HERE.md` remains the consolidated cross-run report index.

## 🚀 Features

### Core Functionality
- **PDB Download**: Automatically fetch PDB files from RCSB PDB database
- **Ligand Extraction**: Identify and extract HETATM residues (ligands) from structures
  - **Grouped Display**: HETATMs grouped by type with counts
  - **Smart Selection**: Choose by type, then specific instance
- **Chain Management**: Filter and select specific protein chains
- **Structure Cleaning**: Remove water molecules, ions, and other unwanted residues
  - **Smart Cleaning Options**: Remove specific residues, ALL HETATMs, or common residues only
  - **Removal Summary**: Shows exactly what will be removed before confirmation

### Enhanced Analysis (v3.0)
- **🆕 PLIP Integration**: Optional protein-ligand interaction analysis
  - **Text-based PLIP parsing**: Reliable interaction detection using PLIP's official report format
  - **Comprehensive interaction types**: Hydrophobic, hydrogen bonds, halogen bonds, π-stacking, water bridges, salt bridges, metal complexes, π-cation interactions
  - Results depend on the installed PLIP version and structure preparation
  - **Future-proof design**: Automatically detects unknown interaction types
- **🆕 Residue-Level Coordinate Extraction**: Enhanced `extract_residue_level_coordinates()` function
  - Individual residue averages for detailed binding site analysis
  - Overall binding site center calculation
  - Comprehensive statistics (residue count, atom count)
  - PLIP-enhanced interaction detection with detailed breakdowns
- **Active Site Analysis**: Extract binding site coordinates using PLIP-enhanced or distance-based methods
- **Pocket Analysis**: Contact and residue descriptors from the selected structure.
  Pocket volume, electrostatic potential, and druggability are reported as
  unevaluated where no validated calculation is implemented.

### Multiple Interface Options (v2.1.0)
- **🆕 Interactive Mode**: User-friendly interactive pipeline with guided prompts
- **🆕 CLI Mode**: Command-line interface with configuration file support
- **🆕 Batch Processing**: Enhanced batch processing with configuration-driven automation
- **🆕 Unified Entry Point**: Single `main.py` entry point for all modes
- **🆕 Docking Project Builder**: Build canonical multi-engine docking projects from prepared receptors, prepared ligands, Excel site metadata, and pair-intent rules
- **🆕 Pairlist Automation Modes**: Support manual pairlists, protein-based pair generation, and optional all-proteins-to-all-ligands expansion

### Reporting & Analysis
- **Multi-PDB Analysis**: Analyze multiple PDB structures in a single session
- **Excel Integration**: Generate comprehensive Excel reports with all results
- **Report Generation**: Generate detailed CSV and Excel reports with all analysis results
- **🆕 Post-Docking Analysis**: Comprehensive analysis of molecular docking results
  - **Docking Score Analysis**: Parse engine scores with explicit score directions
  - **Best Pose Selection**: Rank poses using the selected score or consensus policy
  - **PDB Extraction**: Extract best poses as complete receptor-ligand complex PDB files
  - **Statistical Analysis**: Generate comprehensive statistics and rankings
  - **Visualization**: Create binding affinity distributions and top performer plots
  - **Multi-format Reports**: CSV, Excel, and text summary reports
- **🆕 Multi-Engine Docking Execution**: Run GNINA via container and Vina/Smina via conda or direct binaries from one canonical project
- **🆕 Engine-Aware Analysis**: Compare GNINA, Vina, and Smina together or continue downstream from a selected favorite engine

### 🆕 AutoDock Preparation (v3.0.1)
- **Enhanced AutoDock Preparation**: Comprehensive preparation for AutoDock Vina docking
  - **Multi-format Support**: Handles PDB, SDF, MOL2 for ligands; PDB for receptors
  - **Engine-aware ligand preparation**: Open Babel/Meeko/AutoDockTools orchestration
  - **Flexible Input/Output**: Same-folder or separate-folder scenarios
  - **PDB Ligand Extraction**: Extract ligands directly from PDB files
  - **Quality Control**: Comprehensive validation of prepared files
  - **Configuration Support**: JSON/YAML configuration files
  - **Progress Tracking**: Real-time progress indicators and detailed logging
  - **Error Recovery**: Robust error handling and graceful failure recovery
  - **Chemistry preservation**: Require authoritative ligand graphs and retain atom mappings
  - **Validation status**: Regression-tested code; target-specific scientific benchmarking is still required

## 📋 Requirements

### System Requirements
- Python 3.10 or higher
- Internet connection (for PDB downloads)
- 2GB+ RAM recommended for large structures

### Dependencies
- **Core**: numpy, pandas, biopython
- **Optional (post-docking only)**: plip (interaction analysis plugin path)
- **Visualization**: matplotlib, seaborn
- **Excel Support**: openpyxl (for Excel report generation)
- **Chemical mapping and reference RMSD**: RDKit; Meeko mappings for PDBQT poses
- **External preparation**: Open Babel and the chosen receptor/ligand backend
- **Development**: pytest and build; Jupyter is optional
- **Interactive Workflow**: questionary

## 🛠️ Installation

### Option 1: Using Conda (Recommended)

```bash
# Create and activate conda environment
conda env create -f environment.yml
conda activate pdb-prepare-wizard

# Install the package
pip install -e '.[chemistry]'

# Verify installation
python main.py --help
```

### Option 2: Using pip

```bash
# Create virtual environment
python -m venv pdb-wizard-env
source pdb-wizard-env/bin/activate  # On Windows: pdb-wizard-env\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install the package
pip install -e '.[chemistry]'
```

### Option 3: Manual Installation

```bash
# Install core dependencies
pip install numpy pandas biopython

# Install optional PLIP for post-docking interaction analysis
pip install plip

# Optional downstream visualization extras
pip install py3Dmol prolif

# Install visualization tools
pip install matplotlib seaborn

# Install Excel support
pip install openpyxl

# Install interactive workflow prompts
pip install questionary
```

Post-docking dependency contract:
- Core post-docking workflow: required
- `py3Dmol`: optional
- `ProLIF`: optional
- `LigPlot+`: optional external tool
- `PoseView`: optional network-backed stage

The pipeline must complete when optional stages are unavailable. Missing optional
dependencies are reported as skipped stages, not hard failures.

### 🆕 Recent Updates (v3.0.1)

**🔧 Bug Fixes:**
- **Fixed PDB→SDF→PDBQT Conversion**: Resolved explicit hydrogens requirement for Meeko
- **Enhanced Error Handling**: Improved ligand preparation with proper error recovery
- **Preparation/analysis separation**: preparation no longer depends on interaction-analysis modules

**✅ Test Results:**
- **Ligand Preparation**: 6/6 files successfully converted (100% success rate)
- **Receptor Preparation**: 4/5 files successfully prepared (80% success rate)
- **Quality Control**: All prepared files validated and ready for AutoDock Vina

**🎯 Production Status:**
- **Fully Tested**: Validated with real research data from PDB files
- **Production Ready**: All systems working and ready for molecular docking studies
- **Comprehensive Documentation**: Updated guides and examples available

## 🎯 Usage

### Unified Workflow CLI

`main.py` now has one grouped command surface plus backward-compatible legacy aliases.

```bash
# Guided workflow shell (default)
python main.py
python main.py interactive
python main.py workflow interactive
python main.py workflow init --project-dir docking_project/ --layout-profile docking_legacy
python main.py workflow status --project-dir docking_project/
python main.py workflow jump --project-dir docking_project/ --target analyze.stage.rmsd --engine vina

# PDB functionality
python main.py pdb collect --project-dir docking_project/ -p 7CMD
python main.py pdb fetch -p 7CMD -o results/
python main.py pdb run -p "7CMD,6WX4" -o results/
python main.py pdb batch -c pdb_batch_config.yaml -o batch_results/
python main.py pdb prepare-both \
  --receptors-input receptors_raw/ \
  --ligands-input ligands_raw/ \
  --receptors-output receptors_prep/ \
  --ligands-output ligands_prep/

# Legacy aliases still work
python main.py cli -p 7CMD
python main.py batch -c pdb_batch_config.yaml
python main.py prepare-docking \
  --prepared-proteins receptors_prep/ \
  --prepared-ligands ligands_prep/ \
  --excel multi_pdb_analysis.xlsx \
  --output docking_project/

# Build pairlists explicitly before materialization
python main.py prep pairlist \
  --project-dir docking_project/ \
  --mode cocrystal_plus_all

# Materialize the staged docking project
python main.py prep project \
  --project-dir docking_project/

# Docking
python main.py dock run \
  --project-dir docking_project/ \
  --engines gnina,vina,smina \
  --favorite-engine gnina

# Post-docking analysis
python main.py analyze comparative --project-dir docking_project/
python main.py analyze favorite-engine --project-dir docking_project/ --favorite-engine vina
python main.py analyze stage rmsd --project-dir docking_project/ --engine vina
python main.py analyze interactions pandamap --project-dir docking_project/ --engine smina
python main.py analyze visuals pymol --project-dir docking_project/ --engine gnina

# Installed commands (after pip install -e .)
pdb-prepare-wizard
pdb-wizard-workflow workflow interactive
pdb-wizard-interactive
```

The new workflow shell can resume or inspect state from `.workflow/state.json`, jump directly to specific functions, and continue from docking preparation into docking execution and post-docking analysis without switching entrypoints.

`workflow init --layout-profile docking_legacy` bootstraps a staged docking project with:

- `1-Raw_Ligand/`
- `2-Raw_Protien/`
- `3-Preparation/`
- `4-Docking/`
- `5-Analysis/`

For compatibility, DockForge may create `5-Post_Docking_Analysis -> 5-Analysis`
as a symlink when no prior legacy directory exists.

`1-Raw_Ligand/` is now the single raw-ligand staging directory for mixed input
formats. When the workflow can derive an SDF-normalized copy, it is written
into the same folder beside the source ligand instead of a separate
`1-Raw_Ligand_SDF/` tree.

Within that layout, GNINA-compatible paths such as `4-Docking/gnina_out`,
`4-Docking/logs`, `4-Docking/results`, `4-Docking/scripts`, and
`4-Docking/scripts_hpc` are created alongside the Vina and Smina output roots.

### Batch Processing Configuration

The batch processing mode supports two configuration file formats:

#### YAML Configuration Format

Create a YAML configuration file (e.g., `pdb_batch_config.yaml`):

```yaml
# Global settings applied to all PDB entries unless overridden
global_settings:
  output_directory: "batch_docking_preparation"
  
  # Default cleaning strategy
  cleaning:
    strategy: "common"  # Options: "all", "common", "none", or list of residues
    common_residues: ["HOH", "NA", "CL", "SO4", "CA", "MG", "ZN", "FE", "CU", "MN"]
  
  # Active site analysis settings
  analysis:
    use_enhanced_coordinates: true
    distance_cutoff: 5.0
    method: "plip"  # Options: "plip", "distance"
  
  # AutoDock preparation settings
  autodock:
    force_field: "AMBER"
    ph: 7.4
    allow_bad_res: true
    default_altloc: "A"

# List of PDB entries to process
pdb_entries:
  - pdb_id: "7CMD"
    ligand_selection:
      # Specific ligand to extract, if not provided the first non-common ligand is used
      ligand_name: "TTT"
      # chain_id: "A"  # Optional: specify chain
      # res_id: 1     # Optional: specify residue ID
    cleaning:
      # Override global cleaning settings for this PDB
      strategy: ["HOH", "NA", "CL"]
    analysis:
      # Override global analysis settings for this PDB
      use_enhanced_coordinates: true
    autodock:
      # Override global AutoDock settings for this PDB
      prepare_as: "ligand"  # Options: "ligand", "receptor", "both"
  
  - pdb_id: "6WX4"
    ligand_selection:
      ligand_name: "LIG"
    cleaning:
      strategy: "common"
    analysis:
      use_enhanced_coordinates: true
    autodock:
      prepare_as: "receptor"

# Processing options
processing:
  # Continue on error or stop
  continue_on_error: true
```

#### TXT Configuration Format

Create a TXT configuration file (e.g., `pdb_batch_config.txt`):

```txt
# PDB Batch Processing Configuration (TXT Format)
# This file defines multiple PDB entries for batch processing
# Format: PDBID|LIGAND|CLEANING_STRATEGY|PREPARE_AS
# 
# PDBID: 4-character PDB identifier
# LIGAND: Specific ligand to extract (or AUTO for auto-selection)
# CLEANING_STRATEGY: all, common, none, or comma-separated list
# PREPARE_AS: ligand, receptor, both

7CMD|TTT|common|ligand
6WX4|LIG|common|receptor
1ABC|AUTO|none|both
```

### Batch Processing Usage

```bash
# Process PDBs using YAML configuration
python main.py batch -c pdb_batch_config.yaml

# Process PDBs using TXT configuration
python main.py batch -c pdb_batch_config.txt -o batch_results/

# Direct module usage
python batch_pdb_preparation.py -c config.yaml -o output_dir
```

### Python API Usage

```python
from core_pipeline import MolecularDockingPipeline, extract_residue_level_coordinates

# Initialize pipeline
pipeline = MolecularDockingPipeline(output_dir="my_analysis")

# Download and process PDB
pdb_file = pipeline.fetch_pdb("1ABC")
hetatm_details, unique_hetatms = pipeline.enumerate_hetatms(pdb_file)

# Extract specific ligand
ligand_pdb = pipeline.save_hetatm_as_pdb(pdb_file, "LIG", "A", 1)

# Clean structure
cleaned_pdb = pipeline.clean_pdb(pdb_file, to_remove_list=['HOH', 'NA', 'CL'])

# Analyze binding site with enhanced coordinates
residue_analysis = extract_residue_level_coordinates(cleaned_pdb, "LIG", "A", 1)
if residue_analysis:
    coords = residue_analysis['overall_center']
    pocket_results = pipeline.analyze_pocket_properties(cleaned_pdb, coords)
```

### CLI Configuration (v2.1.0)

Create a configuration file for automated processing:

```bash
# Generate sample config
python cli_pipeline.py --create-config

# Use configuration file
python main.py cli -p "7CMD,6WX4" -c config.json
```

Example configuration:
```json
{
  "description": "Sample configuration for Omni-DockForge legacy PDB CLI",
  "ligand_selection": {
    "7cmd": "TTT",
    "6wx4": "LIG",
    "1abc": "GDP"
  },
  "cleaning": {
    "default": "common",
    "7cmd": ["HOH", "NA", "CL"],
    "custom_pdb": "all"
  },
  "analysis": {
    "use_enhanced_coordinates": true,
    "distance_cutoff": 5.0
  }
}
```

### Advanced Usage - Residue-Level Analysis (v2.1.0)

The new version includes enhanced coordinate extraction at the residue level:

```python
from core_pipeline import MolecularDockingPipeline, extract_residue_level_coordinates

# Initialize pipeline
pipeline = MolecularDockingPipeline()

# Process PDB
pdb_file = pipeline.fetch_pdb("1ABC")
hetatm_details, unique_hetatms = pipeline.enumerate_hetatms(pdb_file)
cleaned_pdb = pipeline.clean_pdb(pdb_file, to_remove_list=['HOH', 'NA', 'CL'])

# Enhanced residue-level coordinate extraction
residue_analysis = extract_residue_level_coordinates(cleaned_pdb, "LIG", "A", 1)

if residue_analysis:
    print(f"Binding site center: {residue_analysis['overall_center']}")
    print(f"Interacting residues: {residue_analysis['num_interacting_residues']}")
    print(f"Total atoms: {residue_analysis['num_interacting_atoms']}")
    
    # Individual residue averages
    for residue, coord in residue_analysis['residue_averages'].items():
        print(f"  {residue}: {coord}")
    
    # Analyze pocket properties
    pocket_results = pipeline.analyze_pocket_properties(
        cleaned_pdb, residue_analysis['overall_center']
    )
```

### AutoDock Preparation Usage

The enhanced AutoDock preparation system is profile-driven and engine-aware.
Preparation generates docking-ready assets and does not run PLIP.

#### Command Line Usage

# Create configuration file
python autodock_preparation.py --create-config

# Edit configuration file
nano autodock_config.json

# Run preparation with custom settings
python autodock_preparation.py \
    --ligands-input ./ligands_raw \
    --receptors-input ./receptors_raw \
    --ligands-output ./ligands_prep \
    --receptors-output ./receptors_prep \
    --force-field AMBER \
    --ph 7.4 \
    --ligand-profile engine_aware_full \
    --selected-engines gnina,vina,smina,autodock4

# Run with bash script
./prep_autodock_enhanced.sh autodock_config.json

#### Python API Usage

from autodock_preparation import AutoDockPreparationPipeline, PreparationConfig

# Create configuration
config = PreparationConfig(
    ligands_input="./ligands_raw",
    receptors_input="./receptors_raw",
    ligands_output="./ligands_prep",
    receptors_output="./receptors_prep",
    force_field="AMBER",
    ph=7.4,
    ligand_preparation_profile="engine_aware_full",
    selected_engines=["gnina", "vina", "smina", "autodock4"],
)

# Initialize and run pipeline
pipeline = AutoDockPreparationPipeline(config)
success = pipeline.run_enhanced_preparation()

if success:
    results = pipeline.analyze_preparation_results("./receptors_prep")
    print(f"Ligands prepared: {results['ligands']['count']}")
    print(f"Receptors prepared: {results['receptors']['count']}")

#### Expected Input Structure

```
project_directory/
├── ligands_raw/                 # Input ligands (mixed formats; normalized SDF copies can live beside the source file)
│   ├── ligand1.sdf
│   ├── ligand1.mol
│   ├── ligand2.mol2
│   └── ligand3.pdb
├── receptors_raw/               # Input receptors (PDB)
│   ├── receptor1.pdb
│   └── receptor2.pdb
└── autodock_config.json         # Configuration file
```

#### Output Structure

```
project_directory/
├── ligands_prep/                # Prepared ligands (PDBQT)
│   ├── ligand1.pdbqt
│   ├── ligand2.pdbqt
│   └── ligand3.pdbqt
├── receptors_prep/              # Prepared receptors (PDBQT)
│   ├── receptor1.pdbqt
│   ├── receptor2.pdbqt
│   └── preparation_summary.txt  # Summary report
└── logs/                        # Log files
```

Ligand preparation normalizes each ligand to an explicit-hydrogen 3D SDF
intermediate before generating PDBQT output. Existing 3D coordinates are kept
when present; 2D inputs such as `.mol` or flat `.sdf` ligands are embedded to
3D first so the Meeko/Vina-family preparation path receives chemically usable
coordinates.
When `autodock4` is selected, the engine-aware path includes AutoDockTools-compatible branching automatically.

### Post-Docking Analysis Usage

Two workflows are available:

1. Recommended for GNINA projects: `post_docking_analysis.simplified_cli`
2. Legacy/config-driven workflow: `post_docking_analysis` (full pipeline)

#### Recommended GNINA Workflow

```bash
# Auto-detects local vs HPC GNINA layout (gnina_out + logs + receptors)
python -m post_docking_analysis.simplified_cli \
  --project-dir /path/to/GNINA_project \
  --output /path/to/post_docking_output

# Explicit folder mode
python -m post_docking_analysis.simplified_cli \
  --sdf-folder /path/to/gnina_out \
  --log-folder /path/to/logs \
  --receptors-folder /path/to/receptors \
  --output /path/to/post_docking_output \
  --pairlist /path/to/pairlist.csv

# Optional runtime switches
# --no-rmsd
# --no-visualizations
# --ligplus-root /path/to/LigPlus
```

#### Legacy/Config Workflow

```bash
python -m post_docking_analysis -i /path/to/docking/results -o /path/to/output
python -m post_docking_analysis --config my_config.yaml -i /path/to/docking/results
```

#### Canonical Multi-Engine Workflow

```bash
# Compare all engines in one canonical project
python -m post_docking_analysis \
  --project-dir /path/to/docking_project \
  --analysis-mode comparative_all_engines \
  -o /path/to/analysis_output

# Continue downstream from one selected engine
python -m post_docking_analysis \
  --project-dir /path/to/docking_project \
  --analysis-mode favorite_engine_continue \
  --favorite-engine vina \
  -o /path/to/analysis_output
```

`favorite_engine_continue` behavior:
- `gnina`: continues into the existing simplified structural workflow
- `vina` / `smina`: continue into a unified structural bridge that reuses the simplified downstream stack where safe. It extracts best-pose complex PDBs, runs hierarchical affinity analysis, polypharmacology, general visualizations, PandaMap, report generation, structural quality, output consolidation, and PDB-based RMSD analysis under `rmsd_analysis/`. Per-protein and global best-pose RMSD always run when comparable complexes exist, and per-complex RMSD is added when multiple poses exist for the same tag. Bridge notes now report stage states explicitly such as `completed`, `disabled`, `missing_dependency`, and `missing_configuration`.

#### Simplified Pipeline Input Layouts

```text
# Local layout
project/
├── gnina_out/        # SDF + log files
├── receptors/
└── pairlist.csv      # optional but recommended

# HPC layout
project/
├── gnina_out/        # SDF files
├── logs/             # log files
├── receptors/
└── pairlist.csv      # optional but recommended
```

#### Visualization Kits (Simplified Pipeline)

1. Affinity overview plots (`matplotlib` / `seaborn`)
2. Hierarchical analysis visualizations
3. RMSD clustering/diversity visualizations
4. PandaMap publication-quality 2D/3D maps
5. py3Dmol interactive 3D HTML views
6. ProLIF interaction maps
7. LigPlot+ interaction diagrams

#### Output Structure (Current)

```text
post_docking_output/
├── analysis/                         # hierarchical analysis outputs
├── complexes/                        # receptor+ligand PDB complexes
├── best_poses/                       # categorized strongest poses
├── reports/                          # tabular summaries
├── rmsd_analysis/                    # scoped RMSD outputs
│   ├── per_complex_all_poses/
│   ├── per_protein_best_poses/
│   └── global_best_poses/
├── interactions/                     # canonical interaction outputs
│   ├── pandamap/
│   ├── prolif/
│   ├── ligplot/
│   └── poseview/
├── 3d_visualizations/                # native py3Dmol outputs
├── visualizations/                   # consolidated visual assets by type
└── raw_data/                         # consolidated raw artifacts + manifest
    ├── visualization_manifest.csv
    └── visualization_manifest.json
```

## 📊 Output Files

The pipeline generates several output files:

- `{PDB_ID}.pdb`: Original downloaded PDB file
- `ligand_{LIGAND}_{CHAIN}_{RESID}.pdb`: Extracted ligand structure
- `cleaned.pdb`: Structure with unwanted residues removed
- `chain_filtered.pdb`: Structure with selected chains only
- `{PDB_ID}_pipeline_results.csv`: Individual PDB analysis report
- `multi_pdb_analysis.xlsx`: **Comprehensive Excel report** (when analyzing multiple PDBs)

## 🔬 Analysis Methods

### Active Site Detection
- **Distance-based**: Uses 5Å cutoff to identify interacting residues
- **PLIP-based**: Uses PLIP library for advanced interaction analysis (if available)

### Pocket Analysis
- **Electrostatic Analysis**: Calculates charge distribution around binding site
- **Hydrophobic Analysis**: Identifies hydrophobic residues within 8Å of pocket center
- **Volume Estimation**: Estimates pocket volume based on interaction sphere
- **Druggability Scoring**: Combines multiple factors to predict druggability (0-1 scale)

## 📈 Example Results

```
📋 Pipeline Summary:
----------------------------------------
PDB_ID                    : 1ABC
Selected_Ligand           : LIG_A_1
Active_Site_Center_X      : 12.345
Active_Site_Center_Y      : 23.456
Active_Site_Center_Z      : 34.567
Interacting_Atoms_Count   : 156
Pocket_Volume_A3          : 523.6
Electrostatic_Score       : 2.5
Nearby_Charged_Residues   : 8
Hydrophobic_Score         : 6
Nearby_Hydrophobic_Residues: 12
Druggability_Score        : 0.85
```

## 🧪 Testing

### Test with Example PDB

```python
# Test with a well-known structure
pipeline = MolecularDockingPipeline()
results = pipeline.run_complete_pipeline("1ABC", interactive=False)
```

### Validation

The pipeline includes several validation steps:
- PDB ID format validation
- File existence checks
- Dependency availability checks
- Error handling for missing ligands

## 🔧 Troubleshooting

### Common Issues

1. **PLIP Installation Issues**
   ```bash
   # Try installing from conda-forge
   conda install -c conda-forge plip
   ```

2. **PDB Download Failures**
   - Check internet connection
   - Verify PDB ID exists in database
   - Try alternative PDB IDs

3. **Memory Issues**
   - Use smaller structures for testing
   - Close other applications
   - Consider using cloud resources for large structures

### Error Messages

- `❌ Error importing Biopython`: Install biopython with `pip install biopython`
- `❌ Missing dependency`: Install missing package from requirements.txt
- `❌ Could not find HETATM`: Structure may not contain ligands

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for details.

### Quick Start for Contributors

```bash
# Fork and clone the repository
git clone https://github.com/yourusername/pdb-prepare-wizard.git
cd pdb-prepare-wizard

# Set up development environment
conda env create -f environment.yml
conda activate pdb-prepare-wizard
pip install -e .

# Run tests
python test_pipeline.py
```

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Biopython team for the excellent structural biology library
- PLIP developers for protein-ligand interaction analysis
- RCSB PDB for providing structural data

## 📞 Support

For issues and questions:
1. Check the troubleshooting section
2. Review error messages carefully
3. Try with different PDB structures
4. Open an issue with detailed error information

## 📚 Documentation

- [Installation Guide](INSTALLATION_GUIDE.md) - Complete setup and installation instructions
- [PDBQT Preparation Guide](PDBQT_PREPARATION_EXPLAINED.md) - Canonical AutoDock/PDBQT preparation guide
- [Post-Docking Analysis Guide](POST_DOCKING_ANALYSIS_GUIDE.md) - Complete guide for analyzing docking results
- [PLIP Interaction Types Guide](PLIP_INTERACTION_TYPES_GUIDE.md) - Understanding PLIP interaction analysis
- [Changelog](CHANGELOG.md) - Version history and changes

## 📁 Project Structure

```
pdb-prepare-wizard/
├── core_pipeline.py              # Core PDB processing pipeline
├── main.py                       # Main entry point (interactive/CLI/batch modes)
├── interactive_pipeline.py       # Interactive user interface
├── cli_pipeline.py              # Command-line interface
├── batch_pdb_preparation.py     # Batch processing module
├── autodock_preparation.py      # AutoDock Vina preparation
├── post_docking_analysis/        # Post-docking analysis module
│   ├── pipeline.py              # Main analysis pipeline
│   ├── affinity_analyzer.py    # Binding affinity analysis
│   ├── rmsd_analyzer.py         # RMSD and clustering analysis
│   └── ...
├── prep_autodock_enhanced.sh    # Supported AutoDock/PDBQT preparation adapter
├── requirements.txt             # Python dependencies
├── environment.yml              # Conda environment
└── setup.py                     # Package installation
```

---

**Version**: 3.0.1  
**Last Updated**: 2025-01-15  
**Python Compatibility**: 3.8+ 
