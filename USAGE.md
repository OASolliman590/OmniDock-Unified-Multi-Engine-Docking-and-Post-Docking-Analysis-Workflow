# CLI recipes

These examples use options shown by the live `python main.py ... --help`
surfaces. Run the relevant help command in your checkout before extending them.

## Initialize and inspect a project

```bash
python main.py workflow init --project-dir docking_project --engines gnina,vina,smina,autodock4 --layout-profile canonical
python main.py workflow status --project-dir docking_project
```

Canonical recipes in this guide always pass the layout explicitly. Omitting
`--layout-profile` selects the `docking_legacy` compatibility default.

The guided shell is `python main.py workflow interactive`. Resume an existing
workflow with:

```bash
python main.py workflow resume --project-dir docking_project
```

## Fetch or prepare structures

```bash
python main.py pdb run -p 7CMD -o pdb_results
python main.py pdb prepare-both \
  --receptors-input proteins_raw --ligands-input ligands_raw \
  --receptors-output prepared_proteins --ligands-output prepared_ligands \
  --ligand-backend engine_aware_full \
  --selected-engines gnina,vina,smina,autodock4
```

Inspect the preparation/QC records and read
[PDB_PREPARATION_USAGE.md](PDB_PREPARATION_USAGE.md) before docking.

## Stage a project and pairlist

```bash
python main.py prep project \
  --prepared-proteins prepared_proteins \
  --prepared-ligands prepared_ligands \
  --output docking_project \
  --engines gnina,vina,smina,autodock4 \
  --layout-profile canonical \
  --pair-mode protein-based

python main.py prep pairlist \
  --project-dir docking_project \
  --mode curated_cartesian \
  --curated-receptors receptor_a \
  --curated-ligands ligand_a,ligand_b \
  --freeze
```

Review `project_manifest.json`, `pairlist.csv`, sites, boxes, aliases, and pair
intent. Freezing materializes the selected curation round; it does not validate
the scientific choices.

## Dry-run, then dock locally

```bash
python main.py dock dry-run --project-dir docking_project --engines gnina,vina
python main.py dock run --project-dir docking_project --engines gnina,vina --skip-completed
```

Use `dock run --dry-run` or the dedicated `dock dry-run` verb before consuming
compute. Review generated commands, executable identities, QC, box geometry,
seeds, resources, and output paths.

## Generate, sync, and submit a deployment

```bash
python main.py dock deploy \
  --project-dir docking_project --engines gnina,vina \
  --mode screen --hpc-profile-file local-profile.json

python main.py dock sync \
  --project-dir docking_project --hpc-profile-file local-profile.json --dry-run

python main.py dock submit \
  --project-dir docking_project --mode screen \
  --hpc-profile-file local-profile.json --dry-run
```

Remove each `--dry-run` only after reviewing its plan. Use the public-safe
templates and [HPC_DEPLOYMENT_GUIDE.md](HPC_DEPLOYMENT_GUIDE.md); do not commit
filled profiles.

## Analyze

```bash
python main.py analyze comparative \
  --project-dir docking_project \
  --analysis-scope full \
  --normalization-method per_engine_rank

python main.py analyze favorite-engine \
  --project-dir docking_project \
  --favorite-engine gnina \
  --analysis-scope full

python main.py analyze stage rmsd --project-dir docking_project
python main.py analyze clean --project-dir docking_project --dry-run
```

`comparative` is for compatible multi-engine evidence;
`favorite-engine` continues one selected engine. Stage commands include
hierarchical, polypharmacology, RMSD, reports, visualizations, and
structure-quality. Interaction commands include the clean PLIP/ProLIF path and
optional specialized tools. See [POST_DOCKING_ANALYSIS_GUIDE.md](POST_DOCKING_ANALYSIS_GUIDE.md).

## Safe exhaustive promotion

Create a reviewed rerun manifest from comparative analysis before exhaustive
deployment:

```bash
python main.py analyze comparative \
  --project-dir docking_project \
  --promote-exhaustive --rerun-engine gnina --top-per-protein 3

python main.py dock deploy \
  --project-dir docking_project --mode exhaustive \
  --from-rerun-manifest path/to/rerun_manifest.csv \
  --hpc-profile-file local-profile.json
```

The CLI requires an explicit `--allow-full-exhaustive` override to bypass the
manifest safeguard. Treat that override as a deliberate, reviewed exception.
