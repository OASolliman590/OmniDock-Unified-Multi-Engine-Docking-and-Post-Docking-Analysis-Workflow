# Receptor and ligand preparation

Preparation creates engine inputs plus the provenance and QC evidence needed by
docking and redocking analysis. It does not establish that a target model or
protonation protocol is scientifically appropriate.

## Commands

Prepare both asset classes:

```bash
python main.py pdb prepare-both \
  --receptors-input proteins_raw \
  --ligands-input ligands_raw \
  --receptors-output prepared_proteins \
  --ligands-output prepared_ligands \
  --ph 7.4 \
  --ligand-backend engine_aware_full \
  --selected-engines gnina,vina,smina,autodock4
```

Use `prepare-protein` or `prepare-ligand` for one side. Available ligand backend
profiles are listed by `python main.py pdb prepare-both --help`; some require
Open Babel or explicit AutoDockTools paths.

## Authoritative chemistry and mappings

Ligand preparation requires an authoritative chemical graph from a source such
as a verified SDF/SMILES record. Coordinates alone—especially PDB/PDBQT
coordinates—do not prove bond orders, protonation, stereochemistry, or atom
identity. Preserve source identifiers, input hashes, graph/atom mappings, and
the prepared-output sidecars.

For a native redocking reference, the graph and heavy-atom coordinates must
match the selected experimental ligand instance. Idealized or newly embedded
conformers are not native references. Instance identity includes chain, residue
number, insertion code, and residue name. Insertion-code cases may require an
explicit verified SDF.

## pH and backend choices

`--ph` selects the requested protonation pH, but adding hydrogens is not itself
pH normalization. Record the backend, version, pH, tautomer/protomer choices,
charges, and warnings. A selected normalization/backend failure is reported;
the workflow must not silently substitute a chemically different procedure.

The engine-aware profile produces only representations supported by selected
engines. Macrocycle closure/glue atom types are preserved. An engine that cannot
represent them rejects the input instead of silently deleting or retyping them.

## Fail-closed QC

Preparation reports invalid/nonfinite coordinates and unsupported chemistry.
Receptor preparation preserves heavy-atom coordinates and surfaces backend
errors rather than retrying a permissive residue-deleting path. Docking QC is a
gate by default; disabling it is an explicit protocol change that must be
justified and recorded.

Review at least:

- biological assembly, chains, alternate locations, missing atoms/residues;
- retained waters, metals, cofactors, covalent components, and selected ligand;
- protonation/tautomer state, formal charge, stereochemistry, and atom mapping;
- receptor and ligand coordinate frames;
- box center/size and the intended site;
- per-engine exported formats, warnings, rejected structures, and sidecars.

Pocket volume, electrostatic potential, or druggability is unavailable unless a
validated estimator actually ran. Missing values are not evidence of a good
pocket.

## Stage the prepared project

```bash
python main.py prep project \
  --prepared-proteins prepared_proteins \
  --prepared-ligands prepared_ligands \
  --output docking_project \
  --engines gnina,vina,smina,autodock4 \
  --layout-profile canonical \
  --pair-mode protein-based
```

Omitting the layout flag selects the `docking_legacy` compatibility default.

Then review and, if needed, curate/freeze `pairlist.csv` with `prep pairlist`.
See [PDBQT preparation](PDBQT_PREPARATION_EXPLAINED.md),
[CLI recipes](USAGE.md), and [correctness behavior](docs/correctness-update.md).
