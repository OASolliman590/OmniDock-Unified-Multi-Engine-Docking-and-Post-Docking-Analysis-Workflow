# PDBQT preparation and atom identity

PDBQT is an engine input format, not an authoritative chemical graph. It adds
partial charges and docking atom types to a coordinate representation, but a
legacy or coordinate-only PDBQT may not preserve enough bond-order,
stereochemistry, macrocycle, or atom-identity information for scientific
geometry comparisons.

## Data contract

Begin with a verified chemical graph (for example SDF/SMILES plus provenance),
then export an engine-aware PDBQT while retaining a reversible atom map. Meeko
SMILES/IDX metadata is an accepted mapping source. Preserve:

- original and prepared structure hashes;
- canonical atom identifiers and output atom order;
- bond orders, formal charges, stereochemistry, protonation, and tautomer state;
- preparation tool/version/options and pH;
- macrocycle closure/glue metadata;
- receptor/reference coordinate-frame provenance.

AutoDockTools paths can be supplied through the `pdb prepare-*` options when the
selected backend uses them. Do not assume output from one preparation tool is
interchangeable with another.

## Engine-aware outputs

GNINA, Vina, and Smina use PDBQT-like ligand/receptor inputs but may differ in
supported atom types and runtime requirements. AutoDock4 additionally requires
AutoGrid4/grid maps and a consistent parameter file. Unsupported macrocycle or
atom-type features must fail rather than be erased to make a file parse.

## Geometry and redocking

RMSD requires graph identity/symmetry and an authoritative atom correspondence.
It is measured after placing predicted and reference ligands in their common
native receptor frame; fitting the ligand onto the reference would hide spatial
displacement and is not used. Multi-pose records must select the requested pose
exactly, without substituting pose one, merging models, or truncating atoms.

SDF/MOL and mapped Meeko macrocycles can be evaluable. Coordinate-only
PDB/PDBQT and legacy AD4 DLG geometry remain `not_evaluable` when the required
graph/mapping is absent. That outcome must reduce coverage; it is not a failed
RMSD threshold and not a pass.

## User verification

Inspect the prepared molecule in a chemically aware viewer, compare heavy-atom
coordinates and graph/charges to the authoritative input, confirm expected
rotatable bonds and macrocycle treatment, and run an engine dry-run. Retain all
sidecars with the project. No format conversion alone establishes protonation
accuracy, pose validity, or affinity prediction quality.

See [preparation usage](PDB_PREPARATION_USAGE.md) and
[post-docking analysis](POST_DOCKING_ANALYSIS_GUIDE.md).
