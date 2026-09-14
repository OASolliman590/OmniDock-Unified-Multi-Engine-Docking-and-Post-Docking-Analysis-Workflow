# Scientific correctness update

This update follows the source audit of commit
`cfe33056c77d817d19439a9e9affe225d00ca488`. The original findings are preserved in
`audit/`; they describe the baseline, not the corrected working tree.

## Changed scientific behavior

### Structure and ligand preparation

- Extraction identifies one ligand instance by chain, residue number, insertion
  code, and residue name, resolves alternate locations consistently, and records
  atom and coordinate provenance. Default cleaning removes that selected
  ligand from the receptor.
- Native reference chemistry is retrieved from the RCSB ligand-instance service
  and accepted only when its graph and heavy-atom coordinates match the selected
  experimental instance. Idealized or newly embedded conformers cannot become
  redocking references. Insertion-code instances currently require an explicit
  verified SDF instead of an automatic service lookup.
- Ligand preparation requires an authoritative graph. Failure of the selected
  pH-normalization/backend step is reported instead of silently substituting
  another chemical procedure. Adding hydrogens alone is not pH normalization.
- Macrocycle closure and glue atom types are retained. An engine that cannot
  represent them rejects the input instead of deleting or retyping atoms.
- Receptor preparation preserves heavy-atom coordinates and reports backend
  errors. It does not silently retry a permissive procedure that deletes
  unsupported residues. Metals, cofactors, missing residues, protonation, and
  biological assembly selection still require target-specific decisions.
- Contact geometry uses actual atom distances and residue identities. Nonfinite
  coordinates/box parameters fail validation. Pocket volume, electrostatic
  potential, and druggability are unavailable when no validated estimator ran.

### Docking execution

- Resuming requires matching input contents, executable identities, protocol,
  and validated outputs. An existing pose filename alone is insufficient.
- New attempts preserve old outputs as history; a failed or empty new attempt
  cannot appear successful through old scores or poses. CLI failures return a
  nonzero exit status.
- Local and deployment jobs use the same execution contract, reviewed parameter
  schema, random seed, and box geometry. Deployment paths are portable.
- AutoDock4 grids round outward, and grid generation and docking use the same
  custom parameter file. GNINA resource options remain independent of seed
  selection. Required QC exceptions stop execution.

### Post-docking analysis

- Energy scores sort ascending; GNINA CNNscore and CNNaffinity sort descending.
  All supported normalizations and consensus modes preserve these directions.
  Consensus score is always higher-is-better.
- Comparisons normalize within engine/scoring function and target/site. Raw
  energy offsets between scoring functions must not decide a consensus winner.
  Reference anchoring requires the same engine and scoring function.
- Redocking RMSD measures displacement in the common receptor frame. It never
  fits the ligand onto the reference to conceal a misplaced pose. RDKit graph
  identity and symmetry determine atom correspondence; mismatched graphs,
  malformed records, absent mappings, and unavailable poses are unevaluable.
- Multi-pose SDF/PDBQT parsing selects the requested record exactly. It cannot
  substitute pose one, merge models, or truncate atoms to make a comparison fit.
- PDBQT geometry requires authoritative atom mappings such as Meeko SMILES/IDX
  metadata. Coordinate-only PDB/PDBQT files cannot establish chemical identity.
  SDF/MOL and mapped Meeko macrocycles are supported; legacy AD4 DLG geometry
  remains unevaluable without an authoritative graph/atom mapping.
- Reference flags are explicitly parsed; blank values and the string `False`
  are not evidence. Reference provenance and matching receptor/reference frame
  identifiers are required. Missing comparisons reduce coverage and cannot
  become successful geometric agreement or high validation confidence.
- Duplicate biological measurement keys require explicit aggregation before
  joining, preventing accidental duplication of docking candidates.

## Reusing existing projects

1. Preserve the original inputs and results before recomputation. The baseline
   reports remain useful historical records, but rankings or validation labels
   affected by these defects need to be regenerated.
2. Prepare ligands from authoritative SDF/SMILES chemistry and retain preparation
   sidecars. For redocking, provide the experimental reference structure,
   accession/source, atom mapping, and common coordinate-frame provenance.
3. Review receptor composition, protonation conditions, target assembly, and box
   geometry. Unsupported chemistry must be resolved explicitly.
4. Rerun preparation/docking as necessary. Legacy output files without the new
   completion evidence are not automatically certified as resumable.
5. Regenerate analysis. Do not interpret `not_evaluable` or missing validation
   as a pass. Inspect comparison coverage and reasons alongside rankings.

For a deliberate score-only import, place `normalized_scores.import.json`
beside `scores/normalized_scores.csv` with this content (replace the digest):

```json
{"kind": "imported_scores", "sha256": "SHA256_OF_THE_CSV_BYTES"}
```

The CSV requires `engine`, `tag`, `protein`, `ligand`, `site_id`, `pose`, and
`affinity_kcal_mol`; engine/tag/pose rows must be unique. Imports are labeled
`explicit_import_unverified_execution`. Raw engine outputs take precedence, and
an unmarked old derived CSV is not accepted as fresh execution evidence. Missing
pose geometry remains unevaluable even when the imported score can be ranked.

These are intentional compatibility changes. The previous permissive behavior
could produce reassuring but unsupported scientific conclusions.

## Reproducibility and installation

Python 3.10+ is supported by the declared package contract. Install the chemistry
dependencies explicitly:

```bash
python -m pip install -e '.[chemistry,dev]'
python -m pytest test -q
python -m build
python scripts/check_wheel.py
```

Docking executables, Open Babel, and the selected preparation backends remain
external requirements. Optional interaction tools are available through the
`interactions` extra; installing an extra does not install every external binary.
Shell preparation and cluster submission require an appropriate Unix environment
(for example Linux or WSL). Core Windows support does not imply native support
for every scientific backend.
Deployment workers require Python 3.10+. Resume is disabled when an engine
wrapper or container hides an executable identity that cannot be verified.
After abrupt process termination, inspect a leftover job lock before removing it.

Artifact caching hashes file contents, nested input files, outputs, pipeline
source, and relevant dependency versions. Missing required inputs, failed
callbacks, missing outputs, and input changes during computation invalidate
results. Optional input absence is explicitly declared. Workflow updates use
process-safe transactions, and state/cache writes use atomic replacement.
Old and partial generated artifacts are preserved under `.artifact_history`
beside the graph cache instead of remaining on published output paths after a
failed recomputation. Consumers cannot read a failed optional branch's stale
output as new evidence.

The wheel includes the root CLI modules, shell adapter, package resources, and
console entrypoints. CI is configured to exercise Python 3.10/3.12 on Linux and
Windows and build the distribution; these remote jobs have not been run as part
of this local correction.

## Limits of this correction

### Local verification record — 2026-09-12

- Windows, Python 3.12.10: **225 passed, 1 skipped, 0 failures** across the full
  226-test suite (209.78 seconds). See the [historical verification record](../audit/verification-history.md).
  The raw test log and JUnit XML were generated artifacts; they are recoverable
  from the fixed base `8d4434d3a33eb83b1e12cad82944b02c83270e47`. The run emitted
  586 warnings; their details are retained in those recoverable artifacts.
- The skipped real Meeko macrocycle export test could not load RDKit's
  `rdDetermineBonds` DLL because Windows Application Control blocked it. This
  integration remains unverified in this environment.
- All four docking command/configuration dry runs passed. Source and wheel
  distributions built successfully, and every console entrypoint imported from
  the extracted wheel outside the source checkout. Shell syntax, Python
  compilation, and `git diff --check` passed.
- Real docking engines, cluster submission, CUDA execution, preparation backend
  integrations, and ProLIF execution were not validated here. The configured
  Linux/Windows CI matrix has not run remotely.

### Scientific validation still required

Regression tests exercise scientific invariants and workflow failures with small
offline fixtures. They do not establish prospective enrichment, affinity
accuracy, clinical usefulness, or robustness across all protein classes.
Physical clash/strain validation is not implemented in this correction:
`physical_validity_status` remains `not_evaluated`. A reference gate labeled
validated describes its configured RMSD/coverage checks only. It does not
validate force-field parameters, physical pose quality, or predicted affinity.

Before a scientific release, run a versioned benchmark with real engine binaries:
native redocking (including deliberately displaced negative controls), repeated
seeds, cross-docking, known actives/decoys with leakage controls, and targets with
metals/cofactors and difficult protonation states. Report pose validity, RMSD,
coverage, enrichment with uncertainty, resource use, and all excluded cases.
Keep engine versions, prepared structures, atom maps, parameters, and results in
the benchmark manifest. New assembly modeling, protonation ensembles, induced
fit, and calibrated affinity prediction remain separate scientific work.
