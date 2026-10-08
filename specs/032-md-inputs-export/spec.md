# Feature Specification: Strict MD-Consumer Input Export

**Feature ID**: `032-md-inputs-export`  
**Created**: 2026-08-21  
**Status**: Gate 1 specification and Gate 2 scientific baseline approved by the Scientific Lead on 2026-08-21; engineering implementation and dependency-backed chemistry verification complete, named-consumer validation and Scientific Lead artifact review pending  
**Source**: Scientific Lead memo “`md_inputs` export stage — engineering spec (v2, validated against CHARMM-GUI)”  

## 1. Outcome

Add a selection-driven post-docking export stage that regenerates molecular-dynamics preparation inputs from the authoritative selected pose and explicit chemistry/receptor lineage. The stage must fail closed when pose identity, bond topology, receptor composition, protonation policy, or consumer-readiness cannot be proved.

This feature prepares inputs for a named downstream consumer profile. It does not parameterize a force field, build a solvated system, prove a protonation state, or claim that a docking pose is experimentally correct.

## 2. Corrected Architecture Decisions

The source memo is accepted with these engineering corrections:

1. `4-Working/scores/unified/best_pose_per_tag_by_engine.csv`, not a derived complex and not the global top-pose atlas, is the pose-selection authority.
2. The selected row’s one-based `pose` and `pose_file` identify the coordinate record. Its normalized affinity is checked against the native record when that native score is present.
3. GNINA SDF is a connectivity-bearing pose source. Arbitrary Vina/Smina PDBQT is not: PDBQT lacks bond orders and carbon-bound hydrogens. A Vina/Smina export therefore also requires either preserved Meeko preparation metadata or an explicit SDF/MOL/MOL2 topology plus proven atom mapping. Otherwise the row is `not_comparable` and no MD files are emitted.
4. Receptor PDBQT is not an MD receptor source. Export requires an explicit prepared PDB/mmCIF lineage. PDBQT may be used only to verify docking-coordinate identity.
5. Open Babel pH correction is a model prediction. Its output charge and protonation state are recorded as predicted unless the Scientific Lead supplies an approved per-ligand state/charge record.
6. `CONECT` preserves PDB connectivity for the consumer but does not replace the bond-order-bearing ligand MOL2/SDF.
7. `system.pdb` is the canonical combined filename. The older `complex.pdb` name is not emitted unless a later compatibility requirement is approved.
8. Output directories use stable identifiers and include the engine to prevent collisions: `5-Analysis/md_inputs/{protein_key}/{ligand_key}/{engine}/`. Display names are metadata, not path identity.
9. The initial strict profile is named `charmm_gui_cgenff_v1`. Passing its validation means “conforms to this local profile,” not universal CHARMM-GUI, CGenFF, or MD readiness.
10. The node is registered in the artifact graph but is not a dependency of reports or the ordinary `full` analysis target. It runs only when explicitly requested.

## 3. Inputs

Required for every export row:

- `4-Working/scores/unified/best_pose_per_tag_by_engine.csv` with `engine`, `tag`, `protein`, `ligand`, `pose`, `pose_file`, and normalized score fields.
- Exact engine and tag/allowlist selection supplied to `analyze md-inputs`.
- The current selected pose file and its SHA-256.
- A connectivity-bearing ligand topology source and a deterministic topology-to-pose atom mapping.
- An explicit receptor PDB/mmCIF source and its SHA-256.
- A named consumer profile.
- An explicit pH and protonation policy.
- An explicit charge authority: Scientific Lead-provided charge/state record or an approved predicted-charge workflow.

Optional metadata include protein/ligand display-name mappings, native engine version and executable path, preparation-tool versions, container digest, docking seed, box definition, and force-field intent. Missing provenance fields are recorded as `unknown`; required chemistry fields cannot be downgraded to `unknown` while still emitting an export.

## 4. Outputs

For each completed row:

```text
5-Analysis/md_inputs/{protein_key}/{ligand_key}/{engine}/
├── receptor.pdb
├── ligand.mol2
├── ligand.sdf                 # optional requested fan-out
├── ligand.pdb                 # optional requested fan-out
├── system.pdb
└── provenance.json
```

The node also writes:

- `5-Analysis/md_inputs/md_inputs_manifest.json` — every requested row, including failures/skips.
- `.meta/md_inputs_config.json` — normalized invocation settings and dependency hashes used for cache invalidation.

No partially validated molecular files are published into a completed row directory. Failed work may use a temporary directory and is removed or quarantined before return; the manifest retains the failure.

## 5. Validation Gates

### Tier A — data correctness

- **G1 Pose selection**: select the exact one-based pose from the authoritative best-by-engine row. Compare the normalized selected score to the native pose score when available, using an engine-specific tolerance and field mapping. Missing native score is explicit, not fabricated.
- **G2 Staleness/provenance**: hash the current pose, topology, receptor, selection CSV, and configuration. Cache reuse requires matching content hashes and tool/method versions, not only mtime/size.
- **G3 Bond fidelity**: require a connectivity-bearing topology and proven atom mapping. Never infer bond order from PDB/PDBQT or by distance.
- **G4 Assembly integrity**: parsed `system.pdb` atom count equals receptor atoms plus the exact protonated ligand atom count; atom serials are unique and coordinates are finite.
- **G5 Receptor identity and composition**: validate the expected polymer class and configured minimum residue/atom evidence. A protein-required profile must contain recognized amino-acid residues. Nucleic-acid or mixed systems require a different explicit profile; they are not silently quarantined as garbage.

### Tier B — `charmm_gui_cgenff_v1` consumer readiness

- **G6 Protonation and charge**: apply the approved ligand protonation policy at the explicit target pH, add all explicit hydrogens, preserve heavy-atom coordinates within tolerance, validate valence/sanitization with the selected chemistry backend, and record net formal charge plus its authority (`human_approved` or `predicted`). Zero implicit-H counters alone are not sufficient evidence.
- **G7 Naming and connectivity**: ligand atom names are unique within the residue, at most four characters, deterministic, and stable across MOL2/PDB outputs. `CONECT` records cover every ligand bond edge representable in PDB; MOL2/SDF retains bond orders.
- **G8 Strict PDB contract**: write fixed-width PDB records with deterministic altloc resolution, nonblank chain ID, correct residue columns, element columns, serial/residue bounds, `TER`, `CONECT`, and `END`. Parse back with an independent strict validator and compare identity/counts.

Every gate emits one of: `completed`, `failed`, `skipped_disabled`, `skipped_missing_dependency`, `skipped_missing_configuration`, or `not_comparable`, with a reason and evidence fields.

## 6. CLI Contract

Proposed command:

```text
python main.py analyze md-inputs \
  --project-dir PROJECT \
  --engine ENGINE \
  --tags-file TAGS.csv \
  --consumer-profile charmm_gui_cgenff_v1 \
  --ph PH \
  --protonation-policy POLICY \
  --charge-map CHARGES.csv
```

Additional bounded flags:

- `--ligand-format {mol2,sdf,pdb}` may be repeated; MOL2 is always emitted for the CGenFF profile.
- `--topology-map CSV` supplies explicit per-ligand topology lineage when it is not already in project metadata.
- `--receptor-map CSV` supplies explicit per-receptor PDB/mmCIF lineage.
- `--no-protonate` is allowed only with a supplied, already-explicit-H ligand topology whose charge/state is approved and validates under G6.
- `--force` invalidates the requested export cache.

`--receptor-clean` and unconstrained `--resname-scheme` are rejected from v1 because they hide chemically meaningful mutations. Cleaning and residue naming must be deterministic named policies with provenance.

The direct command returns nonzero if any requested row fails a required gate. The optional DAG node may be recorded as `skipped_optional`, but the row manifest must retain the underlying scientific status rather than calling the run successful.

## 7. Dependency Policy

- No new hard dependency is silently added to the base analysis environment.
- The chemistry backend is capability-detected and version-recorded.
- Open Babel may implement the first approved protonation/export backend, but absence produces `skipped_missing_dependency` and no molecular outputs.
- Meeko metadata may support Vina/Smina topology recovery only when its preparation remarks and mappings are present and validated.
- Tests use backend-independent contract fixtures plus dependency-gated integration tests. A lenient Open Babel round-trip is necessary but not sufficient.

## 8. Scientific Decisions Required Before Implementation

The Scientific Lead must approve:

1. Target pH and whether it is required per invocation or may default to 7.4.
2. Charge/protonation authority: require a per-ligand approved map, or allow Open Babel prediction with `human_review_required` status.
3. Expected receptor polymer class and the receptor source/cleaning policy.
4. Initial export scope: exact tag allowlist, all classified hits, or another explicit rule.
5. Whether `charmm_gui_cgenff_v1` is the only v1 consumer profile and whether receptor protonation is in or out of this stage.

## 9. Acceptance Criteria

1. A multi-pose GNINA SDF fixture selects a non-first pose by the manifest row and validates its native score.
2. Changing the pose file without changing its name invalidates the cache and changes the recorded hash.
3. A PDBQT-only Vina/Smina row is `not_comparable` and emits no molecular files.
4. A Vina/Smina row with proven topology lineage preserves graph, atom mapping, pose coordinates, and selected score.
5. A DNA-only receptor fails a protein-required profile with an explicit G5 reason.
6. Ligand MOL2 and `system.pdb` contain the same protonated atoms, coordinates, names, and charge authority.
7. Ligand names are unique; every topology bond is covered by MOL2/SDF and every representable edge by PDB `CONECT`.
8. Deliberately shifted PDB columns fail G8 even if Open Babel accepts the file.
9. Heavy-atom coordinates remain unchanged within the approved tolerance after protonation.
10. Missing chemistry backend, missing pH/charge policy, ambiguous topology mapping, and missing receptor lineage each produce distinct non-success statuses.
11. The aggregate manifest and every `provenance.json` validate against a versioned schema and contain SHA-256 hashes, method/tool versions, input lineage, pose number, score, pH, net charge, charge authority, and per-gate evidence.
12. A retained consumer-validation fixture loads under the named profile; the user-reported 2OV5/4ZBE/6LL5 cases remain `reported_unverified` until sanitized inputs and consumer evidence are added to the repository.

## 10. Non-Goals

- No force-field parameter generation, solvation, minimization, MD launch, CHARMM-GUI upload, or remote execution.
- No receptor protonation, histidine-state choice, termini capping, alternate-location choice, water/ion retention, or residue mutation without an approved named policy.
- No bond-order inference from PDB/PDBQT.
- No claim that Open Babel’s predicted microspecies is authoritative for unusual functional groups.
- No reuse of `5-Analysis/complexes/` as an MD chemistry source.
