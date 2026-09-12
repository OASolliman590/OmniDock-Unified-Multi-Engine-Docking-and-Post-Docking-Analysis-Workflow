# Scientific preparation audit

Repository: OmniDock, commit `cfe33056c77d817d19439a9e9affe225d00ca488`. Inspected 12 September 2026. Scope: protein/ligand collection, cleanup, parameterization, preparation quality checks, binding-site metadata, and project materialization. Pipeline source was not modified.

The preparation layer has useful building blocks, but its current success criteria do not establish that the molecular system being docked is the intended chemical system. Highest priorities are explicit receptor composition, chemical identity preservation, and removal of permissive chemistry-changing fallback behavior. The problems below are traceable defects or documented design gaps; they are not evidence that every past docking result is invalid.

## Execution paths actually traced

- Current collection: `workflow/execution.py:626-751` → `cli_pipeline.run_single_pdb_cli` → `core_pipeline` download/extraction/cleanup. The default collector uses `cleaning.default=common`, preserves metals/cofactors, and saves an Excel summary. Interactive collection calls a different preparation function.
- Current parameterization: `workflow/execution.py:1094-1135` constructs `PreparationConfig`, then `AutoDockPreparationPipeline.run_enhanced_preparation` launches `prep_autodock_enhanced.sh`. Receptor handling is in that shell script, not a PDBFixer module. No PDBFixer repair stage was found in these active paths.
- Ligand normalization: `docking/preparation/ligand_preparation.py:428-519`; most profiles use Open Babel first, then Meeko or AutoDockTools. Explicit `meeko_only` and `autodocktools_only` bypass normalization.
- QC: ligand/receptor validators feed `docking/preflight.py` and `docking/cli.py`; these checks affect execution, not just reports. `repair_project_ligands` is a separately callable repair utility; I did not find it automatically invoked by the normal prepare path.
- Legacy batch: `batch_pdb_preparation.py` is independently callable and has different cleanup/parameterization behavior. Do not assume a fix in one interface fixes all interfaces.

No applicable AGENTS.md was present in the audited repository or workspace root. Research used primary documentation. Current documentation sometimes describes newer APIs than an unpinned installation; such version differences are called out instead of treated as reproduced installation failures.

## Findings, ranked by practical impact

### PREP-01 — The default collected receptor retains the selected ligand, and the source complex competes with the cleaned receptor

**Priority P1; high confidence in code/probe; scientific consequence depends on the selected chemistry and backend.**

Evidence: [default cleanup](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/cli_pipeline.py#L263-L295), [collector defaults](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/workflow/execution.py#L662-L671), [receptor enumeration](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/prep_autodock_enhanced.sh#L819-L843), and [project asset indexing](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/project_builder.py#L260-L307). `common` removes solvent names but explicitly does not remove the selected ligand. Both the downloaded `1ABC.pdb` and `1ABC_cleaned.pdb` remain in the raw-protein directory and both qualify for shell preparation. The project builder's PDB-ID alias uses the first sorted match; in the probe, `1ABC` resolves to `1ABC.pdbqt`, not `1ABC_cleaned.pdbqt`. Exact `1ABC_cleaned` resolution is correct; this is not a claim that exact filenames always fail.

**Reproduction:** the actual CLI function, isolated from its unavailable Biopython import, selected `LIG` and returned success after removing only `HOH`. The docking engine was not run. Whether the retained ligand reaches receptor PDBQT depends on the backend: Meeko can reject/drop an unmatched ligand, whereas Open Babel or an applicable template can retain it. Backend-dependent deletion is not a receptor-composition policy.

**Consequence:** the intended search site may remain physically occupied by the crystallographic ligand, or a backend may remove compounds that were intended as cofactors. Docking a new ligand into that occupied pocket asks a different scientific question. The canonical collector has no explicit selected-anchor removal step.

**Remedy:** maintain immutable `source_complex`, `reference_ligand`, and explicitly composed `docking_receptor` artifacts; remove the selected ligand by full instance ID for ordinary redocking/competitive screening. Allow intentional co-docking or retained ligands only through a named experiment mode. Stage exact artifact IDs, reject ambiguous PDB-only lookups, and assert that the excluded ligand's atom mapping does not survive in the receptor. Vina's preparation tutorial explicitly separates receptor and ligand and provides residue-deletion controls; this supports the separation, not a blanket rule to remove every heterogen. [Vina basic docking](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html), [Meeko residue selection](https://meeko.readthedocs.io/en/develop/cli_rec_prep.html).

### PREP-02 — Extraction discards authoritative ligand chemistry before inferring it again from coordinates

**Priority P1; high confidence; metadata loss reproduced, downstream chemical error not benchmarked.**

Evidence: [PDB extraction](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/core_pipeline.py#L433-L491), [ligand sanitizer](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/ligand_quality.py#L323-L407), [PDB-to-SDF conversion](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/workflow/execution.py#L691-L699). Extraction copies one residue through Bio.PDB/PDBIO, and sanitation reads only atom records, renames atoms, and writes neither CONECT nor formal-charge fields. The source CCD identity and source-to-output atom map are not used to establish bond orders or stereochemistry. Converting this PDB to SDF merely asks Open Babel to infer missing chemistry.

**Reproduction:** a two-atom PDB with CONECT records and a `1+` formal charge lost both after the actual sanitizer. A valid file can therefore become a chemically less informative file while reporting success. This does not mean every resulting bond order is wrong.

**Consequence:** ambiguous aromaticity, bond orders, charged groups, stereochemistry, or covalent attachments can be guessed incorrectly. Downstream interaction analysis cannot repair the experiment if the initial molecule was wrong. Vina specifically discourages preparing small molecules from PDB, and Meeko prefers SDF with hydrogens and 3D coordinates. [Vina preparation guidance](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html), [Meeko ligand input requirements](https://meeko.readthedocs.io/en/develop/lig_prep_basic.html).

**Remedy:** retrieve the chemical component's bond/charge definition or an instance SDF; map authoritative chemistry onto the deposited coordinates. RCSB provides CCD definitions and chemical-component-instance SDF endpoints. Preserve the crystallographic reference coordinates and atom mapping separately from screening conformers. Validate heavy-atom connectivity, isomeric identity, formal charge, and missing atoms before export. An ideal-coordinate CCD SDF supplies chemistry but must not silently replace the observed reference pose. [RCSB file download services](https://www.rcsb.org/docs/programmatic-access/file-download-services).

### PREP-03 — Valid flexible-macrocycle atom types are rejected, and the repair removes the ring-closure mechanism

**Priority P1; high confidence; reproduced at the actual validator/repair functions.**

Evidence: [allowed types](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/ligand_quality.py#L19-L44), [CG0/G0 repair](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/ligand_quality.py#L178-L183), [atom removal/retyping](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/ligand_quality.py#L231-L259), [repair invocation](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/ligand_quality.py#L922-L924), [preparation validation](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/autodock_preparation.py#L435-L457).

The validator rejects `CG0` and `G0`. The repair changes `CG0` to carbon and drops `G0`. Those are meaningful macrocycle types, not corrupt atom labels. The rest of the torsion tree is retained without reconstructing a closed-ring molecule.

**Reproduction:** a minimal token fixture containing these two types produced `invalid_atom_types`; repair dropped one atom and retyped one; the repaired file then passed the same validator. This is a lexical/contract reproduction, not a complete physically validated macrocycle docking fixture.

**Consequence:** ordinary preparation can reject valid Meeko macrocycles. Invoking repair can then turn a ring-closure representation into an effectively unconstrained open topology or an otherwise inconsistent torsion tree. Vina describes why the dummy atoms and closure potential are necessary, and recommends Meeko export to restore correct chemical connectivity. [Vina macrocycle protocol](https://autodock-vina.readthedocs.io/en/latest/docking_macrocycle.html).

**Remedy:** use engine/version-specific PDBQT validation, preserve recognized closure types, and never infer a chemical atom type from its name to make a failing structure pass. An unsupported engine should receive an explicitly prepared rigid-ring conformer ensemble or a clear unsupported-capability result. Regenerate from authoritative chemistry; do not delete atoms in an already parameterized torsion tree.

### PREP-04 — Receptor preparation can delete incomplete residues or change chemistry through fallback without a molecular difference check

**Priority P1; high confidence in configuration/control flow; no actual Meeko/Open Babel execution.**

Evidence: [default permissive settings](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/autodock_preparation.py#L33-L44), [workflow configuration](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/workflow/execution.py#L1094-L1109), [receptor backend chain](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/prep_autodock_enhanced.sh#L893-L932), [limited shell validation](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/prep_autodock_enhanced.sh#L649-L680).

`allow_bad_res=True` is default. The shell passes `--allow_bad_res`, suppresses stderr, and accepts success if a minimally valid output file exists. On failure it runs `obabel input -O output -xr`; no atom/residue set comparison checks what disappeared or survived. The optional PDB2PQR path is off by default and also falls back when it fails.

Meeko templates require an exact heavy-atom match; missing atoms can prevent matching. Its permissive option ignores/deletes those residues. **Meeko can add missing hydrogens from templates; it would be incorrect to claim the default Meeko path always lacks hydrogens.** Current development docs rename the old option to `--delete_bad_res`; this is an API compatibility risk, not proof that the user's installed release rejects the flag. [Meeko receptor chemistry](https://meeko.readthedocs.io/en/develop/rec_overview.html), [Meeko receptor options](https://meeko.readthedocs.io/en/develop/rec_cli_options.html).

**Consequence:** a missing side-chain atom near the pocket can remove the whole residue, changing sterics and contacts. A different installed backend can change the receptor while using the same nominal preparation configuration.

**Remedy:** default to strict template matching; retain stdout/stderr and a residue-level change ledger. Fail if pocket residues, required cofactors, or coordinating atoms are removed. Diagnose missing residues/atoms first; repair only justified cases and record modeled regions. PDBFixer exposes missing-residue/atom detection and selective repair, but wholesale loop reconstruction or heterogen removal is not automatically scientifically justified. [PDBFixer manual](https://github.com/openmm/pdbfixer/blob/master/Manual.html).

### PREP-05 — The pH setting is not a consistent molecular-state contract

**Priority P1/P2 depending on target; high confidence in code, state-change effects not simulated.**

Evidence: [ligand protonation fallback](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/ligand_preparation.py#L303-L320), [RDKit path](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/ligand_preparation.py#L244-L275), [fallback metadata](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/ligand_preparation.py#L394-L403), and receptor commands cited in PREP-04.

If Open Babel pH treatment fails, ligand normalization retries simple `-h`. If Open Babel normalization fails altogether, RDKit `AddHs` is used without pH-dependent ionization, yet `protonation_ph` is still attached to the returned metadata. Direct Meeko/ADT profiles likewise report the requested pH without performing this normalization. In receptor preparation, the pH is supplied only to optional PDB2PQR; neither the default Meeko command nor the Open Babel fallback encodes that pH-specific state assignment.

**Consequence:** identical requested pH can yield different protonation/tautomer states across machines and profiles, affecting donor/acceptor typing and, for applicable engines, electrostatics. `AddHs` and ionization-state enumeration are different operations. Open Babel documents distinct `-h` and `-p` behavior, including limitations of its atomwise ionization model. [Open Babel options](https://openbabel.org/docs/Command-line_tools/babel.html).

**Remedy:** save requested pH separately from applied treatment and resulting molecular state. A mandatory pH treatment failure should fail or explicitly produce an untreated artifact. Make protonation, tautomerism, histidine states, terminal groups, and relevant metal coordination configurable molecular states. Enumerate a bounded set of justified ligand states and retain parent/state IDs; report sensitivity across states instead of pretending one uncertain state is exact. This is a targeted improvement, not a requirement to enumerate every conceivable state.

### PREP-06 — Reported pocket volume and druggability are synthetic heuristics presented as measurements

**Priority P1 for scientific reporting; high confidence; reproduced.**

Evidence: [constant volume](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/core_pipeline.py#L742-L749), [residue-count electrostatics](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/core_pipeline.py#L751-L776), [druggability formula/labels](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/core_pipeline.py#L812-L848).

`pocket_volume_A3` is always the volume of an assumed 5 Å sphere: **523.5988 Å³**, independent of receptor geometry. The so-called electrostatic score is a signed residue count with every histidine assigned `+0.5`, not a calculated potential. The arbitrary druggability formula combines these with hydrophobic-residue counts and attaches qualitative labels.

**Reproduction:** the actual pocket method with an empty structural iterator reports volume `523.5987755982989`, druggability `0.508`, and **“Good”**. This directly demonstrates that the label does not require a pocket.

**Remedy:** remove the volume and druggability claims from scientific output until an actual method is selected and validated. Descriptive counts can remain with honest names, units, and assumptions. If cavity volume or electrostatic potential is needed, calculate it from a specified geometric or electrostatic method and record its parameters; do not promote arbitrary formulas to validated predictions. This finding follows directly from the code and does not require an external authority to show that a constant is not a geometry measurement.

### PREP-07 — Site extraction depends on atom order, loses chain identity, and can miss actual ligand contacts

**Priority P1/P2; high confidence in source; mathematical counterexamples, not full PLIP execution.**

Evidence: [PLIP coordinate reconstruction](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/core_pipeline.py#L1100-L1140), [distance fallback](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/core_pipeline.py#L1202-L1236), [fixed box dimensions](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/project_builder.py#L446-L456).

The distance fallback measures protein-atom distances to the ligand **centroid**, not to ligand atoms, and breaks after the first atom satisfying the cutoff. Its “residue average” consequently contains one selected atom per residue. Residue keys omit chain and insertion code. The PLIP path keys by residue name/number, reconstructs only on the ligand's chain, and re-adds all residue atoms for each reported interaction, weighting frequently interacting residues multiple times.

**Counterexamples:** ligand ends at x=−10 and +10 have centroid 0; a protein atom at x=11 is 1 Å from a ligand atom but is excluded by a 5 Å centroid cutoff. A residue with atoms x=1 and x=4 contributes x=1 or x=4 solely according to file atom order. Same-number residues on two chains collide. An interface ligand assigned chain A may contact chain B, which the PLIP reconstruction omits or misidentifies.

**Consequence:** binding-site centers and reported interaction counts can change with file order, backend availability, and chain labels. A fixed 20 Å box has no demonstrated containment of an elongated reference ligand.

**Remedy:** derive a redocking box from the reference ligand's heavy-atom bounding box plus explicit padding, then assert containment. Define residue contact as the minimum heavy-atom distance; key each residue with model, chain, residue number, and insertion code; deduplicate before centroid/count calculations. Use structured PLIP output with protein-side chain identifiers. Meeko already supports an enveloping box and padding. [Meeko box options](https://meeko.readthedocs.io/en/develop/rec_cli_options.html).

### PREP-08 — Alternate-location sanitation can construct a hybrid conformer and erases its provenance

**Priority P2; high confidence; reproduced.**

Evidence: [atomwise selection](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/ligand_quality.py#L365-L380), [renaming and forced occupancy](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/ligand_quality.py#L384-L402).

The preferred altloc is selected separately for each atom. A unique B-only atom is retained even when the rest of the ligand is selected from A; a missing A alternative falls back to whichever location has greatest occupancy. Every output occupancy becomes 1.0, altloc labels disappear, and names change without a mapping.

**Reproduction:** C1(A) at x=0, C1(B) at x=10, and C2(B) at x=11 become two unlabeled atoms at x=0 and x=11. No source conformer contained that pairing. The probe shows the selection failure; it does not assert that every partial-altloc structure is malformed.

**Remedy:** choose a consistent residue/ligand conformer, include only compatible shared atoms, and mark an incomplete conformer for repair/rejection. Save occupancy, original names, altloc selection, and atom mapping. Respect a user-specified altloc; occupancy-based preference alone is an optional policy, not an absolute scientific rule.

### PREP-09 — “Receptor quality” currently accepts nonfinite and otherwise chemically uninterpretable coordinates

**Priority P2; high confidence; reproduced.**

Evidence: [coordinate parsing](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/receptor_quality.py#L29-L51), [QC calculations](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/receptor_quality.py#L92-L181). The validator checks counts and maximum coordinate span, but `float('nan')` is accepted and there is no finiteness check. A malformed line's fallback collects its first three numeric tokens, which can be serial number, residue number, and x rather than x/y/z. A receptor with no valid coordinates has span 0 and can still pass count thresholds.

**Reproduction:** a 100-atom receptor with all x coordinates NaN passes `validate_receptor_file` with no issues. Ligand validation checks atom-type/root/TORSDOF token counts, not finite coordinates, valence, connectedness, branch references, or matching atom sets.

**Remedy:** reject nonfinite/invalid coordinates and duplicate atom identifiers; use a format-aware parser and engine-specific syntax validation; compare source/target heavy-atom sets and ligand identity. Add pocket completeness/clash checks rather than presenting global atom counts as a scientific receptor quality certificate.

### PREP-10 — Biological assembly, model identity, and instance identity are not first-class inputs

**Design gap, generally P2; high confidence; target-dependent scientific risk.**

Evidence: [PDB download](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/core_pipeline.py#L248-L293), [REMARK350 chain-group hints](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/core_pipeline.py#L365-L430), [extraction searches all models](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/core_pipeline.py#L462-L472), [Excel pivot by PDB alone](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/preparation/excel_sites.py#L46-L49). Collection fetches the entry/asymmetric unit; parsing REMARK350 chain group hints does not apply assembly transformations. There is no explicit model selector in extraction. Site metadata collapses to one row per PDB entry.

RCSB explains that the asymmetric unit can contain part, all, or several biological assemblies; selecting chains cannot generate a missing symmetry mate. This matters particularly for interface pockets. [RCSB biological assembly guide](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies).

**Remedy:** identify receptor artifacts by accession + assembly + model + selected chains + molecular state; record symmetry transformations and source coordinates. Preserve multiple binding sites and receptor states instead of collapsing them to PDB ID. Assess a biological assembly for the intended experiment; do not blindly assume every assembly or every crystallographic contact is biologically relevant.

## Additional technical/behavioral concerns

- **Legacy batch AutoDock is a stub (P1 if used):** `batch_pdb_preparation.py:350-353,384-429` marks overall success, calls a routine that never invokes the enhanced preparer, and returns “preparation requested.” Its nested `ligands_dir.mkdir(exist_ok=True)`/`receptors_dir.mkdir` lack `parents=True` while the parent is not created; failure is caught and the caller ignores the return. The file passed for the receptor is the original `pdb_file`, despite the “cleaned PDB” comment. This is not the current unified parameterization path.
- **Legacy helper has a deterministic type mismatch:** `autodock_preparation.py:340-356` indexes `hetatm_details[ligand_name]`, but `enumerate_hetatms` returns a list of tuples. The helper catches the error and returns None. No active caller was found, so this is a dormant exported-API defect rather than a demonstrated normal workflow failure.
- **Resume is not chemistry aware:** `prep_autodock_enhanced.sh:877-881` reuses any existing minimally valid receptor file; it does not compare source hash, pH, backend, selected atoms, or tool version. `project_builder.py:464-474` similarly skips existing destinations. Changing inputs/settings can leave scientifically stale artifacts. Fix with content/config/version hashes and an immutable artifact manifest.
- **ADMET heuristics are conflated with structural QC:** `ligand_quality.py:47-65,779-805` enables drug-likeness, PAINS/Brenk, and reactive SMARTS filters by default, including a generic nitro-group/thiol alert. Those can be useful project-specific screening policies, but they do not establish ADME/toxicity outcomes or invalid molecular geometry. A general docking pipeline should separate these policy flags from broken-input failures, expose exclusions, and support fragment, cofactor, natural-product, or covalent projects explicitly. RDKit-unavailable ADMET checks can also be skipped with a warning, so results can vary by installation.
- **Source geometry and conformer geometry should be distinct:** the Open Babel path always calls `--gen3d` and minimization (`ligand_preparation.py:297-361`) and reports `had_3d_input=False` even when a 3D ligand was supplied. Generating a docking starting conformer is legitimate, but it must not become the experimental reference. The RDKit fallback preserves existing 3D coordinates, so these paths are not equivalent. Preserve both artifacts, record generation seed/settings, and compare intended molecular identity.
- **Positive foundations:** exact engine profiles, explicit optional cofactor/metal preservation, fixed-seed RDKit embedding, JSON step reports, and early QC hooks provide useful seams for improvement. The central issue is the meaning and enforcement of their contracts.

## A concrete redesign

1. **Molecular identity layer:** immutable source mmCIF + CCD/instance ligand chemistry; stable IDs for receptor assembly/model/state, ligand parent/state/conformer, site, and atoms. Engine files are exports, not authoritative chemistry.
2. **Preparation decisions layer:** explicit selected/excluded residue instances, alternate conformers, protonation/tautomer states, modeled atoms, cofactors/metals/waters, and known unsupported chemistry. Emit a human-readable molecular difference report and a machine-readable decision manifest.
3. **Strict engine adapters:** pinned tool versions and capability tables; identical intended chemistry across exports; unsupported macrocycle/metal/covalent features fail explicitly or invoke a documented alternative protocol. Any fallback records a changed protocol and reruns chemical equivalence checks.
4. **Site-aware validity gates:** finite coordinates, graph/charge/stereochemistry identity, heavy-atom completeness, pocket residue preservation, ligand exclusion/intentional retention, steric clashes, reference containment in the search box, and round-trip export checks. Validation reports should distinguish syntax, chemistry, scientific assumptions, and optional screening policies.
5. **Scientific validation suite:** a small curated set of ordinary ligands, charged ligands, tautomers, macrocycles, metal/cofactor systems, alternate conformers, interface pockets, incomplete structures, and multiple models. Include crystallographic redocking with symmetry-aware heavy-atom RMSD, controlled starting conformers/seeds, and preparation-variant sensitivity. Set acceptance criteria before choosing defaults. This suite is future work, not testing completed in this audit.

Near-term order: stop synthetic pocket claims and destructive macrocycle repair; enforce an explicit receptor/reference split and exact asset selection; make bad-residue/protonation fallbacks visible and fail when required chemistry is lost; preserve authoritative ligand chemistry; then add assembly/state modeling and scientific benchmark coverage.

Proposed first release criteria are concrete: all source-to-prepared heavy atoms must be mapped or listed as intentional additions/deletions; all coordinates must be finite; selected alternate locations must be internally consistent; each ordinary redocking receptor must exclude the exact reference ligand; every reference heavy atom must lie within the configured box; source/config/tool changes must invalidate cached preparation; supported macrocycle closure types and their atom references must survive export. For a predeclared redocking benchmark, report symmetry-aware heavy-atom RMSD against the unmodified reference, top-1 and top-5 success at a predeclared cutoff (for example 2 Å), and results across several fixed search seeds. Set target-specific pass-rate goals from the baseline and a held-out set; do not tune a success cutoff after seeing the results or treat redocking success as proof of affinity prediction.

## Verification performed and limits

Ran [preparation-probes.py](preparation-probes.py) successfully on local Python 3.12 with NumPy/pandas and the actual lightweight preparation/QC modules. Seven assertion groups reproduced the behaviors described above: macrocycle tokens; nonfinite receptor QC; ligand charge/connectivity loss; mixed alternate conformers; PDB-only source selection; default selected-ligand retention; empty-pocket “Good” scoring. CLI/pocket methods were extracted with Python AST to isolate unavailable Bio imports; external operations were replaced with deliberately small fixtures. Initial execution failed because the Windows sandbox denied temporary-directory access; the approved rerun completed.

Biopython, RDKit, Meeko, Open Babel, and Bash were unavailable in the local tested Python/PATH environment. No packages were installed; no molecular docking, PLIP, PDBFixer, PDB2PQR, quantum chemistry, or real protein repair was run. Metadata/protocol defects are established; quantitative effects on existing research results require rerunning representative systems with corrected preparation. No pipeline source code was changed, and no scientific result was fabricated to fill a missing measurement.
