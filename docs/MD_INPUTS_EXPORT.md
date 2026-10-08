# Strict MD-input export

`analyze md-inputs` exports exact selected docking poses for the
`charmm_gui_cgenff_v1` consumer profile. It does not run CHARMM-GUI, generate
force-field parameters, protonate the receptor, or start MD.

## Command

```powershell
python main.py analyze md-inputs `
  --project-dir D:\path\to\project `
  --engine gnina `
  --tags-file D:\path\to\tags.csv `
  --receptor-map D:\path\to\receptor_map.csv `
  --ph 7.4 `
  --protonation-policy openbabel_predicted `
  --ligand-format sdf
```

The selected-pose table must already exist at
`4-Working/scores/unified/best_pose_per_tag_by_engine.csv`. Open Babel 3.2.1 is
required to generate new molecular files. When it is absent, the command exits
nonzero and records `skipped_missing_dependency`; it does not publish a partial
row.

## Selection and mapping files

`tags.csv` contains exact tags:

```csv
tag
2OV5_site_1_OX11
4ZBE_site_1_T2Z14
```

`receptor_map.csv` supplies the unchanged prepared PDB or mmCIF lineage. A
tag-specific row is preferred; an unambiguous protein/receptor row is also
accepted.

```csv
tag,protein,receptor_file
2OV5_site_1_OX11,2OV5,D:\prepared\2OV5.pdb
4ZBE_site_1_T2Z14,4ZBE,D:\prepared\4ZBE.cif
```

GNINA uses the selected native SDF record for connectivity. Vina/Smina require
`--topology-map`; arbitrary PDBQT alone is rejected. `atom_map` is a JSON array
mapping each selected PDBQT heavy-atom index to its topology heavy-atom index.
It must be a complete zero-based permutation. One-based permutations are also
accepted and normalized explicitly.

```csv
tag,topology_file,topology_pose,atom_map
2OV5_site_1_OX11,D:\topologies\OX11.sdf,1,"[0,1,2,3,4]"
```

Without `--charge-map`, Open Babel's result is recorded as `predicted` and
`human_review_required`. To supply Scientific Lead-approved integer charges:

```csv
tag,net_charge
2OV5_site_1_OX11,0
```

The stage fails G6 if the generated charge disagrees with an approved map.

## Outputs and interpretation

Completed rows are written to:

```text
5-Analysis/md_inputs/{protein}/{ligand}/{engine}/
├── receptor.pdb
├── ligand.mol2
├── ligand.sdf or ligand.pdb   # only when requested
├── system.pdb
└── provenance.json
```

The aggregate status is in `5-Analysis/md_inputs/md_inputs_manifest.json`, its
run-tracking copy is `5-Analysis/md_inputs/run_tracking/manifest.json`, and the
normalized request/dependency hashes are in `.meta/md_inputs_config.json`.

`completed` means the local G1–G8 profile passed. It does not prove the predicted
microspecies, charge state, force-field assignment, docking hypothesis, or MD
stability. Review every `predicted` charge/protonation state before using a
force-field front-end.
