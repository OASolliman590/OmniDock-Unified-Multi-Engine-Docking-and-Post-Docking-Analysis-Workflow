# USAGE: one target, one ligand

This is the single supported path (Spec 036 R6). Run every command from the repository root with `python main.py ...`. Replace the placeholders in angle brackets.

## 1. Initialise the project
```bash
python main.py workflow init --project-dir <project> --layout-profile canonical \
    --engines vina --project-name <name>
```

## 2. Set the protonation policy (one project-wide policy)
```bash
python main.py workflow protonation-policy --project-dir <project> ...
```
Records the receptor pH/force field and the ligand policy (`ph_model`, `explicit_state` or `as_input`) in the project manifest. Every later stage reads it. Conflicting values fail with `protonation_policy_conflict`. The subcommand is being added by Spec 036 R1; until it exists, do not skip this step.

## 3. Prepare the receptor and the ligand
```bash
python main.py pdb prepare-protein --receptors-input <raw_receptors_dir> \
    --receptors-output <project>/prepared_proteins --force-field AMBER --ph 7.4 \
    --ligand-backend engine_aware_full

python main.py pdb prepare-ligand --ligands-input <raw_ligands_dir> \
    --ligands-output <project>/prepared_ligands --force-field AMBER --ph 7.4 \
    --ligand-backend engine_aware_full
```
`autodock/prep_autodock.sh` is retired and only prints a pointer here.

## 4. Build the pairlist (box centre and box size)
```bash
python main.py prep pairlist --project-dir <project> --mode cocrystal_only \
    --prepared-proteins <project>/prepared_proteins --prepared-ligands <project>/prepared_ligands \
    --excel <project>/proteins_raw/multi_pdb_analysis.xlsx \
    [--raw-ligands <raw_ligands_dir>]   # crystal-frame reference for the containment check
```
- Box centre: from the Excel Summary. Each centre must carry `binding_site_center_method` (`binding_site_center_v1` or `user_explicit`), otherwise the build is refused.
- Box size (default `rg_scaled_v1`): cubic edge = 2.9 x radius of gyration of the ligand's heavy atoms, from the prepared ligand. Recorded per pair as `box_method`, `ligand_rg_angstrom`, `edge_angstrom`.
- Fixed size only when you ask for it: add `--box-size <Å>` (recorded as `user_fixed`).
- Warnings: the box does not contain the reference ligand's heavy atoms, or the edge exceeds 30 Å.

## 5. Materialise the docking project
```bash
python main.py prep project --project-dir <project> --pairlist-file <project>/pairlist.csv \
    --prepared-proteins <project>/prepared_proteins --prepared-ligands <project>/prepared_ligands \
    --excel <project>/proteins_raw/multi_pdb_analysis.xlsx --engines vina --asset-mode symlink
```

## 6. Dock
```bash
python main.py dock run --project-dir <project> --engines vina \
    --vina-binary <path/to/vina> --vina-cpu 4 --seed <base_seed> \
    [--replicates 3] [--exhaustiveness 32] [--energy-range 3] [--dry-run]
```
- `--seed` is required. Replicate k uses seed `base + k - 1`.
- Defaults: exhaustiveness 32, 3 replicates, `--energy_range 3` (Vina/Smina), parameter preset `exhaustive`.
- Basic mode takes exhaustiveness and num_modes from `--parameter-preset` (`exhaustive` gives exhaustiveness 32 and num_modes 40) and refuses `--exhaustiveness` and `--num-modes`; pass `--parameter-mode advanced --exhaustiveness <n> --num-modes <n>` to set them. The effective values and their sources are printed and written to `run_manifest.json`.
- Every replicate writes `<pair_tag>__repNN.<ext>` and is listed in `run_manifest.json` (see `specs/036-scientific-consistency/replicate_contract.md`).
- Use `--dry-run` to write commands and manifests without running the engine.

## 7. Analyse
```bash
python main.py analyze comparative --project-dir <project>

python main.py analyze md-inputs --project-dir <project> --engine vina --tags-file <tags.csv> \
    --receptor-map <receptor_map.csv> [--topology-map <map.csv>] [--charge-map <charges.csv>]
```
`analyze md-inputs` takes the protonation policy from step 2 (the project manifest), so the ligand state is the docked microspecies by default. The `--ph` flag is legacy and only applies to the `openbabel_predicted` policy. Final flag set: Spec 036 R1, being added by Agent A.
