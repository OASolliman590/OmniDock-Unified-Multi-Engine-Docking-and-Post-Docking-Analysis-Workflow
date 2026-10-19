# USAGE

One entry point: `python main.py <group> <command>` from the repository root. The console script `pdb-prepare-wizard` runs the same CLI. Replace the placeholders in angle brackets. Every command and flag below is checked against `python main.py <command> --help`.

Scientific values (pH, force field, pair mode, engines, seeds, box policy) are the Scientific Lead's decisions. The examples show one explicit choice; nothing here is a recommendation.

## One target, one ligand

Run these in order. Each step reads the project created in step 1.

### 1. Initialise the project
```bash
python main.py workflow init --project-dir <project> --layout-profile canonical \
    --engines vina --project-name <name>
```

### 2. Set the protonation policy (project default)
```bash
python main.py workflow protonation-policy --project-dir <project> \
    --receptor-ph <pH> --receptor-force-field <FF> \
    --ligand-policy <ph_model|explicit_state|as_input> --ligand-ph <pH>
```
Records the receptor pH and force field and the ligand policy in the project manifest. Every later stage reads it. Conflicting values fail with `protonation_policy_conflict`; pass `--replace` only to change a stored value on purpose. `--ligand-state-map <csv>` is required for `explicit_state`.

### 3. Prepare the receptor
```bash
python main.py pdb prepare-protein --receptors-input <raw_receptors_dir> \
    --receptors-output <project>/prepared_proteins --force-field <FF> --ph <pH> \
    --ligand-backend engine_aware_full
```

### 4. Prepare the ligand
```bash
python main.py pdb prepare-ligand --ligands-input <raw_ligands_dir> \
    --ligands-output <project>/prepared_ligands --force-field <FF> --ph <pH> \
    --ligand-backend engine_aware_full
```
`pdb prepare-both` does steps 3 and 4 in one call with the same flags.

**Prepare once (Spec 037 R4a).** An existing prepared receptor or ligand is reused, and the step reports it as `reused_unchanged` (in `receptor_preparation_skipped.csv` / `ligand_preparation_skipped.csv` and in `preparation_summary.txt`), when its provenance matches all of: the resolved protonation policy hash, the request (pH, force field, profile, engines, state map), the input file SHA-256, and the output file SHA-256. A changed input, a changed policy or request, or an edited output is prepared again. Add `--force` to prepare every entity again:
```bash
python main.py pdb prepare-ligand ... --force
```

### 5. Build the pairlist (box centre and box size)
```bash
python main.py prep pairlist --project-dir <project> --mode cocrystal_only \
    --prepared-proteins <project>/prepared_proteins --prepared-ligands <project>/prepared_ligands \
    --excel <project>/proteins_raw/multi_pdb_analysis.xlsx \
    [--raw-ligands <raw_ligands_dir>]   # crystal-frame reference for the containment check
```
- `--mode` is one of `cocrystal_only`, `cocrystal_plus_all`, `cocrystal_plus_nonreference`, `curated_cartesian`, `curated_per_protein`.
- Box centre: from the Excel Summary. Each centre must carry `binding_site_center_method` (`binding_site_center_v1` or `user_explicit`), otherwise the build is refused.
- Box size (default `rg_scaled_v1`): cubic edge = 2.9 x radius of gyration of the ligand's heavy atoms, from the prepared ligand. Recorded per pair as `box_method`, `ligand_rg_angstrom`, `edge_angstrom`.
- Fixed size only when you ask for it: add `--box-size <Å>` (recorded as `user_fixed`).
- Warnings: the box does not contain the reference ligand's heavy atoms, or the edge exceeds 30 Å.

### 6. Materialise the docking project
```bash
python main.py prep project --project-dir <project> --pairlist-file <project>/pairlist.csv \
    --prepared-proteins <project>/prepared_proteins --prepared-ligands <project>/prepared_ligands \
    --excel <project>/proteins_raw/multi_pdb_analysis.xlsx --engines vina --asset-mode symlink
```

### 7. Dock
```bash
python main.py dock run --project-dir <project> --engines vina \
    --vina-binary <path/to/vina> --vina-cpu 4 --seed <base_seed> \
    [--replicates 3] [--exhaustiveness 32] [--energy-range 3] [--dry-run]
```
- `--seed` is required. Replicate k uses seed `base + k - 1`.
- Defaults: exhaustiveness 32, 3 replicates, `--energy-range 3` (Vina/Smina), parameter preset `exhaustive`.
- Basic mode takes exhaustiveness and num_modes from `--parameter-preset` (`screening_fast`, `balanced` or `exhaustive`) and refuses `--exhaustiveness` and `--num-modes`; pass `--parameter-mode advanced --exhaustiveness <n> --num-modes <n>` to set them. The effective values and their sources are printed and written to `run_manifest.json`.
- Every replicate writes `<pair_tag>__repNN.<ext>` and is listed in `run_manifest.json` (see `specs/036-scientific-consistency/replicate_contract.md`).
- Use `--dry-run` to write commands and manifests without running the engine.

#### 7a. Deploy (HPC export, Spec 037 R4c)
```bash
python main.py dock deploy --project-dir <project> --engines vina --seed <base_seed> \
    [--replicates 3] [--energy-range 3] [--parameter-mode basic --parameter-preset exhaustive] \
    [--box-scale 1] [--box-padding 0] [--round <round>] [--mode screen|exhaustive]
```
`dock deploy` applies the same parameter schema as `dock run`:
- `--seed` is required. `--replicates` must be at least 1. `--energy-range` is always passed to Vina/Smina, and it must be positive. Its value is recorded.
- Basic mode refuses `--exhaustiveness` and `--num-modes`, as `dock run` does. Use `--parameter-mode advanced` with both to set them.
- `--box-scale` and `--box-padding` apply to the deployed pairlist boxes, as they do in `dock run`.
- The bundle's `deployment_manifest.json` records `parameter_schema` and `engine_effective_parameters` for each engine: the values, and their sources (`basic_preset:<name>`, `user_explicit` or `advanced_default`).

### 8. Analyse
```bash
python main.py analyze comparative --project-dir <project>

python main.py analyze md-inputs --project-dir <project> --engine vina --tags-file <tags.csv> \
    --receptor-map <receptor_map.csv> [--topology-map <map.csv>] [--charge-map <charges.csv>]
```
`analyze md-inputs` takes the protonation policy from step 2 (the project manifest), so the ligand state is the docked microspecies by default. `--engine` is one of `gnina`, `vina`, `smina`. The `--ph` flag is legacy and only applies to the `openbabel_predicted` policy.

**Top-pose policy and cross-engine comparison (Spec 037 R4b).** `--top-pose-policy` defaults to `best_consensus`, which ranks by the v2 consensus rank (`consensus_rank_geometry_qc_v2`, or `single_engine_native_v1` for one engine). `best_affinity` compares raw affinities, so it is accepted for one engine only and is refused for two or more engines with an error. `--min-affinity-advantage` gated on a raw cross-engine difference and was removed; it must stay `0`. Rerun promotion takes the winning engine of each tag from the consensus. Cross-engine affinity tables, plots and summaries stay, and they are labelled `descriptive_not_comparable_across_engines` in their column names, headers and manifests.

## Many proteins and many ligands

The sequence is the same. What changes:

1. **Raw inputs.** Put the receptors and the ligands in their directories (`--receptors-input`, `--ligands-input`). `pdb prepare-protein`, `pdb prepare-ligand` and `pdb prepare-both` prepare every entity they find.
2. **Protonation overrides per entity (Spec 037 R1).** The project default from step 2 applies to every entity unless overridden. Set an override for one receptor or one ligand, by its file stem:
   ```bash
   python main.py workflow protonation-policy --project-dir <project> \
       --receptor <1IEP> --receptor-ph <pH> --receptor-force-field <FF>

   python main.py workflow protonation-policy --project-dir <project> \
       --ligand <STI> --ligand-policy <ph_model|explicit_state|as_input> --ligand-ph <pH>

   python main.py workflow protonation-policy --project-dir <project> \
       --ligand <STI> --ligand-policy explicit_state --ligand-smiles "<SMILES>"
   ```
   An override that differs from a stored value is refused unless you pass `--replace`.
   To remove an override, the entity inherits the project default again. A prepared entity whose policy hash changed is prepared again on the next run:
   ```bash
   python main.py workflow protonation-policy --project-dir <project> --clear-receptor <1IEP>
   python main.py workflow protonation-policy --project-dir <project> --clear-ligand <STI>
   ```
3. **Pairing.** Choose the pair mode with `prep pairlist --mode <mode>` (or `prep project --pair-mode <mode>`). The modes are listed in step 5. Pairing is a scientific decision.
4. **Prepare once.** Entities whose outputs match their provenance are reused (step 3). Use `--force` only for a full re-run.
5. Steps 5 to 8 then run once for the whole project, on the pairs that the pairlist selects.

## Deprecated commands (kept for one release)

These entry points still run. Each prints a one-line notice naming the `main.py` command it replaces, then forwards to it, or refuses with exit code 2 when no equivalent exists.

| Legacy entry point | Action | Equivalent |
|---|---|---|
| `pdb-wizard-legacy-interactive`, `python interactive_pipeline.py` | forwards (takes no arguments) | `python main.py pdb legacy-interactive [-o OUTPUT]` |
| `pdb-wizard-cli`, `python cli_pipeline.py` | forwards; `--create-config` refused; `-v` ignored | `python main.py pdb run -p <pdbs> -o <output> [-c <config>]` |
| `pdb-wizard-batch`, `python batch_pdb_preparation.py` | forwards; `-v` ignored | `python main.py pdb batch -c <config> [-o <output>]` |
| `python autodock_preparation.py` | forwards; `--config` and `--create-config` refused; file inputs refused (use a directory); `--ligand-profile` becomes `--ligand-backend` | `python main.py pdb prepare-both --receptors-input ... --ligands-input ... --force-field ... --ph ...` |
| `pdb-wizard-prepare-docking` | forwards (same flags) | `python main.py prep project ...` |
| `pdb-wizard-dock` | forwards; 14 legacy-only flags refused (`--ligand-qc-*`, `--no-ligand-admet-filters`, `--autodocktools-*`) | `python main.py dock run ...` |

Also accepted by `main.py` itself, mapped without a notice: `interactive`, `cli`, `batch`, `prepare-docking`, and `dock` without a subcommand.

Not converted yet: `post-docking-analysis` (a second post-docking path; its routing is decided under Spec 038 R3). `webui_cli.py` is still a separate script; `python main.py webui` is the supported launcher.

`prep_autodock_enhanced.sh` is the internal shell backend that the preparation code calls. It is not a user entry point. Use `pdb prepare-*`.

Known differences in forwarded legacy commands:
- `pdb-wizard-batch` exits 0 when the batch status is `partial` (for example, when every entry failed). The legacy script exited 1 when no entry succeeded.
