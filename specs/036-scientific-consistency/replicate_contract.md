# Replicate contract (Spec 036 R5b)

Producer: `docking/runners/base.py` (dock run). Consumers: post-docking pooling and the
`pose_reproducibility` metric (not implemented in the docking layer).

## Seeds
- `--seed <base>` is required for every dock run. Without it the run fails.
- `--replicates N` (default 3; 1 is allowed) runs N independent seeded replicates per pair.
- Replicate k (1-based) uses seed `base + k - 1`, so seeds are `base, base+1, ..., base+N-1`.
- Engines: `vina`, `smina`, `gnina`. AutoDock4 always runs one replicate and records `replicates: 1`.

## Pose file naming
- Pair tag: `<receptor stem>_<site_id>_<ligand file name>`. Receptor stem has no `.pdbqt`
  suffix (Spec 036 R6), for example `1IEP_A_protein_site_1_STI_A_201.pdbqt`.
- Replicate job tag: `<pair tag>__rep<NN>` with NN two-digit, 1-based.
- Pose file: `engines/<engine>/poses/<pair tag>__rep<NN><ext>` (`.pdbqt` for Vina/Smina, `.sdf` for GNINA, `.dlg` for AutoDock4).
- Log file: `engines/<engine>/logs/<pair tag>__rep<NN>.log`.
- Regex for parsing: `^(?P<pair>.+)__rep(?P<rep>\d{2})$` applied to the file stem. A stem that does not match is a legacy single pose with no replicate id.
- Legacy pose names (`<receptor file name>_site_..._<ligand>`, pre-Spec 036) are still read as their pair.

## Run manifest (`engines/<engine>/run_manifest.json`)
- `effective`: `{engine, exhaustiveness, num_modes, energy_range, replicates, seeds, replicate_contract}`.
- `jobs[]`: one entry per (pair, replicate), with:
  - `tag` (replicate job tag), `pair_tag`, `replicate_id`, `seed`, `pose_file` (the output path);
  - `replicates`: the full list for this pair, `[{replicate_id, seed, output_path}, ...]`;
  - `box_method`, `ligand_rg_angstrom`, `edge_angstrom`;
  - `exhaustiveness`, `num_modes`, `energy_range`;
  - `status`, `returncode`, `poses_returned`, `warnings` (for example `poses_returned=2 is fewer than num_modes=5`).
- `poses_returned` is the number of `MODEL` records in a Vina/Smina pose file. It is null for dry runs and for engines without this count.

## Normalised scores (`engines/<engine>/scores/normalized_scores.csv`)
- `tag` is the pair tag (replicate suffix removed), so rows of different replicates share one `tag`.
- `replicate_id` and `seed` identify the replicate. `pose_file` gives the exact file.

## Project manifest (`project_manifest.json`)
- `engine_effective_parameters`: the `effective` block per engine.
- `box_warnings`: box warnings from the pairlist and the dock run.

## Not in this contract
- Pooling of replicates for pose selection, and the `pose_reproducibility` metric (in-place RMSD between replicate top poses, and the fraction within 2 Å). These belong to post-docking analysis and must read the naming above.
