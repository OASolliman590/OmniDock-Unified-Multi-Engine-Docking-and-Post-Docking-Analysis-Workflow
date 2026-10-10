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

## Consumers (post-docking analysis)
- Naming helpers: `post_docking_analysis/replicate_names.py` (shared by `docking/runners` and post-docking analysis).
- Pooling: the replicates of one `(engine, protein, tag)` are pooled by the shared selector
  (`post_docking_analysis/pose_selection.py`). The chosen row keeps `replicate_id`, `seed` and `pose_file`. Exact
  ties across replicates are broken by `replicate_id` (lower first), so the choice does not depend on row order.
- Seeds: read from `engines/<engine>/run_manifest.json` (`jobs[].seed`, falling back to `effective.seeds`).
- `pose_reproducibility` (`post_docking_analysis/replicates.py`, Spec 036 R5b): the top pose of each replicate is
  compared pairwise with `compare_graph_poses_in_place` (Spec 031/034 mapping). PDBQT poses use the Meeko REMARK
  lineage reader, SDF poses the SDF graph loader. Columns in the best-pose table:
  `pose_reproducibility_status` (`completed`, `partial_not_comparable`, `not_comparable`, `single_replicate`,
  `no_poses`), `pose_reproducibility_replicates`, `pose_reproducibility_pairs`,
  `pose_reproducibility_pairs_comparable`, `pose_reproducibility_max_rmsd_angstrom`,
  `pose_reproducibility_median_rmsd_angstrom`, `pose_reproducibility_fraction_within_2A`
  (share of comparable pairs with RMSD <= 2.0 Å), `pose_reproducibility_method`, `pose_reproducibility_reasons`.
  One replicate gives `single_replicate` and no RMSD values.
- Coverage: `engine_detector` counts pairs (`pose_pairs_found`), not replicate files.
- GNINA pose export (`pose_extractor.py`): the replicate suffix is split off the log/score tag. The
  export reports the pair tag and records `source_tag`, `source_replicate_id` and `source_seed`.
