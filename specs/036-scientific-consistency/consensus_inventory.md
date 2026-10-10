# Spec 036 R3a — Consensus mode inventory

**Scope:** every consensus mode name and every ranking path that writes a consensus table in
`post_docking_analysis/` and the CLIs (checked at HEAD 2e08308, before Spec 036 R3 changes).
**Reference policy:** `specs/031-appraisal-remediation-p0/decisions/consensus-v2.md` (`consensus_rank_geometry_qc_v2`)
and `atom-mapping.md` (`spec031-atom-mapping-v1`).
**Forbidden methods (Spec 036 R3a / Spec 031 atom-mapping):** centroid "soft alignment",
atom matching by file order or element order, mapping-free RMSD.

## Mode-by-mode verdicts

| Mode | What it computes | RMSD / mapping used | Forbidden method? | Duplicates v2? | Verdict |
|---|---|---|---|---|---|
| `consensus_rank_geometry_qc_v2` | Per-engine equal-weight rank percentiles over the complete-case ligand set per target; composite = mean of percentiles; GNINA pose = highest `cnn_score` (tie: `cnn_affinity`, Vina affinity, pose, tag); Vina/Smina/AD4 pose = lowest affinity. Geometry is a QC tie-breaker only. | Spec 031 graph/symmetry mapping (`compute_geometric_consensus_v2`, `compare_graph_poses`). | No. | n/a (the reference) | **keep, default everywhere** |
| `dockbox_geometric` | Composite = geometric agreement flag + tiny mean rank term. Best row per (engine, tag) by minimum raw affinity across engines. | `compute_geometric_consensus` -> `_pairwise_rmsd` + `_kabsch_rmsd`. Atoms paired by **element sequence and file order**; soft alignment by **centroid-distance sorting** when atom counts differ by <= 2; no bond-graph check. | **Yes**: file-order / element-order matching, centroid soft alignment, mapping-free RMSD. | Yes (geometry is what v2 does with a proof) | **remove** |
| `weighted_hybrid` | 0.55 mean rank percentile + 0.30 best-affinity rank percentile + 0.10 agreement fraction + 0.05 (1 - spread). Raw affinities from different engines are ranked on one scale. Also the **forced fallback for single-engine projects** before this spec. | none | No RMSD. Mixes raw engine scores, which Spec 031 forbids (score incommensurability). | Yes: a cross-engine composite without a proof or equal weights. | **remove** |
| `strict_consensus` | Keeps tags with `agreement_count >= 2`, then 0.65 mean rank + 0.25 best-affinity rank + 0.10 (1 - spread). Single-engine case only warns. | none | No RMSD. Same raw-score mixing as `weighted_hybrid`. | Yes: v2's complete-case rule already is strict agreement, with no imputation or renormalisation. | **remove** |
| `favorite_guardrails` | Favorite engine's rank percentile minus penalties (agreement < 2, spread, winner != favorite). | none | No RMSD. Gives one engine a privileged weight. | Yes, and the purpose is engine-specific weighting. consensus-v2.md says engine-specific weights need separate scientific approval. | **remove** |
| `single_engine_native_v1` (new, R3c) | One engine only. Uses the shared v2 pose selector. Ranks ligands by the engine's native metric (GNINA `cnn_affinity` as the v2 ranking metric, higher is better; Vina/Smina/AD4 affinity, lower is better) as a within-target percentile, `consensus_status = single_engine`. Raises if more than one engine is present. | none (no cross-engine comparison, no geometry) | No. | No: v2 gives a degenerate ranking for one engine (every row 1.0) and is not a single-engine ranking. | **keep, purpose: explicit single-engine ranking** |

Unknown or removed mode names used to fall back silently to `dockbox_geometric` (`consensus.py` `normalize_consensus_mode`). That fallback is removed. Removed names now raise an error that names `consensus_rank_geometry_qc_v2`.

## Other ranking and geometry paths

| Path | Verdict | Reason |
|---|---|---|
| `geometric_consensus.compute_geometric_consensus` (used only by `dockbox_geometric`) | **removed** | Calls `_pairwise_rmsd`. |
| `geometric_consensus._pairwise_rmsd`, `_kabsch_rmsd` on element-ordered coordinates, `_centroid_sorted` | **retired (R3d)** | Forbidden methods. |
| `geometric_consensus.compute_geometric_consensus_v2`, `_strict_graph_pose`, `_extract_pose_signature` (used by redocking) | kept | Spec 031 graph path. |
| `consensus.select_representative_poses_v2` | kept as a thin wrapper over `pose_selection.select_best_pose_rows` (R4) | One selector. |
| `enhanced_rmsd_analyzer.py` (pairwise RMSD clustering, plots) | kept, **labelled `diagnostic_not_spec031`** (R3d) | Plots only; not Spec 031. |
| Best-pose writers using `_best_rows_by_group(..., "affinity_kcal_mol")` (min raw affinity) or `_select_best_pose_rows` (GNINA min `cnn_affinity`) | replaced by `pose_selection.select_best_pose_rows` (R4) | The GNINA worst-pose bug. See the fixture in the R4 test file. |
| `simplified_pipeline_impl.py` `mean_pairwise_rmsd` | out of scope | Best-pose RMSD summaries, not consensus. |
| `webui/targets.py` mode list | updated to v2 default; legacy names removed from the choice list | Same CLI option. |
| `best_engine_per_complex` (`reports/best_engine_per_complex.csv`, `multi_engine_pipeline_impl`) | **changed (Spec 036 follow-up)** | Picked one "best engine" per tag by the minimum raw affinity across engines (`_best_rows_by_group(best_by_engine, ["tag"], "affinity_kcal_mol")`). That compares Vina, Smina and AD4 kcal/mol values and GNINA's Vina-style affinity directly, which Spec 031 forbids. Now one row per tag with each engine's native selected-pose score in its own column (see below). |

### Follow-up change: `best_engine_per_complex.csv` (Spec 036, T008)

- **Before:** one row per tag. Columns are the winning row (`engine`, `affinity_kcal_mol`, `score_primary`, `pose_file`, ...).
  The winner was the lowest raw affinity across engines.
- **After:** one row per tag, with no winner.
  - `<engine>_affinity_kcal_mol`, `<engine>_pose`, `<engine>_replicate_id`, `<engine>_seed` and
    `<engine>_pose_reproducibility_status` for every engine present. GNINA also gives `gnina_cnn_score` and `gnina_cnn_affinity`.
  - `engines_present`, `engines_missing`, `comparison_status`.
  - `comparison_status = not_comparable_across_engines` when two or more engines are present for the tag.
    `single_engine` when only one is.
  - Kept columns: `tag`, `protein`, `ligand`, `site_id`, and the pair metadata and scope columns.
  - Dropped: `engine` and the raw `affinity_kcal_mol` winner columns.
- **Composite ranking:** not added to this file. The v2 composite (`consensus_rank_geometry_qc_v2`) is the only
  cross-engine ranking and is in `consensus_ranked_hits.csv`.
- **Changed (Spec 037 R3b, closes the T008 item):** `MultiEngineAnalysisPipeline._best_by_tag` no longer takes the minimum
  raw affinity across engines. See the next section.

## Spec 037 R3b: summaries use the consensus rank

Lead decision 6 (2026-10-10). The ranking is decided by `build_consensus_rankings` on the same scores that feed the
summary, so no DAG ordering is needed:

| Project | `ranking_method` | Source |
|---|---|---|
| Two or more engines | `consensus_rank_geometry_qc_v2` | `consensus_rank_global`, `consensus_rank_within_protein`, `consensus_score` |
| One engine | `single_engine_native_v1` | the same columns, from the single-engine native rule |
| A tag whose v2 case is incomplete (`consensus_status = consensus_incomplete`) | `ranking_unavailable_consensus_incomplete` | no rank and no score. It is never ordered by a raw affinity. |

- **Winner engine per tag.** A ranked tag keeps the pose of its `winner_engine` (the engine with the highest per-engine
  rank percentile in its complete case). The pose inside each engine is still chosen by `select_best_pose_rows`.
  `winner_engine` is written to the best-pose rows.
- **Changed outputs** (column names kept, `ranking_method` and the consensus columns added):
  - `_best_by_tag`: one row per ranked tag, ordered by rank. A tag without a complete v2 case keeps every engine row
    (one per engine, no winner invented), marked `ranking_unavailable_consensus_incomplete`. Feeds
    `affinity_analysis.best_poses` in the simplified bridge.
  - `_build_downstream_results`: `top_overall` (top 10 ranked rows), `best_per_protein` and `best_per_ligand` (the
    top-ranked row of each group), `protein_summary` and `ligand_summary` (ordered by `top_consensus_rank_global`;
    `best_affinity` is now the affinity of the top-ranked row, not a minimum across engines).
  - Simplified bridge `SimplifiedPostDockingPipeline._best_pose_per_protein_ligand` (picks the row by
    `consensus_rank_global` when present; the old `idxmin` on `vina_affinity` remains only for inputs without a rank).
  - Single-engine writer `_write_single_engine_reports`: `best_poses.csv` (ordered by rank), `protein_summary.csv`,
    and `summary.txt` (method, the pose selection key and the native metric, the top-ranked complex).
  - Comparative consensus call: one engine uses `single_engine_native_v1`, two or more use v2.
- **Not changed by R3b (remaining raw-affinity uses, for the Scientific Lead to decide):**
  - `top_pose_selector.build_top_pose_atlas`, default `top_pose_selection_policy = best_affinity`. It sorts by
    `affinity_kcal_mol` first across engines (`top_pose_per_ligand_global.csv`). The `best_consensus` and `hybrid`
    policies keep the rank as the primary key. Changing the default is a policy decision.
  - `_build_pair_competition` (rerun promotion: `is_target_engine_winner`, `affinity_advantage_kcal_mol`,
    `min_affinity_advantage`). It compares raw affinities across engines per tag.
  - `_write_cross_engine_visualizations`: the wide table `best_affinity` (min over engines) and `affinity_spread`,
    plus `_plot_cross_engine_disagreement` and `_plot_cross_engine_per_protein_batches`. These are diagnostics.
  - `consensus.py` `best_affinity_kcal_mol` (min raw across engines per tag). It is a tie-break in
    `select_rescoring_candidates` (after `consensus_score`) and a sort key in the polypharmacology summaries.
  - `_build_downstream_results` `summary_stats` (min, max and mean of `vina_affinity` across engines per complex) and the
    `mean_affinity` columns of `protein_summary` and `ligand_summary`. These are descriptive aggregates.
  - `_run_non_gnina_rmsd_bridge`, per-protein RMSD tables sorted by `vina_affinity`.
- **Single-engine GNINA label (Spec 037 R3c):** `summary.txt` names `cnn_score` as the pose selection key. The
  single-engine ligand ranking metric is `cnn_affinity`, which the v2 ranking already used. The two are reported on
  separate lines, and the choice between them is open for the Lead.

## Spec 031 conflict (for the Scientific Lead)

`consensus-v2.md` (Spec 031, approved 2026-08-21) says: "Preserve `dockbox_geometric` and existing rank modes under their legacy names for historical reproduction." Spec 036 D3 (2026-10-10) says to remove legacy modes that do not serve a valid purpose. This change follows Spec 036 D3 and removes the four legacy names from new runs. Historical reproduction stays possible by checking out the commit before this change. **Please confirm this reading.**
