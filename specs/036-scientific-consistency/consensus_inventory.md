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
- **Not changed (remaining, see tasks.md T008):** `MultiEngineAnalysisPipeline._best_by_tag` still takes the minimum
  raw affinity across engines. It feeds the legacy downstream tables (`affinity_analysis.best_poses`, `top_overall`,
  `best_per_protein`) in the simplified bridge.

## Spec 031 conflict (for the Scientific Lead)

`consensus-v2.md` (Spec 031, approved 2026-08-21) says: "Preserve `dockbox_geometric` and existing rank modes under their legacy names for historical reproduction." Spec 036 D3 (2026-10-10) says to remove legacy modes that do not serve a valid purpose. This change follows Spec 036 D3 and removes the four legacy names from new runs. Historical reproduction stays possible by checking out the commit before this change. **Please confirm this reading.**
