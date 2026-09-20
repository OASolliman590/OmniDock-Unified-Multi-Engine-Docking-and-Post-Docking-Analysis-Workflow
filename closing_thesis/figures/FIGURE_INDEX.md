# Closing Thesis Visualization Scheme (Option 2 + SertHCl)

## Generated Figures
- `fig01_similarity_heatmap_triptych.(png|svg)`: mean-similarity heatmaps for biology/docking/SertHCl comparators.
- `fig02_delta_heatmaps.(png|svg)`: comparator delta maps from summary table.
- `fig03_enzyme_profiles_per_compound.(png|svg)`: per-compound enzyme profiles across comparators.
- `fig04_mcs_coverage_triptych.(png|svg)`: structural MCS coverage maps.
- `fig05_metric_distributions.(png|svg)`: fingerprint metric distribution by comparator type.
- `fig06_enzyme_comparator_means.(png|svg)`: enzyme-wise average similarity bars.

## Input Tables
- `pair_mean_similarity_option2_serthcl.csv`
- `similarity_metrics_option2_serthcl.csv`
- `mcs_coverage_option2_serthcl.csv`
- `enzyme_query_summary_option2_serthcl.csv`

## Relocated multiscale overlay

Moved from the repository root into this directory without changing bytes:

- `multiscale_overlapped_experimental_vs_triple_engine.(png|svg)`: overlapped experimental versus triple-engine comparison figure pair.

Original fixed-base blob hashes: png `1666b829fa38e5d44392bb52c6308f90c8b2711e`, svg `858514b5471c66e96a3e4fb64ae33d1d304889dd`.

## Added Comparator-Name + MCS Figures
- `fig00_comparator_name_key.(png|svg)`: enzyme-level comparator naming key (biology/docking/anchor).
- `fig01b_similarity_heatmap_named_comparators.(png|svg)`: similarity heatmaps with explicit comparator names in axis labels.
- `fig07_mcs_top_pairs_biology.(png|svg)`: RDKit MCS-highlighted top biology-comparator pair per enzyme.
- `fig08_mcs_top_pairs_docking.(png|svg)`: RDKit MCS-highlighted top docking-comparator pair per enzyme.
- `fig09_mcs_top_pairs_serthcl_anchor.(png|svg)`: RDKit MCS-highlighted top SertHCl-anchor pair per enzyme.

## Added Tables
- `comparator_name_key_option2_serthcl.csv`
- `enzyme_query_summary_with_comparator_names_option2_serthcl.csv`
- `fig07_mcs_top_pairs_biology_table.csv`
- `fig08_mcs_top_pairs_docking_table.csv`
- `fig09_mcs_top_pairs_serthcl_anchor_table.csv`
