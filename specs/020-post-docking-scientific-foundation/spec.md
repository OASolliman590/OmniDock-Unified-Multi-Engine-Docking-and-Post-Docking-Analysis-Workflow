# Feature Specification: Post-Docking Scientific Foundation

**Feature Branch**: `020-post-docking-scientific-foundation`
**Created**: 2026-03-31
**Status**: Draft
**Input**: Deep audit of the post-docking analysis module — two rounds of code review, run log inspection, and DockBox source analysis — revealing systematic scientific design flaws, implementation bugs, and statistical errors across every major analysis module: broken RMSD semantics, unused CNN scores, split-brain pipeline architecture, inverted score normalization, bidirectional pairlist matching, absent redocking validation, no multiple testing correction, corrupt complex PDB construction, and miscalibrated consensus scoring diverging from the DockBox methodology it claims to follow.

---

## Problem Statement

The post-docking analysis module produces a large volume of outputs but the scientific validity of the most critical ranking and classification outputs is compromised by the following root causes, identified through code audit and run log inspection:

### Root Cause 1 — Split-brain pipeline architecture

Two independent pipelines exist with incompatible scientific logic:

| Pipeline | Entry point | Normalization | CNN scores | Consensus |
|---|---|---|---|---|
| Single-engine | `simplified_pipeline.py` | None — raw vina_affinity | Ignored | None |
| Multi-engine | `multi_engine_pipeline.py` | Per-engine rank percentile | Used in consensus | Weighted hybrid |

In practice, the overwhelming majority of runs use GNINA as a single engine and therefore always execute `simplified_pipeline.py`. All per-protein normalization, CNN score utilization, and consensus logic in `consensus.py` and `multi_engine_pipeline.py` is dead code for single-engine runs. Cross-target comparisons in the simplified pipeline are made on raw kcal/mol values, which is scientifically invalid because different proteins have structurally different affinity ranges.

### Root Cause 2 — Redocking validation never executed

The pipeline has 13 co-crystal ligand pairs in the pairlist (CEL, LIA, AXI, AQ4, 7GI, SRE, Y0R, ANP, RCX, STI, 748, BCT, and others). These are the only ground-truth reference points available to validate whether the GNINA docking setup reproduces known binding geometry. No step in either pipeline computes RMSD between the docked top pose and the co-crystal pose for these pairs. The docking setup is therefore never validated, and all downstream prospective rankings are built on unverified confidence.

### Root Cause 3 — GNINA CNN scores collected and discarded

GNINA produces three scores per pose: `vina_affinity`, `cnn_affinity`, `cnn_score`. The pipeline collects all three in `all_scores.csv` but ranks exclusively on `vina_affinity`. GNINA's scientific advantage over plain Vina is precisely the CNN-based scoring:

- `cnn_score`: CNN pose quality metric — should be used for **pose selection** (which of the 20 sampled poses to use)
- `cnn_affinity`: CNN-predicted binding affinity — should be used for **compound ranking** (more accurate than Vina for potency ordering)
- `vina_affinity`: Physics-based, fast — useful as secondary validation

Not using `cnn_affinity` and `cnn_score` means the CNN compute budget at docking time was wasted and the scientifically superior scoring is ignored.

### Root Cause 4 — RMSD per_protein scope computes cross-ligand RMSD

The RMSD module runs in three scopes: `per_complex`, `per_protein`, and `global`. The `per_protein` scope passes the top pose of each different ligand against the same receptor to the Kabsch RMSD calculator. Different ligands have different atom counts and atom types. The Kabsch algorithm cannot align molecules with different atom counts — it produces NaN for all cross-ligand pairs. The NaN filter then removes all valid indices, causing clustering to fail with "Not enough valid poses for N clusters" for every single protein in every run. All RMSD per_protein and global outputs are therefore invalid.

The `per_complex` scope (RMSD across multiple poses of the same ligand) is scientifically valid and works correctly.

### Root Cause 5 — Consensus implementation diverges from DockBox methodology

The current `consensus.py` implements a weighted rank aggregation formula:

```
consensus_score = mean_rank_pct * 0.55 + best_affinity_rank_pct * 0.30 +
                  (1 - agreement_fraction) * 0.10 + spread_norm * 0.05
```

The DockBox methodology (jp43/DockBox) that this module claims to follow does not use rank aggregation at all. DockBox uses:

1. **Independent best-pose selection per engine** using native scoring
2. **Geometric consensus validation**: `consensus = all(RMSD_ij ≤ 2.0 Å)` for all pairs of selected poses across engines

This is a pose-selection and geometric agreement approach, not a score aggregation approach. The difference is fundamental:
- Rank aggregation assumes scores are comparable across engines — they are not
- Geometric consensus asks "do multiple engines agree on where the ligand binds?" — engine-independent
- For single-engine runs, the current formula degenerates to `0.85 * rank_percentile`, making the entire consensus module a no-op that returns the sort order

The weights [0.55, 0.30, 0.10, 0.05] are not cited from any publication and are not calibrated against any benchmark dataset.

### Root Cause 6 — Polypharmacology scoring rewards promiscuity not selectivity

The polypharmacology scorer sums or counts binding events using absolute kcal/mol thresholds (-7.0 strong, -5.0 moderate) applied identically across all targets. This produces a score that correlates with molecular size and non-specific binding, not genuine polypharmacology. PNBA (p-nitrobenzoic acid, MW 167 Da) ranking #1 with "7 targets" is the expected result of this design — small flexible molecules fit into many sites via van der Waals contacts with no biological specificity.

Reference ligands (`exclude_reference_ligands: True`) are excluded from polypharmacology scoring. This is backwards: co-crystal ligands are the per-target gold standard and should anchor the scoring scale.

### Root Cause 7 — Hit classification has no per-target reference anchor

Hit classes (Strong/Moderate/Weak) are assigned by percentile rank within each target's compound pool. The pool includes co-crystal reference ligands without explicit curation. When a reference ligand docks poorly in redocking (e.g., BCT at -3.45 kcal/mol in CAII), the entire percentile distribution is compressed at the low end, inflating the classes of test compounds relative to a weak baseline.

ADMET status, which the spec and explainability fields promise as an independent classification dimension, is a hardcoded stub returning `"not_assessed"` everywhere.

Ligand efficiency (LE = affinity / heavy atom count) is not computed anywhere, resulting in size-biased rankings where large compounds are systematically overvalued.

### Root Cause 8 — MD receptor ensemble results are collapsed without treatment

3KK6 (COX1) and 3LN1 (COX2) each have 6 conformations (crystal + 5 MD snapshots at 600ps, 700ps, 800ps, 900ps, 1000ps). The pairlist maps all to the same protein name. The hierarchical analyzer produces duplicate entries for COX1 and COX2 in the "Best Ligand per Protein" table — once from the crystal structure, once from the MD snapshot aggregate — without distinguishing them. The scientific value of the MD ensemble (receptor flexibility sampling) is lost.

### Root Cause 9 — Score normalization is mathematically inverted in `consensus.py`

The min-max normalization in `consensus.py` computes:
```python
(values - min) / span  →  [0, 1]
```
In docking, **lower affinity = better binding**. The most negative score (strongest binder) maps to `0` and the weakest binder maps to `1`. Every downstream consensus formula then treats `0` as worst and `1` as best — inverting the entire ranking. The z-score normalization has the same problem. The per-engine rank percentile accidentally escapes this because it assigns rank ordinal, not a directional value — but only if the rank direction is correctly set. This inversion means every "Strong Binder" label produced by multi-engine consensus mode may identify the weakest compounds in the pool.

### Root Cause 10 — Pairlist filename matching uses bidirectional substring, causing false assignments

`generate_scores_csv.py` matches log filenames to pairlist rows with:
```python
if pattern in filename or filename in pattern
```
This bidirectional check causes false positives. If the pairlist contains `"3KK6_site_1_CEL"` and a log file is named `"3KK6_site_1_CEL_A_701"`, both halves of the condition match. In datasets with overlapping ligand names (e.g., `MS1` and `MS1_modified`, or PDB codes like `1S1` inside `1S1B`), log files are silently assigned to the wrong protein/ligand mapping. The resulting `all_scores.csv` carries wrong labels with no warning or error.

### Root Cause 11 — Correlation p-values reported without multiple testing correction

`correlation_analyzer.py` computes 6 p-values (3 Pearson + 3 Spearman) at α=0.05 with no Bonferroni or Benjamini-Hochberg correction. Expected false positives from 6 uncorrected tests = 0.3 per analysis run. With per-target sample sizes of n=3–5 (common after exclusion filtering on a 37-ligand dataset), statistical power is near zero. Most "significant" correlations in the output are false positives by construction.

Additionally, Pearson and Spearman are called without enforcing a minimum sample size. With n=2, Pearson always returns r=±1.0, p≈0 — mathematically defined but scientifically meaningless. No guard exists.

### Root Cause 12 — Single-engine QC downgrade fires universally

The QC downgrade rule `agreement_count < 2 → warn_low_engine_agreement → downgrade Strong to Moderate` was designed for multi-engine disagreement. For single-engine runs, `agreement_count = 1 < 2` always. Every compound in every GNINA-only run is therefore downgraded from whatever classification it earned — regardless of pose quality. This is not a QC signal; it is noise applied uniformly.

The downgrade chain then cascades: `warn_low_engine_agreement` drops Strong→Moderate, and a separate `fail_positive_affinity` check can then drop Moderate→Weak. The same final class can be reached via two different rule paths, making the classification semantically inconsistent.

### Root Cause 13 — Complex PDB files lose ligand residue name

`pose_extractor.py` overwrites every ligand atom's residue name with `"UNK"`:
```python
new_line = new_line[:17] + "UNK" + new_line[20:]
```
ProLIF, LigPlot+, and PyMOL all use the residue name to identify the ligand molecule in a complex. When all ligands are `"UNK"`, tools that distinguish multiple ligands or that rely on residue name for interaction detection fail silently or produce incomplete maps. This is the direct cause of the `'ligand'` KeyError warnings in ProLIF for certain complexes.

### Root Cause 14 — `structure_processor.py` is an empty stub

`structure_processor.py` is imported and appears in the module tree, but all its core functions (`split_complexes`, `extract_apo_proteins`, validate/merge logic) only create directories — they contain no structure processing code. The actual complex merging logic is scattered in `pose_extractor.py`. The stub file is maintenance debt with no function, and its presence creates false confidence that complex structure handling has a dedicated, tested module.

### Root Cause 15 — PDBQT fallback parser misdetects chain ID

The fallback PDBQT parser in `pose_extractor.py` infers chain ID with:
```python
if len(parts[4]) == 1 and parts[4].isalpha():  # → chain ID
```
When the residue number is single-digit (residue 1–9, common for small molecules from GNINA) and the next column is a single letter, the parser misinterprets the residue number as a chain ID. Atom coordinates are then written to the wrong chain in the complex PDB. This silently corrupts the 3D structure of small-molecule complexes.

### Root Cause 16 — Pose tie-breaking is non-reproducible across machines

`pose_extractor.py` selects the best pose by iterating `all_scores.csv` rows sequentially. When two poses have identical `vina_affinity`, whichever row appears later wins. The CSV row order depends on the filesystem glob order used when log files were processed — which varies across operating systems and filesystems. The same docking results can produce a different "best pose" on macOS vs. Linux, breaking reproducibility.

---

## Goals

1. Implement redocking RMSD validation as a mandatory pre-analysis gate with pass/fail reporting
2. Unify the two-pipeline architecture so single-engine runs receive per-protein normalization
3. Activate CNN score utilization for GNINA runs: `cnn_score` for pose selection, `cnn_affinity` for ranking
4. Replace RMSD per_protein scope with scientifically valid intra-run pose diversity RMSD
5. Implement DockBox-style geometric consensus for multi-engine runs replacing rank aggregation
6. Replace raw affinity polypharmacology with per-target normalized polypharmacology
7. Implement reference-anchored hit classification using co-crystal ligands as per-target baselines
8. Compute ligand efficiency alongside affinity in all ranking outputs
9. Implement MD ensemble aggregation for conformational receptor sets
10. Fix inverted score normalization in consensus.py so ranking direction is correct
11. Fix bidirectional pairlist filename matching to prevent false label assignments
12. Add multiple testing correction and minimum sample size enforcement to correlation analysis
13. Fix single-engine QC downgrade to not universally penalize single-engine runs
14. Fix QC downgrade cascade to be atomic, not sequential
15. Preserve ligand residue names in complex PDB files
16. Remove `structure_processor.py` stub and consolidate complex creation logic
17. Fix PDBQT chain ID detection bug in fallback parser
18. Make best-pose tie-breaking deterministic across machines and filesystems

## Non-Goals

- Changing the docking preparation or execution pipeline
- Adding new visualization types beyond what already exists
- Replacing ProLIF, LigPlot+, or PandaMap integrations
- Implementing ADMET prediction (stub removal is in scope; new ADMET engine is not)
- Changing the folder topology or session structure

---

## Architectural Contract — Pipeline Unification

### Current state (broken)

```
Entry point
    │
    ├── 1 engine ──→ simplified_pipeline.py
    │                   raw vina_affinity ranking
    │                   no per-protein normalization
    │                   CNN scores ignored
    │                   no consensus
    │                   no redocking validation
    │
    └── N engines ──→ multi_engine_pipeline.py
                        per-engine rank percentile
                        CNN scores used
                        weighted rank consensus
                        top_pose_selector
```

Two pipelines with incompatible scientific logic. Adding a feature to one does not propagate to the other. The split is the root cause of every "scientifically weaker path" issue identified in the audit.

### Target state (unified)

```
Entry point (1..N engines)
    │
    └──→ unified_pipeline.py
              │
              ├── [Step 1]  Find input files
              ├── [Step 2]  Parse scores → all_scores.csv
              │               vina_affinity + cnn_affinity + cnn_score (GNINA)
              │               or engine native score (other engines)
              │
              ├── [Step 3]  Pose selection
              │               cnn_score → best pose (GNINA)
              │               native score → best pose (other engines)
              │
              ├── [Step 4]  Match poses to receptors
              ├── [Step 5]  Create complex PDB files
              │
              ├── [Step 6]  Redocking validation (if co-crystal pairs present)
              │               RMSD docked pose vs co-crystal pose
              │               pass / warn / fail per pair
              │
              ├── [Step 7]  Per-protein normalization (always, 1 or N engines)
              │               rank percentile within each protein
              │               applied to primary scoring metric
              │
              ├── [Step 8]  Binding affinity analysis
              │               Level 1: best pose per pair
              │               Level 2: best ligands per protein (normalized)
              │               Level 3: cross-protein (normalized ranks, not raw)
              │               Level 4: comparative / redocking
              │
              ├── [Step 9]  MD ensemble aggregation
              │               detect multi-snapshot proteins from pairlist
              │               ensemble_min / mean / std per biological protein
              │
              ├── [Step 10] Consensus
              │               N=1: single_engine annotation, no consensus computation
              │               N≥2: DockBox geometric consensus (RMSD ≤ 2.0 Å)
              │               optional: rank_aggregation mode (preserved, not default)
              │
              ├── [Step 11] Ligand efficiency
              │               LE = |affinity| / HAC from SDF
              │
              ├── [Step 12] Reference-anchored hit classification
              │               co-crystal anchor if available + redocking pass
              │               percentile fallback otherwise
              │
              ├── [Step 13] Polypharmacology (normalized)
              │               per-target rank percentile contribution
              │               aggregate by count of targets above threshold
              │
              ├── [Step 14] RMSD analysis (per_complex only)
              │               intra-run pose diversity per complex
              │               redocking RMSD handled in Step 6
              │
              ├── [Step 15] Top pose atlas (top_pose_selector)
              ├── [Step 16] Reports + VALIDATION_SUMMARY.md
              ├── [Step 17] Visualizations (ProLIF, LigPlot+, py3Dmol, PandaMap)
              └── [Step 18] Consolidate output artifacts
```

### Migration rules

| Old file | Fate |
|---|---|
| `simplified_pipeline.py` | Retired. Replaced by `unified_pipeline.py`. |
| `multi_engine_pipeline.py` | Retired. Logic absorbed into `unified_pipeline.py`. |
| `structure_processor.py` | Retired. Stub removed. Complex creation logic consolidated in `pose_extractor.py`. |
| `consensus.py` | Retained. Normalization direction fixed. Geometric consensus added as default. Rank aggregation preserved as optional. Single-engine QC handling added. Atomic downgrade chain. |
| `hierarchical_analyzer.py` | Retained. Called from Step 8 with normalized scores. Co-crystal typo fixed (`"compartive"` → `"comparative"`). |
| `enhanced_rmsd_analyzer.py` | Retained. `per_protein` and `global` scopes removed. `per_complex` unchanged. |
| `top_pose_selector.py` | Retained. Tie-breaking made deterministic. Single-engine confidence annotation added. |
| `pose_extractor.py` | Retained. Residue name preservation added. PDBQT chain ID fix. Deterministic tie-breaking. |
| `generate_scores_csv.py` | Retained. Bidirectional matching replaced. Empty log reporting added. |
| `correlation_analyzer.py` | Retained. BH correction added. Minimum N enforced. Rank normalization replaces min-max. |
| `binding_affinity_analyzer.py` | Retained. Extended with LE and MW columns. |
| `report_generator.py` | Retained. String escape bug fixed. |
| All visualization modules | Retained unchanged. |

### Compatibility rules

- All existing CLI flags and entry points that currently invoke `simplified_pipeline` or `multi_engine_pipeline` MUST be re-routed to `unified_pipeline` transparently. No user-facing flag changes.
- The `simplified_cli.py` and `cli.py` entry points are retained as wrappers. Internally they construct a `unified_pipeline` with `engines=["gnina"]` or the user-supplied engine list.
- Session output paths and folder structure are unchanged.
- All existing CSV output filenames are unchanged so downstream consumers are not broken.

---

## User Stories

### Story 1 — Redocking Validation Gate (P0)

As a computational chemist, I need to know whether my GNINA docking setup reproduces known co-crystal binding geometries before I can trust any prospective ranking.

**Acceptance criteria:**
1. Given a pairlist with co-crystal ligand entries, when the analysis runs, then a redocking RMSD table is generated for all co-crystal pairs showing docked-pose vs. co-crystal-pose RMSD.
2. Given redocking RMSDs, when the validation report is generated, then each pair is classified as pass (RMSD ≤ 2.0 Å), warn (2.0 < RMSD ≤ 3.5 Å), or fail (RMSD > 3.5 Å).
3. Given one or more redocking failures, when the analysis continues, then a validation warning is emitted in the session summary indicating prospective results should be interpreted cautiously.
4. Given no co-crystal pairs in the pairlist, when analysis runs, then a clear notice is emitted that docking validation could not be performed.

---

### Story 2 — Unified Normalization for Single-Engine Runs (P0)

As an analyst running GNINA-only docking, I need cross-target ligand comparisons to use per-target normalized scores rather than raw kcal/mol so that ligand rankings are not biased by the intrinsic affinity range of each target.

**Acceptance criteria:**
1. Given a single-engine GNINA run, when binding affinity analysis runs, then cross-target rankings use per-protein rank percentile normalization.
2. Given per-protein rank percentiles, when the cross-protein matrix is generated, then raw affinities and normalized ranks are both present as separate columns.
3. Given per-protein normalization, when the top ligands globally are reported, then the ranking is based on normalized scores, not mean raw kcal/mol.

---

### Story 3 — CNN Score Utilization for GNINA (P0)

As a GNINA user, I want the CNN-predicted affinity to be used for compound ranking and the CNN pose quality score to be used for pose selection so that GNINA's deep learning advantage is not discarded after docking.

**Acceptance criteria:**
1. Given GNINA output with `cnn_score` and `cnn_affinity` columns, when pose selection runs, then the pose with the highest `cnn_score` per complex is selected as the representative pose.
2. Given pose selection by `cnn_score`, when per-target rankings are generated, then `cnn_affinity` is used as the primary ranking metric with `vina_affinity` as secondary.
3. Given outputs, when ranking tables are written, then `cnn_affinity`, `cnn_score`, and `vina_affinity` are all present as columns with clear column headers indicating their role.
4. Given a run where CNN scores are absent (non-GNINA engine), when analysis runs, then `vina_affinity` is used as the primary ranking metric without error.

---

### Story 4 — Corrected RMSD Semantics (P0)

As a scientist, I need RMSD analysis to compute meaningful conformational metrics rather than cross-ligand atom-count mismatches that produce all-NaN matrices.

**Acceptance criteria:**
1. Given per_complex scope, when RMSD runs, then RMSD is computed across all N poses of the same ligand against the same receptor (intra-run diversity) — unchanged, already valid.
2. Given per_protein scope, when RMSD runs, then RMSD is no longer computed — this scope is removed. Users are directed to per_complex scope or redocking validation for per-target comparisons.
3. Given global scope, when RMSD runs, then RMSD is no longer computed — this scope is removed.
4. Given redocking pairs in the pairlist, when redocking RMSD runs, then RMSD between docked top pose and co-crystal pose is computed for each pair using the same Kabsch heavy-atom RMSD as per_complex.
5. Given RMSD clustering, when it runs on per_complex results, then clustering uses the per-complex RMSD matrix — same molecule, same receptor — so NaN rows cannot appear due to atom count mismatch.

---

### Story 5 — DockBox-Style Geometric Consensus for Multi-Engine (P1)

As a user running multiple docking engines, I want consensus to be based on whether engines agree on the binding pose location, not on a weighted rank formula with arbitrary uncalibrated weights.

**Acceptance criteria:**
1. Given poses from N engines for the same ligand-protein pair, when consensus is computed, then the best pose per engine (by that engine's primary metric) is selected independently.
2. Given N selected poses, when geometric consensus is evaluated, then pairwise RMSD between all selected poses is computed using heavy-atom Kabsch RMSD.
3. Given pairwise RMSDs, when consensus assignment runs, then `geometric_consensus = True` if all pairwise RMSDs ≤ 2.0 Å, `False` otherwise. The cutoff is configurable.
4. Given geometric consensus result, when hit classification runs, then compounds with `geometric_consensus = True` receive a confidence boost; compounds with `geometric_consensus = False` receive the existing QC downgrade.
5. Given a single-engine run, when consensus runs, then `geometric_consensus = True` by definition (one pose, no inter-engine disagreement possible) and the field is present but annotated as `single_engine`.
6. Given the existing rank-aggregation formula in `consensus.py`, it is preserved as an optional `rank_aggregation` mode selectable at runtime but is not the default.

---

### Story 6 — Per-Target Normalized Polypharmacology (P1)

As a scientist, I need the polypharmacology score to reflect relative potency per target rather than the sum of absolute affinities so that small promiscuous molecules are not falsely promoted over genuinely selective multi-target compounds.

**Acceptance criteria:**
1. Given binding results for a ligand across multiple targets, when polypharmacology scoring runs, then the score for each target contribution is the percentile rank of that ligand within that target's compound distribution (same normalization as cross-target ranking).
2. Given per-target normalized contributions, when a polypharmacology aggregate score is computed, then it is the count of targets where the ligand exceeds a configurable normalized rank threshold (default: top 35% per target).
3. Given reference ligands in the dataset, when polypharmacology runs, then co-crystal ligands are included in the per-target pool (not excluded) so they anchor the percentile scale.
4. Given the top polypharmacology table, when it is written, then raw affinity per target, normalized rank per target, and the aggregate score are all present as separate columns.

---

### Story 7 — Reference-Anchored Hit Classification (P1)

As a scientist, I need hit classification to be calibrated against the co-crystal reference ligand for each target, not just a floating percentile of the test compound pool.

**Acceptance criteria:**
1. Given a co-crystal ligand identified for a target (from pairlist), when hit classification runs, then the reference ligand's affinity is used as the per-target calibration anchor.
2. Given a reference anchor, when strong/moderate/weak thresholds are applied, then:
   - Strong: ligand affinity ≤ reference affinity (equal or better than co-crystal)
   - Moderate: reference affinity < ligand affinity ≤ reference affinity + 2.0 kcal/mol
   - Weak: ligand affinity > reference affinity + 2.0 kcal/mol
3. Given no reference ligand for a target, when classification runs, then the existing percentile fallback (top 10% = strong, top 35% = moderate) is used with a flag `classification_basis = "internal_percentile"`.
4. Given a reference ligand, when classification outputs are written, then `classification_basis = "reference_anchored"` and `reference_ligand` and `reference_affinity_kcal_mol` are included as columns.
5. Given CNN scores available, when classification runs, then both `vina_affinity` and `cnn_affinity` thresholds are evaluated and reported separately.

---

### Story 8 — Ligand Efficiency in Rankings (P1)

As a medicinal chemist, I need ligand efficiency (LE) computed alongside raw affinity in all ranking outputs so that large compounds are not systematically overvalued relative to efficient smaller compounds.

**Acceptance criteria:**
1. Given affinity and ligand SMILES or SDF, when LE is computed, then `LE = |affinity_kcal_mol| / heavy_atom_count`.
2. Given LE, when ranking tables are written, then `ligand_efficiency`, `heavy_atom_count`, and `molecular_weight` (if available) are present as additional columns.
3. Given LE, when it is not computable (missing atom count), then `ligand_efficiency = NaN` without error.
4. Given ranking outputs, when a user sorts by LE, then the sort order differs from affinity sort — this discrepancy is documented as intentional (LE-based lead optimization vs. raw potency screening).

---

### Story 9 — MD Ensemble Aggregation (P1)

As a scientist using MD receptor snapshots for ensemble docking, I need results across conformational snapshots to be aggregated into a single per-protein result that preserves the ensemble signal rather than producing duplicate protein entries.

**Acceptance criteria:**
1. Given a pairlist mapping multiple receptor files to the same protein name (e.g., 3KK6, 3KK6_receptor_600ps ... all → COX1), when analysis runs, then a snapshot-level table is generated showing affinity per snapshot per ligand.
2. Given snapshot-level results, when ensemble aggregation runs, then an ensemble-level table is generated with the following metrics per protein per ligand: `ensemble_min_affinity` (best across snapshots), `ensemble_mean_affinity`, `ensemble_std_affinity`, `n_snapshots_where_best_binder`.
3. Given ensemble aggregation, when the "Best Ligand per Protein" table is generated, then each protein appears exactly once, using ensemble_min_affinity as the primary sort key.
4. Given snapshot-level results, they are retained as a separate `snapshot_affinity_matrix.csv` for users who want per-conformation detail.
5. Given a protein with only one receptor file, when analysis runs, then it is treated identically — no ensemble aggregation, snapshot table and ensemble table are the same row.

---

### Story 10 — Validation Summary and Analysis Trust Signal (P0)

As a scientist, I need a single human-readable validation report at the top of each analysis session that tells me how much to trust the downstream outputs.

**Acceptance criteria:**
1. Given an analysis run, when it completes, then a `VALIDATION_SUMMARY.md` is written to the session root with the following sections:
   - Docking setup validation: redocking RMSD pass/warn/fail per co-crystal pair
   - Score distribution sanity: fraction of complexes with positive affinity, missing poses
   - Normalization method applied: raw vs per-protein rank percentile
   - Scoring metric used: vina_affinity or cnn_affinity (for GNINA)
   - Ensemble treatment: which proteins were aggregated
   - Overall confidence signal: `HIGH` (≥80% redocking pass, normalization applied), `MEDIUM`, or `LOW` (no validation possible or >50% redocking fail)
2. Given a `LOW` confidence signal, when START_HERE.md is generated, then it includes a prominent warning before the results navigation section.

### Story 11 — Correct Score Normalization Direction (P0)

As a scientist using consensus rankings, I need score normalization to preserve the docking convention that lower affinity = better binding so that consensus scores identify strong binders, not weak ones.

**Acceptance criteria:**
1. Given vina_affinity values where lower = better, when min-max normalization is applied, then the most negative value maps to `1.0` (best) and the least negative maps to `0.0` (worst).
2. Given z-score normalization, when applied to affinities, then the result is negated so that more negative z-scores (stronger binders) yield higher normalized values.
3. Given rank percentile normalization, when applied, then rank 1 (best affinity) yields percentile 1.0 and rank N (worst) yields 0.0 (ascending rank = descending score value).
4. Given normalized scores, when a consensus formula weights them, then a ligand with the most negative raw affinity produces the highest consensus score — verified on a synthetic dataset where the expected winner is known.

---

### Story 12 — Directional Pairlist Filename Matching (P0)

As a user with a complex naming convention in my docking project, I need log files to be matched to pairlist entries without false positives so that protein and ligand labels in all_scores.csv are correct.

**Acceptance criteria:**
1. Given a pairlist pattern `"3KK6_site_1_CEL"` and a log file `"3KK6_site_1_CEL_A_701"`, when matching runs, then the log is NOT matched to this pattern (directional: pattern must fully match the log stem, not be a substring of it).
2. Given a pairlist pattern `"3KK6_site_1_CEL_A_701"` and a log file `"3KK6_site_1_CEL_A_701.log"`, when matching runs, then the log IS matched.
3. Given overlapping patterns like `"MS1"` and `"MS1_modified"`, when a log named `"receptor_site_1_MS1.log"` is matched, then it matches `"MS1"` exactly and does NOT match `"MS1_modified"`.
4. Given a log file that matches no pairlist pattern, when parsing runs, then the file is processed with filename-derived labels and a warning is emitted listing unmatched files.

---

### Story 13 — Statistically Rigorous Correlation Analysis (P1)

As a scientist reviewing cross-engine or score-biology correlations, I need reported p-values to reflect corrected significance levels and I need correlations suppressed when sample size is insufficient to support them.

**Acceptance criteria:**
1. Given N correlation tests computed per analysis run, when p-values are reported, then Benjamini-Hochberg FDR correction is applied and both raw p-value and adjusted q-value are present as columns.
2. Given fewer than 5 paired non-null samples for a correlation, when correlation is requested, then no correlation is computed and the result row is annotated `insufficient_n` with the actual n reported.
3. Given a per-target correlation with n between 5 and 10, when the result is written, then a `low_power` flag is set alongside the correlation value.
4. Given correlation normalization across engines, when it is applied, then rank-based normalization (not min-max) is used so that wide-range engines do not dominate narrow-range ones.

---

### Story 14 — Single-Engine QC Behavior (P0)

As a GNINA-only user, I need QC flags and hit class downgrades to reflect actual pose quality signals, not universally fire because I used one engine instead of two.

**Acceptance criteria:**
1. Given a single-engine run, when QC status is assigned, then `warn_low_engine_agreement` is NOT emitted — this warning is only meaningful when 2+ engines exist.
2. Given a single-engine run, when `agreement_fraction` is computed, then it is set to `NaN` or `"single_engine"` — not `1.0` — to explicitly signal that inter-engine agreement was not evaluated.
3. Given a hit classification, when the QC downgrade chain runs, then it is evaluated atomically: all applicable QC flags are collected first, then a single classification decision is made. Strong→Moderate and then →Weak on the same compound in the same run is not possible.
4. Given a `fail_positive_affinity` QC flag, when classification runs, then the class is forced to Weak regardless of percentile rank, as before — this rule is retained.

---

### Story 15 — Complex PDB Structure Integrity (P1)

As a scientist using the complex PDB files for visualization and interaction analysis downstream, I need complex structures to preserve ligand identity and have correct chain assignments so that ProLIF, LigPlot+, and PyMOL can identify the ligand correctly.

**Acceptance criteria:**
1. Given a ligand with residue name `CEL` in the source SDF, when the complex PDB is written, then the HETATM records use `CEL` as the residue name — not `UNK`.
2. Given a PDBQT file with residue number `5` and chain ID `A`, when the fallback PDBQT parser runs, then residue number `5` is parsed as the residue number and `A` is parsed as the chain ID — not one misidentified as the other.
3. Given a merged complex PDB with receptor ATOM records and ligand HETATM records, when atom serial numbers are written, then they are renumbered sequentially from 1 with no conflicts between receptor and ligand atoms.
4. Given `structure_processor.py`, when the module is loaded, then it contains no stub functions that only create directories. All complex creation logic is consolidated in `pose_extractor.py` or a dedicated `complex_builder.py`.

---

### Story 16 — Deterministic Pose Selection (P1)

As a scientist running the same analysis on two different machines, I need the same docking output to produce the same best-pose selection so that results are reproducible regardless of filesystem or operating system.

**Acceptance criteria:**
1. Given two poses with identical `vina_affinity` for the same complex, when best pose is selected, then the tie is broken by `cnn_score` (higher wins) if available, then by `pose_index` (lower wins), in that order — independent of CSV row order.
2. Given the same `all_scores.csv` run on macOS and Linux, when best pose selection runs, then the same pose is selected for every complex.
3. Given a complex where all 20 poses have identical scores (degenerate case), when best pose is selected, then pose index 1 (first pose) is always selected.

---

## Scientific Design Decisions

### Decision 1 — DockBox geometric consensus over rank aggregation

**Chosen approach**: For multi-engine runs, consensus is the set of ligands where all engines independently select a pose within 2.0 Å RMSD of each other. This is the jp43/DockBox approach.

**Why not rank aggregation**: Rank aggregation assumes score scales are comparable across engines — they are not. AutoDock4, Vina, GNINA, and DOCK 6 have fundamentally different scoring functions with different units and distributions. Converting to percentile rank within each engine partially corrects this, but then aggregating ranks introduces implicit assumptions about how to weight each engine. The geometric approach avoids all of this: it asks only "do engines agree on where the ligand is?" which is engine-independent.

**Why preserve rank_aggregation as an option**: Some users may want a continuous score for ranking rather than a binary consensus flag. The existing formula is preserved as `consensus_mode = rank_aggregation` for backward compatibility and for cases where pose structures are not available.

### Decision 2 — cnn_affinity as primary ranking metric for GNINA

**Chosen approach**: `cnn_score` selects the representative pose from the 20 sampled; `cnn_affinity` ranks compounds per target.

**Why not vina_affinity**: GNINA's published validation (McNutt et al. 2021, J. Cheminformatics) shows CNN scoring functions consistently outperform the Vina scoring function on pose selection accuracy and affinity rank correlation with experimental data, particularly for congeneric series. Using vina_affinity with GNINA discards the primary reason to choose GNINA over Vina.

**Why keep vina_affinity**: As a secondary metric and for users who need compatibility with Vina-based comparisons or who distrust CNN scoring for their specific target class.

### Decision 3 — 2.0 Å RMSD redocking pass threshold

**Chosen approach**: RMSD ≤ 2.0 Å = pass, 2.0–3.5 Å = warn, > 3.5 Å = fail.

**Why 2.0 Å**: The 2.0 Å heavy-atom RMSD threshold is the standard criterion in the docking literature for acceptable pose reproduction (Pagadala et al. 2017, Biophysical Reviews; Spitzer & Jain 2012, Drug Discovery Today). It is also the default cutoff in DockBox.

**Why 3.5 Å for warn rather than hard fail**: For large flexible ligands (high torsion count) or targets with flat binding sites, 2.0–3.5 Å may still represent the correct binding region with a slightly different conformation. A warn allows downstream review rather than hard blocking.

### Decision 4 — Reference-anchored thresholds in kcal/mol deltas

**Chosen approach**: Strong = affinity ≤ reference, Moderate = reference to reference + 2.0, Weak = worse than reference + 2.0 kcal/mol.

**Why 2.0 kcal/mol margin**: A 2.0 kcal/mol difference in docking score corresponds to approximately 10-fold difference in estimated binding constant (ΔΔG = RT ln(10) ≈ 1.36 kcal/mol per order at 298K). This is a physically meaningful margin commonly used to distinguish meaningful binding differences from docking noise.

**Why not tighter (e.g., 1.0 kcal/mol)**: Docking scoring functions have typical errors of 1–2 kcal/mol compared to experimental values. A 1.0 kcal/mol threshold would classify most differences as noise.

### Decision 5 — Ligand efficiency with heavy atom count from SDF

**Chosen approach**: LE = |affinity| / HAC where HAC (heavy atom count) is computed from the SDF molecule.

**Why**: LE normalizes for molecular complexity, a fundamental metric in fragment-based drug discovery and lead optimization. A compound with LE ≥ 0.3 kcal/mol/HA is generally considered efficient; LE < 0.2 suggests size-driven binding.

---

## Requirements

### FR-000 — Pipeline unification
`simplified_pipeline.py` and `multi_engine_pipeline.py` MUST be retired and replaced by a single `unified_pipeline.py`. All single-engine and multi-engine runs MUST flow through `unified_pipeline`. Engine count (1 or N) determines which consensus mode applies — it does NOT determine which pipeline is invoked. Existing CLI entry points MUST transparently re-route to `unified_pipeline` with no user-facing flag changes. No scientific logic may remain exclusively in either retired file.

### FR-001 — Redocking RMSD validation
System MUST identify co-crystal ligand pairs from pairlist `cocrystal_ligand_pdb` or `cocrystal_ligand_name` columns and compute heavy-atom Kabsch RMSD between the docked top pose and the co-crystal SDF entry for each pair.

### FR-002 — Redocking RMSD classification
System MUST classify each redocking result as `pass` (≤2.0 Å), `warn` (2.0–3.5 Å), or `fail` (>3.5 Å) and report pass rate and per-pair detail.

### FR-003 — Per-protein rank normalization in simplified pipeline
System MUST apply per-protein rank percentile normalization to vina_affinity and cnn_affinity before generating any cross-protein comparison tables, regardless of whether one or multiple engines are used.

### FR-004 — cnn_score for pose selection
System MUST select the pose with maximum `cnn_score` as the representative pose for each complex when `cnn_score` is available.

### FR-005 — cnn_affinity as primary ranking metric
System MUST use `cnn_affinity` as the primary sort key for per-target rankings when available, with `vina_affinity` as secondary. Both columns MUST be present in all ranking outputs.

### FR-006 — RMSD per_protein and global scopes removed
System MUST NOT compute RMSD between top poses of different ligands against the same receptor. The `per_protein` and `global` RMSD scopes MUST be removed. `per_complex` scope is unchanged.

### FR-007 — DockBox geometric consensus for multi-engine
System MUST compute geometric consensus as `all(pairwise_RMSD_ij ≤ cutoff)` across engine-selected poses when 2+ engines provide results for the same complex.

### FR-008 — Single-engine consensus annotation
System MUST annotate `geometric_consensus = single_engine` when only one engine provides results, distinguishing this from multi-engine true/false consensus.

### FR-009 — Per-target normalized polypharmacology
System MUST compute per-target rank percentile contribution for polypharmacology scoring and aggregate by counting targets where normalized rank exceeds threshold.

### FR-010 — Reference ligands included in polypharmacology pool
System MUST include co-crystal reference ligands in the per-target compound pool for polypharmacology normalization. `exclude_reference_ligands` MUST default to `False`.

### FR-011 — Reference-anchored hit classification
System MUST use co-crystal ligand affinity as the per-target classification anchor when a co-crystal pair exists for that target.

### FR-012 — Percentile fallback when no reference
System MUST fall back to top-10%/top-35% percentile classification when no co-crystal reference exists for a target, annotating `classification_basis = "internal_percentile"`.

### FR-013 — Ligand efficiency computation
System MUST compute `ligand_efficiency = |affinity_kcal_mol| / heavy_atom_count` and include it in all per-pair ranking tables.

### FR-014 — MD ensemble aggregation
System MUST detect when multiple receptor files map to the same protein name and produce ensemble aggregate metrics (`ensemble_min`, `ensemble_mean`, `ensemble_std`, `n_snapshots_best`) alongside per-snapshot detail.

### FR-015 — One protein per row in top-ligand tables
System MUST ensure each biological protein appears exactly once in "Best Ligand per Protein" tables, using ensemble aggregation for proteins with multiple receptor snapshots.

### FR-016 — VALIDATION_SUMMARY.md
System MUST write a `VALIDATION_SUMMARY.md` to the session root summarizing redocking validation, normalization method, scoring metric, and overall confidence signal.

### FR-017 — ADMET stub removal
System MUST remove or clearly mark `admet_status = "not_assessed"` placeholder. Fields MUST be absent from outputs unless an ADMET source is available. The spec explicitly defers ADMET computation to a future feature.

### FR-018 — Correct normalization direction
All normalization methods in `consensus.py` MUST map the most favorable docking score (most negative affinity) to the highest normalized value (1.0 for min-max, highest rank for percentile). A synthetic test with known winner MUST be included to verify direction.

### FR-019 — Directional pairlist matching only
`generate_scores_csv.py` MUST match log filenames to pairlist patterns using exact stem equality or strict prefix matching only. Bidirectional substring matching (`pattern in filename or filename in pattern`) MUST be replaced. Unmatched log files MUST be reported in a `unmatched_log_files.txt` in the session output.

### FR-020 — Multiple testing correction in correlation analysis
`correlation_analyzer.py` MUST apply Benjamini-Hochberg FDR correction to all reported p-values. Both raw `p_value` and corrected `q_value` MUST be present in output tables. Correlation MUST be suppressed (result annotated `insufficient_n`) when sample size < 5.

### FR-021 — Single-engine QC explicit handling
`consensus.py` MUST NOT emit `warn_low_engine_agreement` for single-engine runs. `agreement_fraction` MUST be set to `NaN` when engine count is 1. A dedicated `single_engine` flag column MUST be added to consensus outputs.

### FR-022 — Atomic QC downgrade
The QC downgrade chain MUST be evaluated atomically: collect all applicable QC flags first, then apply a single classification decision. A compound MUST NOT be downgraded twice in the same run via sequential rule application.

### FR-023 — Ligand residue name preservation
`pose_extractor.py` MUST preserve the ligand's original residue name from the source SDF/PDBQT when writing HETATM records to the complex PDB. Overwriting with `"UNK"` is not permitted.

### FR-024 — structure_processor.py retired
`structure_processor.py` MUST be retired. All complex structure creation logic MUST be consolidated in `pose_extractor.py` or a dedicated `complex_builder.py`. No stub-only functions that merely create directories may remain in the active module tree.

### FR-025 — PDBQT chain ID detection fix
The fallback PDBQT parser in `pose_extractor.py` MUST use column-based parsing for chain ID rather than the single-character inference heuristic. If column-based parsing fails, the chain ID MUST default to `'A'` with a warning rather than misinterpreting residue numbers.

### FR-026 — Deterministic pose tie-breaking
Best pose selection MUST use a deterministic, filesystem-independent tie-breaking sequence: primary = primary scoring metric (e.g., `cnn_affinity`), secondary = `vina_affinity`, tertiary = `pose_index` (integer, ascending). CSV row order MUST NOT influence the outcome.

### FR-027 — Correlation normalization method
`correlation_analyzer.py` MUST use rank-based normalization (not min-max) when comparing score distributions across engines, so that wide-range engines do not dominate narrow-range ones.

### FR-028 — Empty/malformed log file reporting
`generate_scores_csv.py` MUST track and report log files that produced zero valid score rows (empty output, malformed header, or no scores section found). These MUST be listed in `empty_log_files.txt` in the session output.

### FR-029 — Minimum N per protein for percentile classification
`consensus.py` MUST suppress hit classification (output `classification_basis = "insufficient_n"`) for any target with fewer than 3 ligands after exclusion filtering. Classification on N<3 produces meaningless percentile ranks.

### FR-030 — String escape fix in report generator
`report_generator.py` MUST use `"\n".join(lines)` not `"\\n".join(lines)` for all multi-line string assembly. Dashboard export contract file MUST contain actual newlines.

---

## Key Entities

- **RedockingValidationRecord**: co-crystal pair identity, docked pose path, co-crystal pose path, RMSD value, classification (pass/warn/fail)
- **NormalizedScoreRecord**: raw vina_affinity, raw cnn_affinity, per-protein rank percentile for each, HAC, ligand_efficiency, molecular_weight
- **GeometricConsensusRecord**: ligand-protein pair, per-engine selected pose, pairwise RMSD matrix, geometric_consensus bool, consensus_cutoff_angstrom, single_engine flag
- **EnsembleAffinity**: protein name, snapshot_id, snapshot_affinity, ensemble_min, ensemble_mean, ensemble_std, n_snapshots
- **ReferenceAnchoredClass**: ligand, protein, reference_ligand, reference_affinity, delta_from_reference, hit_class, classification_basis
- **CorrelationRecord**: metric_pair, n_samples, pearson_r, pearson_p_raw, pearson_q_bh, spearman_r, spearman_p_raw, spearman_q_bh, low_power_flag, insufficient_n_flag
- **ComplexStructureRecord**: tag, receptor_file, ligand_file, ligand_residue_name, chain_id, atom_count_receptor, atom_count_ligand, serial_conflict (bool)
- **PoseSelectionRecord**: tag, selected_pose_index, selection_metric, primary_value, secondary_value, tie_broken (bool), tie_break_field

---

## Success Criteria

- **SC-000**: `simplified_pipeline.py` and `multi_engine_pipeline.py` no longer exist as runnable entry points. All analysis runs invoke `unified_pipeline.py`. A single-engine GNINA run and a two-engine run produce outputs from the same code path, differing only in consensus mode annotation.
- **SC-001**: Redocking RMSD is computed for all co-crystal pairs and output to `redocking_validation.csv` in every run with a pairlist.
- **SC-002**: `VALIDATION_SUMMARY.md` is present in every session root and contains a machine-readable confidence signal (`HIGH`, `MEDIUM`, or `LOW`).
- **SC-003**: Cross-protein ranking tables contain both `raw_affinity_kcal_mol` and `per_protein_rank_pct` columns. Sorting by `per_protein_rank_pct` produces a different order than sorting by `raw_affinity_kcal_mol` for the test dataset.
- **SC-004**: GNINA runs use `cnn_score` for pose selection and `cnn_affinity` for ranking. `vina_affinity` is present as a secondary column. All three are in `all_scores.csv`.
- **SC-005**: RMSD clustering no longer emits "Not enough valid poses" warnings. Per_complex scope runs without NaN rows for same-ligand intra-run comparisons.
- **SC-006**: "Best Ligand per Protein" table has exactly N rows where N = number of distinct biological proteins, not number of receptor files.
- **SC-007**: Hit classification outputs include `reference_ligand`, `reference_affinity_kcal_mol`, and `classification_basis` columns.
- **SC-008**: Polypharmacology top table ranks PNBA lower than or equal to its per-target normalized performance. A ligand that scores in the top 35% on ≥3 targets by normalized rank is classified as a polypharmacology candidate.
- **SC-009**: `ligand_efficiency` column is present in all per-pair ranking outputs. No NaN values for complexes where SDF files are available.
- **SC-010**: MD ensemble proteins (3KK6 family → COX1, 3LN1 family → COX2) produce `ensemble_affinity_matrix.csv` with one row per snapshot and `ensemble_summary.csv` with one row per biological protein.
- **SC-011**: A synthetic test with 3 ligands of known affinity order verifies that consensus normalization ranks strongest binder first. The same test fails against the old inverted implementation.
- **SC-012**: Given log files named `"3KK6_site_1_CEL_A_701"` and a pairlist entry `"3KK6_site_1_CEL"`, running `generate_scores_csv` produces no false-positive matches. `unmatched_log_files.txt` is written when any log has no exact pairlist match.
- **SC-013**: Correlation output tables contain `q_value` (BH-corrected) columns. Any run with fewer than 5 paired samples per correlation produces `insufficient_n` annotation rather than a correlation value.
- **SC-014**: A single-engine GNINA run produces no `warn_low_engine_agreement` entries in any output. `agreement_fraction` column contains `NaN` or `"single_engine"` string for all rows.
- **SC-015**: Complex PDB files for co-crystal ligands (e.g., CEL, LIA) contain HETATM records with the correct 3-letter residue code — not `UNK`. ProLIF runs on these complexes without `'ligand'` KeyError.
- **SC-016**: `structure_processor.py` does not exist in the active module tree. `from post_docking_analysis import structure_processor` raises ImportError.
- **SC-017**: Running best-pose selection on the same `all_scores.csv` on macOS and Linux produces byte-identical per-complex pose selections. A test fixture with intentional affinity ties verifies deterministic tie-breaking by pose_index.
- **SC-018**: The dashboard export contract file in `report_generator.py` output contains actual newline characters, not literal `\n` strings.

---

## Risk Areas

- **Redocking SDF extraction**: The co-crystal SDF pose must be extracted from the multi-pose SDF file; indexing may differ by docking run.
- **CNN score availability**: Some GNINA versions or output formats may not export `cnn_affinity`/`cnn_score`; graceful fallback to `vina_affinity` is required.
- **HAC from SDF**: Hydrogens may or may not be explicit in SDF; HAC computation must count only heavy atoms regardless of hydrogen representation.
- **Ensemble detection heuristic**: Detecting "multiple receptor files → same protein" relies on pairlist protein_name column being consistent; if pairlist is absent, ensemble detection must fall back to filename-based grouping with a warning.
- **Reference-anchored thresholds too tight**: For targets where the reference ligand itself docked poorly (high RMSD, bad geometry), anchoring to it will produce misleading Strong classifications. Reference ligands with `redocking_classification = fail` MUST NOT be used as anchors; percentile fallback applies instead.
- **Normalization direction verification**: The inverted normalization fix requires a synthetic unit test. Without it, a regression to inverted logic is silent — the pipeline runs without error but produces wrong rankings.
- **Pairlist matching strictness**: Stricter matching may break projects where log filenames were intentionally longer than pairlist patterns. A migration warning should be emitted on first run after the fix so users can review unmatched files before concluding something is broken.
- **BH correction reduces reported significance**: Users who previously relied on raw p<0.05 correlations from `correlation_analyzer.py` will see fewer "significant" results after BH correction. This is scientifically correct but may feel like a regression. The raw p-value column is retained for reference.
- **Residue name preservation breaks downstream tools relying on UNK**: If any existing downstream script or workflow explicitly filters HETATM records by residue name `"UNK"`, preserving the real residue name will break it. This is a necessary breaking change.
- **structure_processor.py removal**: Any import of `structure_processor` in other modules must be found and removed or redirected as part of FR-024.

---

## Implementation Phases

### Phase 1 — Pipeline Unification (structural prerequisite for everything else)
- **FR-000**: Create `unified_pipeline.py` with the 18-step contract above
- Retire `simplified_pipeline.py` and `multi_engine_pipeline.py`
- **FR-024**: Retire `structure_processor.py` stub; consolidate complex creation in `pose_extractor.py`
- Re-route all CLI entry points to `unified_pipeline`
- At this stage, `unified_pipeline` may internally delegate to existing step implementations — the goal is a single entrypoint, not rewriting every step
- All existing tests must pass through the new entrypoint

### Phase 2 — Critical Bug Fixes (correctness before science)
- **FR-018**: Fix inverted normalization direction in `consensus.py` — add synthetic test
- **FR-019**: Fix bidirectional pairlist matching in `generate_scores_csv.py`; add `unmatched_log_files.txt`
- **FR-021, FR-022**: Fix single-engine QC downgrade and atomic downgrade chain in `consensus.py`
- **FR-023**: Preserve ligand residue name in complex PDB files (`pose_extractor.py`)
- **FR-025**: Fix PDBQT chain ID detection bug in fallback parser
- **FR-026**: Make best-pose tie-breaking deterministic (`pose_extractor.py`)
- **FR-028**: Add empty/malformed log file reporting to `generate_scores_csv.py`
- **FR-030**: Fix string escape bug in `report_generator.py`

### Phase 3 — Validation Gate (unblock scientific trust)
- **FR-001, FR-002, FR-016**: Redocking RMSD + VALIDATION_SUMMARY.md (Step 6)
- **FR-006**: Remove RMSD `per_protein` and `global` scopes; fix clustering warnings

### Phase 4 — Scoring Correctness (fix the ranking science)
- **FR-003**: Per-protein normalization in Step 7 (always, not just multi-engine)
- **FR-004, FR-005**: CNN score utilization for GNINA in Steps 3 and 8
- **FR-013**: Ligand efficiency in Step 11
- **FR-029**: Suppress hit classification for N<3 ligands per target

### Phase 5 — Statistical Rigour
- **FR-020, FR-027**: BH correction + rank normalization in `correlation_analyzer.py`

### Phase 6 — Classification Calibration (fix what rankings mean)
- **FR-011, FR-012**: Reference-anchored hit classification in Step 12
- **FR-009, FR-010**: Per-target normalized polypharmacology in Step 13
- **FR-014, FR-015**: MD ensemble aggregation in Step 9

### Phase 7 — Consensus Overhaul (align with DockBox)
- **FR-007, FR-008**: DockBox geometric consensus in Step 10
- **FR-017**: Remove ADMET stub

---

## References

- McNutt, A.T. et al. (2021). GNINA 1.0: molecular docking with deep learning. *Journal of Cheminformatics*, 13, 43.
- Spitzer, R. & Jain, A.N. (2012). Surflex-Dock: Docking benchmarks and real-world application. *Journal of Computer-Aided Molecular Design*, 26, 687–699.
- Pagadala, N.S. et al. (2017). Software for molecular docking: a review. *Biophysical Reviews*, 9, 91–102.
- Trott, O. & Olson, A.J. (2010). AutoDock Vina. *Journal of Computational Chemistry*, 31, 455–461.
- jp43/DockBox (https://github.com/jp43/DockBox): geometric consensus scoring reference implementation.
