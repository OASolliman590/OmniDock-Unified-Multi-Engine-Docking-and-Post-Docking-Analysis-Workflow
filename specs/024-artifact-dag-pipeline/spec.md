# Feature Specification: Artifact DAG Pipeline

**Feature Branch**: `024-artifact-dag-pipeline`
**Created**: 2026-04-01
**Status**: Draft
**Input**: Architectural review identifying that the current 16-step linear pipeline conflates independent branches (interactions, polypharmacology, comparative analysis) into a forced sequential order, preventing parallelism, incremental re-runs, and clean failure isolation.

---

## Summary

This spec replaces the 16-step linear post-docking pipeline with a **declarative artifact DAG** — a directed acyclic graph where nodes are data artifacts (files) and edges are computations (transforms). The executor derives execution order automatically from declared dependencies, runs independent branches in parallel, skips artifacts whose inputs have not changed, and isolates failures to individual nodes without blocking unrelated work.

The `consensus` node is the **single divergence point** between single-engine and multi-engine workflows. Every artifact upstream and downstream uses the same schema and the same execution path regardless of engine mode. The engine mode is a strategy injected at the `consensus` node only.

All existing specs (017, 019, 020, 021, 022, 023) remain valid — this spec defines the **execution contract** that ties them together.

---

## Problem with the Current Model

The current `unified_pipeline.py` runs 16 steps in a fixed sequence:

```
parse → score → match → complex → affinity → polypharmacology →
rmsd → extract → report → visualize → prolif → pandamap →
poseview → pymol → top_pose → consolidate
```

**Structural failures of this model:**

1. **False sequencing**: ProLIF (step 11) runs before polypharmacology (step 6 is already wrong — it runs before structures exist). Steps have no declared dependency contract, only position.

2. **No parallelism**: ProLIF on 169 complexes and polypharmacology scoring are independent but run serially.

3. **Cascade failures**: One ProLIF crash for one ligand can abort the entire remaining pipeline.

4. **No incrementalism**: Re-running with a different normalization method re-parses all raw scores, re-builds all complexes, re-runs all ProLIF.

5. **Undeclared inputs**: No step declares what files it needs. A step fails with a KeyError at runtime rather than "input artifact missing" at startup.

6. **Single/multi-engine divergence is spread across 8+ steps** rather than isolated to one node.

---

## Goals

1. Replace the linear step sequence with an artifact DAG with explicit input/output contracts per node.
2. Run independent nodes in parallel using a thread pool.
3. Skip nodes whose output artifacts are up-to-date (hash-based cache).
4. Isolate failures per node: one node failing marks it `failed` but does not block sibling nodes.
5. Support scope-based execution: requesting any artifact triggers only its dependency subgraph.
6. Confine single-engine vs. multi-engine divergence to the `consensus` node only.
7. Provide a machine-readable `dag_execution_report.json` after every run.

## Non-Goals

- Adopting an external workflow engine (Snakemake, Prefect, Airflow). The DAG executor is a lightweight internal class (~200 lines).
- Distributed/cloud execution (handled in spec 014 HPC deployment).
- Changing what is computed at each node (covered in upstream specs 017–023).

---

## The Artifact DAG

### Full Graph

```
[project_manifest.json]  [pairlist.csv]  [raw_docking_outputs/]
          │                    │                    │
          └──────────┬─────────┘                    │
                     ▼                              │
           engine_detection_report ◄────────────────┘  (spec 022)
                     │
                     ▼
           engine_scope_config                         (spec 023)
                     │
                     ▼
              raw_scores.csv                           (N-001)
                     │
                     ▼
         normalized_scores.csv                         (N-002)
                     │
                     ▼
    ┌────── consensus_ranked.csv ───────┐              (N-003)
    │      [STRATEGY NODE — see below]  │
    │                                   │
    ▼                                   ▼
validation_gate.json               engine_agreement.csv  (multi only)
    │
    ▼
classified_hits.csv                                    (N-004)
    │
    ├──────────────────┬───────────────────────────────┐
    ▼                  ▼                               ▼
best_poses/       polypharmacology/           comparative/           (parallel)
    │             biology_correlation/        structure_quality/
    ▼
complexes/                                             (N-005)
    │
    ├────────────┬────────────┬───────────────┐
    ▼            ▼            ▼               ▼
prolif/      pandamap/    poseview/        pymol/      (parallel, N-006)
    │
    └──────────────────────────────────────────────────┐
                                                       ▼
                                                  reports/           (N-007)
                                         (waits for all branches)
```

### The Consensus Node Strategy

The `consensus` node is parameterized by `engine_mode` (resolved from `engine_scope_config`):

```
engine_mode = multi  →  DockBox geometric consensus (spec 017)
                         cross-engine rank aggregation
                         engine_agreement.csv produced

engine_mode = single/gnina   →  cnn_affinity primary rank
                                 cnn_confidence annotation
                                 engine_agreement.csv = null

engine_mode = single/vina    →  vina_affinity primary rank
                                 rmsd_lb pose diversity
                                 engine_agreement.csv = null

engine_mode = single/smina   →  vina_affinity primary rank
                                 scoring_function provenance
                                 engine_agreement.csv = null
```

**Output contract of `consensus_ranked.csv` is IDENTICAL for all strategies:**

| Column | Type | Notes |
|--------|------|-------|
| `ligand` | str | |
| `protein` | str | |
| `tag` | str | |
| `engine_mode` | str | `multi \| single` |
| `engines_in_scope` | str | comma-separated |
| `rank` | int | 1 = best |
| `affinity_kcal_mol` | float | primary engine's affinity |
| `primary_score_name` | str | `cnn_affinity \| vina_affinity` |
| `normalized_score` | float | rank percentile 0–1 (1 = best binder) |
| `consensus_score` | float | weighted aggregate or identity for solo |
| `engine_support_count` | int | N engines agreeing; 1 for solo |
| `single_engine_mode` | bool | |
| `scoped_engine_count` | int | engines used this session |
| `total_detected_engines` | int | from detection report |

Everything downstream reads `consensus_ranked.csv` and is blind to engine mode.

---

## Artifact Contracts

### N-001: `raw_scores.csv`

**Inputs**: `{engine}_out/` directories, `pairlist.csv`, `engine_scope_config`
**Computation**: Parse scores from each scoped engine's output (log or scores CSV). Merge into one frame. Attach pairlist metadata.
**Output columns**: `engine, tag, protein, ligand, site_id, pose, affinity_kcal_mol, score_name_primary, score_primary, cnn_score, cnn_affinity, intramol_energy, rmsd_lb, rmsd_ub, pose_file, log_file`
**Failure modes**: `missing_engine_output`, `parse_failed`, `pairlist_mismatch`
**Cache key**: hash of all input file mtimes

### N-002: `normalized_scores.csv`

**Inputs**: `raw_scores.csv`
**Computation**: Per-engine rank percentile (default), z-score, min-max — all direction-correct (most negative affinity → highest score). Spec 020 FR-018.
**Output columns**: all N-001 columns + `rank_pct, z_score, minmax_norm`
**Failure modes**: `insufficient_scores_per_engine` (N < 2)
**Cache key**: hash of `raw_scores.csv`

### N-003: `consensus_ranked.csv` + `engine_agreement.csv`

**Inputs**: `normalized_scores.csv`, `engine_scope_config`, `validation_gate.json`
**Computation**: Strategy node — see consensus strategies above.
**Output**: `consensus_ranked.csv` (schema above), `engine_agreement.csv` (multi only)
**Failure modes**: `no_valid_scores`, `consensus_strategy_failed`
**Cache key**: hash of `normalized_scores.csv` + `engine_scope_config`

### N-004: `classified_hits.csv`

**Inputs**: `consensus_ranked.csv`, `validation_gate.json`, `pairlist.csv`
**Computation**: Hit classification (reference-anchor or target percentile, spec 020 FR-011/012). Ligand efficiency. CNN confidence (GNINA solo). QC gates. ADMET annotation if available.
**Output columns**: all N-003 columns + `hit_class, ligand_efficiency, cnn_confidence, qc_status, admet_status, classification_basis, classification_policy`
**Failure modes**: `no_consensus_scores`, `classification_policy_mismatch`
**Cache key**: hash of `consensus_ranked.csv` + `validation_gate.json`

### N-005: `complexes/`

**Inputs**: `classified_hits.csv`, pose files from engine-specific output directories, top-pose atlas (spec 019)
**Computation**: Select best pose per (ligand, protein) per top-pose policy. Extract pose. Merge receptor PDBQT + ligand pose → complex PDB (ATOM + HETATM). Preserve residue identity (spec 020 FR-026).
**Output**: one `{protein}__{ligand}__complex.pdb` per selected pair
**Failure modes**: `pose_file_missing`, `pdb_merge_failed`, `residue_identity_lost`
**Cache key**: hash of `classified_hits.csv` + pose file mtimes

### N-006: Interaction nodes (parallel)

All four nodes share the same input (`complexes/`, `classified_hits.csv`) and run fully in parallel:

| Node | Output | Tool | Failure isolation |
|------|--------|------|------------------|
| `prolif/` | frequency tables, barcode, network | ProLIF / RDKit | per-complex; others continue |
| `pandamap/` | 2D interaction maps | PandaMap | per-complex |
| `poseview/` | 2D diagrams | PoseView CLI | per-complex |
| `pymol/` | 3D sessions + scripts | PyMOL | per-complex; skip if binary absent |

Each interaction node runs a sub-DAG: one task per complex, parallelized via `ThreadPoolExecutor(max_workers=rmsd_workers)`.

### N-007: `reports/`

**Inputs**: `classified_hits.csv`, `engine_agreement.csv`, `polypharmacology/`, `comparative/`, `biology_correlation/`, `interactions/`, `complexes/`
**Computation**: Generate `START_HERE.md` (spec 021), HTML/CSV hit reports, consolidated run summary, dashboard export contract, `dag_execution_report.json`.
**Failure modes**: `missing_upstream_artifact` (non-fatal — reports generated with available inputs, missing sections noted)
**Cache key**: hash of all input artifact mtimes

---

## The ArtifactGraph Executor

A lightweight internal executor in `post_docking_analysis/artifact_graph.py`:

```python
class ArtifactNode:
    name: str
    inputs: List[str]           # artifact names this node depends on
    outputs: List[str]          # artifact paths this node produces
    compute: Callable           # function to run
    optional: bool = False      # if True, failure doesn't block dependents

class ArtifactGraph:
    def register(self, node: ArtifactNode): ...
    def request(self, artifact: str, force: bool = False): ...
    def _is_cached(self, node: ArtifactNode) -> bool: ...   # hash check
    def _topo_sort(self, target: str) -> List[List[ArtifactNode]]: ...
    def _execute(self, plan: List[List[ArtifactNode]]): ...  # ThreadPoolExecutor per tier
```

**Execution model:**
1. `request(artifact)` → trace all dependencies recursively
2. Filter out nodes whose outputs are up-to-date (hash unchanged) unless `force=True`
3. Topological sort → group into tiers (nodes in the same tier have no inter-dependencies)
4. Execute each tier with `ThreadPoolExecutor`; advance to next tier when all complete (or fail)
5. Failed nodes are marked `failed`; their dependents are skipped with `blocked_by_failure`
6. Optional nodes (`prolif`, `pymol`) that fail do not mark dependents as blocked

**Cache invalidation:**
- Each node computes `cache_key = sha256(sorted(input_file_mtimes + input_hashes))`
- Cache key stored in `4-Working/dag_cache.json`
- If `cache_key` matches stored key and all outputs exist → skip

---

## Scope System

The `--scope` flag maps to an artifact request:

```python
SCOPE_MAP = {
    "full":            "reports",
    "comparison_only": "classified_hits.csv",
    "structures_only": "complexes",
    "interactions":    "interactions",
    "report_only":     "reports",          # deps cached; only N-007 runs
    "top_pose_only":   "best_poses",
}
```

Requesting `classified_hits.csv` triggers: ingestion → normalization → consensus → classification. No structures, no interactions, no reports. Exactly the minimum needed.

---

## Execution Report

`dag_execution_report.json` written after every run:

```json
{
  "run_id": "interactive_20260401_050015",
  "engine_mode": "multi",
  "engines_in_scope": ["gnina", "vina", "smina"],
  "scope_requested": "full",
  "artifact_requested": "reports",
  "nodes": {
    "raw_scores": {"status": "completed", "duration_s": 1.2, "cache_hit": false},
    "normalized_scores": {"status": "completed", "duration_s": 0.3, "cache_hit": true},
    "consensus_ranked": {"status": "completed", "duration_s": 2.1, "cache_hit": false},
    "classified_hits": {"status": "completed", "duration_s": 0.8, "cache_hit": false},
    "complexes": {"status": "completed", "duration_s": 14.3, "cache_hit": false},
    "prolif": {"status": "completed", "duration_s": 187.2, "cache_hit": false},
    "pandamap": {"status": "completed", "duration_s": 42.1, "cache_hit": false},
    "poseview": {"status": "failed", "error": "binary_not_found", "cache_hit": false},
    "pymol": {"status": "skipped_optional", "cache_hit": false},
    "polypharmacology": {"status": "completed", "duration_s": 3.4, "cache_hit": false},
    "comparative": {"status": "completed", "duration_s": 1.1, "cache_hit": false},
    "reports": {"status": "completed_with_warnings", "duration_s": 2.2, "cache_hit": false}
  },
  "total_duration_s": 201.4,
  "parallel_time_saved_s": 143.6
}
```

---

## Migration from 16-Step Linear Pipeline

| Old step | New node | Notes |
|---|---|---|
| 1: find files | ingestion (pre-N-001) | |
| 2: generate scores | N-001: raw_scores | |
| 3: match poses | N-001: raw_scores | merged into ingestion |
| 4: create complexes | N-005: complexes | now after classification |
| 5: binding affinity | N-002+N-003: normalize+consensus | |
| 6: polypharmacology | parallel branch from N-004 | no longer after RMSD |
| 7: RMSD | sub-task inside N-005 | per-complex only |
| 8: extract poses | N-005: complexes | |
| 9: generate reports | N-007: reports | now last |
| 10: visualizations | N-006: interactions | parallel |
| 11: ProLIF | N-006: prolif | parallel |
| 12: PandaMap | N-006: pandamap | parallel |
| 13: PoseView | N-006: poseview | optional |
| 14: PyMOL | N-006: pymol | optional |
| 15: top pose | N-004 + N-005 | top-pose atlas between classify and structure |
| 16: consolidate | N-007: reports | |

**Key ordering fix**: complex creation (old step 4) now happens after classification (N-004), not before. This means only hits above QC threshold get complex PDBs built — saving significant I/O for large ligand sets.

---

## Integration with Prior Specs

| Spec | Integration point |
|------|------------------|
| 017: DockBox Consensus | Multi-engine strategy inside N-003 |
| 019: Top Pose Atlas | Between N-004 and N-005; top_pose_atlas artifact |
| 020: Scientific Foundation | N-001 (score parsing), N-002 (normalization), N-003 (consensus bugs), N-005 (complex PDB) |
| 021: Output Topology | All artifact write paths; `5-Analysis/` root enforced |
| 022: Engine Detection | Pre-graph; populates `engine_detection_report` before N-001 |
| 023: Engine Scope Filter | Pre-N-003; populates `engine_scope_config` |

---

## Acceptance Criteria

1. A 3-engine full-scope run completes with ProLIF and polypharmacology running in parallel; `dag_execution_report.json` shows non-zero `parallel_time_saved_s`.
2. Re-running with unchanged inputs after a full run: only `reports` node re-executes (all others `cache_hit: true`).
3. Re-running after changing normalization method: N-002 cache miss triggers N-003, N-004, N-005, N-006, N-007 re-execution; N-001 remains cached.
4. ProLIF parse failure for one complex marks `prolif` node as `completed_with_warnings`; `pandamap`, `polypharmacology`, `comparative`, and `reports` nodes all complete normally.
5. `--scope comparison_only` runs nodes N-001 through N-004 only; no pose files or complexes written.
6. GNINA solo run: `consensus_ranked.csv` `primary_score_name = cnn_affinity`; same schema as multi-engine run.
7. `classified_hits.csv` schema is identical for single-engine and multi-engine runs.
8. Optional node failure (`poseview`, `pymol`) does not block `reports` node.
9. `dag_execution_report.json` written on every run; contains status for every node in the requested subgraph.
10. Requesting `--scope report_only` on a project where all upstream artifacts are cached generates reports in < 5 seconds.

---

## Risks

- Hash-based caching requires stable artifact serialization (CSV column order, float precision). Must use deterministic writers everywhere.
- ThreadPoolExecutor for I/O-bound interaction nodes (ProLIF, PandaMap) is appropriate. CPU-bound consensus computation should not be parallelized with GIL-bound threads — use `ProcessPoolExecutor` or serial execution for N-003.
- The migration reorders step 4 (complex creation) to after classification. Existing users who run `--scope structures_only` expecting complexes for all ligands will now only get complexes for classified hits. This is intentional but must be documented.

---

## Implementation Phases

### Phase 1: ArtifactGraph Executor
Implement `post_docking_analysis/artifact_graph.py` with `ArtifactNode`, `ArtifactGraph`, topological sort, tier-based `ThreadPoolExecutor`, hash cache (FR: all executor contracts above).

### Phase 2: Node Registration
Register all nodes (N-001 through N-007 + sub-nodes) with their input/output contracts. Wire existing pipeline functions as `compute` callables — no logic changes yet.

### Phase 3: Consensus Strategy Injection
Implement the strategy pattern at N-003: inject engine mode from `engine_scope_config`; confirm `consensus_ranked.csv` schema is identical for all strategies.

### Phase 4: Parallel Interaction Nodes
Convert N-006 sub-nodes (prolif, pandamap, poseview, pymol) to run as parallel tasks within the `ThreadPoolExecutor`; add per-complex failure isolation.

### Phase 5: Scope System
Wire `SCOPE_MAP` to `graph.request(artifact)`; validate all scope values produce correct subgraph.

### Phase 6: Execution Report + Cache
Emit `dag_execution_report.json`; implement hash-based cache with `4-Working/dag_cache.json`.

### Phase 7: Migration + Tests
Remove old step-sequence loop from `unified_pipeline.py`; replace with `graph.request(scope_artifact)`. Smoke tests for all acceptance criteria.
