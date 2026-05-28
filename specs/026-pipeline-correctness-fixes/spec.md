# Spec 026: Pipeline Correctness Fixes

## Problem

A cross-code audit revealed 12 correctness bugs and silent-failure paths across the post-docking pipeline. They fall into five categories:

1. **Top-pose selector direction bug** — `consensus_score` sort direction and summary aggregation are correct for rank-pct consensus modes but directionally wrong for agreement-based modes (`dockbox_geometric`, `favorite_guardrails`), silently picking the weakest pose instead of the strongest.
2. **`reference_anchor` classification not wired in DAG** — `_dag_compute_classified_hits_node` always passes `reference_baselines=pd.DataFrame()`, so the `reference_anchor` hit-class policy silently degrades to percentile behavior even when a validated baseline table exists on disk.
3. **Engine detection cache ignores `override_engine`** — `detect_engines()` returns the cached manifest entry before the override logic runs, so `--engine gnina` on a cached multi-engine project has no effect without `--redetect`.
4. **RMSD scope interface dead code** — `_normalize_rmsd_scopes()` only accepts values already in `_DEFAULT_RMSD_SCOPES = ("per_complex",)`, making `per_protein` and `global` scopes permanently unreachable despite being handled in downstream code.
5. **`strict_consensus` silently no-ops in single-engine mode** — the agreement-count filter (`agreement_count >= 2`) is skipped when `single_engine=True`, so the user's policy choice has no effect and no warning is emitted.

---

## Goal

Fix all 12 defects (enumerated below) with minimal, targeted changes. No new features, no refactoring beyond what each fix requires.

---

## Affected Files

| File | Issues |
|---|---|
| `post_docking_analysis/top_pose_selector.py` | 1a, 1b |
| `post_docking_analysis/_legacy_multi_engine_pipeline_impl.py` | 2a, 2b |
| `post_docking_analysis/engine_detector.py` | 3 |
| `post_docking_analysis/_legacy_simplified_pipeline_impl.py` | 4 |
| `post_docking_analysis/consensus.py` | 5 |

---

## Fix Definitions

### Fix 1a — `consensus_score` sort direction in `build_top_pose_atlas()`

**File:** `top_pose_selector.py`

**Root cause:** `consensus_score` is computed differently per consensus mode:
- `strict_consensus` / `weighted_hybrid`: formula is `mean_rank_pct * w1 + best_affinity_rank_pct * w2 + ...` — a rank-percentile blend where **lower = better**.
- `dockbox_geometric`: formula is `geometric_agreement (0/1) + mean_rank_pct * 1e-3` — an agreement score where **higher = better** (geometric_agreement dominates).
- `favorite_guardrails`: formula is `base_favorite - guardrail_penalty` — an explicit score where **higher = better**.

The atlas sort is always `ascending=True`, which is correct for rank-pct modes but picks the worst pose for agreement-based modes.

**Fix:**
- Add a `consensus_mode: str = ""` parameter to `build_top_pose_atlas()`.
- Define a helper `_consensus_score_ascending(mode: str) -> bool` that returns `True` for rank-pct modes (`strict_consensus`, `weighted_hybrid`, default/empty) and `False` for agreement-based modes (`dockbox_geometric`, `favorite_guardrails`).
- Pass the boolean into the three sort calls at lines ~204–216 and ~225 so ascending direction varies per mode.
- Propagate `consensus_mode` from `_dag_compute_top_pose_atlas_node()` using `self.consensus_mode`.

### Fix 1b — `best_consensus_score` summary aggregation

**File:** `top_pose_selector.py`

**Root cause:** Summary aggregation at line ~232:
```python
best_consensus_score=("consensus_score", "min"),
```
`min` is correct for rank-pct modes (lower = better) but wrong for agreement-based modes (should be `max`).

**Fix:**
- Compute `best_consensus_score` conditionally:
  - If `_consensus_score_ascending(consensus_mode)` → use `min`
  - Else → use `max`
- Since pandas `agg` doesn't accept a conditional aggregation function inline, compute this as a separate step after the main `agg()` call by replacing the column.

---

### Fix 2a — `reference_baselines` always empty in DAG classified-hits node

**File:** `_legacy_multi_engine_pipeline_impl.py`

**Root cause:** `_dag_compute_classified_hits_node()` (line ~685–693) hardcodes:
```python
reference_baselines=pd.DataFrame(),
```
The validation gate node writes `reference_baselines.csv` to disk but `_dag_compute_classified_hits_node` never reads it.

**Fix:**
- In `_dag_compute_classified_hits_node()`, load `reference_baselines.csv` from the validation gate output path before calling `classify_hits_target_aware()`:
  ```python
  baselines_path = paths["validation_gate"].parent / "reference_baselines.csv"
  reference_baselines = pd.read_csv(baselines_path) if baselines_path.exists() else pd.DataFrame()
  ```
- Pass the loaded frame (or empty DataFrame if absent) to `reference_baselines=`.

### Fix 2b — `reference_baselines` path must be added to `_dag_artifact_paths()`

**File:** `_legacy_multi_engine_pipeline_impl.py`

**Root cause:** `paths["validation_gate"]` points to `validation_gate_status.json`. `reference_baselines.csv` lives in the same directory but has no canonical key in `_dag_artifact_paths()`.

**Fix:**
- Add `"reference_baselines"` to `_dag_artifact_paths()` pointing to `<validation_gate_dir>/reference_baselines.csv`.
- Use `paths["reference_baselines"]` in Fix 2a instead of a manual path construction.

---

### Fix 3 — Engine detection cache ignores `override_engine`

**File:** `engine_detector.py`

**Root cause:** Lines 199–202:
```python
if not redetect and isinstance(manifest.get("detected_engines"), dict):
    cached = dict(manifest["detected_engines"])
    cached["cache_hit"] = True
    return cached   # ← returns before override logic at lines 227–233
```
`override_engine` is never applied to the cached result.

**Fix:**
- After retrieving the cached dict, apply the override before returning:
  ```python
  if not redetect and isinstance(manifest.get("detected_engines"), dict):
      cached = dict(manifest["detected_engines"])
      cached["cache_hit"] = True
      override = str(override_engine or "").strip().lower()
      if override:
          cached["detection_override"] = True
          cached["override_reason"] = f"--engine {override}"
          cached["routed_to_engine"] = override
          cached["routing_decision"] = "single_engine"
      return cached
  ```

---

### Fix 4 — RMSD scope normalizer blocks `per_protein` and `global`

**File:** `_legacy_simplified_pipeline_impl.py`

**Root cause:** `_normalize_rmsd_scopes()` (lines 70–88) only accepts tokens already in `_DEFAULT_RMSD_SCOPES = ("per_complex",)`. The guard at line 86:
```python
if value in _DEFAULT_RMSD_SCOPES and value not in normalized:
```
filters out `"per_protein"` and `"global"` because they are not in `_DEFAULT_RMSD_SCOPES`. The downstream code at lines 2117–2119 that reads these scopes is therefore dead.

**Fix:**
- Define the full valid set of scope tokens:
  ```python
  _VALID_RMSD_SCOPES = ("per_complex", "per_protein", "global")
  _DEFAULT_RMSD_SCOPES = ("per_complex",)
  ```
- Change the guard in `_normalize_rmsd_scopes()` to check against `_VALID_RMSD_SCOPES`:
  ```python
  if value in _VALID_RMSD_SCOPES and value not in normalized:
      normalized.append(value)
  ```
- The `"all"` / `"*"` expansion should return all valid scopes: `return _VALID_RMSD_SCOPES`.

---

### Fix 5 — `strict_consensus` silently no-ops in single-engine mode

**File:** `consensus.py`

**Root cause:** Line 328:
```python
if mode == "strict_consensus" and not single_engine:
    by_tag = by_tag[by_tag["agreement_count"] >= 2].copy()
```
In single-engine mode, `strict_consensus` has no filtering effect and there is no warning.

**Fix:**
- When `mode == "strict_consensus"` and `single_engine=True`, add a `logger.warning()` (or populate a `warnings` list that is returned in the result metadata) to explicitly state that strict consensus filtering cannot be applied in single-engine mode.
- Document the behavior in the function docstring: *"In single-engine mode, strict_consensus cannot enforce multi-engine agreement and behaves identically to weighted_hybrid. A warning is emitted."*
- Do **not** change the filtering logic — skipping the filter in single-engine mode is scientifically correct since there is only one engine. The fix is transparency, not behavior change.

---

## Inputs

No new data inputs. All fixes operate on existing DAG artifacts and code paths.

## Outputs

No new artifact files. Fix 2a/2b reads an existing file that was already being written.

---

## Constraints

- No new public API surface beyond `consensus_mode` parameter on `build_top_pose_atlas()`.
- All fixes must be backward-compatible: default arguments preserve existing behavior for callers that do not pass the new parameter.
- No changes to the DAG topology (no new nodes, no new edges) beyond adding the `reference_baselines` path key.
- Tests must cover the directional fix: a `best_consensus` policy with `dockbox_geometric` mode must select the row with the **highest** `consensus_score`, not lowest.

---

## Non-Goals

- Fixing the geometric consensus `pairwise_mapping_incomplete` bug (GNINA SDF ↔ Vina PDBQT atom index mismatch) — tracked separately.
- Fixing GNINA multi-pose SDF atom count mismatch in redocking validation — tracked separately.
- Fixing `favorite_engine_continue` upstream validation — tracked separately.
- Any UI or CLI changes.
