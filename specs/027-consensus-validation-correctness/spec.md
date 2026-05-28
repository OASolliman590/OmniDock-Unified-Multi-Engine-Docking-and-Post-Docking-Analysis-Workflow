# Spec 027: Consensus and Validation Correctness

## Problem

Three bugs deferred from Spec 026 silently corrupt scientific results for multi-engine GNINA/Vina runs:

1. **Geometric consensus `pairwise_mapping_incomplete`** — `_pairwise_rmsd()` fails when GNINA SDF and Vina PDBQT atoms don't map one-to-one (element sequence mismatch or atom count difference). Sets `geometric_agreement=False` for 185 pairs despite actual RMSD of 0.02 Å. Silently downgrades all 185 complexes to Weak/Uncategorized in hit classification.

2. **GNINA SDF multi-pose parse in redocking validation** — `_compute_pose_rmsd()` in `redocking_validation.py` calls `_parse_pdb_like_heavy_atoms()`, which delegates to `_parse_sdf_heavy_atoms()`. That function reads the full file's atom block from a fixed line-offset approach. For multi-pose GNINA SDF, line-offset parsing can misread atom counts (88 = 2×44 reported vs 44 expected), producing `atom_count_mismatch:88!=44` for all GNINA complexes. This collapses the validation gate to `needs_review` universally.

3. **`favorite_engine_continue` — no upfront engine membership check** — `run()` checks `if not favorite_engine` (empty string) but never verifies the named engine is in `self.engines_in_scope`. An invalid engine name passes silently and fails late with an unhelpful downstream error.

---

## Goal

Fix all three bugs with targeted, minimal changes. No new public API beyond small parameter additions.

---

## Affected Files

| File | Issue |
|---|---|
| `post_docking_analysis/geometric_consensus.py` | 1 |
| `post_docking_analysis/redocking_validation.py` | 2 |
| `post_docking_analysis/_legacy_multi_engine_pipeline_impl.py` | 3 |

---

## Fix Definitions

### Fix 1 — Geometric consensus: tolerate cross-format atom mapping

**File:** `geometric_consensus.py`

**Root cause:** `_pairwise_rmsd()` (line ~243) immediately rejects pairs with different `coords_a.shape[0] != coords_b.shape[0]` as `atom_count_mismatch`. For GNINA SDF vs Vina PDBQT this happens because:
- GNINA SDF: element parsed from fixed columns `[31:34]` of each atom line
- Vina PDBQT: element parsed from `[76:78]` (standard PDB element column) or derived from atom name
- Different parsers may include/exclude certain atoms (e.g., atoms at block boundaries, aromatic disambiguation)

Result: even when both engines bind to the same location (RMSD 0.02 Å), `pairwise_mapping_incomplete` is recorded and `geometric_agreement = False`.

**Fix strategy — two-tier fallback in `_pairwise_rmsd()`:**

**Tier 1 (existing):** Require equal atom count AND equal element sequence → compute RMSD exactly. No change.

**Tier 2 (new fallback when Tier 1 fails):** When `atom_count_mismatch` or `element_sequence_mismatch` occurs between two engines with different pose formats (`.sdf` vs `.pdbqt`), attempt a *soft alignment*:
- Take the minimum atom count of the two poses as `n_common`
- Sort both coordinate arrays by distance from their own centroid
- Compute RMSD on the first `n_common` atoms (centroid-sorted order)
- Record result with `reason = "soft_aligned_count_mismatch"` and a new flag `soft_alignment_used = True`
- Only mark `geometric_agreement = True` if soft-aligned RMSD ≤ cutoff **and** `|count_a - count_b| <= 2` (i.e., discrepancy is small — likely hydrogen handling difference, not a different molecule)

**Tier 3 (existing fallback when Tier 2 still fails):** Record `pairwise_mapping_incomplete` as before.

**Changes to `_pairwise_rmsd()`:**
- Accept optional `pose_formats: Dict[str, str]` parameter (engine → suffix, e.g. `{"gnina": ".sdf", "vina": ".pdbqt"}`)
- When `atom_count_mismatch` AND formats differ AND `|count_a - count_b| <= 2`: attempt Tier 2 soft alignment
- Add `soft_alignment_used` field to each `pair_row`

**Changes to `compute_geometric_consensus()`:**
- Collect each engine's pose file suffix from the `pose_file` column
- Pass `pose_formats` to `_pairwise_rmsd()`

**Changes to `geometric_agreement` final condition (line ~332):**
```python
# Old: valid_pairs == expected_pairs (strict)
# New: valid_pairs == expected_pairs OR all failing pairs used soft alignment with RMSD <= cutoff
```
Specifically: a tag's `geometric_agreement = True` when all pairwise RMSDs (including soft-aligned) are ≤ cutoff and no pair had a hard failure (missing file, wrong molecule).

---

### Fix 2 — Redocking validation: use `$$$$`-split SDF parser for GNINA poses

**File:** `redocking_validation.py`

**Root cause:** `_parse_sdf_heavy_atoms()` (line ~78) reads heavy atoms from a flat `lines` list using `natoms = int(lines[3][0:3])`. For GNINA multi-pose SDF files this can fail in two ways:
1. If the SDF header has a blank first line or non-standard spacing, `natoms` is mis-parsed → reads too many or too few atoms
2. If `natoms` is parsed correctly but the atom table is followed immediately by bond lines that accidentally parse as atom coordinates, extra atoms are counted

`geometric_consensus.py` already has a robust `_extract_sdf_pose(path, pose_index)` that splits on `$$$$` first, then reads each block independently. This is the correct approach for multi-pose SDF.

**Fix:**
- In `_parse_pdb_like_heavy_atoms()` (line ~110), when `file_path.suffix.lower() in (".sdf", ".mol")`, call `_extract_sdf_pose(file_path, pose_index=1)` from `geometric_consensus.py` instead of `_parse_sdf_heavy_atoms(lines)`.
- `_extract_sdf_pose` returns `(coords, elements, error)` — same signature as the existing parse functions.
- Remove `_parse_sdf_heavy_atoms` from `redocking_validation.py` — it is now replaced by the shared implementation.
- Add `pose_index: int = 1` parameter to `_parse_pdb_like_heavy_atoms()` so callers can specify which pose to extract (defaults to 1 = best pose).
- Update `_compute_pose_rmsd()` to pass `pose_index=1` explicitly.

**Import:** `from post_docking_analysis.geometric_consensus import _extract_sdf_pose` at the top of `redocking_validation.py`. Since this creates a dependency between two internal modules, add a note in `geometric_consensus.py` that `_extract_sdf_pose` is part of the semi-public parser API.

---

### Fix 3 — `favorite_engine_continue`: validate engine membership before running

**File:** `_legacy_multi_engine_pipeline_impl.py`

**Root cause:** `run()` at line ~1842:
```python
favorite_engine = self.favorite_engine or self.manifest.get("favorite_engine") or self.engine
if not favorite_engine:
    raise ValueError("favorite_engine_continue requires ...")
```
The check only guards against empty string. If `favorite_engine = "gnina"` but the project only has Vina outputs, the code continues, calls `self._write_single_engine_reports(scores, "gnina", ...)`, finds no gnina rows in `scores`, and produces an empty report — no error, no warning.

**Fix:**
- After the empty-string check, add a validation guard:
```python
if self.engines_in_scope and favorite_engine not in self.engines_in_scope:
    raise ValueError(
        f"favorite_engine '{favorite_engine}' is not in engines_in_scope "
        f"{self.engines_in_scope}. Use --engine to specify a valid engine."
    )
```
- If `self.engines_in_scope` is empty (not yet populated), fall back to checking `self.valid_engines` from the detection report if available.
- The check must be at the start of the `favorite_engine_continue` branch, before any file writes or report generation.

---

## Inputs

No new data files. All fixes operate on existing pose files and code paths.

## Outputs

No new artifact files. Fix 1 adds `soft_alignment_used` and `soft_alignment_rmsd` fields to `geometric_pairwise_rmsd_json` entries (additive, backward-compatible).

---

## Constraints

- Fix 1: `geometric_agreement` must remain `False` when `|count_a - count_b| > 2` (different molecules or gross parse failure — do not soft-align these)
- Fix 1: `soft_alignment_used = True` must be recorded in the pairwise JSON so the result is auditable
- Fix 2: `_extract_sdf_pose` must not be renamed or made private-only; it is now shared by two modules
- Fix 3: the guard must only fire when `engines_in_scope` is non-empty; an empty scope list (all engines allowed) must not block execution
- All fixes: backward-compatible defaults; existing callers with no new arguments retain identical behavior

---

## Non-Goals

- Fixing SMINA/Vina 20×20×20 box size for 3LN1 MS-series — those are correct runs for the CEL pocket; the alternative binding seen by GNINA is a separate phenomenon
- Geometric consensus for 3+ engine runs using full Hungarian matching — Tier 2 soft alignment (centroid sort) is sufficient for 2-engine mismatches
- Changing the redocking validation RMSD thresholds or gate logic
