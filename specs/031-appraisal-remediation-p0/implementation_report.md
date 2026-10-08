# Spec 031 — Gate 1 and Gate 2 Implementation Report

**Date:** 2026-08-19  
**Status:** Gate 1 engineering and Gate 2 scientific-method scopes implemented, locally verified, and accepted by the Scientific Lead on 2026-08-21.  
**Baseline:** `9d7400c` with preserved pre-existing working-tree changes.  
**External execution:** `skipped_disabled` — no docking engine execution, HPC submission, upload, cloud action, or remote mutation was performed. The all-engine workflow check used dry-run mode only.

## Outcome

The evidence-backed engineering corrections from the revised plan are implemented. The original assertion that Vina output was dispatched to the GNINA parser was disproved in both the assessed and current code. The retained zero-row observation was therefore treated as an ingestion/layout defect: artifacts are now resolved read-only and per type, mixed layouts are diagnosed, and parse-to-zero outcomes fail explicitly.

No scientific policy was guessed. The Phase B implementations follow the four dated Scientific Lead decision records and remain opt-in or additive where required. No external docking, HPC/cloud execution, upload, or remote mutation occurred.

## Before and After

| Area | Before | After |
|---|---|---|
| Vina failure classification | Plan claimed Vina was routed through the GNINA parser. | Evidence shows Vina already used `parse_vina_pdbqt`; the superseded claim is retained as such and remediation targets artifact resolution. |
| Ingestion/layout | Readers depended on fixed directories and could silently produce zero rows. | Scores, poses, logs, and metadata resolve independently across canonical, flat, numbered, legacy, and mixed layouts without migration writes; ambiguity and selected paths are reported. |
| Empty ingestion | Missing or mismatched artifacts could look like successful empty analysis. | Zero parsed rows raise an actionable diagnostic with engine, parser, expected artifacts, resolved paths, and ambiguity information. |
| Cache/resume | GNINA raw aggregation, normalized caches, DAG dependencies, and detection cache could disagree or remain stale. | Raw artifact and pair-intent signatures participate in invalidation; GNINA `all_scores.csv` refreshes; `force` propagates through detection, aggregation, and DAG execution. |
| Pose selection | Redocking could ignore the selected row/model or fall back to model 1. | Shared strict one-based PDB/PDBQT model extraction honors the selected pose, accepts a final EOF-closed model, records its index, and rejects malformed/out-of-range selection. |
| Co-crystal metadata | False-like strings/NaN could become truthy and pairlist metadata could be lost. | Central boolean normalization and atomic pairlist propagation preserve explicit co-crystal/source fields without inventing lineage. |
| Plot compatibility | The affected DataFrame operation was incompatible with newer pandas and lacked a real plot regression. | Modern `DataFrame.map` is used with compatibility fallback; the actual hit-class matrix path is executed and verified to write a non-empty figure. |
| Manifest paths | Relative paths could contain platform-native backslashes. | Serialized manifest/report paths use POSIX `/` while runtime filesystem resolution remains native. |
| Exact deployment | Deployment could not consume a one-off exact pairlist without affecting canonical state. | `--pairlist-file` is mutually exclusive with rerun manifests, deploys only those rows, preserves canonical pairlist/curation state, and records source path, kind, and SHA-256 in project and per-engine Slurm/Condor manifests. |
| Packaging/runtime | The installed console target was not the packaged CLI callable; direct imports were absent from dependency manifests. | Entry point targets `workflow.cli:main`; demonstrated direct dependencies include PyYAML, SciPy, and scikit-learn. |
| Scientific wording | Pocket outputs could read as calibrated druggability measurements. | Outputs disclose `uncalibrated_pocket_heuristics_v1`, method fields, non-calibrated status, and limitations; legacy numeric behavior is preserved. |
| Windows/runtime robustness | Symlink-only smoke assertions, a system-Python subprocess, checkpoint timestamp ties, and concurrent job JSON reads caused nondeterminism. | Smoke accepts the documented README fallback, subprocesses use the active interpreter, checkpoint ties are stable, and job-record I/O is synchronized/retried. |
| RMSD correspondence | Scientific RMSD could rely on atom order or centroid-sorted soft alignment without proving molecular correspondence. | `spec031-atom-mapping-v1` requires complete heavy-atom, bond-order/charge/stereo-labelled graph isomorphism, enumerates symmetry, reports minimum Kabsch RMSD, and emits no numeric value for `not_comparable`. |
| Consensus semantics | Legacy modes mixed ranking/geometry behavior and GNINA pose/ranking metrics were not governed by the approved v2 contract. | Opt-in `consensus_rank_geometry_qc_v2` uses engine-specific pose selection, within-target complete cases, equal-weight percentiles, no imputation/renormalization, missing-engine disclosure, and comparable-only geometry QC/ties; legacy modes remain available unchanged. |
| Binding-site center | Enhanced and fallback preparation could replace the ligand centroid with a PLIP/contact-residue atom centroid. | `binding_site_center_v1` uses the explicitly selected holo ligand's heavy-atom centroid (or explicit user coordinates), separates contact centroids, records source/model/altloc/checksum provenance, and returns `skipped_missing_configuration` for unconfigured apo cases. |
| Scientific reference status | Tags, site labels, and filename/display-name patterns could create scientific reference/redocking status. | `explicit_reference_metadata_v2` requires an explicit normalized flag plus identity/lineage. Misleading names remain non-reference; legacy inference is opt-in, warned, report-only, and never anchor-eligible without confirmed metadata. |
| Atom-mapping dependency | The declared Open Babel environment did not provide a reliable in-process, edge-labelled automorphism contract for the approved fixtures. | NetworkX `>=3.2,<4` is the single bounded direct backend in all applicable manifests; missing backend yields `skipped_missing_dependency`/`not_comparable` with no heuristic fallback. |

## Verification

| Check | Result |
|---|---|
| Spec 031 focused regression files | `completed` — 31 passed |
| Gate 2 focused suites | `completed` — 30 passed (mapping, consensus v2/top-pose direction, center, reference policy, and selected-pose integration) |
| All Spec 031 regression files | `completed` — 54 passed, 16 warnings |
| Full pytest suite | `completed` — 161 passed, 1 skipped, 16 warnings |
| Complete DockForge smoke suite | `completed` — all checks passed |
| All-engine CLI workflow | `completed` — GNINA, Vina, Smina, and AutoDock4 dry-run manifests/assets verified |
| Standalone `test_pipeline.py` | `completed` |
| Changed-tree `compileall` | `completed` |
| `git diff --check` | `completed` |
| NMRbox retained-project reproduction | `skipped_missing_configuration` |
| External HPC/cloud/docking execution | `skipped_disabled` |

The 16 pytest warnings are non-failing Windows symlink-fallback notices plus one Matplotlib `tight_layout` warning. The complete smoke run also reports unavailable optional tools (Open Babel and `jq`) where tests intentionally exercise fallback/skip behavior. Some pre-existing logging handlers emit CP1252 Unicode diagnostics during the Windows smoke run; the suite still completes successfully, but that console-noise issue is outside the approved remediation semantics.

## Retained Scientific Fixture Assumptions

- Receptor: PDB `1IEP`, chain `A`.
- Ligand: `STI`, residue `A:201`.
- Retained Vina protocol: Vina `1.2.7` Python mode, `sf=vina`, center `[16.757738, 56.582794, 13.778150]` Å, box `20 × 20 × 20` Å, exhaustiveness `8`, seed `12345`, CPU `2`, requested poses `5`.
- Protonation pH, tautomer states, charge model, hydrogen policy, and exact preparation-tool version remain unknown because the retained evidence does not establish them.
- The fixture is an offline regression artifact, not a newly approved biological protocol and not evidence of affinity, selectivity, mechanism, or efficacy.

## Gate 2 Approval Addendum — 2026-08-21

The Scientific Lead approved the recommended baseline, including the bounded dependency policy. Normative records:

- `decisions/atom-mapping.md`
- `decisions/consensus-v2.md`
- `decisions/grid-center.md`
- `decisions/reference-classification.md`

The bounded implementation sequence is `phase_b_plan.md`. T013, T015, T017, and T019 are implemented and locally verified. T032 was accepted and closed on 2026-08-21 under `decisions/gate-2-implementation-acceptance.md`.

## Gate 2 Scientific Assumptions and Limits

- The code does not select or alter receptor accession, ligand identity, chain/residue, pH, protonation, tautomer, charge model, box dimensions, exhaustiveness, seeds, engine version, or scoring mode.
- PDB/PDBQT coordinates without an explicit SDF/MOL bond topology are `not_comparable`; no atom-order or centroid fallback is used in current scientific outputs.
- V2 consensus is a relative within-target docking rank, not an experimental affinity, selectivity, mechanism, or efficacy claim. Raw engine scores are not equated.
- Reference-affinity anchors require explicit source correspondence and a passing graph-mapped redocking result. Names and legacy candidates cannot satisfy that requirement.
