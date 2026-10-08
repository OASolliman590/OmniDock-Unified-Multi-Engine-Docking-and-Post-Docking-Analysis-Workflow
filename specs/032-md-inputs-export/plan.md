# Implementation Plan: Strict MD-Consumer Input Export

## 1. Gate State

**Current gate**: Gate 4 implementation. The Scientific Lead approved the complete recommended baseline in `decisions/scientific-baseline-approval.md` on 2026-08-21.

No external service, CHARMM-GUI upload, HPC job, cloud resource, docking run, or MD run is authorized.

## 2. Smallest Safe Slice

The first implementation targets `charmm_gui_cgenff_v1` ligand-plus-system export. It adds no general MD workflow framework and does not modify the existing visualization `complexes` stage.

Architecture:

```text
best_pose_per_tag_by_engine.csv + explicit selection/config
                         |
                         v
              resolve exact pose record
                         |
     connectivity topology + proven atom mapping
                         |
                         v
       approved protonation/charge backend
                         |
      explicit prepared receptor PDB/mmCIF
                         |
                         v
        deterministic writer -> strict validator
                         |
                         v
     atomic publish + row/aggregate provenance
```

## 3. Modules

- Add `post_docking_analysis/md_inputs.py` for request/config models, source resolution, export orchestration, atomic publication, and manifest writing.
- Add `post_docking_analysis/md_chemistry.py` for a narrow backend protocol and the approved Open Babel implementation. Backend absence is a typed capability result.
- Extend `post_docking_analysis/complex_validation.py` with separate MD-profile validators. Keep the existing lightweight visualization validator behavior compatible.
- Extend `post_docking_analysis/multi_engine_pipeline_impl.py` with `best_poses_by_engine`, `md_inputs_config`, and `md_inputs` artifacts plus `_dag_compute_md_inputs_node()`.
- Extend `workflow/cli.py` and its execution dispatch with `analyze md-inputs`. Keep `post_docking_analysis.cli` compatibility only where the existing entry point requires it.
- Add versioned JSON schemas under `post_docking_analysis/schemas/` for the aggregate manifest and row provenance.

## 4. DAG and Cache Contract

The graph node declares the authoritative selection CSV, normalized configuration artifact, pair/project metadata, selected pose files, topology files, receptor files, and backend/method version material as inputs. Because the current generic `ArtifactGraph` cache hashes mtime and size, `md_inputs` performs its own SHA-256 provenance comparison before reuse. A later generic content-hash refactor is out of scope.

The node is registered as optional and is addressable as scope `md_inputs`, but it is not added to `reports` or `full`. Direct CLI dispatch inspects the node and row manifests and returns nonzero for any required failure.

## 5. Implementation Sequence

1. Freeze request, row-manifest, provenance, and validation schemas.
2. Implement exact selection and source-lineage resolution with no chemistry writes.
3. Implement strict PDB parsing/writing and Tier A validation.
4. Implement the approved chemistry backend and Tier B validation.
5. Implement atomic per-row publication and aggregate manifest.
6. Register DAG paths/node/scope and CLI command.
7. Add contract, fixture, cache, dependency, and consumer-profile tests.
8. Run Spec Kit gates 5–10 and present evidence for Gate 11 review.

## 6. Verification Commands

Planned commands, adjusted to the local environment when implementation begins:

```text
python -m compileall -q post_docking_analysis workflow test
pytest -q test/test_spec032_md_inputs.py
pytest -q test/test_spec031_*.py
python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix
pytest -q
git diff --check
```

Dependency-backed consumer fixtures are additionally run in an environment containing the pinned approved chemistry backend. A skipped dependency test is reported as `skipped_missing_dependency`, never as a pass.

## 7. Rollback and Compatibility

All new output is confined to `5-Analysis/md_inputs/` and `.meta/md_inputs_config.json`. Existing complex, interaction, consensus, and top-pose outputs are unchanged. Removing the node, CLI route, and new output folders restores prior behavior without data migration.
