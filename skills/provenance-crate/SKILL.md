---
name: provenance-crate
description: Generate and verify RO-Crate and run-tracking provenance manifests containing code revisions, input hashes, parameters, and tool metadata.
---

# Provenance Crate Skill

## Trigger Conditions
Activate upon completion of docking or post-docking analysis to produce publication-grade provenance.

## Explicit Inputs and Outputs
- **Inputs:** `run_tracking/manifest.json`, input files, environment metadata, git commit hash.
- **Outputs:** Provenance manifest, SHA-256 hash table, RO-Crate metadata JSON.

## Refusal-to-Guess Rules
- **NEVER** omit tool versions, engine commit IDs, random seeds, or box dimensions.
- **NEVER** forge or backfill checksums for missing files.

## Verification Commands
- Inspect `run_tracking/manifest.json` and verify SHA-256 hashes against input files.

## Known Limitations & Data-Use Notes
- Records exact digital provenance for local and HPC executions.
