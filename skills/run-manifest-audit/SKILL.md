---
name: run-manifest-audit
description: Audit run manifests, stage statuses, and reproducibility metadata in `run_tracking/` and `.meta/`.
---

# Run Manifest Audit Skill

## Trigger Conditions
Activate before and after workflow execution to inspect or verify run tracking records.

## Explicit Inputs and Outputs
- **Inputs:** `run_tracking/manifest.json`, `run_tracking/stage_status.json`, `.meta/`.
- **Outputs:** Completeness audit, stage timing summary, failure diagnostics.

## Refusal-to-Guess Rules
- **NEVER** leave a failed stage unrecorded in `stage_status.json`.
- **NEVER** omit environment or container digests from the manifest.

## Verification Commands
- `python -c "import json; m=json.load(open('run_tracking/manifest.json')); print(m.keys())"`
