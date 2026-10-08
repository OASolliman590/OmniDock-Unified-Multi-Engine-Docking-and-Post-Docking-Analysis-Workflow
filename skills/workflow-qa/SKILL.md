---
name: workflow-qa
description: Enforce quality assurance, stage dependency checks, cache validation, and output topology across pipeline runs.
---

# Workflow QA Skill

## Trigger Conditions
Activate after running any preparation, docking, or post-docking analysis stage to verify outputs.

## Explicit Inputs and Outputs
- **Inputs:** Project directory, `.workflow/state.json`, `run_tracking/`, generated artifacts.
- **Outputs:** Stage status report, missing artifact list, DAG cache audit.

## Refusal-to-Guess Rules
- **NEVER** mark a stage as `completed` if the required report, summary CSV, or visualization is absent.
- **NEVER** treat an error or empty table as a successful skip.

## Verification Commands
- `pytest test/test_webui.py`
- `python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix`

## Known Limitations & Data-Use Notes
- QA confirms contract adherence; does not replace biological evaluation.
