---
name: docking-fixture-test
description: Execute deterministic test fixtures and redocking validation suites against known holo complexes (e.g., 1IEP / STI-571).
---

# Docking Fixture Test Skill

## Trigger Conditions
Activate during test runs, continuous integration, or validation of preparation and parser changes.

## Explicit Inputs and Outputs
- **Inputs:** Test fixture structures (`test/fixtures/`, RCSB 1IEP), expected binding modes, reference RMSDs.
- **Outputs:** Fixture execution report, calculated RMSDs, score extraction metrics.

## Refusal-to-Guess Rules
- **NEVER** use unredistributable or undocumented data for fixtures.
- **NEVER** loosen RMSD tolerance thresholds merely to make a test pass.

## Verification Commands
- `python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix`
- `pytest test/test_webui.py`
