---
name: engine-parity-review
description: Verify multi-engine feature parity, parser compatibility, and execution contracts across GNINA, Vina, Smina, and AutoDock4.
---

# Engine Parity Review Skill

## Trigger Conditions
Activate when modifying engine runners, log parsers, pose extractors, or HPC adapters to ensure all 4 engines maintain parity.

## Explicit Inputs and Outputs
- **Inputs:** Engine runner modules (`docking/runners/`), docking output files (PDBQT/SDF/LOG).
- **Outputs:** Parity matrix report, parsed score tables, engine failure classifications.

## Refusal-to-Guess Rules
- **NEVER** route non-GNINA outputs through a GNINA-specific parser (e.g. Vina PDBQT vs GNINA log).
- **NEVER** assume an engine succeeded if its log or output pose file is missing.

## Verification Commands
- `python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix`
