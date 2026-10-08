---
name: docking-parameter-audit
description: Audit docking parameters (box coordinates, padding, exhaustiveness, num_modes, seed, engine extensions) for consistency and correctness.
---

# Docking Parameter Audit Skill

## Trigger Conditions
Activate before launching docking jobs to verify parameter resolution and prevent unrecorded overrides.

## Explicit Inputs and Outputs
- **Inputs:** User parameter dict, preset name (`screening_fast`, `balanced`, `exhaustive`), engine name.
- **Outputs:** Resolved parameter schema JSON, preflight validation status.

## Refusal-to-Guess Rules
- **NEVER** silently swap advanced settings for basic presets.
- **NEVER** omit the random seed or engine version from resolved parameter metadata.

## Verification Commands
- `python -c "from docking.parameter_schema import resolve_parameter_schema; print(resolve_parameter_schema('balanced'))"`
