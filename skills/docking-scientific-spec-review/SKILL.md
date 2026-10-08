---
name: docking-scientific-spec-review
description: Review molecular docking specifications for binding pocket definitions, protonation policies, scoring function assumptions, and validation gates.
---

# Docking Scientific Spec Review Skill

## Trigger Conditions
Activate when drafting or modifying docking parameter presets, scoring algorithms, or redocking validation criteria.

## Explicit Inputs and Outputs
- **Inputs:** Docking parameter schema, pocket coordinates, reference ligands.
- **Outputs:** Scientific audit checklist, scoring compatibility notes, parameter validation.

## Refusal-to-Guess Rules
- **NEVER** treat GNINA CNN scores, Vina affinities, and AD4 energies as directly comparable without explicit normalization.
- **NEVER** silently alter exhaustiveness, seed, or pH defaults.

## Verification Commands
- Check against `docking/parameter_schema.py` and `docking/models.py`.
