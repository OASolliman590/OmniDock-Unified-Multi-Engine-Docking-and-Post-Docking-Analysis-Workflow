---
name: statistical-review
description: Review statistical metrics, normalization models, ranking algorithms, and consensus scoring for mathematical correctness.
---

# Statistical Review Skill

## Trigger Conditions
Activate when reviewing consensus score aggregation, z-score transformations, percentile rankings, or correlation metrics.

## Explicit Inputs and Outputs
- **Inputs:** Score distributions, ranking tables, consensus configurations.
- **Outputs:** Statistical validity report, distribution diagnostics, tie-handling verification.

## Refusal-to-Guess Rules
- **NEVER** pool disparate energy scales without explicit normalization.
- **NEVER** present target percentile rankings as calibrated binding free energies.

## Verification Commands
- Validate score distributions using `scipy.stats` / `numpy` in Python.

## Known Limitations & Data-Use Notes
- Assesses statistical validity of the computational output.
