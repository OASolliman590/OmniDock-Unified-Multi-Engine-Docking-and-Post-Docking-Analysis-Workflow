---
name: scientific-spec-review
description: Review specifications for scientific rigor, domain validity, biological assumptions, and explicit parameter declarations.
---

# Scientific Specification Review Skill

## Trigger Conditions
Activate before implementing any new algorithm, post-docking analysis stage, score transformation, or biological metric calculation.

## Explicit Inputs and Outputs
- **Inputs:** `spec.md`, proposed algorithm/formula, reference publications, target biological questions.
- **Outputs:** Scientific review checklist, identified ambiguities, explicit approval recommendation.

## Refusal-to-Guess Rules
- **NEVER** approve a specification that treats heuristic descriptors as physical binding affinities without explicit disclosure.
- **NEVER** allow cross-engine score averaging without a formal normalization model.
- **FLAG** any hardcoded thresholds lacking scientific citations or empirical justification.

## Verification Commands
- Check mathematical formulation against authoritative literature (e.g., Vina, GNINA, PLIP manuals).

## Known Limitations & Data-Use Notes
- Scientific review evaluates methodological defensibility, not experimental truth.
