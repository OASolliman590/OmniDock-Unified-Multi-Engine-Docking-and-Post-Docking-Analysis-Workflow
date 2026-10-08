# Spec 031 Scientific Interpretation Contract

Docking outputs in OmniDock are computational hypotheses and relative ranking signals. They are not experimental proof of affinity, kinetics, selectivity, mechanism, safety, or clinical efficacy.

## Scores and classifications

- Vina/Smina empirical affinity, AutoDock4 estimated binding energy, and GNINA CNN outputs are not directly commensurate across engines.
- Raw scores must not be compared across different protein targets as though they share one calibrated physical scale.
- Any consensus output must name the engines included, failed, or missing and record its normalization, formula version, direction, weights/components, missing-data policy, tie rule, and pose-selection rule.
- `classification_basis` is part of the result contract. `target_percentile_fallback` means the label is cohort-relative, not a calibrated activity threshold. `insufficient_n` is unclassified rather than successful.

## Pocket descriptors

The following outputs are explicitly heuristic descriptors:

- fixed 5 Å sphere volume;
- charged-residue counts and empirical charge sums;
- hydrophobic-residue counts;
- the combined legacy “druggability” score.

These values are not calibrated physical pocket volumes, electrostatic potentials, binding free energies, or validated druggability probabilities. New outputs carry `descriptor_contract`, method fields, `druggability_score_calibrated=false`, and `scientific_limitation`.

## Deployment mode

Deployment asset generation records the resolved scheduler `submit_mode`, selected pair source, and pair-source SHA-256. `--pairlist-file` scopes generation to exactly the supplied rows and does not rewrite the canonical pairlist or pair-curation state. Generating assets does not itself authorize remote submission.

## Statuses

Required stages distinguish `completed`, `failed`, `skipped_disabled`, `skipped_missing_dependency`, `skipped_network`, `skipped_missing_configuration`, and `not_comparable`. A missing required artifact is never treated as success.
