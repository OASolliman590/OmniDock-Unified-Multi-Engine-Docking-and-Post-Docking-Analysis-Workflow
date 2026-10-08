# Gate 2 Decision — Reference-Ligand Classification and Migration

**Decision ID:** `explicit_reference_metadata_v2`  
**Status:** Approved  
**Approved by:** Scientific Lead  
**Approval date:** 2026-08-21  
**Implements:** SR-004 / T018

## Default Scientific Policy

1. Explicit pair/extraction metadata is authoritative. Filename, tag, display-name, and broad text heuristics cannot establish scientific reference or self-redocking status.
2. A scientific reference requires a normalized explicit reference/co-crystal flag plus explicit ligand identity or source lineage. A self-redocking benchmark additionally requires explicit receptor/ligand source-structure correspondence.
3. Imported or manual pairlists may establish reference status only through those explicit fields and schema validation.
4. Missing or incomplete lineage is `unclassified` and non-reference, with a machine-readable warning. It cannot become a redocking validation row or a reference-affinity anchor.
5. `"False"`, `0`, empty, null, and NaN remain false. Names such as `novel_ligand_1` never establish reference status.
6. Current scientific outputs use method/version `explicit_reference_metadata_v2`.

## Legacy Reproduction Mode

- Legacy inference is opt-in and never the default.
- It emits method/version `legacy_inferred_reference_v1`, a prominent warning, the exact heuristic and source token, and a migration candidate report.
- Legacy inference may reproduce/display historical classification separately, but it cannot establish observed provenance, a current validation benchmark, or a reference-affinity anchor until the user confirms the candidate and writes explicit metadata.
- No historical project is silently rewritten. Migration is additive and reviewable.

## Consolidation Boundary

Reference/scientific-source classification must be implemented behind one shared policy API. Ligand display-name extraction for visualization or structure labeling is a separate non-scientific concern and may retain name parsing only when it cannot alter reference, validation, ranking, or confidence status.

## Required Provenance

Record classification method/version, normalized flag, identity and lineage fields used, source record/path and checksum when available, validation-anchor eligibility, legacy mode and heuristic when used, warnings, user-confirmation/migration state, and final classification reason.

## Acceptance Cases

- Fully explicit co-crystal metadata becomes a reference and is eligible for validation subject to source correspondence.
- Explicit flag without required identity/lineage becomes `unclassified` with an incomplete-metadata reason.
- `novel_ligand_1`, `_ligand_`, `"False"`, empty, and NaN remain non-reference.
- Legacy mode emits candidates/warnings but cannot create a current validation anchor.
- Default and legacy outputs are distinguishable and reproducible by method/version.
