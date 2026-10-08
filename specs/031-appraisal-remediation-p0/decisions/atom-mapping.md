# Gate 2 Decision — Atom Mapping and RMSD Comparability

**Decision ID:** `spec031-atom-mapping-v1`  
**Status:** Approved  
**Approved by:** Scientific Lead  
**Approval date:** 2026-08-21  
**Implements:** SR-001 / T012

## Approved Policy

1. RMSD uses heavy atoms only. Explicit and implicit hydrogens are excluded from both mapping and the denominator.
2. Atom correspondence must be proven from a shared canonical ligand topology/lineage or from a chemistry-aware, bond-order-aware graph isomorphism.
3. The graph contract includes element, normalized bond order/aromaticity, formal charge, and stereochemical labels when the source format preserves them. If required labels are absent and shared lineage does not prove identity, the comparison is `not_comparable`.
4. All symmetry-equivalent graph mappings are eligible. The reported value is the minimum heavy-atom RMSD after Kabsch alignment across valid mappings.
5. Exact atom order is only a fast path when atom identity and topology lineage are proven. Element sorting, nearest-neighbor matching, and centroid-radius sorting cannot establish scientific correspondence.
6. Full heavy-atom coverage is required. Atom-count mismatch, constitutional isomers, topology/protonation mismatch, ambiguous incomplete mapping, or missing bond topology are `not_comparable`. No common-substructure RMSD is approved in this version.
7. Alternate locations use highest occupancy, then altloc `A`, then lexical order; the selected model, conformer, and alternate-location rule are recorded.
8. Inter-engine geometric agreement is `RMSD <= 2.0 Å`. Redocking is `pass` at `<= 2.0 Å`, `warn` at `> 2.0 Å and <= 3.5 Å`, and `fail` at `> 3.5 Å`.
9. A numeric RMSD must never be emitted for `not_comparable` mappings. Centroid-based values may remain only as explicitly labeled non-scientific diagnostics and cannot affect agreement, validation, ranking, or confidence.

## Topology Source Priority

1. Explicit prepared-ligand topology plus retained atom lineage and checksum.
2. Explicit reference/docked chemistry files whose graphs can be parsed and proven isomorphic.
3. Otherwise, `not_comparable` with a machine-readable reason.

PDBQT coordinates without a provable external/shared topology do not establish bond order by themselves.

## Required Provenance

Each comparison records method/version, topology sources and SHA-256 values, parser/backend version, atom policy, mapped and total heavy atoms, coverage, valid mapping count, selected mapping identifier, alignment method, cutoff, result status, RMSD when valid, and failure reason when invalid.

## Dependency Policy

A narrowly scoped and tested chemistry dependency is authorized only if an implementation spike demonstrates that the declared Open Babel environment cannot reliably enumerate the approved graph/symmetry mappings. The selected dependency must be version-bounded and documented in every applicable environment manifest. If the backend is unavailable, the stage reports `skipped_missing_dependency` and comparisons report `not_comparable`; heuristic fallback is prohibited.

## Compatibility

This decision supersedes Spec 027's centroid-sorted soft-agreement semantics for scientific outputs. Legacy fields may be retained for reproduction only when labeled with their legacy method and excluded from current agreement decisions.

## Acceptance Cases

- Permuted atom order maps and produces the expected RMSD.
- Symmetry-equivalent atoms produce the minimum valid RMSD deterministically.
- Same-formula constitutional isomers are `not_comparable`.
- Atom-count mismatch is `not_comparable`.
- Missing topology/backend produces no numeric RMSD.
- Proven atom-order identity takes the fast path with the same result as graph mapping.
