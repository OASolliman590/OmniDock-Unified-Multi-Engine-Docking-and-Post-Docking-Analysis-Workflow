# Gate 2 Decision — Binding-Site Center and Provenance

**Decision ID:** `binding_site_center_v1`  
**Status:** Approved  
**Approved by:** Scientific Lead  
**Approval date:** 2026-08-21  
**Implements:** SR-003 / T016

## Approved Policy

1. For a holo structure with an explicitly selected co-crystal ligand, `binding_site_center` is the arithmetic centroid of that ligand's heavy-atom Cartesian coordinates in Å.
2. The ligand is selected from explicit accession/model/chain/residue metadata. Filename resemblance does not select it.
3. Hydrogens are excluded. Alternate locations use highest occupancy, then altloc `A`, then lexical order. The selected model is explicit; model 1 is used only for a single-model structure or a documented compatibility default.
4. `ligand_centroid`, `contact_residue_centroid`, and `binding_site_center` remain distinct fields. Contact-residue centroids may be reported as descriptors but do not replace the approved center.
5. For apo or no-reference cases, the workflow requires explicit user coordinates and their provenance. Without them, the center stage is `skipped_missing_configuration` and no docking box is invented.
6. No cavity detector or PLIP contact-atom center is approved in this scope. Contact-atom centering remains deferred until atom/ring identifiers have a versioned parser contract.
7. Existing box dimensions, exhaustiveness, seeds, scoring modes, and all other docking parameters are unchanged by this decision.
8. `overall_center` and `ligand_center` remain compatibility aliases for one migration window, but both callers must consume `binding_site_center` as the authoritative field.

## Required Provenance

Record center method/version, coordinates and units, coordinate frame, source structure path and SHA-256, accession, model, chain, residue name/number/insertion code, heavy-atom count, hydrogen/altloc policy, parser/tool version, explicit-user source when applicable, compatibility aliases written, and status/reason.

## Acceptance Cases

- Holo fixture centroid equals the selected ligand heavy-atom centroid.
- Hydrogens and unselected altloc coordinates do not affect the center.
- Two ligands require explicit selection; no automatic first-ligand choice is allowed.
- Apo/no-reference input without user coordinates returns `skipped_missing_configuration`.
- Compatibility aliases equal `binding_site_center` during the migration window.
- Box dimensions remain byte-for-byte/field-for-field unchanged except for added provenance.
