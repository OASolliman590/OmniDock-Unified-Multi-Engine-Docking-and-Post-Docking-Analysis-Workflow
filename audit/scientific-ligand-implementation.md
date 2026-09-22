# Scientific ligand implementation record

Date: 2026-09-14

Baseline: `8d4434d3a33eb83b1e12cad82944b02c83270e47`

This bounded change covers ordinary single-component ligand preparation,
direct preparation configuration, normalization provenance, and the final
PDBQT chemistry contract.

## Implemented guards

- `validate_preparation_ph` rejects nonfinite values and JSON/Python boolean
  values. The direct environment seam raises for invalid pH instead of using
  the default. Profile aliases and omitted defaults remain supported, while an
  unknown public profile token is reported as invalid with the accepted list.
- Every RDKit-loaded ordinary ligand is required to contain one connected
  component before an external backend is invoked. No salt stripping,
  component selection, or enumeration is performed.
- Normalization records a source-to-normalized heavy-atom mapping using RDKit's
  stereochemistry-aware substructure matcher. Element, isotope, radical,
  aromaticity, heavy-atom connectivity, and bond order are checked explicitly.
  RDKit's `useChirality=True` match plus mapping-aware defined-stereo state
  protects stereochemistry; undefined source stereo cannot be silently
  assigned.
- Hydrogen/formal-charge differences are recorded per mapped atom. A selected
  pH procedure may permit them, but every formal-charge change must have the
  same per-atom hydrogen delta. A charge-only or otherwise unpaired change is
  rejected. Mapping ambiguity above the safety cap fails closed.
- Source 3D coordinates are compared at three-decimal precision for mapped
  heavy atoms. Existing coordinates are reported as `source_3d_preserved`;
  Open Babel-generated conformers are explicitly marked
  `openbabel_generated_nondeterministic` and are not treated as native. Every
  normalized and direct-preparation input artifact must also have one finite
  3D conformer; graph-only identity checks remain available to callers that
  do not yet have coordinates.
- Normalized SDF candidates and final PDBQT candidates are staged and checked
  before publication. Successful normalized artifacts are retained beside the
  prepared output, with input, normalized-SDF, prepared-input, and output
  SHA-256 provenance in the preparation summary.
- Meeko-mapped PDBQT output must match the prepared input graph, hydrogen state,
  defined stereochemistry, and mapped coordinates at PDBQT precision. Missing
  authoritative mapping is a limited legacy status for Open Babel/ADT output;
  malformed/conflicting mapping is invalid. Open Babel/ADT output therefore
  cannot be labeled exact chemistry merely because it is syntactically valid.

## Verification

Local execution record (the final chemistry command could not collect tests):

```text
python -m py_compile docking/models.py docking/preparation/ligand_identity.py docking/preparation/ligand_preparation.py test/test_scientific_preparation_config.py test/test_scientific_ligand_contract.py
python -m pytest test/test_scientific_preparation_config.py -q  # 15 passed
python -m pytest test/test_dockforge_smoke.py -q -k "profile_validator or ligand_output_contract_validator or prepare_ph_guard or prepare_guard"  # 4 passed
python -m pytest test/test_docking_correctness.py -q  # 32 passed
python -m pytest test/test_scientific_ligand_contract.py -q  # collection blocked by host RDKit DLL policy
```

The chemistry module intentionally imports RDKit strictly. On this Windows
host, `rdkit.Chem` cannot load `rdCIPLabeler` because of Application Control;
the new chemistry contract tests therefore require the existing chemistry CI
environment. No chemistry package was installed or bypassed locally. No real
Open Babel, Meeko, ADT, docking engine, remote benchmark, or NMRBox execution
was performed.
