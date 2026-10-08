# Decision needed — PDB2PQR heavy-atom changes vs the conservation gate

**Status:** decided 2026-10-08 by Scientific Lead — option 1. Implemented in `docking/preparation/receptor_preparation.py` (Spec 033 R2c).

## Observation (2026-10-08, PDB2PQR 3.4.1, AMBER, pH 7.4, Meeko 0.7.1)

The input was 1IEP chain A, protein only (2229 heavy atoms, residues 225–498). With default PDB2PQR settings, the strict receptor module refuses the result because heavy atoms changed:

- **Added:** C-terminal carboxylate oxygen (OXT) on GLN 498. It is absent from the crystal chain because the chain is truncated.
- **Moved (hydrogen-bond optimisation flips):** His 295, His 375, His 396 (ring atoms swapped), Asn 414 (OD1/ND2), Gln 252 (OE1/NE2), Thr 272 (OG1/CG2).

With `--noopt --nodebump`, every flip disappears and only the terminal-OXT addition remains. With PDB2PQR disabled, Meeko's residue templates give 470 HD atoms, conservation passes, and the protonation is template-based rather than pH-titrated.

## Options

1. **(Agent recommendation)** Run PDB2PQR with `--noopt --nodebump` (titration at the explicit pH, no side-chain flips). Have the conservation gate accept **only** terminal atoms that PDB2PQR adds to complete a chain terminus, listing each one in provenance. Any other added, missing or moved heavy atom still fails. This keeps crystal heavy-atom geometry and gets pH-dependent protonation.
2. Allow PDB2PQR flip optimisation as well, recording each flip as an accepted heavy-atom move. This is closer to standard PDB2PQR practice, but it changes crystal side-chain atom identities.
3. Keep PDB2PQR off by default (Meeko templates, `template_selected_not_ph_titrated`), as in `main@7fc130a`. Enable PDB2PQR per project.

Until this decision is made, the branch default (`receptor_use_pdb2pqr: true`) refuses most crystal receptors whose chain ends without OXT. That is the fail-closed behaviour.

## Decision

Option 1 (Scientific Lead, 2026-10-08). The implemented rule is:

1. **PDB2PQR invocation** (receptor_use_pdb2pqr=true): `pdb2pqr30 --ff <force_field> --with-ph <ph> --noopt --nodebump --keep-chain <in> <out.pqr>`. The flags are recorded in provenance as `pdb2pqr_options`. `--keep-chain` is an addition to the decided pair. It only labels chains so that terminal additions can be attributed. It was verified not to move atoms: heavy-atom count and positions are unchanged.
2. **Removed or moved heavy atoms** always fail (`heavy_atom_conservation_failed`).
3. **Added heavy atoms** are accepted only if all of these hold. The atom is an oxygen named `O` or `OXT` (after PDB2PQR, OpenBabel and Meeko). It sits on the highest-numbered residue of its chain in the input ATOM records. Internal chain breaks are not termini. There is at most one such addition per terminus. Accepted additions are listed as `accepted_terminal_additions` and set `heavy_atom_conservation` to `passed_with_terminal_additions`. Any other addition fails and is named in the details.
4. **receptor_use_pdb2pqr=false** keeps the strict rule: no additions are accepted.

## Verified result (1IEP chain A, protein only, default config)

- `protonation_status`: `pdb2pqr_applied`; `requested_ph`: 7.4; force field AMBER.
- Accepted addition: chain A, GLN 498, atom O (the C-terminal carboxylate oxygen, PDBQT type OA).
- Heavy atoms: 2229 in the input, 2230 in the output.
- HD atoms: 472 in total, 422 within 1.1 Å of N or NA.
- Output sha256 `1062c9940d194ba90e778eb7269a5f75fe865bd875159ec93cfd3e9058ce6387` (scratch run; sha256 depends on the Meeko and PDB2PQR versions).

Software evidence and scientific evidence are reported separately. The His, Asn, Gln and Thr side chains keep their crystal orientation, because PDB2PQR optimisation is off. Histidine tautomers are not experimentally confirmed. The pKa and protonation states are predictions.
