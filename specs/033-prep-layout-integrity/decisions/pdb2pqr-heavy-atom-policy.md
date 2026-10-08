# Decision needed — PDB2PQR heavy-atom changes vs the conservation gate

**Status:** open; Scientific Lead decision required.

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
