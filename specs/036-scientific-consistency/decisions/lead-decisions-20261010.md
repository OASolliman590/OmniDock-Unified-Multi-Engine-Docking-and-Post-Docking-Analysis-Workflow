# Scientific Lead decisions — 2026-10-10 (after the 1IEP/STI re-dock)

1. **Spec 032 T021: accepted.** Basis: the 1IEP/STI re-dock evidence. CHARMM-GUI/CGenFF consumer evidence remains `reported_unverified`.
2. **Charge-blind parent mapping: approved.** For atom correspondence only, in redocking and MD-export G3. Heavy-atom positions are charge-independent. Every result records `charge_blind_parent_mapping`.
3. **Legacy consensus modes: removal stands.** The Lead's direction is "the pipeline must be seamless and coherent". A single consensus method (`consensus_rank_geometry_qc_v2`, plus `single_engine_native_v1` for one-engine projects) serves that. This supersedes the Spec 031 `consensus-v2.md` note about keeping `dockbox_geometric`. Past results remain reproducible from commits before `a2154ef`.
4. **pH: one governing value.** The project protonation policy's pH governs every stage: receptor preparation, the ligand `ph_model`, and the MD export. The default MD path reuses the docked microspecies, which was produced under that policy, so no second pH is applied.
5. **Basic mode refuses explicit `--exhaustiveness`**, as it already refuses `--num-modes`. The effective value and its source are recorded.
6. **"Top overall" and best-per-protein summaries use the v2 consensus rank.** Raw cross-engine affinities are no longer compared.
7. **Merge into `main` once everything works**, through the Spec 035 per-change reconciliation.
8. **Design direction:** a layered protonation policy for multi-protein, multi-ligand projects, and the ligand radius of gyration computed once at ligand preparation. Both are specified in Spec 037.
