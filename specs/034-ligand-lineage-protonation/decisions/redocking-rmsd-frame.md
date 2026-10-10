# Decision needed — redocking RMSD frame (in place vs Kabsch-superposed)

**Status:** open; Scientific Lead decision required. No code change is made under Spec 034.

The accepted Spec 031 method (`spec031-atom-mapping-v1`) computes RMSD after Kabsch superposition. Redocking validation conventionally reports **in-place** heavy-atom RMSD: the docked pose and the crystal ligand already share the receptor frame, and superposition can hide a misplaced or rotated pose.

On the Phase 2 1IEP/STI poses, a diagnostic in-place RMSD placed the top pose within the conventional 2 Å threshold and the other poses more than 12 Å away. Superposed values would be smaller and would not separate correctly from misplaced poses.

## Options
1. **(Agent recommendation)** For redocking validation only, report in-place RMSD as the primary value and the Kabsch value as secondary. Inter-engine geometric consensus keeps Kabsch, because there the frames are also shared but pose-shape agreement is the question being asked.
2. Keep Kabsch everywhere (status quo) and document the limitation.
