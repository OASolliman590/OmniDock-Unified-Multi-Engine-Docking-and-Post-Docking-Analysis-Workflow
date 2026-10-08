---
name: interaction-analysis-review
description: Verify protein-ligand interaction profiling contracts (PLIP, ProLIF, LigPlot+, PoseView, PandaMap).
---

# Interaction Analysis Review Skill

## Trigger Conditions
Activate when analyzing 2D/3D interaction profiles, hydrogen bonds, hydrophobic contacts, and pi-stacking interactions.

## Explicit Inputs and Outputs
- **Inputs:** Receptor-ligand complex PDB/PDBQT files, interaction engine configs.
- **Outputs:** Interaction fingerpints, 2D diagrams, PandaMap publication graphics.

## Refusal-to-Guess Rules
- **NEVER** claim an interaction exists without tool-verified geometric criteria (distance/angle).
- **NEVER** run only half of the joint interaction contract (ProLIF + LigPlot+) when joint analysis is requested.

## Verification Commands
- Inspect interaction outputs in `5-Analysis/` and `6-Visualizations/`.
