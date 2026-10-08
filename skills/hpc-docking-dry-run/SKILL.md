---
name: hpc-docking-dry-run
description: Generate and validate HPC submission scripts (Slurm, PBS, LSF) in dry-run mode without submitting jobs.
---

# HPC Docking Dry-Run Skill

## Trigger Conditions
Activate when generating HPC job scripts, configuring cluster profiles, or preparing batch docking submissions.

## Explicit Inputs and Outputs
- **Inputs:** HPC profile configuration, scheduler type, resource specifications (CPUs, GPUs, walltime, partitions).
- **Outputs:** Generated job scripts, staging commands, submission dry-run report.

## Refusal-to-Guess Rules
- **NEVER** submit live jobs to an HPC cluster without explicit human authorization.
- **NEVER** hardcode cluster credentials, tokens, or private scratch paths.

## Verification Commands
- `python main.py dock dry-run --help`
- `python main.py hpc profile --help`
