---
name: docking-provenance-report
description: Generate comprehensive scientific provenance reports detailing parameters, software versions, input checksums, and consensus methods.
---

# Docking Provenance Report Skill

## Trigger Conditions
Activate when generating final study reports, thesis chapters, or publication supplementary materials.

## Explicit Inputs and Outputs
- **Inputs:** Project outputs, run manifests, stage logs, consensus results.
- **Outputs:** Markdown/HTML provenance report, methods text with citation placeholders.

## Refusal-to-Guess Rules
- **NEVER** claim experimental binding affinity from computational docking scores alone.
- **NEVER** omit docking box geometry or engine scoring functions from methods descriptions.

## Verification Commands
- Inspect generated reports in `7-Reports/`.
