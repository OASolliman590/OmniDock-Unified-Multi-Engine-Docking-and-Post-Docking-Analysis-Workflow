# OmniDock Cleanup Progress

**Fixed base:** `8d4434d3a33eb83b1e12cad82944b02c83270e47`
**Branch:** `codex/cleanup-simplify-docs`

| Batch | Status | Commit | Verification |
| --- | --- | --- | --- |
| Repository checker | queued | — | — |
| Evidence-preserving hygiene | queued | — | — |
| Docking deployment refactor | queued | — | — |
| CLI parser extraction | queued | — | — |
| Documentation consolidation | queued | — | — |
| Final review/distribution | queued | — | — |

## Needs your eyes

- Code Simplify is not callable or installable in this session; Grok plus the installed codebase-design and code-review skills are the explicit fallback.
- The separate scientific-integrity task owns behavior changes in preparation and post-docking execution. This cleanup owns shared public documentation and non-scientific orchestration.
- Windows Application Control currently blocks RDKit's `rdinchi` DLL, so a clean Linux/Windows CI run is required before the draft PR can be considered fully verified.

## Baseline

- 300 tracked files; 13,327,496 tracked bytes.
- Five broken relative Markdown links.
- Five generated audit log/XML artifacts totaling 184,736 bytes.
- One tracked dangling private absolute symlink.
- Host-usable regression: 145 passed; two RDKit-dependent smoke failures. Three RDKit-focused modules fail collection because the same DLL is blocked.

## End-of-run checklist

- Full diff reviewed against the design specification.
- Repository, test, build, wheel, shell, and link gates pass or exact host blockers are recorded.
- Independent Standards and Spec reviews have no unresolved blocking findings.
- Draft PR is open with CI results and scientific-integrity sequencing; branch is not merged.
