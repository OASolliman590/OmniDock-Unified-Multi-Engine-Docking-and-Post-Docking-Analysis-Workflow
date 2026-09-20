# OmniDock Cleanup Progress

**Fixed base:** `8d4434d3a33eb83b1e12cad82944b02c83270e47`
**Branch:** `codex/cleanup-simplify-docs`

| Batch | Status | Commit | Verification |
| --- | --- | --- | --- |
| Repository checker | complete | `968b36d` | 9 passed, 1 native-symlink capability skip; 11 expected base findings; see [exact final evidence](2026-09-14-repository-cleanup-verification.md) |
| Evidence-preserving hygiene | complete | `3c1cb00` | checker/static gates clean; two shell scripts parse; five research assets are exact renames; see [evidence](2026-09-14-repository-cleanup-verification.md) |
| Docking deployment refactor | complete | `b7cb01f` | deployment/correctness contracts pass with expected compatibility warnings; see [evidence](2026-09-14-repository-cleanup-verification.md) |
| CLI parser extraction | complete | `5c9ff96` | parser/workflow and representative CLI/deployment smoke checks pass; see [evidence](2026-09-14-repository-cleanup-verification.md) |
| Documentation consolidation | complete | `391bd66` | 125 Markdown files, 0 broken relative links, 7/7 help surfaces, and five exact archive renames; see [evidence](2026-09-14-repository-cleanup-verification.md) |
| Final review/distribution | active | — | local static/build/wheel gates complete; independent fix re-review and draft PR/CI pending; see [evidence](2026-09-14-repository-cleanup-verification.md) |

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
