# Verification history

This file is a curated **historical** record. It is not a current test, build,
engine, or scientific-validation run.

Two layers are kept distinct:

| Layer | What it is | What it is not |
| --- | --- | --- |
| Immutable historical evidence | Outcomes recorded on **2026-09-12** from commit `cfe33056c77d817d19439a9e9affe225d00ca488` / fixed review base `8d4434d3a33eb83b1e12cad82944b02c83270e47`. Raw generated logs and XML that used to live under `audit/` were removed from the working tree in this hygiene batch; their Git blob hashes remain recoverable from that fixed base. | A rerun. These numbers were **not** reproduced as part of this cleanup. |
| Current validation | `python scripts/check_repository.py`, focused contract tests, Bash `bash -n` on remaining tracked shell scripts, `python -m compileall`, and `git diff --check` for this hygiene batch. CI is configured to run the repository checker after pytest and before the distribution build. | Proof that docking engines, cluster submission, CUDA, preparation backends, or ProLIF work. Not a scientific benchmark. |

Do not treat the 2026-09-12 counts as current CI status. Do not treat deletion of
legacy helpers or generated logs as scientific validation of preparation,
docking, or post-docking results.

Fixed review base: `8d4434d3a33eb83b1e12cad82944b02c83270e47`.

## Historical local verification record — 2026-09-12

These facts are copied from the preserved historical claims; they were **not**
rerun for this batch.

- Host/runtime: Windows, Python 3.12.10.
- Full pytest suite: **225 passed, 1 skipped, 0 failures**, 586 warnings, 209.78 seconds (226 collected tests).
- The skip was the real Meeko macrocycle export test: Windows Application Control blocked RDKit `rdDetermineBonds`. That integration remains unverified on that host.
- Four docking command/configuration dry runs passed.
- Source and wheel distributions built successfully.
- Every console entrypoint imported from the extracted `pdb_prepare_wizard-3.0.1-py3-none-any.whl` outside the checkout.
- Real engines, cluster submission, CUDA, preparation backends, and ProLIF were **not** validated by that historical record.
- The configured Linux/Windows CI matrix had not been run remotely as part of that record.

Raw artifacts that previously held those details:

- `audit/correctness-tests.log` (`54fe4a9cb7e83e39d80f6ed43699d76d357d4066`)
- `audit/correctness-tests.xml` (`a001145a26077519937918237bfa2b6dc2807bd8`)
- `audit/build-correctness.log` (`1096b5d64e9c37a1f59698c8dbce21392ab5d076`)
- `audit/docking-dry-run.log` (`83974e86779979cf05d1c9c6617ff6226b0599c8`)
- `audit/wheel-verification.log` (`eb19abe30872cefeffdbaa256d628db68d78da00`)

## Binary-safe recovery

The only verified binary-safe materialization command in this checkout is
restore of the **original** path. Git writes the blob bytes itself:

```text
git restore --source=8d4434d3a33eb83b1e12cad82944b02c83270e47 -- <original-path>
```

Inspect blob identity (hash / contents on the terminal). This is **not** a
PowerShell binary-write recipe:

```text
git cat-file blob <hash>
```

If the working tree should then hold the file at a **new** location, restore
the original path first and use an ordinary byte-preserving filesystem move
(for example Explorer, `Move-Item`, or `shutil.move`). Do not recover
`.docx`, `.png`, `.svg`, `.xml`, or other binary blobs with PowerShell text
redirection (`>`, `Out-File`, or `Set-Content` without a byte API). That
recodes newlines and corrupts the file.

## Removed and moved paths

Every path below is recorded with its original Git mode and blob hash at the
fixed base. Recovery always uses the **original** path and hash, even when the
working tree now has a new location.

| Original path | Mode | Blob hash | Disposition | Rationale | Recovery |
| --- | --- | --- | --- | --- | --- |
| `4-K_pneumoniae_/5-Post_Docking_Analysis` | `120000` | `4f1f79df14c6d15fc7c4e7ac4db10e383492579b` | deleted | Tracked symlink whose blob is the absolute private workstation target `/Users/omara.soliman/Desktop/Projects /My Projects/8-PDB-Prepare-Wizard/4-K_pneumoniae_/5-Analysis`. Unreferenced by code or documentation. Dangling in this checkout. Not a scientific result. | Restore original path: `git restore --source=8d4434d3a33eb83b1e12cad82944b02c83270e47 -- 4-K_pneumoniae_/5-Post_Docking_Analysis`. Inspect identity: `git cat-file blob 4f1f79df14c6d15fc7c4e7ac4db10e383492579b` |
| `audit/build-correctness.log` | `100644` | `1096b5d64e9c37a1f59698c8dbce21392ab5d076` | deleted (generated) | Generated local build log. Historical outcome summarized above. | Restore original path: `git restore --source=8d4434d3a33eb83b1e12cad82944b02c83270e47 -- audit/build-correctness.log`. Inspect identity: `git cat-file blob 1096b5d64e9c37a1f59698c8dbce21392ab5d076` |
| `audit/correctness-tests.log` | `100644` | `54fe4a9cb7e83e39d80f6ed43699d76d357d4066` | deleted (generated) | Generated pytest log for the 2026-09-12 Windows run. | Restore original path: `git restore --source=8d4434d3a33eb83b1e12cad82944b02c83270e47 -- audit/correctness-tests.log`. Inspect identity: `git cat-file blob 54fe4a9cb7e83e39d80f6ed43699d76d357d4066` |
| `audit/correctness-tests.xml` | `100644` | `a001145a26077519937918237bfa2b6dc2807bd8` | deleted (generated) | Generated JUnit XML for the same historical run. | Restore original path: `git restore --source=8d4434d3a33eb83b1e12cad82944b02c83270e47 -- audit/correctness-tests.xml`. Inspect identity: `git cat-file blob a001145a26077519937918237bfa2b6dc2807bd8` |
| `audit/docking-dry-run.log` | `100644` | `83974e86779979cf05d1c9c6617ff6226b0599c8` | deleted (generated) | Generated docking dry-run log. | Restore original path: `git restore --source=8d4434d3a33eb83b1e12cad82944b02c83270e47 -- audit/docking-dry-run.log`. Inspect identity: `git cat-file blob 83974e86779979cf05d1c9c6617ff6226b0599c8` |
| `audit/wheel-verification.log` | `100644` | `eb19abe30872cefeffdbaa256d628db68d78da00` | deleted (generated) | Generated wheel/entrypoint verification log. | Restore original path: `git restore --source=8d4434d3a33eb83b1e12cad82944b02c83270e47 -- audit/wheel-verification.log`. Inspect identity: `git cat-file blob eb19abe30872cefeffdbaa256d628db68d78da00` |
| `autodock/prep_autodock.sh` | `100644` | `d09d4fd3e86088cef6cf91851a9e84dc22c4cc12` | deleted (retired helper) | Failed Bash parsing (`bash -n`: syntax error near unexpected token `2` at `for mol in "$RAW_LIG"/*.{mol2,sdf} 2>/dev/null`). Neither packaged (`MANIFEST.in` / wheel) nor called by active Python code. Predated and bypassed the supported authoritative-graph / fail-closed preparation path. **Deletion is repository hygiene, not a scientific validation.** | Restore original path: `git restore --source=8d4434d3a33eb83b1e12cad82944b02c83270e47 -- autodock/prep_autodock.sh`. Inspect identity: `git cat-file blob d09d4fd3e86088cef6cf91851a9e84dc22c4cc12` |
| `autodock/prep_ligands_custom.sh` | `100755` | `e04ba8fd6939a325e16efa08ac21ea8def213ca4` | deleted (retired helper) | Failed Bash parsing at the same brace-expansion redirect. Neither packaged nor called by active code. Predated and bypassed the supported preparation path. **Deletion is repository hygiene, not a scientific validation.** | Restore original path: `git restore --source=8d4434d3a33eb83b1e12cad82944b02c83270e47 -- autodock/prep_ligands_custom.sh`. Inspect identity: `git cat-file blob e04ba8fd6939a325e16efa08ac21ea8def213ca4` |
| `Sertraline_Docking_Methods_Numbered_With_Citation_Placeholders.docx` | `100644` | `1ca6ff470de8355ba6ad3de544a93a5a65ec81b9` | moved to `closing_thesis/manuscripts/Sertraline_Docking_Methods_Numbered_With_Citation_Placeholders.docx` | Research manuscript. Bytes unchanged (destination `git hash-object` matches). | Restore original path: `git restore --source=8d4434d3a33eb83b1e12cad82944b02c83270e47 -- Sertraline_Docking_Methods_Numbered_With_Citation_Placeholders.docx`, then a byte-preserving filesystem move if relocating. Inspect identity: `git cat-file blob 1ca6ff470de8355ba6ad3de544a93a5a65ec81b9` |
| `Sertraline_Docking_Methods_Omni_DockForge.docx` | `100644` | `59d9b909cfc91a3b993c78610c30f247628c19e8` | moved to `closing_thesis/manuscripts/Sertraline_Docking_Methods_Omni_DockForge.docx` | Research manuscript. Bytes unchanged. | Restore original path: `git restore --source=8d4434d3a33eb83b1e12cad82944b02c83270e47 -- Sertraline_Docking_Methods_Omni_DockForge.docx`, then a byte-preserving filesystem move if relocating. Inspect identity: `git cat-file blob 59d9b909cfc91a3b993c78610c30f247628c19e8` |
| `Sertraline_Docking_Methods_Thesis_Chapter.docx` | `100644` | `707148a0c70dd727dd799d4712394876656f69ba` | moved to `closing_thesis/manuscripts/Sertraline_Docking_Methods_Thesis_Chapter.docx` | Research manuscript. Bytes unchanged. | Restore original path: `git restore --source=8d4434d3a33eb83b1e12cad82944b02c83270e47 -- Sertraline_Docking_Methods_Thesis_Chapter.docx`, then a byte-preserving filesystem move if relocating. Inspect identity: `git cat-file blob 707148a0c70dd727dd799d4712394876656f69ba` |
| `multiscale_overlapped_experimental_vs_triple_engine.png` | `100644` | `1666b829fa38e5d44392bb52c6308f90c8b2711e` | moved to `closing_thesis/figures/multiscale_overlapped_experimental_vs_triple_engine.png` | Research figure. Bytes unchanged. | Restore original path: `git restore --source=8d4434d3a33eb83b1e12cad82944b02c83270e47 -- multiscale_overlapped_experimental_vs_triple_engine.png`, then a byte-preserving filesystem move if relocating. Inspect identity: `git cat-file blob 1666b829fa38e5d44392bb52c6308f90c8b2711e` |
| `multiscale_overlapped_experimental_vs_triple_engine.svg` | `100644` | `858514b5471c66e96a3e4fb64ae33d1d304889dd` | moved to `closing_thesis/figures/multiscale_overlapped_experimental_vs_triple_engine.svg` | Research figure. Bytes unchanged. | Restore original path: `git restore --source=8d4434d3a33eb83b1e12cad82944b02c83270e47 -- multiscale_overlapped_experimental_vs_triple_engine.svg`, then a byte-preserving filesystem move if relocating. Inspect identity: `git cat-file blob 858514b5471c66e96a3e4fb64ae33d1d304889dd` |

Supported replacement for the retired shell helpers is `prep_autodock_enhanced.sh`
(packaged, syntax-valid) or the Python preparation workflow. See
`docs/dockforge_migration_notes.md`.

## Intentionally retained

Curated audit reports, executable probes, and tracked structured probe JSON
under `audit/` were **not** deleted. Closing-thesis manuscripts and figures were
**not** deleted. `audit/postdocking_probes.py` remains tracked; its
`postdocking-probe-artifacts/results.json` output was never tracked at the
audited or fixed base and is not invented here.
