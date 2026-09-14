# Independent scientific/code review

Reviewed the scientific-integrity working diff against baseline `8d4434d3a33eb83b1e12cad82944b02c83270e47` on 2026-09-14. Scope: preparation configuration, structure selection, receptor preparation, ligand preparation/identity, pose geometry, consensus and redocking validation; associated scientific fixtures and methods documentation.

## Result

No unresolved material P1/P2 findings in the inspected final implementation. This is a bounded source review, not a chemistry-runtime validation or approval to merge without CI. Windows Application Control prevents local RDKit execution; chemistry-required fixtures must run in the declared CI environment. Independent tester results are recorded separately in `scientific-verification.md`.

## Findings resolved during review

- RDKit normalization staged a `.candidate` file that its suffix-dispatch loader could not read. Staging now retains `.sdf`, with a positive 2D-to-3D normalization fixture.
- A source without 3D coordinates could bypass conditional coordinate comparison and receive an exact mapped-output claim. Prepared normalized artifacts, direct preparation inputs and the mapped-output validator now require finite 3D input coordinates.
- Non-Meeko output summaries could declare `status=exact` with `is_valid=false` and pass the public contract. Exact validation now requires a true validity result for every backend.
- Durable normalized artifacts were published before backend/output validation. Their publication now follows graph, coordinate, syntax and mapped-output checks.

Also inspected the corrected source/query mapping direction, mapped-label stereo ordering, paired per-atom proton/charge changes, explicit heavy element/isotope/radical/bond comparisons, absence-versus-malformed mapping status, and removal of native-provenance inference from arbitrary 3D coordinates.

## Evidence and limits

Positive and negative fixtures cover source identity, a real chiral protonation example, non-self-inverse atom reordering, stereo/graph mutations, mapped explicit-H output, final state/coordinate changes, receptor atom identity and frame preservation, exact pose record selection, coherent PDBQT metadata scopes and consensus input ambiguity. Inspection supports the intended assertions; chemistry tests were not executed by this reviewer.

Receptor checks conserve represented atom identities/elements and coordinates while distinguishing upstream titration from final pH validation. Methods documentation does not claim biological protonation correctness, physical pose validation, calibrated affinity, or a complete replicate policy. Symmetric mapping choice remains conservative and may reject otherwise equivalent ambiguous inputs. Public summary validation checks supplied assertions and artifact hashes; it is not an independent re-execution of chemistry from an untrusted summary.

Final handoff inspection confirmed the backend-failure fixture preserves both preexisting PDBQT and normalized SDF bytes after successful staged normalization and a backend exception. git diff --check returned exit code 0 (line-ending warnings only). No further material findings.

## Full-suite compatibility follow-up

The subsequent full filtered suite exposed three existing smoke cases that pass `expected_engines=[]` to request inferred engine scope. The initial independent source review missed this established default-call compatibility; the new empty-scope rejection was a regression.

Reviewed the bounded correction: `None` and an empty engine collection infer the observed normalized engine scope, while every nonempty declared scope still rejects malformed tokens, duplicate normalized engines and engines outside that scope. The added two-engine positive fixture checks agreement count 2 and fraction 1.0. The preceding row-identity, scoring-function and duplicate-row guards are still executed before inference. No new scientific denominator claim is introduced, and no remaining material finding was identified in this correction. Execution results belong to the independent tester's verification record; chemistry CI remains required.
