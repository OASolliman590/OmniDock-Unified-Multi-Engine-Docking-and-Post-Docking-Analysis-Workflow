# Scientific methods and limits

This document describes the scope of the scientific-integrity correction. The
current-code findings, reproductions and acceptance tests are recorded in
`audit/scientific-*-assessment.md`; those assessments distinguish new findings
from the earlier baseline audit. The pinned starting point and verification
sequence are in [the correction plan](../audit/scientific-integrity-plan.md).

## Separate scientific questions

OmniDock's evidence must be interpreted separately for each question:

| Question | Evidence required | What it does not establish |
|---|---|---|
| Chemical identity | Authoritative graph, atom correspondence, charge/isotope/stereochemistry checks | Correct biological protonation or tautomer state |
| Preparation execution | Input/output hashes, backend/configuration, successful format and chemistry checks | Universal suitability for all targets |
| Reference-pose reproduction | Exact selected pose, matching chemistry and receptor frame, symmetry-aware RMSD | Physical plausibility or experimental affinity |
| Physical pose quality | An explicitly executed chemical/geometric/steric validator with recorded settings | Experimental activity or selectivity |
| Relative prioritization | Compatible scoring context and explicit candidate/engine coverage | Calibrated binding free energy |
| Screening performance | A specified external evaluation protocol and experimental labels | Prospective performance on unrelated targets |

This table defines the interpretation policy for this correction. RMSD and
physical plausibility require separate assessment: PoseBusters evaluates
chemical consistency and intra/intermolecular geometry in addition to
reference-pose displacement. [Buttenschoen et al. (2024/02), methods and abstract]

## Preparation policy

The target-specific protocol must identify the structure model and assembly,
chains, alternate conformers, waters, metals, cofactors and covalent chemistry
to retain. It must document any repair, missing region or terminal treatment.
These are required protocol decisions, not biological conclusions inferred by
OmniDock. Unsupported cases should be reported with their reason rather than
resolved by silently deleting atoms or choosing a chemical state.

Hydrogen addition, titration-state assignment and partial-charge assignment are
separate operations. PDB2PQR selects its PROPKA calculation through
`--titration-state-method propka`; `--with-ph` supplies the pH used by the selected
calculation. Its source adds hydrogens separately and applies force-field
parameters later. [PDB2PQR (accessed 2026/09), `build_main_parser` and `non_trivial`]

The final docking backend must be recorded separately from any upstream
protonation procedure. A successful upstream calculation alone is insufficient
evidence that downstream templates preserved each assigned state. A methods
report should identify that uncertainty rather than describe the final receptor
as pH-validated.

The preparation contract distinguishes hydrogen changes from conserved heavy
atoms. A receptor atom map may permit serial renumbering while requiring each
chain/residue/atom identity, element and represented coordinate to remain the
same. This check does not establish a valid histidine, terminal, disulfide,
cofactor or metal-coordination state. Missing or ambiguous chemistry needs a
target-specific protocol and backend support.

Source-scoped frame identifiers are deliberately conservative. Selections made
from the same original PDB share its frame identifier; different original source
bytes are not automatically asserted to share a frame, even if their coordinate
strings happen to match. Reprepare related receptor/ligand artifacts together
when moving between frame-identifier versions. Existing files are not certified
merely because a historical sidecar has a matching text label.

Ordinary ligand preparation requires one connected component. Its normalization
ledger compares mapped heavy-atom identity, bonds and defined stereochemistry,
and records paired hydrogen/formal-charge changes under the requested pH
procedure. A retained normalized SDF and byte hashes identify the actual input
to the export backend. Mapped Meeko output must match that prepared chemical
state and its represented coordinates. An output without authoritative mapping
has limited validation scope; successful PDBQT syntax checks do not establish
exact output chemistry.

Preserved source 3D coordinates do not by themselves establish experimental or
native provenance. Open Babel conformer generation is not asserted to be
deterministic; the retained generated artifact identifies the realized input.
Symmetric mappings that cannot be resolved within the bounded search fail
conservatively. These checks do not enumerate tautomers or prove the selected
protonation state is biologically correct.

## Pose and score interpretation

The strict RMSD path requires identical isomeric molecular graphs and computes
RMSD in the existing coordinate frame across stereochemistry-compatible graph
mappings. It performs no ligand fitting. Missing graphs, mappings, frames or
exact records are not evaluable. [OmniDock baseline (2026/09),
`post_docking_analysis/pose_geometry.py:fixed_frame_rmsd` and
`post_docking_analysis/redocking_validation.py:run_redocking_validation`]

Protonation-equivalent or tautomer-equivalent comparison is a different metric
and must not be silently substituted for strict chemical identity. Any future
such metric needs a named/versioned definition, recorded atom mapping and
state differences, and controls that reject stereoisomers or unrelated graphs.

Docking scores are model outputs. Report their metric, direction, scoring
function, target/site and comparison cohort. A consensus is relative
prioritization among the included candidates. Its evidence should retain
failed/excluded attempts and the intended denominator; repeated seeds must not
become independent engine votes by accident. These are interpretation and input
contract requirements of this correction, not an affinity-calibration claim.

The consensus boundary expects one already-selected pose per engine/candidate.
Duplicate rows require an explicit upstream run/replicate policy. They must not
be silently treated as additional votes or collapsed through a pivot aggregate.
A missing expected engine remains in the declared denominator. This input guard
does not reconstruct seed/attempt identity already lost upstream. In particular,
the existing best-row selection can collapse multiple attempts before this
boundary. Preserving run and seed identity through the full pipeline remains
deferred; the guard is not an end-to-end replicate policy.

## Deferred scientific validation

Physical pose-quality validation remains `not_evaluated` unless a real validator
is executed. A future integration must pin the validator and settings, preserve
individual check results, and test deliberately impossible geometries and
receptor clashes. This correction does not introduce unvalidated universal
clash/strain thresholds. [Buttenschoen et al. (2024/02), methods]

Assembly modeling, missing-residue reconstruction, target-specific protonation
ensembles, metal/cofactor parameterization, covalent docking, calibrated affinity
prediction and prospective screening validation remain separate work. Existing
benchmark artifacts are historical evidence from another checkout; no benchmark
or remote execution is performed by this correction.

## Sources

[Buttenschoen et al., 2024/02] Martin Buttenschoen, Garrett M. Morris and Charlotte
M. Deane. "PoseBusters: AI-based docking methods fail to generate physically valid
poses or generalise to novel sequences." Chemical Science 15, 3130–3139.
https://doi.org/10.1039/D3SC04185A

[PDB2PQR, accessed 2026/09] PDB2PQR developers. "pdb2pqr.main." PDB2PQR 3.7.1
source documentation.
https://pdb2pqr.readthedocs.io/en/latest/_modules/pdb2pqr/main.html

[OmniDock baseline, 2026/09] OmniDock contributors. Source at
`8d4434d3a33eb83b1e12cad82944b02c83270e47`.
https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/tree/8d4434d3a33eb83b1e12cad82944b02c83270e47
