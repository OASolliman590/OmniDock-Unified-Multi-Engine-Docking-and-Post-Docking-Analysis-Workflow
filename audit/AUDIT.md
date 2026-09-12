# OmniDock: scientific and engineering assessment

**Audited commit:** `cfe33056c77d817d19439a9e9affe225d00ca488` · **Date:** 12 September 2026

**Verdict:** OmniDock has enough useful infrastructure to become a strong docking platform, but its automated rankings and validation labels are not yet scientifically dependable. The immediate opportunity is to make the existing calculations and experimental records trustworthy. Expanding the number of engines or report figures before doing that would amplify existing errors.

This audit inspected protein and ligand preparation, engine execution, post-docking analysis, and development architecture. Three scientific subagents worked alongside an engineering audit. We traced active routes, checked primary scientific documentation, and executed small adversarial fixtures. This establishes concrete algorithmic and orchestration defects; it does not establish how frequently they affected previous experiments. No tracked pipeline source was changed, and no remote review, issue, commit, or job was submitted.

## Read the findings by area

| Area | Detailed report | Overall appraisal |
|---|---|---|
| Protein/ligand preparation | [Preparation audit](preparation-audit.md) | Composition, chemical identity, protonation treatment, and site geometry need stronger guarantees. |
| Docking and HPC | [Docking audit](docking-audit.md) | Useful adapters; run identity, output acceptance and local/HPC protocol parity are unreliable. |
| Post-docking science | [Post-docking audit](postdocking-audit.md) | Highest-risk area: incorrect pose RMSD, score directions, missingness, reference validation, and cross-engine comparisons. |
| Engineering and development | [Technical audit](technical-audit.md) | Preserve the existing infrastructure, but centralize scientific decisions and enforce artifact/state contracts. |
| Implementation sequence | [Revolution plan](REVOLUTION_PLAN.md) | A phased redesign with explicit acceptance criteria and a scientific benchmark. |

## Findings that should determine the next development cycle

P1 means repair before relying on affected scientific conclusions or scaling campaigns. P2 means a substantive reliability or capability issue. These priorities are audit judgments, not claims that every route or historical result is affected.

| Priority | Finding | Evidence and affected behavior | Required correction |
|---|---|---|---|
| P1 | **Pose RMSD removes binding-site displacement.** | A ligand shifted **100 Å** receives ≈ **0 Å** in the redocking and geometric RMSD routines. Both fit the ligand independently. Active redocking and default geometric consensus are affected. [PD-01](postdocking-audit.md) | Preserve the receptor coordinate frame; map chemically equivalent atoms; never fit the ligand to hide pose error. |
| P1 | **Score direction changes between modules.** | GNINA selector chooses CNNaffinity **4 over 8**; zscore normalization reverses −10/−7/−4 ordering; weighted/strict atlas chooses consensus **0.1 over 0.9**. These affect different supported paths, not every default option. [PD-02](postdocking-audit.md) | Define metric direction, units and purpose once; every selector, classifier and export must consume the same definition. |
| P1 | **Validation can certify almost no evidence.** | **1 pass + 99 missing references → HIGH / validated**. NaN and text `False` can become true benchmark flags; reference lookup can accept ordinary conformers. [PD-05](postdocking-audit.md) | Require explicit experimental reference provenance; report coverage separately; enforce target/engine-specific validation eligibility. |
| P1 | **The scored pose may not be the evaluated pose.** | PDBQT models are concatenated for redocking; selected SDF pose is not consistently respected; an ordinary second SDF record fails geometric parsing. [PD-03](postdocking-audit.md) | One exact pose identifier across score, geometry, interaction and export records. |
| P1 | **Missing or chemically mismatched engine poses can count as agreement.** | Two readable poses plus a missing third still give `all_pairs_within_cutoff`; some atom-count mismatches are truncated and accepted. [PD-04](postdocking-audit.md) | Retain expected denominators and failures; require chemical graph/atom-map consistency for same-ligand comparison. |
| P1 | **Cross-engine energy minima still decide winners.** | Winner engine, some reference deltas and exports use raw score minima/spread across engines. [PD-06](postdocking-audit.md) | Compare raw values only within compatible scoring conditions; aggregate calibrated results or within-condition ranks. |
| P1 | **A failed or changed experiment can inherit an old result.** | Failed Vina rerun retains an earlier −9 score; empty files count as resumable; failed campaign returns exit 0. [D1–D2](docking-audit.md) | Immutable attempt identity, explicit output validation and manifest-based result collection. |
| P1 | **Analysis can ignore changed docking output.** | Actual pipeline: −7.2 changed to −11.2 in raw PDBQT; ordinary rerun remains −7.2 with `cache_hit`; forced parse yields −11.2. [T1](technical-audit.md) | Hash all scientific inputs and effective settings; connect them to the cache key. |
| P1 | **Local and HPC execution can change the protocol.** | Local seed 42 and 40 Å box become no seed flag and 20 Å box in deployment. AD4 fallback files retain local paths after relocation. [D3–D4](docking-audit.md) | Compile a single resolved experiment before selecting an executor; stage all transitive assets portably. |
| P1 | **The default collected receptor retains the selected ligand.** | `common` cleaning removes solvent but not the selected ligand; original and cleaned complexes both enter preparation, and an ambiguous PDB alias can favor the original. Final retention depends on the backend. [PREP-01](preparation-audit.md) | Separate source complex, experimental reference and intended receptor composition; remove the selected instance explicitly for ordinary redocking. |
| P1 | **Preparation can lose authoritative chemistry.** | Sanitizer drops connectivity/formal-charge fields. Valid macrocycle `CG0/G0` types are rejected; separately invoked repair removes/retypes ring-closure atoms. [PREP-02–03](preparation-audit.md) | Carry a chemical graph and atom map; use engine-aware validation; regenerate unsupported inputs from chemistry. |
| P1 | **Pocket properties are presented as measurements without measuring the pocket.** | Volume is always **523.6 Å³**; an empty structural iterator gets **Good druggability, 0.508**. [PREP-06](preparation-audit.md) | Remove/rename heuristic claims; use a validated cavity method only if that output is needed. |
| P1/P2 | **QC can pass structures it did not meaningfully validate.** | NaN coordinates pass receptor QC; QC exceptions become warnings with `ok=True`; permissive receptor conversion has no residue-change enforcement. [PREP-04/09](preparation-audit.md), [D5/D10](docking-audit.md) | Explicit pass/fail/not-evaluable states, finite-coordinate checks and composition change ledgers. |
| P1/P2 | **Interaction maps can lose chemical and receptor context.** | In an explicitly requested interaction branch, PDB conversion can discard retained heterogens/chain identity and then reconstruct ligand chemistry without a template. Source evidence; ProLIF stack not run here. [PD-07](postdocking-audit.md) | Analyze the selected mapped SDF plus preserved receptor/cofactor state; distinguish failed computation from an empty fingerprint. |
| P2 | **Workflow, distribution and test contracts need hardening.** | Concurrent state updates lose one writer; failed graph prerequisites do not block all descendants; nonexistent outputs can be marked completed. Default pytest collects none, fresh Windows setup requires symlink privileges, and packaging omits standalone entrypoint modules. [T2–T6](technical-audit.md) | Transactional state, enforced stage outputs, discoverable tests/CI, clean-wheel verification and explicit platform support. |

The detailed reports include additional site-centroid/chain errors, altloc mixing, pH fallback semantics, GNINA resource flags, AD4 grid/parameter consistency, duplicate annotation joins and unsupported polypharmacology interpretations.

## Scientific interpretation

Three different questions need separate answers:

1. **Was a plausible pose sampled?** Assess chemical/physical validity and best-of-N pose recovery.
2. **Did the ranking policy select the correct pose?** Assess the actual Top1/TopN policy in a preserved receptor frame.
3. **Does the protocol prioritize experimentally active compounds?** Assess held-out enrichment or affinity ranking using appropriate experimental endpoints.

Success at one does not establish the others. RMSD alone also cannot establish physical plausibility; checks of stereochemistry, geometry and protein–ligand clashes are justified by the [PoseBusters primary study](https://pubs.rsc.org/en/content/articlehtml/2024/sc/d3sc04185a). For screening evaluation, use measured actives/inactives and account for benchmark bias; [LIT-PCBA](https://pubmed.ncbi.nlm.nih.gov/32282202/) is one primary-source-backed starting point, alongside target-specific assays.

Keep CNNscore, CNNaffinity and empirical docking energy distinct. The GNINA paper describes different pose and affinity prediction roles; three related Vina/Smina/GNINA executions are not three independent biological confirmations. This independence judgment follows from their shared methods, not from a measured correlation in this audit. [GNINA primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8191141/).

AD4 and Vina raw energies are explicitly not interchangeable according to the [Vina tutorial](https://autodock-vina.readthedocs.io/en/stable/docking_basic.html). Likewise, a top percentile within a target is an ordering statistic, not a calibrated binding probability, potency estimate, or selectivity measurement. Keep polypharmacology output as a hypothesis matrix until target-specific experimental calibration supports stronger labels.

## Engineering appraisal

**Preserve:** explicit pairlists; real engine adapters; project manifests; preparation and preflight hooks; scope selection; the artifact-graph idea; checkpoint/provenance scaffolding; modular parsers; optional report integrations.

**Change:** scientific decisions spread across several large pipeline implementations; filename-derived identity; broad success/fallback semantics; raw directory scanning; mutable shared results; inconsistent dependency declarations; route-dependent chemistry. A clean module interface should make it difficult to select a score without knowing its direction, compare poses without knowing their identity/frame, or reuse a result without knowing its exact inputs.

**Development priority:** a discoverable scientific regression suite before another broad refactor. Several existing smoke expectations bless incorrect behavior, so simply making the existing suite green is insufficient. The proposed architecture and acceptance tests are in [REVOLUTION_PLAN.md](REVOLUTION_PLAN.md).

## What was actually verified

| Verification | Outcome |
|---|---|
| Preparation probes | Seven synthetic assertion groups passed, reproducing the reported defects; no conversion backend run. |
| Post-docking probes | Reproduced RMSD, selected-pose parsing, missingness, ranking and annotation counterexamples. GNINA selection uses unchanged AST method bodies; see script for isolation details. |
| Docking probes | Reproduced failure/resume, CLI status, command/GPF/DPF generation, HPC parity and QC defects with mocked processes; no engine or scheduler execution. |
| Technical probes | Reproduced stale cache, graph contract failures and lost state updates. Integrated cache test used the real pipeline registration and parser. |
| Existing custom smoke command | First three groups passed; stopped at checkpoint ordering. Full suite did not complete. |
| Eight selected existing smoke checks | Seven passed; one blocked by Windows symlink privileges. Passing does not imply scientific validity. |
| Default pytest discovery | No tests collected after disabling unrelated host plugins. Custom `_smoke_*` functions are not registered pytest tests. |

Evidence is preserved in this directory: `preparation-probes.py`, `postdocking_probes.py`, `docking_probe.py`, `technical_probes.py`, `integrated_technical_probes.py`, their JSON results and smoke logs. Chemistry/engine environments were not installed, no real docking benchmark was run, and historical result sets were not validated. Those are the next validation activities after correcting the established defects.

**Recommended next implementation:** one correctness release covering chemical/pose identity, score semantics, receptor-frame RMSD, validation denominators and immutable runs/cache invalidation, followed by a curated target benchmark. Preserve useful orchestration while replacing the unsafe calculations behind it.
