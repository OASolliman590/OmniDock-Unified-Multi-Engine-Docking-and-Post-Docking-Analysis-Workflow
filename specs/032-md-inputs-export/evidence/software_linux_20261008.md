# Spec 032 software evidence — Linux re-verification (2026-10-08)

**Evidence type:** software only. This record is not part of the T021 scientific acceptance packet and does not change the named-consumer status, which remains `reported_unverified`.

**Tested commit:** `46e0f6a` on branch `chore/install-agentic-bioinformatics-stack`.

**Full environment and reproduction commands:** `specs/031-appraisal-remediation-p0/evidence/software_linux_20261008.md`. In summary: Ubuntu 24.04.5, Python 3.10.22 from conda-forge, and Open Babel 3.2.1 (conda-forge CLI, PyPI 3.2.1 Python module).

## Results

| Command | Result |
|---|---|
| `python -m pytest test/test_spec032_md_inputs.py test/test_spec031_*.py -q` | 69 passed, 1 skipped. The skip is the private Spec 031 1IEP fixture; no Spec 032 test was skipped. |
| `python -m pytest test/ -q -rs` | 177 passed, 1 skipped, 0 failed |
| `python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix` | exit 0 |

The previous record (Windows, CPython 3.12 Open Babel 3.2.1 wheel) is now joined by a Linux pass. The Open Babel-backed test exercises `CorrectForPH(7.4)`, explicit hydrogens, formal-charge capture, MOL2 parse-back, strict `system.pdb` assembly and provenance publication. It **ran and passed** on Linux.

No Spec 032 source file changed. The fixes needed for a clean clone touched only test fixtures and test guards (`62e1f27`, `46e0f6a`).

## Not established

- The scientific acceptability of the protonation or charge state Open Babel predicts for any real ligand. Such states stay `human_review_required`.
- Named-consumer (CHARMM-GUI/CGenFF) acceptance remains `reported_unverified`.
- T021 is still open.
