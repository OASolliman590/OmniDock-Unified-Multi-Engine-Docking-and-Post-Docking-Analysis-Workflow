# Spec 033 software verification (2026-10-08)

**Evidence type:** software only.
**Commits:** `93b4192` (R1 layout), `a4cdf1e` (R2/R3 strict receptor preparation), `d4974cf` (option 1 PDB2PQR policy).
**Environment:** as in `specs/031-appraisal-remediation-p0/evidence/software_linux_20261008.md`, plus conda-forge PDB2PQR 3.4.1, PROPKA 3.5.1, Meeko 0.7.1 and RDKit 2026.03.6. The environment's `bin` directory must be on `PATH` so the tests find `pdb2pqr30` and `mk_prepare_receptor.py`.

| Command (at `d4974cf`) | Result |
|---|---|
| `python -m pytest test/ -q -p no:cacheprovider` | 210 passed, 1 skipped (private Spec 031 1IEP package absent) |
| `python -m pytest test/test_spec033_layout.py -q` | 6 passed (2 fail against pre-fix code; negative control) |
| `python -m pytest test/test_spec033_receptor_prep.py -q` | 27 passed, including the real 1IEP chain-A PDB2PQR → Meeko acceptance test |
| `python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix` | exit 0 |
| `python -m compileall -q docking workflow`; `git diff --check` | clean |

When PDB2PQR or Meeko is absent, the dependency-backed receptor tests skip with an explicit reason, and receptor preparation itself fails with `backend_unavailable`.
