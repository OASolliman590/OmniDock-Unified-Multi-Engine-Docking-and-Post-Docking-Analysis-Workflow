#!/usr/bin/env bash
# RETIRED (Spec 036 R6).
#
# This script hard-coded PDB2PQR with protonation flips, which bypasses the
# project-wide protonation policy (Spec 036 R1) and the strict receptor module.
# It does nothing now. Use the strict module instead:
#
#   python main.py pdb prepare-protein --receptors-input <raw_dir> --receptors-output <prep_dir> \
#       --force-field AMBER --ph 7.4 --ligand-backend engine_aware_full
#
# See USAGE.md for the single supported path.

echo "autodock/prep_autodock.sh is retired (Spec 036 R6)." >&2
echo "Use: python main.py pdb prepare-protein (strict receptor module). See USAGE.md." >&2
exit 1
