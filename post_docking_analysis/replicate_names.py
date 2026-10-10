"""Pair-tag and replicate naming shared by the docking runners and post-docking analysis.

Contract: ``specs/036-scientific-consistency/replicate_contract.md`` (Spec 036 R5b) and the
receptor-stem tag rule (Spec 036 R6).

- Pair tag: ``<receptor stem>_<site_id>_<ligand file name>``. The receptor stem has no
  ``.pdbqt``/``.pdb`` suffix, the ligand file name keeps its suffix. Example:
  ``1IEP_A_protein_site_1_STI_A_201.pdbqt``.
- Legacy pair tag (pre-R6): ``<receptor file name>_<site_id>_<ligand file name>``.
- Replicate job tag: ``<pair tag>__rep<NN>``, NN two-digit and 1-based.
"""
from __future__ import annotations

import re
from typing import Optional, Tuple

REPLICATE_TOKEN = "__rep"
REPLICATE_STEM_PATTERN = re.compile(r"^(?P<pair>.+)__rep(?P<rep>\d{2})$")


def receptor_stem_name(receptor: str) -> str:
    """Receptor file name without its .pdbqt or .pdb suffix (Spec 036 R6)."""
    name = str(receptor)
    for suffix in (".pdbqt", ".pdb"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return name


def pair_tag(receptor: str, site_id: str, ligand: str) -> str:
    return f"{receptor_stem_name(receptor)}_{site_id}_{ligand}"


def legacy_pair_tag(receptor: str, site_id: str, ligand: str) -> str:
    return f"{receptor}_{site_id}_{ligand}"


def replicate_job_tag(pair_tag_value: str, replicate_id: int) -> str:
    return f"{pair_tag_value}{REPLICATE_TOKEN}{int(replicate_id):02d}"


def split_replicate_stem(stem: str) -> Tuple[str, Optional[int]]:
    """Return (pair_tag, replicate_id). Legacy stems without a replicate suffix return None."""
    match = REPLICATE_STEM_PATTERN.match(str(stem))
    if not match:
        return str(stem), None
    return match.group("pair"), int(match.group("rep"))
