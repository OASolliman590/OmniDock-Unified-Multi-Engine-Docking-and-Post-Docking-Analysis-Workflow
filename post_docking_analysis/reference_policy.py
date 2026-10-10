"""Shared explicit scientific-reference classification policy for Spec 031."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Dict, Mapping, Optional

import pandas as pd

from post_docking_analysis.value_normalization import normalize_boolean


REFERENCE_METHOD = "explicit_reference_metadata_v2"
LEGACY_METHOD = "legacy_inferred_reference_v1"

_IDENTITY_FIELDS = (
    "cocrystal_ligand_name",
    "reference_ligand_name",
    "reference_ligand_id",
)
_LINEAGE_FIELDS = (
    "reference_pose_file",
    "reference_ligand_file",
    "cocrystal_pose_file",
    "cocrystal_ligand_file",
    "reference_source_structure_id",
)
_RECEPTOR_LINEAGE_FIELDS = (
    "reference_receptor_file",
    "reference_receptor_id",
    "reference_source_structure_id",
    "pdb_id",
)
_LEGACY_TOKENS = ("reference", "cocrystal", "co-crystal", "native", "redocking", "benchmark", "control")


def _cell_text(value: object) -> str:
    """Text of one cell; a missing value (None or NaN after a merge) is the empty string, never 'nan'."""
    if value is None:
        return ""
    try:
        if bool(pd.isna(value)):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _text(row: Mapping[str, object], fields: tuple[str, ...]) -> Dict[str, str]:
    values = {field: _cell_text(row.get(field)) for field in fields}
    return {field: value for field, value in values.items() if value}


@dataclass(frozen=True)
class ReferenceClassification:
    classification: str
    is_reference: bool
    validation_anchor_eligible: bool
    method: str
    reason: str
    normalized_flag: bool
    identity_fields: Dict[str, str]
    lineage_fields: Dict[str, str]
    source_checksums: Dict[str, str]
    warning: str = ""
    legacy_candidate: bool = False
    legacy_heuristic: str = ""

    def to_dict(self) -> Dict[str, object]:
        return dict(self.__dict__)


def classify_reference_row(
    row: Mapping[str, object], *, legacy_mode: bool = False
) -> ReferenceClassification:
    flag = normalize_boolean(row.get("is_cocrystal_benchmark")) or normalize_boolean(row.get("is_reference"))
    identity = _text(row, _IDENTITY_FIELDS)
    lineage = _text(row, _LINEAGE_FIELDS)
    receptor_lineage = _text(row, _RECEPTOR_LINEAGE_FIELDS)
    source_checksums: Dict[str, str] = {}
    for field, value in {**lineage, **receptor_lineage}.items():
        path = Path(value).expanduser()
        if path.is_file():
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            source_checksums[field] = digest.hexdigest()
    if flag and (identity or lineage):
        source_correspondence = bool(lineage and receptor_lineage)
        return ReferenceClassification(
            classification="reference",
            is_reference=True,
            validation_anchor_eligible=source_correspondence,
            method=REFERENCE_METHOD,
            reason=("explicit_flag_identity_and_source_lineage" if source_correspondence else "explicit_flag_and_identity_or_lineage"),
            normalized_flag=True,
            identity_fields=identity,
            lineage_fields={**lineage, **receptor_lineage},
            source_checksums=source_checksums,
            warning="" if source_correspondence else "reference_source_correspondence_incomplete",
        )
    if flag:
        return ReferenceClassification(
            classification="unclassified",
            is_reference=False,
            validation_anchor_eligible=False,
            method=REFERENCE_METHOD,
            reason="explicit_flag_missing_identity_and_lineage",
            normalized_flag=True,
            identity_fields=identity,
            lineage_fields=lineage,
            source_checksums=source_checksums,
            warning="incomplete_reference_metadata",
        )
    if legacy_mode:
        for field in ("pair_source", "site_id", "tag", "ligand", "ligand_display_name"):
            value = str(row.get(field) or "").strip().lower()
            token = next((candidate for candidate in _LEGACY_TOKENS if candidate in value), "")
            if token:
                return ReferenceClassification(
                    classification="non_reference",
                    is_reference=False,
                    validation_anchor_eligible=False,
                    method=LEGACY_METHOD,
                    reason="legacy_inference_candidate_requires_confirmation",
                    normalized_flag=False,
                    identity_fields=identity,
                    lineage_fields=lineage,
                    source_checksums=source_checksums,
                    warning="legacy_candidate_cannot_establish_scientific_reference",
                    legacy_candidate=True,
                    legacy_heuristic=f"{field}:contains:{token}",
                )
    return ReferenceClassification(
        classification="non_reference",
        is_reference=False,
        validation_anchor_eligible=False,
        method=REFERENCE_METHOD,
        reason="explicit_reference_flag_false_or_missing",
        normalized_flag=False,
        identity_fields=identity,
        lineage_fields=lineage,
        source_checksums=source_checksums,
    )


def reference_mask(frame: Optional[pd.DataFrame]) -> pd.Series:
    if frame is None or frame.empty:
        return pd.Series(dtype=bool)
    return frame.apply(lambda row: classify_reference_row(row).is_reference, axis=1).astype(bool)


def classify_reference_frame(frame: pd.DataFrame, *, legacy_mode: bool = False) -> pd.DataFrame:
    output = frame.copy()
    classifications = [classify_reference_row(row, legacy_mode=legacy_mode) for _, row in output.iterrows()]
    output["reference_classification"] = [item.classification for item in classifications]
    output["reference_classification_method"] = [item.method for item in classifications]
    output["reference_classification_reason"] = [item.reason for item in classifications]
    output["reference_validation_anchor_eligible"] = [item.validation_anchor_eligible for item in classifications]
    output["reference_classification_warning"] = [item.warning for item in classifications]
    output["reference_source_checksums_json"] = [
        json.dumps(item.source_checksums, sort_keys=True) for item in classifications
    ]
    output["legacy_reference_candidate"] = [item.legacy_candidate for item in classifications]
    output["legacy_reference_heuristic"] = [item.legacy_heuristic for item in classifications]
    return output
