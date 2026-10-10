"""
Redocking validation helpers for reference-aware post-docking analysis.

This module provides a conservative validation workflow:
- detects benchmark/reference rows from pair metadata
- attempts RMSD validation when both docked and reference pose files are available
- never fakes validation when pose geometry is unavailable
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from docking.project_layout import shared_ligands_dir
from post_docking_analysis.atom_mapping import (
    IN_PLACE_METHOD,
    MAPPING_METHOD,
    LineageError,
    attach_coordinates_to_topology,
    compare_pose_pair,
    load_pdbqt_lineage_pose,
    load_sdf_graph_pose,
    not_comparable_result,
)
from post_docking_analysis.geometric_consensus import _extract_pdb_like_pose, _extract_sdf_pose
from post_docking_analysis.reference_policy import (
    classify_reference_frame,
    classify_reference_row,
)

logger = logging.getLogger(__name__)

_REFERENCE_SITE_TOKENS = {
    "reference",
    "comparative",
    "compartive",
    "native",
    "redocking",
    "control",
    "benchmark",
    "known",
    "lapi",
}
_REFERENCE_TEXT_TOKENS = (
    "reference",
    "co-crystal",
    "cocrystal",
    "native",
    "control",
    "benchmark",
    "redocking",
)


def _is_reference_candidate_row(row: pd.Series) -> bool:
    return classify_reference_row(row).is_reference


def _parse_pdb_like_heavy_atoms(
    file_path: Path,
    *,
    pose_index: int = 1,
) -> Tuple[np.ndarray, List[str], str]:
    if file_path.suffix.lower() in (".sdf", ".mol"):
        return _extract_sdf_pose(file_path, pose_index=pose_index)
    return _extract_pdb_like_pose(file_path, pose_index=pose_index)


def _compute_pose_rmsd(
    docked_pose_file: Path,
    reference_pose_file: Path,
    *,
    docked_pose_index: int = 1,
    reference_pose_index: int = 1,
    docked_topology_file: Optional[Path] = None,
    reference_topology_file: Optional[Path] = None,
) -> Tuple[Optional[float], str, Dict[str, object]]:
    """Redocking RMSD of a docked pose against the reference ligand.

    Spec 034 T002a: the returned value is the in-place heavy-atom RMSD (``in_place_v1``,
    no superposition) minimised over the valid symmetry mappings; it drives classification.
    The Kabsch-superposed value over the same mapping set is in the details as
    ``kabsch_rmsd_angstrom`` (secondary).
    """
    def _load_pose(coordinate_file: Path, pose_index: int, topology_file: Optional[Path]):
        if coordinate_file.suffix.lower() in {".sdf", ".mol"}:
            return load_sdf_graph_pose(coordinate_file, pose_index=pose_index)
        if topology_file is None or topology_file.suffix.lower() not in {".sdf", ".mol"}:
            raise ValueError("missing_explicit_sdf_topology")
        coords, elements, error = _parse_pdb_like_heavy_atoms(coordinate_file, pose_index=pose_index)
        if error or coords.size == 0:
            raise ValueError(f"pose_parse_failed:{error or 'no_heavy_atoms'}")
        return attach_coordinates_to_topology(coords, elements, topology_file)

    # Spec 034 R1d: a PDBQT pose with Meeko lineage remarks uses its lineage graph.
    # Without remarks the existing explicit-topology path is unchanged.
    try:
        lineage_docked = (
            load_pdbqt_lineage_pose(docked_pose_file, docked_pose_index)
            if docked_pose_file.suffix.lower() == ".pdbqt"
            else None
        )
    except (LineageError, TypeError, ValueError) as exc:
        reason = exc.reason if isinstance(exc, LineageError) else f"invalid_pose_index:{exc}"
        if reason != "no_lineage_remarks":
            result = not_comparable_result(f"docked_pose_not_comparable:{reason}")
            return None, result.reason, result.to_dict()
        lineage_docked = None
    if lineage_docked is not None:
        try:
            reference = _load_pose(reference_pose_file, reference_pose_index, reference_topology_file)
        except (OSError, RuntimeError, ValueError) as exc:
            result = not_comparable_result(f"reference_pose_not_comparable:{exc}")
            return None, result.reason, result.to_dict()
        # Spec 036 R1b/R2: redocking compares neutral parent graphs (charge-blind atom correspondence).
        lineage_result = compare_pose_pair(
            lineage_docked, reference, rmsd_frame=IN_PLACE_METHOD, charge_blind_parent_mapping=True
        )
        if not lineage_result.comparable:
            return None, f"not_comparable:{lineage_result.reason}", lineage_result.to_dict()
        return float(lineage_result.rmsd_angstrom), "", lineage_result.to_dict()

    try:
        docked = _load_pose(docked_pose_file, docked_pose_index, docked_topology_file)
    except (OSError, RuntimeError, ValueError) as exc:
        result = not_comparable_result(f"docked_pose_not_comparable:{exc}")
        return None, result.reason, result.to_dict()
    try:
        reference = _load_pose(reference_pose_file, reference_pose_index, reference_topology_file)
    except (OSError, RuntimeError, ValueError) as exc:
        result = not_comparable_result(f"reference_pose_not_comparable:{exc}")
        return None, result.reason, result.to_dict()
    result = compare_pose_pair(docked, reference, rmsd_frame=IN_PLACE_METHOD, charge_blind_parent_mapping=True)
    if not result.comparable:
        return None, f"not_comparable:{result.reason}", result.to_dict()
    return float(result.rmsd_angstrom), "", result.to_dict()


def _candidate_reference_roots(project_dir: Path) -> List[Path]:
    roots: List[Path] = []
    for candidate in (
        project_dir / "0-Input" / "references",
        project_dir / "0-Input" / "ligands",
        project_dir / "1-Raw_Ligand",
        project_dir / "ligands_raw",
        shared_ligands_dir(project_dir),
    ):
        if candidate.exists() and candidate.is_dir():
            roots.append(candidate)
    return roots


def _resolve_existing_path(project_dir: Path, value: object) -> Optional[Path]:
    token = str(value or "").strip()
    if not token:
        return None
    path = Path(token).expanduser()
    if path.exists():
        return path.resolve()
    if not path.is_absolute():
        candidate = (project_dir / path).resolve()
        if candidate.exists():
            return candidate
    return None


def _iter_ligand_name_tokens(row: pd.Series) -> Iterable[str]:
    seen: set[str] = set()
    for key in ("cocrystal_ligand_name", "ligand", "ligand_display_name", "reference_ligand_name"):
        token = str(row.get(key) or "").strip()
        if not token:
            continue
        stem = Path(token).stem
        if stem and stem.lower() not in seen:
            seen.add(stem.lower())
            yield stem
        if token.lower() not in seen:
            seen.add(token.lower())
            yield token


def _find_reference_pose_file(project_dir: Path, row: pd.Series, *, legacy_mode: bool = False) -> Optional[Path]:
    explicit_fields = (
        "reference_pose_file",
        "reference_ligand_file",
        "cocrystal_pose_file",
        "cocrystal_ligand_file",
        "native_pose_file",
        "native_ligand_file",
    )
    for field in explicit_fields:
        resolved = _resolve_existing_path(project_dir, row.get(field))
        if resolved is not None and resolved.is_file():
            return resolved

    if not legacy_mode:
        return None

    roots = _candidate_reference_roots(project_dir)
    if not roots:
        return None

    extensions = (".pdb", ".pdbqt")
    for token in _iter_ligand_name_tokens(row):
        for root in roots:
            exact_candidates: List[Path] = []
            fuzzy_candidates: List[Path] = []
            for ext in extensions:
                exact = root / f"{token}{ext}"
                if exact.exists() and exact.is_file():
                    exact_candidates.append(exact.resolve())
            if exact_candidates:
                return sorted(exact_candidates)[0]
            for ext in extensions:
                fuzzy_candidates.extend(sorted(root.glob(f"*{token}*{ext}")))
            if fuzzy_candidates:
                return Path(fuzzy_candidates[0]).resolve()
    return None


KABSCH_SECONDARY_METHOD = "kabsch_secondary"


def _rmsd_frame_summary(pass_threshold_angstrom: float, warn_threshold_angstrom: float) -> Dict[str, object]:
    return {
        "redocking_rmsd_frame": IN_PLACE_METHOD,
        "decision_rmsd_column": "redocking_rmsd_angstrom",
        "secondary_rmsd_method": KABSCH_SECONDARY_METHOD,
        "secondary_rmsd_column": "kabsch_secondary_rmsd_angstrom",
        "pass_threshold_angstrom": float(pass_threshold_angstrom),
        "warn_threshold_angstrom": float(warn_threshold_angstrom),
    }


def _build_validation_summary_markdown(summary: Dict[str, object]) -> str:
    lines = [
        "# Validation Summary",
        "",
        f"- total_reference_rows: {int(summary.get('total_reference_rows', 0) or 0)}",
        f"- evaluable_rows: {int(summary.get('evaluable_rows', 0) or 0)}",
        f"- pass_rows: {int(summary.get('pass_rows', 0) or 0)}",
        f"- warn_rows: {int(summary.get('warn_rows', 0) or 0)}",
        f"- fail_rows: {int(summary.get('fail_rows', 0) or 0)}",
        f"- not_evaluable_rows: {int(summary.get('not_evaluable_rows', 0) or 0)}",
        f"- pass_rate: {summary.get('pass_rate', '0.000')}",
        f"- confidence_signal: {summary.get('confidence_signal', 'LOW')}",
        "",
        "## Interpretation",
        "",
    ]
    confidence = str(summary.get("confidence_signal", "LOW")).upper()
    if confidence == "HIGH":
        lines.append("Validation confidence is high: most evaluable redocking rows passed.")
    elif confidence == "MEDIUM":
        lines.append("Validation confidence is moderate: mixed redocking pass/fail outcomes.")
    else:
        lines.append("Validation confidence is low: no/limited evaluable rows or weak redocking reproduction.")
    if summary.get("redocking_rmsd_frame"):
        lines.extend(
            [
                "",
                "## RMSD frame",
                "",
                f"- decision_rmsd: {summary.get('decision_rmsd_column', 'redocking_rmsd_angstrom')} "
                f"({summary.get('redocking_rmsd_frame')}, in-place heavy-atom RMSD, no superposition). "
                "This value drives pass/warn/fail.",
                f"- thresholds: pass <= {float(summary.get('pass_threshold_angstrom', 2.0)):.3f} A, "
                f"warn <= {float(summary.get('warn_threshold_angstrom', 3.5)):.3f} A, otherwise fail.",
                f"- secondary_rmsd: {summary.get('secondary_rmsd_column', 'kabsch_secondary_rmsd_angstrom')} "
                f"({summary.get('secondary_rmsd_method', KABSCH_SECONDARY_METHOD)}, Kabsch-superposed over the same "
                "valid mappings). Reported only; it does not affect classification.",
            ]
        )
    return "\n".join(lines) + "\n"


def run_redocking_validation(
    *,
    project_dir: Path,
    best_by_engine: pd.DataFrame,
    output_dir: Path,
    pass_threshold_angstrom: float = 2.0,
    warn_threshold_angstrom: float = 3.5,
    legacy_reference_mode: bool = False,
) -> Dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    validation_file = output_dir / "redocking_validation.csv"
    baselines_file = output_dir / "reference_baselines.csv"
    summary_file = output_dir / "VALIDATION_SUMMARY.md"
    gate_file = output_dir / "validation_gate_status.json"

    validation_columns = [
        "protein",
        "ligand",
        "tag",
        "engine",
        "docked_pose_file",
        "docked_pose_index",
        "reference_pose_file",
        "redocking_rmsd_angstrom",
        "redocking_classification",
        "validation_reason",
        "mapping_method",
        "mapping_status",
        "mapping_reason",
        "alignment_method",
        "mapped_heavy_atoms",
        "total_heavy_atoms_docked",
        "total_heavy_atoms_reference",
        "mapping_coverage",
        "valid_mapping_count",
        "mapping_backend",
        "mapping_backend_version",
        "docked_topology_sha256",
        "reference_topology_sha256",
        "lineage_source",
        "redocking_rmsd_frame",
        "kabsch_secondary_rmsd_angstrom",
        "rmsd_secondary_method",
        "charge_blind_parent_mapping",
    ]
    baseline_columns = [
        "protein",
        "reference_tag",
        "reference_ligand",
        "reference_engine",
        "reference_affinity",
        "redocking_classification",
        "anchor_status",
        "anchor_reason",
    ]

    if best_by_engine is None or best_by_engine.empty:
        empty = pd.DataFrame(
            columns=validation_columns
        )
        empty.to_csv(validation_file, index=False)
        pd.DataFrame(columns=baseline_columns).to_csv(baselines_file, index=False)
        summary_payload = {
            "total_reference_rows": 0,
            "evaluable_rows": 0,
            "pass_rows": 0,
            "warn_rows": 0,
            "fail_rows": 0,
            "not_evaluable_rows": 0,
            "pass_rate": "0.000",
            "confidence_signal": "LOW",
        }
        summary_payload["validation_gate_state"] = "needs_review"
        summary_payload["allow_reference_anchor"] = False
        summary_payload["validation_gate_reason"] = "no_reference_rows"
        summary_payload.update(_rmsd_frame_summary(pass_threshold_angstrom, warn_threshold_angstrom))
        summary_file.write_text(_build_validation_summary_markdown(summary_payload), encoding="utf-8")
        gate_file.write_text(json.dumps(summary_payload, indent=2), encoding="utf-8")
        return {
            "validation_file": str(validation_file),
            "summary_file": str(summary_file),
            "validation_gate_file": str(gate_file),
            "reference_baselines_file": str(baselines_file),
            "summary": summary_payload,
            "validation_gate": summary_payload,
            "validation_df": empty,
            "reference_baselines_df": pd.DataFrame(),
        }

    frame = classify_reference_frame(best_by_engine, legacy_mode=legacy_reference_mode)
    if legacy_reference_mode:
        frame[frame["legacy_reference_candidate"]].to_csv(
            output_dir / "legacy_reference_migration_candidates.csv", index=False
        )
    frame["is_reference_candidate"] = frame["reference_classification"].eq("reference")
    references = frame[frame["is_reference_candidate"]].copy()

    # Deduplicate: one RMSD computation per unique (protein, ligand, tag, engine) tuple.
    # best_by_engine may contain multiple score-row entries for the same pose file when
    # the upstream frame carries per-mode rows; take one representative row per pose.
    dedup_keys = [k for k in ("protein", "ligand", "tag", "engine") if k in references.columns]
    if dedup_keys:
        references = references.drop_duplicates(subset=dedup_keys, keep="first")

    records: List[Dict[str, object]] = []
    for _, row in references.iterrows():
        docked_pose_path = _resolve_existing_path(project_dir, row.get("pose_file"))
        reference_pose_path = _find_reference_pose_file(project_dir, row, legacy_mode=False)
        raw_pose_index = row.get("pose", 1)
        try:
            docked_pose_index = 1 if pd.isna(raw_pose_index) else int(raw_pose_index)
        except (TypeError, ValueError, OverflowError):
            docked_pose_index = raw_pose_index
        rmsd: Optional[float] = None
        mapping_details: Dict[str, object] = not_comparable_result("not_attempted").to_dict()
        classification = "not_evaluable"
        reason = ""

        if docked_pose_path is None or not docked_pose_path.exists():
            reason = "missing_docked_pose_file"
        elif reference_pose_path is None or not reference_pose_path.exists():
            reason = "missing_reference_pose_file"
        else:
            rmsd, error, mapping_details = _compute_pose_rmsd(
                docked_pose_path,
                reference_pose_path,
                docked_pose_index=docked_pose_index,
                docked_topology_file=_resolve_existing_path(
                    project_dir, row.get("topology_file", row.get("ligand_topology_file"))
                ),
                reference_topology_file=_resolve_existing_path(
                    project_dir, row.get("reference_ligand_file", row.get("reference_topology_file"))
                ),
            )
            if rmsd is None:
                reason = error or "rmsd_computation_failed"
            else:
                if rmsd <= float(pass_threshold_angstrom):
                    classification = "pass"
                elif rmsd <= float(warn_threshold_angstrom):
                    classification = "warn"
                else:
                    classification = "fail"

        records.append(
            {
                "protein": str(row.get("protein", "")),
                "ligand": str(row.get("ligand", "")),
                "tag": str(row.get("tag", "")),
                "engine": str(row.get("engine", "")),
                "docked_pose_file": str(docked_pose_path) if docked_pose_path else "",
                "docked_pose_index": docked_pose_index,
                "reference_pose_file": str(reference_pose_path) if reference_pose_path else "",
                "redocking_rmsd_angstrom": rmsd,
                "redocking_classification": classification,
                "validation_reason": reason,
                "mapping_method": MAPPING_METHOD,
                "mapping_status": mapping_details.get("status", "not_comparable"),
                "mapping_reason": mapping_details.get("reason", reason),
                "alignment_method": mapping_details.get("alignment_method", "kabsch"),
                "mapped_heavy_atoms": mapping_details.get("mapped_heavy_atoms", 0),
                "total_heavy_atoms_docked": mapping_details.get("total_heavy_atoms_a", 0),
                "total_heavy_atoms_reference": mapping_details.get("total_heavy_atoms_b", 0),
                "mapping_coverage": mapping_details.get("mapping_coverage", 0.0),
                "valid_mapping_count": mapping_details.get("valid_mapping_count", 0),
                "mapping_backend": mapping_details.get("backend", "networkx"),
                "mapping_backend_version": mapping_details.get("backend_version", ""),
                "docked_topology_sha256": mapping_details.get("topology_sha256_a", ""),
                "reference_topology_sha256": mapping_details.get("topology_sha256_b", ""),
                "lineage_source": mapping_details.get("lineage_source", ""),
                "redocking_rmsd_frame": IN_PLACE_METHOD if rmsd is not None else "",
                "kabsch_secondary_rmsd_angstrom": (
                    mapping_details.get("kabsch_rmsd_angstrom") if rmsd is not None else None
                ),
                "rmsd_secondary_method": (
                    KABSCH_SECONDARY_METHOD if rmsd is not None and mapping_details.get("kabsch_rmsd_angstrom") is not None else ""
                ),
                "charge_blind_parent_mapping": bool(mapping_details.get("charge_blind_parent_mapping", False)),
            }
        )

    validation_df = pd.DataFrame(records, columns=validation_columns)
    validation_df.to_csv(validation_file, index=False)

    merged = references.merge(
        validation_df[["tag", "engine", "redocking_classification"]],
        on=["tag", "engine"],
        how="left",
    )

    baseline_rows: List[Dict[str, object]] = []
    for protein, group in merged.groupby("protein", dropna=False):
        anchor_eligible = group.get("reference_validation_anchor_eligible", False)
        if not isinstance(anchor_eligible, pd.Series):
            anchor_eligible = pd.Series(bool(anchor_eligible), index=group.index)
        valid = group[(group["redocking_classification"] == "pass") & anchor_eligible.astype(bool)].copy()
        if not valid.empty:
            valid["affinity_kcal_mol"] = pd.to_numeric(valid.get("affinity_kcal_mol"), errors="coerce")
            valid = valid[np.isfinite(valid["affinity_kcal_mol"])].copy()
        if not valid.empty:
            best_ref = valid.sort_values(["affinity_kcal_mol", "tag"]).iloc[0]
            baseline_rows.append(
                {
                    "protein": str(protein),
                    "reference_tag": str(best_ref.get("tag", "")),
                    "reference_ligand": str(best_ref.get("ligand", "")),
                    "reference_engine": str(best_ref.get("engine", "")),
                    "reference_affinity": float(best_ref.get("affinity_kcal_mol")),
                    "redocking_classification": "pass",
                    "anchor_status": "eligible",
                    "anchor_reason": "validated_reference",
                }
            )
        else:
            baseline_rows.append(
                {
                    "protein": str(protein),
                    "reference_tag": "",
                    "reference_ligand": "",
                    "reference_engine": "",
                    "reference_affinity": np.nan,
                    "redocking_classification": "not_evaluable",
                    "anchor_status": "fallback_percentile",
                    "anchor_reason": "no_passed_redocking_reference",
                }
            )
    reference_baselines_df = pd.DataFrame(baseline_rows, columns=baseline_columns)
    reference_baselines_df.to_csv(baselines_file, index=False)

    evaluable_mask = validation_df["redocking_classification"].isin({"pass", "warn", "fail"})
    evaluable_count = int(evaluable_mask.sum())
    pass_count = int((validation_df["redocking_classification"] == "pass").sum())
    warn_count = int((validation_df["redocking_classification"] == "warn").sum())
    fail_count = int((validation_df["redocking_classification"] == "fail").sum())
    total_reference_rows = int(len(validation_df))
    not_evaluable_count = int((validation_df["redocking_classification"] == "not_evaluable").sum())
    pass_rate = float(pass_count / evaluable_count) if evaluable_count else 0.0
    if total_reference_rows == 0 or evaluable_count == 0:
        confidence = "LOW"
    elif pass_rate >= 0.80:
        confidence = "HIGH"
    elif pass_rate >= 0.50:
        confidence = "MEDIUM"
    else:
        confidence = "LOW"

    summary_payload = {
        "total_reference_rows": total_reference_rows,
        "evaluable_rows": evaluable_count,
        "pass_rows": pass_count,
        "warn_rows": warn_count,
        "fail_rows": fail_count,
        "not_evaluable_rows": not_evaluable_count,
        "pass_rate": f"{pass_rate:.3f}",
        "confidence_signal": confidence,
    }
    gate_state = "validated" if confidence == "HIGH" else "needs_review"
    allow_reference_anchor = gate_state == "validated"
    gate_reason = (
        "high_confidence_redocking"
        if gate_state == "validated"
        else "redocking_confidence_not_high"
    )
    summary_payload.update(
        {
            "validation_gate_state": gate_state,
            "allow_reference_anchor": bool(allow_reference_anchor),
            "validation_gate_reason": gate_reason,
        }
    )
    summary_payload.update(_rmsd_frame_summary(pass_threshold_angstrom, warn_threshold_angstrom))
    summary_file.write_text(_build_validation_summary_markdown(summary_payload), encoding="utf-8")
    gate_file.write_text(json.dumps(summary_payload, indent=2), encoding="utf-8")
    logger.info(
        "✅ Redocking validation completed: refs=%s evaluable=%s pass=%s warn=%s fail=%s",
        total_reference_rows,
        evaluable_count,
        pass_count,
        warn_count,
        fail_count,
    )
    return {
        "validation_file": str(validation_file),
        "summary_file": str(summary_file),
        "validation_gate_file": str(gate_file),
        "reference_baselines_file": str(baselines_file),
        "summary": summary_payload,
        "validation_gate": summary_payload,
        "validation_df": validation_df,
        "reference_baselines_df": reference_baselines_df,
    }


# ---------------------------------------------------------------------------
# Spec 036 R2: per-target redocking status (marks downstream outputs; never blocks).
# ---------------------------------------------------------------------------

REDOCKING_STATUS_VALIDATED = "validated"
REDOCKING_STATUS_FAILED = "failed"
REDOCKING_STATUS_NOT_EVALUATED = "not_evaluated"
REDOCKING_STATUS_NOT_COMPARABLE = "not_comparable"
REDOCKING_STATUS_COLUMNS = (
    "redocking_validation_status",
    "redocking_validation_reason",
    "redocking_reference",
)
_NOT_COMPARABLE_PREFIXES = ("not_comparable", "docked_pose_not_comparable", "reference_pose_not_comparable")
_TARGET_KEY = "protein"


def _reference_label(row: pd.Series) -> str:
    tag = str(row.get("tag", "") or "")
    engine = str(row.get("engine", "") or "")
    return f"{tag};engine={engine}" if engine else tag


def _target_status_for_group(group: pd.DataFrame, pass_threshold: float, warn_threshold: float) -> Dict[str, str]:
    if group.empty:
        return {
            "redocking_validation_status": REDOCKING_STATUS_NOT_EVALUATED,
            "redocking_validation_reason": "no_reference_rows_for_target",
            "redocking_reference": "",
        }
    classification = group["redocking_classification"].astype(str)
    rmsd = pd.to_numeric(group["redocking_rmsd_angstrom"], errors="coerce")
    evaluable = group[classification.isin({"pass", "warn", "fail"}) & rmsd.notna()].copy()
    if not evaluable.empty:
        evaluable["_rmsd"] = pd.to_numeric(evaluable["redocking_rmsd_angstrom"], errors="coerce")
        passed = evaluable[evaluable["redocking_classification"] == "pass"]
        if not passed.empty:
            best = passed.sort_values(["_rmsd", "tag", "engine"]).iloc[0]
            return {
                "redocking_validation_status": REDOCKING_STATUS_VALIDATED,
                "redocking_validation_reason": (
                    f"in_place_rmsd_{float(best['_rmsd']):.3f}A_le_{float(pass_threshold):.1f}A"
                ),
                "redocking_reference": _reference_label(best),
            }
        best = evaluable.sort_values(["_rmsd", "tag", "engine"]).iloc[0]
        band = str(best["redocking_classification"])
        if band == "warn":
            reason = (
                f"in_place_rmsd_{float(best['_rmsd']):.3f}A_in_warn_band_"
                f"{float(pass_threshold):.1f}-{float(warn_threshold):.1f}A_not_validated"
            )
        else:
            reason = f"in_place_rmsd_{float(best['_rmsd']):.3f}A_above_{float(warn_threshold):.1f}A"
        return {
            "redocking_validation_status": REDOCKING_STATUS_FAILED,
            "redocking_validation_reason": reason,
            "redocking_reference": _reference_label(best),
        }
    reasons = group["validation_reason"].fillna("").astype(str)
    comparability = reasons[reasons.str.startswith(_NOT_COMPARABLE_PREFIXES)]
    if not comparability.empty:
        reason = comparability.value_counts().index[0]
        return {
            "redocking_validation_status": REDOCKING_STATUS_NOT_COMPARABLE,
            "redocking_validation_reason": str(reason),
            "redocking_reference": "",
        }
    reason = reasons[reasons != ""]
    return {
        "redocking_validation_status": REDOCKING_STATUS_NOT_EVALUATED,
        "redocking_validation_reason": str(reason.value_counts().index[0]) if not reason.empty else "not_evaluated",
        "redocking_reference": "",
    }


def build_target_redocking_status(
    validation_df: Optional[pd.DataFrame],
    targets: Optional[Iterable[str]] = None,
    *,
    pass_threshold_angstrom: float = 2.0,
    warn_threshold_angstrom: float = 3.5,
) -> pd.DataFrame:
    """Derive one redocking status per target (receptor), from ``run_redocking_validation`` rows.

    Rules, in order:
    - ``validated``: at least one reference row passes in-place RMSD <= 2.0 A. The reference is the
      passing row with the lowest in-place RMSD.
    - ``failed``: reference rows are evaluable but none passes. The reference is the lowest-RMSD row.
      A warn-band result (2.0-3.5 A) is not validated, and is reported as ``failed`` with a warn reason.
    - ``not_comparable``: no evaluable row, and the reason is a Spec 031 mapping or comparability reason.
    - ``not_evaluated``: no reference rows for the target, or missing pose files.

    Targets with no rows are returned as ``not_evaluated`` when ``targets`` is given.
    """
    columns = [_TARGET_KEY, *REDOCKING_STATUS_COLUMNS]
    frame = validation_df.copy() if validation_df is not None else pd.DataFrame()
    if frame.empty or _TARGET_KEY not in frame.columns:
        frame = pd.DataFrame(columns=[_TARGET_KEY, "tag", "engine", "redocking_classification",
                                      "redocking_rmsd_angstrom", "validation_reason"])
    frame[_TARGET_KEY] = frame[_TARGET_KEY].astype(str)
    names = sorted(set(frame[_TARGET_KEY].tolist()) | {str(item) for item in (targets or [])})
    records: List[Dict[str, object]] = []
    for name in names:
        group = frame[frame[_TARGET_KEY] == name]
        status = _target_status_for_group(group, pass_threshold_angstrom, warn_threshold_angstrom)
        records.append({_TARGET_KEY: name, **status})
    return pd.DataFrame(records, columns=columns)


def attach_redocking_status(frame: pd.DataFrame, status_df: pd.DataFrame) -> pd.DataFrame:
    """Left-join the per-target redocking status onto a table keyed by ``protein``.

    Rows whose target has no status row get ``not_evaluated`` with ``no_reference_rows_for_target``.
    Existing columns with the same names are replaced.
    """
    working = frame.copy()
    for column in REDOCKING_STATUS_COLUMNS:
        if column in working.columns:
            working = working.drop(columns=[column])
    if working.empty:
        for column in REDOCKING_STATUS_COLUMNS:
            working[column] = pd.Series(dtype=object)
        return working
    lookup = status_df if status_df is not None and not status_df.empty else pd.DataFrame(
        columns=[_TARGET_KEY, *REDOCKING_STATUS_COLUMNS]
    )
    lookup = lookup[[_TARGET_KEY, *REDOCKING_STATUS_COLUMNS]].drop_duplicates(_TARGET_KEY)
    keys = working[_TARGET_KEY].astype(str) if _TARGET_KEY in working.columns else pd.Series("", index=working.index)
    mapping = lookup.set_index(_TARGET_KEY)
    for column in REDOCKING_STATUS_COLUMNS:
        values = keys.map(mapping[column]) if not mapping.empty else pd.Series(index=working.index, dtype=object)
        working[column] = values.where(values.notna(), "")
    missing = working["redocking_validation_status"].astype(str) == ""
    working.loc[missing, "redocking_validation_status"] = REDOCKING_STATUS_NOT_EVALUATED
    working.loc[missing, "redocking_validation_reason"] = "no_reference_rows_for_target"
    return working
