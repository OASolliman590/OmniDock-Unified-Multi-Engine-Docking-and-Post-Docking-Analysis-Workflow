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
from post_docking_analysis.geometric_consensus import _extract_sdf_pose
from post_docking_analysis.pose_geometry import pose_rmsd, selected_record
from post_docking_analysis.score_semantics import explicit_true
from post_docking_analysis.pose_geometry import content_hash

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
    if explicit_true(row.get("is_cocrystal_benchmark")):
        return True

    pair_source = str(row.get("pair_source") or "").strip().lower()
    if pair_source in {"cocrystal", "reference", "native", "redocking", "comparative", "compartive"}:
        return True

    return False


def _atom_element_from_pdb_line(line: str) -> str:
    element = line[76:78].strip()
    if element:
        return element.upper()
    atom_name = line[12:16].strip()
    letters = "".join(ch for ch in atom_name if ch.isalpha())
    if not letters:
        return ""
    if len(letters) >= 2 and letters[0].isalpha() and letters[1].islower():
        return letters[:2].upper()
    return letters[0].upper()


def _parse_pdb_like_heavy_atoms(
    file_path: Path,
    *,
    pose_index: int = 1,
) -> Tuple[np.ndarray, List[str], str]:
    coords: List[List[float]] = []
    elements: List[str] = []

    if file_path.suffix.lower() in (".sdf", ".mol"):
        return _extract_sdf_pose(file_path, pose_index=pose_index)

    try:
        lines = selected_record(file_path, pose_index).splitlines()
    except Exception as exc:
        return np.empty((0, 3), dtype=float), [], f"read_error:{exc}"

    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if len(line) < 54:
            continue
        element = _atom_element_from_pdb_line(line)
        if not element or element.upper().startswith("H"):
            continue
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except Exception:
            continue
        coords.append([x, y, z])
        elements.append(element.upper())
    if not coords:
        return np.empty((0, 3), dtype=float), [], "no_heavy_atoms"
    return np.array(coords, dtype=float), elements, ""


def _kabsch_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    if coords_a.shape != coords_b.shape or coords_a.shape[0] <= 0:
        return float("nan")
    # Legacy API name: docking validation must remain in the receptor frame.
    diff = coords_a - coords_b
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def _positive_pose_index(value, default=None) -> int:
    if pd.isna(value):
        if default is None:
            raise ValueError("missing_pose_index")
        return default
    number = float(value)
    if not np.isfinite(number) or number < 1 or number != int(number):
        raise ValueError("invalid_pose_index")
    return int(number)


def _experimental_provenance(row: pd.Series) -> bool:
    source = row.get("reference_source")
    accession = row.get("reference_pdb_id")
    if pd.isna(source) or pd.isna(accession) or not str(accession).strip():
        return False
    label, separator, identifier = str(source).strip().lower().partition(":")
    return label in {"cocrystal", "experimental"} and (
        not separator or identifier.upper() == str(accession).strip().upper()
    )


def _compute_pose_rmsd(docked_pose_file: Path, reference_pose_file: Path, *, pose_index: int = 1, reference_pose_index: int = 1) -> Tuple[Optional[float], str]:
    try:
        value = pose_rmsd(docked_pose_file, reference_pose_file, pose_a=pose_index, pose_b=reference_pose_index)
        return (value, "") if np.isfinite(value) else (None, "rmsd_not_finite")
    except Exception as exc:
        return None, str(exc)


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


def _find_reference_pose_file(project_dir: Path, row: pd.Series) -> Optional[Path]:
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

    # Raw ligands and fuzzy file-name matches are not experimental references.
    return None


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
    if confidence == "SUFFICIENT_PROTOCOL_CHECKS":
        lines.append("The configured reference-count, coverage and pose-reproduction checks passed. This is a protocol check, not evidence of predictive affinity or physical pose validity.")
    else:
        lines.append("Reference checks are limited by coverage, independent reference count, or pose reproduction. Reference anchoring is disabled.")
    return "\n".join(lines) + "\n"


def run_redocking_validation(
    *,
    project_dir: Path,
    best_by_engine: pd.DataFrame,
    output_dir: Path,
    pass_threshold_angstrom: float = 2.0,
    warn_threshold_angstrom: float = 3.5,
    minimum_reference_complexes: int = 3,
    minimum_coverage: float = 1.0,
    expected_reference_rows: Optional[pd.DataFrame] = None,
    expected_engines: Optional[List[str]] = None,
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
        "reference_pose_file",
        "redocking_rmsd_angstrom",
        "redocking_classification",
        "validation_reason",
        "pose", "scoring_function",
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
        "engine", "scoring_function",
    ]

    best_by_engine = best_by_engine.copy() if best_by_engine is not None else pd.DataFrame()
    if expected_reference_rows is not None and not expected_reference_rows.empty:
        known = set(zip(best_by_engine.get("tag", []), best_by_engine.get("engine", [])))
        missing = []
        engines = expected_engines or sorted(set(best_by_engine.get("engine", [])))
        for _, candidate in expected_reference_rows.iterrows():
            if not _is_reference_candidate_row(candidate):
                continue
            identity = candidate.to_dict()
            identity["protein"] = identity.get("protein", identity.get("receptor", ""))
            identity["tag"] = identity.get("tag", f"{identity['protein']}_{identity.get('site_id', '')}_{identity.get('ligand', '')}")
            for engine in engines:
                if (identity["tag"], engine) not in known:
                    missing.append({**identity, "engine":engine, "pose":1, "pose_file":"", "affinity_kcal_mol":np.nan})
        if missing:
            best_by_engine = pd.concat([best_by_engine, pd.DataFrame(missing)], ignore_index=True)
    if best_by_engine.empty:
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

    frame = best_by_engine.copy()
    frame["is_reference_candidate"] = frame.apply(_is_reference_candidate_row, axis=1)
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
        reference_pose_path = _find_reference_pose_file(project_dir, row)
        rmsd: Optional[float] = None
        classification = "not_evaluable"
        reason = ""

        if docked_pose_path is None or not docked_pose_path.exists():
            reason = "missing_docked_pose_file"
        elif reference_pose_path is None or not reference_pose_path.exists():
            reason = "missing_reference_pose_file"
        elif not _experimental_provenance(row):
            reason = "missing_experimental_reference_provenance"
        elif not str(row.get("reference_pdb_id", "") or "").strip() or pd.isna(row.get("reference_pdb_id")):
            reason = "missing_reference_source_identifier"
        elif pd.isna(row.get("receptor_frame_id")) or pd.isna(row.get("reference_frame_id")) or not str(row.get("receptor_frame_id", "")).strip() or str(row.get("receptor_frame_id")) != str(row.get("reference_frame_id")):
            reason = "unverified_common_receptor_frame"
        else:
            try:
                pose_index = _positive_pose_index(row.get("pose", 1))
                reference_index = _positive_pose_index(row.get("reference_pose_index", 1), default=1)
                rmsd, error = _compute_pose_rmsd(docked_pose_path, reference_pose_path, pose_index=pose_index, reference_pose_index=reference_index)
            except (TypeError, ValueError, OverflowError):
                rmsd, error = None, "invalid_pose_index"
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
                "reference_pose_file": str(reference_pose_path) if reference_pose_path else "",
                "redocking_rmsd_angstrom": rmsd,
                "redocking_classification": classification,
                "validation_reason": reason,
                "pose": row.get("pose", 1),
                "scoring_function": row.get("scoring_function", row.get("engine", "")),
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
    if "scoring_function" not in merged:
        merged["scoring_function"] = merged["engine"]
    for (protein, engine, scoring_function), group in merged.groupby(["protein", "engine", "scoring_function"], dropna=False):
        valid = group[group["redocking_classification"] == "pass"].copy()
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
                    "engine": str(engine), "scoring_function": str(scoring_function),
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
                    "engine": str(engine), "scoring_function": str(scoring_function),
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
    coverage = float(evaluable_count / total_reference_rows) if total_reference_rows else 0.0
    # Aliases and additional engines do not create independent experimental complexes.
    unique_reference_count = len({
        (str(row.get("reference_pdb_id", "")).upper(), str(row.get("reference_frame_id", "")),
         content_hash(reference_path))
        for _, row in references.iterrows() if _experimental_provenance(row)
        if (reference_path := _find_reference_pose_file(project_dir, row)) is not None
    })
    sufficient = coverage >= minimum_coverage and unique_reference_count >= minimum_reference_complexes and pass_rate >= 0.80
    # This is protocol-check evidence, not a statistical confidence estimate.
    confidence = "SUFFICIENT_PROTOCOL_CHECKS" if sufficient else "LIMITED"

    summary_payload = {
        "total_reference_rows": total_reference_rows,
        "evaluable_rows": evaluable_count,
        "pass_rows": pass_count,
        "warn_rows": warn_count,
        "fail_rows": fail_count,
        "not_evaluable_rows": not_evaluable_count,
        "pass_rate": f"{pass_rate:.3f}",
        "confidence_signal": confidence,
        "coverage": coverage,
        "unique_reference_complexes": unique_reference_count,
        "minimum_reference_complexes": int(minimum_reference_complexes),
        "minimum_coverage": float(minimum_coverage),
    }
    gate_state = "validated" if sufficient else "needs_review"
    allow_reference_anchor = gate_state == "validated"
    gate_reason = (
        "sufficient_protocol_checks"
        if gate_state == "validated"
        else "insufficient_protocol_checks_or_coverage"
    )
    summary_payload.update(
        {
            "validation_gate_state": gate_state,
            "allow_reference_anchor": bool(allow_reference_anchor),
            "validation_gate_reason": gate_reason,
        }
    )
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
