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
    if bool(row.get("is_cocrystal_benchmark")):
        return True

    pair_source = str(row.get("pair_source") or "").strip().lower()
    if pair_source in {"cocrystal", "reference", "native", "redocking", "comparative", "compartive"}:
        return True

    site_id = str(row.get("site_id") or "").strip().lower()
    if site_id in _REFERENCE_SITE_TOKENS:
        return True

    text = " ".join(
        str(row.get(key) or "")
        for key in ("tag", "ligand", "cocrystal_ligand_name", "ligand_display_name")
    ).lower()
    return any(token in text for token in _REFERENCE_TEXT_TOKENS)


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
        lines = file_path.read_text(encoding="utf-8", errors="ignore").splitlines()
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
    a = coords_a - np.mean(coords_a, axis=0)
    b = coords_b - np.mean(coords_b, axis=0)
    covariance = np.dot(a.T, b)
    u, _, vt = np.linalg.svd(covariance)
    rot = np.dot(vt.T, u.T)
    if np.linalg.det(rot) < 0:
        vt[-1, :] *= -1
        rot = np.dot(vt.T, u.T)
    aligned = np.dot(a, rot)
    diff = aligned - b
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def _compute_pose_rmsd(docked_pose_file: Path, reference_pose_file: Path) -> Tuple[Optional[float], str]:
    dock_coords, dock_elements, dock_err = _parse_pdb_like_heavy_atoms(docked_pose_file, pose_index=1)
    if dock_coords.size == 0:
        return None, f"docked_pose_parse_failed:{dock_err or 'unknown'}"
    ref_coords, ref_elements, ref_err = _parse_pdb_like_heavy_atoms(reference_pose_file, pose_index=1)
    if ref_coords.size == 0:
        return None, f"reference_pose_parse_failed:{ref_err or 'unknown'}"
    if dock_coords.shape[0] != ref_coords.shape[0]:
        return None, f"atom_count_mismatch:{dock_coords.shape[0]}!={ref_coords.shape[0]}"
    if dock_elements != ref_elements:
        return None, "element_sequence_mismatch"
    rmsd = _kabsch_rmsd(dock_coords, ref_coords)
    if not np.isfinite(rmsd):
        return None, "rmsd_not_finite"
    return float(rmsd), ""


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
    return "\n".join(lines) + "\n"


def run_redocking_validation(
    *,
    project_dir: Path,
    best_by_engine: pd.DataFrame,
    output_dir: Path,
    pass_threshold_angstrom: float = 2.0,
    warn_threshold_angstrom: float = 3.5,
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
        else:
            rmsd, error = _compute_pose_rmsd(docked_pose_path, reference_pose_path)
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
