"""
Canonical visualization report suite for DAG-first post-docking analysis.

This module intentionally isolates per-figure failures so suite generation never raises.
"""
from __future__ import annotations

import json
import logging
import re
import textwrap
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from post_docking_analysis.protein_naming import format_protein_label

try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import seaborn as sns

    _PLOTTING_AVAILABLE = True
except Exception:  # pragma: no cover - dependency guard
    _PLOTTING_AVAILABLE = False
    plt = None  # type: ignore
    sns = None  # type: ignore
    mcolors = None  # type: ignore


logger = logging.getLogger(__name__)

_SUITE_STYLE: Dict[str, object] = {
    "dpi": 300,
    "figsize": (12, 7),
    "small_figsize": (10, 6),
    "large_figsize": (15, 9),
    "font_scale": 1.0,
    "context": "paper",
    "style": "whitegrid",
    "palette": "deep",
}

# DockForge co-crystal ligand naming convention: {PDB_CODE}_ligand_{resname}_{chain}_{seqnum}
_COCRYSTAL_LIGAND_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}_ligand_", re.IGNORECASE)
# MD snapshot suffix pattern, e.g. "3LN1_receptor_600ps.pdbqt"
_MD_SNAPSHOT_RE = re.compile(r"_(\d+ps)\b", re.IGNORECASE)

_HIT_CLASS_COLORS = {
    "Strong": "#2e7d32",
    "Moderate": "#9ccc65",
    "Weak": "#f9a825",
    "Inactive": "#ef5350",
    "Unclassified": "#bdbdbd",
    "Uncategorized": "#bdbdbd",
}


@dataclass
class FigureResult:
    group: str
    name: str
    path: str
    title: str
    conditional: bool = False
    generated: bool = False
    reason: str = ""
    error: str = ""


@dataclass
class PlotContext:
    output_root: Path
    classified_hits: pd.DataFrame
    consensus_ranked: pd.DataFrame
    engine_agreement: pd.DataFrame
    normalized_scores: pd.DataFrame
    top_pose_global: pd.DataFrame
    engine_scope_config: Dict[str, object]
    validation_gate: Dict[str, object]
    validation_table: pd.DataFrame
    reference_baselines: pd.DataFrame
    is_multi_engine: bool
    engine_name: Optional[str]
    primary_score_col: str
    primary_score_label: str
    has_references: bool
    validation_passed: bool
    validation_status: str
    engine_count: int
    protein_name_map: Optional[Dict[str, Dict[str, str]]] = field(default=None)
    engine_affinity_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    style: Dict[str, object] = field(default_factory=lambda: dict(_SUITE_STYLE))


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _safe_read_csv(path: Path) -> pd.DataFrame:
    try:
        if path.exists() and path.is_file():
            return pd.read_csv(path)
    except Exception as exc:
        logger.warning("Visualization suite: failed to read %s (%s)", path, exc)
    return pd.DataFrame()


def _safe_read_json(path: Path) -> Dict[str, object]:
    try:
        if path.exists() and path.is_file():
            payload = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(payload, dict):
                return payload
    except Exception as exc:
        logger.warning("Visualization suite: failed to read %s (%s)", path, exc)
    return {}


def _ensure_columns(frame: pd.DataFrame, columns: Sequence[str], default: object = np.nan) -> pd.DataFrame:
    working = frame.copy()
    for column in columns:
        if column not in working.columns:
            working[column] = default
    return working


def _best_rows(frame: pd.DataFrame, group_cols: Sequence[str], score_col: str) -> pd.DataFrame:
    if frame is None or frame.empty:
        return pd.DataFrame(columns=list(frame.columns) if frame is not None else [])
    working = frame.copy()
    working[score_col] = pd.to_numeric(working.get(score_col), errors="coerce")
    working = working[np.isfinite(working[score_col])].copy()
    if working.empty:
        return working
    ordered = working.sort_values(list(group_cols) + [score_col], ascending=[True] * len(group_cols) + [True])
    return ordered.groupby(list(group_cols), dropna=False, sort=False).head(1).copy()


def _resolve_engine_mode(scope_config: Mapping[str, object], normalized_scores: pd.DataFrame) -> Tuple[bool, Optional[str], str, str]:
    engines = [str(item).strip().lower() for item in (scope_config.get("engines_in_scope") or []) if str(item).strip()]
    is_multi = len(engines) >= 2
    engine_name = None if is_multi else (engines[0] if engines else "")
    primary_col = "affinity_kcal_mol"
    primary_label = "vina_affinity"
    if not is_multi and engine_name == "gnina":
        cnn_aff = pd.to_numeric(normalized_scores.get("cnn_affinity"), errors="coerce")
        if cnn_aff.notna().any():
            primary_col = "cnn_affinity"
            primary_label = "cnn_affinity"
    elif not is_multi and engine_name in {"vina", "smina"}:
        primary_col = "affinity_kcal_mol"
        primary_label = "vina_affinity"
    return is_multi, (engine_name or None), primary_col, primary_label


def _detect_references(classified_hits: pd.DataFrame) -> bool:
    if classified_hits is None or classified_hits.empty:
        return False
    frame = classified_hits.copy()
    if "ligand_type" in frame.columns:
        series = frame["ligand_type"].astype(str).str.lower()
        if series.eq("reference").any():
            return True
    if "is_cocrystal_benchmark" in frame.columns:
        series = frame["is_cocrystal_benchmark"].astype(str).str.lower()
        if series.isin({"true", "1", "yes"}).any():
            return True
    if "tag" in frame.columns:
        if frame["tag"].astype(str).str.contains("reference_", case=False, na=False).any():
            return True
    # DockForge co-crystal naming convention: {PDB_CODE}_ligand_{resname}_{chain}_{seqnum}
    if "ligand" in frame.columns:
        if frame["ligand"].astype(str).str.contains(_COCRYSTAL_LIGAND_RE).any():
            return True
    return False


def _validation_state(validation_gate: Mapping[str, object]) -> Tuple[str, bool]:
    state = str(validation_gate.get("status") or validation_gate.get("validation_gate_state") or "").strip().lower()
    passed = state in {"validated", "pass", "high_confidence_redocking"}
    return state or "unknown", passed


def _apply_suite_style() -> None:
    if not _PLOTTING_AVAILABLE:
        return
    sns.set_theme(
        context=str(_SUITE_STYLE["context"]),
        style=str(_SUITE_STYLE["style"]),
        palette=str(_SUITE_STYLE["palette"]),
        font_scale=float(_SUITE_STYLE["font_scale"]),
    )


def _save_placeholder(path: Path, title: str, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not _PLOTTING_AVAILABLE:
        path.write_text(message + "\n", encoding="utf-8")
        return
    fig, ax = plt.subplots(figsize=_SUITE_STYLE["small_figsize"])
    ax.axis("off")
    ax.set_title(title)
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _as_relative(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except Exception:
        return str(path.resolve())


def _pname(ctx: "PlotContext", raw: str) -> str:
    """Return a human-readable protein display label using the protein_name_map."""
    mapping = getattr(ctx, "protein_name_map", None)
    text = str(raw).strip()
    clean = text[:-6] if text.lower().endswith(".pdbqt") else text
    snap = _MD_SNAPSHOT_RE.search(clean)
    display = format_protein_label(clean, mapping)

    # Append MD snapshot suffix when present and not already in display name
    if snap and snap.group(1).lower() not in display.lower():
        display = f"{display} {snap.group(1)}"

    return display


def _is_reference_ligand(value: object) -> bool:
    text = str(value or "").strip()
    if not text:
        return False
    return bool(_COCRYSTAL_LIGAND_RE.search(text) or re.search(r"\bref(?:erence)?\b", text, flags=re.IGNORECASE))


def _ligand_label(value: object, *, reference: Optional[bool] = None) -> str:
    text = str(value or "").strip()
    if not text:
        text = "Unknown Ligand"
    label = text[:-6] if text.lower().endswith(".pdbqt") else text
    if reference is None:
        reference = _is_reference_ligand(label)
    if reference and "[Ref]" not in label:
        label = f"{label} [Ref]"
    return label


def _ligand_reference_mask(frame: pd.DataFrame) -> pd.Series:
    if frame is None or frame.empty:
        return pd.Series(dtype=bool)
    mask = pd.Series(False, index=frame.index)
    if "ligand_type" in frame.columns:
        mask |= frame["ligand_type"].astype(str).str.lower().eq("reference")
    if "is_cocrystal_benchmark" in frame.columns:
        raw = frame["is_cocrystal_benchmark"]
        mask |= raw.astype(str).str.lower().isin({"true", "1", "yes"})
    if "tag" in frame.columns:
        mask |= frame["tag"].astype(str).str.contains("reference_", case=False, na=False)
    if "ligand" in frame.columns:
        mask |= frame["ligand"].astype(str).apply(_is_reference_ligand)
    return mask


def _group_ligand_columns(columns: Sequence[object]) -> Tuple[List[str], List[str]]:
    labels = [str(col) for col in columns]
    references = [col for col in labels if _is_reference_ligand(col)]
    series = [col for col in labels if col not in references]
    return references, series


def _apply_category_tick_layout(ax, *, axis: str = "x", rotation: int = 55, wrap_width: int = 20, fontsize: int = 8) -> None:
    ticklabels = ax.get_xticklabels() if axis == "x" else ax.get_yticklabels()
    wrapped = [textwrap.fill(str(t.get_text()), width=wrap_width) for t in ticklabels]
    if axis == "x":
        ax.set_xticklabels(wrapped, rotation=rotation, ha="right", fontsize=fontsize)
    else:
        ax.set_yticklabels(wrapped, fontsize=fontsize)


def _figure_width_for_categories(n_items: int, *, minimum: float = 12.0, per_item: float = 0.42) -> float:
    return max(minimum, per_item * max(n_items, 1))


def _build_engine_rank_matrix(ctx: PlotContext) -> pd.DataFrame:
    best = _best_scored_frame(ctx, include_engine=True)
    if best.empty:
        return pd.DataFrame()
    matrix = best.pivot_table(index="tag", columns="engine", values="affinity_kcal_mol", aggfunc="min")
    if matrix.shape[1] < 2:
        return pd.DataFrame()
    matrix = matrix.dropna(axis=0, how="any")
    if len(matrix) < 3:
        return pd.DataFrame()
    return matrix


def _build_engine_agreement_summary(ctx: PlotContext) -> pd.DataFrame:
    denominator = float(max(ctx.engine_count, 1))
    source = pd.DataFrame()
    if not ctx.consensus_ranked.empty and "tag" in ctx.consensus_ranked.columns:
        source = ctx.consensus_ranked.copy()
        source["agreement_fraction"] = pd.to_numeric(source.get("agreement_fraction"), errors="coerce")
        if "agreement_count" in source.columns and source["agreement_fraction"].isna().all() and denominator > 1:
            source["agreement_fraction"] = pd.to_numeric(source.get("agreement_count"), errors="coerce") / denominator
        source["ligand"] = source.get("ligand", source.get("tag", ""))
        source["tag"] = source["tag"].astype(str)
    if source.empty or source["agreement_fraction"].isna().all():
        agreement = _ensure_columns(ctx.engine_agreement, ["tag", "engine_support_count"], default=np.nan)
        if agreement.empty:
            return pd.DataFrame()
        agreement["engine_support_count"] = pd.to_numeric(agreement.get("engine_support_count"), errors="coerce")
        agreement["agreement_fraction"] = agreement["engine_support_count"] / denominator
        mapping = pd.DataFrame()
        for frame in (ctx.classified_hits, ctx.normalized_scores):
            if not frame.empty and "tag" in frame.columns and "ligand" in frame.columns:
                mapping = frame[["tag", "ligand"]].drop_duplicates("tag")
                break
        source = agreement.merge(mapping, on="tag", how="left")
        if "ligand" not in source.columns:
            source["ligand"] = source["tag"]
    source["ligand"] = source["ligand"].astype(str)
    source["display_ligand"] = source["ligand"].map(_ligand_label)
    source["agreement_fraction"] = pd.to_numeric(source.get("agreement_fraction"), errors="coerce")
    source = source[source["agreement_fraction"].notna()].copy()
    if source.empty:
        return source
    source["full_support"] = source["agreement_fraction"] >= 0.999
    return (
        source.groupby(["ligand", "display_ligand"], dropna=False)
        .agg(
            matched_complexes=("tag", "nunique"),
            full_support_fraction=("full_support", "mean"),
            mean_agreement_fraction=("agreement_fraction", "mean"),
        )
        .reset_index()
        .sort_values(["full_support_fraction", "mean_agreement_fraction", "ligand"], ascending=[False, False, True])
    )


def _pose_diversity_diagnostic(values: pd.Series) -> Optional[str]:
    if values.empty:
        return "No rmsd_lb values available."
    finite = pd.to_numeric(values, errors="coerce").dropna()
    if finite.empty:
        return "No finite rmsd_lb values available."
    if (finite == 0.0).all():
        return "All rmsd_lb values are 0.0; this suggests pose-collapse or missing diversity metadata rather than informative spread."
    if finite.nunique() <= 1:
        return "rmsd_lb is constant across all Vina poses; pose diversity is not informative for this run."
    return None


def _normalize_class(value: object) -> str:
    token = str(value or "").strip()
    if not token:
        return "Unclassified"
    normalized = token.capitalize()
    if normalized == "Strong":
        return "Strong"
    if normalized == "Moderate":
        return "Moderate"
    if normalized == "Weak":
        return "Weak"
    if normalized == "Inactive":
        return "Inactive"
    return "Unclassified"


def _load_protein_name_map(mapping_file: Optional[Path]) -> Optional[Dict[str, Dict[str, str]]]:
    """Load the protein name mapping CSV produced by protein_naming.py into a mapping dict."""
    if not mapping_file or not mapping_file.exists():
        return None
    try:
        df = pd.read_csv(mapping_file)
        pdb_to_name: Dict[str, str] = {}
        rec_to_name: Dict[str, str] = {}
        for _, row in df.iterrows():
            rt = str(row.get("record_type", "")).strip()
            key = str(row.get("key", "")).strip()
            display = str(row.get("display_name", "")).strip()
            if not key or not display:
                continue
            if rt == "pdb_code":
                pdb_to_name[key.upper()] = display
            else:
                rec_to_name[key] = display
        return {"pdb_code_to_name": pdb_to_name, "receptor_to_name": rec_to_name}
    except Exception as exc:
        logger.warning("Visualization suite: failed to load protein name mapping %s (%s)", mapping_file, exc)
        return None


def _build_plot_context(
    *,
    classified_hits_file: Path,
    consensus_ranked_file: Path,
    engine_agreement_file: Path,
    normalized_scores_file: Path,
    engine_scope_config_file: Path,
    validation_gate_file: Path,
    top_pose_global_file: Path,
    output_dir: Path,
    protein_name_mapping_file: Optional[Path] = None,
) -> PlotContext:
    classified_hits = _safe_read_csv(classified_hits_file)
    consensus_ranked = _safe_read_csv(consensus_ranked_file)
    engine_agreement = _safe_read_csv(engine_agreement_file)
    normalized_scores = _safe_read_csv(normalized_scores_file)
    top_pose_global = _safe_read_csv(top_pose_global_file)
    scope_config = _safe_read_json(engine_scope_config_file)
    validation_gate = _safe_read_json(validation_gate_file)
    _rdv_primary = validation_gate_file.parent / "redocking_validation.csv"
    _rdv_raw = validation_gate_file.parent.parent / "raw" / "redocking_validation.csv"
    validation_table = _safe_read_csv(_rdv_primary)
    # Fallback: the consensus-level copy may be empty; use the raw copy if so.
    if (validation_table is None or validation_table.empty) and _rdv_raw.exists():
        validation_table = _safe_read_csv(_rdv_raw)
    reference_baselines = _safe_read_csv(validation_gate_file.parent / "reference_baselines.csv")
    engine_affinity_matrix = _safe_read_csv(validation_gate_file.parent / "engine_affinity_matrix.csv")
    protein_name_map = _load_protein_name_map(protein_name_mapping_file)

    is_multi, engine_name, primary_col, primary_label = _resolve_engine_mode(scope_config, normalized_scores)
    has_references = _detect_references(classified_hits)
    validation_status, validation_passed = _validation_state(validation_gate)
    engine_count = int(len([e for e in (scope_config.get("engines_in_scope") or []) if str(e).strip()]))

    return PlotContext(
        output_root=output_dir.resolve(),
        classified_hits=classified_hits,
        consensus_ranked=consensus_ranked,
        engine_agreement=engine_agreement,
        normalized_scores=normalized_scores,
        top_pose_global=top_pose_global,
        engine_scope_config=scope_config,
        validation_gate=validation_gate,
        validation_table=validation_table,
        reference_baselines=reference_baselines,
        is_multi_engine=is_multi,
        engine_name=engine_name,
        primary_score_col=primary_col,
        primary_score_label=primary_label,
        has_references=has_references,
        validation_passed=validation_passed,
        validation_status=validation_status,
        engine_count=engine_count,
        protein_name_map=protein_name_map,
        engine_affinity_matrix=engine_affinity_matrix,
    )


def _run_plot_function(
    ctx: PlotContext,
    *,
    group: str,
    name: str,
    title: str,
    fn: Callable[[PlotContext, Path], None],
    figure_overrides: Optional[Mapping[str, Callable[[PlotContext, Path], None]]] = None,
    conditional: bool = False,
    skip_reason: str = "",
) -> FigureResult:
    path = ctx.output_root / group / name
    if skip_reason:
        return FigureResult(
            group=group,
            name=name,
            path=str(path),
            title=title,
            conditional=conditional,
            generated=False,
            reason=skip_reason,
        )
    runner = figure_overrides.get(name) if figure_overrides else None
    if runner is None:
        runner = fn
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        runner(ctx, path)
        generated = bool(path.exists())
        return FigureResult(
            group=group,
            name=name,
            path=str(path),
            title=title,
            conditional=conditional,
            generated=generated,
            reason="" if generated else "figure_not_written",
        )
    except Exception as exc:  # no-raise contract
        logger.warning("Visualization suite figure failed: %s (%s)", name, exc)
        return FigureResult(
            group=group,
            name=name,
            path=str(path),
            title=title,
            conditional=conditional,
            generated=False,
            reason="figure_exception",
            error=str(exc),
        )


def _best_scored_frame(ctx: PlotContext, *, include_engine: bool = True) -> pd.DataFrame:
    scores = _ensure_columns(
        ctx.normalized_scores,
        ["engine", "tag", "protein", "ligand", "affinity_kcal_mol", "cnn_affinity", "cnn_score", "rmsd_lb", "rmsd_ub"],
        default=np.nan,
    )
    if scores.empty:
        return scores
    score_col = ctx.primary_score_col if ctx.primary_score_col in scores.columns else "affinity_kcal_mol"
    grouped = ["engine", "tag"] if include_engine else ["tag"]
    best = _best_rows(scores, grouped, score_col)
    if best.empty:
        best = _best_rows(scores, ["engine", "tag"], "affinity_kcal_mol")
    best["affinity_kcal_mol"] = pd.to_numeric(best.get("affinity_kcal_mol"), errors="coerce")
    return best


def _build_affinity_matrix(ctx: PlotContext) -> Tuple[pd.DataFrame, List[str], List[str]]:
    frame = _ensure_columns(ctx.classified_hits, ["protein", "ligand", "best_affinity_kcal_mol", "affinity_kcal_mol"], default=np.nan)
    if frame.empty:
        return pd.DataFrame(), [], []
    frame["protein"] = frame["protein"].astype(str)
    frame["ligand"] = frame["ligand"].astype(str)
    frame["best_affinity_kcal_mol"] = pd.to_numeric(frame.get("best_affinity_kcal_mol"), errors="coerce")
    frame["affinity_kcal_mol"] = pd.to_numeric(frame.get("affinity_kcal_mol"), errors="coerce")
    frame["selected_affinity"] = frame["best_affinity_kcal_mol"].where(frame["best_affinity_kcal_mol"].notna(), frame["affinity_kcal_mol"])
    best = _best_rows(frame, ["protein", "ligand"], "selected_affinity")
    matrix = best.pivot(index="protein", columns="ligand", values="selected_affinity")

    references: List[str] = []
    if ctx.has_references:
        ref_mask = pd.Series(False, index=best.index)
        if "ligand_type" in best.columns:
            ref_mask |= best["ligand_type"].astype(str).str.lower().eq("reference")
        if "is_cocrystal_benchmark" in best.columns:
            ref_mask |= best["is_cocrystal_benchmark"].astype(str).str.lower().isin({"true", "1", "yes"})
        if "tag" in best.columns:
            ref_mask |= best["tag"].astype(str).str.contains("reference_", case=False, na=False)
        if "ligand" in best.columns:
            ref_mask |= best["ligand"].astype(str).str.contains(_COCRYSTAL_LIGAND_RE)
        references = sorted(best.loc[ref_mask, "ligand"].dropna().astype(str).unique().tolist())
    all_ligands = [str(col) for col in matrix.columns]
    series = [lig for lig in all_ligands if lig not in references]
    return matrix, references, series


def _plot_score_distribution(ctx: PlotContext, output: Path) -> None:
    best = _best_scored_frame(ctx, include_engine=True)
    if best.empty:
        _save_placeholder(output, "Score Distribution", "No score rows available.")
        return
    score_col = ctx.primary_score_col if ctx.primary_score_col in best.columns else "affinity_kcal_mol"
    best[score_col] = pd.to_numeric(best.get(score_col), errors="coerce")
    best = best[np.isfinite(best[score_col])].copy()
    if best.empty:
        _save_placeholder(output, "Score Distribution", "No finite scores available.")
        return

    if ctx.is_multi_engine:
        engines = sorted(best["engine"].astype(str).dropna().unique().tolist())
        n = max(1, len(engines))
        fig, axes = plt.subplots(1, n, figsize=(6 * n, 5), squeeze=False)
        for idx, engine in enumerate(engines):
            ax = axes[0, idx]
            sub = best[best["engine"].astype(str) == engine]
            favorable = sub[sub[score_col] <= 0][score_col]
            unfavorable = sub[sub[score_col] > 0][score_col]
            if not favorable.empty:
                ax.hist(favorable, bins=20, color="#1f77b4", alpha=0.75, label="Favorable (<=0)")
            if not unfavorable.empty:
                ax.hist(unfavorable, bins=20, color="#d62728", alpha=0.75, label="Unfavorable (>0)")
            mean_val = float(sub[score_col].mean())
            median_val = float(sub[score_col].median())
            ax.axvline(mean_val, color="#d62728", linestyle="--", linewidth=1.5, label="Mean")
            ax.axvline(median_val, color="#2ca02c", linestyle="--", linewidth=1.5, label="Median")
            ax.axvline(0.0, color="#ff7f0e", linestyle=":", linewidth=1.5, label="Zero")
            un_count = int((sub[score_col] > 0).sum())
            total = int(len(sub))
            pct = (100.0 * un_count / total) if total else 0.0
            ax.text(
                0.98,
                0.95,
                f"Unfavorable: {un_count}/{total} ({pct:.1f}%)",
                ha="right",
                va="top",
                transform=ax.transAxes,
                fontsize=9,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7),
            )
            ax.set_title(f"{engine.upper()} ({ctx.primary_score_label})")
            ax.set_xlabel("Score (kcal/mol)")
            ax.set_ylabel("Count")
            ax.legend(loc="best", fontsize=8)
    else:
        fig, ax = plt.subplots(figsize=_SUITE_STYLE["figsize"])
        favorable = best[best[score_col] <= 0][score_col]
        unfavorable = best[best[score_col] > 0][score_col]
        if not favorable.empty:
            ax.hist(favorable, bins=20, color="#1f77b4", alpha=0.75, label="Favorable (<=0)")
        if not unfavorable.empty:
            ax.hist(unfavorable, bins=20, color="#d62728", alpha=0.75, label="Unfavorable (>0)")
        mean_val = float(best[score_col].mean())
        median_val = float(best[score_col].median())
        ax.axvline(mean_val, color="#d62728", linestyle="--", linewidth=1.5, label="Mean")
        ax.axvline(median_val, color="#2ca02c", linestyle="--", linewidth=1.5, label="Median")
        ax.axvline(0.0, color="#ff7f0e", linestyle=":", linewidth=1.5, label="Zero")
        un_count = int((best[score_col] > 0).sum())
        total = int(len(best))
        pct = (100.0 * un_count / total) if total else 0.0
        ax.text(
            0.98,
            0.95,
            f"Unfavorable: {un_count}/{total} ({pct:.1f}%)",
            ha="right",
            va="top",
            transform=ax.transAxes,
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7),
        )
        label = (ctx.engine_name or "single_engine").upper()
        ax.set_title(f"{label} Score Distribution ({ctx.primary_score_label})")
        ax.set_xlabel("Score (kcal/mol)")
        ax.set_ylabel("Count")
        ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _polypharmacology_table(ctx: PlotContext) -> pd.DataFrame:
    frame = _ensure_columns(ctx.classified_hits, ["ligand", "protein", "docking_quality_class", "best_affinity_kcal_mol", "affinity_kcal_mol"], default=np.nan)
    if frame.empty:
        return frame
    frame["ligand"] = frame["ligand"].astype(str)
    frame["protein"] = frame["protein"].astype(str)
    frame["docking_quality_class"] = frame["docking_quality_class"].map(_normalize_class)
    frame["best_affinity_kcal_mol"] = pd.to_numeric(frame.get("best_affinity_kcal_mol"), errors="coerce")
    frame["affinity_kcal_mol"] = pd.to_numeric(frame.get("affinity_kcal_mol"), errors="coerce")
    frame["selected_affinity"] = frame["best_affinity_kcal_mol"].where(frame["best_affinity_kcal_mol"].notna(), frame["affinity_kcal_mol"])
    best = _best_rows(frame, ["ligand", "protein"], "selected_affinity")
    best["moderate_or_better"] = best["docking_quality_class"].isin({"Strong", "Moderate"})
    return best


def _plot_consensus_leaderboard(ctx: PlotContext, output: Path) -> None:
    table = _polypharmacology_table(ctx)
    if table.empty:
        _save_placeholder(output, "Consensus Leaderboard", "No classified hits available.")
        return
    table["reference_ligand"] = _ligand_reference_mask(table)
    grouped = (
        table.groupby("ligand", dropna=False)
        .agg(
            targets=("protein", "nunique"),
            poly_score=("moderate_or_better", "sum"),
            mean_affinity=("selected_affinity", "mean"),
            reference_ligand=("reference_ligand", "max"),
        )
        .reset_index()
    )
    grouped.sort_values(["poly_score", "mean_affinity", "ligand"], ascending=[False, True, True], inplace=True)
    grouped["display_ligand"] = grouped.apply(
        lambda row: _ligand_label(row.get("ligand", ""), reference=bool(row.get("reference_ligand", False))),
        axis=1,
    )

    better_ref_map: Dict[str, int] = {}
    if ctx.has_references:
        ref_rows = table.copy()
        ref_mask = pd.Series(False, index=ref_rows.index)
        if "ligand_type" in ref_rows.columns:
            ref_mask |= ref_rows["ligand_type"].astype(str).str.lower().eq("reference")
        if "tag" in ref_rows.columns:
            ref_mask |= ref_rows["tag"].astype(str).str.contains("reference_", case=False, na=False)
        if "ligand" in ref_rows.columns:
            ref_mask |= ref_rows["ligand"].astype(str).str.contains(_COCRYSTAL_LIGAND_RE)
        refs = (
            ref_rows.loc[ref_mask]
            .groupby("protein", dropna=False)["selected_affinity"]
            .min()
            .rename("_ref_affinity")
            .reset_index()
        )
        # Drop any pre-existing reference_affinity column before merging to avoid _x/_y collision
        merge_table = table.drop(columns=["reference_affinity"], errors="ignore")
        merged = merge_table.merge(refs, on="protein", how="left")
        merged["better_than_reference"] = merged["selected_affinity"] < merged["_ref_affinity"]
        better_ref = merged.groupby("ligand", dropna=False)["better_than_reference"].sum().astype(int)
        better_ref_map = better_ref.to_dict()

    fig, ax = plt.subplots(figsize=(max(12, _figure_width_for_categories(len(grouped), minimum=12.0, per_item=0.28)), max(6, 0.42 * len(grouped))))
    y = np.arange(len(grouped))
    cmap = plt.cm.Greens
    norm = mcolors.Normalize(vmin=float(grouped["poly_score"].min()), vmax=float(grouped["poly_score"].max() or 1.0))
    colors = [cmap(norm(v)) for v in grouped["poly_score"]]
    ax.barh(y, grouped["poly_score"], color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels(grouped["display_ligand"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Polypharmacology score (Moderate-or-better targets)")
    ax.set_title("Consensus Leaderboard")
    for idx, row in grouped.reset_index(drop=True).iterrows():
        extra = f", better_ref={better_ref_map.get(str(row['ligand']), 0)}" if ctx.has_references else ""
        ax.text(
            float(row["poly_score"]) + 0.05,
            idx,
            f"targets={int(row['targets'])}{extra}",
            va="center",
            fontsize=8,
        )
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_affinity_per_ligand(ctx: PlotContext, output: Path) -> None:
    table = _polypharmacology_table(ctx)
    if table.empty:
        _save_placeholder(output, "Affinity Per Ligand", "No ligand affinities available.")
        return
    table["reference_ligand"] = _ligand_reference_mask(table)
    score_rank = {"Strong": 0, "Moderate": 1, "Weak": 2, "Inactive": 3, "Unclassified": 4, "Uncategorized": 4}
    grouped = (
        table.groupby("ligand", dropna=False)
        .agg(
            best_affinity=("selected_affinity", "min"),
            hit_class=("docking_quality_class", lambda s: sorted([_normalize_class(v) for v in s], key=lambda k: score_rank.get(k, 9))[0]),
            reference_ligand=("reference_ligand", "max"),
        )
        .reset_index()
        .sort_values(["best_affinity", "ligand"], ascending=[True, True])
    )
    grouped["display_ligand"] = grouped.apply(
        lambda row: _ligand_label(row.get("ligand", ""), reference=bool(row.get("reference_ligand", False))),
        axis=1,
    )
    fig, ax = plt.subplots(figsize=(_figure_width_for_categories(len(grouped), minimum=15.0, per_item=0.48), _SUITE_STYLE["large_figsize"][1]))
    colors = [_HIT_CLASS_COLORS.get(str(v), "#bdbdbd") for v in grouped["hit_class"]]
    ax.bar(grouped["display_ligand"], grouped["best_affinity"], color=colors)
    ax.axhline(-7.0, color="#d62728", linestyle="--", linewidth=1.2, label="-7 kcal/mol")
    ax.set_ylabel("Best affinity (kcal/mol)")
    ax.set_xlabel("Ligand")
    ax.set_title("Affinity Per Ligand")
    _apply_category_tick_layout(ax, axis="x", rotation=70, wrap_width=18, fontsize=8)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_full_affinity_heatmap(ctx: PlotContext, output: Path) -> None:
    matrix, references, series = _build_affinity_matrix(ctx)
    if matrix.empty:
        _save_placeholder(output, "Full Affinity Heatmap", "No affinity matrix available.")
        return
    if references:
        diagonal_map: Dict[str, str] = {}
        ref_rows = ctx.classified_hits.copy()
        ref_mask_h = pd.Series(False, index=ref_rows.index)
        if "ligand_type" in ref_rows.columns:
            ref_mask_h |= ref_rows["ligand_type"].astype(str).str.lower().eq("reference")
        if "tag" in ref_rows.columns:
            ref_mask_h |= ref_rows["tag"].astype(str).str.contains("reference_", case=False, na=False)
        if "ligand" in ref_rows.columns:
            ref_mask_h |= ref_rows["ligand"].astype(str).str.contains(_COCRYSTAL_LIGAND_RE)
        ref_rows = ref_rows[ref_mask_h].copy()
        if not ref_rows.empty:
            ref_rows = _ensure_columns(ref_rows, ["protein", "ligand", "best_affinity_kcal_mol", "affinity_kcal_mol"], default=np.nan)
            ref_rows["best_affinity_kcal_mol"] = pd.to_numeric(ref_rows.get("best_affinity_kcal_mol"), errors="coerce")
            ref_rows["affinity_kcal_mol"] = pd.to_numeric(ref_rows.get("affinity_kcal_mol"), errors="coerce")
            ref_rows["selected_affinity"] = ref_rows["best_affinity_kcal_mol"].where(ref_rows["best_affinity_kcal_mol"].notna(), ref_rows["affinity_kcal_mol"])
            best_ref = _best_rows(ref_rows, ["protein"], "selected_affinity")
            diagonal_map = dict(zip(best_ref["protein"].astype(str), best_ref["ligand"].astype(str)))
        ref_matrix = matrix[references].copy()
        for protein in ref_matrix.index.astype(str):
            expected = diagonal_map.get(protein, "")
            for ligand in ref_matrix.columns.astype(str):
                if ligand != expected:
                    ref_matrix.loc[protein, ligand] = np.nan
        parts: List[pd.DataFrame] = [ref_matrix]
        spacer = pd.DataFrame(index=ref_matrix.index, data={" ": np.nan})
        parts.append(spacer)
        if series:
            parts.append(matrix[series].copy())
        combined = pd.concat(parts, axis=1)
    else:
        combined = matrix.copy()

    ref_cols, series_cols = _group_ligand_columns(combined.columns)
    if ref_cols and series_cols:
        spacer = pd.DataFrame(index=combined.index, data={" ": np.nan})
        combined = pd.concat([combined[ref_cols], spacer, combined[series_cols]], axis=1)
    combined = combined.loc[combined.mean(axis=1, skipna=True).sort_values().index]
    combined.index = [_pname(ctx, p) for p in combined.index]
    combined.columns = [_ligand_label(col) if str(col).strip() != "" else " " for col in combined.columns]
    vmin = float(min(-10.0, np.nanmin(combined.values))) if np.isfinite(np.nanmin(combined.values)) else -10.0
    vmax = float(max(2.0, np.nanmax(combined.values))) if np.isfinite(np.nanmax(combined.values)) else 2.0
    fig, ax = plt.subplots(figsize=(max(12, 0.7 * len(combined.columns)), max(6, 0.55 * len(combined.index))))
    sns.heatmap(
        combined,
        cmap="RdYlGn_r",
        center=0.0,
        vmin=vmin,
        vmax=vmax,
        annot=True,
        fmt=".2f",
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Affinity (kcal/mol)"},
        ax=ax,
    )
    ax.set_title("Full Affinity Heatmap")
    ax.set_xlabel("Ligands")
    ax.set_ylabel("Proteins")
    _apply_category_tick_layout(ax, axis="x", rotation=70, wrap_width=16, fontsize=8)
    _apply_category_tick_layout(ax, axis="y", wrap_width=24, fontsize=9)
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_series_affinity_heatmap(ctx: PlotContext, output: Path) -> None:
    matrix, references, series = _build_affinity_matrix(ctx)
    if matrix.empty:
        _save_placeholder(output, "Series Affinity Heatmap", "No affinity matrix available.")
        return
    if series:
        series_matrix = matrix[series].copy()
    else:
        series_matrix = matrix.copy()
    series_matrix = series_matrix.transpose()
    if series_matrix.empty:
        _save_placeholder(output, "Series Affinity Heatmap", "No series ligands available.")
        return
    series_matrix = series_matrix.loc[series_matrix.mean(axis=1, skipna=True).sort_values().index]
    # After transpose, columns are proteins — rename them for display
    series_matrix.columns = [_pname(ctx, p) for p in series_matrix.columns]
    fig, ax = plt.subplots(figsize=(max(8, 0.6 * len(series_matrix.columns)), max(6, 0.45 * len(series_matrix.index))))
    sns.heatmap(
        series_matrix,
        cmap="RdYlGn_r",
        center=0.0,
        annot=True,
        fmt=".2f",
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Affinity (kcal/mol)"},
        ax=ax,
    )
    ax.set_title("Series Affinity Heatmap")
    ax.set_xlabel("Proteins")
    ax.set_ylabel("Series Ligands")
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_top_ligand_profiles(ctx: PlotContext, output: Path) -> None:
    table = _polypharmacology_table(ctx)
    if table.empty:
        _save_placeholder(output, "Top Ligand Profiles", "No ligand-protein rows available.")
        return
    ligand_summary = (
        table.groupby("ligand", dropna=False)
        .agg(
            poly_score=("moderate_or_better", "sum"),
            mean_affinity=("selected_affinity", "mean"),
        )
        .reset_index()
    )
    ligand_summary.sort_values(["poly_score", "mean_affinity"], ascending=[False, True], inplace=True)
    top_n = int(min(10, len(ligand_summary)))
    selected_ligands = ligand_summary.head(top_n)["ligand"].astype(str).tolist()
    prof = table[table["ligand"].astype(str).isin(selected_ligands)].copy()
    if prof.empty:
        _save_placeholder(output, "Top Ligand Profiles", "No selected top ligands available.")
        return
    proteins_order = (
        prof.groupby("protein", dropna=False)["selected_affinity"]
        .mean()
        .sort_values()
        .index.astype(str)
        .tolist()
    )
    proteins_display = [_pname(ctx, p) for p in proteins_order]
    fig, ax = plt.subplots(figsize=_SUITE_STYLE["large_figsize"])
    for ligand in selected_ligands:
        sub = prof[prof["ligand"].astype(str) == ligand].copy()
        sub = sub.groupby("protein", dropna=False)["selected_affinity"].min().reindex(proteins_order)
        ax.plot(proteins_display, sub.values, marker="o", linewidth=1.8, label=ligand)
    ax.axhline(-7.0, color="#d62728", linestyle="--", linewidth=1.2, label="-7 threshold")
    ax.axhline(0.0, color="#ff9896", linestyle=":", linewidth=1.2, label="0 threshold")
    ax.set_title("Top Ligand Profiles Across Proteins")
    ax.set_xlabel("Protein")
    ax.set_ylabel("Best affinity (kcal/mol)")
    ax.tick_params(axis="x", rotation=55)
    ax.legend(loc="best", fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_polypharmacology_delta_heatmap(ctx: PlotContext, output: Path) -> None:
    matrix, references, series = _build_affinity_matrix(ctx)
    if matrix.empty or not references or not series:
        _save_placeholder(output, "Polypharmacology Delta Heatmap", "No reference/series matrix available.")
        return
    ref_rows = ctx.classified_hits.copy()
    mask = pd.Series(False, index=ref_rows.index)
    if "ligand_type" in ref_rows.columns:
        mask |= ref_rows["ligand_type"].astype(str).str.lower().eq("reference")
    if "tag" in ref_rows.columns:
        mask |= ref_rows["tag"].astype(str).str.contains("reference_", case=False, na=False)
    if "ligand" in ref_rows.columns:
        mask |= ref_rows["ligand"].astype(str).str.contains(_COCRYSTAL_LIGAND_RE)
    ref_rows = ref_rows[mask].copy()
    ref_rows = _ensure_columns(ref_rows, ["protein", "ligand", "best_affinity_kcal_mol", "affinity_kcal_mol"], default=np.nan)
    if ref_rows.empty:
        _save_placeholder(output, "Polypharmacology Delta Heatmap", "No reference rows detected.")
        return
    ref_rows["best_affinity_kcal_mol"] = pd.to_numeric(ref_rows.get("best_affinity_kcal_mol"), errors="coerce")
    ref_rows["affinity_kcal_mol"] = pd.to_numeric(ref_rows.get("affinity_kcal_mol"), errors="coerce")
    ref_rows["selected_affinity"] = ref_rows["best_affinity_kcal_mol"].where(ref_rows["best_affinity_kcal_mol"].notna(), ref_rows["affinity_kcal_mol"])
    protein_reference = _best_rows(ref_rows, ["protein"], "selected_affinity")[["protein", "selected_affinity"]]
    protein_reference.rename(columns={"selected_affinity": "reference_affinity"}, inplace=True)

    series_matrix = matrix[series].copy()
    delta = series_matrix.sub(
        protein_reference.set_index("protein")["reference_affinity"],
        axis=0,
    )
    delta = delta.transpose()
    if delta.empty:
        _save_placeholder(output, "Polypharmacology Delta Heatmap", "No delta values available.")
        return
    delta.columns = [_pname(ctx, p) for p in delta.columns]
    fig, ax = plt.subplots(figsize=(max(8, 0.55 * len(delta.columns)), max(6, 0.45 * len(delta.index))))
    vmax = np.nanmax(np.abs(delta.values)) if np.isfinite(np.nanmax(np.abs(delta.values))) else 5.0
    vmax = float(max(vmax, 1.0))
    sns.heatmap(
        delta,
        cmap="RdYlGn_r",
        center=0.0,
        vmin=-vmax,
        vmax=vmax,
        annot=True,
        fmt=".2f",
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Delta vs reference (kcal/mol)"},
        ax=ax,
    )
    ax.set_title("Polypharmacology Delta Heatmap (series - reference)")
    ax.set_xlabel("Protein")
    ax.set_ylabel("Series Ligand")
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_hit_class_matrix(ctx: PlotContext, output: Path) -> None:
    table = _polypharmacology_table(ctx)
    if table.empty:
        _save_placeholder(output, "Hit Class Matrix", "No class table available.")
        return
    pivot = (
        table.groupby(["ligand", "protein"], dropna=False)["docking_quality_class"]
        .first()
        .unstack(fill_value="Unclassified")
    )
    class_order = {"Unclassified": 0, "Inactive": 1, "Weak": 2, "Moderate": 3, "Strong": 4}
    numeric = pivot.apply(lambda column: column.map(lambda value: class_order.get(_normalize_class(value), 0)))
    selectivity = (pivot.apply(lambda column: column.map(lambda value: _normalize_class(value) in {"Strong", "Moderate"})).sum(axis=1)).astype(int)
    numeric = numeric.loc[selectivity.sort_values(ascending=False).index]
    strong_counts = (pivot.apply(lambda column: column.map(lambda value: _normalize_class(value) == "Strong")).sum(axis=0)).astype(int)
    numeric = numeric[strong_counts.sort_values(ascending=False).index]

    fig = plt.figure(figsize=(max(9, 0.55 * numeric.shape[1]), max(6, 0.45 * numeric.shape[0] + 1.5)))
    gs = fig.add_gridspec(2, 1, height_ratios=[max(1, numeric.shape[0]), 1], hspace=0.12)
    ax_main = fig.add_subplot(gs[0, 0])
    cmap = mcolors.ListedColormap(["#ffffff", "#ef5350", "#f9a825", "#9ccc65", "#2e7d32"])
    sns.heatmap(
        numeric,
        cmap=cmap,
        vmin=0,
        vmax=4,
        cbar=False,
        linewidths=0.4,
        linecolor="white",
        ax=ax_main,
    )
    ax_main.set_title("Hit Class Matrix")
    ax_main.set_xlabel("")
    ax_main.set_ylabel("Ligand")

    ax_bottom = fig.add_subplot(gs[1, 0])
    sel_df = pd.DataFrame([selectivity.reindex(numeric.index).values], columns=numeric.index, index=["Exploratory target count"])
    sns.heatmap(
        sel_df,
        cmap="Greens",
        annot=True,
        fmt="d",
        cbar=False,
        linewidths=0.4,
        linecolor="white",
        ax=ax_bottom,
    )
    ax_bottom.set_xlabel("Ligand")
    ax_bottom.set_ylabel("")
    ax_bottom.tick_params(axis="x", rotation=75)
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_affinity_by_protein(ctx: PlotContext, output: Path) -> None:
    best = _best_scored_frame(ctx, include_engine=True)
    if best.empty:
        _save_placeholder(output, "Affinity by Protein", "No score rows available.")
        return
    best["protein"] = best["protein"].astype(str)
    best["affinity_kcal_mol"] = pd.to_numeric(best.get("affinity_kcal_mol"), errors="coerce")
    best = best[np.isfinite(best["affinity_kcal_mol"])].copy()
    if best.empty:
        _save_placeholder(output, "Affinity by Protein", "No finite affinity values available.")
        return

    best["protein"] = best["protein"].map(lambda p: _pname(ctx, p))
    n_proteins = int(best["protein"].nunique())
    fig, ax = plt.subplots(figsize=(_figure_width_for_categories(n_proteins, minimum=15.0, per_item=0.9), _SUITE_STYLE["large_figsize"][1]))
    if ctx.is_multi_engine and "engine" in best.columns:
        sns.boxplot(data=best, x="protein", y="affinity_kcal_mol", hue="engine", ax=ax)
        sns.stripplot(
            data=best,
            x="protein",
            y="affinity_kcal_mol",
            hue="engine",
            dodge=True,
            alpha=0.28,
            size=2.5,
            ax=ax,
        )
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[: len(set(best["engine"]))], labels[: len(set(best["engine"]))], loc="best", fontsize=8)
    else:
        sns.boxplot(data=best, x="protein", y="affinity_kcal_mol", color="#4c72b0", ax=ax)
        sns.stripplot(data=best, x="protein", y="affinity_kcal_mol", color="#1f77b4", alpha=0.4, size=3, ax=ax)

    ann = _ensure_columns(ctx.classified_hits, ["protein", "ligand", "best_affinity_kcal_mol", "affinity_kcal_mol"], default=np.nan)
    if not ann.empty:
        ann["best_affinity_kcal_mol"] = pd.to_numeric(ann.get("best_affinity_kcal_mol"), errors="coerce")
        ann["affinity_kcal_mol"] = pd.to_numeric(ann.get("affinity_kcal_mol"), errors="coerce")
        ann["selected_affinity"] = ann["best_affinity_kcal_mol"].where(ann["best_affinity_kcal_mol"].notna(), ann["affinity_kcal_mol"])
        ann["reference_ligand"] = _ligand_reference_mask(ann)
        star_rows = _best_rows(ann, ["protein"], "selected_affinity")
        for _, row in star_rows.iterrows():
            raw_protein = str(row.get("protein", ""))
            protein = _pname(ctx, raw_protein)
            if protein not in list(best["protein"].astype(str).unique()):
                continue
            x_pos = list(best["protein"].astype(str).unique()).index(protein)
            y_val = float(row.get("selected_affinity", np.nan))
            if np.isfinite(y_val):
                ax.scatter([x_pos], [y_val], marker="*", s=120, color="#ffd700", edgecolor="black", linewidth=0.6, zorder=5)
                ax.annotate(
                    _ligand_label(row.get("ligand", ""), reference=bool(row.get("reference_ligand", False))),
                    xy=(x_pos, y_val),
                    xytext=(0, 8),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="white", alpha=0.7, linewidth=0.0),
                )
    ax.axhline(-7.0, color="#d62728", linestyle="--", linewidth=1.2)
    ax.axhline(0.0, color="#ff7f0e", linestyle=":", linewidth=1.2)
    ax.set_title("Affinity by Protein")
    ax.set_xlabel("Protein")
    ax.set_ylabel("Affinity (kcal/mol)")
    _apply_category_tick_layout(ax, axis="x", rotation=45, wrap_width=18, fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_best_ligand_per_protein(ctx: PlotContext, output: Path) -> None:
    frame = _ensure_columns(ctx.classified_hits, ["protein", "ligand", "best_affinity_kcal_mol", "affinity_kcal_mol"], default=np.nan)
    if frame.empty:
        _save_placeholder(output, "Best Ligand Per Protein", "No classified hits available.")
        return
    frame["best_affinity_kcal_mol"] = pd.to_numeric(frame.get("best_affinity_kcal_mol"), errors="coerce")
    frame["affinity_kcal_mol"] = pd.to_numeric(frame.get("affinity_kcal_mol"), errors="coerce")
    frame["selected_affinity"] = frame["best_affinity_kcal_mol"].where(frame["best_affinity_kcal_mol"].notna(), frame["affinity_kcal_mol"])
    frame["reference_ligand"] = _ligand_reference_mask(frame)
    best = _best_rows(frame, ["protein"], "selected_affinity")
    if best.empty:
        _save_placeholder(output, "Best Ligand Per Protein", "No best-per-protein rows available.")
        return
    best.sort_values(["selected_affinity", "protein"], inplace=True)
    best["protein"] = best["protein"].map(lambda p: _pname(ctx, p))
    fig, ax = plt.subplots(figsize=(_SUITE_STYLE["large_figsize"][0], max(6, 0.55 * len(best))))
    ax.barh(best["protein"], best["selected_affinity"], color="#4c72b0")
    for idx, row in best.reset_index(drop=True).iterrows():
        x_text = float(row["selected_affinity"]) + 0.08
        ax.text(
            x_text,
            idx,
            _ligand_label(row.get("ligand", ""), reference=bool(row.get("reference_ligand", False))),
            va="center",
            fontsize=8,
        )
    ax.set_title("Best Ligand Per Protein")
    ax.set_xlabel("Best affinity (kcal/mol)")
    ax.set_ylabel("Protein")
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_hit_class_distribution_per_protein(ctx: PlotContext, output: Path) -> None:
    frame = _ensure_columns(ctx.classified_hits, ["protein", "docking_quality_class"], default="")
    if frame.empty:
        _save_placeholder(output, "Hit Class Distribution Per Protein", "No class distribution rows available.")
        return
    frame["protein"] = frame["protein"].astype(str)
    frame["docking_quality_class"] = frame["docking_quality_class"].map(_normalize_class)
    counts = (
        frame.groupby(["protein", "docking_quality_class"], dropna=False)
        .size()
        .unstack(fill_value=0)
    )
    class_order = ["Strong", "Moderate", "Weak", "Inactive", "Unclassified"]
    for col in class_order:
        if col not in counts.columns:
            counts[col] = 0
    counts = counts[class_order]
    strong_order = counts["Strong"].sort_values(ascending=False).index
    counts = counts.loc[strong_order]
    fig, ax = plt.subplots(figsize=_SUITE_STYLE["large_figsize"])
    bottom = np.zeros(len(counts), dtype=float)
    x = np.arange(len(counts))
    for cls in class_order:
        values = counts[cls].to_numpy(dtype=float)
        ax.bar(x, values, bottom=bottom, color=_HIT_CLASS_COLORS.get(cls, "#bdbdbd"), label=cls)
        bottom += values
    ax.set_xticks(x)
    ax.set_xticklabels([_pname(ctx, p) for p in counts.index.tolist()], rotation=55)
    ax.set_title("Hit Class Distribution Per Protein")
    ax.set_xlabel("Protein")
    ax.set_ylabel("Ligand count")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_engine_rank_correlation(ctx: PlotContext, output: Path) -> None:
    matrix = _build_engine_rank_matrix(ctx)
    if matrix.empty:
        _save_placeholder(output, "Engine Rank Correlation", "Need at least three matched multi-engine complexes for a defensible rank-correlation heatmap.")
        return
    corr = matrix.corr(method="spearman")
    fig, ax = plt.subplots(figsize=(max(6, 1.2 * corr.shape[0]), max(5, 1.0 * corr.shape[0])))
    sns.heatmap(
        corr,
        cmap="RdBu",
        center=0.0,
        vmin=-1.0,
        vmax=1.0,
        annot=True,
        fmt=".2f",
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Spearman rank correlation"},
        ax=ax,
    )
    ax.set_title("Engine Rank Correlation")
    ax.text(
        0.98,
        0.02,
        f"Matched complexes: {len(matrix)}",
        ha="right",
        va="bottom",
        transform=ax.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7),
    )
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_engine_agreement_per_complex(ctx: PlotContext, output: Path) -> None:
    ligand_fraction = _build_engine_agreement_summary(ctx)
    if ligand_fraction.empty:
        _save_placeholder(output, "Engine Agreement Per Complex", "No ligand-level engine agreement rows were available.")
        return
    colors = []
    for value in ligand_fraction["full_support_fraction"]:
        if value >= 0.67:
            colors.append("#2e7d32")
        elif value >= 0.33:
            colors.append("#f9a825")
        else:
            colors.append("#d32f2f")
    fig, ax = plt.subplots(figsize=(_figure_width_for_categories(len(ligand_fraction), minimum=15.0, per_item=0.45), _SUITE_STYLE["large_figsize"][1]))
    ax.bar(ligand_fraction["display_ligand"], ligand_fraction["full_support_fraction"], color=colors)
    ax.set_ylim(0.0, 1.05)
    ax.set_title("Engine Agreement Per Ligand")
    ax.set_ylabel("Full-support fraction")
    ax.set_xlabel("Ligand")
    _apply_category_tick_layout(ax, axis="x", rotation=70, wrap_width=18, fontsize=8)
    for idx, row in ligand_fraction.reset_index(drop=True).iterrows():
        ax.text(
            idx,
            float(row["full_support_fraction"]) + 0.02,
            f"n={int(row['matched_complexes'])}\nmean={float(row['mean_agreement_fraction']):.2f}",
            ha="center",
            va="bottom",
            fontsize=7,
        )
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_per_protein_cross_engine_batches(ctx: PlotContext, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    best = _best_scored_frame(ctx, include_engine=True)
    if best.empty:
        _save_placeholder(output_dir / "cross_engine_placeholder.png", "Per-Protein Cross-Engine", "No score rows available.")
        return
    best["engine"] = best["engine"].astype(str)
    best["protein"] = best["protein"].astype(str)
    proteins = sorted(best["protein"].dropna().unique().tolist())
    generated = 0
    for protein in proteins:
        sub = best[best["protein"] == protein].copy()
        if sub["engine"].nunique() < 2:
            continue
        pivot = sub.pivot_table(index="ligand", columns="engine", values="affinity_kcal_mol", aggfunc="min")
        if pivot.empty:
            continue
        fig, ax = plt.subplots(figsize=(max(8, 0.4 * len(pivot.index)), 5))
        for engine in sorted(pivot.columns.astype(str).tolist()):
            ax.plot(pivot.index.astype(str), pivot[engine], marker="o", linewidth=1.4, label=engine.upper())
        ax.set_title(f"Cross-Engine Batch: {protein}")
        ax.set_xlabel("Ligand")
        ax.set_ylabel("Affinity (kcal/mol)")
        ax.tick_params(axis="x", rotation=70)
        ax.legend(loc="best", fontsize=8)
        fig.tight_layout()
        fig.savefig(output_dir / f"{protein}_cross_engine.png", dpi=int(_SUITE_STYLE["dpi"]))
        plt.close(fig)
        generated += 1
    if generated == 0:
        _save_placeholder(output_dir / "cross_engine_placeholder.png", "Per-Protein Cross-Engine", "No proteins had >=2 engines.")


def _plot_cnn_confidence_distribution(ctx: PlotContext, output: Path) -> None:
    best = _best_scored_frame(ctx, include_engine=True)
    if best.empty:
        _save_placeholder(output, "CNN Confidence Distribution", "No GNINA rows available.")
        return
    gnina = best[best["engine"].astype(str).str.lower() == "gnina"].copy()
    if gnina.empty:
        _save_placeholder(output, "CNN Confidence Distribution", "No GNINA rows available.")
        return
    score = pd.to_numeric(gnina.get("cnn_score"), errors="coerce")
    gnina["cnn_confidence"] = pd.cut(
        score,
        bins=[-np.inf, 0.3, 0.5, np.inf],
        labels=["low", "moderate", "high"],
        right=False,
    ).astype(str)
    gnina["cnn_confidence"] = gnina["cnn_confidence"].replace("nan", "low")
    counts = (
        gnina.groupby(["protein", "cnn_confidence"], dropna=False)
        .size()
        .unstack(fill_value=0)
    )
    for col in ("low", "moderate", "high"):
        if col not in counts.columns:
            counts[col] = 0
    counts = counts[["high", "moderate", "low"]]
    fig, ax = plt.subplots(figsize=(_figure_width_for_categories(len(counts), minimum=15.0, per_item=0.9), _SUITE_STYLE["large_figsize"][1]))
    bottom = np.zeros(len(counts), dtype=float)
    x = np.arange(len(counts))
    colors = {"high": "#2e7d32", "moderate": "#f9a825", "low": "#d32f2f"}
    for level in ("high", "moderate", "low"):
        values = counts[level].to_numpy(dtype=float)
        ax.bar(x, values, bottom=bottom, color=colors[level], label=level)
        bottom += values
    ax.set_xticks(x)
    ax.set_xticklabels([_pname(ctx, p) for p in counts.index.astype(str).tolist()])
    ax.set_title("GNINA CNN Confidence Distribution")
    ax.set_xlabel("Protein")
    ax.set_ylabel("Ligand count")
    ax.legend(loc="best", fontsize=8)
    _apply_category_tick_layout(ax, axis="x", rotation=45, wrap_width=18, fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_pose_diversity(ctx: PlotContext, output: Path) -> None:
    best = _best_scored_frame(ctx, include_engine=True)
    if best.empty:
        _save_placeholder(output, "Pose Diversity", "No Vina rows available.")
        return
    vina = best[best["engine"].astype(str).str.lower() == "vina"].copy()
    if vina.empty:
        _save_placeholder(output, "Pose Diversity", "No Vina rows available.")
        return
    vina["rmsd_lb"] = pd.to_numeric(vina.get("rmsd_lb"), errors="coerce")
    values = vina["rmsd_lb"].dropna()
    diagnostic = _pose_diversity_diagnostic(values)
    if diagnostic:
        _save_placeholder(output, "Pose Diversity", diagnostic)
        return
    fig, ax = plt.subplots(figsize=_SUITE_STYLE["figsize"])
    ax.hist(values, bins=20, color="#4c72b0", alpha=0.8)
    ax.set_title("Vina Pose Diversity (rmsd_lb)")
    ax.set_xlabel("rmsd_lb")
    ax.set_ylabel("Count")
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _has_rmsd_data(ctx: PlotContext) -> bool:
    """Return True if there are any finite RMSD values available to plot."""
    if ctx.validation_table is None or ctx.validation_table.empty:
        return False
    vals = pd.to_numeric(ctx.validation_table.get("redocking_rmsd_angstrom"), errors="coerce").dropna()
    return len(vals) > 0


def _gate_annotation(ctx: PlotContext) -> str:
    """Return a short annotation string describing the validation gate status."""
    state = ctx.validation_status
    return f"Validation gate: {state.upper()}"


def _validation_joined_table(ctx: PlotContext) -> pd.DataFrame:
    validation = _ensure_columns(
        ctx.validation_table,
        ["protein", "ligand", "tag", "engine", "redocking_rmsd_angstrom", "redocking_classification", "reference_affinity"],
        default=np.nan,
    )
    if validation.empty:
        return validation
    validation["redocking_rmsd_angstrom"] = pd.to_numeric(validation.get("redocking_rmsd_angstrom"), errors="coerce")
    scores = _ensure_columns(ctx.normalized_scores, ["tag", "engine", "affinity_kcal_mol"], default=np.nan)
    keys = ["tag", "engine"]
    if "pose" in validation and "pose" in scores:
        keys.append("pose")
        scores = scores.drop_duplicates(keys)
    else:
        scores = _best_rows(scores, keys, "affinity_kcal_mol")
    merged = validation.merge(scores[keys + ["affinity_kcal_mol"]], on=keys, how="left", validate="many_to_one")
    merged["affinity_kcal_mol"] = pd.to_numeric(merged["affinity_kcal_mol"], errors="coerce")
    baseline_keys = ["protein", "engine", "scoring_function"]
    if not ctx.reference_baselines.empty and all(key in ctx.reference_baselines and key in merged for key in baseline_keys):
        refs = ctx.reference_baselines[baseline_keys + ["reference_affinity"]].copy()
        if refs.duplicated(baseline_keys).any():
            raise ValueError("Ambiguous engine/profile reference baselines")
        merged = merged.merge(refs, on=baseline_keys, how="left", suffixes=("", "_baseline"), validate="many_to_one")
        merged["reference_affinity"] = pd.to_numeric(merged["reference_affinity_baseline"], errors="coerce")
        merged.drop(columns=["reference_affinity_baseline"], inplace=True)
    else:
        merged["reference_affinity"] = np.nan
    return merged


def _plot_rmsd_distribution(ctx: PlotContext, output: Path) -> None:
    data = _validation_joined_table(ctx)
    values = pd.to_numeric(data.get("redocking_rmsd_angstrom"), errors="coerce").dropna()
    if values.empty:
        _save_placeholder(output, "RMSD Distribution", "No RMSD values available.")
        return
    bins = np.linspace(0.0, max(float(values.max()) + 0.2, 4.5), 18)
    fig, ax = plt.subplots(figsize=_SUITE_STYLE["figsize"])
    ax.hist(values[values < 2.0], bins=bins, color="#2e7d32", alpha=0.8, label="< 2 A")
    ax.hist(values[(values >= 2.0) & (values <= 4.0)], bins=bins, color="#f9a825", alpha=0.8, label="2-4 A")
    ax.hist(values[values > 4.0], bins=bins, color="#d32f2f", alpha=0.8, label="> 4 A")
    ax.axvline(2.0, color="#2e7d32", linestyle="--", linewidth=1.2)
    ax.axvline(4.0, color="#d32f2f", linestyle="--", linewidth=1.2)
    pass_rate = float((values < 2.0).sum() / max(len(values), 1))
    ax.text(
        0.98,
        0.95,
        f"Empirical pass rate (<2 A): {pass_rate * 100:.1f}%\nN={len(values)}",
        ha="right",
        va="top",
        transform=ax.transAxes,
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7),
    )
    gate_label = _gate_annotation(ctx)
    color = "#2e7d32" if ctx.validation_passed else "#d32f2f"
    ax.text(
        0.02, 0.95, gate_label, ha="left", va="top", transform=ax.transAxes,
        fontsize=8, color=color,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7),
    )
    ax.set_title("Redocking RMSD Distribution")
    ax.set_xlabel("RMSD (angstrom)")
    ax.set_ylabel("Count")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_reference_vs_docked(ctx: PlotContext, output: Path) -> None:
    data = _validation_joined_table(ctx)
    x = pd.to_numeric(data.get("reference_affinity"), errors="coerce")
    y = pd.to_numeric(data.get("affinity_kcal_mol"), errors="coerce")
    if x.notna().sum() == 0:
        _save_placeholder(output, "Reference vs Docked", "Reference baselines are unavailable, so reference-vs-docked affinity cannot be plotted.")
        return
    if y.notna().sum() == 0:
        _save_placeholder(output, "Reference vs Docked", "Docked affinity values are unavailable for the validated reference rows.")
        return
    valid = data[x.notna() & y.notna()].copy()
    if valid.empty:
        _save_placeholder(output, "Reference vs Docked", "No reference/docked affinity pairs available.")
        return
    valid["reference_affinity"] = pd.to_numeric(valid.get("reference_affinity"), errors="coerce")
    valid["affinity_kcal_mol"] = pd.to_numeric(valid.get("affinity_kcal_mol"), errors="coerce")
    fig, ax = plt.subplots(figsize=_SUITE_STYLE["figsize"])
    ax.scatter(valid["reference_affinity"], valid["affinity_kcal_mol"], color="#4c72b0", alpha=0.8)
    if len(valid) <= 18:
        for _, row in valid.iterrows():
            ax.text(
                float(row["reference_affinity"]),
                float(row["affinity_kcal_mol"]),
                _pname(ctx, str(row.get("protein", ""))),
                fontsize=7,
                alpha=0.8,
            )
    min_lim = float(min(valid["reference_affinity"].min(), valid["affinity_kcal_mol"].min()))
    max_lim = float(max(valid["reference_affinity"].max(), valid["affinity_kcal_mol"].max()))
    ax.plot([min_lim, max_lim], [min_lim, max_lim], linestyle="--", color="#d62728", linewidth=1.2)
    gate_label = _gate_annotation(ctx)
    color = "#2e7d32" if ctx.validation_passed else "#d32f2f"
    ax.text(
        0.02, 0.98, gate_label, ha="left", va="top", transform=ax.transAxes,
        fontsize=8, color=color,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7),
    )
    ax.set_title("Reference vs Docked Affinity")
    ax.set_xlabel("Reference affinity (kcal/mol)")
    ax.set_ylabel("Docked affinity (kcal/mol)")
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_rmsd_per_complex(ctx: PlotContext, output: Path) -> None:
    data = _validation_joined_table(ctx)
    data["redocking_rmsd_angstrom"] = pd.to_numeric(data.get("redocking_rmsd_angstrom"), errors="coerce")
    data["affinity_kcal_mol"] = pd.to_numeric(data.get("affinity_kcal_mol"), errors="coerce")
    if data["redocking_rmsd_angstrom"].notna().sum() == 0:
        _save_placeholder(output, "RMSD Per Complex", "No redocking RMSD values are available for plotting.")
        return
    if data["affinity_kcal_mol"].notna().sum() == 0:
        _save_placeholder(output, "RMSD Per Complex", "Docked affinity values are unavailable for the validated rows, so RMSD-versus-affinity cannot be plotted.")
        return
    valid = data[data["redocking_rmsd_angstrom"].notna() & data["affinity_kcal_mol"].notna()].copy()
    if valid.empty:
        _save_placeholder(output, "RMSD Per Complex", "No RMSD-affinity rows available.")
        return
    colors = []
    for value in valid["redocking_rmsd_angstrom"]:
        if value < 2.0:
            colors.append("#2e7d32")
        elif value <= 4.0:
            colors.append("#f9a825")
        else:
            colors.append("#d32f2f")
    fig, ax = plt.subplots(figsize=_SUITE_STYLE["figsize"])
    ax.scatter(valid["affinity_kcal_mol"], valid["redocking_rmsd_angstrom"], color=colors, alpha=0.85)
    ax.axhline(2.0, color="#2e7d32", linestyle="--", linewidth=1.0, label="2 Å threshold")
    ax.axhline(4.0, color="#f9a825", linestyle="--", linewidth=1.0, label="4 Å threshold")
    ax.axvline(-7.0, color="#d62728", linestyle="--", linewidth=1.2, label="-7 kcal/mol")
    gate_label = _gate_annotation(ctx)
    color = "#2e7d32" if ctx.validation_passed else "#d32f2f"
    ax.text(
        0.02, 0.98, gate_label, ha="left", va="top", transform=ax.transAxes,
        fontsize=8, color=color,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7),
    )
    ax.set_title("RMSD Per Complex")
    ax.set_xlabel("Docked affinity (kcal/mol)")
    ax.set_ylabel("Redocking RMSD (angstrom)")
    if len(valid) <= 18:
        for _, row in valid.iterrows():
            ax.annotate(
                _pname(ctx, str(row.get("protein", ""))),
                xy=(float(row["affinity_kcal_mol"]), float(row["redocking_rmsd_angstrom"])),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=7,
            )
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=int(_SUITE_STYLE["dpi"]))
    plt.close(fig)


def _plot_per_engine_affinity_heatmaps(ctx: PlotContext, output_dir: Path) -> List[FigureResult]:
    """
    Generate one affinity heatmap per engine (protein × ligand).

    Returns a list of FigureResult entries, one per engine, for manifest registration.
    Uses engine_affinity_matrix which has columns: tag, {engine1}, {engine2}, ...
    """
    results: List[FigureResult] = []
    em = ctx.engine_affinity_matrix
    if em is None or em.empty or "tag" not in em.columns:
        return results

    engines = [str(e).strip().lower() for e in (ctx.engine_scope_config.get("engines_in_scope") or []) if str(e).strip()]
    engine_cols = [e for e in engines if e in em.columns]
    if not engine_cols:
        return results

    def _parse_tag(tag: str):
        parts = tag.split("_site_")
        if len(parts) >= 2:
            protein = parts[0]
            rest = "_".join(parts[1:])
            ligand = "_".join(rest.split("_")[1:]) if "_" in rest else rest
            return protein, ligand
        return None, None

    em = em.copy()
    em["_protein"] = em["tag"].map(lambda t: _parse_tag(str(t))[0])
    em["_ligand"] = em["tag"].map(lambda t: _parse_tag(str(t))[1])
    em = em.dropna(subset=["_protein", "_ligand"])

    # Determine vmin/vmax shared across all engines for a consistent color scale
    all_vals = pd.to_numeric(em[engine_cols].values.flatten(), errors="coerce")
    finite_vals = all_vals[np.isfinite(all_vals)]
    if len(finite_vals) == 0:
        return results
    vmin_shared = float(min(-10.0, np.nanmin(finite_vals)))
    vmax_shared = float(max(2.0, np.nanmax(finite_vals)))

    output_dir.mkdir(parents=True, exist_ok=True)
    for engine in engine_cols:
        sub = em[["_protein", "_ligand", engine]].copy()
        sub[engine] = pd.to_numeric(sub[engine], errors="coerce")
        try:
            pivot = sub.pivot_table(index="_protein", columns="_ligand", values=engine, aggfunc="min")
        except Exception:
            continue
        if pivot.empty:
            continue
        pivot = pivot.loc[pivot.mean(axis=1, skipna=True).sort_values().index]
        pivot.index = [_pname(ctx, p) for p in pivot.index]

        # Separate co-crystal reference columns from series
        ref_cols, series_cols = _group_ligand_columns(pivot.columns)
        if ref_cols and series_cols:
            spacer = pd.DataFrame(index=pivot.index, data={" ": np.nan})
            pivot = pd.concat([pivot[ref_cols], spacer, pivot[series_cols]], axis=1)
        pivot.columns = [_ligand_label(col) if str(col).strip() != "" else " " for col in pivot.columns]

        fig, ax = plt.subplots(figsize=(max(12, 0.7 * len(pivot.columns)), max(6, 0.55 * len(pivot.index))))
        sns.heatmap(
            pivot,
            cmap="RdYlGn_r",
            center=0.0,
            vmin=vmin_shared,
            vmax=vmax_shared,
            annot=True,
            fmt=".1f",
            linewidths=0.4,
            linecolor="white",
            cbar_kws={"label": "Affinity (kcal/mol)"},
            ax=ax,
        )
        ax.set_title(f"{engine.upper()} Affinity Heatmap")
        ax.set_xlabel("Ligands")
        ax.set_ylabel("Proteins")
        _apply_category_tick_layout(ax, axis="x", rotation=70, wrap_width=16, fontsize=8)
        _apply_category_tick_layout(ax, axis="y", wrap_width=24, fontsize=9)
        fig.tight_layout()
        out_path = output_dir / f"{engine}_affinity_heatmap.png"
        fig.savefig(out_path, dpi=int(_SUITE_STYLE["dpi"]))
        plt.close(fig)
        results.append(
            FigureResult(
                group="multi_target",
                name=out_path.name,
                path=str(out_path),
                title=f"{engine.upper()} Affinity Heatmap",
                conditional=True,
                generated=out_path.exists(),
            )
        )
    return results


def _write_figures_manifest(
    output_root: Path,
    *,
    engine_mode: str,
    engine_name: Optional[str],
    figures: List[FigureResult],
) -> Dict[str, object]:
    figures_manifest_path = output_root / "figures_manifest.json"
    generated_rows = [row for row in figures if row.generated]
    skipped_rows = [row for row in figures if not row.generated]
    payload = {
        "generated_at": _utc_now_iso(),
        "engine_mode": engine_mode,
        "engine": engine_name if engine_name else None,
        "figures": [
            {
                "group": row.group,
                "name": row.name,
                "path": row.path,
                "title": row.title,
                "conditional": bool(row.conditional),
                "generated": bool(row.generated),
            }
            for row in figures
        ],
        "skipped": [
            {
                "name": row.name,
                "reason": row.reason or row.error or "not_generated",
            }
            for row in skipped_rows
        ],
        "summary": {
            "total": int(len(figures)),
            "generated": int(len(generated_rows)),
            "skipped": int(len(skipped_rows)),
            "groups": {
                group: int(sum(1 for row in generated_rows if row.group == group))
                for group in sorted({row.group for row in figures})
            },
        },
    }
    output_root.mkdir(parents=True, exist_ok=True)
    figures_manifest_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def generate_visualization_suite(
    *,
    classified_hits_file: Path,
    consensus_ranked_file: Path,
    engine_agreement_file: Path,
    normalized_scores_file: Path,
    engine_scope_config_file: Path,
    validation_gate_file: Path,
    top_pose_global_file: Path,
    output_dir: Path,
    protein_name_mapping_file: Optional[Path] = None,
    figure_overrides: Optional[Mapping[str, Callable[[PlotContext, Path], None]]] = None,
) -> Dict[str, object]:
    """
    Generate the canonical visualization suite.

    This function never raises. Any per-figure exception is captured in the manifest.
    """
    output_root = Path(output_dir).expanduser().resolve()
    ctx = _build_plot_context(
        classified_hits_file=Path(classified_hits_file),
        consensus_ranked_file=Path(consensus_ranked_file),
        engine_agreement_file=Path(engine_agreement_file),
        normalized_scores_file=Path(normalized_scores_file),
        engine_scope_config_file=Path(engine_scope_config_file),
        validation_gate_file=Path(validation_gate_file),
        top_pose_global_file=Path(top_pose_global_file),
        output_dir=output_root,
        protein_name_mapping_file=Path(protein_name_mapping_file) if protein_name_mapping_file else None,
    )

    figures: List[FigureResult] = []
    if not _PLOTTING_AVAILABLE:
        for group, name, title in (
            ("overview", "score_distribution.png", "Score Distribution"),
            ("ranking", "consensus_leaderboard.png", "Consensus Leaderboard"),
        ):
            figures.append(
                FigureResult(
                    group=group,
                    name=name,
                    path=str(output_root / group / name),
                    title=title,
                    conditional=False,
                    generated=False,
                    reason="plotting_backend_unavailable",
                )
            )
        return _write_figures_manifest(
            output_root,
            engine_mode="multi_engine" if ctx.is_multi_engine else "single_engine",
            engine_name=ctx.engine_name,
            figures=figures,
        )

    _apply_suite_style()

    figures.append(
        _run_plot_function(
            ctx,
            group="overview",
            name="score_distribution.png",
            title="Score Distribution",
            fn=_plot_score_distribution,
            figure_overrides=figure_overrides,
        )
    )
    figures.append(
        _run_plot_function(
            ctx,
            group="ranking",
            name="consensus_leaderboard.png",
            title="Consensus Leaderboard",
            fn=_plot_consensus_leaderboard,
            figure_overrides=figure_overrides,
        )
    )
    figures.append(
        _run_plot_function(
            ctx,
            group="ranking",
            name="affinity_per_ligand.png",
            title="Affinity Per Ligand",
            fn=_plot_affinity_per_ligand,
            figure_overrides=figure_overrides,
        )
    )
    figures.append(
        _run_plot_function(
            ctx,
            group="multi_target",
            name="full_affinity_heatmap.png",
            title="Full Affinity Heatmap",
            fn=_plot_full_affinity_heatmap,
            figure_overrides=figure_overrides,
        )
    )
    figures.append(
        _run_plot_function(
            ctx,
            group="multi_target",
            name="series_affinity_heatmap.png",
            title="Series Affinity Heatmap",
            fn=_plot_series_affinity_heatmap,
            figure_overrides=figure_overrides,
        )
    )
    # Per-engine affinity heatmaps (multi-engine only, placed in multi_target group)
    if ctx.is_multi_engine and not ctx.engine_affinity_matrix.empty:
        try:
            per_engine_dir = output_root / "multi_target"
            engine_fig_results = _plot_per_engine_affinity_heatmaps(ctx, per_engine_dir)
            figures.extend(engine_fig_results)
        except Exception as exc:
            logger.warning("Visualization suite: per-engine heatmaps failed (%s)", exc)

    figures.append(
        _run_plot_function(
            ctx,
            group="multi_target",
            name="top_ligand_profiles.png",
            title="Top Ligand Profiles",
            fn=_plot_top_ligand_profiles,
            figure_overrides=figure_overrides,
        )
    )
    figures.append(
        _run_plot_function(
            ctx,
            group="multi_target",
            name="polypharmacology_delta_heatmap.png",
            title="Polypharmacology Delta Heatmap",
            fn=_plot_polypharmacology_delta_heatmap,
            figure_overrides=figure_overrides,
            conditional=True,
            skip_reason="" if ctx.has_references else "no_references",
        )
    )
    figures.append(
        _run_plot_function(
            ctx,
            group="multi_target",
            name="hit_class_matrix.png",
            title="Hit Class Matrix",
            fn=_plot_hit_class_matrix,
            figure_overrides=figure_overrides,
        )
    )
    figures.append(
        _run_plot_function(
            ctx,
            group="per_protein",
            name="affinity_by_protein.png",
            title="Affinity by Protein",
            fn=_plot_affinity_by_protein,
            figure_overrides=figure_overrides,
        )
    )
    figures.append(
        _run_plot_function(
            ctx,
            group="per_protein",
            name="best_ligand_per_protein.png",
            title="Best Ligand Per Protein",
            fn=_plot_best_ligand_per_protein,
            figure_overrides=figure_overrides,
        )
    )
    figures.append(
        _run_plot_function(
            ctx,
            group="per_protein",
            name="hit_class_distribution_per_protein.png",
            title="Hit Class Distribution Per Protein",
            fn=_plot_hit_class_distribution_per_protein,
            figure_overrides=figure_overrides,
        )
    )

    if ctx.is_multi_engine:
        figures.append(
            _run_plot_function(
                ctx,
                group="engine_agreement",
                name="engine_rank_correlation.png",
                title="Engine Rank Correlation",
                fn=_plot_engine_rank_correlation,
                figure_overrides=figure_overrides,
                conditional=True,
            )
        )
        figures.append(
            _run_plot_function(
                ctx,
                group="engine_agreement",
                name="engine_agreement_per_complex.png",
                title="Engine Agreement Per Ligand",
                fn=_plot_engine_agreement_per_complex,
                figure_overrides=figure_overrides,
                conditional=True,
            )
        )
        # directory-level output for cross-engine batches
        cross_dir = output_root / "engine_agreement" / "per_protein_cross_engine"
        try:
            _plot_per_protein_cross_engine_batches(ctx, cross_dir)
            generated = bool(cross_dir.exists())
            figures.append(
                FigureResult(
                    group="engine_agreement",
                    name="per_protein_cross_engine",
                    path=str(cross_dir),
                    title="Per-Protein Cross-Engine Batches",
                    conditional=True,
                    generated=generated,
                    reason="" if generated else "directory_not_written",
                )
            )
        except Exception as exc:
            figures.append(
                FigureResult(
                    group="engine_agreement",
                    name="per_protein_cross_engine",
                    path=str(cross_dir),
                    title="Per-Protein Cross-Engine Batches",
                    conditional=True,
                    generated=False,
                    reason="figure_exception",
                    error=str(exc),
                )
            )
    # CNN confidence and pose diversity: generate whenever the relevant engine was used,
    # regardless of single- vs multi-engine mode.
    _scope_engines = [str(e).strip().lower() for e in (ctx.engine_scope_config.get("engines_in_scope") or [])]
    figures.append(
        _run_plot_function(
            ctx,
            group="engine_agreement",
            name="cnn_confidence_distribution.png",
            title="CNN Confidence Distribution",
            fn=_plot_cnn_confidence_distribution,
            figure_overrides=figure_overrides,
            conditional=True,
            skip_reason="" if "gnina" in _scope_engines else "gnina_not_in_scope",
        )
    )
    figures.append(
        _run_plot_function(
            ctx,
            group="engine_agreement",
            name="pose_diversity.png",
            title="Pose Diversity",
            fn=_plot_pose_diversity,
            figure_overrides=figure_overrides,
            conditional=True,
            skip_reason="" if "vina" in _scope_engines else "vina_not_in_scope",
        )
    )

    figures.append(
        _run_plot_function(
            ctx,
            group="validation",
            name="rmsd_distribution.png",
            title="RMSD Distribution",
            fn=_plot_rmsd_distribution,
            figure_overrides=figure_overrides,
            conditional=True,
            skip_reason="" if _has_rmsd_data(ctx) else "no_rmsd_data",
        )
    )
    figures.append(
        _run_plot_function(
            ctx,
            group="validation",
            name="reference_vs_docked.png",
            title="Reference vs Docked",
            fn=_plot_reference_vs_docked,
            figure_overrides=figure_overrides,
            conditional=True,
            skip_reason="" if _has_rmsd_data(ctx) else "no_rmsd_data",
        )
    )
    figures.append(
        _run_plot_function(
            ctx,
            group="validation",
            name="rmsd_per_complex.png",
            title="RMSD Per Complex",
            fn=_plot_rmsd_per_complex,
            figure_overrides=figure_overrides,
            conditional=True,
            skip_reason="" if _has_rmsd_data(ctx) else "no_rmsd_data",
        )
    )

    return _write_figures_manifest(
        output_root,
        engine_mode="multi_engine" if ctx.is_multi_engine else "single_engine",
        engine_name=ctx.engine_name,
        figures=figures,
    )
