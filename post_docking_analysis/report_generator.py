"""
Report generator for post-docking analysis pipeline.

This module handles the generation of summary reports in various formats.
"""
import pandas as pd
from pathlib import Path
from typing import Dict, List, Mapping, Optional
import json
import csv
from datetime import datetime, timezone


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _safe_relative(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except Exception:
        return str(path)


def generate_hit_classification_reports(
    consensus_df: pd.DataFrame,
    output_dir: Path,
    classifier_config: Dict[str, object],
    output_stem: str = "consensus_hit_classification",
) -> Dict[str, str]:
    """
    Write classifier configuration and class-distribution reports.
    """
    reports_dir = Path(output_dir).expanduser().resolve()
    reports_dir.mkdir(parents=True, exist_ok=True)

    config_payload = {
        "classifier_config": dict(classifier_config or {}),
        "input_rows": int(len(consensus_df)) if isinstance(consensus_df, pd.DataFrame) else 0,
    }
    config_file = reports_dir / f"{output_stem}_config.json"
    config_file.write_text(json.dumps(config_payload, indent=2), encoding="utf-8")

    if not isinstance(consensus_df, pd.DataFrame) or consensus_df.empty:
        distribution_df = pd.DataFrame(
            columns=["docking_quality_class", "row_count", "fraction_of_total", "qc_pass_fraction", "admet_pass_fraction"]
        )
        qc_breakdown_df = pd.DataFrame(columns=["docking_quality_class", "qc_status", "admet_status", "row_count"])
    else:
        frame = consensus_df.copy()
        frame["docking_quality_class"] = frame.get("docking_quality_class", pd.Series(["unclassified"] * len(frame))).astype(str)
        frame["qc_status"] = frame.get("qc_status", pd.Series(["unknown"] * len(frame))).astype(str)
        frame["admet_status"] = frame.get("admet_status", pd.Series(["unknown"] * len(frame))).astype(str)
        total_rows = max(int(len(frame)), 1)

        grouped = frame.groupby("docking_quality_class", dropna=False)
        distribution_df = grouped.agg(
            row_count=("docking_quality_class", "size"),
            qc_pass_fraction=("qc_status", lambda col: float((col == "pass").sum()) / max(len(col), 1)),
            admet_pass_fraction=("admet_status", lambda col: float((col == "pass").sum()) / max(len(col), 1)),
        ).reset_index()
        distribution_df["fraction_of_total"] = distribution_df["row_count"].astype(float) / float(total_rows)
        distribution_df = distribution_df[
            ["docking_quality_class", "row_count", "fraction_of_total", "qc_pass_fraction", "admet_pass_fraction"]
        ].sort_values(["row_count", "docking_quality_class"], ascending=[False, True])

        qc_breakdown_df = (
            frame.groupby(["docking_quality_class", "qc_status", "admet_status"], dropna=False)
            .size()
            .reset_index(name="row_count")
            .sort_values(["docking_quality_class", "row_count"], ascending=[True, False])
        )

    distribution_file = reports_dir / f"{output_stem}_distribution.csv"
    qc_breakdown_file = reports_dir / f"{output_stem}_qc_breakdown.csv"
    distribution_df.to_csv(distribution_file, index=False)
    qc_breakdown_df.to_csv(qc_breakdown_file, index=False)

    return {
        "classification_config_file": str(config_file),
        "classification_distribution_file": str(distribution_file),
        "classification_qc_breakdown_file": str(qc_breakdown_file),
    }


def generate_consolidated_run_summary(
    output_dir: Path,
    *,
    run_context: Mapping[str, object],
    engine_summary: Optional[pd.DataFrame] = None,
    classification_distribution: Optional[pd.DataFrame] = None,
    qc_breakdown: Optional[pd.DataFrame] = None,
    output_stem: str = "consolidated_run_summary",
) -> Dict[str, str]:
    """
    Emit consolidated run-level summary artifacts (JSON + CSV + optional table mirrors).
    """
    reports_dir = Path(output_dir).expanduser().resolve()
    reports_dir.mkdir(parents=True, exist_ok=True)

    payload: Dict[str, object] = {
        "generated_at": _utc_now_iso(),
        "run_context": dict(run_context or {}),
        "engine_count": int(len(engine_summary)) if isinstance(engine_summary, pd.DataFrame) else 0,
        "classification_rows": int(len(classification_distribution))
        if isinstance(classification_distribution, pd.DataFrame)
        else 0,
        "qc_breakdown_rows": int(len(qc_breakdown)) if isinstance(qc_breakdown, pd.DataFrame) else 0,
    }

    json_file = reports_dir / f"{output_stem}.json"
    csv_file = reports_dir / f"{output_stem}.csv"
    json_file.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    rows = [{"field": str(key), "value": json.dumps(value) if isinstance(value, (dict, list)) else str(value)} for key, value in payload.items()]
    pd.DataFrame(rows, columns=["field", "value"]).to_csv(csv_file, index=False)

    artifacts: Dict[str, str] = {
        "consolidated_summary_json_file": str(json_file),
        "consolidated_summary_csv_file": str(csv_file),
    }
    if isinstance(engine_summary, pd.DataFrame) and not engine_summary.empty:
        engine_file = reports_dir / "consolidated_engine_summary.csv"
        engine_summary.to_csv(engine_file, index=False)
        artifacts["consolidated_engine_summary_file"] = str(engine_file)
    if isinstance(classification_distribution, pd.DataFrame) and not classification_distribution.empty:
        dist_file = reports_dir / "consolidated_hit_class_distribution.csv"
        classification_distribution.to_csv(dist_file, index=False)
        artifacts["consolidated_hit_class_distribution_file"] = str(dist_file)
    if isinstance(qc_breakdown, pd.DataFrame) and not qc_breakdown.empty:
        qc_file = reports_dir / "consolidated_qc_breakdown.csv"
        qc_breakdown.to_csv(qc_file, index=False)
        artifacts["consolidated_qc_breakdown_file"] = str(qc_file)
    return artifacts


def generate_dashboard_index(
    output_dir: Path,
    *,
    run_context: Mapping[str, object],
    artifact_paths: Mapping[str, object],
    output_stem: str = "dashboard_export_index",
) -> Dict[str, str]:
    """
    Emit dashboard-oriented JSON index and a simple contract markdown.
    """
    reports_dir = Path(output_dir).expanduser().resolve()
    reports_dir.mkdir(parents=True, exist_ok=True)

    entries = []
    for key, value in artifact_paths.items():
        if value in (None, ""):
            continue
        path = Path(str(value)).expanduser()
        entries.append(
            {
                "key": str(key),
                "path": str(path),
                "exists": bool(path.exists()),
            }
        )

    index_payload = {
        "generated_at": _utc_now_iso(),
        "run_context": dict(run_context or {}),
        "entries": entries,
    }
    index_file = reports_dir / f"{output_stem}.json"
    index_file.write_text(json.dumps(index_payload, indent=2), encoding="utf-8")

    contract_file = reports_dir / "DASHBOARD_EXPORT_CONTRACT.md"
    lines = [
        "# Dashboard Export Contract",
        "",
        "This file describes the JSON index consumed by dashboard workflows.",
        "",
        "## JSON File",
        f"- `{index_file.name}`",
        "",
        "## Schema",
        "- `generated_at`: ISO-8601 UTC timestamp.",
        "- `run_context`: key/value metadata for run scope and scoring modes.",
        "- `entries[]`: array of export pointers:",
        "  - `key`: logical artifact identifier.",
        "  - `path`: absolute path to produced artifact.",
        "  - `exists`: boolean file existence check at generation time.",
        "",
        "## Notes",
        "- Consumers should tolerate missing optional artifacts (`exists=false`).",
        "- Paths are absolute to support interactive/local dashboard loading.",
        "",
    ]
    contract_file.write_text("\n".join(lines), encoding="utf-8")
    return {
        "dashboard_index_file": str(index_file),
        "dashboard_contract_file": str(contract_file),
    }


def generate_start_here_index(
    reports_root: Path,
    *,
    project_root: Path,
    run_context: Mapping[str, object],
    artifact_paths: Mapping[str, object],
    output_name: str = "START_HERE.md",
) -> str:
    """
    Generate a final report navigation index in `7-Reports/START_HERE.md`.
    """
    reports_dir = Path(reports_root).expanduser().resolve()
    reports_dir.mkdir(parents=True, exist_ok=True)
    project_dir = Path(project_root).expanduser().resolve()

    key_rows = []
    for key, value in artifact_paths.items():
        if value in (None, ""):
            continue
        path = Path(str(value)).expanduser()
        key_rows.append((str(key), path, bool(path.exists())))
    key_rows.sort(key=lambda row: row[0])

    lines = [
        "# DockForge Reports Index",
        "",
        "## Run Context",
        f"- generated_at: `{_utc_now_iso()}`",
        f"- analysis_mode: `{run_context.get('analysis_mode', '')}`",
        f"- analysis_scope: `{run_context.get('analysis_scope', '')}`",
        f"- consensus_mode: `{run_context.get('consensus_mode', '')}`",
        f"- normalization_method: `{run_context.get('normalization_method', '')}`",
        f"- engines_in_scope: `{', '.join(run_context.get('engines_in_scope', []) or [])}`",
        f"- scope_source: `{run_context.get('scope_source', '')}`",
        f"- scoped_engine_count: `{run_context.get('scoped_engine_count', '')}`",
        f"- detected_engine_count: `{run_context.get('detected_engine_count', '')}`",
        "",
        "## Engine Scope",
    ]
    excluded = run_context.get("excluded_engines") or []
    if excluded:
        for entry in excluded:
            if isinstance(entry, dict):
                lines.append(
                    f"- excluded `{entry.get('engine', '')}` ({entry.get('exclusion_reason', 'user_excluded')})"
                )
            else:
                lines.append(f"- excluded `{entry}`")
    else:
        lines.append("- no excluded engines")
    notice = str(run_context.get("scope_discrepancy_notice", "") or "").strip()
    if notice:
        lines.extend(
            [
                "",
                "## Scope Notice",
                f"- {notice}",
            ]
        )
    lines.extend(
        [
            "",
        "## Key Artifacts",
        ]
    )
    for key, path, exists in key_rows:
        rel = _safe_relative(path, project_dir)
        status = "exists" if exists else "missing"
        lines.append(f"- `{key}` -> `{rel}` ({status})")
    lines.extend(
        [
            "",
            "## Notes",
            "- This index is generated automatically for project-level navigation.",
            "- Missing artifacts may indicate optional features were disabled for this run.",
            "",
        ]
    )

    output_file = reports_dir / output_name
    output_file.write_text("\n".join(lines), encoding="utf-8")
    return str(output_file)


def generate_analysis_root_index(
    analysis_root: Path,
    *,
    project_root: Path,
    session_root: Path,
    reports_root: Path,
    output_name: str = "START_HERE.md",
) -> str:
    """
    Generate the canonical analysis navigation file in `5-Analysis/START_HERE.md`.
    """
    analysis_dir = Path(analysis_root).expanduser().resolve()
    analysis_dir.mkdir(parents=True, exist_ok=True)
    project_dir = Path(project_root).expanduser().resolve()
    latest_session = Path(session_root).expanduser().resolve()
    reports_dir = Path(reports_root).expanduser().resolve()
    legacy_root = project_dir / "5-Post_Docking_Analysis"

    key_paths = [
        ("latest_session", latest_session, "Most recent analysis session."),
        ("session_complexes", latest_session / "complexes", "MD-ready receptor-pose complex PDB files."),
        ("session_scores", latest_session / "scores", "Per-session score tables and normalized outputs."),
        ("session_best_poses", latest_session / "best_poses", "Best-pose exports and selection tables."),
        ("session_interactions", latest_session / "interactions", "Interaction analyses (ProLIF, PandaMap, PoseView, PyMOL)."),
        ("session_visualizations", latest_session / "visualizations", "Mirrored 2D/3D visual outputs for the latest session."),
        ("session_reports", latest_session / "reports", "Per-session reports and summaries."),
        ("comparative_root", analysis_dir / "comparative", "Cross-run comparative hit ranking tables."),
        ("top_pose_root", analysis_dir / "top_pose_ligand_performance", "Top-pose atlas outputs for hit triage."),
        ("polypharmacology_root", analysis_dir / "polypharmacology", "Multi-target docking-biology summaries."),
        ("structure_quality_root", analysis_dir / "structure_quality", "Validation and structure-quality outputs."),
        ("reports_root", reports_dir, "Consolidated cross-run reports and dashboard exports."),
    ]

    lines = [
        "# DockForge Analysis Index",
        "",
        f"- generated_at: `{_utc_now_iso()}`",
        f"- latest_session_id: `{latest_session.name}`",
        "",
        "## Open This First",
        f"- Latest session folder: `{_safe_relative(latest_session, project_dir)}`",
        f"- MD-ready complexes: `{_safe_relative(latest_session / 'complexes', project_dir)}`",
        f"- Top-pose atlas: `{_safe_relative(analysis_dir / 'top_pose_ligand_performance', project_dir)}`",
        f"- Consolidated reports: `{_safe_relative(reports_dir, project_dir)}`",
        "",
        "## Key Artifact Map",
    ]
    for key, path, description in key_paths:
        lines.append(f"- `{key}` -> `{_safe_relative(path, project_dir)}`")
        lines.append(f"  {description}")

    if legacy_root.exists() and not legacy_root.is_symlink():
        lines.extend(
            [
                "",
                "## Legacy Warning",
                f"- Existing legacy directory detected at `{_safe_relative(legacy_root, project_dir)}`.",
                "- It was left untouched for safety. New runs write to `5-Analysis/`.",
            ]
        )

    output_file = analysis_dir / output_name
    output_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    raw_data_root = analysis_dir / "raw_data"
    raw_data_root.mkdir(parents=True, exist_ok=True)
    (raw_data_root / output_name).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(output_file)

def generate_csv_reports(analysis_results: Dict[str, pd.DataFrame], output_dir: Path):
    """
    Generate CSV reports from analysis results.
    
    Parameters
    ----------
    analysis_results : Dict[str, pd.DataFrame]
        Dictionary containing analysis results
    output_dir : Path
        Output directory for reports
    """
    print("📝 Generating CSV reports...")
    
    # Create output directory
    reports_dir = output_dir / "reports"
    reports_dir.mkdir(exist_ok=True)
    
    # Generate individual CSV files for each result type
    for name, df in analysis_results.items():
        if isinstance(df, pd.DataFrame):
            csv_file = reports_dir / f"{name}.csv"
            df.to_csv(csv_file, index=False)
            print(f"✅ {name} report saved to: {csv_file}")
    
    print("✅ CSV reports generated successfully!")

def generate_excel_report(analysis_results: Dict[str, pd.DataFrame], output_dir: Path):
    """
    Generate an Excel report from analysis results.
    
    Parameters
    ----------
    analysis_results : Dict[str, pd.DataFrame]
        Dictionary containing analysis results
    output_dir : Path
        Output directory for reports
    """
    print("📝 Generating Excel report...")
    
    # Check if openpyxl is available
    try:
        import openpyxl
    except ImportError:
        print("⚠️  openpyxl not available - Excel report generation skipped")
        return
    
    # Create output directory
    reports_dir = output_dir / "reports"
    reports_dir.mkdir(exist_ok=True)
    
    # Generate Excel file with multiple sheets
    excel_file = reports_dir / "docking_analysis_results.xlsx"
    
    with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
        for name, df in analysis_results.items():
            if isinstance(df, pd.DataFrame):
                # Limit sheet name to 31 characters (Excel limit)
                sheet_name = name[:31]
                df.to_excel(writer, sheet_name=sheet_name, index=False)
    
    print(f"✅ Excel report saved to: {excel_file}")

def generate_summary_report(analysis_results: Dict[str, pd.DataFrame], output_dir: Path):
    """
    Generate a summary report with key findings.
    
    Parameters
    ----------
    analysis_results : Dict[str, pd.DataFrame]
        Dictionary containing analysis results
    output_dir : Path
        Output directory for reports
    """
    print("📝 Generating summary report...")
    
    # Create output directory
    reports_dir = output_dir / "reports"
    reports_dir.mkdir(exist_ok=True)
    
    # Get key results
    best_poses = analysis_results['best_poses']
    summary_stats = analysis_results['summary_stats']
    
    # Generate summary text
    summary_lines = [
        "Post-Docking Analysis Summary Report",
        "==================================",
        "",
        f"Total complexes analyzed: {len(best_poses)}",
        f"Average binding affinity: {best_poses['vina_affinity'].mean():.2f} kcal/mol",
        f"Best binding affinity: {best_poses['vina_affinity'].min():.2f} kcal/mol",
        f"Worst binding affinity: {best_poses['vina_affinity'].max():.2f} kcal/mol",
        "",
        "Top 5 Performers:",
    ]
    
    # Add top 5 performers
    top_5 = best_poses.head(5)
    for idx, (_, row) in enumerate(top_5.iterrows(), 1):
        summary_lines.append(
            f"  {idx}. {row['complex_name']}: {row['vina_affinity']:.2f} kcal/mol"
        )
    
    # Save summary report
    summary_file = reports_dir / "summary_report.txt"
    with open(summary_file, 'w') as f:
        f.write('\n'.join(summary_lines))
    
    print(f"✅ Summary report saved to: {summary_file}")

def generate_json_report(analysis_results: Dict[str, pd.DataFrame], output_dir: Path):
    """
    Generate a JSON report from analysis results.
    
    Parameters
    ----------
    analysis_results : Dict[str, pd.DataFrame]
        Dictionary containing analysis results
    output_dir : Path
        Output directory for reports
    """
    print("📝 Generating JSON report...")
    
    # Create output directory
    reports_dir = output_dir / "reports"
    reports_dir.mkdir(exist_ok=True)
    
    # Convert DataFrames to dictionaries for JSON serialization
    json_data = {}
    for name, df in analysis_results.items():
        if isinstance(df, pd.DataFrame):
            json_data[name] = df.to_dict(orient='records')
    
    # Save JSON report
    json_file = reports_dir / "analysis_results.json"
    with open(json_file, 'w') as f:
        json.dump(json_data, f, indent=2)
    
    print(f"✅ JSON report saved to: {json_file}")

def generate_all_reports(analysis_results: Dict[str, pd.DataFrame], output_dir: Path):
    """
    Generate all reports.
    
    Parameters
    ----------
    analysis_results : Dict[str, pd.DataFrame]
        Dictionary containing analysis results
    output_dir : Path
        Output directory for reports
    """
    print("📝 Generating all reports...")
    
    generate_csv_reports(analysis_results, output_dir)
    generate_excel_report(analysis_results, output_dir)
    generate_summary_report(analysis_results, output_dir)
    generate_json_report(analysis_results, output_dir)
    
    print("✅ All reports generated successfully!")
