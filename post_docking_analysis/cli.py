"""
Command-line interface for post-docking analysis pipeline.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _prompt_cli_engine_fallback(project_dir: str | Path) -> str | None:
    try:
        from docking.project_layout import load_manifest
        from .engine_detector import SUPPORTED_ENGINES
    except ImportError:
        from docking.project_layout import load_manifest
        from engine_detector import SUPPORTED_ENGINES

    if not sys.stdin.isatty():
        return None

    manifest = load_manifest(Path(project_dir).expanduser().resolve())
    declared = [
        str(engine or "").strip().lower()
        for engine in (manifest.get("engines") or [])
        if str(engine or "").strip().lower()
    ]
    candidates = []
    for engine in declared + list(SUPPORTED_ENGINES):
        if engine and engine not in candidates:
            candidates.append(engine)
    if not candidates:
        return None

    print("No valid docking-engine outputs were auto-detected.")
    print("Manual fallback engine selection is available because --interactive-fallback was passed.")
    for idx, engine in enumerate(candidates, start=1):
        print(f"  {idx}. {engine}")

    while True:
        raw = input("Select engine number (blank to cancel): ").strip()
        if not raw:
            return None
        try:
            index = int(raw)
        except ValueError:
            print("Please enter a valid number.")
            continue
        if 1 <= index <= len(candidates):
            return candidates[index - 1]
        print("Selection out of range.")

def main():
    """
    Main function to run the pipeline from command line.
    """
    parser = argparse.ArgumentParser(
        description="Post-Docking Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default settings
  python -m post_docking_analysis
  
  # Specify input and output directories
  python -m post_docking_analysis -i /path/to/docking/results -o /path/to/output
  
  # Use a configuration file
  python -m post_docking_analysis --config my_config.json
  
  # Run preprocessing to generate all_scores.csv
  python -m post_docking_analysis.preprocess /path/to/docking/results
  
  # Run with specific options
  python -m post_docking_analysis -i /path/to/results --no-visualizations
        """
    )
    
    # Input/Output arguments
    parser.add_argument("-i", "--input", 
                        help="Input directory containing docking results")
    parser.add_argument("-o", "--output", 
                        help="Output directory for results")
    parser.add_argument("--project-dir",
                        help="Canonical multi-engine docking project directory")
    parser.add_argument("--config", 
                        help="Configuration file path")
    
    # Processing options
    parser.add_argument("--no-split", action="store_true",
                        help="Skip complex splitting step")
    parser.add_argument("--no-apo", action="store_true",
                        help="Skip apo protein extraction")
    parser.add_argument("--no-ligands", action="store_true",
                        help="Skip ligand extraction")
    parser.add_argument("--no-analysis", action="store_true",
                        help="Skip binding affinity analysis")
    parser.add_argument("--no-visualizations", action="store_true",
                        help="Skip visualization generation")
    parser.add_argument("--no-reports", action="store_true",
                        help="Skip report generation")
    parser.add_argument("--fix-chains", action="store_true",
                        help="Fix chain issues in structures")
    
    # Directory structure
    parser.add_argument("--structure", choices=["auto", "single", "multi"],
                        help="Directory structure type (auto/single/multi)")
    
    # Scoring options
    parser.add_argument("--no-cnn", action="store_true",
                        help="Disable CNN scoring")
    
    # Output format options
    parser.add_argument("--no-csv", action="store_true",
                        help="Disable CSV output")
    parser.add_argument("--no-excel", action="store_true",
                        help="Disable Excel output")
    parser.add_argument("--no-pdb", action="store_true",
                        help="Disable PDB output")
    parser.add_argument("--no-mol2", action="store_true",
                        help="Disable MOL2 output")
    
    # Preprocessing options
    parser.add_argument("--preprocess", action="store_true",
                        help="Run preprocessing to generate all_scores.csv and identify pairs")
    parser.add_argument("--force", action="store_true",
                        help="Force regeneration during preprocessing or bypass DAG cache for artifact-scope runs")
    parser.add_argument("--pairlist", 
                        help="Path to pairlist.csv for accurate receptor-ligand mapping")
    parser.add_argument("--analysis-mode",
                        choices=["auto", "single_engine", "comparative_all_engines", "favorite_engine_continue"],
                        help="Run engine-aware analysis for a canonical multi-engine project")
    parser.add_argument(
        "--analysis-scope",
        choices=["full", "comparison_only", "rescoring_only", "top_pose_only", "qc_only", "report_only"],
        default="full",
        help="Project-centric analysis scope (full/comparison/rescoring/top-pose/QC/report)",
    )
    parser.add_argument(
        "--scope",
        choices=["full", "comparison_only", "structures_only", "interactions", "report_only", "top_pose_only"],
        help="Experimental artifact-DAG scope request; executes the minimum subgraph for the requested artifact family",
    )
    parser.add_argument("--engine",
                        help="Engine override; when provided it forces single_engine mode for the named engine")
    parser.add_argument("--engines",
                        help="Comma-separated engine scope filter for this analysis session, for example: gnina,vina")
    parser.add_argument("--engine-preset",
                        help="Named engine scope preset from project_manifest.json[\"engine_presets\"]")
    parser.add_argument("--favorite-engine",
                        help="Favorite engine used by favorite_engine_continue")
    parser.add_argument("--redetect", action="store_true",
                        help="Force re-probing engine outputs even when detection is cached in the manifest")
    parser.add_argument("--interactive-fallback", action="store_true",
                        help="Allow fallback interaction when no valid engines are detected")
    parser.add_argument("--promote-exhaustive", action="store_true",
                        help="Emit an exhaustive rerun manifest from comparative engine results")
    parser.add_argument("--rerun-engine",
                        help="Engine to promote into exhaustive reruns")
    parser.add_argument("--top-per-protein", type=int, default=1,
                        help="Number of best ligands per protein to promote into the rerun manifest")
    parser.add_argument("--winner-only", action="store_true",
                        help="Only promote pairs where the target engine wins across engines")
    parser.add_argument("--min-affinity-advantage", type=float, default=0.0,
                        help="Require the target engine to beat the best competing engine by at least this many kcal/mol")
    parser.add_argument("--max-rerun-pairs", type=int, default=0,
                        help="Optional global cap on promoted exhaustive rerun pairs")
    parser.add_argument("--pair-allowlist",
                        help="Optional TXT/CSV file restricting promotion to explicit pair tags or receptor/site_id/ligand rows")
    parser.add_argument(
        "--consensus-mode",
        choices=["dockbox_geometric", "weighted_hybrid", "strict_consensus", "favorite_guardrails"],
        default="dockbox_geometric",
        help="Consensus policy for comparative hit ranking and rerun promotion",
    )
    parser.add_argument(
        "--rescoring-scope",
        choices=["top_n_per_protein", "top_n_global"],
        default="top_n_per_protein",
        help="Scope used to select rescoring candidate shortlist",
    )
    parser.add_argument(
        "--rescoring-top-n",
        type=int,
        default=3,
        help="Top-N candidates retained by the selected rescoring scope",
    )
    parser.add_argument(
        "--normalization-method",
        choices=["per_engine_rank", "per_engine_minmax", "per_engine_zscore"],
        default="per_engine_rank",
        help="Normalization method used before cross-engine consensus",
    )
    parser.add_argument(
        "--biology-file",
        help="Optional CSV/TSV/JSON table with biological annotations for correlation",
    )
    parser.add_argument(
        "--biology-mapping-mode",
        choices=["auto", "tag", "protein_ligand", "protein", "ligand"],
        default="auto",
        help="How to map biology annotations to docking rows",
    )
    parser.add_argument(
        "--hit-class-policy",
        choices=["target_percentile", "reference_anchor"],
        default="target_percentile",
        help="Target-aware hit classification policy",
    )
    parser.add_argument(
        "--hit-class-strong-percentile",
        type=float,
        default=0.10,
        help="Per-target percentile threshold for Strong hits",
    )
    parser.add_argument(
        "--hit-class-moderate-percentile",
        type=float,
        default=0.35,
        help="Per-target percentile threshold for Moderate hits",
    )
    parser.add_argument(
        "--top-pose-policy",
        choices=["hybrid", "best_affinity", "best_consensus"],
        default="best_affinity",
        help="Top-pose atlas selection policy",
    )
    parser.add_argument(
        "--top-pose-aggregation",
        choices=["best_target"],
        default="best_target",
        help="Top-pose atlas global aggregation mode",
    )
    parser.add_argument(
        "--best-pose-selection-metric",
        choices=["auto", "vina_affinity", "cnn_affinity"],
        default="auto",
        help="Metric used for best-pose complex PDB extraction",
    )
    parser.add_argument(
        "--complex-query",
        help="Optional complex filter query, for example: protein=2FVD;ligand=Sorafenib or exclude_ligand=Acetazolamide. Use 0/skip/all to disable.",
    )
    parser.add_argument("--prompt-protein-names", action="store_true",
                        help="Prompt for protein display names before analysis starts")
    parser.add_argument("--prompt-ligand-names", action="store_true",
                        help="Prompt for ligand display names before analysis starts")
    parser.add_argument(
        "--rmsd-scopes",
        default="per_complex",
        help="Comma-separated RMSD scopes to run (current default and recommended: per_complex)",
    )
    parser.add_argument(
        "--no-resume-rmsd",
        action="store_true",
        help="Disable RMSD checkpoint reuse and recompute selected scopes from scratch",
    )
    
    # Help and info
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose output")
    parser.add_argument("--version", action="version", version="Post-Docking Analysis Pipeline 1.0")
    
    args = parser.parse_args()

    if args.project_dir or args.analysis_mode:
        project_dir = args.project_dir or args.input
        if not project_dir:
            parser.error("--project-dir or --input is required for engine-aware analysis")
        try:
            from .unified_pipeline import UnifiedPostDockingPipeline
        except ImportError:
            from unified_pipeline import UnifiedPostDockingPipeline
        pipeline_kwargs = dict(
            project_dir=project_dir,
            output_dir=args.output,
            analysis_mode=args.analysis_mode or "auto",
            engine=args.engine,
            engines=args.engines,
            engine_preset=args.engine_preset,
            favorite_engine=args.favorite_engine,
            redetect=args.redetect,
            interactive_fallback=args.interactive_fallback,
            promote_exhaustive=args.promote_exhaustive,
            rerun_engine=args.rerun_engine,
            top_per_protein=args.top_per_protein,
            winner_only=args.winner_only,
            min_affinity_advantage=args.min_affinity_advantage,
            max_rerun_pairs=args.max_rerun_pairs,
            pair_allowlist=args.pair_allowlist,
            consensus_mode=args.consensus_mode,
            rescoring_scope=args.rescoring_scope,
            rescoring_top_n=args.rescoring_top_n,
            analysis_scope=args.analysis_scope,
            normalization_method=args.normalization_method,
            biology_file=args.biology_file,
            biology_mapping_mode=args.biology_mapping_mode,
            hit_class_policy=args.hit_class_policy,
            hit_class_strong_percentile=args.hit_class_strong_percentile,
            hit_class_moderate_percentile=args.hit_class_moderate_percentile,
            top_pose_selection_policy=args.top_pose_policy,
            top_pose_global_aggregation=args.top_pose_aggregation,
            best_pose_selection_metric=args.best_pose_selection_metric,
            complex_query=args.complex_query,
            prompt_protein_names=args.prompt_protein_names,
            prompt_ligand_names=args.prompt_ligand_names,
            rmsd_scopes=args.rmsd_scopes,
            resume_rmsd=not args.no_resume_rmsd,
            dag_scope=args.scope,
            dag_force=bool(args.force and args.scope),
        )
        try:
            engine_pipeline = UnifiedPostDockingPipeline(**pipeline_kwargs)
        except ValueError as exc:
            if args.interactive_fallback and "No valid docking-engine outputs were detected" in str(exc):
                selected_engine = _prompt_cli_engine_fallback(project_dir)
                if not selected_engine:
                    print("Interactive fallback cancelled.")
                    sys.exit(1)
                pipeline_kwargs["engine"] = selected_engine
                pipeline_kwargs["favorite_engine"] = selected_engine
                engine_pipeline = UnifiedPostDockingPipeline(**pipeline_kwargs)
            else:
                raise
        success = engine_pipeline.run()
        sys.exit(0 if success else 1)
    
    # If preprocessing is requested, run the preprocessing script
    if args.preprocess:
        from .preprocess import preprocess_analysis
        success = preprocess_analysis(
            args.input or ".",
            args.output,
            args.pairlist,
            args.force
        )
        sys.exit(0 if success else 1)
    
    parser.error(
        "Legacy config-driven CLI execution is disabled because it routes to deprecated stub analyzers.\n"
        "Use either:\n"
        "  1) python -m post_docking_analysis --project-dir <project_root> [--analysis-mode ...]\n"
        "  2) python -m post_docking_analysis.simplified_cli --project-dir <gnina_root> --output <out_dir>\n"
        "  3) python -m post_docking_analysis --preprocess -i <input_dir> [-o <output_dir>]\n"
    )

if __name__ == "__main__":
    main()
