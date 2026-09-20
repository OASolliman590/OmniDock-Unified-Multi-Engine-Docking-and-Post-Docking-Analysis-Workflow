from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

from .cli_parser import (
    ALL_ANALYSIS_TARGETS,
    DEFAULT_CLEAN_INTERACTION_DATASET_ROOT,
    JUMPABLE_ANALYSIS_TARGETS,
    _parse_rmsd_workers,
    build_parser,
)


def _has_tty() -> bool:
    return sys.stdin.isatty() and sys.stdout.isatty()


def _tty_required(command_label: str) -> int:
    print(f"❌ {command_label} requires an interactive terminal (TTY).")
    print("   Use `python main.py --help` to see the non-interactive command surfaces.")
    return 2


def _build_dock_argv(args: argparse.Namespace, forced_engine: Optional[str] = None, force_dry_run: bool = False) -> List[str]:
    argv: List[str] = ["--project-dir", args.project_dir]
    engines = forced_engine or args.engines
    if engines:
        argv.extend(["--engines", engines])
    if args.mode:
        argv.extend(["--mode", args.mode])
    if args.favorite_engine:
        argv.extend(["--favorite-engine", args.favorite_engine])
    if args.dry_run or force_dry_run:
        argv.append("--dry-run")
    if args.skip_completed:
        argv.append("--skip-completed")
    if args.no_ligand_qc_gate:
        argv.append("--no-ligand-qc-gate")
    if args.no_receptor_qc_gate:
        argv.append("--no-receptor-qc-gate")
    for flag in [
        "--gnina-image",
        "--gnina-binary",
        "--gnina-device",
        "--gnina-cnn-scoring",
        "--vina-conda-env",
        "--vina-binary",
        "--smina-conda-env",
        "--smina-binary",
        "--smina-scoring",
        "--autodock4-binary",
        "--autogrid4-binary",
        "--autodock4-parameter-file",
        "--parameter-mode",
        "--parameter-preset",
        "--execution-environment",
        "--execution-workdir",
        "--prerequisites-dir",
        "--required-files-dir",
        "--shared-conda-env",
        "--container-image",
        "--receptor-qc-min-atom-count",
        "--receptor-qc-min-heavy-atom-count",
        "--receptor-qc-min-chain-count",
        "--receptor-qc-max-coordinate-span",
    ]:
        value = getattr(args, flag[2:].replace("-", "_"))
        if value:
            argv.extend([flag, str(value)])
    for flag in [
        "--gnina-cpu",
        "--vina-cpu",
        "--smina-cpu",
        "--autodock4-spacing",
        "--autodock4-ga-pop-size",
        "--autodock4-ga-num-evals",
        "--autodock4-ga-num-generations",
        "--autodock4-ga-run",
        "--autodock4-ls-search-freq",
        "--autodock4-torsdof",
        "--exhaustiveness",
        "--num-modes",
        "--seed",
        "--box-scale",
        "--box-padding",
    ]:
        value = getattr(args, flag[2:].replace("-", "_"))
        if value is not None:
            argv.extend([flag, str(value)])
    return argv


def _build_deploy_argv(args: argparse.Namespace) -> List[str]:
    argv: List[str] = ["--project-dir", args.project_dir, "--mode", args.mode]
    if args.engines:
        argv.extend(["--engines", args.engines])
    if args.round:
        argv.extend(["--round", args.round])
    if args.from_rerun_manifest:
        argv.extend(["--from-rerun-manifest", args.from_rerun_manifest])
    if args.allow_full_exhaustive:
        argv.append("--allow-full-exhaustive")
    if args.skip_completed:
        argv.append("--skip-completed")
    if args.no_ligand_qc_gate:
        argv.append("--no-ligand-qc-gate")
    if args.no_ligand_admet_filters:
        argv.append("--no-ligand-admet-filters")
    if args.no_receptor_qc_gate:
        argv.append("--no-receptor-qc-gate")
    for flag in [
        "--gnina-image",
        "--gnina-binary",
        "--gnina-device",
        "--gnina-cnn-scoring",
        "--vina-conda-env",
        "--vina-binary",
        "--smina-conda-env",
        "--smina-binary",
        "--smina-scoring",
        "--autodock4-binary",
        "--autogrid4-binary",
        "--autodock4-parameter-file",
        "--slurm-time",
        "--slurm-mem",
        "--slurm-partition",
        "--slurm-account",
        "--slurm-job-name-prefix",
        "--slurm-extra-args",
        "--hpc-profile",
        "--hpc-profile-file",
        "--remote-project-dir",
        "--ligand-qc-max-lipinski-violations",
        "--ligand-qc-max-molecular-weight",
        "--ligand-qc-max-logp",
        "--ligand-qc-max-tpsa",
        "--ligand-qc-max-rotatable-bonds",
        "--ligand-qc-max-formal-charge-abs",
        "--ligand-qc-min-heavy-atom-count",
    ]:
        value = getattr(args, flag[2:].replace("-", "_"))
        if value:
            argv.extend([flag, str(value)])
    for flag in [
        "--ligand-qc-allow-pains",
        "--ligand-qc-allow-brenk",
        "--ligand-qc-allow-reactive",
    ]:
        if getattr(args, flag[2:].replace("-", "_")):
            argv.append(flag)
    for flag in [
        "--gnina-cpu",
        "--vina-cpu",
        "--smina-cpu",
        "--autodock4-spacing",
        "--autodock4-ga-pop-size",
        "--autodock4-ga-num-evals",
        "--autodock4-ga-num-generations",
        "--autodock4-ga-run",
        "--autodock4-ls-search-freq",
        "--autodock4-torsdof",
        "--exhaustiveness",
        "--num-modes",
        "--slurm-cpus-per-task",
        "--slurm-gpus",
        "--slurm-array-parallelism",
        "--receptor-qc-min-atom-count",
        "--receptor-qc-min-heavy-atom-count",
        "--receptor-qc-min-chain-count",
        "--receptor-qc-max-coordinate-span",
    ]:
        value = getattr(args, flag[2:].replace("-", "_"))
        if value is not None:
            argv.extend([flag, str(value)])
    return argv


def _build_sync_argv(args: argparse.Namespace) -> List[str]:
    argv: List[str] = ["--project-dir", args.project_dir]
    for flag in [
        "--hpc-profile",
        "--hpc-profile-file",
        "--ssh-target",
        "--remote-project-dir",
    ]:
        value = getattr(args, flag[2:].replace("-", "_"))
        if value:
            argv.extend([flag, str(value)])
    if args.delete:
        argv.append("--delete")
    if args.dry_run:
        argv.append("--dry-run")
    return argv


def _build_submit_argv(args: argparse.Namespace) -> List[str]:
    argv: List[str] = ["--project-dir", args.project_dir]
    for flag in [
        "--engines",
        "--round",
        "--mode",
        "--hpc-profile",
        "--hpc-profile-file",
        "--ssh-target",
        "--remote-project-dir",
    ]:
        value = getattr(args, flag[2:].replace("-", "_"))
        if value:
            argv.extend([flag, str(value)])
    if args.dry_run:
        argv.append("--dry-run")
    if args.allow_full_exhaustive:
        argv.append("--allow-full-exhaustive")
    if args.no_ligand_qc_gate:
        argv.append("--no-ligand-qc-gate")
    if args.no_ligand_admet_filters:
        argv.append("--no-ligand-admet-filters")
    if args.no_receptor_qc_gate:
        argv.append("--no-receptor-qc-gate")
    for flag in [
        "--ligand-qc-max-lipinski-violations",
        "--ligand-qc-max-molecular-weight",
        "--ligand-qc-max-logp",
        "--ligand-qc-max-tpsa",
        "--ligand-qc-max-rotatable-bonds",
        "--ligand-qc-max-formal-charge-abs",
        "--ligand-qc-min-heavy-atom-count",
        "--receptor-qc-min-atom-count",
        "--receptor-qc-min-heavy-atom-count",
        "--receptor-qc-min-chain-count",
        "--receptor-qc-max-coordinate-span",
    ]:
        value = getattr(args, flag[2:].replace("-", "_"))
        if value is not None:
            argv.extend([flag, str(value)])
    for flag in [
        "--ligand-qc-allow-pains",
        "--ligand-qc-allow-brenk",
        "--ligand-qc-allow-reactive",
    ]:
        if getattr(args, flag[2:].replace("-", "_")):
            argv.append(flag)
    return argv


def _run_analysis_dispatch(args: argparse.Namespace, target: str) -> int:
    from .execution import run_analysis_comparative, run_analysis_favorite, run_analysis_target

    try:
        rmsd_workers = _parse_rmsd_workers(getattr(args, "rmsd_workers", 0))
    except Exception as exc:
        raise SystemExit(f"Invalid --rmsd-workers value: {exc}") from exc

    if target == "analyze.comparative":
        result = run_analysis_comparative(
            args.project_dir,
            args.output,
            config_file=getattr(args, "config_file", None),
            promote_exhaustive=bool(getattr(args, "promote_exhaustive", False)),
            rerun_engine=getattr(args, "rerun_engine", None),
            top_per_protein=int(getattr(args, "top_per_protein", 1) or 1),
            winner_only=bool(getattr(args, "winner_only", False)),
            min_affinity_advantage=float(getattr(args, "min_affinity_advantage", 0.0) or 0.0),
            max_rerun_pairs=int(getattr(args, "max_rerun_pairs", 0) or 0),
            pair_allowlist=getattr(args, "pair_allowlist", None),
            consensus_mode=str(getattr(args, "consensus_mode", "dockbox_geometric") or "dockbox_geometric"),
            rescoring_scope=str(getattr(args, "rescoring_scope", "top_n_per_protein") or "top_n_per_protein"),
            rescoring_top_n=int(getattr(args, "rescoring_top_n", 3) or 3),
            prompt_protein_names=bool(getattr(args, "prompt_protein_names", False)),
            prompt_ligand_names=bool(getattr(args, "prompt_ligand_names", False)),
            complex_query=getattr(args, "complex_query", None),
            exclude_problematic_ligands=bool(getattr(args, "exclude_problematic_ligands", False)),
            positive_affinity_threshold=float(getattr(args, "positive_affinity_threshold", 0.0) or 0.0),
            minimum_pose_count=int(getattr(args, "minimum_pose_count", 1) or 1),
            rmsd_workers=rmsd_workers,
            analysis_scope=str(getattr(args, "analysis_scope", "full") or "full"),
            normalization_method=str(getattr(args, "normalization_method", "per_engine_rank") or "per_engine_rank"),
            biology_file=getattr(args, "biology_file", None),
            biology_mapping_mode=str(getattr(args, "biology_mapping_mode", "auto") or "auto"),
            hit_class_policy=str(getattr(args, "hit_class_policy", "target_percentile") or "target_percentile"),
            hit_class_strong_percentile=float(getattr(args, "hit_class_strong_percentile", 0.10) or 0.10),
            hit_class_moderate_percentile=float(getattr(args, "hit_class_moderate_percentile", 0.35) or 0.35),
            top_pose_selection_policy=str(getattr(args, "top_pose_policy", "best_affinity") or "best_affinity"),
            top_pose_global_aggregation=str(getattr(args, "top_pose_aggregation", "best_target") or "best_target"),
            best_pose_selection_metric=str(getattr(args, "best_pose_selection_metric", "auto") or "auto"),
        )
        if out_dir := result.outputs.get("output_dir"):
            print(f"Comparative analysis output: {out_dir}")
        if consensus_file := result.outputs.get("consensus_ranked_hits_file"):
            print(f"Consensus ranked hits: {consensus_file}")
        if rescoring_file := result.outputs.get("rescoring_candidates_file"):
            print(f"Rescoring candidates: {rescoring_file}")
        if explain_file := result.outputs.get("consensus_explainability_file"):
            print(f"Consensus explainability: {explain_file}")
        if corr_file := result.outputs.get("engine_rank_correlation_per_protein_file"):
            print(f"Per-protein rank correlations: {corr_file}")
        if corr_global_file := result.outputs.get("engine_rank_correlation_global_file"):
            print(f"Global rank correlations: {corr_global_file}")
        if canonical_scores := result.outputs.get("canonical_scores_root"):
            print(f"Canonical score outputs: {canonical_scores}")
        if class_file := result.outputs.get("consensus_ranked_hits_with_classes_file"):
            print(f"Hit classes file: {class_file}")
        if bio_map := result.outputs.get("biology_mapping_report_file"):
            print(f"Biology mapping report: {bio_map}")
        if bio_corr := result.outputs.get("biology_correlation_global_file"):
            print(f"Biology global correlations: {bio_corr}")
        if top_pose_global := result.outputs.get("top_pose_per_ligand_global_file"):
            print(f"Top-pose global atlas: {top_pose_global}")
        if top_pose_summary := result.outputs.get("top_pose_summary_file"):
            print(f"Top-pose summary: {top_pose_summary}")
        if top_pose_canonical := result.outputs.get("top_pose_canonical_root"):
            print(f"Top-pose canonical mirror: {top_pose_canonical}")
        if rerun_manifest := result.outputs.get("rerun_manifest_file"):
            print(f"Exhaustive rerun manifest: {rerun_manifest}")
        return 0 if result.status == "completed" else 1
    if target == "analyze.favorite_engine":
        if not args.favorite_engine:
            raise SystemExit("--favorite-engine is required for analyze.favorite_engine")
        result = run_analysis_favorite(
            args.project_dir,
            args.favorite_engine,
            args.output,
            config_file=getattr(args, "config_file", None),
            prompt_protein_names=bool(getattr(args, "prompt_protein_names", False)),
            prompt_ligand_names=bool(getattr(args, "prompt_ligand_names", False)),
            complex_query=getattr(args, "complex_query", None),
            rmsd_scopes=getattr(args, "rmsd_scopes", "per_complex"),
            resume_rmsd=not bool(getattr(args, "no_resume_rmsd", False)),
            exclude_problematic_ligands=bool(getattr(args, "exclude_problematic_ligands", False)),
            positive_affinity_threshold=float(getattr(args, "positive_affinity_threshold", 0.0) or 0.0),
            minimum_pose_count=int(getattr(args, "minimum_pose_count", 1) or 1),
            rmsd_workers=rmsd_workers,
            speed_profile=str(getattr(args, "speed_profile", "standard") or "standard"),
            analysis_scope=str(getattr(args, "analysis_scope", "full") or "full"),
            normalization_method=str(getattr(args, "normalization_method", "per_engine_rank") or "per_engine_rank"),
            biology_file=getattr(args, "biology_file", None),
            biology_mapping_mode=str(getattr(args, "biology_mapping_mode", "auto") or "auto"),
            hit_class_policy=str(getattr(args, "hit_class_policy", "target_percentile") or "target_percentile"),
            hit_class_strong_percentile=float(getattr(args, "hit_class_strong_percentile", 0.10) or 0.10),
            hit_class_moderate_percentile=float(getattr(args, "hit_class_moderate_percentile", 0.35) or 0.35),
            top_pose_selection_policy=str(getattr(args, "top_pose_policy", "best_affinity") or "best_affinity"),
            top_pose_global_aggregation=str(getattr(args, "top_pose_aggregation", "best_target") or "best_target"),
            best_pose_selection_metric=str(getattr(args, "best_pose_selection_metric", "auto") or "auto"),
        )
        if top_pose_global := result.outputs.get("top_pose_per_ligand_global_file"):
            print(f"Top-pose global atlas: {top_pose_global}")
        if top_pose_summary := result.outputs.get("top_pose_summary_file"):
            print(f"Top-pose summary: {top_pose_summary}")
        if top_pose_canonical := result.outputs.get("top_pose_canonical_root"):
            print(f"Top-pose canonical mirror: {top_pose_canonical}")
        return 0 if result.status == "completed" else 1
    result = run_analysis_target(
        target,
        project_dir=getattr(args, "project_dir", None),
        config_file=getattr(args, "config_file", None),
        output_dir=getattr(args, "output", None),
        engine=getattr(args, "engine", None),
        favorite_engine=getattr(args, "favorite_engine", None),
        sdf_folder=getattr(args, "sdf_folder", None),
        log_folder=getattr(args, "log_folder", None),
        receptors_folder=getattr(args, "receptors_folder", None),
        pairlist=getattr(args, "pairlist", None),
        ligplus_root=getattr(args, "ligplus_root", None),
        enable_poseview=bool(getattr(args, "enable_poseview", False)),
        prompt_protein_names=bool(getattr(args, "prompt_protein_names", False)),
        prompt_ligand_names=bool(getattr(args, "prompt_ligand_names", False)),
        complex_query=getattr(args, "complex_query", None),
        rmsd_scopes=getattr(args, "rmsd_scopes", "per_complex"),
        resume_rmsd=not bool(getattr(args, "no_resume_rmsd", False)),
        exclude_problematic_ligands=bool(getattr(args, "exclude_problematic_ligands", False)),
        positive_affinity_threshold=float(getattr(args, "positive_affinity_threshold", 0.0) or 0.0),
        minimum_pose_count=int(getattr(args, "minimum_pose_count", 1) or 1),
        rmsd_workers=rmsd_workers,
        speed_profile=str(getattr(args, "speed_profile", "standard") or "standard"),
        dataset_root=getattr(args, "dataset_root", None),
        python_bin=getattr(args, "python_bin", None),
        pipeline_execution_mode=str(getattr(args, "execution_mode", "local") or "local"),
        target_filter=getattr(args, "target_filter", None),
        ligand_filter=getattr(args, "ligand_filter", None),
        max_targets=int(getattr(args, "max_targets", 0) or 0),
        max_ligands=int(getattr(args, "max_ligands", 0) or 0),
        max_poses=int(getattr(args, "max_poses", 0) or 0),
        run_chemistry_fixes=not bool(getattr(args, "skip_chem_fixes", False)),
        run_layered_plip=not bool(getattr(args, "skip_layered_plip", False)),
        run_plip=not bool(getattr(args, "no_plip", False)),
        run_prolif=not bool(getattr(args, "no_prolif", False)),
        dry_run=bool(getattr(args, "dry_run", False)),
    )
    return 0 if result.status == "completed" else 1


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.group == "workflow":
        if args.workflow_command == "init":
            from .execution import run_workflow_init

            engines = [engine.strip() for engine in args.engines.split(",") if engine.strip()]
            result = run_workflow_init(
                args.project_dir,
                engines=engines,
                project_name=args.project_name or "",
                favorite_engine=args.favorite_engine or "",
                layout_profile=args.layout_profile,
            )
            print(f"Initialized workflow root: {result.root}")
            if layout_type := result.outputs.get("layout_type"):
                print(f"Layout type: {layout_type}")
            if registered_as := result.outputs.get("registered_as"):
                print(f"Registered as: {registered_as}")
            compatibility_profiles = result.outputs.get("compatibility_profiles", [])
            if compatibility_profiles:
                print(f"Compatibility profiles: {', '.join(str(p) for p in compatibility_profiles)}")
            return 0 if result.status == "completed" else 1
        if args.workflow_command == "interactive":
            if not _has_tty():
                return _tty_required("workflow interactive")
            from .interactive import run_interactive_workflow
            return run_interactive_workflow(args.project_dir or args.output)
        if args.workflow_command == "resume":
            if not _has_tty():
                return _tty_required("workflow resume")
            from .interactive import run_interactive_workflow
            return run_interactive_workflow(args.project_dir)
        if args.workflow_command in {"clone-checkpoint", "clone-maturation"}:
            from .execution import run_workflow_clone_checkpoint

            result = run_workflow_clone_checkpoint(
                source_project_dir=args.source_project_dir,
                target_project_dir=args.target_project_dir,
                layout_profile=getattr(args, "layout_profile", None),
            )
            print(f"Checkpoint workspace clone: {result.outputs.get('target_project_dir', '')}")
            print(f"Source project: {result.outputs.get('source_project_dir', '')}")
            if layout := result.outputs.get("layout_profile"):
                print(f"Layout profile: {layout}")
            return 0 if result.status == "completed" else 1
        if args.workflow_command == "status":
            from .state import summarize_state

            for line in summarize_state(Path(args.project_dir).expanduser().resolve()):
                print(line)
            return 0
        if args.workflow_command == "jump":
            return _run_analysis_dispatch(args, args.target)
        parser.error("A workflow subcommand is required")

    if args.group == "pdb":
        if args.pdb_command == "fetch":
            from .execution import run_pdb_fetch

            run_pdb_fetch(args.pdbs, args.output)
            return 0
        if args.pdb_command == "collect":
            if args.selection_mode == "interactive" and not _has_tty():
                return _tty_required("pdb collect --selection-mode interactive")
            from .execution import run_pdb_collect

            result = run_pdb_collect(
                args.project_dir,
                args.pdbs,
                selection_mode=args.selection_mode,
                preserve_metals=not bool(args.drop_metals),
                preserve_cofactors=not bool(args.drop_cofactors),
                preserve_full_receptor=bool(args.preserve_full_receptor),
                preserve_residues=_split_csv_option(args.preserve_residues),
                preferred_ligands=_split_csv_option(args.preferred_ligands),
            )
            return 0 if result.status in {"completed", "partial"} else 1
        if args.pdb_command == "run":
            from .execution import run_pdb_cli

            result = run_pdb_cli(args.pdbs, args.output, args.config)
            return 0 if result.status == "completed" else 1
        if args.pdb_command == "batch":
            from .execution import run_pdb_batch

            result = run_pdb_batch(args.config, args.output)
            return 0 if result.status in {"completed", "partial"} else 1
        if args.pdb_command == "legacy-interactive":
            if not _has_tty():
                return _tty_required("pdb legacy-interactive")
            from .execution import run_legacy_pdb_interactive

            run_legacy_pdb_interactive(args.output)
            return 0
        if args.pdb_command in {"prepare-protein", "prepare-ligand", "prepare-both"}:
            from .execution import run_autodock_prepare

            target = {
                "prepare-protein": "pdb.prepare_protein",
                "prepare-ligand": "pdb.prepare_ligand",
                "prepare-both": "pdb.prepare_both",
            }[args.pdb_command]
            result = run_autodock_prepare(
                target,
                args.receptors_input,
                args.ligands_input,
                args.receptors_output,
                args.ligands_output,
                force_field=args.force_field,
                ph=args.ph,
                ligand_preparation_backend=args.ligand_backend,
                selected_engines=[token.strip().lower() for token in str(getattr(args, "selected_engines", "") or "").split(",") if token.strip()],
                autodocktools_prepare_ligand4=args.autodocktools_prepare_ligand4,
                autodocktools_prepare_receptor4=args.autodocktools_prepare_receptor4,
                autodocktools_python=args.autodocktools_python,
            )
            return 0 if result.status == "completed" else 1
        parser.error("A pdb subcommand is required")

    if args.group == "prep":
        if args.prep_command == "pairlist":
            from .execution import run_prepare_pairlist

            if not args.mode and not args.freeze:
                raise SystemExit("--mode is required unless --freeze is used to materialize an existing round")
            curated_receptors = _split_csv_option(args.curated_receptors)
            curated_ligands = _split_csv_option(args.curated_ligands)
            curated_mapping = _parse_curated_mapping(args.curated_mapping)
            result = run_prepare_pairlist(
                project_dir=args.project_dir,
                mode=args.mode,
                prepared_proteins=args.prepared_proteins,
                prepared_ligands=args.prepared_ligands,
                excel_path=args.excel,
                default_site_id=args.default_site_id,
                default_box_size=args.default_box_size,
                curated_receptors=curated_receptors,
                curated_ligands=curated_ligands,
                curated_mapping=curated_mapping,
                iterative=bool(args.iterative),
                round_id=args.round,
                freeze=bool(args.freeze),
                prompt_protein_aliases=bool(args.prompt_protein_aliases),
                prompt_ligand_aliases=bool(args.prompt_ligand_aliases),
            )
            print(f"Pairlist mode: {result.outputs.get('pair_mode', args.mode or '')}")
            if round_id := result.outputs.get("round_id"):
                print(f"Pair-curation round: {round_id}")
            if state_file := result.outputs.get("pair_curation_state_file"):
                print(f"Pair-curation state: {state_file}")
            if pairlist_file := result.outputs.get("pairlist_file"):
                print(f"Pairlist file: {pairlist_file}")
            if pair_intent_file := result.outputs.get("pair_intent_file"):
                print(f"Pair intent file: {pair_intent_file}")
            if protein_alias_file := result.outputs.get("protein_alias_file"):
                print(f"Protein alias file: {protein_alias_file}")
            if ligand_alias_file := result.outputs.get("ligand_alias_file"):
                print(f"Ligand alias file: {ligand_alias_file}")
            if pair_count := result.outputs.get("pair_count"):
                print(f"Pair count: {pair_count}")
            return 0 if result.status == "completed" else 1
        if args.prep_command == "project":
            from .execution import run_prepare_project

            _hydrate_project_prep_defaults(args)
            prep_args = _namespace_to_argv(args, skip={"group", "prep_command", "project_dir"})
            result = run_prepare_project(prep_args)
            return 0 if result.status == "completed" else 1
        parser.error("A prep subcommand is required")

    if args.group == "dock":
        if args.dock_command == "run":
            from .execution import run_docking

            result = run_docking(_build_dock_argv(args))
            return 0 if result.status == "completed" else 1
        if args.dock_command == "dry-run":
            from .execution import run_docking

            result = run_docking(_build_dock_argv(args, force_dry_run=True))
            return 0 if result.status == "completed" else 1
        if args.dock_command == "engine":
            from .execution import run_docking

            result = run_docking(_build_dock_argv(args, forced_engine=args.engine_name), engine_target=f"dock.engine.{args.engine_name}")
            return 0 if result.status == "completed" else 1
        if args.dock_command == "deploy":
            from .execution import run_docking_deploy

            result = run_docking_deploy(_build_deploy_argv(args))
            return 0 if result.status == "completed" else 1
        if args.dock_command == "sync":
            from .execution import run_docking_sync

            result = run_docking_sync(_build_sync_argv(args))
            return 0 if result.status in {"completed", "partial"} else 1
        if args.dock_command == "submit":
            from .execution import run_docking_submit

            result = run_docking_submit(_build_submit_argv(args))
            return 0 if result.status in {"completed", "partial"} else 1
        parser.error("A dock subcommand is required")

    if args.group == "analyze":
        if args.analyze_command == "comparative":
            return _run_analysis_dispatch(args, "analyze.comparative")
        if args.analyze_command == "favorite-engine":
            return _run_analysis_dispatch(args, "analyze.favorite_engine")
        if args.analyze_command == "clean":
            return _run_analysis_dispatch(args, "analyze.interactions.clean")
        if args.analyze_command == "stage":
            target = {
                "hierarchical": "analyze.stage.hierarchical",
                "polypharmacology": "analyze.stage.polypharmacology",
                "rmsd": "analyze.stage.rmsd",
                "reports": "analyze.stage.reports",
                "visualizations": "analyze.stage.visualizations",
                "structure-quality": "analyze.stage.structure_quality",
            }[args.stage_command]
            return _run_analysis_dispatch(args, target)
        if args.analyze_command == "interactions":
            target = f"analyze.interactions.{args.interactions_command}"
            return _run_analysis_dispatch(args, target)
        if args.analyze_command == "visuals":
            target = "analyze.visuals.py3dmol" if args.visuals_command == "py3dmol" else "analyze.visuals.pymol"
            return _run_analysis_dispatch(args, target)
        parser.error("An analyze subcommand is required")

    parser.print_help()
    return 0


def _namespace_to_argv(args: argparse.Namespace, skip: Optional[set[str]] = None) -> List[str]:
    skip = skip or set()
    argv: List[str] = []
    for key, value in vars(args).items():
        if key in skip or value in (None, False):
            continue
        flag = f"--{key.replace('_', '-')}"
        if value is True:
            argv.append(flag)
        else:
            argv.extend([flag, str(value)])
    return argv


def _split_csv_option(value: Optional[str]) -> Optional[List[str]]:
    if not value:
        return None
    return [item.strip() for item in value.split(",") if item.strip()]


def _parse_curated_mapping(value: Optional[str]) -> Optional[dict[str, list[str]]]:
    if not value:
        return None
    mapping: dict[str, list[str]] = {}
    for chunk in value.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        receptor, _, ligands = chunk.partition("=")
        mapping[receptor.strip()] = [item.strip() for item in ligands.split(",") if item.strip()]
    return mapping or None


def _hydrate_project_prep_defaults(args: argparse.Namespace) -> None:
    if not args.project_dir:
        required = {
            "prepared_proteins": args.prepared_proteins,
            "prepared_ligands": args.prepared_ligands,
            "output": args.output,
        }
        if not args.pairlist_file:
            required["excel"] = args.excel
        missing = [name for name, value in required.items() if not value]
        if missing:
            raise SystemExit(
                f"--{' --'.join(missing)} is required when --project-dir is not provided"
            )
        return

    from docking.project_layout import ensure_project_layout, load_manifest, pair_intent_path, pairlist_path

    project_root = Path(args.project_dir).expanduser().resolve()
    layout = ensure_project_layout(project_root, args.layout_profile)
    manifest = load_manifest(project_root)
    args.output = args.output or str(project_root)
    args.prepared_proteins = args.prepared_proteins or manifest.get("prepared_proteins_dir") or str(layout["prepared_proteins"])
    args.prepared_ligands = args.prepared_ligands or manifest.get("prepared_ligands_dir") or str(layout["prepared_ligands"])
    args.raw_proteins = args.raw_proteins or manifest.get("raw_proteins_dir") or str(layout["raw_proteins"])
    args.raw_ligands = args.raw_ligands or manifest.get("raw_ligands_dir") or str(layout["raw_ligands"])
    default_pair_intent = pair_intent_path(project_root, manifest.get("layout_profile", args.layout_profile))
    default_pairlist = pairlist_path(project_root, manifest.get("layout_profile", args.layout_profile))
    if not args.pair_intent and default_pair_intent.exists():
        args.pair_intent = str(default_pair_intent)
    if not args.pairlist_file and default_pairlist.exists():
        args.pairlist_file = str(default_pairlist)
    excel_candidate = args.excel or manifest.get("source_paths", {}).get("excel_path") or str(layout["raw_proteins"] / "multi_pdb_analysis.xlsx")
    if excel_candidate:
        excel_path = Path(str(excel_candidate)).expanduser()
        if args.excel:
            args.excel = str(excel_path)
        elif excel_path.exists():
            args.excel = str(excel_path)
        elif args.pairlist_file:
            args.excel = None
        else:
            args.excel = str(excel_path)
