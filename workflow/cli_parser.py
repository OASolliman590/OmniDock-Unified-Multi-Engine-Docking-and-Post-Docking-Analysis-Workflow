from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from docking.models import MAX_PREPARATION_PH, MIN_PREPARATION_PH


ALL_ANALYSIS_TARGETS = [
    "analyze.stage.hierarchical",
    "analyze.stage.polypharmacology",
    "analyze.stage.rmsd",
    "analyze.stage.reports",
    "analyze.stage.visualizations",
    "analyze.stage.structure_quality",
    "analyze.interactions.pandamap",
    "analyze.interactions.prolif",
    "analyze.interactions.ligplot",
    "analyze.interactions.poseview",
    "analyze.interactions.clean",
    "analyze.visuals.py3dmol",
    "analyze.visuals.pymol",
]

DEFAULT_CLEAN_INTERACTION_DATASET_ROOT = str(
    Path(
        os.environ.get(
            "DOCKFORGE_CLEAN_INTERACTION_DATASET_ROOT",
            str(Path.cwd()),
        )
    )
    .expanduser()
    .resolve()
)


JUMPABLE_ANALYSIS_TARGETS = [
    "analyze.comparative",
    "analyze.favorite_engine",
    *ALL_ANALYSIS_TARGETS,
]


def _add_rmsd_args(parser: argparse.ArgumentParser) -> None:
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


def _add_analysis_common_args(parser: argparse.ArgumentParser, require_project: bool = False) -> None:
    parser.add_argument("--project-dir", required=require_project, help="Canonical project root or GNINA project root")
    parser.add_argument("--config-file", help="Optional post-docking analysis YAML/JSON config file")
    parser.add_argument("-o", "--output", help="Output directory")
    parser.add_argument("--engine", help="Engine for stage-level analysis in canonical projects")
    parser.add_argument("--favorite-engine", help="Favorite engine for favorite-engine continuation")
    parser.add_argument("--sdf-folder", help="Explicit GNINA SDF folder")
    parser.add_argument("--log-folder", help="Explicit GNINA log folder")
    parser.add_argument("--receptors-folder", help="Explicit GNINA receptors folder")
    parser.add_argument("--pairlist", help="Optional pairlist.csv")
    parser.add_argument("--ligplus-root", help="Optional LigPlus root")
    parser.add_argument(
        "--complex-query",
        help="Optional complex filter query, for example: protein=2FVD;ligand=Sorafenib or exclude_ligand=Acetazolamide. Use 0/skip/all to disable.",
    )
    parser.add_argument("--enable-poseview", action="store_true", help="Enable PoseView stage")
    parser.add_argument("--prompt-protein-names", action="store_true", help="Prompt for protein display names before analysis starts")
    parser.add_argument("--prompt-ligand-names", action="store_true", help="Prompt for ligand display names before analysis starts")
    parser.add_argument(
        "--exclude-problematic-ligands",
        action="store_true",
        help="Exclude ligands with positive/weak affinity or missing pose support before analysis",
    )
    parser.add_argument(
        "--positive-affinity-threshold",
        type=float,
        default=0.0,
        help="Exclude ligands when best affinity is greater than this threshold (kcal/mol)",
    )
    parser.add_argument(
        "--minimum-pose-count",
        type=int,
        default=1,
        help="Minimum required poses per ligand when problematic-ligand filtering is enabled",
    )
    parser.add_argument(
        "--rmsd-workers",
        default="0",
        help="RMSD workers (0=auto, 1=single-thread, >1 parallel; accepts auto/single/parallel)",
    )
    parser.add_argument(
        "--speed-profile",
        choices=["standard", "fast"],
        default="standard",
        help="Analysis speed profile (fast skips optional heavy stages where supported)",
    )
    _add_rmsd_args(parser)


def _add_clean_interaction_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--project-dir", help="Optional workflow project root for state tracking")
    parser.add_argument("--config-file", help="Optional post-docking analysis YAML/JSON config file")
    parser.add_argument(
        "--dataset-root",
        default=DEFAULT_CLEAN_INTERACTION_DATASET_ROOT,
        help="Dataset root containing pose_ensemble_appraisal and chemistry-fix scripts",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output directory for clean-pipeline command manifests/logs (defaults under dataset root)",
    )
    parser.add_argument(
        "--python-bin",
        default=sys.executable,
        help="Python executable used to run external scripts/modules",
    )
    parser.add_argument(
        "--execution-mode",
        choices=["local", "local_dry_run", "hpc"],
        default="local",
        help="pose_ensemble_appraisal execution mode",
    )
    parser.add_argument("--target-filter", help="Optional target filter passed to pose_ensemble_appraisal")
    parser.add_argument("--ligand-filter", help="Optional ligand filter passed to pose_ensemble_appraisal")
    parser.add_argument(
        "--max-targets",
        type=int,
        default=3,
        help="Safe default cap for targets (<=0 disables cap)",
    )
    parser.add_argument(
        "--max-ligands",
        type=int,
        default=0,
        help="Optional cap for ligands (<=0 disables cap)",
    )
    parser.add_argument(
        "--max-poses",
        type=int,
        default=3,
        help="Safe default cap for poses (<=0 disables cap)",
    )
    parser.add_argument(
        "--skip-chem-fixes",
        action="store_true",
        help="Skip chemistry/format-fix scripts and run interaction extraction only",
    )
    parser.add_argument(
        "--skip-layered-plip",
        action="store_true",
        help="Skip layered PLIP extraction wrapper stage",
    )
    parser.add_argument(
        "--no-plip",
        action="store_true",
        help="Disable PLIP stage in pose_ensemble_appraisal",
    )
    parser.add_argument(
        "--no-prolif",
        action="store_true",
        help="Disable ProLIF stage in pose_ensemble_appraisal",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Plan command set without executing external scripts",
    )


def _add_workflow_group(subparsers) -> None:
    workflow = subparsers.add_parser("workflow", help="Workflow shell and state utilities")
    workflow_sub = workflow.add_subparsers(dest="workflow_command")
    workflow_interactive = workflow_sub.add_parser("interactive", help="Launch the guided workflow shell")
    workflow_interactive.add_argument("-o", "--output", help="Initial workflow root")
    workflow_interactive.add_argument("--project-dir", help="Existing project root to resume immediately")
    workflow_init = workflow_sub.add_parser("init", help="Initialize a project directory for workflow use")
    workflow_init.add_argument("--project-dir", required=True, help="Project root to initialize")
    workflow_init.add_argument("--engines", default="gnina,vina,smina,autodock4", help="Comma-separated engine list for canonical projects")
    workflow_init.add_argument("--project-name", help="Optional project name for the manifest")
    workflow_init.add_argument("--favorite-engine", help="Optional favorite engine to record")
    workflow_init.add_argument("--layout-profile", choices=["canonical", "docking_legacy"], default="docking_legacy")
    workflow_status = workflow_sub.add_parser("status", help="Show workflow state summary")
    workflow_status.add_argument("--project-dir", required=True, help="Workflow/project root")
    workflow_resume = workflow_sub.add_parser("resume", help="Resume using the interactive shell")
    workflow_resume.add_argument("--project-dir", required=True, help="Workflow/project root")
    workflow_clone_checkpoint = workflow_sub.add_parser(
        "clone-checkpoint",
        help="Clone a project into a Checkpoint & Revise workspace",
    )
    workflow_clone_checkpoint.add_argument("--source-project-dir", required=True, help="Source project root to clone")
    workflow_clone_checkpoint.add_argument(
        "--target-project-dir",
        help="Target checkpoint workspace root (default: sibling <source>_checkpoint)",
    )
    workflow_clone_checkpoint.add_argument(
        "--layout-profile",
        choices=["canonical", "docking_legacy"],
        help="Optional layout profile override for re-initializing the clone",
    )
    # Backward-compatible alias
    workflow_clone_legacy = workflow_sub.add_parser(
        "clone-maturation",
        help="(deprecated) Alias for clone-checkpoint",
    )
    workflow_clone_legacy.add_argument("--source-project-dir", required=True, help="Source project root to clone")
    workflow_clone_legacy.add_argument(
        "--target-project-dir",
        help="Target checkpoint workspace root (default: sibling <source>_checkpoint)",
    )
    workflow_clone_legacy.add_argument(
        "--layout-profile",
        choices=["canonical", "docking_legacy"],
        help="Optional layout profile override for re-initializing the clone",
    )
    workflow_jump = workflow_sub.add_parser("jump", help="Jump directly to an analysis target")
    workflow_jump.add_argument("--target", required=True, choices=JUMPABLE_ANALYSIS_TARGETS, help="Analysis target to execute")
    _add_analysis_common_args(workflow_jump, require_project=True)
    workflow_jump.add_argument("--promote-exhaustive", action="store_true", help="Emit an exhaustive rerun manifest from comparative results")
    workflow_jump.add_argument("--rerun-engine", help="Engine to promote into exhaustive reruns")
    workflow_jump.add_argument("--top-per-protein", type=int, default=1, help="Number of best ligands per protein to promote into the rerun manifest")
    workflow_jump.add_argument("--winner-only", action="store_true", help="Only promote pairs where the target engine wins against the other engines")
    workflow_jump.add_argument("--min-affinity-advantage", type=float, default=0.0, help="Require the target engine to beat the best competing engine by at least this many kcal/mol")
    workflow_jump.add_argument("--max-rerun-pairs", type=int, default=0, help="Optional global cap on promoted exhaustive rerun pairs")
    workflow_jump.add_argument("--pair-allowlist", help="Optional TXT/CSV file restricting promotion to explicit pair tags or receptor/site_id/ligand rows")
    workflow_jump.add_argument(
        "--consensus-mode",
        choices=["dockbox_geometric", "weighted_hybrid", "strict_consensus", "favorite_guardrails"],
        default="dockbox_geometric",
        help="Consensus policy for comparative hit ranking and rerun promotion",
    )
    workflow_jump.add_argument(
        "--rescoring-scope",
        choices=["top_n_per_protein", "top_n_global"],
        default="top_n_per_protein",
        help="Scope used to select rescoring candidate shortlist",
    )
    workflow_jump.add_argument(
        "--rescoring-top-n",
        type=int,
        default=3,
        help="Top-N candidates retained by the selected rescoring scope",
    )


def _add_pdb_group(subparsers) -> None:
    pdb = subparsers.add_parser("pdb", help="PDB preparation surfaces")
    pdb_sub = pdb.add_subparsers(dest="pdb_command")
    pdb_fetch = pdb_sub.add_parser("fetch", help="Fetch/download PDB structures")
    pdb_fetch.add_argument("-p", "--pdbs", required=True, help="PDB ID(s) or file path")
    pdb_fetch.add_argument("-o", "--output", required=True, help="Output directory")
    pdb_collect = pdb_sub.add_parser("collect", help="Fetch PDBs into a staged docking project")
    pdb_collect.add_argument("--project-dir", required=True, help="Docking project root")
    pdb_collect.add_argument("-p", "--pdbs", required=True, help="PDB ID(s) or file path")
    pdb_collect.add_argument(
        "--selection-mode",
        choices=["auto", "heuristic", "interactive"],
        default="heuristic",
        help="Ligand selection strategy during PDB collection",
    )
    pdb_collect.add_argument(
        "--preserve-residues",
        help="Comma-separated residues to keep during default cleaning (e.g. ZN,MN,NAD,NAP,NDP)",
    )
    pdb_collect.add_argument(
        "--preferred-ligands",
        help="Comma-separated preferred ligand residue names for heuristic mode",
    )
    pdb_collect.add_argument(
        "--drop-metals",
        action="store_true",
        help="Remove common metals during default cleaning instead of preserving them",
    )
    pdb_collect.add_argument(
        "--drop-cofactors",
        action="store_true",
        help="Remove common cofactors during default cleaning instead of preserving them",
    )
    pdb_collect.add_argument(
        "--preserve-full-receptor",
        action="store_true",
        help="Skip cleaning entirely and keep full receptor unchanged",
    )
    pdb_run = pdb_sub.add_parser("run", help="Run CLI-style PDB preparation")
    pdb_run.add_argument("-p", "--pdbs", required=True, help="PDB ID(s) or file path")
    pdb_run.add_argument("-o", "--output", required=True, help="Output directory")
    pdb_run.add_argument("-c", "--config", help="Optional JSON config")
    pdb_batch = pdb_sub.add_parser("batch", help="Run batch PDB preparation")
    pdb_batch.add_argument("-c", "--config", required=True, help="Batch config file")
    pdb_batch.add_argument("-o", "--output", default="batch_docking_preparation", help="Output directory")
    pdb_legacy = pdb_sub.add_parser("legacy-interactive", help="Run the legacy interactive PDB flow")
    pdb_legacy.add_argument("-o", "--output", help="Output directory")
    for command_name, help_text in [
        ("prepare-protein", "Prepare receptors/proteins for docking"),
        ("prepare-ligand", "Prepare ligands for docking"),
        ("prepare-both", "Prepare receptors and ligands for docking"),
    ]:
        command = pdb_sub.add_parser(command_name, help=help_text)
        command.add_argument("--receptors-input", help="Raw receptors input directory")
        command.add_argument("--ligands-input", help="Raw ligands input directory")
        command.add_argument("--receptors-output", help="Prepared receptors output directory")
        command.add_argument("--ligands-output", help="Prepared ligands output directory")
        command.add_argument("--force-field", default="AMBER", help="Force field")
        command.add_argument(
            "--ph",
            type=float,
            default=7.4,
            help=f"Protonation pH ({MIN_PREPARATION_PH:.1f}-{MAX_PREPARATION_PH:.1f})",
        )
        command.add_argument(
            "--ligand-backend",
            choices=[
                "engine_aware_full",
                "openbabel_only",
                "meeko_only",
                "autodocktools_only",
                "openbabel_meeko",
                "openbabel_meeko_autodock",
                "openbabel_autodocktools",
            ],
            default="engine_aware_full",
            help="Ligand preparation profile",
        )
        command.add_argument(
            "--selected-engines",
            default="",
            help="Optional comma-separated selected engines for engine-aware profile (example: gnina,vina,smina,autodock4)",
        )
        command.add_argument("--autodocktools-prepare-ligand4", help="Optional path to prepare_ligand4.py")
        command.add_argument("--autodocktools-prepare-receptor4", help="Optional path to prepare_receptor4.py")
        command.add_argument("--autodocktools-python", help="Optional python executable for prepare_ligand4.py")


def _add_prep_group(subparsers) -> None:
    prep = subparsers.add_parser("prep", help="Docking project preparation")
    prep_sub = prep.add_subparsers(dest="prep_command")
    prep_pairlist = prep_sub.add_parser("pairlist", help="Build pair intent and pairlist.csv from prepared project assets")
    prep_pairlist.add_argument("--project-dir", required=True, help="Docking project root")
    prep_pairlist.add_argument(
        "--mode",
        choices=["cocrystal_only", "cocrystal_plus_all", "cocrystal_plus_nonreference", "curated_cartesian", "curated_per_protein"],
        help="Pairlist generation mode",
    )
    prep_pairlist.add_argument("--prepared-proteins", help="Prepared proteins directory override")
    prep_pairlist.add_argument("--prepared-ligands", help="Prepared ligands directory override")
    prep_pairlist.add_argument("--excel", help="multi_pdb_analysis.xlsx override")
    prep_pairlist.add_argument("--default-site-id", default="site_1")
    prep_pairlist.add_argument("--default-box-size", type=float, default=20.0)
    prep_pairlist.add_argument("--curated-receptors", help="Comma-separated receptor selection for curated_cartesian")
    prep_pairlist.add_argument("--curated-ligands", help="Comma-separated ligand selection for curated_cartesian")
    prep_pairlist.add_argument(
        "--curated-mapping",
        help="Semicolon-separated receptor=lig1,lig2 mappings for curated_per_protein",
    )
    prep_pairlist.add_argument("--iterative", action="store_true", help="Persist the selection as a reusable pair-curation round")
    prep_pairlist.add_argument("--round", help="Optional pair-curation round identifier")
    prep_pairlist.add_argument("--freeze", action="store_true", help="Materialize the selected or existing round into pairlist.csv and pair_intent.csv")
    prep_pairlist.add_argument("--prompt-protein-aliases", action="store_true", help="Prompt for persistent protein display aliases before writing pair metadata")
    prep_pairlist.add_argument("--prompt-ligand-aliases", action="store_true", help="Prompt for persistent ligand display aliases before writing pair metadata")
    prep_project = prep_sub.add_parser("project", help="Build a canonical docking project")
    prep_project.add_argument("--project-dir", help="Existing staged docking project root")
    prep_project.add_argument("--prepared-proteins", help="Prepared proteins directory")
    prep_project.add_argument("--prepared-ligands", help="Prepared ligands directory")
    prep_project.add_argument("--pair-intent", help="Manual pair-intent CSV")
    prep_project.add_argument("--pairlist-file", help="Prebuilt pairlist.csv to materialize")
    prep_project.add_argument("--excel", help="multi_pdb_analysis.xlsx path")
    prep_project.add_argument("--output", help="Output project directory")
    prep_project.add_argument("--raw-proteins", help="Optional raw proteins directory")
    prep_project.add_argument("--raw-ligands", help="Optional raw ligands directory")
    prep_project.add_argument("--default-site-id", default="site_1")
    prep_project.add_argument("--default-box-size", type=float, default=20.0)
    prep_project.add_argument("--asset-mode", choices=["symlink", "copy"], default="symlink")
    prep_project.add_argument("--engines", default="gnina,vina,smina,autodock4")
    prep_project.add_argument("--project-name")
    prep_project.add_argument("--layout-profile", choices=["canonical", "docking_legacy"], default="docking_legacy")
    prep_project.add_argument(
        "--pair-mode",
        choices=["manual", "protein-based", "cocrystal_only", "cocrystal_plus_all", "curated_cartesian", "curated_per_protein"],
        default="manual",
    )
    prep_project.add_argument("--dock-all-proteins-to-all-ligands", action="store_true")


def _add_dock_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--project-dir", required=True, help="Docking project root")
    parser.add_argument("--engines", help="Comma-separated engine list")
    parser.add_argument("--mode", choices=["standard", "score-only"], default="standard")
    parser.add_argument("--favorite-engine", help="Optional favorite engine to record")
    parser.add_argument("--dry-run", action="store_true", help="Dry-run without execution")
    parser.add_argument("--skip-completed", action="store_true", help="Skip completed rows")
    parser.add_argument("--gnina-image")
    parser.add_argument("--gnina-binary")
    parser.add_argument("--gnina-device")
    parser.add_argument("--gnina-cpu", type=int)
    parser.add_argument("--gnina-cnn-scoring", default="rescore")
    parser.add_argument("--vina-conda-env")
    parser.add_argument("--vina-binary")
    parser.add_argument("--vina-cpu", type=int)
    parser.add_argument("--smina-conda-env")
    parser.add_argument("--smina-binary")
    parser.add_argument("--smina-cpu", type=int)
    parser.add_argument("--smina-scoring")
    parser.add_argument("--autodock4-binary")
    parser.add_argument("--autogrid4-binary")
    parser.add_argument("--autodock4-parameter-file")
    parser.add_argument("--autodock4-spacing", type=float)
    parser.add_argument("--autodock4-ga-pop-size", type=int)
    parser.add_argument("--autodock4-ga-num-evals", type=int)
    parser.add_argument("--autodock4-ga-num-generations", type=int)
    parser.add_argument("--autodock4-ga-run", type=int)
    parser.add_argument("--autodock4-ls-search-freq", type=float)
    parser.add_argument("--autodock4-torsdof", type=int)
    parser.add_argument("--exhaustiveness", type=int, default=16)
    parser.add_argument("--num-modes", type=int, default=20)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--box-scale", type=float, default=1.0)
    parser.add_argument("--box-padding", type=float, default=0.0)
    parser.add_argument("--parameter-mode", choices=["basic", "advanced"], default="basic")
    parser.add_argument("--parameter-preset", choices=["screening_fast", "balanced", "exhaustive"], default="balanced")
    parser.add_argument(
        "--execution-environment",
        choices=["local_cpu", "local_gpu", "conda_env", "container", "remote_hpc"],
        default="local_cpu",
    )
    parser.add_argument("--execution-workdir")
    parser.add_argument("--prerequisites-dir")
    parser.add_argument("--required-files-dir")
    parser.add_argument("--shared-conda-env")
    parser.add_argument("--container-image")
    parser.add_argument("--no-ligand-qc-gate", action="store_true")
    parser.add_argument("--no-receptor-qc-gate", action="store_true")
    parser.add_argument("--receptor-qc-min-atom-count", type=int, default=100)
    parser.add_argument("--receptor-qc-min-heavy-atom-count", type=int, default=60)
    parser.add_argument("--receptor-qc-min-chain-count", type=int, default=1)
    parser.add_argument("--receptor-qc-max-coordinate-span", type=float, default=500.0)


def _parse_rmsd_workers(value: object) -> int:
    token = str(value or "").strip().lower()
    if not token:
        return 0
    if token in {"auto", "default"}:
        return 0
    if token in {"single", "single-thread", "single_thread"}:
        return 1
    if token in {"parallel", "multi", "multi-thread", "multi_thread", ">1", ">=2"}:
        return 2
    if token.startswith(">"):
        parsed = int(token[1:])
        return 2 if parsed <= 1 else parsed
    parsed = int(token)
    if parsed < 0:
        raise ValueError("RMSD worker count must be >= 0")
    return parsed


def _add_deploy_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--project-dir", required=True, help="Docking project root")
    parser.add_argument("--engines", help="Comma-separated engine list")
    parser.add_argument("--round", help="Optional logical round identifier")
    parser.add_argument("--mode", choices=["screen", "exhaustive"], default="screen")
    parser.add_argument("--from-rerun-manifest", help="Optional rerun manifest CSV from analyze comparative")
    parser.add_argument(
        "--allow-full-exhaustive",
        action="store_true",
        help="Allow exhaustive deployment without a rerun manifest. Use only for deliberate full-pair exhaustive campaigns.",
    )
    parser.add_argument("--skip-completed", action="store_true", help="Skip rows with existing pose output")
    parser.add_argument("--gnina-image")
    parser.add_argument("--gnina-binary")
    parser.add_argument("--gnina-device")
    parser.add_argument("--gnina-cpu", type=int)
    parser.add_argument("--gnina-cnn-scoring", default="rescore")
    parser.add_argument("--vina-conda-env")
    parser.add_argument("--vina-binary")
    parser.add_argument("--vina-cpu", type=int)
    parser.add_argument("--smina-conda-env")
    parser.add_argument("--smina-binary")
    parser.add_argument("--smina-cpu", type=int)
    parser.add_argument("--smina-scoring")
    parser.add_argument("--autodock4-binary")
    parser.add_argument("--autogrid4-binary")
    parser.add_argument("--autodock4-parameter-file")
    parser.add_argument("--autodock4-spacing", type=float)
    parser.add_argument("--autodock4-ga-pop-size", type=int)
    parser.add_argument("--autodock4-ga-num-evals", type=int)
    parser.add_argument("--autodock4-ga-num-generations", type=int)
    parser.add_argument("--autodock4-ga-run", type=int)
    parser.add_argument("--autodock4-ls-search-freq", type=float)
    parser.add_argument("--autodock4-torsdof", type=int)
    parser.add_argument("--exhaustiveness", type=int, default=16)
    parser.add_argument("--num-modes", type=int, default=20)
    parser.add_argument("--slurm-time")
    parser.add_argument("--slurm-mem")
    parser.add_argument("--slurm-cpus-per-task", type=int)
    parser.add_argument("--slurm-gpus", type=int)
    parser.add_argument("--slurm-partition")
    parser.add_argument("--slurm-account")
    parser.add_argument("--slurm-job-name-prefix")
    parser.add_argument("--slurm-array-parallelism", type=int)
    parser.add_argument("--slurm-extra-args")
    parser.add_argument("--hpc-profile")
    parser.add_argument("--hpc-profile-file")
    parser.add_argument("--remote-project-dir")
    parser.add_argument("--no-ligand-qc-gate", action="store_true")
    parser.add_argument("--no-ligand-admet-filters", action="store_true")
    parser.add_argument("--ligand-qc-max-lipinski-violations", type=int, default=1)
    parser.add_argument("--ligand-qc-max-molecular-weight", type=float, default=650.0)
    parser.add_argument("--ligand-qc-max-logp", type=float, default=6.0)
    parser.add_argument("--ligand-qc-max-tpsa", type=float, default=180.0)
    parser.add_argument("--ligand-qc-max-rotatable-bonds", type=int, default=15)
    parser.add_argument("--ligand-qc-max-formal-charge-abs", type=int, default=2)
    parser.add_argument("--ligand-qc-min-heavy-atom-count", type=int, default=6)
    parser.add_argument("--ligand-qc-allow-pains", action="store_true")
    parser.add_argument("--ligand-qc-allow-brenk", action="store_true")
    parser.add_argument("--ligand-qc-allow-reactive", action="store_true")
    parser.add_argument("--no-receptor-qc-gate", action="store_true")
    parser.add_argument("--receptor-qc-min-atom-count", type=int, default=100)
    parser.add_argument("--receptor-qc-min-heavy-atom-count", type=int, default=60)
    parser.add_argument("--receptor-qc-min-chain-count", type=int, default=1)
    parser.add_argument("--receptor-qc-max-coordinate-span", type=float, default=500.0)


def _add_sync_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--project-dir", required=True, help="Docking project root")
    parser.add_argument("--hpc-profile")
    parser.add_argument("--hpc-profile-file")
    parser.add_argument("--ssh-target")
    parser.add_argument("--remote-project-dir")
    parser.add_argument("--delete", action="store_true")
    parser.add_argument("--dry-run", action="store_true")


def _add_submit_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--project-dir", required=True, help="Docking project root")
    parser.add_argument("--engines", help="Comma-separated engine list")
    parser.add_argument("--round", help="Deployment round identifier")
    parser.add_argument("--mode", choices=["screen", "exhaustive"])
    parser.add_argument("--hpc-profile")
    parser.add_argument("--hpc-profile-file")
    parser.add_argument("--ssh-target")
    parser.add_argument("--remote-project-dir")
    parser.add_argument(
        "--allow-full-exhaustive",
        action="store_true",
        help="Allow exhaustive submission even when the bundle does not record a rerun manifest.",
    )
    parser.add_argument("--no-ligand-qc-gate", action="store_true")
    parser.add_argument("--no-ligand-admet-filters", action="store_true")
    parser.add_argument("--ligand-qc-max-lipinski-violations", type=int, default=1)
    parser.add_argument("--ligand-qc-max-molecular-weight", type=float, default=650.0)
    parser.add_argument("--ligand-qc-max-logp", type=float, default=6.0)
    parser.add_argument("--ligand-qc-max-tpsa", type=float, default=180.0)
    parser.add_argument("--ligand-qc-max-rotatable-bonds", type=int, default=15)
    parser.add_argument("--ligand-qc-max-formal-charge-abs", type=int, default=2)
    parser.add_argument("--ligand-qc-min-heavy-atom-count", type=int, default=6)
    parser.add_argument("--ligand-qc-allow-pains", action="store_true")
    parser.add_argument("--ligand-qc-allow-brenk", action="store_true")
    parser.add_argument("--ligand-qc-allow-reactive", action="store_true")
    parser.add_argument("--no-receptor-qc-gate", action="store_true")
    parser.add_argument("--receptor-qc-min-atom-count", type=int, default=100)
    parser.add_argument("--receptor-qc-min-heavy-atom-count", type=int, default=60)
    parser.add_argument("--receptor-qc-min-chain-count", type=int, default=1)
    parser.add_argument("--receptor-qc-max-coordinate-span", type=float, default=500.0)
    parser.add_argument("--dry-run", action="store_true")


def _add_dock_group(subparsers) -> None:
    dock = subparsers.add_parser("dock", help="Run docking engines")
    dock_sub = dock.add_subparsers(dest="dock_command")
    dock_run = dock_sub.add_parser("run", help="Run docking")
    _add_dock_args(dock_run)
    dock_dry = dock_sub.add_parser("dry-run", help="Generate docking commands without execution")
    _add_dock_args(dock_dry)
    dock_engine = dock_sub.add_parser("engine", help="Run one docking engine")
    dock_engine.add_argument("engine_name", choices=["gnina", "vina", "smina", "autodock4"], help="Engine name")
    _add_dock_args(dock_engine)
    dock_deploy = dock_sub.add_parser("deploy", help="Generate Slurm deployment assets")
    _add_deploy_args(dock_deploy)
    dock_sync = dock_sub.add_parser("sync", help="Sync a prepared docking project to the HPC")
    _add_sync_args(dock_sync)
    dock_submit = dock_sub.add_parser("submit", help="Submit a previously synced deployment on the HPC")
    _add_submit_args(dock_submit)


def _add_analysis_group(subparsers) -> None:
    analyze = subparsers.add_parser("analyze", help="Post-docking analysis")
    analyze_sub = analyze.add_subparsers(dest="analyze_command")
    analyze_comparative = analyze_sub.add_parser("comparative", help="Comparative multi-engine analysis")
    analyze_comparative.add_argument("--project-dir", required=True, help="Canonical project root")
    analyze_comparative.add_argument("--config-file", help="Optional post-docking analysis YAML/JSON config file")
    analyze_comparative.add_argument("-o", "--output", help="Output directory")
    analyze_comparative.add_argument("--promote-exhaustive", action="store_true", help="Emit an exhaustive rerun manifest from comparative results")
    analyze_comparative.add_argument("--rerun-engine", help="Engine to promote into exhaustive reruns")
    analyze_comparative.add_argument("--top-per-protein", type=int, default=1, help="Number of best ligands per protein to promote into the rerun manifest")
    analyze_comparative.add_argument("--winner-only", action="store_true", help="Only promote pairs where the target engine wins against the other engines")
    analyze_comparative.add_argument("--min-affinity-advantage", type=float, default=0.0, help="Require the target engine to beat the best competing engine by at least this many kcal/mol")
    analyze_comparative.add_argument("--max-rerun-pairs", type=int, default=0, help="Optional global cap on promoted exhaustive rerun pairs")
    analyze_comparative.add_argument("--pair-allowlist", help="Optional TXT/CSV file restricting promotion to explicit pair tags or receptor/site_id/ligand rows")
    analyze_comparative.add_argument(
        "--consensus-mode",
        choices=["dockbox_geometric", "weighted_hybrid", "strict_consensus", "favorite_guardrails"],
        default="dockbox_geometric",
        help="Consensus policy for comparative hit ranking and rerun promotion",
    )
    analyze_comparative.add_argument(
        "--rescoring-scope",
        choices=["top_n_per_protein", "top_n_global"],
        default="top_n_per_protein",
        help="Scope used to select rescoring candidate shortlist",
    )
    analyze_comparative.add_argument(
        "--rescoring-top-n",
        type=int,
        default=3,
        help="Top-N candidates retained by the selected rescoring scope",
    )
    analyze_comparative.add_argument(
        "--analysis-scope",
        choices=["full", "comparison_only", "rescoring_only", "top_pose_only", "qc_only", "report_only"],
        default="full",
        help="Project-centric analysis scope",
    )
    analyze_comparative.add_argument(
        "--normalization-method",
        choices=["per_engine_rank", "per_engine_minmax", "per_engine_zscore"],
        default="per_engine_rank",
        help="Normalization method used before cross-engine consensus",
    )
    analyze_comparative.add_argument("--biology-file", help="Optional biology annotation file (CSV/TSV/JSON)")
    analyze_comparative.add_argument(
        "--biology-mapping-mode",
        choices=["auto", "tag", "protein_ligand", "protein", "ligand"],
        default="auto",
        help="How to map biology rows to docking rows",
    )
    analyze_comparative.add_argument(
        "--hit-class-policy",
        choices=["target_percentile", "reference_anchor"],
        default="target_percentile",
        help="Target-aware hit classification policy",
    )
    analyze_comparative.add_argument("--hit-class-strong-percentile", type=float, default=0.10)
    analyze_comparative.add_argument("--hit-class-moderate-percentile", type=float, default=0.35)
    analyze_comparative.add_argument(
        "--top-pose-policy",
        choices=["hybrid", "best_affinity", "best_consensus"],
        default="best_affinity",
        help="Top-pose atlas selection policy",
    )
    analyze_comparative.add_argument(
        "--top-pose-aggregation",
        choices=["best_target"],
        default="best_target",
        help="Top-pose atlas global aggregation mode",
    )
    analyze_comparative.add_argument(
        "--best-pose-selection-metric",
        choices=["auto", "vina_affinity", "cnn_affinity"],
        default="auto",
        help="Metric used for best-pose complex PDB extraction",
    )
    analyze_comparative.add_argument(
        "--complex-query",
        help="Optional complex filter query, for example: protein=2FVD;ligand=Sorafenib or exclude_ligand=Acetazolamide. Use 0/skip/all to disable.",
    )
    analyze_comparative.add_argument("--prompt-protein-names", action="store_true", help="Prompt for protein display names before analysis starts")
    analyze_comparative.add_argument("--prompt-ligand-names", action="store_true", help="Prompt for ligand display names before analysis starts")
    analyze_comparative.add_argument(
        "--exclude-problematic-ligands",
        action="store_true",
        help="Exclude ligands with positive/weak affinity or missing pose support before comparative analysis",
    )
    analyze_comparative.add_argument(
        "--positive-affinity-threshold",
        type=float,
        default=0.0,
        help="Exclude ligands when best affinity is greater than this threshold (kcal/mol)",
    )
    analyze_comparative.add_argument(
        "--minimum-pose-count",
        type=int,
        default=1,
        help="Minimum required poses per ligand when problematic-ligand filtering is enabled",
    )
    analyze_comparative.add_argument(
        "--rmsd-workers",
        default="0",
        help="RMSD workers (0=auto, 1=single-thread, >1 parallel; accepts auto/single/parallel)",
    )
    analyze_favorite = analyze_sub.add_parser("favorite-engine", help="Favorite-engine continuation")
    analyze_favorite.add_argument("--project-dir", required=True, help="Canonical project root")
    analyze_favorite.add_argument("--config-file", help="Optional post-docking analysis YAML/JSON config file")
    analyze_favorite.add_argument("--favorite-engine", required=True, help="Favorite engine")
    analyze_favorite.add_argument("-o", "--output", help="Output directory")
    analyze_favorite.add_argument(
        "--complex-query",
        help="Optional complex filter query, for example: protein=2FVD;ligand=Sorafenib or exclude_ligand=Acetazolamide. Use 0/skip/all to disable.",
    )
    analyze_favorite.add_argument("--prompt-protein-names", action="store_true", help="Prompt for protein display names before analysis starts")
    analyze_favorite.add_argument("--prompt-ligand-names", action="store_true", help="Prompt for ligand display names before analysis starts")
    analyze_favorite.add_argument(
        "--exclude-problematic-ligands",
        action="store_true",
        help="Exclude ligands with positive/weak affinity or missing pose support before favorite-engine continuation",
    )
    analyze_favorite.add_argument(
        "--positive-affinity-threshold",
        type=float,
        default=0.0,
        help="Exclude ligands when best affinity is greater than this threshold (kcal/mol)",
    )
    analyze_favorite.add_argument(
        "--minimum-pose-count",
        type=int,
        default=1,
        help="Minimum required poses per ligand when problematic-ligand filtering is enabled",
    )
    analyze_favorite.add_argument(
        "--rmsd-workers",
        default="0",
        help="RMSD workers (0=auto, 1=single-thread, >1 parallel; accepts auto/single/parallel)",
    )
    analyze_favorite.add_argument(
        "--speed-profile",
        choices=["standard", "fast"],
        default="standard",
        help="Analysis speed profile (fast skips optional heavy stages where supported)",
    )
    analyze_favorite.add_argument(
        "--analysis-scope",
        choices=["full", "comparison_only", "rescoring_only", "top_pose_only", "qc_only", "report_only"],
        default="full",
        help="Project-centric analysis scope",
    )
    analyze_favorite.add_argument(
        "--normalization-method",
        choices=["per_engine_rank", "per_engine_minmax", "per_engine_zscore"],
        default="per_engine_rank",
        help="Normalization method used before cross-engine consensus",
    )
    analyze_favorite.add_argument("--biology-file", help="Optional biology annotation file (CSV/TSV/JSON)")
    analyze_favorite.add_argument(
        "--biology-mapping-mode",
        choices=["auto", "tag", "protein_ligand", "protein", "ligand"],
        default="auto",
        help="How to map biology rows to docking rows",
    )
    analyze_favorite.add_argument(
        "--hit-class-policy",
        choices=["target_percentile", "reference_anchor"],
        default="target_percentile",
        help="Target-aware hit classification policy",
    )
    analyze_favorite.add_argument("--hit-class-strong-percentile", type=float, default=0.10)
    analyze_favorite.add_argument("--hit-class-moderate-percentile", type=float, default=0.35)
    analyze_favorite.add_argument(
        "--top-pose-policy",
        choices=["hybrid", "best_affinity", "best_consensus"],
        default="best_affinity",
        help="Top-pose atlas selection policy",
    )
    analyze_favorite.add_argument(
        "--top-pose-aggregation",
        choices=["best_target"],
        default="best_target",
        help="Top-pose atlas global aggregation mode",
    )
    analyze_favorite.add_argument(
        "--best-pose-selection-metric",
        choices=["auto", "vina_affinity", "cnn_affinity"],
        default="auto",
        help="Metric used for best-pose complex PDB extraction",
    )
    _add_rmsd_args(analyze_favorite)

    stage = analyze_sub.add_parser("stage", help="Stage-specific analysis commands")
    stage_sub = stage.add_subparsers(dest="stage_command")
    for name in ["hierarchical", "polypharmacology", "rmsd", "reports", "visualizations", "structure-quality"]:
        cmd = stage_sub.add_parser(name, help=f"Run {name} only")
        _add_analysis_common_args(cmd)

    interactions = analyze_sub.add_parser("interactions", help="Interaction-only tools")
    interactions_sub = interactions.add_subparsers(dest="interactions_command")
    for name in ["pandamap", "prolif", "ligplot", "poseview"]:
        if name in {"pandamap", "prolif", "ligplot"}:
            help_text = f"Route {name} target through clean interaction pipeline"
        else:
            help_text = f"Run {name} only"
        cmd = interactions_sub.add_parser(name, help=help_text)
        _add_analysis_common_args(cmd)
    interactions_clean = interactions_sub.add_parser(
        "clean",
        help="Run clean interaction pipeline (complex-fix + PLIP + ProLIF)",
    )
    _add_clean_interaction_args(interactions_clean)

    analyze_clean = analyze_sub.add_parser(
        "clean",
        help="Shortcut for the clean interaction pipeline (same as: analyze interactions clean)",
    )
    _add_clean_interaction_args(analyze_clean)

    visuals = analyze_sub.add_parser("visuals", help="Visualization-only tools")
    visuals_sub = visuals.add_subparsers(dest="visuals_command")
    for name in ["py3dmol", "pymol"]:
        cmd = visuals_sub.add_parser(name, help=f"Run {name} only")
        _add_analysis_common_args(cmd)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Omni-DockForge unified workflow CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py workflow init --project-dir docking_project
  python main.py workflow interactive
  python main.py pdb run -p 7CMD -o results/
  python main.py prep project --prepared-proteins receptors_prep --prepared-ligands ligands_prep --excel multi_pdb_analysis.xlsx --pair-mode protein-based --output docking_project/
  python main.py dock run --project-dir docking_project --engines gnina,vina,smina,autodock4 --dry-run
  python main.py analyze comparative --project-dir docking_project
  python main.py analyze clean
  python main.py analyze stage rmsd --project-dir docking_project --engine vina
        """,
    )
    subparsers = parser.add_subparsers(dest="group")
    _add_workflow_group(subparsers)
    _add_pdb_group(subparsers)
    _add_prep_group(subparsers)
    _add_dock_group(subparsers)
    _add_analysis_group(subparsers)
    return parser
