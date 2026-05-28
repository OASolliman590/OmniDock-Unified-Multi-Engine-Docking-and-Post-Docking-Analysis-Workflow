"""
Simplified CLI for Post-Docking Analysis Pipeline.

Supports two input modes:

1. **--project-dir** (recommended for GNINA HPC workflow):
   Auto-detects directory layout (local or HPC) and resolves
   sdf/log/receptors/pairlist paths automatically.

2. **Explicit folders** (--sdf-folder, --log-folder, --receptors-folder):
   For full manual control over each input path.
"""
import argparse
import sys
from pathlib import Path
from typing import Dict, Optional


def _has_tty() -> bool:
    return sys.stdin.isatty() and sys.stdout.isatty()


def _expand_path(raw_value: str) -> Path:
    return Path(str(raw_value).strip()).expanduser()


def _detect_layout(project_dir: Path) -> Dict[str, Optional[Path]]:
    from .gnina_hpc_adapter import detect_gnina_layout

    return detect_gnina_layout(project_dir)


def _is_dockforge_manifest_project(project_dir: Optional[Path]) -> bool:
    if project_dir is None:
        return False
    try:
        from docking.project_layout import manifest_path

        return manifest_path(project_dir).exists()
    except Exception:
        return False


def _resolve_project_favorite_engine(project_dir: Path) -> Optional[str]:
    try:
        from docking.project_layout import load_manifest

        manifest = load_manifest(project_dir)
        favorite = str(manifest.get("favorite_engine", "") or "").strip().lower()
        if favorite:
            return favorite
        engines = [str(item).strip().lower() for item in (manifest.get("engines") or []) if str(item).strip()]
        if "gnina" in engines:
            return "gnina"
        return engines[0] if engines else None
    except Exception:
        return None


def _prompt_text(question: str, default: Optional[str] = None, required: bool = False) -> str:
    while True:
        suffix = f" [{default}]" if default else ""
        value = input(f"{question}{suffix}: ").strip()
        if value:
            return value
        if default is not None:
            return str(default)
        if not required:
            return ""
        print("❌ This value is required.")


def _path_matches_kind(path: Path, path_type: str) -> bool:
    if path_type == "dir":
        return path.is_dir()
    if path_type == "file":
        return path.is_file()
    return path.exists()


def _prompt_path(
    question: str,
    default: Optional[str] = None,
    path_type: str = "dir",
    must_exist: bool = True,
    required: bool = True,
) -> Optional[Path]:
    while True:
        value = _prompt_text(question, default=default, required=required)
        if not value:
            return None

        candidate = _expand_path(value)
        if not must_exist:
            return candidate
        if _path_matches_kind(candidate, path_type):
            return candidate

        expected = "directory" if path_type == "dir" else "file"
        print(f"❌ Path not found ({expected} required): {candidate}")


def _prompt_optional_path(
    question: str,
    default: Optional[Path] = None,
    path_type: str = "file",
) -> Optional[Path]:
    default_text = str(default) if default else None
    while True:
        suffix = f" [{default_text}]" if default_text else ""
        value = input(f"{question}{suffix}: ").strip()

        if value.lower() in {"skip", "none", "-"}:
            return None
        if not value:
            if default is None:
                return None
            if _path_matches_kind(default, path_type):
                return default
            print(f"❌ Default path is not valid: {default}")
            continue

        candidate = _expand_path(value)
        if _path_matches_kind(candidate, path_type):
            return candidate

        expected = "directory" if path_type == "dir" else "file"
        print(f"❌ Path not found ({expected} required): {candidate}")


def _collect_inputs_interactively(args: argparse.Namespace) -> Dict[str, Optional[Path]]:
    print("🧭 Interactive input mode")
    print("   Type paths directly. For optional values, press Enter or type 'skip'.")

    if args.project_dir:
        use_project_mode = True
    elif args.sdf_folder or args.log_folder or args.receptors_folder:
        use_project_mode = False
    else:
        mode = _prompt_text(
            "Use GNINA project-root auto detection? (y/n)",
            default="y",
            required=True,
        ).lower()
        use_project_mode = mode in {"y", "yes"}

    project_dir = None
    sdf_folder = None
    log_folder = None
    receptors_folder = None
    pairlist_file = None

    if use_project_mode:
        project_default = str(args.project_dir).strip() if args.project_dir else None
        project_dir = _prompt_path(
            "GNINA project root directory",
            default=project_default,
            path_type="dir",
            must_exist=True,
            required=True,
        )
        layout = _detect_layout(project_dir)
        print(f"📂 Auto-detected GNINA layout: {layout['layout'].upper()}")

        sdf_folder = layout.get("sdf_folder")
        log_folder = layout.get("log_folder")
        receptors_folder = layout.get("receptors_folder")
        pairlist_file = layout.get("pairlist_file")

        if not sdf_folder:
            sdf_folder = _prompt_path("SDF folder", path_type="dir", must_exist=True, required=True)
        if not log_folder:
            log_folder = _prompt_path("Log folder", path_type="dir", must_exist=True, required=True)
        if not receptors_folder:
            receptors_folder = _prompt_path("Receptors folder", path_type="dir", must_exist=True, required=True)

        if args.pairlist:
            pairlist_file = _expand_path(args.pairlist)
        pairlist_file = _prompt_optional_path(
            "pairlist.csv path (optional, Enter keeps default, 'skip' ignores)",
            default=pairlist_file,
            path_type="file",
        )
    else:
        sdf_folder = _prompt_path(
            "SDF folder",
            default=args.sdf_folder,
            path_type="dir",
            must_exist=True,
            required=True,
        )
        log_folder = _prompt_path(
            "Log folder",
            default=args.log_folder,
            path_type="dir",
            must_exist=True,
            required=True,
        )
        receptors_folder = _prompt_path(
            "Receptors folder",
            default=args.receptors_folder,
            path_type="dir",
            must_exist=True,
            required=True,
        )

        pairlist_default = _expand_path(args.pairlist) if args.pairlist else None
        pairlist_file = _prompt_optional_path(
            "pairlist.csv path (optional)",
            default=pairlist_default,
            path_type="file",
        )

    output_dir = _prompt_path(
        "Output directory",
        default=args.output,
        path_type="dir",
        must_exist=False,
        required=True,
    )

    if args.prompt_protein_names:
        prompt_protein_names = True
    else:
        prompt_answer = _prompt_text(
            "Prompt protein names for each detected PDB/receptor? (y/n)",
            default="y",
            required=True,
        ).lower()
        prompt_protein_names = prompt_answer in {"y", "yes"}

    if args.prompt_ligand_names:
        prompt_ligand_names = True
    else:
        ligand_prompt_answer = _prompt_text(
            "Prompt ligand names for detected ligands? (y/n)",
            default="n",
            required=True,
        ).lower()
        prompt_ligand_names = ligand_prompt_answer in {"y", "yes"}

    complex_query = _prompt_text(
        "Optional complex filter query (press Enter or type 0/skip to include all)",
        default=args.complex_query or "",
        required=False,
    )

    return {
        "project_dir": project_dir,
        "sdf_folder": sdf_folder,
        "log_folder": log_folder,
        "receptors_folder": receptors_folder,
        "pairlist_file": pairlist_file,
        "output_dir": output_dir,
        "prompt_protein_names": prompt_protein_names,
        "prompt_ligand_names": prompt_ligand_names,
        "complex_query": complex_query,
    }


def main():
    """Main CLI function."""
    parser = argparse.ArgumentParser(
        description="Simplified Post-Docking Analysis Pipeline (GNINA Focus)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect layout from GNINA HPC project root
  python -m post_docking_analysis.simplified_cli \\
    --project-dir /path/to/Sertaline_docking \\
    --output /path/to/output

  # Basic usage with 3 folders (explicit)
  python -m post_docking_analysis.simplified_cli \\
    --sdf-folder /path/to/sdf \\
    --log-folder /path/to/logs \\
    --receptors-folder /path/to/receptors \\
    --output /path/to/output

  # Interactive prompting for missing fields
  python -m post_docking_analysis.simplified_cli --interactive
        """
    )

    # Auto-detect mode
    parser.add_argument(
        "--project-dir",
        help="GNINA project root directory (auto-detects HPC or local layout). "
             "Overrides --sdf-folder, --log-folder, --receptors-folder, and --pairlist."
    )

    # Explicit folder mode
    parser.add_argument("--sdf-folder", help="Folder containing SDF pose files")
    parser.add_argument("--log-folder", help="Folder containing log files (can be same as sdf-folder)")
    parser.add_argument("--receptors-folder", help="Folder containing receptor PDBQT files")
    parser.add_argument("--output", help="Output directory for results")

    # Optional arguments
    parser.add_argument("--pairlist", help="Path to pairlist.csv for receptor-ligand mapping")
    parser.add_argument("--config-file", help="Optional post-docking analysis YAML/JSON config file")
    parser.add_argument(
        "--complex-query",
        help="Optional complex filter query, for example: protein=2FVD;ligand=Sorafenib or exclude_ligand=Acetazolamide. Use 0/skip/all to disable.",
    )
    parser.add_argument(
        "--ligplus-root",
        help="Optional LigPlus root path for auto LigPlot generation. "
             "If omitted, pipeline checks LIGPLUS_ROOT / LIGPLUS_HOME."
    )
    parser.add_argument(
        "--prompt-protein-names",
        action="store_true",
        help="Prompt for protein display names before analysis starts"
    )
    parser.add_argument(
        "--prompt-ligand-names",
        action="store_true",
        help="Prompt for ligand display names before analysis starts"
    )
    parser.add_argument(
        "--enable-poseview",
        action="store_true",
        help="Enable PoseView REST API interaction diagrams (requires internet access)"
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Force interactive prompts for input paths and naming settings"
    )

    # Analysis options
    parser.add_argument(
        "--no-rmsd",
        action="store_true",
        help="Deprecated compatibility flag. RMSD is always required in the clean interaction pipeline.",
    )
    parser.add_argument(
        "--no-visualizations",
        action="store_true",
        help="Deprecated compatibility flag. Core visualization/interaction stages are always required.",
    )
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

    args = parser.parse_args()

    if args.no_rmsd or args.no_visualizations:
        print(
            "❌ --no-rmsd/--no-visualizations are no longer supported.\n"
            "   Clean post-docking interaction runs require RMSD and core interaction visual stages."
        )
        return 1

    interactive_tty = _has_tty()
    has_explicit_folders = bool(args.sdf_folder and args.log_folder and args.receptors_folder)
    missing_required = (not args.output) or (not args.project_dir and not has_explicit_folders)
    should_prompt = bool(args.interactive or (interactive_tty and missing_required))

    if args.interactive and not interactive_tty:
        print("❌ --interactive requires a terminal TTY session.")
        return 1

    prompted_session = False
    prompt_protein_names = bool(args.prompt_protein_names)
    prompt_ligand_names = bool(args.prompt_ligand_names)
    resolved_project_dir: Optional[Path] = None
    sdf_folder: Optional[Path] = None
    log_folder: Optional[Path] = None
    receptors_folder: Optional[Path] = None
    output_dir: Optional[Path] = None
    pairlist_file: Optional[Path] = None
    complex_query = str(args.complex_query or "").strip()

    if should_prompt:
        if not interactive_tty:
            print(
                "❌ Missing required arguments and no interactive terminal detected.\n"
                "   Provide required flags explicitly (--output and input paths), or run in a TTY."
            )
            return 1
        try:
            resolved = _collect_inputs_interactively(args)
        except KeyboardInterrupt:
            print("\n⏹️ Interactive input cancelled.")
            return 1

        prompted_session = True
        sdf_folder = resolved["sdf_folder"]
        log_folder = resolved["log_folder"]
        receptors_folder = resolved["receptors_folder"]
        output_dir = resolved["output_dir"]
        resolved_project_dir = resolved["project_dir"]
        pairlist_file = resolved["pairlist_file"]
        prompt_protein_names = bool(resolved["prompt_protein_names"])
        prompt_ligand_names = bool(resolved["prompt_ligand_names"])
        complex_query = str(resolved.get("complex_query") or "").strip()
    else:
        if not args.output:
            print("❌ --output is required in non-interactive mode")
            return 1

        output_dir = _expand_path(args.output)
        pairlist_file = _expand_path(args.pairlist) if args.pairlist else None

        if args.project_dir:
            project_dir = _expand_path(args.project_dir)
            if not project_dir.exists():
                print(f"❌ Project directory does not exist: {project_dir}")
                return 1
            resolved_project_dir = project_dir

            layout = _detect_layout(project_dir)
            print(f"📂 Auto-detected GNINA layout: {layout['layout'].upper()}")

            sdf_folder = layout["sdf_folder"]
            log_folder = layout["log_folder"]
            receptors_folder = layout["receptors_folder"]
            if not pairlist_file:
                pairlist_file = layout["pairlist_file"]

            if not sdf_folder:
                print("❌ Could not find SDF/gnina_out folder in project directory")
                return 1
            if not log_folder:
                print("❌ Could not find log files in project directory")
                return 1
            if not receptors_folder:
                print("❌ Could not find receptors folder in project directory")
                return 1
        else:
            if not has_explicit_folders:
                print(
                    "❌ Must provide either --project-dir OR all of "
                    "--sdf-folder, --log-folder, --receptors-folder"
                )
                return 1

            sdf_folder = _expand_path(args.sdf_folder)
            log_folder = _expand_path(args.log_folder)
            receptors_folder = _expand_path(args.receptors_folder)

            if not sdf_folder.exists():
                print(f"❌ SDF folder does not exist: {sdf_folder}")
                return 1
            if not log_folder.exists():
                print(f"❌ Log folder does not exist: {log_folder}")
                return 1
            if not receptors_folder.exists():
                print(f"❌ Receptors folder does not exist: {receptors_folder}")
                return 1

    if prompted_session and not args.prompt_protein_names:
        # Interactive sessions default to explicit per-target naming.
        prompt_protein_names = True

    should_delegate_unified = _is_dockforge_manifest_project(resolved_project_dir)
    config_overrides = {}
    resolved_config_file = str(args.config_file or "").strip()
    if resolved_config_file:
        try:
            from .config_manager import load_config_overrides

            config_overrides = load_config_overrides(resolved_config_file)
            print(f"🧩 Loaded analysis config overrides: {resolved_config_file}")
        except Exception as exc:
            print(f"❌ Could not load analysis config file {resolved_config_file}: {exc}")
            return 1
    if should_delegate_unified:
        favorite_engine = _resolve_project_favorite_engine(resolved_project_dir) or "gnina"
        print(
            "ℹ️ Legacy simplified CLI wrapper: delegating execution to unified post-docking pipeline."
        )
        from workflow.execution import run_analysis_favorite

        result = run_analysis_favorite(
            project_dir=str(resolved_project_dir),
            favorite_engine=favorite_engine,
            output_dir=str(output_dir),
            enable_poseview=args.enable_poseview,
            prompt_protein_names=prompt_protein_names,
            prompt_ligand_names=prompt_ligand_names,
            complex_query=complex_query or None,
            rmsd_scopes=args.rmsd_scopes,
            resume_rmsd=not args.no_resume_rmsd,
            speed_profile="standard",
            analysis_scope="full",
            config_file=resolved_config_file or None,
        )
        if str(result.status).lower() == "completed":
            print("\n✅ Post-docking analysis completed successfully!")
            print(f"📂 Results saved to: {output_dir}")
            return 0
        print("\n❌ Post-docking analysis failed!")
        return 1

    from .simplified_pipeline import SimplifiedPostDockingPipeline

    pipeline = SimplifiedPostDockingPipeline(
        sdf_folder=str(sdf_folder),
        log_folder=str(log_folder),
        receptors_folder=str(receptors_folder),
        output_dir=str(output_dir),
        pairlist_file=str(pairlist_file) if pairlist_file else None,
        run_rmsd=True,
        run_visualizations=True,
        ligplus_root=args.ligplus_root,
        prompt_protein_names=prompt_protein_names,
        prompt_ligand_names=prompt_ligand_names,
        enable_poseview=args.enable_poseview,
        complex_query=complex_query or None,
        rmsd_scopes=args.rmsd_scopes,
        resume_rmsd=not args.no_resume_rmsd,
        analysis_config=config_overrides,
    )

    success = pipeline.run()

    if success:
        print("\n✅ Post-docking analysis completed successfully!")
        print(f"📂 Results saved to: {output_dir}")
        return 0

    print("\n❌ Post-docking analysis failed!")
    return 1


if __name__ == "__main__":
    sys.exit(main())
