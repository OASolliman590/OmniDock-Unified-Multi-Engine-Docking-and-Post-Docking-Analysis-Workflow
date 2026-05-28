"""
CLI for batch LigPlot+ generation using LigPlus executables.
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import List, Optional

from .ligplot_integration import (
    LigPlotRunner,
    export_ligplot_outputs,
    iter_pdb_files,
    list_ligands_from_pairlist,
)


def _parse_filters(choice: str, ligands: List[str]) -> List[str]:
    if not choice:
        return []
    filters: List[str] = []
    parts = [p.strip() for p in choice.split(",") if p.strip()]
    for part in parts:
        if part.isdigit() and ligands:
            idx = int(part) - 1
            if 0 <= idx < len(ligands):
                filters.append(ligands[idx])
        else:
            filters.append(part)
    return filters


def _prompt_ligand_filters(ligands: List[str]) -> List[str]:
    if not ligands:
        choice = input("Ligand filter (substring, empty for all): ").strip()
        return _parse_filters(choice, [])

    preview = ligands[:30]
    print("Available ligands (showing first 30):")
    for i, ligand in enumerate(preview, 1):
        print(f"  {i}. {ligand}")
    if len(ligands) > len(preview):
        print(f"  ... and {len(ligands) - len(preview)} more")

    choice = input(
        "Select ligand by number or substring (comma-separated, empty for all): "
    ).strip()
    return _parse_filters(choice, ligands)


def _filter_files(pdb_files: List[Path], filters: List[str]) -> List[Path]:
    if not filters:
        return pdb_files
    lowered = [f.lower() for f in filters]
    selected = [
        pdb_file
        for pdb_file in pdb_files
        if any(f in pdb_file.name.lower() for f in lowered)
    ]
    return selected


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Batch LigPlot+ generation for PDB complexes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-dir", required=True, help="Root folder with PDB files")
    parser.add_argument(
        "--output-dir",
        help="Output folder for LigPlot results (default: sibling ligplot_output)",
    )
    parser.add_argument(
        "--ligplus-root",
        help="LigPlus root folder (contains LigPlus.jar and lib/)",
    )
    parser.add_argument("--exe-dir", help="Override LigPlus executables folder")
    parser.add_argument("--param-dir", help="Override LigPlus params folder")
    parser.add_argument("--pairlist", help="Optional pairlist.csv for ligand list")
    parser.add_argument("--ligand-filter", help="Substring filter for ligand selection")
    parser.add_argument("--no-interactive", action="store_true", help="Skip prompts")
    parser.add_argument("--no-recursive", action="store_true", help="Do not recurse")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite outputs")
    parser.add_argument("--dry-run", action="store_true", help="List files only")
    parser.add_argument("--limit", type=int, default=0, help="Process only N files")
    parser.add_argument("--skip-hbplus", action="store_true", help="Skip HBADD/HBPLUS runs")
    parser.add_argument("--skip-hbadd", action="store_true", help="Skip HBADD run")
    parser.add_argument(
        "--hydrogenate",
        choices=["none", "obabel", "auto"],
        default="auto",
        help="Add hydrogens via OpenBabel (auto retries if initial plot fails)",
    )
    parser.add_argument(
        "--strip-metals",
        choices=["none", "auto", "always"],
        default="auto",
        help="Remove metal ATOM records if LigPlot crashes (auto aliases metals first, then strips on failure)",
    )
    parser.add_argument(
        "--hbplus-hparam",
        type=float,
        default=2.70,
        help="HBPLUS hydrogen distance cutoff (H-A)",
    )
    parser.add_argument(
        "--hbplus-dparam",
        type=float,
        default=3.35,
        help="HBPLUS donor-acceptor distance cutoff (D-A)",
    )
    parser.add_argument(
        "--hbplus-nb-hparam",
        type=float,
        default=2.90,
        help="HBPLUS nonbond hydrogen cutoff (H-A)",
    )
    parser.add_argument(
        "--hbplus-nb-dparam",
        type=float,
        default=3.90,
        help="HBPLUS nonbond distance cutoff (D-A)",
    )
    parser.add_argument(
        "--contact-type",
        choices=["0", "1", "2"],
        default="1",
        help="Non-bonded contact mode (0=hydrophobic only, 1=carbon/sulfur vs any, 2=all nonbonded)",
    )
    parser.add_argument(
        "--export",
        choices=["none", "png", "pdf", "both"],
        default="none",
        help="Export ligplot.ps to png/pdf each run",
    )
    parser.add_argument(
        "--export-dpi",
        type=int,
        default=300,
        help="DPI for PNG export (ghostscript/magick)",
    )

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    input_dir = Path(args.input_dir).expanduser().resolve()
    if not input_dir.exists():
        print(f"Input directory not found: {input_dir}")
        return 1

    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else input_dir.parent / "ligplot_output"
    )

    default_root = os.environ.get("LIGPLUS_ROOT") or os.environ.get("LIGPLUS_HOME")
    if not default_root:
        repo_root = Path(__file__).resolve().parents[1]
        for candidate in (
            repo_root / ".workflow" / "tools" / "LigPlus",
            repo_root / "tools" / "LigPlus",
            repo_root / "ligplus",
        ):
            if candidate.exists():
                default_root = str(candidate.resolve())
                break

    ligplus_root = (
        Path(args.ligplus_root or default_root).expanduser().resolve()
        if (args.ligplus_root or default_root)
        else None
    )
    if ligplus_root is None:
        print("LigPlus root is required. Use --ligplus-root.")
        return 1

    pairlist_file: Optional[Path] = None
    if args.pairlist:
        pairlist_file = Path(args.pairlist).expanduser().resolve()
    else:
        candidate = input_dir.parent / "pairlist.csv"
        if candidate.exists():
            pairlist_file = candidate

    pdb_files = iter_pdb_files(input_dir, recursive=not args.no_recursive)
    if not pdb_files:
        print(f"No PDB files found in {input_dir}")
        return 1

    ligands: List[str] = []
    if pairlist_file:
        ligands = list_ligands_from_pairlist(pairlist_file)

    filters: List[str] = []
    if args.ligand_filter:
        filters = _parse_filters(args.ligand_filter, ligands)
    elif not args.no_interactive:
        filters = _prompt_ligand_filters(ligands)

    pdb_files = _filter_files(pdb_files, filters)
    if args.limit > 0:
        pdb_files = pdb_files[: args.limit]

    if not pdb_files:
        print("No files matched the ligand filter.")
        return 1

    if args.dry_run:
        print(f"Would process {len(pdb_files)} PDB files")
        for pdb_file in pdb_files:
            print(pdb_file)
        return 0

    runner = LigPlotRunner(
        ligplus_root=ligplus_root,
        exe_dir=Path(args.exe_dir).expanduser().resolve() if args.exe_dir else None,
        param_dir=Path(args.param_dir).expanduser().resolve() if args.param_dir else None,
        contact_type=args.contact_type,
        hbplus_hparam=f"{args.hbplus_hparam:.2f}",
        hbplus_dparam=f"{args.hbplus_dparam:.2f}",
        hbplus_nb_hparam=f"{args.hbplus_nb_hparam:.2f}",
        hbplus_nb_dparam=f"{args.hbplus_nb_dparam:.2f}",
        use_hbplus=not args.skip_hbplus,
        use_hbadd=not args.skip_hbadd,
        hydrogenate=args.hydrogenate,
        strip_metals=args.strip_metals,
    )

    try:
        runner.validate()
    except Exception as exc:
        print(f"LigPlot validation failed: {exc}")
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)
    total = len(pdb_files)
    success = 0

    for idx, pdb_file in enumerate(pdb_files, 1):
        rel = pdb_file.relative_to(input_dir)
        target_dir = output_dir / rel.parent / pdb_file.stem
        target_dir.mkdir(parents=True, exist_ok=True)

        logging.info("[%d/%d] %s", idx, total, pdb_file.name)
        ok = runner.run_single(
            pdb_file,
            target_dir,
            overwrite=args.overwrite,
        )
        if ok:
            if args.export != "none":
                export_ok = export_ligplot_outputs(
                    target_dir,
                    [args.export],
                    overwrite=args.overwrite,
                    dpi=args.export_dpi,
                )
                if not export_ok:
                    logging.warning("Export failed for %s", pdb_file.name)
            success += 1

    logging.info("Completed: %d/%d succeeded", success, total)
    logging.info("Output folder: %s", output_dir)
    return 0 if success == total else 1


if __name__ == "__main__":
    sys.exit(main())
