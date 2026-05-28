"""
LigPlot+ batch integration using LigPlus executables.

This module runs HBADD + HBPLUS + LIGPLOT in a temporary working directory
and copies LigPlot outputs to a mirrored output structure.
"""
from __future__ import annotations

import csv
import logging
import os
import platform
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

logger = logging.getLogger(__name__)

WATER_RESNAMES = {"HOH", "WAT", "SOL", "H2O", "DOD"}
COMMON_NONLIGAND_RESNAMES = {
    "SO4",
    "PO4",
    "DMS",
    "EDO",
    "GOL",
    "ACT",
}
METAL_ELEMENTS = {
    "ZN",
    "FE",
    "MG",
    "MN",
    "CU",
    "CO",
    "NI",
    "CA",
    "CD",
    "HG",
    "NA",
    "K",
}
METAL_ALIAS_ELEMENT = "S"


@dataclass
class LigPlotOptions:
    # "2" includes all non-bonded contacts, giving richer interaction diagrams.
    contact_type: str = "2"
    plot_metal: bool = False
    include_waters: bool = False
    filter_waters: bool = False
    # Keep batch runs resilient to per-file edge cases.
    no_abort: bool = True
    use_hbplus: bool = True
    use_hbadd: bool = True
    # Auto hydrogenation tends to improve interaction completeness.
    hydrogenate: str = "auto"
    strip_metals: str = "auto"
    hbplus_hparam: str = "2.70"
    hbplus_dparam: str = "3.35"
    hbplus_nb_hparam: str = "2.90"
    hbplus_nb_dparam: str = "3.90"


@dataclass
class LigPlotPaths:
    ligplus_root: Path
    exe_dir: Path
    param_dir: Path
    ligplot_bin: Path
    hbplus_bin: Path
    hbadd_bin: Path
    het_dictionary: Path


def _strip_quotes(value: str) -> str:
    value = value.strip()
    if (value.startswith("'") and value.endswith("'")) or (
        value.startswith('"') and value.endswith('"')
    ):
        return value[1:-1]
    return value


def read_ligplus_params(ligplus_root: Path) -> Dict[str, str]:
    params_path = ligplus_root / "lib" / "params" / "ligplus.par"
    params: Dict[str, str] = {}
    if not params_path.exists():
        return params
    with params_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            params[key.strip()] = value.strip()
    return params


def resolve_paths(
    ligplus_root: Path,
    exe_dir: Optional[Path] = None,
    param_dir: Optional[Path] = None,
) -> LigPlotPaths:
    ligplus_root = ligplus_root.resolve()
    params = read_ligplus_params(ligplus_root)

    if exe_dir is None:
        system = platform.system().lower()
        machine = platform.machine().lower()
        candidates: List[Path] = []
        if system == "darwin":
            if machine == "arm64":
                candidates = [ligplus_root / "lib" / "exe_mac64"]
            else:
                candidates = [ligplus_root / "lib" / "exe_mac"]
            candidates.append(ligplus_root / "lib" / "exe_mac64")
        elif system == "linux":
            candidates = [
                ligplus_root / "lib" / "exe_linux64",
                ligplus_root / "lib" / "exe_linux",
            ]
        elif system.startswith("win"):
            candidates = [
                ligplus_root / "lib" / "exe_win",
                ligplus_root / "lib" / "exe_win32",
            ]
        for candidate in candidates:
            if candidate.exists():
                exe_dir = candidate
                break
        if exe_dir is None:
            exe_dir = ligplus_root / "lib" / "exe"

    if param_dir is None:
        param_dir = ligplus_root / "lib" / "params"

    ligplot_bin = exe_dir / "ligplot"
    hbplus_bin = exe_dir / "hbplus"
    hbadd_bin = exe_dir / "hbadd"

    het_dictionary = params.get("HET_GROUP_DICTIONARY", "")
    if not het_dictionary:
        het_dictionary = str(ligplus_root / "components.cif")
    het_dictionary = Path(_strip_quotes(het_dictionary)).expanduser()

    return LigPlotPaths(
        ligplus_root=ligplus_root,
        exe_dir=exe_dir,
        param_dir=param_dir,
        ligplot_bin=ligplot_bin,
        hbplus_bin=hbplus_bin,
        hbadd_bin=hbadd_bin,
        het_dictionary=het_dictionary,
    )


def build_options(
    params: Dict[str, str],
    contact_type: Optional[str] = None,
    hbplus_hparam: Optional[str] = None,
    hbplus_dparam: Optional[str] = None,
    hbplus_nb_hparam: Optional[str] = None,
    hbplus_nb_dparam: Optional[str] = None,
    plot_metal: bool = False,
    include_waters: bool = False,
    filter_waters: bool = False,
    no_abort: Optional[bool] = None,
    use_hbplus: bool = True,
    use_hbadd: bool = True,
    hydrogenate: Optional[str] = None,
    strip_metals: Optional[str] = None,
) -> LigPlotOptions:
    options = LigPlotOptions()
    if contact_type:
        options.contact_type = contact_type
    else:
        options.contact_type = params.get("CONTACT_TYPE", options.contact_type)
    options.plot_metal = plot_metal
    options.include_waters = include_waters
    options.filter_waters = filter_waters
    if no_abort is None:
        options.no_abort = params.get("NO_ABORT", "FALSE").upper() == "TRUE"
    else:
        options.no_abort = no_abort
    options.use_hbplus = use_hbplus
    options.use_hbadd = use_hbadd
    options.hydrogenate = (hydrogenate or "auto").lower()
    options.strip_metals = (strip_metals or "auto").lower()
    options.hbplus_hparam = hbplus_hparam or params.get(
        "HBPLUS_HPARAM", options.hbplus_hparam
    )
    options.hbplus_dparam = hbplus_dparam or params.get(
        "HBPLUS_DPARAM", options.hbplus_dparam
    )
    options.hbplus_nb_hparam = hbplus_nb_hparam or params.get(
        "HBPLUS_NB_HPARAM", options.hbplus_nb_hparam
    )
    options.hbplus_nb_dparam = hbplus_nb_dparam or params.get(
        "HBPLUS_NB_DPARAM", options.hbplus_nb_dparam
    )
    return options


def detect_ligand_range(
    pdb_file: Path, ligand_resname_hint: Optional[str] = None
) -> Optional[Tuple[str, str, str, str, str]]:
    residues: List[Tuple[str, str, str, str]] = []
    atom_counts: Dict[Tuple[str, str, str, str], int] = {}
    hinted_counts: Dict[Tuple[str, str, str, str], int] = {}
    seen = set()
    hint = (ligand_resname_hint or "").strip().upper()

    with pdb_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip()
            if not resname or resname in WATER_RESNAMES or resname in COMMON_NONLIGAND_RESNAMES:
                continue
            chain = line[21].strip() or " "
            resnum = line[22:26].strip()
            icode = line[26].strip()
            if not resnum:
                continue
            key = (resname, chain, resnum, icode)
            atom_counts[key] = atom_counts.get(key, 0) + 1
            if hint and resname.upper() == hint:
                hinted_counts[key] = hinted_counts.get(key, 0) + 1
            if key not in seen:
                residues.append(key)
                seen.add(key)

    if not atom_counts:
        return None

    selected_pool = hinted_counts if hinted_counts else atom_counts
    selected = max(selected_pool, key=selected_pool.get)
    resname, chain, _resnum, _icode = selected

    same_res = [key for key in residues if key[0] == resname and key[1] == chain]
    first_resnum = same_res[0][2]
    last_resnum = same_res[-1][2]
    return resname, first_resnum, resname, last_resnum, chain


def _run_command(cmd: List[str], cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        check=False,
    )


def _run_tool(cmd: List[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )


def _hydrogenate_with_obabel(input_pdb: Path, output_pdb: Path) -> bool:
    cmd = ["obabel", "-ipdb", str(input_pdb), "-opdb", "-O", str(output_pdb), "-h"]
    result = _run_tool(cmd)
    if result.returncode != 0 or not output_pdb.exists():
        logger.warning(
            "OpenBabel hydrogenation failed: %s%s",
            result.stdout.strip(),
            f" | {result.stderr.strip()}" if result.stderr.strip() else "",
        )
        return False
    return True


def _is_metal_atom(line: str) -> bool:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return False
    element = line[76:78].strip().upper()
    if not element:
        element = line[12:16].strip().upper()
    return element in METAL_ELEMENTS


def _has_metal_atoms(pdb_file: Path) -> bool:
    with pdb_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("ATOM") and _is_metal_atom(line):
                return True
    return False


def _strip_metals_from_pdb(input_pdb: Path, output_pdb: Path) -> int:
    removed = 0
    with input_pdb.open("r", encoding="utf-8") as inp, output_pdb.open(
        "w", encoding="utf-8"
    ) as out:
        for line in inp:
            if line.startswith("ATOM") and _is_metal_atom(line):
                removed += 1
                continue
            out.write(line)
    return removed


def _format_atom_name(element: str) -> str:
    element = element.strip().upper()
    if not element:
        return "    "
    if len(element) == 1:
        return f" {element}  "
    return f"{element[:2]:<4}"


def _alias_metals_for_ligplot(
    input_pdb: Path, output_pdb: Path, alias_element: str = METAL_ALIAS_ELEMENT
) -> int:
    replaced = 0
    atom_field = _format_atom_name(alias_element)
    element_field = alias_element.strip().upper().rjust(2)
    with input_pdb.open("r", encoding="utf-8") as inp, output_pdb.open(
        "w", encoding="utf-8"
    ) as out:
        for line in inp:
            if line.startswith("ATOM") and _is_metal_atom(line):
                raw = line.rstrip("\n")
                newline = "\n" if line.endswith("\n") else ""
                if len(raw) < 80:
                    raw = raw.ljust(80)
                raw = raw[:12] + atom_field + raw[16:]
                raw = raw[:76] + element_field + raw[78:]
                line = raw + newline
                replaced += 1
            out.write(line)
    return replaced


def _clear_generated_files(tmp_dir: Path) -> None:
    for path in tmp_dir.glob("ligplot.*"):
        path.unlink(missing_ok=True)
    for name in ("ligplus.hhb", "ligplus.nnb", "hbdebug.dat", "hbadd.bonds"):
        (tmp_dir / name).unlink(missing_ok=True)


def _ligplot_outputs_ok(tmp_dir: Path) -> bool:
    drw_file = tmp_dir / "ligplot.drw"
    ps_file = tmp_dir / "ligplot.ps"
    sum_file = tmp_dir / "ligplot.sum"
    if not drw_file.exists() or drw_file.stat().st_size == 0:
        return False
    if not ps_file.exists() or ps_file.stat().st_size == 0:
        return False
    if sum_file.exists() and sum_file.stat().st_size == 0:
        return False
    return True


def _describe_binary(path: Path) -> str:
    try:
        result = subprocess.run(
            ["file", str(path)], capture_output=True, text=True, check=False
        )
        return result.stdout.strip()
    except FileNotFoundError:
        return ""


class LigPlotRunner:
    def __init__(
        self,
        ligplus_root: Path,
        exe_dir: Optional[Path] = None,
        param_dir: Optional[Path] = None,
        contact_type: Optional[str] = None,
        hbplus_hparam: Optional[str] = None,
        hbplus_dparam: Optional[str] = None,
        hbplus_nb_hparam: Optional[str] = None,
        hbplus_nb_dparam: Optional[str] = None,
        plot_metal: bool = False,
        include_waters: bool = False,
        filter_waters: bool = False,
        no_abort: Optional[bool] = None,
        use_hbplus: bool = True,
        use_hbadd: bool = True,
        hydrogenate: str = "none",
        strip_metals: str = "auto",
    ) -> None:
        self.params = read_ligplus_params(ligplus_root)
        self.paths = resolve_paths(ligplus_root, exe_dir=exe_dir, param_dir=param_dir)
        self.options = build_options(
            self.params,
            contact_type=contact_type,
            hbplus_hparam=hbplus_hparam,
            hbplus_dparam=hbplus_dparam,
            hbplus_nb_hparam=hbplus_nb_hparam,
            hbplus_nb_dparam=hbplus_nb_dparam,
            plot_metal=plot_metal,
            include_waters=include_waters,
            filter_waters=filter_waters,
            no_abort=no_abort,
            use_hbplus=use_hbplus,
            use_hbadd=use_hbadd,
            hydrogenate=hydrogenate,
            strip_metals=strip_metals,
        )

        self.last_run_report: Dict[str, object] = {}

    def validate(self) -> None:
        missing = []
        for path in [self.paths.ligplot_bin, self.paths.hbplus_bin, self.paths.hbadd_bin]:
            if not path.exists():
                missing.append(str(path))
        if missing:
            raise FileNotFoundError(f"Missing LigPlus executables: {', '.join(missing)}")

        if platform.system().lower() == "darwin":
            desc = _describe_binary(self.paths.ligplot_bin)
            machine = platform.machine().lower()
            if machine == "x86_64" and "arm64" in desc:
                raise RuntimeError(
                    "LigPlot binary is arm64 but this Python is x86_64. "
                    "Install the x86_64 LigPlus build or run in an arm64 shell."
                )
        if self.options.hydrogenate == "obabel" and not shutil.which("obabel"):
            raise RuntimeError("Open Babel (obabel) not found in PATH; cannot hydrogenate.")

    def run_single(
        self,
        pdb_file: Path,
        output_dir: Path,
        ligand_resname_hint: Optional[str] = None,
        overwrite: bool = False,
    ) -> bool:
        self.last_run_report = {}
        output_dir = output_dir.resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        output_ps = output_dir / "ligplot.ps"
        if output_ps.exists() and not overwrite:
            logger.info("Skipping %s (output exists)", pdb_file.name)
            self.last_run_report = {
                "pdb_file": str(pdb_file),
                "success": True,
                "skipped_existing": True,
                "failure_reason": "",
                "attempts": [],
            }
            return True

        ligand_range = detect_ligand_range(pdb_file, ligand_resname_hint)
        if not ligand_range:
            reason = "No ligand residues detected"
            logger.warning("%s in %s", reason, pdb_file.name)
            self.last_run_report = {
                "pdb_file": str(pdb_file),
                "success": False,
                "skipped_existing": False,
                "failure_reason": reason,
                "attempts": [],
            }
            return False

        resname_start, resnum_start, resname_end, resnum_end, chain = ligand_range
        resname_arg_start = ""
        resname_arg_end = ""
        if chain == " ":
            resname_arg_start = f"-n{resname_start}"
            resname_arg_end = f"-n{resname_end}"

        ligplot_prm = self.paths.param_dir / "ligplot.prm"
        if not ligplot_prm.exists():
            raise FileNotFoundError(f"ligplot.prm not found: {ligplot_prm}")

        with tempfile.TemporaryDirectory(prefix="ligplot_") as tmpdir:
            tmp_dir = Path(tmpdir)
            orig_pdb = tmp_dir / "orig.pdb"
            work_pdb = tmp_dir / "ligplus.pdb"
            shutil.copy2(pdb_file, orig_pdb)

            attempt_reports: List[Dict[str, object]] = []

            def _first_line(text_value: str) -> str:
                compact = str(text_value or "").strip()
                if not compact:
                    return ""
                return compact.splitlines()[0][:220]

            def attempt_run(hydrogenate: bool, metal_mode: str) -> Dict[str, object]:
                _clear_generated_files(tmp_dir)
                shutil.copy2(orig_pdb, work_pdb)

                report: Dict[str, object] = {
                    "hydrogenate": bool(hydrogenate),
                    "metal_mode": str(metal_mode),
                    "hbadd_rc": None,
                    "hbplus_nb_rc": None,
                    "hbplus_rc": None,
                    "ligplot_rc": None,
                    "output_ok": False,
                    "errors": [],
                    "generated_outputs": [],
                }

                if hydrogenate:
                    hyd_pdb = tmp_dir / "ligplus_h.pdb"
                    if _hydrogenate_with_obabel(work_pdb, hyd_pdb):
                        hyd_pdb.replace(work_pdb)
                    else:
                        report["errors"].append("OpenBabel hydrogenation failed")

                if metal_mode == "alias":
                    alias_pdb = tmp_dir / "ligplus_alias.pdb"
                    replaced = _alias_metals_for_ligplot(work_pdb, alias_pdb)
                    report["aliased_metals"] = int(replaced)
                    alias_pdb.replace(work_pdb)
                elif metal_mode == "strip":
                    stripped_pdb = tmp_dir / "ligplus_strip.pdb"
                    removed = _strip_metals_from_pdb(work_pdb, stripped_pdb)
                    report["stripped_metals"] = int(removed)
                    stripped_pdb.replace(work_pdb)

                if self.options.use_hbplus:
                    hbadd_cmd = [
                        str(self.paths.hbadd_bin),
                        str(work_pdb),
                        str(self.paths.het_dictionary),
                        "-wkdir",
                        str(tmp_dir) + os.sep,
                    ]
                    if self.options.use_hbadd:
                        hbadd_result = _run_command(hbadd_cmd, tmp_dir)
                        report["hbadd_rc"] = int(hbadd_result.returncode)
                        if hbadd_result.returncode != 0:
                            msg = _first_line(hbadd_result.stderr) or _first_line(hbadd_result.stdout)
                            report["errors"].append(
                                f"HBADD exit {hbadd_result.returncode}{': ' + msg if msg else ''}"
                            )

                    hbplus_nb_cmd = [
                        str(self.paths.hbplus_bin),
                        "-L",
                        "-h",
                        self.options.hbplus_nb_hparam,
                        "-d",
                        self.options.hbplus_nb_dparam,
                        "-N",
                        str(work_pdb),
                        "-wkdir",
                        str(tmp_dir) + os.sep,
                    ]
                    hbplus_nb_result = _run_command(hbplus_nb_cmd, tmp_dir)
                    report["hbplus_nb_rc"] = int(hbplus_nb_result.returncode)
                    if hbplus_nb_result.returncode != 0:
                        msg = _first_line(hbplus_nb_result.stderr) or _first_line(hbplus_nb_result.stdout)
                        report["errors"].append(
                            f"HBPLUS(nonbond) exit {hbplus_nb_result.returncode}{': ' + msg if msg else ''}"
                        )

                    hbplus_cmd = [
                        str(self.paths.hbplus_bin),
                        "-L",
                        "-h",
                        self.options.hbplus_hparam,
                        "-d",
                        self.options.hbplus_dparam,
                        str(work_pdb),
                        "-wkdir",
                        str(tmp_dir) + os.sep,
                    ]
                    hbplus_result = _run_command(hbplus_cmd, tmp_dir)
                    report["hbplus_rc"] = int(hbplus_result.returncode)
                    if hbplus_result.returncode != 0:
                        msg = _first_line(hbplus_result.stderr) or _first_line(hbplus_result.stdout)
                        report["errors"].append(
                            f"HBPLUS exit {hbplus_result.returncode}{': ' + msg if msg else ''}"
                        )

                ligplot_cmd = [
                    str(self.paths.ligplot_bin),
                    str(work_pdb),
                    resname_arg_start,
                    resnum_start,
                    resname_arg_end,
                    resnum_end,
                    chain,
                    "-wkdir",
                    str(tmp_dir) + os.sep,
                    "-prm",
                    str(ligplot_prm),
                    "-ctype",
                    str(self.options.contact_type),
                ]

                if self.options.plot_metal:
                    ligplot_cmd.append("-m")
                if self.options.filter_waters:
                    ligplot_cmd.append("-wcut")
                elif self.options.include_waters:
                    ligplot_cmd.append("-wat")
                if self.options.no_abort:
                    ligplot_cmd.append("-no_abort")

                ligplot_result = _run_command(ligplot_cmd, tmp_dir)
                report["ligplot_rc"] = int(ligplot_result.returncode)
                if ligplot_result.returncode != 0:
                    msg = _first_line(ligplot_result.stderr) or _first_line(ligplot_result.stdout)
                    report["errors"].append(
                        f"LIGPLOT exit {ligplot_result.returncode}{': ' + msg if msg else ''}"
                    )

                report["generated_outputs"] = sorted(path.name for path in tmp_dir.glob("ligplot.*"))
                report["output_ok"] = bool(_ligplot_outputs_ok(tmp_dir))
                if not report["output_ok"]:
                    report["errors"].append("Missing required LigPlot outputs (ligplot.ps and ligplot.drw)")
                return report

            hydrogenate_mode = self.options.hydrogenate
            strip_mode = self.options.strip_metals
            has_metals = _has_metal_atoms(orig_pdb)
            base_metal_mode = "strip" if strip_mode == "always" else "none"
            base_hydrogenate = hydrogenate_mode == "obabel"

            attempted = set()

            def try_attempt(hydrogenate: bool, metal_mode: str) -> bool:
                key = (hydrogenate, metal_mode)
                if key in attempted:
                    return False
                attempted.add(key)
                report = attempt_run(hydrogenate, metal_mode)
                attempt_reports.append(report)
                if report.get("output_ok"):
                    return True

                detail = "; ".join(str(msg) for msg in report.get("errors", [])[:3])
                if detail:
                    logger.info(
                        "LigPlot attempt failed for %s [hydrogenate=%s, metals=%s]: %s",
                        pdb_file.name,
                        hydrogenate,
                        metal_mode,
                        detail,
                    )
                else:
                    logger.info(
                        "LigPlot attempt failed for %s [hydrogenate=%s, metals=%s]",
                        pdb_file.name,
                        hydrogenate,
                        metal_mode,
                    )
                return False

            success = try_attempt(base_hydrogenate, base_metal_mode)

            if not success and hydrogenate_mode in {"auto", "obabel"} and base_hydrogenate:
                logger.info("Retrying LigPlot without hydrogenation for %s", pdb_file.name)
                success = try_attempt(False, base_metal_mode)

            if not success and strip_mode == "auto" and has_metals:
                logger.info("Retrying LigPlot with metal alias for %s", pdb_file.name)
                success = try_attempt(False, "alias")
                if not success:
                    logger.info("Retrying LigPlot with metal removal for %s", pdb_file.name)
                    success = try_attempt(False, "strip")

            if not success and hydrogenate_mode == "auto" and shutil.which("obabel"):
                logger.info("Retrying LigPlot with OpenBabel hydrogenation for %s", pdb_file.name)
                success = try_attempt(True, base_metal_mode)
                if not success and strip_mode == "auto" and has_metals:
                    logger.info("Retrying LigPlot with OpenBabel + metal alias for %s", pdb_file.name)
                    success = try_attempt(True, "alias")
                    if not success:
                        logger.info("Retrying LigPlot with OpenBabel + metal removal for %s", pdb_file.name)
                        success = try_attempt(True, "strip")

            if not success:
                reason_items: List[str] = []
                for idx, report in enumerate(attempt_reports, start=1):
                    errors = report.get("errors", [])
                    first_error = str(errors[0]) if errors else "Missing required outputs"
                    reason_items.append(
                        f"attempt#{idx}[hydrogenate={report.get('hydrogenate')}, metals={report.get('metal_mode')}] {first_error}"
                    )
                failure_reason = " | ".join(reason_items[:6])
                if len(reason_items) > 6:
                    failure_reason += f" | ... (+{len(reason_items) - 6} more attempts)"

                logger.warning(
                    "LigPlot failed for %s after %d attempts: %s",
                    pdb_file.name,
                    len(attempt_reports),
                    failure_reason or "unknown failure",
                )
                self.last_run_report = {
                    "pdb_file": str(pdb_file),
                    "success": False,
                    "skipped_existing": False,
                    "failure_reason": failure_reason or "unknown failure",
                    "attempts": attempt_reports,
                }
                return False

            outputs = list(tmp_dir.glob("ligplot.*"))
            if not outputs:
                reason = "No LigPlot output files produced"
                logger.warning("%s for %s", reason, pdb_file.name)
                self.last_run_report = {
                    "pdb_file": str(pdb_file),
                    "success": False,
                    "skipped_existing": False,
                    "failure_reason": reason,
                    "attempts": attempt_reports,
                }
                return False

            copied_outputs: List[str] = []
            for out_file in outputs:
                shutil.copy2(out_file, output_dir / out_file.name)
                copied_outputs.append(out_file.name)

            if len(attempt_reports) > 1:
                final = attempt_reports[-1]
                logger.info(
                    "LigPlot recovered for %s after %d attempts [hydrogenate=%s, metals=%s]",
                    pdb_file.name,
                    len(attempt_reports),
                    final.get("hydrogenate"),
                    final.get("metal_mode"),
                )

            self.last_run_report = {
                "pdb_file": str(pdb_file),
                "success": True,
                "skipped_existing": False,
                "failure_reason": "",
                "attempts": attempt_reports,
                "copied_outputs": copied_outputs,
            }
            return True


def list_ligands_from_pairlist(pairlist_file: Path) -> List[str]:
    ligands: List[str] = []
    if not pairlist_file.exists():
        return ligands
    with pairlist_file.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            ligand = (row.get("ligand") or "").strip()
            if ligand:
                ligands.append(ligand)
    return sorted(set(ligands))


def iter_pdb_files(input_dir: Path, recursive: bool = True) -> List[Path]:
    if recursive:
        return sorted(input_dir.rglob("*.pdb"))
    return sorted(input_dir.glob("*.pdb"))


def export_ligplot_outputs(
    output_dir: Path,
    formats: Iterable[str],
    overwrite: bool = False,
    dpi: int = 300,
) -> bool:
    output_dir = output_dir.resolve()
    ps_file = output_dir / "ligplot.ps"
    if not ps_file.exists():
        logger.warning("No ligplot.ps found in %s", output_dir)
        return False

    requested = {fmt.lower() for fmt in formats if fmt}
    if not requested or requested == {"none"}:
        return True

    success = True

    def ensure_pdf() -> Optional[Path]:
        pdf_file = output_dir / "ligplot.pdf"
        if pdf_file.exists() and not overwrite:
            return pdf_file

        if shutil.which("pstopdf"):
            cmd = ["pstopdf", str(ps_file), "-o", str(pdf_file)]
            result = _run_tool(cmd)
        elif shutil.which("ps2pdf"):
            cmd = ["ps2pdf", str(ps_file), str(pdf_file)]
            result = _run_tool(cmd)
        elif shutil.which("gs"):
            cmd = [
                "gs",
                "-q",
                "-dNOPAUSE",
                "-dBATCH",
                "-sDEVICE=pdfwrite",
                f"-sOutputFile={pdf_file}",
                str(ps_file),
            ]
            result = _run_tool(cmd)
        else:
            logger.warning("No PS->PDF converter found (pstopdf/ps2pdf/gs). Install Ghostscript for export.")
            return None

        if result.returncode != 0:
            logger.warning(
                "PS->PDF conversion failed: %s%s",
                result.stdout.strip(),
                f" | {result.stderr.strip()}" if result.stderr.strip() else "",
            )
            return None
        return pdf_file

    if "pdf" in requested or "both" in requested:
        pdf_file = ensure_pdf()
        if pdf_file is None:
            success = False

    if "png" in requested or "both" in requested:
        png_file = output_dir / "ligplot.png"
        if png_file.exists() and not overwrite:
            return success

        if shutil.which("gs"):
            cmd = [
                "gs",
                "-q",
                "-dNOPAUSE",
                "-dBATCH",
                "-dTextAlphaBits=4",
                "-dGraphicsAlphaBits=4",
                "-sDEVICE=pngalpha",
                f"-r{dpi}",
                f"-sOutputFile={png_file}",
                str(ps_file),
            ]
            result = _run_tool(cmd)
        elif shutil.which("magick"):
            cmd = [
                "magick",
                "convert",
                "-density",
                str(dpi),
                str(ps_file),
                str(png_file),
            ]
            result = _run_tool(cmd)
        elif shutil.which("convert"):
            cmd = [
                "convert",
                "-density",
                str(dpi),
                str(ps_file),
                str(png_file),
            ]
            result = _run_tool(cmd)
        elif shutil.which("sips"):
            # Try direct PS -> PNG conversion via sips first.
            cmd = ["sips", "-s", "format", "png", str(ps_file), "--out", str(png_file)]
            result = _run_tool(cmd)
            if result.returncode != 0:
                pdf_file = ensure_pdf()
                if pdf_file is None:
                    return False
                cmd = ["sips", "-s", "format", "png", str(pdf_file), "--out", str(png_file)]
                result = _run_tool(cmd)
        else:
            logger.warning("No PS->PNG converter found (gs/magick/convert/sips). Install Ghostscript for export.")
            return False

        if result.returncode != 0:
            logger.warning(
                "PS->PNG conversion failed: %s%s",
                result.stdout.strip(),
                f" | {result.stderr.strip()}" if result.stderr.strip() else "",
            )
            success = False

    return success
