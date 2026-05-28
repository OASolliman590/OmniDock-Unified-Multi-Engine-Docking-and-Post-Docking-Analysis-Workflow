from __future__ import annotations

import os
import shlex
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from post_docking_analysis.docking_parser import parse_autodock4_dlg
from ..models import PairlistRow
from .base import DockingEngineRunner


class AutoDock4Runner(DockingEngineRunner):
    name = "autodock4"
    binary_name = "autodock4"
    pose_extension = ".dlg"

    def _resolve_ligand_path(self, ligand_name: str) -> Path:
        candidates = [
            self.project_root / "4-Docking" / "ligands_engine_specific" / "autodock4" / ligand_name,
            self.project_root / "ligands_engine_specific" / "autodock4" / ligand_name,
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return self.shared_ligands / ligand_name

    def _resolve_receptor_path(self, receptor_name: str) -> Path:
        candidates = [
            self.project_root / "4-Docking" / "receptors_engine_specific" / "autodock4" / receptor_name,
            self.project_root / "receptors_engine_specific" / "autodock4" / receptor_name,
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return self.shared_receptors / receptor_name

    @staticmethod
    def _extract_atom_types(pdbqt_file: Path) -> List[str]:
        atom_types = set()
        with open(pdbqt_file, "r", encoding="utf-8", errors="replace") as handle:
            for raw_line in handle:
                if not raw_line.startswith(("ATOM", "HETATM")):
                    continue
                atom_type = raw_line[77:79].strip() if len(raw_line) >= 79 else raw_line.split()[-1].strip()
                if atom_type:
                    atom_types.add(atom_type)
        return sorted(atom_types)

    @staticmethod
    def _estimate_about_center(ligand_file: Path) -> Tuple[float, float, float]:
        xs: List[float] = []
        ys: List[float] = []
        zs: List[float] = []
        with open(ligand_file, "r", encoding="utf-8", errors="replace") as handle:
            for raw_line in handle:
                if not raw_line.startswith(("ATOM", "HETATM")):
                    continue
                try:
                    xs.append(float(raw_line[30:38].strip()))
                    ys.append(float(raw_line[38:46].strip()))
                    zs.append(float(raw_line[46:54].strip()))
                except ValueError:
                    continue
        if not xs:
            return 0.0, 0.0, 0.0
        return sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)

    @staticmethod
    def _grid_points(size: float, spacing: float) -> int:
        points = int(round(float(size) / float(spacing)))
        return max(2, points)

    @staticmethod
    def _extract_torsdof(pdbqt_file: Path) -> Optional[int]:
        try:
            with open(pdbqt_file, "r", encoding="utf-8", errors="replace") as handle:
                for raw_line in handle:
                    if not raw_line.startswith("TORSDOF"):
                        continue
                    tokens = raw_line.strip().split()
                    if len(tokens) < 2:
                        continue
                    return int(float(tokens[1]))
        except Exception:
            return None
        return None

    def _resolve_adt_param_generators(self) -> Optional[Dict[str, str]]:
        py_exec = str(self.runtime.get("autodocktools_python") or "").strip()
        gpf_script = str(self.runtime.get("autodocktools_prepare_gpf4") or "").strip()
        dpf_script = str(self.runtime.get("autodocktools_prepare_dpf4") or "").strip()

        if not (gpf_script and dpf_script):
            ligand_script = str(
                self.runtime.get("autodocktools_prepare_ligand4")
                or self.runtime.get("autodocktools_prepare_receptor4")
                or ""
            ).strip()
            if ligand_script:
                parent = Path(ligand_script).expanduser().resolve().parent
                if not gpf_script:
                    candidate = parent / "prepare_gpf4.py"
                    if candidate.exists():
                        gpf_script = str(candidate)
                if not dpf_script:
                    candidate = parent / "prepare_dpf4.py"
                    if candidate.exists():
                        dpf_script = str(candidate)

        if not py_exec:
            py_exec = os.environ.get("MGLTOOLS_PYTHON", "").strip()

        if not (py_exec and gpf_script and dpf_script):
            return None
        return {
            "python": py_exec,
            "prepare_gpf4": gpf_script,
            "prepare_dpf4": dpf_script,
        }

    def _resolve_parameter_file(self) -> str:
        return str(
            self.runtime.get("parameter_file_path")
            or self.runtime.get("parameter_file")
            or "AD4.1_bound.dat"
        )

    def _write_gpf(
        self,
        gpf_file: Path,
        receptor_file: Path,
        ligand_file: Path,
        map_prefix: Path,
        center: Tuple[float, float, float],
        box_size: Tuple[float, float, float],
        spacing: float,
    ) -> List[str]:
        receptor_types = self._extract_atom_types(receptor_file)
        ligand_types = self._extract_atom_types(ligand_file)
        if not receptor_types:
            raise ValueError(f"Could not determine receptor atom types from {receptor_file}")
        if not ligand_types:
            raise ValueError(f"Could not determine ligand atom types from {ligand_file}")
        npts = tuple(self._grid_points(size, spacing) for size in box_size)

        lines = [
            f"npts {npts[0]} {npts[1]} {npts[2]}",
            f"gridfld {map_prefix}.maps.fld",
            f"spacing {spacing:.3f}",
            f"receptor_types {' '.join(receptor_types)}",
            f"ligand_types {' '.join(ligand_types)}",
            f"receptor {receptor_file}",
            f"gridcenter {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}",
            "smooth 0.5",
        ]
        for atom_type in ligand_types:
            lines.append(f"map {map_prefix}.{atom_type}.map")
        lines.extend(
            [
                f"elecmap {map_prefix}.e.map",
                f"dsolvmap {map_prefix}.d.map",
                "dielectric -0.1465",
            ]
        )
        gpf_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return ligand_types

    def _write_dpf(
        self,
        dpf_file: Path,
        map_prefix: Path,
        ligand_file: Path,
        ligand_types: List[str],
        about: Tuple[float, float, float],
    ) -> None:
        ga_pop_size = int(self.runtime.get("ga_pop_size", 150))
        ga_num_evals = int(self.runtime.get("ga_num_evals", 25000000))
        ga_num_generations = int(self.runtime.get("ga_num_generations", 27000))
        ga_run = int(self.runtime.get("ga_run", 100))
        ls_search_freq = float(self.runtime.get("ls_search_freq", 0.06))
        torsdof_override = int(self.runtime.get("torsdof", 0))
        torsdof = torsdof_override if torsdof_override > 0 else (self._extract_torsdof(ligand_file) or 0)
        parameter_file = self._resolve_parameter_file()
        seed = self.runtime.get("seed")
        if seed in (None, ""):
            seed_line = "seed pid time"
        else:
            seed_int = int(seed)
            seed_line = f"seed {seed_int} {seed_int}"

        lines = [
            "autodock_parameter_version 4.2",
            f"parameter_file {parameter_file}",
            "intelec",
            seed_line,
            f"ligand_types {' '.join(ligand_types)}",
            f"fld {map_prefix}.maps.fld",
        ]
        for atom_type in ligand_types:
            lines.append(f"map {map_prefix}.{atom_type}.map")
        lines.extend(
            [
                f"elecmap {map_prefix}.e.map",
                f"desolvmap {map_prefix}.d.map",
                f"move {ligand_file}",
                f"about {about[0]:.3f} {about[1]:.3f} {about[2]:.3f}",
                "tran0 random",
                "quaternion0 random",
                "dihe0 random",
                f"torsdof {torsdof}",
                "rmstol 2.0",
                "extnrg 1000.0",
                "e0max 0.0 10000",
                f"ga_pop_size {ga_pop_size}",
                f"ga_num_evals {ga_num_evals}",
                f"ga_num_generations {ga_num_generations}",
                "ga_elitism 1",
                "ga_mutation_rate 0.02",
                "ga_crossover_rate 0.8",
                "ga_window_size 10",
                "ga_cauchy_alpha 0.0",
                "ga_cauchy_beta 1.0",
                "set_ga",
                "sw_max_its 300",
                "sw_max_succ 4",
                "sw_max_fail 4",
                "sw_rho 1.0",
                "sw_lb_rho 0.01",
                f"ls_search_freq {ls_search_freq}",
                "set_psw1",
                "unbound_model bound",
                f"ga_run {ga_run}",
                "analysis",
            ]
        )
        dpf_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def build_command(self, row: PairlistRow, pose_file: Path, log_file: Path) -> List[str]:
        autogrid_binary = str(self.runtime.get("autogrid_binary") or "autogrid4")
        autodock_binary = str(self.runtime.get("binary") or self.runtime.get("autodock_binary") or self.binary_name)
        spacing = float(self.runtime.get("spacing", 0.375))

        receptor_file = self._resolve_receptor_path(row.receptor)
        ligand_file = self._resolve_ligand_path(row.ligand)

        pair_dir = self.layout["root"] / "pairs" / row.tag
        pair_dir.mkdir(parents=True, exist_ok=True)
        map_prefix = pair_dir / "maps"
        gpf_file = pair_dir / "grid.gpf"
        glg_file = pair_dir / "grid.glg"
        dpf_file = pair_dir / "docking.dpf"
        scripts = self._resolve_adt_param_generators()
        if scripts:
            npts = tuple(self._grid_points(size, spacing) for size in (row.size_x, row.size_y, row.size_z))
            receptor_local = receptor_file.name
            ligand_local = ligand_file.name
            parameter_value = self._resolve_parameter_file()
            parameter_local = parameter_value
            parameter_link = ""
            if "/" in parameter_value:
                parameter_local = Path(parameter_value).name
                parameter_link = f"ln -sf {shlex.quote(parameter_value)} {shlex.quote(parameter_local)}; "
            gpf_cmd = " ".join(
                [
                    shlex.quote(scripts["python"]),
                    shlex.quote(scripts["prepare_gpf4"]),
                    "-l",
                    shlex.quote(ligand_local),
                    "-r",
                    shlex.quote(receptor_local),
                    "-o",
                    shlex.quote(str(gpf_file.name)),
                    "-p",
                    shlex.quote(f"npts={npts[0]},{npts[1]},{npts[2]}"),
                    "-p",
                    shlex.quote(f"gridcenter={row.center_x:.3f},{row.center_y:.3f},{row.center_z:.3f}"),
                    "-p",
                    shlex.quote(f"spacing={spacing:.3f}"),
                ]
            )
            dpf_parts: List[str] = [
                shlex.quote(scripts["python"]),
                shlex.quote(scripts["prepare_dpf4"]),
                "-l",
                shlex.quote(ligand_local),
                "-r",
                shlex.quote(receptor_local),
                "-o",
                shlex.quote(str(dpf_file.name)),
                "-p",
                shlex.quote(f"parameter_file={parameter_local}"),
                "-p",
                shlex.quote(f"ga_pop_size={int(self.runtime.get('ga_pop_size', 150))}"),
                "-p",
                shlex.quote(f"ga_num_evals={int(self.runtime.get('ga_num_evals', 25000000))}"),
                "-p",
                shlex.quote(f"ga_num_generations={int(self.runtime.get('ga_num_generations', 27000))}"),
                "-p",
                shlex.quote(f"ga_run={int(self.runtime.get('ga_run', 100))}"),
                "-p",
                shlex.quote(f"ls_search_freq={float(self.runtime.get('ls_search_freq', 0.06))}"),
            ]
            torsdof_override = int(self.runtime.get("torsdof", 0))
            if torsdof_override > 0:
                dpf_parts.extend(["-p", shlex.quote(f"torsdof={torsdof_override}")])
            seed = self.runtime.get("seed")
            if seed not in (None, ""):
                seed_int = int(seed)
                dpf_parts.extend(["-p", shlex.quote(f"seed={seed_int},{seed_int}")])
            dpf_cmd = " ".join(dpf_parts)
            shell_cmd = (
                "set -euo pipefail; "
                f"cd {shlex.quote(str(pair_dir))}; "
                f"ln -sf {shlex.quote(str(receptor_file))} {shlex.quote(receptor_local)}; "
                f"ln -sf {shlex.quote(str(ligand_file))} {shlex.quote(ligand_local)}; "
                f"{parameter_link}"
                f"{gpf_cmd}; "
                f"{dpf_cmd}; "
                f"{shlex.quote(autogrid_binary)} -p {shlex.quote(str(gpf_file.name))} -l {shlex.quote(str(glg_file.name))}; "
                f"{shlex.quote(autodock_binary)} -p {shlex.quote(str(dpf_file.name))} -l {shlex.quote(str(pose_file))}"
            )
        else:
            ligand_types = self._write_gpf(
                gpf_file=gpf_file,
                receptor_file=receptor_file,
                ligand_file=ligand_file,
                map_prefix=map_prefix,
                center=(row.center_x, row.center_y, row.center_z),
                box_size=(row.size_x, row.size_y, row.size_z),
                spacing=spacing,
            )
            about = self._estimate_about_center(ligand_file)
            self._write_dpf(
                dpf_file=dpf_file,
                map_prefix=map_prefix,
                ligand_file=ligand_file,
                ligand_types=ligand_types,
                about=about,
            )
            shell_cmd = (
                "set -euo pipefail; "
                f"{shlex.quote(autogrid_binary)} -p {shlex.quote(str(gpf_file))} -l {shlex.quote(str(glg_file))}; "
                f"{shlex.quote(autodock_binary)} -p {shlex.quote(str(dpf_file))} -l {shlex.quote(str(pose_file))}"
            )
        return ["bash", "-lc", shell_cmd]

    def collect_normalized_scores(self, pairlist_rows: List[PairlistRow]) -> pd.DataFrame:
        pair_index: Dict[str, PairlistRow] = {row.tag: row for row in pairlist_rows}
        rows = []
        for dlg_file in sorted(self.layout["poses"].glob(f"*{self.pose_extension}")):
            tag = dlg_file.stem
            parsed = parse_autodock4_dlg(dlg_file)
            pair = pair_index.get(tag)
            for _, record in parsed.iterrows():
                affinity = float(record.get("autodock4_affinity"))
                rows.append(
                    {
                        "engine": self.name,
                        "tag": tag,
                        "protein": pair.receptor if pair else "",
                        "ligand": pair.ligand if pair else "",
                        "site_id": pair.site_id if pair else "",
                        "pose": int(record.get("pose", 0)),
                        "affinity_kcal_mol": affinity,
                        "score_name_primary": "autodock4_affinity",
                        "score_primary": affinity,
                        "score_name_secondary": "",
                        "score_secondary": None,
                        "rmsd_lb": None,
                        "rmsd_ub": None,
                        "pose_file": str(dlg_file),
                        "log_file": str(self.layout["logs"] / f"{tag}.log"),
                    }
                )
        return pd.DataFrame(rows)
