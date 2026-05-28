from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import pandas as pd

from post_docking_analysis.docking_parser import parse_vina_pdbqt
from ..models import PairlistRow
from .base import DockingEngineRunner


class VinaRunner(DockingEngineRunner):
    name = "vina"
    binary_name = "vina"
    pose_extension = ".pdbqt"

    def build_command(self, row: PairlistRow, pose_file: Path, log_file: Path) -> List[str]:
        runtime = self.runtime
        env_name = str(runtime.get("conda_env") or "").strip()
        binary = str(runtime.get("binary") or self.binary_name)
        exhaustiveness = int(runtime.get("exhaustiveness", 16))
        num_modes = int(runtime.get("num_modes", 20))
        cpu = runtime.get("cpu")
        seed = runtime.get("seed")

        command: List[str] = []
        if env_name:
            command.extend(["conda", "run", "-n", env_name])
        command.extend(
            [
                binary,
                "--receptor",
                str(self.shared_receptors / row.receptor),
                "--ligand",
                str(self.shared_ligands / row.ligand),
                "--center_x",
                str(row.center_x),
                "--center_y",
                str(row.center_y),
                "--center_z",
                str(row.center_z),
                "--size_x",
                str(row.size_x),
                "--size_y",
                str(row.size_y),
                "--size_z",
                str(row.size_z),
                "--exhaustiveness",
                str(exhaustiveness),
                "--num_modes",
                str(num_modes),
                "--out",
                str(pose_file),
            ]
        )
        if cpu not in (None, ""):
            command.extend(["--cpu", str(cpu)])
        if seed not in (None, ""):
            command.extend(["--seed", str(seed)])
        return command

    def collect_normalized_scores(self, pairlist_rows: List[PairlistRow]) -> pd.DataFrame:
        pair_index: Dict[str, PairlistRow] = {row.tag: row for row in pairlist_rows}
        rows = []
        for pose_file in sorted(self.layout["poses"].glob(f"*{self.pose_extension}")):
            tag = pose_file.stem
            parsed = parse_vina_pdbqt(pose_file)
            pair = pair_index.get(tag)
            for _, record in parsed.iterrows():
                rows.append(
                    {
                        "engine": self.name,
                        "tag": tag,
                        "protein": pair.receptor if pair else "",
                        "ligand": pair.ligand if pair else "",
                        "site_id": pair.site_id if pair else "",
                        "pose": int(record.get("pose", 0)),
                        "affinity_kcal_mol": float(record.get("vina_affinity")),
                        "score_name_primary": "vina_affinity",
                        "score_primary": float(record.get("vina_affinity")),
                        "score_name_secondary": "",
                        "score_secondary": None,
                        "rmsd_lb": float(record.get("rmsd_lb")) if pd.notna(record.get("rmsd_lb")) else None,
                        "rmsd_ub": float(record.get("rmsd_ub")) if pd.notna(record.get("rmsd_ub")) else None,
                        "pose_file": str(pose_file),
                        "log_file": str(self.layout["logs"] / f"{tag}.log"),
                    }
                )
        return pd.DataFrame(rows)
