from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import pandas as pd

from post_docking_analysis.docking_parser import parse_vina_pdbqt
from ..models import PairlistRow
from .base import DockingEngineRunner, build_pair_index, split_replicate_stem


class VinaRunner(DockingEngineRunner):
    name = "vina"
    binary_name = "vina"
    pose_extension = ".pdbqt"

    def build_command(
        self,
        row: PairlistRow,
        pose_file: Path,
        log_file: Path,
        seed: Optional[int] = None,
    ) -> List[str]:
        runtime = self.runtime
        env_name = str(runtime.get("conda_env") or "").strip()
        binary = str(runtime.get("binary") or self.binary_name)
        exhaustiveness = int(runtime.get("exhaustiveness", 32))
        num_modes = int(runtime.get("num_modes", 20))
        energy_range = float(runtime.get("energy_range", 3.0))
        cpu = runtime.get("cpu")
        if seed is None:
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
                "--energy_range",
                f"{energy_range:g}",
                "--out",
                str(pose_file),
            ]
        )
        if cpu not in (None, ""):
            command.extend(["--cpu", str(cpu)])
        if seed not in (None, ""):
            command.extend(["--seed", str(seed)])
        return command

    def count_returned_poses(self, pose_file: Path) -> Optional[int]:
        count = 0
        with open(pose_file, "r", encoding="utf-8", errors="replace") as handle:
            for raw_line in handle:
                if raw_line.startswith("MODEL"):
                    count += 1
        return count

    def collect_normalized_scores(self, pairlist_rows: List[PairlistRow]) -> pd.DataFrame:
        pair_index = build_pair_index(pairlist_rows)
        rows = []
        for pose_file in sorted(self.layout["poses"].glob(f"*{self.pose_extension}")):
            pair_tag, replicate_id = split_replicate_stem(pose_file.stem)
            parsed = parse_vina_pdbqt(pose_file)
            pair = pair_index.get(pair_tag)
            for _, record in parsed.iterrows():
                rows.append(
                    {
                        "engine": self.name,
                        "tag": pair_tag,
                        "replicate_id": replicate_id,
                        "seed": self.seed_for_replicate(replicate_id) if replicate_id else None,
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
                        "log_file": str(self.layout["logs"] / f"{pose_file.stem}.log"),
                    }
                )
        return pd.DataFrame(rows)
