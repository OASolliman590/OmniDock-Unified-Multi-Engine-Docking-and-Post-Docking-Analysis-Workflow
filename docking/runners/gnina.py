from __future__ import annotations

from pathlib import Path
from typing import List

import pandas as pd

from post_docking_analysis.generate_scores_csv import generate_all_scores_csv
from typing import Optional

from ..models import PairlistRow
from .base import DockingEngineRunner, build_pair_index, split_replicate_stem


class GninaRunner(DockingEngineRunner):
    name = "gnina"
    binary_name = "gnina"
    pose_extension = ".sdf"

    def build_command(
        self,
        row: PairlistRow,
        pose_file: Path,
        log_file: Path,
        seed: Optional[int] = None,
    ) -> List[str]:
        runtime = self.runtime
        image = str(runtime.get("image") or "")
        binary = str(runtime.get("binary") or self.binary_name)
        raw_use_gpu = runtime.get("use_gpu")
        use_gpu = True if raw_use_gpu in (None, "") else bool(raw_use_gpu)
        score_only = bool(runtime.get("score_only", False))
        exhaustiveness = int(runtime.get("exhaustiveness", 32))
        num_modes = int(runtime.get("num_modes", 20))
        cnn_scoring = str(runtime.get("cnn_scoring", "rescore"))
        device = runtime.get("device")
        cpu = runtime.get("cpu")
        if seed is None:
            seed = runtime.get("seed")

        command: List[str] = []
        if image:
            command = ["apptainer", "exec"]
            if use_gpu:
                command.append("--nv")
            command.extend([image, binary])
        else:
            command = [binary]

        command.extend(
            [
                "-r",
                str(self.shared_receptors / row.receptor),
                "-l",
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
                "--cnn_scoring",
                cnn_scoring,
                "--out",
                str(pose_file),
                "--log",
                str(log_file),
            ]
        )
        if score_only:
            command.append("--score_only")
        if seed not in (None, ""):
            command.extend(["--seed", str(seed)])
        if use_gpu and device not in (None, ""):
            command.extend(["--device", str(device)])
        elif cpu not in (None, ""):
            command.extend(["--cpu", str(cpu)])
        return command

    def collect_normalized_scores(self, pairlist_rows: List[PairlistRow]) -> pd.DataFrame:
        all_scores = self.layout["scores"] / "all_scores.csv"
        if not generate_all_scores_csv(
            self.layout["poses"],
            output_file=all_scores,
            pairlist_file=self.pairlist_file,
            log_dir=self.layout["logs"],
        ):
            return pd.DataFrame()

        df = pd.read_csv(all_scores)
        pair_index = build_pair_index(pairlist_rows)
        rows = []
        for _, record in df.iterrows():
            pose_stem = str(record.get("tag", ""))
            tag, replicate_id = split_replicate_stem(pose_stem)
            pair = pair_index.get(tag)
            rows.append(
                {
                    "engine": self.name,
                    "tag": tag,
                    "replicate_id": replicate_id,
                    "seed": self.seed_for_replicate(replicate_id),
                    "protein": pair.receptor if pair else "",
                    "ligand": pair.ligand if pair else "",
                    "site_id": pair.site_id if pair else "",
                    "pose": int(record.get("mode", 0)),
                    "affinity_kcal_mol": float(record.get("vina_affinity")),
                    "score_name_primary": "cnn_affinity",
                    "score_primary": float(record.get("cnn_affinity")) if pd.notna(record.get("cnn_affinity")) else None,
                    "score_name_secondary": "cnn_score",
                    "score_secondary": float(record.get("cnn_score")) if pd.notna(record.get("cnn_score")) else None,
                    "rmsd_lb": None,
                    "rmsd_ub": None,
                    "pose_file": str(self.layout["poses"] / f"{pose_stem}{self.pose_extension}"),
                    "log_file": str(self.layout["logs"] / f"{pose_stem}.log"),
                }
            )
        return pd.DataFrame(rows)
