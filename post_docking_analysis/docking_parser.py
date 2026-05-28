"""
Docking result parser for post-docking analysis pipeline.

This module handles parsing of PDBQT files to extract binding affinity and RMSD
values from Vina-family engines. Vina emits ``REMARK VINA RESULT`` while smina
often emits ``REMARK minimizedAffinity`` and ``REMARK minimizedRMSD`` per model.
"""

from pathlib import Path
import re
from typing import Dict, List, Optional

import pandas as pd


def parse_autodock4_dlg(dlg_file: Path) -> pd.DataFrame:
    """Parse an AutoDock4 DLG file and extract per-pose binding energies."""
    results: List[Dict[str, float]] = []
    current_pose: Optional[int] = None
    current_run: Optional[int] = None
    affinity_pattern = re.compile(r"=\s*([+-]?\d+(?:\.\d+)?)\s*kcal/mol", re.IGNORECASE)

    with open(dlg_file, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith("DOCKED: MODEL"):
                parts = line.split()
                try:
                    current_pose = int(parts[-1])
                except (IndexError, ValueError):
                    current_pose = None
                continue

            if line.startswith("DOCKED: USER") and "Run =" in line:
                try:
                    current_run = int(line.split("Run =", 1)[1].strip().split()[0])
                except (IndexError, ValueError):
                    current_run = None
                continue

            if "Estimated Free Energy of Binding" not in line:
                continue

            match = affinity_pattern.search(line)
            if not match:
                continue
            try:
                affinity = float(match.group(1))
            except ValueError:
                continue

            pose_number = current_pose or current_run or (len(results) + 1)
            results.append(
                {
                    "pose": int(pose_number),
                    "autodock4_affinity": affinity,
                }
            )

    return pd.DataFrame(results)

def parse_vina_pdbqt(pdbqt_file: Path) -> pd.DataFrame:
    """
    Parse a Vina PDBQT file and extract binding affinity and RMSD values.
    
    Parameters
    ----------
    pdbqt_file : Path
        Path to the PDBQT file
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing pose information
    """
    results: List[Dict[str, float]] = []
    current_pose: Optional[int] = None
    current_affinity: Optional[float] = None
    current_rmsd: Optional[float] = None

    def flush_current() -> None:
        nonlocal current_affinity, current_rmsd, current_pose
        if current_affinity is None:
            return
        pose_number = current_pose or (len(results) + 1)
        rmsd_value = current_rmsd if current_rmsd is not None else 0.0
        results.append(
            {
                "pose": pose_number,
                "vina_affinity": current_affinity,
                "rmsd_lb": rmsd_value,
                "rmsd_ub": rmsd_value,
            }
        )
        current_affinity = None
        current_rmsd = None

    with open(pdbqt_file, "r", encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith("MODEL "):
                flush_current()
                parts = line.split()
                try:
                    current_pose = int(parts[1])
                except (IndexError, ValueError):
                    current_pose = len(results) + 1
                continue

            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        results.append(
                            {
                                "pose": current_pose or (len(results) + 1),
                                "vina_affinity": float(parts[3]),
                                "rmsd_lb": float(parts[4]),
                                "rmsd_ub": float(parts[5]),
                            }
                        )
                    except ValueError:
                        continue
                continue

            if line.startswith("REMARK minimizedAffinity"):
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        current_affinity = float(parts[2])
                    except ValueError:
                        current_affinity = None
                continue

            if line.startswith("REMARK minimizedRMSD"):
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        current_rmsd = float(parts[2])
                    except ValueError:
                        current_rmsd = None
                continue

            if line.startswith("ENDMDL"):
                flush_current()

    flush_current()
    return pd.DataFrame(results)

def parse_all_docking_results(complexes: List[Dict[str, Path]]) -> Dict[str, pd.DataFrame]:
    """
    Parse docking results for all complexes.
    
    Parameters
    ----------
    complexes : List[Dict[str, Path]]
        List of complexes with docking result files
        
    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary mapping complex names to their parsed results
    """
    all_results = {}
    
    for complex_info in complexes:
        complex_name = complex_info["name"]
        if "docking_result" in complex_info:
            try:
                df = parse_vina_pdbqt(complex_info["docking_result"])
                if not df.empty:
                    all_results[complex_name] = df
                else:
                    print(f"⚠️  No poses found in {complex_info['docking_result']}")
            except Exception as e:
                print(f"❌ Error parsing {complex_info['docking_result']}: {e}")
        else:
            print(f"⚠️  No docking result file for {complex_name}")
    
    return all_results
