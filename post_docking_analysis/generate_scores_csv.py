#!/usr/bin/env python3
"""
Script to generate all_scores.csv from GNINA log files.

This script parses GNINA log files to extract docking scores and creates
the required all_scores.csv file for the post-docking analysis pipeline.
It can use pairlist.csv for accurate complex naming if available.
"""

import os
import re
import csv
import argparse
import json
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple

def load_pairlist_mapping(pairlist_file: Path) -> Dict[str, str]:
    """
    Load pairlist mapping for accurate complex naming.
    
    Parameters
    ----------
    pairlist_file : Path
        Path to pairlist.csv file
        
    Returns
    -------
    dict
        Dictionary mapping log filenames to complex names
    """
    try:
        # Read the pairlist CSV
        df = pd.read_csv(pairlist_file)
        
        # Create mapping from receptor_site_ligand to tag names
        mapping = {}
        for _, row in df.iterrows():
            receptor = row['receptor']
            site_id = row['site_id']
            ligand = row['ligand']
            
            # Create expected log filename pattern
            log_pattern = f"{receptor}_{site_id}_{ligand}"
            tag_name = f"{receptor}_{site_id}_{ligand}"
            
            mapping[log_pattern] = tag_name
            
        logging.getLogger(__name__).info("Loaded %s mappings from pairlist.csv", len(mapping))
        return mapping
        
    except Exception as e:
        logging.getLogger(__name__).warning("Could not load pairlist.csv: %s", e)
        return {}


def _resolve_tag_from_pairlist(filename_stem: str, pairlist_mapping: Dict[str, str]) -> Tuple[str, bool]:
    """
    Resolve a log stem to a canonical tag using strict, deterministic matching.

    Returns
    -------
    Tuple[str, bool]
        Resolved tag and whether it came from pairlist mapping.
    """
    if not pairlist_mapping:
        return filename_stem, False
    mapped = pairlist_mapping.get(filename_stem)
    if mapped:
        return mapped, True
    return filename_stem, False


def parse_gnina_log(log_file: Path, pairlist_mapping: Optional[Dict[str, str]] = None):
    """
    Parse a GNINA log file to extract docking scores.
    
    Parameters
    ----------
    log_file : Path
        Path to the GNINA log file
    pairlist_mapping : dict, optional
        Mapping from log patterns to tag names
        
    Returns
    -------
    list
        List of dictionaries containing score information
    """
    scores: List[Dict[str, object]] = []
    
    # Extract tag name from filename
    filename = log_file.stem
    
    # Use pairlist mapping if available
    tag_name, matched_pairlist = _resolve_tag_from_pairlist(filename, pairlist_mapping or {})
    if log_file.name.lower().endswith(".runner.log") or (pairlist_mapping is not None and not matched_pairlist):
        return [], False
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        # Look for the scoring table
        # Pattern: mode |  affinity  |  intramol  |    CNN     |   CNN
        #          | (kcal/mol) | (kcal/mol) | pose score | affinity
        # -----+------------+------------+------------+----------
        #    1      -12.53       -0.67       0.9954      7.774
        
        lines = content.split('\n')
        in_scores_section = False
        mode_line_pattern = re.compile(r'^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')
        mode_line_pattern_3 = re.compile(r'^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')
        
        for line in lines:
            # Detect start of scores section
            if 'mode |' in line and 'affinity' in line:
                in_scores_section = True
                continue
                
            # Process score lines
            if in_scores_section:
                match = mode_line_pattern.match(line)
                if match:
                    mode = int(match.group(1))
                    vina_affinity = float(match.group(2))
                    intramol = float(match.group(3))
                    cnn_score = float(match.group(4))
                    cnn_affinity = float(match.group(5))
                    
                    scores.append({
                        'tag': tag_name,
                        'mode': mode,
                        'vina_affinity': vina_affinity,
                        'cnn_affinity': cnn_affinity,
                        'cnn_score': cnn_score
                    })
                    continue

                match3 = mode_line_pattern_3.match(line)
                if match3:
                    mode = int(match3.group(1))
                    vina_affinity = float(match3.group(2))
                    scores.append({
                        'tag': tag_name,
                        'mode': mode,
                        'vina_affinity': vina_affinity,
                        'cnn_affinity': None,
                        'cnn_score': None
                    })
                    
            # Do not terminate on non-indented lines here. GNINA score rows are
            # typically indented in the raw log, but `line.strip()` removes that
            # indentation and caused the parser to stop after the first pose.
            # Keeping the scan open is safe because `mode_line_pattern` is the
            # only path that records scores.
                
    except Exception as e:
        print(f"⚠️  Error parsing {log_file}: {e}")
        
    return scores, matched_pairlist

def _find_log_files(*search_dirs):
    """
    Search for GNINA log files across multiple candidate directories.

    Parameters
    ----------
    *search_dirs : Path
        Directories to search (in priority order).

    Returns
    -------
    List[Path]
        Sorted list of unique log-file paths.
    """
    found = {}
    for d in search_dirs:
        d = Path(d)
        if not d.exists():
            continue
        for f in d.glob("*.log"):
            if f.name.lower().endswith(".runner.log"):
                continue
            found.setdefault(f.name, f)
    if not found:
        for d in search_dirs:
            d = Path(d)
            if not d.exists():
                continue
            for f in d.glob("*_log"):
                found.setdefault(f.name, f)
    return sorted(found.values())


def generate_all_scores_csv(
    gnina_out_dir,
    output_file=None,
    pairlist_file=None,
    log_dir=None,
    log_files=None,
):
    """
    Generate all_scores.csv from GNINA log files.
    
    Parameters
    ----------
    gnina_out_dir : Path
        Directory containing GNINA output files (SDF poses).
    output_file : Path, optional
        Output CSV file path (default: gnina_out_dir/all_scores.csv).
    pairlist_file : Path, optional
        Path to pairlist.csv for accurate complex naming.
    log_dir : Path, optional
        Separate directory containing log files (GNINA HPC layout).
        When *None*, logs are searched inside *gnina_out_dir* and, as a
        fallback, in a sibling ``logs/`` directory.
    log_files : Iterable[Path], optional
        Explicit log-file list to parse. When provided, directory discovery is
        skipped and only the supplied files are parsed.
        
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    gnina_out_dir = Path(gnina_out_dir)
    
    if not gnina_out_dir.exists():
        print(f"❌ GNINA output directory not found: {gnina_out_dir}")
        return False
        
    if output_file is None:
        output_file = gnina_out_dir / "all_scores.csv"
        
    # Load pairlist mapping if provided
    pairlist_mapping = {}
    if pairlist_file and Path(pairlist_file).exists():
        print(f"🔍 Using pairlist mapping from: {pairlist_file}")
        pairlist_mapping = load_pairlist_mapping(pairlist_file)
        
    # Build list of directories to search for log files.
    # Priority: explicit log_dir > gnina_out_dir > sibling logs/
    search_dirs = []
    if log_dir:
        search_dirs.append(Path(log_dir))
    search_dirs.append(gnina_out_dir)
    sibling_logs = gnina_out_dir.parent / "logs"
    if sibling_logs.is_dir():
        search_dirs.append(sibling_logs)

    if log_files is not None:
        log_files = [Path(path) for path in log_files if Path(path).exists()]
    else:
        log_files = _find_log_files(*search_dirs)
        
    if not log_files:
        print(f"❌ No log files found in {[str(d) for d in search_dirs]}")
        return False
        
    print(f"📊 Found {len(log_files)} log files")
    
    # Parse all log files
    all_scores: List[Dict[str, object]] = []
    unmatched_logs: List[str] = []
    empty_logs: List[str] = []
    for log_file in log_files:
        print(f"🔍 Parsing {log_file.name}...")
        scores, matched_pairlist = parse_gnina_log(log_file, pairlist_mapping if pairlist_file else None)
        pose_file = gnina_out_dir / f"{log_file.stem}.sdf"
        if not pose_file.is_file() or not pose_file.stat().st_size:
            empty_logs.append(log_file.name + ":missing_pose_output")
            continue
        completion = Path(str(pose_file) + ".completion.json")
        if completion.is_file():
            try:
                from post_docking_analysis.pose_geometry import content_hash
                status = json.loads(completion.read_text(encoding="utf-8"))
                if status.get("status") != "completed" or status.get("pose_sha256") != content_hash(pose_file):
                    continue
            except (ValueError, OSError):
                continue
        if pairlist_mapping and not matched_pairlist:
            unmatched_logs.append(log_file.name)
            continue
        if not scores:
            empty_logs.append(log_file.name)
        all_scores.extend(scores)
        
    if not all_scores:
        print("❌ No scores extracted from log files")
        return False
        
    # Sort scores by tag and mode for consistency
    all_scores.sort(key=lambda x: (x['tag'], x['mode']))

    report_dir = Path(output_file).parent
    report_dir.mkdir(parents=True, exist_ok=True)
    if unmatched_logs:
        unmatched_file = report_dir / "unmatched_log_files.txt"
        unmatched_file.write_text(
            "\n".join(sorted(set(unmatched_logs))) + "\n",
            encoding="utf-8",
        )
        print(f"⚠️  Unmatched pairlist logs: {len(set(unmatched_logs))} (report: {unmatched_file})")
    if empty_logs:
        empty_file = report_dir / "empty_log_files.txt"
        empty_file.write_text(
            "\n".join(sorted(set(empty_logs))) + "\n",
            encoding="utf-8",
        )
        print(f"⚠️  Logs with no parsed scores: {len(set(empty_logs))} (report: {empty_file})")
        
    # Write to CSV
    try:
        with open(output_file, 'w', newline='') as csvfile:
            fieldnames = ['tag', 'mode', 'vina_affinity', 'cnn_affinity', 'cnn_score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for score in all_scores:
                writer.writerow(score)
                
        print(f"✅ Successfully generated {output_file}")
        print(f"   Total scores: {len(all_scores)}")
        print(f"   Unique complexes: {len(set(score['tag'] for score in all_scores))}")
        
        return True
        
    except Exception as e:
        print(f"❌ Error writing CSV file: {e}")
        return False

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Generate all_scores.csv from GNINA log files")
    parser.add_argument("gnina_out_dir", help="Directory containing GNINA output files")
    parser.add_argument("-o", "--output", help="Output CSV file (default: gnina_out_dir/all_scores.csv)")
    parser.add_argument("-p", "--pairlist", help="Pairlist CSV file for accurate complex naming")
    
    args = parser.parse_args()
    
    success = generate_all_scores_csv(args.gnina_out_dir, args.output, args.pairlist)
    
    if success:
        print("\n🎉 all_scores.csv generation completed successfully!")
    else:
        print("\n❌ all_scores.csv generation failed!")
        return 1
        
    return 0

if __name__ == "__main__":
    exit(main())
