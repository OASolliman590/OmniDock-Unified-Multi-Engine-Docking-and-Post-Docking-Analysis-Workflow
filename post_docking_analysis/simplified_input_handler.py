"""
Simplified Input Handler for Post-Docking Analysis.

Handles the simplified 3-folder input structure:
- sdf_folder: Docking poses (SDF files)
- log_folder: Docking logs (can be same as sdf_folder)
- receptors_folder: Receptor PDBQT files
"""
from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd


def find_sdf_files(sdf_folder: Path) -> List[Path]:
    """
    Find all SDF pose files in the folder.
    
    Parameters
    ----------
    sdf_folder : Path
        Folder containing SDF files
        
    Returns
    -------
    List[Path]
        List of SDF file paths
    """
    if not sdf_folder.exists():
        return []
    
    # Look for common SDF patterns
    sdf_files = set(sdf_folder.glob("*.sdf"))
    sdf_files.update(sdf_folder.glob("*_top.sdf"))
    
    return sorted([f for f in sdf_files if f.is_file()])


def find_log_files(log_folder: Path) -> List[Path]:
    """
    Find all log files in the folder.
    
    Parameters
    ----------
    log_folder : Path
        Folder containing log files
        
    Returns
    -------
    List[Path]
        List of log file paths
    """
    if not log_folder.exists():
        return []
    
    # Look for common log patterns
    log_files = list(log_folder.glob("*.log"))
    log_files.extend(log_folder.glob("*_log"))
    
    return sorted(log_files)


def find_receptor_files(receptors_folder: Path) -> List[Path]:
    """
    Find all receptor PDBQT files in the folder.
    
    Parameters
    ----------
    receptors_folder : Path
        Folder containing receptor PDBQT files
        
    Returns
    -------
    List[Path]
        List of receptor file paths
    """
    if not receptors_folder.exists():
        return []
    
    receptor_files = list(receptors_folder.glob("*.pdbqt"))
    
    return sorted(receptor_files)


def load_pairlist(pairlist_file: Optional[Path]) -> pd.DataFrame:
    """
    Load pairlist.csv if provided.
    
    Parameters
    ----------
    pairlist_file : Path or None
        Path to pairlist.csv
        
    Returns
    -------
    pd.DataFrame
        Pairlist DataFrame, empty if file doesn't exist
    """
    if pairlist_file is None or not pairlist_file.exists():
        return pd.DataFrame()
    
    try:
        pairlist_df = pd.read_csv(pairlist_file)
        pair_intent_file = pairlist_file.parent / "metadata" / "pair_intent.csv"
        if pair_intent_file.exists():
            try:
                pair_intent_df = pd.read_csv(pair_intent_file)
                required = {"receptor", "site_id", "ligand"}
                if required.issubset(pair_intent_df.columns):
                    return pair_intent_df
            except Exception:
                pass
        return pairlist_df
    except Exception as e:
        print(f"⚠️  Warning: Could not load pairlist.csv: {e}")
        return pd.DataFrame()


def auto_detect_pairlist_file(*paths: Optional[Path]) -> Optional[Path]:
    """
    Attempt deterministic pairlist auto-discovery from input-folder ancestry.

    Search order for each provided path:
    1) path directory
    2) up to four parent levels above that directory
    """
    roots: List[Path] = []
    seen = set()

    for raw_path in paths:
        if raw_path is None:
            continue
        candidate = Path(raw_path).expanduser()
        if candidate.is_file():
            candidate = candidate.parent
        try:
            candidate = candidate.resolve()
        except Exception:
            pass

        ancestry: List[Path] = []
        current = candidate
        for _ in range(5):
            ancestry.append(current)
            parent = current.parent
            if parent == current:
                break
            current = parent
        for root in ancestry:
            key = str(root)
            if key in seen:
                continue
            seen.add(key)
            if root.exists() and root.is_dir():
                roots.append(root)

    for root in roots:
        direct = root / "pairlist.csv"
        if direct.is_file():
            return direct

    for root in roots:
        for csv_path in sorted(root.glob("*pairlist*.csv")):
            if csv_path.is_file():
                return csv_path

    return None


def match_poses_to_receptors(
    sdf_files: List[Path],
    receptor_files: List[Path],
    pairlist_df: pd.DataFrame = None
) -> List[Dict]:
    """
    Match SDF pose files to receptors using pairlist or filename patterns.
    
    Parameters
    ----------
    sdf_files : List[Path]
        List of SDF pose files
    receptor_files : List[Path]
        List of receptor PDBQT files
    pairlist_df : pd.DataFrame, optional
        Pairlist DataFrame for matching
        
    Returns
    -------
    List[Dict]
        List of matched complexes with receptor and pose info
    """
    complexes = []

    def _tokens(name: str) -> List[str]:
        if not name:
            return []
        try:
            stem = Path(str(name)).stem
        except Exception:
            stem = ""
        tokens = [str(name), stem]
        return [t.lower() for t in tokens if t]

    def _match_receptor(receptor_name: str) -> Optional[Path]:
        receptor_tokens = _tokens(receptor_name)
        for rf in receptor_files:
            rf_name = rf.name.lower()
            if any(tok in rf_name for tok in receptor_tokens):
                return rf
        return None

    def _match_sdf(receptor_name: str, ligand_name: str) -> Optional[Path]:
        receptor_tokens = _tokens(receptor_name)
        ligand_tokens = _tokens(ligand_name)

        # Prefer SDFs that contain both receptor and ligand identifiers.
        candidates = []
        for sf in sdf_files:
            sf_name = sf.name.lower()
            if any(tok in sf_name for tok in receptor_tokens) and any(tok in sf_name for tok in ligand_tokens):
                candidates.append(sf)

        if candidates:
            return sorted(candidates)[0]

        # Fallback: ligand-only match (legacy behavior).
        for sf in sdf_files:
            sf_name = sf.name.lower()
            if any(tok in sf_name for tok in ligand_tokens):
                return sf
        return None
    
    # If pairlist is available, use it for matching
    if pairlist_df is not None and not pairlist_df.empty:
        for _, row in pairlist_df.iterrows():
            receptor_name = row.get('receptor', '')
            ligand_name = row.get('ligand', '')
            site_id = row.get('site_id', 'unknown')
            
            # Find matching receptor file
            receptor_file = _match_receptor(str(receptor_name))
            
            # Find matching SDF file
            sdf_file = _match_sdf(str(receptor_name), str(ligand_name))
            
            if receptor_file and sdf_file:
                complexes.append({
                    'complex_name': f"{receptor_name}_{site_id}_{ligand_name}",
                    'receptor_file': receptor_file,
                    'pose_file': sdf_file,
                    'site_id': site_id,
                    'receptor_name': receptor_name,
                    'ligand_name': ligand_name
                })
    else:
        # Fallback: filename pattern matching
        # Extract base names and try to match
        for sdf_file in sdf_files:
            # Strip common suffixes added by GNINA HPC workflow
            sdf_base = sdf_file.stem
            for suffix in ('_poses', '_top', '_out'):
                if sdf_base.endswith(suffix):
                    sdf_base = sdf_base[:-len(suffix)]
            
            # Try to find matching receptor
            for receptor_file in receptor_files:
                receptor_base = receptor_file.stem.replace('_cleaned', '')
                
                # Simple matching: if receptor name contains part of SDF name or vice versa
                if sdf_base in receptor_base or receptor_base in sdf_base:
                    complexes.append({
                        'complex_name': sdf_base,
                        'receptor_file': receptor_file,
                        'pose_file': sdf_file,
                        'site_id': 'unknown',
                        'receptor_name': receptor_base,
                        'ligand_name': sdf_base
                    })
                    break
    
    return complexes
