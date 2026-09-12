"""
Enhanced RMSD Analysis with actual structure-based calculations.

Calculates RMSD from PDB structures, performs pose clustering,
and analyzes conformational diversity.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor
import os
from sklearn.cluster import KMeans, DBSCAN
import matplotlib.pyplot as plt
import seaborn as sns
import logging

try:
    from Bio.PDB import PDBParser
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

logger = logging.getLogger(__name__)


def _rmsd_pair_worker(payload: Tuple[int, int, str, str, bool]) -> Tuple[int, int, float]:
    i, j, file_a, file_b, ligand_only = payload
    rmsd_value = calculate_rmsd_between_structures(Path(file_a), Path(file_b), ligand_only=ligand_only)
    return i, j, rmsd_value


def _should_include_atom(residue, atom, ligand_only: bool) -> bool:
    """Return True if atom should be included in RMSD calculation."""
    # Skip hydrogens for stable heavy-atom RMSD.
    if atom.name.startswith("H"):
        return False

    is_standard_residue = residue.id[0] == " "
    if ligand_only:
        # Ligand atoms are written as HETATM/non-standard residues.
        return not is_standard_residue
    return True


def _kabsch_aligned_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """Legacy helper now preserves the supplied coordinate frame."""
    if coords1.shape != coords2.shape or coords1.size == 0:
        return np.nan
    return float(np.sqrt(np.mean(np.sum((coords1 - coords2) ** 2, axis=1))))



def calculate_rmsd_between_structures(pdb_file1: Path, pdb_file2: Path, ligand_only: bool = True) -> float:
    """Mapped same-ligand RMSD in a declared common receptor frame."""
    return _calculate_rmsd_simple(pdb_file1, pdb_file2, ligand_only)



def _calculate_rmsd_simple(pdb_file1: Path, pdb_file2: Path, ligand_only: bool = True) -> float:
    from .pose_geometry import pose_rmsd
    import json
    if not ligand_only:
        return np.nan
    try:
        files = [Path(pdb_file1), Path(pdb_file2)]
        frames = [json.loads(p.with_suffix(".geometry.json").read_text())["receptor_frame_id"] for p in files]
        if not all(frames) or frames[0] != frames[1]:
            return np.nan
        return pose_rmsd(files[0].with_suffix(".ligand.sdf"), files[1].with_suffix(".ligand.sdf"))
    except Exception as exc:
        logger.warning("RMSD not evaluable: %s", exc)
        return np.nan



def calculate_rmsd_matrix_from_pdbs(
    pdb_files: List[Path],
    ligand_only: bool = True,
    max_pairs: Optional[int] = None,
    num_workers: int = 0,
) -> Tuple[np.ndarray, List[str]]:
    """
    Calculate RMSD matrix from PDB files.
    
    Parameters
    ----------
    pdb_files : List[Path]
        List of PDB file paths
    ligand_only : bool
        Calculate RMSD only for ligand atoms
    max_pairs : int, optional
        Maximum number of pairs to calculate (for performance)
        
    Returns
    -------
    Tuple[np.ndarray, List[str]]
        RMSD matrix and list of filenames
    """
    logger.info(f"📏 Calculating RMSD matrix from {len(pdb_files)} PDB files...")
    
    n = len(pdb_files)
    rmsd_matrix = np.zeros((n, n))
    filenames = [f.stem for f in pdb_files]
    
    total_pairs = n * (n - 1) // 2
    if int(num_workers or 0) <= 0:
        worker_count = max((os.cpu_count() or 2) - 1, 1)
    else:
        worker_count = max(int(num_workers), 1)

    pair_jobs: List[Tuple[int, int, str, str, bool]] = []
    skipped_for_limit = 0
    for i in range(n):
        for j in range(i + 1, n):
            if max_pairs and len(pair_jobs) >= max_pairs:
                rmsd_matrix[i, j] = np.nan
                rmsd_matrix[j, i] = np.nan
                skipped_for_limit += 1
                continue
            pair_jobs.append((i, j, str(pdb_files[i]), str(pdb_files[j]), ligand_only))

    calculated = 0
    if worker_count > 1 and len(pair_jobs) > 1:
        logger.info(f"⚙️ RMSD parallel mode enabled with {worker_count} workers")
        chunksize = max(1, len(pair_jobs) // max(worker_count * 4, 1))
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            for i, j, rmsd in executor.map(_rmsd_pair_worker, pair_jobs, chunksize=chunksize):
                rmsd_matrix[i, j] = rmsd
                rmsd_matrix[j, i] = rmsd
                calculated += 1
                if calculated % 25 == 0:
                    logger.debug(f"  Calculated {calculated}/{len(pair_jobs)} pairs...")
    else:
        for i, j, file_a, file_b, lig_only in pair_jobs:
            rmsd = calculate_rmsd_between_structures(Path(file_a), Path(file_b), ligand_only=lig_only)
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd
            calculated += 1
            if calculated % 10 == 0:
                logger.debug(f"  Calculated {calculated}/{len(pair_jobs)} pairs...")
    
    logger.info(f"✅ RMSD matrix calculated ({calculated} pairs)")
    if skipped_for_limit:
        logger.info(f"ℹ️ Skipped {skipped_for_limit} pairs due to max_pairs limit")
    return rmsd_matrix, filenames


def analyze_pose_clustering_enhanced(
    poses_data: pd.DataFrame,
    rmsd_matrix: np.ndarray,
    pdb_files: List[Path],
    method: str = 'kmeans',
    n_clusters: int = 3,
    eps: float = 2.0,
    min_samples: int = 2
) -> Dict:
    """
    Enhanced pose clustering analysis with actual RMSD data.
    
    Parameters
    ----------
    poses_data : pd.DataFrame
        DataFrame containing pose metadata
    rmsd_matrix : np.ndarray
        RMSD matrix between poses
    pdb_files : List[Path]
        List of PDB file paths
    method : str
        Clustering method ('kmeans' or 'dbscan')
    n_clusters : int
        Number of clusters for K-means
    eps : float
        Epsilon parameter for DBSCAN
    min_samples : int
        Minimum samples for DBSCAN
        
    Returns
    -------
    Dict
        Dictionary containing clustering results
    """
    logger.info(f"🔍 Analyzing pose clustering using {method}...")
    
    # Remove NaN values for clustering
    valid_mask = ~np.isnan(rmsd_matrix).any(axis=1)
    valid_indices = np.where(valid_mask)[0]
    
    if len(valid_indices) < n_clusters:
        logger.warning(f"⚠️  Not enough valid poses for {n_clusters} clusters")
        n_clusters = max(1, len(valid_indices) // 2)  # At least 1, or half of available
    
    if n_clusters == 0 or len(valid_indices) == 0:
        logger.warning("⚠️  No valid poses for clustering")
        return {
            'poses_with_clusters': poses_data.copy(),
            'cluster_summary': pd.DataFrame(),
            'cluster_centroids': pd.DataFrame(),
            'rmsd_matrix': rmsd_matrix,
            'cluster_labels': np.array([]),
            'valid_indices': np.array([])
        }
    
    valid_rmsd = rmsd_matrix[np.ix_(valid_indices, valid_indices)]
    valid_poses = poses_data.iloc[valid_indices].copy()
    valid_poses = valid_poses.reset_index(drop=True)
    
    # Perform clustering
    if method == 'kmeans':
        clusterer = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        cluster_labels = clusterer.fit_predict(valid_rmsd)
    elif method == 'dbscan':
        clusterer = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
        cluster_labels = clusterer.fit_predict(valid_rmsd)
    else:
        raise ValueError(f"Unknown clustering method: {method}")
    
    # Add cluster labels to poses data
    valid_poses['cluster'] = cluster_labels
    
    # Analyze clusters
    cluster_summary = []
    cluster_centroids = []
    
    for cluster_id in sorted(set(cluster_labels)):
        if cluster_id == -1:  # Skip noise points in DBSCAN
            continue
        
        cluster_poses = valid_poses[valid_poses['cluster'] == cluster_id]
        
        if len(cluster_poses) == 0:
            continue
        
        # Cluster statistics
        cluster_summary.append({
            'cluster': cluster_id,
            'size': len(cluster_poses),
            'avg_affinity': cluster_poses['vina_affinity'].mean() if 'vina_affinity' in cluster_poses.columns else np.nan,
            'min_affinity': cluster_poses['vina_affinity'].min() if 'vina_affinity' in cluster_poses.columns else np.nan,
            'max_affinity': cluster_poses['vina_affinity'].max() if 'vina_affinity' in cluster_poses.columns else np.nan,
            'avg_rmsd': valid_rmsd[np.ix_(cluster_poses.index, cluster_poses.index)].mean() if len(cluster_poses) > 1 else 0.0
        })
        
        # Find centroid (pose with lowest average RMSD to others in cluster)
        if len(cluster_poses) > 1:
            cluster_indices = cluster_poses.index
            cluster_rmsd = valid_rmsd[np.ix_(cluster_indices, cluster_indices)]
            avg_rmsd_per_pose = cluster_rmsd.mean(axis=1)
            centroid_idx = cluster_indices[np.argmin(avg_rmsd_per_pose)]
        else:
            centroid_idx = cluster_poses.index[0]
        
        cluster_centroids.append({
            'cluster': cluster_id,
            'centroid_pose': valid_poses.loc[centroid_idx, 'tag'] if 'tag' in valid_poses.columns else f"pose_{centroid_idx}",
            'centroid_affinity': valid_poses.loc[centroid_idx, 'vina_affinity'] if 'vina_affinity' in valid_poses.columns else np.nan,
            'cluster_size': len(cluster_poses),
            'avg_affinity': cluster_poses['vina_affinity'].mean() if 'vina_affinity' in cluster_poses.columns else np.nan
        })
    
    cluster_summary_df = pd.DataFrame(cluster_summary)
    cluster_centroids_df = pd.DataFrame(cluster_centroids)
    
    logger.info(f"✅ Pose clustering completed")
    logger.info(f"   Found {len(cluster_summary_df)} clusters")
    if len(cluster_centroids_df) > 0:
        logger.info(f"   Best cluster affinity: {cluster_centroids_df['avg_affinity'].min():.2f} kcal/mol")
    
    return {
        'poses_with_clusters': valid_poses,
        'cluster_summary': cluster_summary_df,
        'cluster_centroids': cluster_centroids_df,
        'rmsd_matrix': valid_rmsd,
        'cluster_labels': cluster_labels,
        'valid_indices': valid_indices
    }


def analyze_conformational_diversity_enhanced(
    poses_data: pd.DataFrame,
    rmsd_matrix: np.ndarray
) -> Dict:
    """
    Enhanced conformational diversity analysis.
    
    Parameters
    ----------
    poses_data : pd.DataFrame
        DataFrame containing pose metadata
    rmsd_matrix : np.ndarray
        RMSD matrix between poses
        
    Returns
    -------
    Dict
        Dictionary containing diversity analysis results
    """
    logger.info("🌊 Analyzing conformational diversity...")
    
    # Remove NaN values
    valid_mask = ~np.isnan(rmsd_matrix).any(axis=1)
    valid_indices = np.where(valid_mask)[0]
    valid_rmsd = rmsd_matrix[np.ix_(valid_indices, valid_indices)]
    valid_poses = poses_data.iloc[valid_indices].copy()
    valid_poses = valid_poses.reset_index(drop=True)

    if len(valid_poses) < 2:
        logger.warning("⚠️  Not enough valid poses for diversity statistics")
        empty_df = pd.DataFrame()
        overall_stats = {
            'total_poses': len(valid_poses),
            'avg_pairwise_rmsd': np.nan,
            'max_pairwise_rmsd': np.nan,
            'min_pairwise_rmsd': np.nan,
            'rmsd_std': np.nan,
            'median_pairwise_rmsd': np.nan
        }
        return {
            'diversity_metrics': empty_df,
            'overall_stats': overall_stats,
            'most_diverse_poses': empty_df,
            'least_diverse_poses': empty_df,
            'rmsd_matrix': valid_rmsd
        }
    
    # Calculate diversity metrics for each pose
    diversity_metrics = []
    
    for i, (idx, pose) in enumerate(valid_poses.iterrows()):
        # Calculate average RMSD to all other poses
        other_rmsds = np.concatenate([valid_rmsd[i, :i], valid_rmsd[i, i+1:]])
        other_rmsds = other_rmsds[~np.isnan(other_rmsds)]
        
        if len(other_rmsds) == 0:
            continue
        
        diversity_metrics.append({
            'tag': pose.get('tag', f'pose_{idx}'),
            'vina_affinity': pose.get('vina_affinity', np.nan),
            'avg_rmsd_to_others': np.mean(other_rmsds),
            'max_rmsd_to_others': np.max(other_rmsds),
            'min_rmsd_to_others': np.min(other_rmsds),
            'rmsd_std': np.std(other_rmsds),
            'median_rmsd': np.median(other_rmsds)
        })
    
    diversity_df = pd.DataFrame(diversity_metrics)
    
    # Overall diversity statistics
    upper_triangle = valid_rmsd[np.triu_indices_from(valid_rmsd, k=1)]
    upper_triangle = upper_triangle[~np.isnan(upper_triangle)]
    
    if len(upper_triangle) == 0:
        overall_stats = {
            'total_poses': len(valid_poses),
            'avg_pairwise_rmsd': np.nan,
            'max_pairwise_rmsd': np.nan,
            'min_pairwise_rmsd': np.nan,
            'rmsd_std': np.nan,
            'median_pairwise_rmsd': np.nan
        }
    else:
        overall_stats = {
            'total_poses': len(valid_poses),
            'avg_pairwise_rmsd': np.mean(upper_triangle),
            'max_pairwise_rmsd': np.max(upper_triangle),
            'min_pairwise_rmsd': np.min(upper_triangle),
            'rmsd_std': np.std(upper_triangle),
            'median_pairwise_rmsd': np.median(upper_triangle)
        }
    
    # Identify most diverse poses
    if len(diversity_df) > 0:
        most_diverse = diversity_df.nlargest(5, 'avg_rmsd_to_others')
        least_diverse = diversity_df.nsmallest(5, 'avg_rmsd_to_others')
    else:
        most_diverse = pd.DataFrame()
        least_diverse = pd.DataFrame()
    
    logger.info(f"✅ Conformational diversity analysis completed")
    if not np.isnan(overall_stats['avg_pairwise_rmsd']):
        logger.info(f"   Average pairwise RMSD: {overall_stats['avg_pairwise_rmsd']:.2f} Å")
    else:
        logger.info("   Average pairwise RMSD: N/A")
    if len(most_diverse) > 0:
        logger.info(f"   Most diverse pose: {most_diverse.iloc[0]['tag']} ({most_diverse.iloc[0]['avg_rmsd_to_others']:.2f} Å)")
    
    return {
        'diversity_metrics': diversity_df,
        'overall_stats': overall_stats,
        'most_diverse_poses': most_diverse,
        'least_diverse_poses': least_diverse,
        'rmsd_matrix': valid_rmsd
    }


def create_rmsd_visualizations_enhanced(
    clustering_results: Dict,
    diversity_results: Dict,
    output_dir: Path,
    dpi: int = 300
) -> List[Path]:
    """
    Create enhanced RMSD visualizations with heatmaps and 2D plots.
    
    Parameters
    ----------
    clustering_results : Dict
        Results from pose clustering analysis
    diversity_results : Dict
        Results from diversity analysis
    output_dir : Path
        Output directory for visualizations
    dpi : int
        DPI for image output
        
    Returns
    -------
    List[Path]
        List of created visualization files
    """
    logger.info("📊 Creating enhanced RMSD visualizations...")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    created_files = []
    
    rmsd_matrix = np.asarray(clustering_results.get('rmsd_matrix', np.array([])), dtype=float)
    poses_data = clustering_results.get('poses_with_clusters', pd.DataFrame())

    def _pose_label(row: pd.Series, index: int) -> str:
        if isinstance(row, pd.Series):
            if 'tag' in row and pd.notna(row['tag']):
                return str(row['tag'])
            protein = str(row.get('protein_display_name') or row.get('protein') or "").strip()
            ligand = str(row.get('ligand') or "").strip()
            pose = row.get('pose', index + 1)
            if protein or ligand:
                return f"{protein}:{ligand}:p{pose}"
        return f"pose_{index + 1}"

    # 1. RMSD Heatmap (correctly labeled with robust NaN handling)
    if rmsd_matrix.size > 0:
        clean_matrix = np.array(rmsd_matrix, copy=True)
        finite_values = clean_matrix[np.isfinite(clean_matrix)]
        fill_value = float(np.nanmax(finite_values)) if finite_values.size else 0.0
        clean_matrix[~np.isfinite(clean_matrix)] = fill_value

        n_poses = clean_matrix.shape[0]
        # Keep labels readable by downsampling ticks for large matrices.
        max_ticks = 25
        if n_poses <= max_ticks:
            tick_positions = np.arange(n_poses)
        else:
            tick_positions = np.linspace(0, n_poses - 1, max_ticks, dtype=int)
        tick_positions = np.unique(tick_positions)

        labels = []
        if isinstance(poses_data, pd.DataFrame) and not poses_data.empty:
            for idx in range(min(n_poses, len(poses_data))):
                labels.append(_pose_label(poses_data.iloc[idx], idx))
        while len(labels) < n_poses:
            labels.append(f"pose_{len(labels) + 1}")

        short_labels = [labels[i][:24] for i in tick_positions]
        fig, ax = plt.subplots(figsize=(14, 12))
        sns.heatmap(
            clean_matrix,
            cmap='viridis',
            square=True,
            cbar_kws={'label': 'RMSD (Å)'},
            xticklabels=False,
            yticklabels=False,
            ax=ax
        )
        ax.set_xticks(tick_positions + 0.5)
        ax.set_yticks(tick_positions + 0.5)
        ax.set_xticklabels(short_labels, rotation=90, fontsize=8)
        ax.set_yticklabels(short_labels, rotation=0, fontsize=8)
        ax.set_title('RMSD Similarity Matrix', fontsize=16, fontweight='bold')
        ax.set_xlabel('Pose', fontsize=12)
        ax.set_ylabel('Pose', fontsize=12)
        plt.tight_layout()
    else:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, "RMSD heatmap unavailable\n(no valid RMSD matrix)", ha='center', va='center')
        ax.axis('off')
        ax.set_title('RMSD Similarity Matrix', fontsize=16, fontweight='bold')
        plt.tight_layout()

    heatmap_file = output_dir / 'rmsd_heatmap.png'
    plt.savefig(heatmap_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    created_files.append(heatmap_file)
    
    # 2. Cluster analysis plot (always render; placeholder if data missing)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    if isinstance(poses_data, pd.DataFrame) and len(poses_data) > 0 and 'cluster' in poses_data.columns:
        unique_clusters = sorted([c for c in poses_data['cluster'].dropna().unique() if c != -1])
        colors = plt.cm.Set3(np.linspace(0, 1, max(len(unique_clusters), 1)))
        if unique_clusters and 'vina_affinity' in poses_data.columns:
            for i, cluster in enumerate(unique_clusters):
                cluster_data = poses_data[poses_data['cluster'] == cluster]
                if len(cluster_data) > 1:
                    jitter = np.linspace(-0.08, 0.08, len(cluster_data))
                else:
                    jitter = np.zeros(len(cluster_data))
                ax1.scatter(
                    np.full(len(cluster_data), cluster) + jitter,
                    cluster_data['vina_affinity'],
                    c=[colors[i]],
                    label=f'Cluster {cluster}',
                    alpha=0.75,
                    s=50
                )
            ax1.legend()
        else:
            ax1.text(0.5, 0.5, "No clustered affinity data", ha='center', va='center')
        ax1.set_xlabel('Cluster', fontsize=12)
        ax1.set_ylabel('Binding Affinity (kcal/mol)', fontsize=12)
        ax1.set_title('Binding Affinity by Cluster', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        if unique_clusters:
            ax1.set_xticks(unique_clusters)

        cluster_sizes = poses_data['cluster'].value_counts().sort_index()
        cluster_sizes = cluster_sizes[cluster_sizes.index != -1]
        if len(cluster_sizes) > 0:
            x_positions = np.arange(len(cluster_sizes))
            ax2.bar(x_positions, cluster_sizes.values, color=colors[:len(cluster_sizes)])
            ax2.set_xticks(x_positions)
            ax2.set_xticklabels(cluster_sizes.index.astype(str))
            ax2.set_ylabel('Number of Poses', fontsize=12)
        else:
            ax2.text(0.5, 0.5, "No cluster size data", ha='center', va='center')
        ax2.set_xlabel('Cluster', fontsize=12)
        ax2.set_title('Cluster Size Distribution', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
    else:
        ax1.text(0.5, 0.5, "Cluster analysis unavailable", ha='center', va='center')
        ax2.text(0.5, 0.5, "Cluster distribution unavailable", ha='center', va='center')
        ax1.axis('off')
        ax2.axis('off')
    plt.tight_layout()
    cluster_file = output_dir / 'cluster_analysis.png'
    plt.savefig(cluster_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    created_files.append(cluster_file)
    
    # 3. Diversity analysis plot (always render; placeholder if missing)
    diversity_df = diversity_results.get('diversity_metrics', pd.DataFrame())
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    if isinstance(diversity_df, pd.DataFrame) and len(diversity_df) > 0:
        if 'vina_affinity' in diversity_df.columns and 'avg_rmsd_to_others' in diversity_df.columns:
            scatter = ax1.scatter(
                diversity_df['avg_rmsd_to_others'],
                diversity_df['vina_affinity'],
                alpha=0.75,
                c=diversity_df['vina_affinity'],
                cmap='viridis',
                s=50
            )
            ax1.set_xlabel('Average RMSD to Other Poses (Å)', fontsize=12)
            ax1.set_ylabel('Binding Affinity (kcal/mol)', fontsize=12)
            ax1.set_title('Binding Affinity vs Conformational Diversity', fontsize=14, fontweight='bold')
            ax1.grid(True, alpha=0.3)
            plt.colorbar(scatter, ax=ax1, label='Affinity (kcal/mol)')
        else:
            ax1.text(0.5, 0.5, "Diversity/affinity columns missing", ha='center', va='center')
            ax1.axis('off')

        if 'avg_rmsd_to_others' in diversity_df.columns:
            ax2.hist(diversity_df['avg_rmsd_to_others'], bins=20, alpha=0.75, color='skyblue', edgecolor='black')
            mean_rmsd = diversity_df['avg_rmsd_to_others'].mean()
            ax2.axvline(mean_rmsd, color='red', linestyle='--', label=f'Mean: {mean_rmsd:.2f} Å')
            ax2.set_xlabel('Average RMSD to Other Poses (Å)', fontsize=12)
            ax2.set_ylabel('Frequency', fontsize=12)
            ax2.set_title('Distribution of Conformational Diversity', fontsize=14, fontweight='bold')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
        else:
            ax2.text(0.5, 0.5, "No diversity distribution data", ha='center', va='center')
            ax2.axis('off')
    else:
        ax1.text(0.5, 0.5, "Diversity analysis unavailable", ha='center', va='center')
        ax2.text(0.5, 0.5, "Diversity distribution unavailable", ha='center', va='center')
        ax1.axis('off')
        ax2.axis('off')
    plt.tight_layout()
    diversity_file = output_dir / 'diversity_analysis.png'
    plt.savefig(diversity_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    created_files.append(diversity_file)
    
    logger.info(f"✅ Created {len(created_files)} RMSD visualizations")
    return created_files
