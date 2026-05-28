"""
Correlation analysis module for post-docking analysis pipeline.

This module handles correlation analysis between different scoring functions,
particularly Vina vs CNN scores, and other comparative analyses.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')


MIN_CORRELATION_N = 5
LOW_POWER_N = 10


def _bh_fdr_adjust(p_values: List[float]) -> List[float]:
    """
    Benjamini-Hochberg FDR correction.
    """
    if not p_values:
        return []
    indexed = [(idx, float(p)) for idx, p in enumerate(p_values)]
    indexed.sort(key=lambda item: item[1])
    m = len(indexed)
    adjusted = [np.nan] * m
    prev = 1.0
    for rank, (orig_idx, pval) in enumerate(reversed(indexed), start=1):
        bh = pval * m / float(m - rank + 1)
        prev = min(prev, bh)
        adjusted[orig_idx] = min(max(prev, 0.0), 1.0)
    return adjusted


def _rank_normalize(series: pd.Series, *, lower_is_better: bool = True) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    finite = values[np.isfinite(values)]
    if finite.empty:
        return pd.Series(np.nan, index=series.index, dtype=float)
    if len(finite) == 1:
        out = pd.Series(np.nan, index=series.index, dtype=float)
        out.loc[finite.index] = 1.0
        return out
    # Lower affinity is better; convert to [0,1] rank score where 1.0=best.
    ascending = not lower_is_better
    ranks = finite.rank(method="average", ascending=ascending)
    scaled = (ranks - 1.0) / float(len(finite) - 1)
    out = pd.Series(np.nan, index=series.index, dtype=float)
    out.loc[finite.index] = scaled
    return out


def _safe_corr_with_stats(series_a: pd.Series, series_b: pd.Series) -> Dict[str, object]:
    paired = pd.DataFrame({"a": pd.to_numeric(series_a, errors="coerce"), "b": pd.to_numeric(series_b, errors="coerce")}).dropna()
    n = int(len(paired))
    result = {
        "paired_rows": n,
        "insufficient_n": n < MIN_CORRELATION_N,
        "low_power": (MIN_CORRELATION_N <= n < LOW_POWER_N),
        "pearson_corr": np.nan,
        "pearson_p_value": np.nan,
        "spearman_corr": np.nan,
        "spearman_p_value": np.nan,
    }
    if n < MIN_CORRELATION_N:
        return result
    rank_a = _rank_normalize(paired["a"], lower_is_better=True)
    rank_b = _rank_normalize(paired["b"], lower_is_better=True)
    try:
        pearson_corr, pearson_p = pearsonr(rank_a, rank_b)
        result["pearson_corr"] = float(pearson_corr)
        result["pearson_p_value"] = float(pearson_p)
    except Exception:
        pass
    try:
        spearman_corr, spearman_p = spearmanr(rank_a, rank_b)
        result["spearman_corr"] = float(spearman_corr)
        result["spearman_p_value"] = float(spearman_p)
    except Exception:
        pass
    return result

def analyze_vina_cnn_correlation(scores_df: pd.DataFrame) -> Dict:
    """
    Analyze correlation between Vina and CNN scores.
    
    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame containing Vina and CNN scores
        
    Returns
    -------
    Dict
        Dictionary containing correlation analysis results
    """
    print("📊 Analyzing Vina vs CNN correlation...")
    
    # Check if CNN scores are available
    if 'cnn_affinity' not in scores_df.columns or 'cnn_score' not in scores_df.columns:
        print("⚠️ CNN scores not available - skipping correlation analysis")
        return {'error': 'CNN scores not available'}
    
    # Remove rows with missing values
    valid_data = scores_df.dropna(subset=['vina_affinity', 'cnn_affinity', 'cnn_score'])
    
    if len(valid_data) == 0:
        print("⚠️ No valid data for correlation analysis")
        return {'error': 'No valid data'}
    
    # Calculate correlations
    vina_affinity = valid_data['vina_affinity']
    cnn_affinity = valid_data['cnn_affinity']
    cnn_score = valid_data['cnn_score']
    
    corr_pairs = [
        ("vina_cnn_affinity", vina_affinity, cnn_affinity),
        ("vina_cnn_score", vina_affinity, cnn_score),
        ("cnn_affinity_score", cnn_affinity, cnn_score),
    ]
    pair_stats = {name: _safe_corr_with_stats(left, right) for name, left, right in corr_pairs}

    # Collect p-values for BH correction (Pearson + Spearman tests).
    raw_p_values: List[float] = []
    p_slots: List[Tuple[str, str]] = []
    for pair_name, stats_payload in pair_stats.items():
        for metric_name, p_col in (("pearson", "pearson_p_value"), ("spearman", "spearman_p_value")):
            p_value = stats_payload.get(p_col)
            if pd.notna(p_value):
                raw_p_values.append(float(p_value))
                p_slots.append((pair_name, metric_name))
    corrected = _bh_fdr_adjust(raw_p_values)
    q_lookup: Dict[Tuple[str, str], float] = {
        slot: corrected[idx] for idx, slot in enumerate(p_slots)
    }
    
    # Create correlation matrix
    correlation_data = valid_data[['vina_affinity', 'cnn_affinity', 'cnn_score']].corr()
    
    # Statistical significance
    n_samples = len(valid_data)
    
    correlation_results = {
        'n_samples': n_samples,
        'pearson_correlations': {},
        'spearman_correlations': {},
        'correlation_matrix': correlation_data,
        'valid_data': valid_data,
        'minimum_required_n': MIN_CORRELATION_N,
        'low_power_threshold_n': LOW_POWER_N,
    }

    for pair_name, stats_payload in pair_stats.items():
        pearson_p = stats_payload.get("pearson_p_value")
        spearman_p = stats_payload.get("spearman_p_value")
        pearson_q = q_lookup.get((pair_name, "pearson"), np.nan)
        spearman_q = q_lookup.get((pair_name, "spearman"), np.nan)
        correlation_results['pearson_correlations'][pair_name] = {
            'correlation': stats_payload.get("pearson_corr", np.nan),
            'p_value': pearson_p,
            'q_value': pearson_q,
            'significant': bool(pd.notna(pearson_q) and pearson_q < 0.05),
            'paired_rows': int(stats_payload.get("paired_rows", 0)),
            'insufficient_n': bool(stats_payload.get("insufficient_n", False)),
            'low_power': bool(stats_payload.get("low_power", False)),
        }
        correlation_results['spearman_correlations'][pair_name] = {
            'correlation': stats_payload.get("spearman_corr", np.nan),
            'p_value': spearman_p,
            'q_value': spearman_q,
            'significant': bool(pd.notna(spearman_q) and spearman_q < 0.05),
            'paired_rows': int(stats_payload.get("paired_rows", 0)),
            'insufficient_n': bool(stats_payload.get("insufficient_n", False)),
            'low_power': bool(stats_payload.get("low_power", False)),
        }
    
    print(f"✅ Correlation analysis completed for {n_samples} samples")
    vina_affinity_pearson = correlation_results['pearson_correlations'].get('vina_cnn_affinity', {})
    vina_score_pearson = correlation_results['pearson_correlations'].get('vina_cnn_score', {})
    print(
        "   Vina-CNN Affinity Pearson r: "
        f"{float(vina_affinity_pearson.get('correlation', np.nan)):.3f} "
        f"(p={float(vina_affinity_pearson.get('p_value', np.nan)):.3f}, "
        f"q={float(vina_affinity_pearson.get('q_value', np.nan)):.3f})"
    )
    print(
        "   Vina-CNN Score Pearson r: "
        f"{float(vina_score_pearson.get('correlation', np.nan)):.3f} "
        f"(p={float(vina_score_pearson.get('p_value', np.nan)):.3f}, "
        f"q={float(vina_score_pearson.get('q_value', np.nan)):.3f})"
    )
    
    return correlation_results

def analyze_score_distributions(scores_df: pd.DataFrame) -> Dict:
    """
    Analyze distributions of different scoring functions.
    
    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame containing scoring data
        
    Returns
    -------
    Dict
        Dictionary containing distribution analysis results
    """
    print("📈 Analyzing score distributions...")
    
    # Calculate descriptive statistics
    score_columns = ['vina_affinity']
    if 'cnn_affinity' in scores_df.columns:
        score_columns.append('cnn_affinity')
    if 'cnn_score' in scores_df.columns:
        score_columns.append('cnn_score')
    
    distribution_stats = {}
    for col in score_columns:
        if col in scores_df.columns:
            data = scores_df[col].dropna()
            distribution_stats[col] = {
                'mean': data.mean(),
                'std': data.std(),
                'min': data.min(),
                'max': data.max(),
                'median': data.median(),
                'q25': data.quantile(0.25),
                'q75': data.quantile(0.75),
                'skewness': stats.skew(data),
                'kurtosis': stats.kurtosis(data),
                'n_samples': len(data)
            }
    
    # Normality tests
    normality_tests = {}
    for col in score_columns:
        if col in scores_df.columns:
            data = scores_df[col].dropna()
            if len(data) >= 3:  # Minimum sample size for Shapiro-Wilk
                try:
                    shapiro_stat, shapiro_p = stats.shapiro(data)
                    normality_tests[col] = {
                        'shapiro_wilk': {
                            'statistic': shapiro_stat,
                            'p_value': shapiro_p,
                            'normal': shapiro_p > 0.05
                        }
                    }
                except:
                    normality_tests[col] = {'error': 'Could not perform normality test'}
    
    distribution_results = {
        'descriptive_stats': distribution_stats,
        'normality_tests': normality_tests
    }
    
    print("✅ Score distribution analysis completed")
    return distribution_results

def analyze_score_agreement(scores_df: pd.DataFrame, threshold: float = 0.1) -> Dict:
    """
    Analyze agreement between different scoring functions.
    
    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame containing scoring data
    threshold : float
        Threshold for considering scores in agreement
        
    Returns
    -------
    Dict
        Dictionary containing agreement analysis results
    """
    print("🤝 Analyzing score agreement...")
    
    if 'cnn_affinity' not in scores_df.columns:
        print("⚠️ CNN scores not available - skipping agreement analysis")
        return {'error': 'CNN scores not available'}
    
    # Rank-normalize scores to [0,1] with 1.0=stronger binder across engines.
    vina_norm = _rank_normalize(scores_df['vina_affinity'], lower_is_better=True)
    cnn_norm = _rank_normalize(scores_df['cnn_affinity'], lower_is_better=True)
    
    # Calculate agreement
    score_diff = np.abs(vina_norm - cnn_norm)
    agreement_mask = score_diff <= threshold
    
    agreement_results = {
        'threshold': threshold,
        'total_comparisons': len(scores_df),
        'agreements': np.sum(agreement_mask),
        'disagreements': np.sum(~agreement_mask),
        'agreement_percentage': (np.sum(agreement_mask) / len(scores_df)) * 100,
        'mean_score_difference': np.mean(score_diff),
        'std_score_difference': np.std(score_diff),
        'agreement_details': {
            'vina_norm': vina_norm,
            'cnn_norm': cnn_norm,
            'score_differences': score_diff,
            'agreement_mask': agreement_mask
        }
    }
    
    print(f"✅ Score agreement analysis completed")
    print(f"   Agreement rate: {agreement_results['agreement_percentage']:.1f}%")
    print(f"   Mean score difference: {agreement_results['mean_score_difference']:.3f}")
    
    return agreement_results

def create_correlation_visualizations(correlation_results: Dict, distribution_results: Dict,
                                    agreement_results: Dict, output_dir: Path) -> List[Path]:
    """
    Create visualizations for correlation analysis.
    
    Parameters
    ----------
    correlation_results : Dict
        Results from correlation analysis
    distribution_results : Dict
        Results from distribution analysis
    agreement_results : Dict
        Results from agreement analysis
    output_dir : Path
        Output directory for visualizations
        
    Returns
    -------
    List[Path]
        List of created visualization files
    """
    print("📊 Creating correlation visualizations...")
    
    output_dir.mkdir(exist_ok=True)
    created_files = []
    
    if 'error' in correlation_results:
        print("⚠️ Skipping correlation visualizations due to missing data")
        return created_files
    
    valid_data = correlation_results['valid_data']
    
    # 1. Correlation matrix heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    correlation_matrix = correlation_results['correlation_matrix']
    
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    sns.heatmap(correlation_matrix, mask=mask, annot=True, cmap='coolwarm', center=0,
                square=True, linewidths=0.5, cbar_kws={"shrink": 0.8}, ax=ax)
    ax.set_title('Correlation Matrix: Vina vs CNN Scores', fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    heatmap_file = output_dir / 'correlation_matrix.png'
    plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
    plt.close()
    created_files.append(heatmap_file)
    
    # 2. Scatter plots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Score Correlations and Distributions', fontsize=16, fontweight='bold')
    
    # Vina vs CNN Affinity
    ax1 = axes[0, 0]
    ax1.scatter(valid_data['vina_affinity'], valid_data['cnn_affinity'], alpha=0.6)
    ax1.set_xlabel('Vina Affinity (kcal/mol)')
    ax1.set_ylabel('CNN Affinity')
    ax1.set_title('Vina vs CNN Affinity')
    
    # Add correlation info
    pearson_r = correlation_results['pearson_correlations']['vina_cnn_affinity']['correlation']
    pearson_p = correlation_results['pearson_correlations']['vina_cnn_affinity']['p_value']
    ax1.text(0.05, 0.95, f'Pearson r = {pearson_r:.3f}\np = {pearson_p:.3f}', 
             transform=ax1.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax1.grid(True, alpha=0.3)
    
    # Vina vs CNN Score
    ax2 = axes[0, 1]
    ax2.scatter(valid_data['vina_affinity'], valid_data['cnn_score'], alpha=0.6, color='orange')
    ax2.set_xlabel('Vina Affinity (kcal/mol)')
    ax2.set_ylabel('CNN Score')
    ax2.set_title('Vina vs CNN Score')
    
    # Add correlation info
    pearson_r = correlation_results['pearson_correlations']['vina_cnn_score']['correlation']
    pearson_p = correlation_results['pearson_correlations']['vina_cnn_score']['p_value']
    ax2.text(0.05, 0.95, f'Pearson r = {pearson_r:.3f}\np = {pearson_p:.3f}', 
             transform=ax2.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax2.grid(True, alpha=0.3)
    
    # Score distributions
    ax3 = axes[1, 0]
    ax3.hist(valid_data['vina_affinity'], bins=20, alpha=0.7, color='blue', label='Vina', density=True)
    if 'cnn_affinity' in valid_data.columns:
        ax3.hist(valid_data['cnn_affinity'], bins=20, alpha=0.7, color='red', label='CNN', density=True)
    ax3.set_xlabel('Affinity')
    ax3.set_ylabel('Density')
    ax3.set_title('Score Distributions')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Agreement analysis
    ax4 = axes[1, 1]
    if 'error' not in agreement_results:
        score_diffs = agreement_results['agreement_details']['score_differences']
        ax4.hist(score_diffs, bins=20, alpha=0.7, color='green', edgecolor='black')
        ax4.axvline(agreement_results['threshold'], color='red', linestyle='--', 
                   label=f'Threshold: {agreement_results["threshold"]}')
        ax4.set_xlabel('Normalized Score Difference')
        ax4.set_ylabel('Frequency')
        ax4.set_title('Score Agreement Analysis')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    scatter_file = output_dir / 'correlation_scatter_plots.png'
    plt.savefig(scatter_file, dpi=300, bbox_inches='tight')
    plt.close()
    created_files.append(scatter_file)
    
    # 3. Summary statistics plot
    if distribution_results and 'descriptive_stats' in distribution_results:
        fig, ax = plt.subplots(figsize=(12, 8))
        
        stats_data = distribution_results['descriptive_stats']
        score_types = list(stats_data.keys())
        means = [stats_data[score]['mean'] for score in score_types]
        stds = [stats_data[score]['std'] for score in score_types]
        
        x_pos = np.arange(len(score_types))
        bars = ax.bar(x_pos, means, yerr=stds, capsize=5, alpha=0.7)
        ax.set_xlabel('Score Type')
        ax.set_ylabel('Value')
        ax.set_title('Score Statistics Summary')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(score_types, rotation=45)
        ax.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for i, (mean, std) in enumerate(zip(means, stds)):
            ax.text(i, mean + std + 0.1, f'{mean:.2f}±{std:.2f}', 
                   ha='center', va='bottom')
        
        plt.tight_layout()
        stats_file = output_dir / 'score_statistics.png'
        plt.savefig(stats_file, dpi=300, bbox_inches='tight')
        plt.close()
        created_files.append(stats_file)
    
    print(f"✅ Created {len(created_files)} correlation visualizations")
    return created_files


def compute_cross_engine_rank_correlations(best_by_engine: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute per-protein and global cross-engine correlations from one-best-pose rows.

    Parameters
    ----------
    best_by_engine : pd.DataFrame
        Table containing at least: engine, tag, protein, affinity_kcal_mol

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        (per_protein_correlations, global_correlations)
    """
    required = {"engine", "tag", "protein", "affinity_kcal_mol"}
    result_columns = [
        "engine_a",
        "engine_b",
        "paired_tags",
        "spearman_rank_corr",
        "spearman_p_value",
        "spearman_q_value",
        "spearman_significant",
        "pearson_affinity_corr",
        "pearson_p_value",
        "pearson_q_value",
        "pearson_significant",
        "insufficient_n",
        "low_power",
        "minimum_required_n",
        "normalization_method",
    ]
    per_protein_columns = ["protein"] + result_columns
    if best_by_engine is None or best_by_engine.empty or not required.issubset(best_by_engine.columns):
        return pd.DataFrame(columns=per_protein_columns), pd.DataFrame(columns=result_columns)

    frame = best_by_engine.copy()
    frame["engine"] = frame["engine"].astype(str)
    frame["tag"] = frame["tag"].astype(str)
    frame["protein"] = frame["protein"].astype(str)
    frame["affinity_kcal_mol"] = pd.to_numeric(frame["affinity_kcal_mol"], errors="coerce")
    frame = frame[np.isfinite(frame["affinity_kcal_mol"])].copy()
    if frame.empty:
        return pd.DataFrame(columns=per_protein_columns), pd.DataFrame(columns=result_columns)

    def _pairwise_corr_rows(pivot_df: pd.DataFrame, protein: Optional[str] = None) -> List[Dict[str, object]]:
        rows: List[Dict[str, object]] = []
        engines = sorted([str(column) for column in pivot_df.columns])
        for index, engine_a in enumerate(engines):
            for engine_b in engines[index + 1 :]:
                pair_frame = pivot_df[[engine_a, engine_b]].dropna()
                paired_tags = int(len(pair_frame))
                row: Dict[str, object] = {
                    "engine_a": engine_a,
                    "engine_b": engine_b,
                    "paired_tags": paired_tags,
                    "spearman_rank_corr": np.nan,
                    "spearman_p_value": np.nan,
                    "spearman_q_value": np.nan,
                    "spearman_significant": False,
                    "pearson_affinity_corr": np.nan,
                    "pearson_p_value": np.nan,
                    "pearson_q_value": np.nan,
                    "pearson_significant": False,
                    "insufficient_n": paired_tags < MIN_CORRELATION_N,
                    "low_power": MIN_CORRELATION_N <= paired_tags < LOW_POWER_N,
                    "minimum_required_n": MIN_CORRELATION_N,
                    "normalization_method": "rank_normalized",
                }
                if protein is not None:
                    row["protein"] = protein
                if paired_tags >= MIN_CORRELATION_N:
                    norm_a = _rank_normalize(pair_frame[engine_a], lower_is_better=True)
                    norm_b = _rank_normalize(pair_frame[engine_b], lower_is_better=True)
                    norm_frame = pd.DataFrame({"a": norm_a, "b": norm_b}).dropna()
                    try:
                        spearman_val, spearman_p = spearmanr(norm_frame["a"], norm_frame["b"])
                    except Exception:
                        spearman_val, spearman_p = (np.nan, np.nan)
                    try:
                        pearson_val, pearson_p = pearsonr(norm_frame["a"], norm_frame["b"])
                    except Exception:
                        pearson_val, pearson_p = (np.nan, np.nan)
                    row["spearman_rank_corr"] = float(spearman_val) if pd.notna(spearman_val) else np.nan
                    row["spearman_p_value"] = float(spearman_p) if pd.notna(spearman_p) else np.nan
                    row["pearson_affinity_corr"] = float(pearson_val) if pd.notna(pearson_val) else np.nan
                    row["pearson_p_value"] = float(pearson_p) if pd.notna(pearson_p) else np.nan
                rows.append(row)
        return rows

    def _apply_bh(df: pd.DataFrame) -> pd.DataFrame:
        if df.empty:
            return df
        for metric in ("spearman", "pearson"):
            p_col = f"{metric}_p_value"
            q_col = f"{metric}_q_value"
            sig_col = f"{metric}_significant"
            valid_mask = (~df["insufficient_n"].astype(bool)) & pd.to_numeric(df[p_col], errors="coerce").notna()
            pvals = pd.to_numeric(df.loc[valid_mask, p_col], errors="coerce").tolist()
            corrected = _bh_fdr_adjust(pvals)
            df[q_col] = np.nan
            if len(corrected) > 0:
                df.loc[valid_mask, q_col] = corrected
            df[sig_col] = False
            df.loc[valid_mask, sig_col] = pd.to_numeric(df.loc[valid_mask, q_col], errors="coerce") < 0.05
        return df

    per_protein_rows: List[Dict[str, object]] = []
    for protein, group in frame.groupby("protein", dropna=False):
        pivot = group.pivot_table(index="tag", columns="engine", values="affinity_kcal_mol", aggfunc="min")
        per_protein_rows.extend(_pairwise_corr_rows(pivot, protein=str(protein)))

    global_pivot = frame.pivot_table(index="tag", columns="engine", values="affinity_kcal_mol", aggfunc="min")
    global_rows = _pairwise_corr_rows(global_pivot, protein=None)

    per_protein_df = pd.DataFrame(per_protein_rows, columns=per_protein_columns)
    global_df = pd.DataFrame(global_rows, columns=result_columns)
    per_protein_df = _apply_bh(per_protein_df)
    global_df = _apply_bh(global_df)
    return per_protein_df, global_df
