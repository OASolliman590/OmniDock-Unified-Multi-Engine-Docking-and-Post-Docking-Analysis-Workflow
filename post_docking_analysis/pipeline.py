"""
Main pipeline module for post-docking analysis.
"""
import sys
from pathlib import Path
import pandas as pd
import logging

# Add the docking_analysis directory to the path so we can import its scripts
docking_analysis_path = Path(__file__).parent.parent / "docking_analysis"
sys.path.insert(0, str(docking_analysis_path))

# Use relative imports for same-package modules
from .config_manager import load_config
from .input_handler import find_docking_files, validate_complex_files
from .docking_parser import parse_all_docking_results
from .structure_quality import assess_structure_quality, create_quality_visualizations
from .correlation_analyzer import analyze_vina_cnn_correlation, analyze_score_distributions, analyze_score_agreement, create_correlation_visualizations
from .pymol_visualizer import PyMOLVisualizer, create_comparative_analysis
from .pymol_generate import render_pymol_scene
from .pandamap_integration import PandaMapAnalyzer
from .plip_integration import PLIPAnalyzer, run_plip_analysis
from .logging_config import setup_logging, get_logger

class PostDockingAnalysisPipeline:
    """
    Main pipeline for post-docking analysis.
    
    This pipeline handles the complete analysis of molecular docking results,
    from parsing output files to generating comprehensive reports.
    """
    
    def __init__(self, input_dir="", output_dir="", config_file=None):
        """
        Initialize the pipeline.
        
        Parameters
        ----------
        input_dir : str
            Path to the directory containing docking results
        output_dir : str
            Path to the directory where results will be saved
        config_file : str, optional
            Path to configuration file (YAML or JSON)
        """
        # Load configuration
        self.config = load_config(config_file)
        
        # Override config with command line arguments if provided
        if input_dir:
            self.config.set("paths.input_dir", input_dir)
        if output_dir:
            self.config.set("paths.output_dir", output_dir)
        
        # Set paths from configuration
        input_dir_value = str(self.config.get("paths.input_dir", "") or "").strip()
        self.input_dir = Path(input_dir_value).expanduser().resolve() if input_dir_value else None
        self.output_dir = Path(self.config.get("paths.output_dir")).expanduser().resolve()

        receptors_dir_value = str(self.config.get("paths.receptors_dir", "") or "").strip()
        self.receptors_dir = (
            Path(receptors_dir_value).expanduser().resolve() if receptors_dir_value else None
        )

        gnina_out_dir_value = str(self.config.get("paths.gnina_out_dir", "") or "").strip()
        self.gnina_out_dir = (
            Path(gnina_out_dir_value).expanduser().resolve() if gnina_out_dir_value else None
        )
        
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up logging
        log_file = self.output_dir / "post_docking_analysis.log"
        log_level = self.config.get("advanced.log_level", "INFO")
        self.logger = setup_logging(str(log_file), log_level)
        self.logger.info("Post-Docking Analysis Pipeline initialized")
        
        # Initialize results storage
        self.results = {}
        self.complexes = []
        self.docking_results = {}
        
        # Initialize plugin manager
        self.plugin_manager = None

    @staticmethod
    def _infer_strong_binder_threshold(best_poses: pd.DataFrame, configured_threshold="auto") -> float:
        """Infer strong-binder threshold from best-pose affinities with legacy-compatible defaults."""
        values = pd.to_numeric(best_poses.get("vina_affinity"), errors="coerce")
        values = values[values.notna()]
        if values.empty:
            return -8.0

        if configured_threshold == "auto":
            q25 = float(values.quantile(0.25))
            q75 = float(values.quantile(0.75))
            iqr = q75 - q25
            return float(q25 - 1.5 * iqr)

        try:
            return float(configured_threshold)
        except (TypeError, ValueError):
            return -8.0

    @staticmethod
    def _add_protein_ligand_columns(best_poses: pd.DataFrame) -> pd.DataFrame:
        """Attach protein/ligand decomposition columns using the legacy complex-name convention."""
        annotated = best_poses.copy()
        parts = annotated["complex_name"].astype(str).str.split("_")
        annotated["protein"] = parts.str[0].fillna("Unknown")
        annotated["ligand"] = parts.apply(lambda values: "_".join(values[1:]) if isinstance(values, list) and len(values) > 1 else "Unknown")
        return annotated

    def _analyze_affinity_dataframe(self, full_df: pd.DataFrame, comparative_benchmark: str, configured_threshold="auto") -> dict:
        """Legacy-compatible affinity analysis that does not depend on removed stub modules."""
        df = full_df.copy()
        df["vina_affinity"] = pd.to_numeric(df.get("vina_affinity"), errors="coerce")
        df = df[df["vina_affinity"].notna()]
        df = df[df["vina_affinity"] < 0]
        if comparative_benchmark != "*":
            mask = df["complex_name"].astype(str).str.contains(comparative_benchmark, case=False, na=False)
            df = df[mask]
        if df.empty:
            return {
                "full_data": pd.DataFrame(),
                "best_poses": pd.DataFrame(),
                "summary_stats": pd.DataFrame(),
                "top_overall": pd.DataFrame(),
                "best_per_protein": pd.DataFrame(),
                "best_per_ligand": pd.DataFrame(),
                "strong_binder_threshold": -8.0,
            }

        best_poses = df.loc[df.groupby("complex_name")["vina_affinity"].idxmin()].copy()
        best_poses = best_poses.sort_values("vina_affinity")
        best_poses = self._add_protein_ligand_columns(best_poses)
        threshold = self._infer_strong_binder_threshold(best_poses, configured_threshold)
        best_poses["binder_category"] = pd.cut(
            best_poses["vina_affinity"],
            bins=[-float("inf"), threshold, -6.0, 0],
            labels=["Strong binder", "Moderate binder", "Weak binder"],
        )

        summary_stats = df.groupby("complex_name").agg({"vina_affinity": ["min", "max", "mean", "std"]}).round(3)
        summary_stats.columns = ["_".join(col).strip() for col in summary_stats.columns]
        summary_stats = summary_stats.reset_index()

        top_cols = [col for col in ["complex_name", "protein", "ligand", "vina_affinity", "binder_category", "pose"] if col in best_poses.columns]
        top_overall = best_poses.head(10)[top_cols]
        best_per_protein = best_poses.groupby("protein").first().reset_index().sort_values("vina_affinity")
        best_per_ligand = best_poses.groupby("ligand").first().reset_index().sort_values("vina_affinity")

        return {
            "full_data": df,
            "best_poses": best_poses,
            "summary_stats": summary_stats,
            "top_overall": top_overall,
            "best_per_protein": best_per_protein,
            "best_per_ligand": best_per_ligand,
            "strong_binder_threshold": float(threshold),
        }
        
    def validate_input(self):
        """
        Validate the input directory and files.
        
        Returns
        -------
        bool
            True if input is valid, False otherwise
        """
        if self.input_dir is None:
            self.logger.error("❌ Input directory is required and cannot be empty")
            return False

        if not self.input_dir.exists():
            self.logger.error(f"❌ Input directory does not exist: {self.input_dir}")
            return False
            
        if not self.input_dir.is_dir():
            self.logger.error(f"❌ Input path is not a directory: {self.input_dir}")
            return False
            
        # Validate GNINA directories if GNINA analysis is enabled
        if "gnina" in self.config.get("analysis.docking_types", []) or self.config.get("analysis.docking_types") == "gnina":
            if self.receptors_dir and not self.receptors_dir.exists():
                self.logger.warning(f"⚠️  Receptors directory does not exist: {self.receptors_dir}")
            if self.gnina_out_dir and not self.gnina_out_dir.exists():
                self.logger.warning(f"⚠️  GNINA output directory does not exist: {self.gnina_out_dir}")
            
        return True
        
    def run_pipeline(self):
        """
        Run the complete post-docking analysis pipeline.
        """
        self.logger.info("🚀 Starting Post-Docking Analysis Pipeline")
        self.logger.info(f"📂 Input directory: {self.input_dir.absolute() if self.input_dir else '(unset)'}")
        self.logger.info(f"📂 Output directory: {self.output_dir.absolute()}")

        try:
            # Validate input
            if not self.validate_input():
                return False

            docking_types = self.config.get("analysis.docking_types", [])
            if isinstance(docking_types, str):
                docking_types = [docking_types]

            # Check if GNINA analysis is enabled
            if "gnina" in docking_types:
                # Check if all_scores.csv exists, if not, try to generate it
                gnina_scores = self.gnina_out_dir / "all_scores.csv" if self.gnina_out_dir else self.input_dir / "gnina_out" / "all_scores.csv"
                
                if not gnina_scores.exists():
                    self.logger.info("⚠️  all_scores.csv not found, attempting to generate from log files...")
                    if self._generate_all_scores_csv():
                        self.logger.info("✅ Successfully generated all_scores.csv")
                        gnina_scores = self.gnina_out_dir / "all_scores.csv" if self.gnina_out_dir else self.input_dir / "gnina_out" / "all_scores.csv"
                    else:
                        self.logger.warning("⚠️  Failed to generate all_scores.csv, continuing with regular analysis...")
                
                # Fast-path: GNINA results already aggregated to CSV
                if gnina_scores.exists():
                    self.logger.info("⚡ Detected GNINA scores CSV. Using streamlined analysis path.")
                    if not self._analyze_from_gnina_scores(gnina_scores):
                        return False
                    # Reports
                    if not self.generate_reports():
                        return False
                    # Visualizations
                    if self.config.get("analysis.generate_visualizations", True):
                        if not self.generate_visualizations():
                            return False
                    # Best poses extraction
                    if self.config.get("analysis.extract_poses", True):
                        if not self.extract_best_poses_pdb():
                            return False
                    # PandaMap interaction analysis
                    if self.config.get("visualization.generate_2d_interactions", True):
                        if not self.generate_pandamap_interactions():
                            return False
                    # PLIP 2D interaction diagrams
                    if self.config.get("visualization.plip.enabled", True):
                        if not self.generate_plip_diagrams():
                            return False
                    # Execute plugins
                    if self.config.get("advanced.enable_plugins", True):
                        # Initialize plugins
                        self.initialize_plugins()
                        # Execute plugins
                        if not self.execute_plugins():
                            return False
                    self.logger.info("✅ Post-Docking Analysis (GNINA fast-path) completed successfully!")
                    return True

            self.logger.info("🔄 Running standard analysis pipeline...")
            return self._run_standard_pipeline()
                
        except Exception as e:
            self.logger.error(f"❌ Pipeline execution failed: {e}")
            return False

    def _run_standard_pipeline(self):
        """
        Run the legacy multi-step analysis path for non-GNINA or non-fast-path input.
        """
        structure_type = self.config.get("advanced.directory_structure", "AUTO")
        self.complexes = find_docking_files(self.input_dir, structure_type)
        self.complexes = validate_complex_files(self.complexes)
        if not self.complexes:
            self.logger.error("❌ No valid docking complexes found for standard analysis")
            return False

        if not self.parse_docking_results():
            return False
        if not self.analyze_binding_affinities():
            return False

        if (
            self.config.get("binding_affinity.analyze_by_protein", True)
            or self.config.get("binding_affinity.analyze_by_ligand", True)
        ):
            if not self.analyze_protein_ligand_breakdown():
                return False

        if self.config.get("analysis.rmsd_analysis", True):
            if not self.analyze_rmsd_and_clustering():
                return False

        if not self.generate_reports():
            return False

        if self.config.get("analysis.generate_visualizations", True):
            if not self.generate_visualizations():
                return False

        if self.config.get("analysis.extract_poses", True):
            if not self.extract_best_poses_pdb():
                return False

        if self.config.get("visualization.generate_2d_interactions", True):
            if not self.generate_pandamap_interactions():
                return False

        if self.config.get("visualization.plip.enabled", True):
            if not self.generate_plip_diagrams():
                return False

        if self.config.get("advanced.enable_plugins", True):
            self.initialize_plugins()
            if not self.execute_plugins():
                return False

        self.logger.info("✅ Post-Docking Analysis (standard pipeline) completed successfully!")
        return True
    
    def _analyze_from_gnina_scores(self, scores_csv: Path):
        """
        Streamlined analysis when GNINA all_scores.csv is available.
        Populates self.results with full_data, best_poses, summary_stats, top_overall.
        """
        try:
            self.logger.info("🔍 Loading GNINA scores CSV...")
            df = pd.read_csv(scores_csv)
            if df.empty:
                self.logger.error("❌ GNINA scores CSV is empty")
                return False
            # Normalize columns
            required = {'tag', 'mode', 'vina_affinity'}
            if not required.issubset(set(df.columns)):
                self.logger.error("❌ GNINA scores CSV missing required columns")
                return False
            # Create a clean DataFrame with proper column names
            full_df = pd.DataFrame({
                'complex_name': df['tag'],
                'pose': df['mode'],
                'vina_affinity': df['vina_affinity']
            })
            # Best poses per tag
            best_poses = df.loc[df.groupby('tag')['vina_affinity'].idxmin()].copy()
            best_poses = best_poses.sort_values('vina_affinity')
            best_poses.rename(columns={'tag': 'complex_name', 'mode': 'pose'}, inplace=True)
            best_poses = best_poses[['complex_name', 'pose', 'vina_affinity']]
            # Summary
            summary_stats = full_df.groupby('complex_name').agg({
                'vina_affinity': ['min', 'max', 'mean', 'std']
            }).round(3)
            summary_stats.columns = ['_'.join(col).strip() for col in summary_stats.columns]
            summary_stats = summary_stats.reset_index()
            # Top
            top_overall = best_poses.head(10)[['complex_name', 'vina_affinity', 'pose']]
            
            comparative_benchmark = self.config.get("analysis.comparative_benchmark", "*")
            strong_binder_threshold = self.config.get("binding_affinity.strong_binder_threshold", "auto")
            analysis_results = self._analyze_affinity_dataframe(full_df, comparative_benchmark, strong_binder_threshold)
            
            self.results = {
                'full_data': full_df,
                'best_poses': best_poses,
                'summary_stats': summary_stats,
                'top_overall': top_overall,
                'strong_binder_threshold': analysis_results['strong_binder_threshold']
            }
            self.logger.info(f"✅ GNINA scores loaded: {len(full_df)} poses, {len(best_poses)} complexes")
            return True
        except Exception as e:
            self.logger.error(f"❌ Error reading GNINA scores: {e}", exc_info=True)
            return False
            
    def _generate_all_scores_csv(self):
        """
        Generate all_scores.csv from GNINA log files if it doesn't exist.
        
        Returns
        -------
        bool
            True if generation successful, False otherwise
        """
        try:
            self.logger.info("🔄 Attempting to generate all_scores.csv from log files...")
            
            # Import the generate_scores_csv function
            from .generate_scores_csv import generate_all_scores_csv
            
            # Determine GNINA output directory
            gnina_out_dir = self.gnina_out_dir if self.gnina_out_dir else self.input_dir / "gnina_out"
            
            if not gnina_out_dir.exists():
                self.logger.error(f"❌ GNINA output directory not found: {gnina_out_dir}")
                return False
            
            # Detect separate logs/ directory (GNINA HPC layout)
            log_dir = None
            sibling_logs = gnina_out_dir.parent / "logs"
            if sibling_logs.is_dir() and any(sibling_logs.glob("*.log")):
                log_dir = sibling_logs
                self.logger.info(f"📂 Detected HPC layout — logs in {log_dir}")
                
            # Look for pairlist.csv in the project directory
            pairlist_file = self.input_dir / "pairlist.csv"
            if not pairlist_file.exists():
                # Try parent directory
                pairlist_file = self.input_dir.parent / "pairlist.csv"
                
            # Generate the CSV file
            output_file = gnina_out_dir / "all_scores.csv"
            success = generate_all_scores_csv(
                gnina_out_dir, output_file,
                pairlist_file if pairlist_file.exists() else None,
                log_dir=log_dir
            )
            
            if success:
                self.logger.info(f"✅ Successfully generated {output_file}")
                return True
            else:
                self.logger.error("❌ Failed to generate all_scores.csv")
                return False
                
        except Exception as e:
            self.logger.error(f"❌ Error generating all_scores.csv: {e}", exc_info=True)
            return False
        
    def parse_docking_results(self):
        """
        Parse docking result files and extract scores.
        
        Returns
        -------
        bool
            True if parsing successful, False otherwise
        """
        self.logger.info("🔍 Parsing docking results...")
        
        try:
            # Parse all docking results
            self.docking_results = parse_all_docking_results(self.complexes)
            
            if not self.docking_results:
                self.logger.error("❌ No docking results parsed successfully")
                return False
                
            self.logger.info(f"✅ Parsed docking results for {len(self.docking_results)} complexes")
            return True
        except Exception as e:
            self.logger.error(f"❌ Error parsing docking results: {e}", exc_info=True)
            return False
        
    def analyze_binding_affinities(self):
        """
        Analyze binding affinities and identify top performers.
        
        Returns
        -------
        bool
            True if analysis successful, False otherwise
        """
        print("📊 Analyzing binding affinities...")
        
        # Get comparative benchmark and strong binder threshold from configuration
        comparative_benchmark = self.config.get("analysis.comparative_benchmark", "*")
        strong_binder_threshold = self.config.get("binding_affinity.strong_binder_threshold", "auto")
        
        # Create a comprehensive DataFrame with all results
        all_data = []
        for complex_name, df in self.docking_results.items():
            for _, row in df.iterrows():
                all_data.append({
                    'complex_name': complex_name,
                    'pose': row['pose'],
                    'vina_affinity': row['vina_affinity'],
                    'rmsd_lb': row['rmsd_lb'],
                    'rmsd_ub': row['rmsd_ub']
                })
        
        if not all_data:
            print("❌ No data to analyze")
            return False
            
        full_df = pd.DataFrame(all_data)
        
        # Analyze binding affinities with comparative benchmark and dynamic threshold
        analysis_results = self._analyze_affinity_dataframe(full_df, comparative_benchmark, strong_binder_threshold)
        
        self.results = analysis_results
        
        print(f"✅ Binding affinities analyzed for {len(full_df)} poses across {len(analysis_results['best_poses'])} complexes")
        print(f"   Best binding affinity: {analysis_results['best_poses']['vina_affinity'].min():.2f} kcal/mol")
        print(f"   Strong binder threshold: {analysis_results['strong_binder_threshold']:.2f} kcal/mol")
        return True
        
    def generate_reports(self):
        """
        Generate summary reports in various formats.
        
        Returns
        -------
        bool
            True if report generation successful, False otherwise
        """
        print("📝 Generating summary reports...")
        
        # Create output directory
        reports_dir = self.output_dir / "reports"
        reports_dir.mkdir(exist_ok=True)
        
        # Generate individual CSV files for each result type
        for name, df in self.results.items():
            if isinstance(df, pd.DataFrame):
                csv_file = reports_dir / f"{name}.csv"
                df.to_csv(csv_file, index=False)
                print(f"✅ {name} report saved to: {csv_file}")
        
        # Generate Excel report
        try:
            import openpyxl
            excel_file = reports_dir / "docking_analysis_results.xlsx"
            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                for name, df in self.results.items():
                    if isinstance(df, pd.DataFrame):
                        # Limit sheet name to 31 characters (Excel limit)
                        sheet_name = name[:31]
                        df.to_excel(writer, sheet_name=sheet_name, index=False)
            print(f"✅ Excel report saved to: {excel_file}")
        except ImportError:
            print("⚠️  openpyxl not available - Excel report generation skipped")
        
        # Generate summary report
        best_poses = self.results['best_poses'].copy()
        if 'vina_affinity' in best_poses.columns:
            best_poses['vina_affinity'] = pd.to_numeric(best_poses['vina_affinity'], errors='coerce')
            best_poses = best_poses[best_poses['vina_affinity'].notna()].sort_values('vina_affinity')
        summary_lines = [
            "Post-Docking Analysis Summary Report",
            "==================================",
            "",
            f"Total complexes analyzed: {len(best_poses)}",
            f"Average binding affinity: {best_poses['vina_affinity'].mean():.2f} kcal/mol",
            f"Best binding affinity: {best_poses['vina_affinity'].min():.2f} kcal/mol",
            f"Worst binding affinity: {best_poses['vina_affinity'].max():.2f} kcal/mol",
            "",
            "Top 5 Performers:",
        ]
        
        # Add top 5 performers
        top_5 = best_poses.head(5)
        for idx, (_, row) in enumerate(top_5.iterrows(), 1):
            complex_name = row['complex_name'] if isinstance(row['complex_name'], str) else str(row['complex_name'])
            vina_affinity = pd.to_numeric(row.get('vina_affinity'), errors='coerce')
            if pd.isna(vina_affinity):
                continue
            pose = pd.to_numeric(row.get('pose', 1), errors='coerce')
            pose = int(pose) if pd.notna(pose) else 1
            
            summary_lines.append(
                f"  {idx}. {complex_name}: {vina_affinity:.2f} kcal/mol (Pose {pose})"
            )
        
        # Save summary report
        summary_file = reports_dir / "summary_report.txt"
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write('\n'.join(summary_lines))
        print(f"✅ Summary report saved to: {summary_file}")
        
        print("✅ Summary reports generated successfully!")
        return True
        
    def generate_visualizations(self):
        """
        Generate visualizations of the results focused on best poses only.
        
        Returns
        -------
        bool
            True if visualization generation successful, False otherwise
        """
        print("📊 Generating visualizations (best poses only)...")
        
        # Check if matplotlib is available
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError:
            print("⚠️  matplotlib/seaborn not available - visualization generation skipped")
            return True
        
        # Create output directory
        viz_dir = self.output_dir / "visualizations"
        viz_dir.mkdir(exist_ok=True)
        
        # Generate all visualizations using the enhanced module
        try:
            from . import visualizer
            plot_files = visualizer.generate_all_visualizations(self.results, self.output_dir, self.config.config)
            print(f"✅ Generated {len(plot_files)} visualizations successfully!")
            return True
        except Exception as e:
            print(f"❌ Error generating visualizations: {e}")
            return False
        
    def extract_best_poses_pdb(self):
        """
        Extract best poses as PDB files by combining receptor and ligand coordinates.
        
        Returns
        -------
        bool
            True if extraction successful, False otherwise
        """
        print("🔬 Extracting best poses as PDB files...")
        
        # Prefer GNINA SDF-based extraction when a gnina_out folder exists
        try:
            from . import pose_extractor
            # Prefer configured GNINA output directory when provided.
            possible_gnina_dirs = []
            if self.gnina_out_dir is not None:
                possible_gnina_dirs.append(self.gnina_out_dir)
            possible_gnina_dirs.extend(
                [
                    self.input_dir / "gnina_out",
                    self.input_dir / "gnina_out_cox2",
                    self.input_dir / "gnina_out_inha",
                ]
            )
            
            gnina_dir = None
            for possible_dir in possible_gnina_dirs:
                if possible_dir.exists():
                    gnina_dir = possible_dir
                    break
            
            if gnina_dir is not None:
                written = pose_extractor.extract_best_poses_from_gnina(
                    self.input_dir,
                    self.output_dir,
                    self.config.config,
                    gnina_dir=gnina_dir,
                    receptors_dir=self.receptors_dir,
                )
                if written > 0:
                    # Organize poses by affinity
                    best_poses_dir = self.output_dir / "best_poses"
                    if best_poses_dir.exists():
                        # Use the calculated threshold from the analysis results
                        threshold = self.results.get('strong_binder_threshold', -8.0)
                        pose_extractor.organize_poses_by_affinity(best_poses_dir, threshold)
                        pose_extractor.create_pose_summary_report(best_poses_dir, self.output_dir / "reports")
                        
                        # Create best binding poses summary folder
                        self._create_best_binding_poses_summary(best_poses_dir, threshold)
                    return True
                else:
                    print("⚠️  GNINA-based extraction wrote 0 files; falling back to PDBQT extraction")
        except Exception as e:
            print(f"⚠️  GNINA-based extraction unavailable: {e}")
        
        # Fallback: extract from Vina PDBQT outputs if available
        try:
            from openbabel import pybel  # noqa: F401
        except ImportError:
            print("⚠️  Open Babel not available - PDB extraction skipped")
            return True
        
        poses_dir = self.output_dir / "best_poses_pdb"
        poses_dir.mkdir(exist_ok=True)
        
        best_poses = self.results['best_poses']
        extracted_count = 0
        
        for _, row in best_poses.iterrows():
            complex_name = row['complex_name']
            pose_number = int(row['pose'])
            
            complex_info = None
            for comp in self.complexes:
                if comp['name'] == complex_name:
                    complex_info = comp
                    break
            
            if not complex_info:
                print(f"⚠️  Complex info not found for {complex_name}")
                continue
                
            try:
                pdb_content = self._extract_pose_from_pdbqt(
                    complex_info['docking_result'], 
                    pose_number,
                    complex_info.get('receptor'),
                    complex_name
                )
                
                if pdb_content:
                    pdb_file = poses_dir / f"{complex_name}_pose{pose_number}.pdb"
                    with open(pdb_file, 'w', encoding='utf-8') as f:
                        f.write(pdb_content)
                    extracted_count += 1
                    print(f"✅ Extracted {complex_name} pose {pose_number}")
                else:
                    print(f"⚠️  Failed to extract {complex_name} pose {pose_number}")
                    
            except Exception as e:
                print(f"❌ Error extracting {complex_name} pose {pose_number}: {e}")
        
        print(f"✅ Extracted {extracted_count} best poses as PDB files to: {poses_dir}")
        
        # Optional: auto-render PyMOL PNGs for each extracted pose if PyMOL is available
        try:
            rendered_dir = self.output_dir / "pymol_renders"
            for pdb_file in poses_dir.glob("*.pdb"):
                render_pymol_scene(pdb_file, rendered_dir, pdb_file.stem)
        except Exception as e:
            print(f"⚠️  Skipping PyMOL auto-render: {e}")
        return True
        
    def _extract_pose_from_pdbqt(self, pdbqt_file, pose_number, receptor_file, complex_name):
        """
        Extract a specific pose from PDBQT file and combine with receptor.
        
        Parameters
        ----------
        pdbqt_file : Path
            Path to PDBQT file
        pose_number : int
            Pose number to extract (1-based)
        receptor_file : Path
            Path to receptor PDBQT file
        complex_name : str
            Name of the complex
            
        Returns
        -------
        str
            PDB content as string, or None if failed
        """
        try:
            from openbabel import pybel
            
            # Read the PDBQT file and extract the specific pose
            poses = list(pybel.readfile("pdbqt", str(pdbqt_file)))
            
            if pose_number > len(poses):
                print(f"⚠️  Pose {pose_number} not found in {pdbqt_file}")
                return None
                
            ligand_pose = poses[pose_number - 1]  # Convert to 0-based index
            
            # Convert ligand to PDB format
            ligand_pdb = ligand_pose.write("pdb")
            ligand_lines = []
            for line in ligand_pdb.split('\n'):
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # Fix the line format and assign chain B
                    line = line.ljust(80)
                    new_line = f"HETATM{line[6:21]}B{line[22:]}"
                    new_line = new_line[:17] + "UNK" + new_line[20:]
                    ligand_lines.append(new_line)
            
            # Read receptor if available
            receptor_lines = []
            if receptor_file and receptor_file.exists():
                try:
                    receptor_mol = next(pybel.readfile("pdbqt", str(receptor_file)))
                    receptor_pdb = receptor_mol.write("pdb")
                    for line in receptor_pdb.split('\n'):
                        if line.startswith('ATOM'):
                            # Fix the line format and assign chain A
                            line = line.ljust(80)
                            new_line = f"ATOM  {line[6:21]}A{line[22:]}"
                            receptor_lines.append(new_line)
                except Exception as e:
                    print(f"⚠️  Could not read receptor {receptor_file}: {e}")
            
            # Combine receptor and ligand
            all_lines = receptor_lines + ligand_lines + ["END"]
            return '\n'.join(all_lines)
            
        except Exception as e:
            print(f"❌ Error extracting pose from PDBQT: {e}")
            return None
    
    def analyze_protein_ligand_breakdown(self):
        """
        Analyze best performance by protein and by ligand.
        
        Returns
        -------
        bool
            True if analysis successful, False otherwise
        """
        print("🧬 Analyzing protein vs ligand breakdown...")
        
        if not hasattr(self, 'results') or 'best_poses' not in self.results:
            print("❌ No results available for protein-ligand breakdown")
            return False
        
        try:
            best_poses = self.results['best_poses'].copy()
            if best_poses.empty:
                print("⚠️  Best-pose table is empty; protein-ligand breakdown skipped")
                return True
            best_poses = self._add_protein_ligand_columns(best_poses)

            best_per_protein = best_poses.groupby('protein').agg({
                'vina_affinity': 'min',
                'complex_name': 'first',
                'ligand': 'first',
            }).reset_index().sort_values('vina_affinity')
            best_per_protein.columns = ['protein', 'best_affinity', 'best_complex', 'best_ligand']

            best_per_ligand = best_poses.groupby('ligand').agg({
                'vina_affinity': 'min',
                'complex_name': 'first',
                'protein': 'first',
            }).reset_index().sort_values('vina_affinity')
            best_per_ligand.columns = ['ligand', 'best_affinity', 'best_complex', 'best_protein']

            protein_summary = best_poses.groupby('protein').agg({
                'vina_affinity': ['min', 'max', 'mean', 'std', 'count'],
            }).round(3)
            protein_summary.columns = ['min_affinity', 'max_affinity', 'mean_affinity', 'std_affinity', 'pose_count']
            protein_summary = protein_summary.reset_index()

            ligand_summary = best_poses.groupby('ligand').agg({
                'vina_affinity': ['min', 'max', 'mean', 'std', 'count'],
            }).round(3)
            ligand_summary.columns = ['min_affinity', 'max_affinity', 'mean_affinity', 'std_affinity', 'pose_count']
            ligand_summary = ligand_summary.reset_index()

            breakdown_results = {
                'best_per_protein': best_per_protein,
                'best_per_ligand': best_per_ligand,
                'protein_summary': protein_summary,
                'ligand_summary': ligand_summary,
                'scores_with_breakdown': best_poses,
            }
            
            # Store results
            if not hasattr(self, 'results'):
                self.results = {}
            self.results['protein_ligand_breakdown'] = breakdown_results
            
            # Save breakdown reports
            breakdown_dir = self.output_dir / "reports"
            breakdown_dir.mkdir(exist_ok=True)
            
            breakdown_results['best_per_protein'].to_csv(
                breakdown_dir / "best_per_protein.csv", index=False
            )
            breakdown_results['best_per_ligand'].to_csv(
                breakdown_dir / "best_per_ligand.csv", index=False
            )
            breakdown_results['protein_summary'].to_csv(
                breakdown_dir / "protein_summary.csv", index=False
            )
            breakdown_results['ligand_summary'].to_csv(
                breakdown_dir / "ligand_summary.csv", index=False
            )
            
            print("✅ Protein vs ligand breakdown completed")
            return True
            
        except Exception as e:
            print(f"❌ Error in protein-ligand breakdown: {e}")
            return False
    
    def analyze_rmsd_and_clustering(self):
        """
        Perform RMSD analysis and pose clustering with comparative benchmarking.
        
        Returns
        -------
        bool
            True if analysis successful, False otherwise
        """
        print("📏 Analyzing RMSD and pose clustering...")
        
        if not hasattr(self, 'results') or 'full_data' not in self.results:
            print("❌ No results available for RMSD analysis")
            return False
        
        try:
            from .enhanced_rmsd_analyzer import (
                analyze_conformational_diversity_enhanced,
                analyze_pose_clustering_enhanced,
                calculate_rmsd_matrix_from_pdbs,
                create_rmsd_visualizations_enhanced,
            )

            clustering_method = self.config.get("rmsd.clustering_method", "kmeans")
            n_clusters = int(self.config.get("rmsd.kmeans_clusters", 3))
            rmsd_workers = int(self.config.get("rmsd.workers", 0) or 0)
            poses_dir = self.output_dir / "best_poses_pdb"
            pdb_files = sorted(poses_dir.glob("*.pdb")) if poses_dir.exists() else []
            if len(pdb_files) < 2:
                print("⚠️  RMSD analysis skipped: need at least two best-pose PDB files")
                self.results['rmsd_analysis'] = {
                    'state': 'skipped_insufficient_pdbs',
                    'poses_dir': str(poses_dir),
                    'pdb_count': len(pdb_files),
                }
                return True

            poses_df = self.results.get('best_poses', pd.DataFrame()).copy()
            if poses_df.empty:
                tags = [path.stem for path in pdb_files]
                poses_df = pd.DataFrame({'tag': tags, 'vina_affinity': [np.nan] * len(tags)})
            if 'tag' not in poses_df.columns:
                if 'complex_name' in poses_df.columns:
                    poses_df['tag'] = poses_df['complex_name'].astype(str)
                else:
                    poses_df['tag'] = [path.stem for path in pdb_files[: len(poses_df)]]

            rmsd_matrix, filenames = calculate_rmsd_matrix_from_pdbs(
                pdb_files, ligand_only=True, num_workers=rmsd_workers
            )
            matrix_n = rmsd_matrix.shape[0]
            poses_df = poses_df.head(matrix_n).reset_index(drop=True)
            if len(poses_df) < matrix_n:
                missing = matrix_n - len(poses_df)
                extra = pd.DataFrame({
                    'tag': filenames[len(poses_df):matrix_n],
                    'vina_affinity': [np.nan] * missing,
                })
                poses_df = pd.concat([poses_df, extra], ignore_index=True)

            clustering_results = analyze_pose_clustering_enhanced(
                poses_df,
                rmsd_matrix,
                pdb_files,
                method=clustering_method,
                n_clusters=n_clusters,
            )
            diversity_results = analyze_conformational_diversity_enhanced(
                clustering_results.get('poses_with_clusters', poses_df),
                clustering_results.get('rmsd_matrix', rmsd_matrix),
            )
            
            # Store results
            if not hasattr(self, 'results'):
                self.results = {}
            self.results['rmsd_analysis'] = {
                'clustering': clustering_results,
                'diversity': diversity_results,
                'rmsd_matrix': rmsd_matrix
            }
            
            # Create RMSD visualizations
            viz_dir = self.output_dir / "visualizations"
            viz_dir.mkdir(exist_ok=True)
            create_rmsd_visualizations_enhanced(
                clustering_results, diversity_results, viz_dir,
                dpi=self.config.get("visualization.dpi", 300)
            )
            
            # Save RMSD reports
            reports_dir = self.output_dir / "reports"
            reports_dir.mkdir(exist_ok=True)
            
            clustering_results['poses_with_clusters'].to_csv(
                reports_dir / "poses_with_clusters.csv", index=False
            )
            clustering_results['cluster_summary'].to_csv(
                reports_dir / "cluster_summary.csv", index=False
            )
            clustering_results['cluster_centroids'].to_csv(
                reports_dir / "cluster_centroids.csv", index=False
            )
            diversity_results['diversity_metrics'].to_csv(
                reports_dir / "diversity_metrics.csv", index=False
            )
            
            print("✅ RMSD analysis and clustering completed")
            return True
            
        except Exception as e:
            print(f"❌ Error in RMSD analysis: {e}")
            return False
    
    def assess_structure_quality(self):
        """
        Assess structure quality of best poses.
        
        Returns
        -------
        bool
            True if analysis successful, False otherwise
        """
        print("🔍 Assessing structure quality...")
        
        if not hasattr(self, 'results') or 'best_poses' not in self.results:
            print("❌ No results available for structure quality assessment")
            return False
        
        try:
            # Get best poses PDB files
            poses_dir = self.output_dir / "best_poses_pdb"
            if not poses_dir.exists():
                print("⚠️ Best poses PDB files not found - skipping quality assessment")
                return True
            
            quality_results = []
            pdb_files = list(poses_dir.glob("*.pdb"))
            
            for pdb_file in pdb_files:
                quality_assessment = assess_structure_quality(pdb_file)
                quality_results.append(quality_assessment)
            
            # Store results
            if not hasattr(self, 'results'):
                self.results = {}
            self.results['structure_quality'] = quality_results
            
            # Create quality visualizations
            viz_dir = self.output_dir / "visualizations"
            viz_dir.mkdir(exist_ok=True)
            create_quality_visualizations(quality_results, viz_dir)
            
            # Save quality reports
            reports_dir = self.output_dir / "reports"
            reports_dir.mkdir(exist_ok=True)
            
            quality_summary = []
            for result in quality_results:
                quality_summary.append({
                    'pdb_file': Path(result['pdb_file']).name,
                    'quality_score': result['quality_score'],
                    'overall_quality': result['overall_quality'],
                    'total_clashes': result['clash_data']['total_clashes'],
                    'allowed_residues_pct': result['ramachandran_quality']['allowed_percentage']
                })
            
            quality_df = pd.DataFrame(quality_summary)
            quality_df.to_csv(reports_dir / "structure_quality_summary.csv", index=False)
            
            print("✅ Structure quality assessment completed")
            return True
            
        except Exception as e:
            print(f"❌ Error in structure quality assessment: {e}")
            return False
    
    def analyze_correlations(self):
        """
        Analyze correlations between different scoring functions.
        
        Returns
        -------
        bool
            True if analysis successful, False otherwise
        """
        print("📊 Analyzing score correlations...")
        
        if not hasattr(self, 'results') or 'full_data' not in self.results:
            print("❌ No results available for correlation analysis")
            return False
        
        try:
            # Analyze Vina-CNN correlations
            correlation_results = analyze_vina_cnn_correlation(self.results['full_data'])
            
            # Analyze score distributions
            distribution_results = analyze_score_distributions(self.results['full_data'])
            
            # Analyze score agreement
            agreement_results = analyze_score_agreement(self.results['full_data'])
            
            # Store results
            if not hasattr(self, 'results'):
                self.results = {}
            self.results['correlation_analysis'] = {
                'correlations': correlation_results,
                'distributions': distribution_results,
                'agreement': agreement_results
            }
            
            # Create correlation visualizations
            viz_dir = self.output_dir / "visualizations"
            viz_dir.mkdir(exist_ok=True)
            create_correlation_visualizations(
                correlation_results, distribution_results, agreement_results, viz_dir
            )
            
            # Save correlation reports
            reports_dir = self.output_dir / "reports"
            reports_dir.mkdir(exist_ok=True)
            
            if 'error' not in correlation_results:
                correlation_summary = {
                    'n_samples': correlation_results['n_samples'],
                    'vina_cnn_affinity_pearson': correlation_results['pearson_correlations']['vina_cnn_affinity']['correlation'],
                    'vina_cnn_affinity_p_value': correlation_results['pearson_correlations']['vina_cnn_affinity']['p_value'],
                    'vina_cnn_score_pearson': correlation_results['pearson_correlations']['vina_cnn_score']['correlation'],
                    'vina_cnn_score_p_value': correlation_results['pearson_correlations']['vina_cnn_score']['p_value']
                }
                
                correlation_df = pd.DataFrame([correlation_summary])
                correlation_df.to_csv(reports_dir / "correlation_summary.csv", index=False)
            
            print("✅ Correlation analysis completed")
            return True
            
        except Exception as e:
            print(f"❌ Error in correlation analysis: {e}")
            return False
    
    def create_pymol_visualizations(self):
        """
        Create 3D visualizations using PyMOL.
        
        Returns
        -------
        bool
            True if analysis successful, False otherwise
        """
        print("🎬 Creating PyMOL visualizations...")
        
        try:
            # Create PyMOL visualizer with configuration
            pymol_dir = self.output_dir / "pymol_visualizations"
            visualizer = PyMOLVisualizer(pymol_dir, self.config.config)
            
            # Get best poses PDB files
            poses_dir = self.output_dir / "best_poses_pdb"
            if poses_dir.exists():
                pdb_files = self._rank_pose_files_by_affinity(list(poses_dir.glob("*.pdb")))
                
                if len(pdb_files) >= 2:
                    # Create comparative analysis between the two strongest ranked poses
                    reference_pdb = pdb_files[0]
                    novel_pdb = pdb_files[1]
                    
                    comparative_results = create_comparative_analysis(
                        reference_pdb, novel_pdb, pymol_dir, 
                        highlight_residues=[212, 213, 214],
                        config=self.config.config
                    )
                    
                    # Create best poses gallery
                    gallery_session = visualizer.create_best_poses_gallery(
                        pdb_files[:5], "best_poses_gallery"  # Limit to first 5 poses
                    )
                    
                    # Store results
                    if not hasattr(self, 'results'):
                        self.results = {}
                    self.results['pymol_visualizations'] = {
                        'comparative_analysis': comparative_results,
                        'gallery_session': gallery_session,
                        'output_directory': pymol_dir
                    }
                    
                    print("✅ PyMOL visualizations created")
                    return True
                else:
                    print("⚠️ Need at least 2 PDB files for comparative analysis")
                    return True
            else:
                print("⚠️ Best poses PDB files not found - skipping PyMOL visualizations")
                return True
                
        except Exception as e:
            print(f"❌ Error creating PyMOL visualizations: {e}")
            return False

    def generate_pandamap_interactions(self):
        """
        Generate protein-ligand interaction visualizations using PandaMap.
        
        Returns
        -------
        bool
            True if analysis successful, False otherwise
        """
        print("🐼 Generating PandaMap interaction visualizations...")
        
        try:
            # Check if best poses PDB files exist
            poses_dir = self.output_dir / "best_poses_pdb"
            if not poses_dir.exists():
                print("⚠️ Best poses PDB files not found - skipping PandaMap analysis")
                return True
            
            # Create PandaMap output directory
            pandamap_dir = self.output_dir / "pandamap_analysis"
            
            # Initialize PandaMap analyzer with configuration
            conda_env = self.config.get("advanced", {}).get("pandamap_conda_env", "pandamap")
            analyzer = PandaMapAnalyzer(conda_env=conda_env, config=self.config.config)
            
            # Generate comprehensive analysis
            ligand_name = self._infer_pandamap_ligand_name()
            summary = analyzer.generate_comprehensive_analysis(
                poses_dir=poses_dir,
                output_dir=pandamap_dir,
                ligand_name=ligand_name
            )
            
            # Store results
            if not hasattr(self, 'results'):
                self.results = {}
            self.results['pandamap_analysis'] = summary
            
            print(f"✅ PandaMap analysis completed:")
            print(f"   📊 Generated {summary['generated_2d_maps']} 2D interaction maps")
            print(f"   🌐 Generated {summary['generated_3d_visualizations']} 3D HTML visualizations")
            print(f"   📄 Generated {summary['generated_reports']} detailed reports")
            print(f"   📂 Output directory: {pandamap_dir}")
            
            return True
            
        except Exception as e:
            print(f"❌ Error generating PandaMap interactions: {e}")
            return False

    def generate_plip_diagrams(self):
        """
        Generate 2D protein-ligand interaction diagrams using PLIP.
        
        PLIP (Protein-Ligand Interaction Profiler) generates detailed 2D
        diagrams showing hydrogen bonds, hydrophobic contacts, π-stacking,
        salt bridges, and other interactions.
        
        Returns
        -------
        bool
            True if analysis successful, False otherwise
        """
        print("🔬 Generating PLIP 2D interaction diagrams...")
        
        # Check if PLIP is enabled in configuration
        plip_config = self.config.get("visualization.plip", {})
        if not plip_config.get("enabled", True):
            print("ℹ️  PLIP analysis disabled in configuration")
            return True
        
        try:
            # Determine poses directory - check multiple locations
            poses_dir = None
            
            # Priority 1: best_poses with subdirectories (strong/moderate/weak)
            best_poses_dir = self.output_dir / "best_poses"
            if best_poses_dir.exists():
                # Check if it has the expected subdirectory structure
                subdirs = ['strong_binders', 'moderate_binders', 'weak_binders']
                if any((best_poses_dir / d).exists() for d in subdirs):
                    poses_dir = best_poses_dir
                    print(f"📁 Using categorized poses from: {poses_dir}")
            
            # Priority 2: best_poses_pdb flat directory
            if poses_dir is None:
                best_poses_pdb_dir = self.output_dir / "best_poses_pdb"
                if best_poses_pdb_dir.exists():
                    poses_dir = best_poses_pdb_dir
                    print(f"📁 Using poses from: {poses_dir}")
            
            if poses_dir is None or not poses_dir.exists():
                print("⚠️ No poses directory found - skipping PLIP analysis")
                return True
            
            # Create PLIP output directory
            plip_output_dir = self.output_dir / "plip_diagrams"
            plip_output_dir.mkdir(parents=True, exist_ok=True)
            
            # Initialize PLIP analyzer with configuration
            analyzer = PLIPAnalyzer(config=plip_config)
            try:
                max_workers = max(1, int(plip_config.get("max_workers", 4) or 1))
            except (TypeError, ValueError):
                self.logger.warning(
                    "⚠️ Invalid PLIP max_workers value '%s'; falling back to 1",
                    plip_config.get("max_workers"),
                )
                max_workers = 1
            
            if not analyzer.plip_available:
                print("⚠️ PLIP not installed - skipping 2D diagram generation")
                print("   Install with: pip install plip")
                return True
            
            # Check if this is a categorized structure or flat directory
            subdirs = ['strong_binders', 'moderate_binders', 'weak_binders']
            is_categorized = any((poses_dir / d).exists() for d in subdirs)
            
            if is_categorized:
                # Analyze the complete best_poses structure
                all_results = analyzer.analyze_best_poses_directory(
                    poses_dir, plip_output_dir, max_workers=max_workers
                )
                
                # Calculate totals
                total_analyzed = sum(len(r) for r in all_results.values())
                total_successful = sum(
                    sum(1 for res in r.values() if res.success)
                    for r in all_results.values()
                )
            else:
                # Single directory
                results = analyzer.analyze_directory(poses_dir, plip_output_dir, max_workers=max_workers)
                total_analyzed = len(results)
                total_successful = sum(1 for r in results.values() if r.success)
            
            # Store results
            if not hasattr(self, 'results'):
                self.results = {}
            self.results['plip_diagrams'] = {
                'total_analyzed': total_analyzed,
                'successful': total_successful,
                'output_dir': str(plip_output_dir)
            }
            
            print(f"✅ PLIP analysis completed:")
            print(f"   📊 Analyzed {total_analyzed} complexes")
            print(f"   ✅ Generated {total_successful} 2D interaction diagrams")
            print(f"   📂 Output directory: {plip_output_dir}")
            
            return True
            
        except Exception as e:
            print(f"❌ Error generating PLIP diagrams: {e}")
            import traceback
            traceback.print_exc()
            return False

    def initialize_plugins(self):
        """
        Initialize plugin manager and load available plugins.
        """
        print("ℹ️  Plugin system has been retired; skipping plugin initialization.")
        self.plugin_manager = None
        return True

    def execute_plugins(self):
        """
        Execute all loaded plugins with current results.
        
        Returns
        -------
        bool
            True if plugins executed successfully, False otherwise
        """
        print("ℹ️  Plugin execution skipped: plugin system retired.")
        self.results['plugin_results'] = {
            "status": "skipped_retired_plugin_system",
            "reason": "Plugin manager and plugin bundles were removed during post-docking remediation.",
        }
        return True

    def _load_pose_summary_dataframe(self):
        """Load pose_summary.csv when available for affinity-based ranking."""
        summary_file = self.output_dir / "reports" / "pose_summary.csv"
        if not summary_file.exists():
            return None
        try:
            return pd.read_csv(summary_file)
        except Exception as exc:
            self.logger.warning(f"⚠️  Could not read pose summary for ranking: {exc}")
            return None

    def _infer_pandamap_ligand_name(self) -> str:
        """Infer a more useful ligand name than the legacy UNK fallback."""
        try:
            best_poses = self.results.get("best_poses")
            if isinstance(best_poses, pd.DataFrame) and not best_poses.empty:
                if "ligand" in best_poses.columns:
                    ligand = str(best_poses.iloc[0].get("ligand") or "").strip()
                    if ligand:
                        return ligand
                complex_name = str(best_poses.iloc[0].get("complex_name") or "").strip()
                if complex_name:
                    parts = complex_name.split("_")
                    if len(parts) >= 2:
                        return parts[-1]
        except Exception:
            pass
        return "UNK"

    def _rank_pose_files_by_affinity(self, files):
        """Rank extracted pose files using pose_summary.csv affinity when available."""
        pose_summary_df = self._load_pose_summary_dataframe()
        if pose_summary_df is None or pose_summary_df.empty:
            return sorted(files, key=lambda path: path.name)

        pose_summary_df = pose_summary_df.copy()
        if "pdb_file" in pose_summary_df.columns:
            pose_summary_df["pdb_file"] = pose_summary_df["pdb_file"].astype(str)

        ranked = []
        for pdb_file in files:
            rel_name = pdb_file.name
            match = (
                pose_summary_df[pose_summary_df["pdb_file"].str.endswith(rel_name)]
                if "pdb_file" in pose_summary_df.columns
                else pd.DataFrame()
            )
            if match.empty:
                affinity = float("inf")
            else:
                affinity = pd.to_numeric(match.iloc[0].get("vina_affinity"), errors="coerce")
                if pd.isna(affinity):
                    affinity = float("inf")
            ranked.append((float(affinity), pdb_file.name, pdb_file))
        ranked.sort(key=lambda item: (item[0], item[1]))
        return [item[2] for item in ranked]

    def _create_best_binding_poses_summary(self, best_poses_dir: Path, threshold: float):
        """
        Create a summary folder with the best binding poses from each category.
        
        Parameters
        ----------
        best_poses_dir : Path
            Directory containing organized poses
        threshold : float
            Strong binder threshold
        """
        import shutil
        
        # Create best binding poses summary directory
        summary_dir = self.output_dir / "best_binding_poses"
        summary_dir.mkdir(exist_ok=True)
        
        # Copy top 5 strong binders
        strong_binders_dir = best_poses_dir / "strong_binders"
        if strong_binders_dir.exists():
            strong_files = list(strong_binders_dir.glob("*.pdb"))
            strong_files = self._rank_pose_files_by_affinity(strong_files)
            
            for i, pdb_file in enumerate(strong_files[:5], 1):
                dest_file = summary_dir / f"top_{i}_strong_binder_{pdb_file.name}"
                shutil.copy2(pdb_file, dest_file)
        
        # Copy top 5 moderate binders
        moderate_binders_dir = best_poses_dir / "moderate_binders"
        if moderate_binders_dir.exists():
            moderate_files = list(moderate_binders_dir.glob("*.pdb"))
            moderate_files = self._rank_pose_files_by_affinity(moderate_files)
            
            for i, pdb_file in enumerate(moderate_files[:5], 1):
                dest_file = summary_dir / f"top_{i}_moderate_binder_{pdb_file.name}"
                shutil.copy2(pdb_file, dest_file)
        
        # Create summary report
        summary_file = summary_dir / "best_binding_poses_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("Best Binding Poses Summary\n")
            f.write("==========================\n\n")
            f.write(f"Strong binder threshold: {threshold:.2f} kcal/mol\n\n")
            
            f.write("Strong Binders (≤{:.2f} kcal/mol):\n".format(threshold))
            if strong_binders_dir.exists():
                strong_files = _rank_pose_files(list(strong_binders_dir.glob("*.pdb")))
                for i, pdb_file in enumerate(strong_files[:5], 1):
                    f.write(f"  {i}. {pdb_file.name}\n")
            else:
                f.write("  No strong binders found\n")
            
            f.write(f"\nModerate Binders (-6.0 to {threshold:.2f} kcal/mol):\n")
            if moderate_binders_dir.exists():
                moderate_files = _rank_pose_files(list(moderate_binders_dir.glob("*.pdb")))
                for i, pdb_file in enumerate(moderate_files[:5], 1):
                    f.write(f"  {i}. {pdb_file.name}\n")
            else:
                f.write("  No moderate binders found\n")
        
        print(f"✅ Best binding poses summary created in: {summary_dir}")


def main():
    """
    Main function to run the pipeline from command line.
    """
    print(
        "❌ Direct legacy pipeline entrypoint is disabled.\n"
        "Use:\n"
        "  1) python -m post_docking_analysis --project-dir <project_root> [--analysis-mode ...]\n"
        "  2) python -m post_docking_analysis.simplified_cli --project-dir <gnina_root> --output <out_dir>\n"
    )
    return 1

if __name__ == "__main__":
    sys.exit(main())
