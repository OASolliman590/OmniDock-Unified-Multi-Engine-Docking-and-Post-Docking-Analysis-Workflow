"""
PLIP (Protein-Ligand Interaction Profiler) integration module for post-docking analysis.

Generates 2D interaction diagrams and interaction reports using PLIP for local PDB files.
This provides PoseView-style 2D diagrams showing hydrogen bonds, hydrophobic contacts,
π-stacking, salt bridges, and other protein-ligand interactions.

Reference: https://github.com/pharmai/plip
"""
import subprocess
import shutil
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from concurrent.futures import ThreadPoolExecutor, as_completed
import tempfile

logger = logging.getLogger(__name__)


def _coerce_worker_count(value: Any, default: int = 1) -> int:
    """Coerce configured worker counts to a safe positive integer."""
    try:
        return max(1, int(value or default))
    except (TypeError, ValueError):
        logger.warning("⚠️ Invalid PLIP max_workers value '%s'; falling back to %s", value, default)
        return max(1, int(default))


@dataclass
class PLIPResult:
    """Results from PLIP analysis of a single complex."""
    complex_name: str
    pdb_file: Path
    success: bool
    png_diagram: Optional[Path] = None
    pymol_session: Optional[Path] = None
    xml_report: Optional[Path] = None
    txt_report: Optional[Path] = None
    binding_site_count: int = 0
    interactions: Dict[str, int] = field(default_factory=dict)
    error_message: Optional[str] = None


class PLIPAnalyzer:
    """
    PLIP analyzer for generating 2D interaction diagrams from local PDB files.
    
    This class provides methods to analyze protein-ligand complexes and generate:
    - 2D interaction diagrams (PNG)
    - PyMOL session files (.pse)
    - XML reports with detailed interaction data
    - Text reports for human-readable output
    """
    
    def __init__(self, config: Dict = None):
        """
        Initialize PLIP analyzer.
        
        Parameters
        ----------
        config : Dict, optional
            Configuration dictionary with PLIP settings
        """
        self.config = config or {}
        self._use_python_module = False
        self.plip_available = self._check_plip_installation()
        
        # Default settings
        self.output_formats = self.config.get('output_formats', ['png', 'xml', 'txt'])
        self.generate_pymol = self.config.get('generate_pymol_session', True)
        self.timeout = self.config.get('timeout_seconds', 120)
        
    def _check_plip_installation(self) -> bool:
        """Check if PLIP is installed and accessible."""
        # First try importing as Python module (most reliable)
        try:
            from plip.structure.preparation import PDBComplex
            logger.info("✅ PLIP Python module available")
            self._use_python_module = True
            return True
        except ImportError:
            pass
        
        # Try command line as fallback
        try:
            # Try direct command
            result = subprocess.run(
                ['plip', '--version'],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                logger.info(f"✅ PLIP CLI found: {result.stdout.strip()}")
                self._use_python_module = False
                self._plip_cmd = ['plip']
                return True
        except FileNotFoundError:
            pass
        
        # Try python -m plip.plipcmd
        try:
            import sys
            result = subprocess.run(
                [sys.executable, '-m', 'plip.plipcmd', '--version'],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                logger.info(f"✅ PLIP module CLI found")
                self._use_python_module = False
                self._plip_cmd = [sys.executable, '-m', 'plip.plipcmd']
                return True
        except Exception:
            pass
        
        logger.warning(
            "⚠️ PLIP not found. Install with:\n"
            "   pip install plip\n"
            "   or: conda install -c conda-forge plip"
        )
        return False
    
    def analyze_complex(
        self,
        pdb_file: Path,
        output_dir: Path,
        complex_name: Optional[str] = None
    ) -> PLIPResult:
        """
        Analyze a protein-ligand complex and generate 2D interaction diagram.
        
        Parameters
        ----------
        pdb_file : Path
            Path to the PDB file containing the complex
        output_dir : Path
            Output directory for results
        complex_name : str, optional
            Name for the complex (defaults to filename stem)
            
        Returns
        -------
        PLIPResult
            Results from the PLIP analysis
        """
        pdb_file = Path(pdb_file)
        output_dir = Path(output_dir)
        
        if complex_name is None:
            complex_name = pdb_file.stem
        
        # Create output directory for this complex
        complex_output_dir = output_dir / complex_name
        complex_output_dir.mkdir(parents=True, exist_ok=True)
        
        if not self.plip_available:
            return PLIPResult(
                complex_name=complex_name,
                pdb_file=pdb_file,
                success=False,
                error_message="PLIP not installed"
            )
        
        if not pdb_file.exists():
            return PLIPResult(
                complex_name=complex_name,
                pdb_file=pdb_file,
                success=False,
                error_message=f"PDB file not found: {pdb_file}"
            )
        
        logger.info(f"🔬 Analyzing: {complex_name}")
        
        try:
            # Determine which method to use
            if hasattr(self, '_use_python_module') and self._use_python_module:
                # Use Python API directly
                return self._analyze_with_python_api(pdb_file, complex_output_dir, complex_name)
            
            # Use command line
            plip_cmd = getattr(self, '_plip_cmd', ['plip'])
            
            # Build PLIP command
            cmd = plip_cmd + [
                '-f', str(pdb_file),
                '-o', str(complex_output_dir),
            ]

            requested_formats = {
                str(value).strip().lower()
                for value in (self.output_formats or [])
                if str(value).strip()
            }
            if "xml" in requested_formats:
                cmd.append('-x')
            if "txt" in requested_formats or "text" in requested_formats:
                cmd.append('-t')
            if "png" in requested_formats:
                cmd.append('-p')
            
            # Add PyMOL session generation if enabled
            if self.generate_pymol:
                cmd.append('-y')  # Generate PyMOL session
            
            # Run PLIP
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout
            )
            
            if result.returncode != 0:
                logger.warning(f"⚠️ PLIP returned non-zero: {result.stderr}")
                # Continue anyway - PLIP sometimes returns non-zero but still produces output
            
            # Find generated files
            png_files = list(complex_output_dir.glob("*.png"))
            pse_files = list(complex_output_dir.glob("*.pse"))
            xml_files = list(complex_output_dir.glob("*.xml"))
            txt_files = list(complex_output_dir.glob("report.txt")) or list(complex_output_dir.glob("*.txt"))
            
            # Parse interaction summary from text report
            interactions = {}
            binding_site_count = 0
            if txt_files:
                interactions, binding_site_count = self._parse_plip_report(txt_files[0])
            
            plip_result = PLIPResult(
                complex_name=complex_name,
                pdb_file=pdb_file,
                success=len(png_files) > 0 or len(xml_files) > 0,
                png_diagram=png_files[0] if png_files else None,
                pymol_session=pse_files[0] if pse_files else None,
                xml_report=xml_files[0] if xml_files else None,
                txt_report=txt_files[0] if txt_files else None,
                binding_site_count=binding_site_count,
                interactions=interactions
            )
            
            if plip_result.success:
                logger.info(f"  ✅ Generated diagram for {complex_name}")
                if interactions:
                    interaction_summary = ", ".join(
                        f"{k}: {v}" for k, v in interactions.items() if v > 0
                    )
                    logger.info(f"     Interactions: {interaction_summary}")
            else:
                logger.warning(f"  ⚠️ No output generated for {complex_name}")
                plip_result.error_message = "No output files generated"
            
            return plip_result
            
        except subprocess.TimeoutExpired:
            logger.error(f"❌ PLIP timeout for {complex_name}")
            return PLIPResult(
                complex_name=complex_name,
                pdb_file=pdb_file,
                success=False,
                error_message="PLIP analysis timed out"
            )
        except Exception as e:
            logger.error(f"❌ Error analyzing {complex_name}: {e}")
            return PLIPResult(
                complex_name=complex_name,
                pdb_file=pdb_file,
                success=False,
                error_message=str(e)
            )
    
    def _analyze_with_python_api(
        self,
        pdb_file: Path,
        output_dir: Path,
        complex_name: str
    ) -> PLIPResult:
        """
        Analyze complex using PLIP Python API directly.
        
        This is used when the command-line tool is not available
        but the Python module is installed.
        """
        try:
            from plip.structure.preparation import PDBComplex
            from plip.exchange.report import BindingSiteReport
            import xml.etree.ElementTree as ET
            
            # Load and analyze the complex
            mol = PDBComplex()
            mol.load_pdb(str(pdb_file))
            mol.analyze()
            
            interactions = {
                'hydrogen_bonds': 0,
                'hydrophobic': 0,
                'pi_stacking': 0,
                'pi_cation': 0,
                'salt_bridges': 0,
                'halogen_bonds': 0,
                'water_bridges': 0,
                'metal_complexes': 0
            }
            
            binding_site_count = 0
            txt_report_lines = [f"PLIP Analysis Report for {complex_name}", "=" * 50, ""]
            
            # Iterate through binding sites
            for site_key, site in mol.interaction_sets.items():
                binding_site_count += 1
                txt_report_lines.append(f"Binding Site: {site_key}")
                txt_report_lines.append("-" * 30)
                
                # Count interactions
                interactions['hydrogen_bonds'] += len(site.hbonds_ldon) + len(site.hbonds_pdon)
                interactions['hydrophobic'] += len(site.hydrophobic_contacts)
                interactions['pi_stacking'] += len(site.pistacking)
                interactions['pi_cation'] += len(site.pication_laro) + len(site.pication_paro)
                interactions['salt_bridges'] += len(site.saltbridge_lneg) + len(site.saltbridge_pneg)
                interactions['halogen_bonds'] += len(site.halogen_bonds)
                interactions['water_bridges'] += len(site.water_bridges)
                interactions['metal_complexes'] += len(site.metal_complexes)
                
                # Add to report
                txt_report_lines.append(f"  Hydrogen Bonds: {len(site.hbonds_ldon) + len(site.hbonds_pdon)}")
                txt_report_lines.append(f"  Hydrophobic: {len(site.hydrophobic_contacts)}")
                txt_report_lines.append(f"  Pi-Stacking: {len(site.pistacking)}")
                txt_report_lines.append(f"  Salt Bridges: {len(site.saltbridge_lneg) + len(site.saltbridge_pneg)}")
                txt_report_lines.append("")
            
            # Save text report
            txt_report_file = output_dir / "report.txt"
            with open(txt_report_file, 'w') as f:
                f.write('\n'.join(txt_report_lines))
            
            # Note: PNG generation requires the command-line tool
            # XML generation could be added here if needed
            
            return PLIPResult(
                complex_name=complex_name,
                pdb_file=pdb_file,
                success=True,
                txt_report=txt_report_file,
                binding_site_count=binding_site_count,
                interactions=interactions
            )
            
        except Exception as e:
            logger.error(f"❌ Python API analysis failed for {complex_name}: {e}")
            return PLIPResult(
                complex_name=complex_name,
                pdb_file=pdb_file,
                success=False,
                error_message=str(e)
            )
    
    def _parse_plip_report(self, report_file: Path) -> tuple:
        """
        Parse PLIP text report to extract interaction counts.
        
        Returns
        -------
        tuple
            (interaction_counts dict, binding_site_count)
        """
        interactions = {
            'hydrogen_bonds': 0,
            'hydrophobic': 0,
            'pi_stacking': 0,
            'pi_cation': 0,
            'salt_bridges': 0,
            'halogen_bonds': 0,
            'water_bridges': 0,
            'metal_complexes': 0
        }
        binding_site_count = 0
        
        try:
            with open(report_file, 'r') as f:
                content = f.read()
            
            # Count binding sites
            binding_site_count = content.count('Interacting chain(s):')
            
            # Parse interaction counts from report
            import re
            
            patterns = {
                'hydrogen_bonds': r'Hydrogen Bonds:\s*(\d+)',
                'hydrophobic': r'Hydrophobic Interactions:\s*(\d+)',
                'pi_stacking': r'pi-Stacking:\s*(\d+)',
                'pi_cation': r'pi-Cation Interactions:\s*(\d+)',
                'salt_bridges': r'Salt Bridges:\s*(\d+)',
                'halogen_bonds': r'Halogen Bonds:\s*(\d+)',
                'water_bridges': r'Water Bridges:\s*(\d+)',
                'metal_complexes': r'Metal Complexes:\s*(\d+)'
            }
            
            for key, pattern in patterns.items():
                matches = re.findall(pattern, content, re.IGNORECASE)
                interactions[key] = sum(int(m) for m in matches)
            
        except Exception as e:
            logger.debug(f"Could not parse PLIP report: {e}")
        
        return interactions, binding_site_count
    
    def analyze_directory(
        self,
        poses_dir: Path,
        output_dir: Path,
        binding_category: Optional[str] = None,
        max_workers: int = 4
    ) -> Dict[str, PLIPResult]:
        """
        Analyze all PDB files in a directory.
        
        Parameters
        ----------
        poses_dir : Path
            Directory containing PDB files
        output_dir : Path
            Output directory for results
        binding_category : str, optional
            Binding category label (e.g., 'strong', 'moderate', 'weak')
        max_workers : int
            Maximum parallel workers for analysis
            
        Returns
        -------
        Dict[str, PLIPResult]
            Dictionary mapping complex names to their results
        """
        poses_dir = Path(poses_dir)
        output_dir = Path(output_dir)
        
        if binding_category:
            output_dir = output_dir / binding_category
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Find all PDB files
        pdb_files = list(poses_dir.glob("*.pdb"))
        
        if not pdb_files:
            logger.warning(f"⚠️ No PDB files found in {poses_dir}")
            return {}
        
        logger.info(f"📁 Processing {len(pdb_files)} PDB files from {poses_dir.name}")
        
        results = {}

        worker_count = max(1, min(_coerce_worker_count(max_workers, default=1), len(pdb_files)))
        if worker_count == 1:
            for pdb_file in pdb_files:
                result = self.analyze_complex(pdb_file, output_dir)
                results[result.complex_name] = result
        else:
            logger.info(f"⚙️ Using PLIP worker pool: {worker_count}")
            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                future_map = {
                    executor.submit(self.analyze_complex, pdb_file, output_dir): pdb_file
                    for pdb_file in pdb_files
                }
                for future in as_completed(future_map):
                    pdb_file = future_map[future]
                    try:
                        result = future.result()
                    except Exception as exc:
                        result = PLIPResult(
                            complex_name=pdb_file.stem,
                            pdb_file=pdb_file,
                            success=False,
                            error_message=str(exc),
                        )
                    results[result.complex_name] = result
        
        # Summary
        successful = sum(1 for r in results.values() if r.success)
        logger.info(f"✅ Completed: {successful}/{len(results)} complexes analyzed successfully")
        
        return results
    
    def analyze_best_poses_directory(
        self,
        best_poses_dir: Path,
        output_dir: Path,
        max_workers: int = 4,
    ) -> Dict[str, Dict[str, PLIPResult]]:
        """
        Analyze the complete best_poses directory structure.
        
        Expects structure:
        best_poses/
        ├── strong_binders/
        ├── moderate_binders/
        └── weak_binders/
        
        Parameters
        ----------
        best_poses_dir : Path
            Path to best_poses directory
        output_dir : Path
            Output directory for results
            
        Returns
        -------
        Dict[str, Dict[str, PLIPResult]]
            Nested dictionary: category -> complex_name -> result
        """
        best_poses_dir = Path(best_poses_dir)
        output_dir = Path(output_dir)
        
        # Expected subdirectories
        categories = ['strong_binders', 'moderate_binders', 'weak_binders']
        
        all_results = {}
        total_processed = 0
        total_successful = 0
        
        for category in categories:
            category_dir = best_poses_dir / category
            
            if category_dir.exists():
                logger.info(f"\n📂 Processing {category}...")
                results = self.analyze_directory(
                    category_dir,
                    output_dir,
                    binding_category=category,
                    max_workers=max_workers,
                )
                all_results[category] = results
                
                total_processed += len(results)
                total_successful += sum(1 for r in results.values() if r.success)
            else:
                logger.info(f"ℹ️ Skipping {category} (directory not found)")
        
        # Generate summary
        self._generate_summary(all_results, output_dir)
        
        logger.info(f"\n🎉 PLIP Analysis Complete!")
        logger.info(f"   Total complexes: {total_processed}")
        logger.info(f"   Successful: {total_successful}")
        logger.info(f"   Output: {output_dir}")
        
        return all_results
    
    def _generate_summary(
        self,
        all_results: Dict[str, Dict[str, PLIPResult]],
        output_dir: Path
    ):
        """Generate a JSON summary of all results."""
        summary = {
            'total_analyzed': 0,
            'total_successful': 0,
            'by_category': {},
            'interaction_totals': {
                'hydrogen_bonds': 0,
                'hydrophobic': 0,
                'pi_stacking': 0,
                'pi_cation': 0,
                'salt_bridges': 0,
                'halogen_bonds': 0,
                'water_bridges': 0,
                'metal_complexes': 0
            },
            'complexes': []
        }
        
        for category, results in all_results.items():
            category_summary = {
                'count': len(results),
                'successful': sum(1 for r in results.values() if r.success),
                'complexes': []
            }
            
            for name, result in results.items():
                summary['total_analyzed'] += 1
                if result.success:
                    summary['total_successful'] += 1
                
                # Aggregate interactions
                for itype, count in result.interactions.items():
                    summary['interaction_totals'][itype] += count
                
                complex_info = {
                    'name': name,
                    'category': category,
                    'success': result.success,
                    'binding_sites': result.binding_site_count,
                    'interactions': result.interactions,
                    'png_diagram': str(result.png_diagram) if result.png_diagram else None,
                    'error': result.error_message
                }
                
                category_summary['complexes'].append(complex_info)
                summary['complexes'].append(complex_info)
            
            summary['by_category'][category] = category_summary
        
        # Save summary
        summary_file = output_dir / 'plip_summary.json'
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"📊 Summary saved to {summary_file}")


def run_plip_analysis(
    poses_dir: Path,
    output_dir: Path,
    config: Dict = None
) -> Dict:
    """
    Main entry point for running PLIP analysis on docking results.
    
    Parameters
    ----------
    poses_dir : Path
        Path to best_poses directory or single poses directory
    output_dir : Path
        Output directory for results
    config : Dict, optional
        Configuration dictionary
        
    Returns
    -------
    Dict
        Analysis results summary
    """
    print("🔬 Running PLIP 2D Interaction Analysis...")
    
    poses_dir = Path(poses_dir)
    output_dir = Path(output_dir)

    if not poses_dir.exists():
        return {
            'success': False,
            'error': f'Poses directory not found: {poses_dir}',
        }
    if not poses_dir.is_dir():
        return {
            'success': False,
            'error': f'Poses path is not a directory: {poses_dir}',
        }
    
    # Initialize analyzer
    analyzer = PLIPAnalyzer(config)
    max_workers = _coerce_worker_count((config or {}).get("max_workers", 4), default=1)
    
    if not analyzer.plip_available:
        print("❌ PLIP is not installed. Please install with:")
        print("   pip install plip")
        print("   or: conda install -c conda-forge plip")
        return {'success': False, 'error': 'PLIP not installed'}
    
    # Check if this is a best_poses directory structure
    subdirs = ['strong_binders', 'moderate_binders', 'weak_binders']
    is_best_poses_structure = any((poses_dir / d).exists() for d in subdirs)
    
    if is_best_poses_structure:
        # Analyze the complete structure
        results = analyzer.analyze_best_poses_directory(poses_dir, output_dir, max_workers=max_workers)
        
        # Flatten for return
        total = sum(len(r) for r in results.values())
        successful = sum(
            sum(1 for res in r.values() if res.success)
            for r in results.values()
        )
        if total == 0:
            return {
                'success': False,
                'error': f'No PDB files found under categorized poses directory: {poses_dir}',
                'output_dir': str(output_dir),
            }
        
        return {
            'success': True,
            'total_analyzed': total,
            'successful': successful,
            'by_category': {
                cat: len(res) for cat, res in results.items()
            },
            'output_dir': str(output_dir)
        }
    else:
        # Single directory
        results = analyzer.analyze_directory(poses_dir, output_dir, max_workers=max_workers)
        if not results:
            return {
                'success': False,
                'error': f'No PDB files found in poses directory: {poses_dir}',
                'output_dir': str(output_dir),
            }
        
        return {
            'success': True,
            'total_analyzed': len(results),
            'successful': sum(1 for r in results.values() if r.success),
            'output_dir': str(output_dir)
        }


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python plip_integration.py <poses_dir> [output_dir]")
        print("\nExample:")
        print("  python plip_integration.py ./best_poses ./plip_output")
        sys.exit(1)
    
    poses_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("./plip_output")
    
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s'
    )
    
    results = run_plip_analysis(poses_dir, output_dir)
    print(f"\nResults: {json.dumps(results, indent=2)}")
