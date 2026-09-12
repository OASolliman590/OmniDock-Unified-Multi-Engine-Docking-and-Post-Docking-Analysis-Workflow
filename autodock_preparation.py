#!/usr/bin/env python3
"""
Enhanced AutoDock Preparation Module
Integrates with PDB Prepare Wizard for comprehensive molecular docking preparation
"""

import os
import json
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
import shutil
import tempfile

from docking.models import (
    normalize_ligand_preparation_profile,
    resolve_effective_ligand_preparation_profile,
    validate_ligand_preparation_profile,
)
from docking.preparation.ligand_quality import validate_prepared_ligand_pdbqt
from docking.preparation.ligand_preparation import prepare_ligand_for_vina_family


def profile_requires_autodocktools(profile: str, selected_engines: Optional[List[str]] = None) -> bool:
    effective = resolve_effective_ligand_preparation_profile(profile, selected_engines)
    if effective in {"autodocktools_only", "openbabel_autodocktools", "openbabel_meeko_autodock"}:
        return True
    return False

@dataclass
class PreparationConfig:
    """Configuration for AutoDock preparation"""
    ligands_input: str
    receptors_input: str
    ligands_output: str
    receptors_output: str
    force_field: str = "AMBER"
    ph: float = 7.4
    allow_bad_res: bool = False
    default_altloc: str = "A"
    receptor_use_pdb2pqr: bool = False
    validate_outputs: bool = True
    min_file_size_kb: int = 1
    ligand_preparation_backend: str = "engine_aware_full"  # backward-compatible field
    ligand_preparation_profile: str = "engine_aware_full"
    selected_engines: List[str] = None
    autodocktools_prepare_ligand4: str = ""
    autodocktools_prepare_receptor4: str = ""
    autodocktools_python: str = ""

class AutoDockPreparationPipeline:
    """
    Enhanced AutoDock preparation pipeline.
    """
    
    def __init__(self, config: Optional[PreparationConfig] = None):
        self.config = config or PreparationConfig(
            ligands_input="./ligands_raw",
            receptors_input="./receptors_raw", 
            ligands_output="./ligands_prep",
            receptors_output="./receptors_prep"
        )
        if self.config.selected_engines is None:
            self.config.selected_engines = []
        # Keep legacy backend field in sync with the new profile field.
        profile = str(self.config.ligand_preparation_profile or self.config.ligand_preparation_backend or "engine_aware_full").strip()
        self.config.ligand_preparation_profile = profile
        self.config.ligand_preparation_backend = profile
        self.logger = self._setup_logging()
        self.script_dir = Path(__file__).parent
        self.enhanced_script = self.script_dir / "prep_autodock_enhanced.sh"

    def _tool_usable(self, tool: str) -> bool:
        if not shutil.which(tool):
            return False
        if tool in {"mk_prepare_ligand.py", "mk_prepare_receptor.py"}:
            completed = subprocess.run(
                [tool, "--help"],
                capture_output=True,
                text=True,
            )
            return completed.returncode == 0
        return True

    def _resolve_autodocktools_prepare_ligand4(self) -> Optional[Path]:
        candidates: List[Path] = []

        script_hint = (self.config.autodocktools_prepare_ligand4 or "").strip()
        if script_hint:
            candidates.append(Path(script_hint).expanduser())

        for env_key in ("AUTODOCKTOOLS_PREPARE_LIGAND4", "ADT_PREPARE_LIGAND4"):
            env_value = str(os.environ.get(env_key, "") or "").strip()
            if env_value:
                candidates.append(Path(env_value).expanduser())

        for command_name in ("prepare_ligand4.py", "prepare_ligand4"):
            resolved = shutil.which(command_name)
            if resolved:
                candidates.append(Path(resolved).expanduser())

        mgltools_root = str(os.environ.get("MGLTOOLS_PATH", "") or "").strip()
        if mgltools_root:
            candidates.append(Path(mgltools_root).expanduser() / "MGLToolsPckgs" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py")

        candidates.extend(
            [
                self.script_dir / ".workflow" / "tools" / "autodocktools-prepare-py3k" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py",
                self.script_dir / "tools" / "autodocktools-prepare-py3k" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py",
                Path.home() / "mgltools_x86_64Linux2_1.5.7" / "MGLToolsPckgs" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py",
                Path("/opt/mgltools/1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"),
                Path("/usr/local/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"),
            ]
        )

        for candidate in candidates:
            if candidate.exists() and candidate.is_file():
                return candidate.resolve()
        return None

    @staticmethod
    def _derive_receptor_script_from_ligand_script(ligand_script: Optional[Path]) -> Optional[Path]:
        if ligand_script is None:
            return None
        ligand_path = Path(ligand_script).expanduser()
        receptor_candidate = ligand_path.with_name("prepare_receptor4.py")
        if receptor_candidate.exists() and receptor_candidate.is_file():
            return receptor_candidate.resolve()
        return None

    def _resolve_autodocktools_prepare_receptor4(self, ligand_script: Optional[Path] = None) -> Optional[Path]:
        candidates: List[Path] = []

        script_hint = (self.config.autodocktools_prepare_receptor4 or "").strip()
        if script_hint:
            candidates.append(Path(script_hint).expanduser())

        for env_key in ("AUTODOCKTOOLS_PREPARE_RECEPTOR4", "ADT_PREPARE_RECEPTOR4"):
            env_value = str(os.environ.get(env_key, "") or "").strip()
            if env_value:
                candidates.append(Path(env_value).expanduser())

        if ligand_script is not None:
            derived = self._derive_receptor_script_from_ligand_script(ligand_script)
            if derived is not None:
                candidates.append(derived)

        configured_ligand_hint = (self.config.autodocktools_prepare_ligand4 or "").strip()
        if configured_ligand_hint:
            derived = self._derive_receptor_script_from_ligand_script(Path(configured_ligand_hint).expanduser())
            if derived is not None:
                candidates.append(derived)

        for env_key in ("AUTODOCKTOOLS_PREPARE_LIGAND4", "ADT_PREPARE_LIGAND4"):
            env_value = str(os.environ.get(env_key, "") or "").strip()
            if env_value:
                derived = self._derive_receptor_script_from_ligand_script(Path(env_value).expanduser())
                if derived is not None:
                    candidates.append(derived)

        for command_name in ("prepare_receptor4.py", "prepare_receptor4"):
            resolved = shutil.which(command_name)
            if resolved:
                candidates.append(Path(resolved).expanduser())

        mgltools_root = str(os.environ.get("MGLTOOLS_PATH", "") or "").strip()
        if mgltools_root:
            candidates.append(Path(mgltools_root).expanduser() / "MGLToolsPckgs" / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py")

        candidates.extend(
            [
                self.script_dir / ".workflow" / "tools" / "autodocktools-prepare-py3k" / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py",
                self.script_dir / "tools" / "autodocktools-prepare-py3k" / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py",
                Path.home() / "mgltools_x86_64Linux2_1.5.7" / "MGLToolsPckgs" / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py",
                Path("/opt/mgltools/1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"),
                Path("/usr/local/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"),
            ]
        )

        for candidate in candidates:
            if candidate.exists() and candidate.is_file():
                return candidate.resolve()
        return None
        
    def _setup_logging(self) -> logging.Logger:
        """Setup logging for the preparation pipeline"""
        logger = logging.getLogger('autodock_preparation')
        logger.setLevel(logging.INFO)
        
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            
        return logger
    
    def create_config_file(self, config_path: str = "autodock_config.json") -> str:
        """
        Create a configuration file for the enhanced bash script
        
        Args:
            config_path: Path to save the configuration file
            
        Returns:
            Path to the created configuration file
        """
        config = {
            "input": {
                "ligands": {
                    "path": self.config.ligands_input,
                    "formats": ["sdf", "mol", "mol2", "pdb"],
                    "in_same_folder": False
                },
                "receptors": {
                    "path": self.config.receptors_input,
                    "formats": ["pdb"],
                    "in_same_folder": False
                }
            },
            "output": {
                "ligands": self.config.ligands_output,
                "receptors": self.config.receptors_output,
                "logs": "./logs"
            },
            "preparation": {
                "force_field": self.config.force_field,
                "ph": self.config.ph,
                "allow_bad_res": self.config.allow_bad_res,
                "default_altloc": self.config.default_altloc,
                "receptor_use_pdb2pqr": bool(self.config.receptor_use_pdb2pqr),
                "ligand_preparation_backend": normalize_ligand_preparation_profile(self.config.ligand_preparation_profile or self.config.ligand_preparation_backend),
                "ligand_preparation_profile": normalize_ligand_preparation_profile(self.config.ligand_preparation_profile or self.config.ligand_preparation_backend),
                "selected_engines": [str(engine).strip().lower() for engine in (self.config.selected_engines or []) if str(engine).strip()],
                "autodocktools_prepare_ligand4": self.config.autodocktools_prepare_ligand4,
                "autodocktools_prepare_receptor4": self.config.autodocktools_prepare_receptor4,
                "autodocktools_python": self.config.autodocktools_python,
            },
            "quality_control": {
                "validate_outputs": self.config.validate_outputs,
                "check_file_sizes": True,
                "min_file_size_kb": self.config.min_file_size_kb
            }
        }
        
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=2)
            
        self.logger.info(f"Configuration file created: {config_path}")
        return config_path
    
    def check_dependencies(self) -> Tuple[bool, List[str]]:
        """
        Check if all required dependencies are available
        
        Returns:
            Tuple of (all_available, missing_dependencies)
        """
        required_tools = ["obabel", "jq"]
        optional_tools = ["pdb2pqr30", "mk_prepare_ligand.py", "mk_prepare_receptor.py"]

        missing = []
        for tool in required_tools:
            if not shutil.which(tool):
                missing.append(tool)

        if missing:
            self.logger.error(f"Missing dependencies: {missing}")
            return False, missing

        missing_optional = [tool for tool in optional_tools if not self._tool_usable(tool)]
        if missing_optional:
            self.logger.warning(
                "Preparation tools missing; stages requiring these tools will fail with diagnostics: %s",
                ", ".join(missing_optional),
            )

        profile = normalize_ligand_preparation_profile(
            self.config.ligand_preparation_profile or self.config.ligand_preparation_backend
        )
        self.config.ligand_preparation_profile = profile
        self.config.ligand_preparation_backend = profile
        compatibility = validate_ligand_preparation_profile(profile, self.config.selected_engines)
        for warning in compatibility.warnings:
            self.logger.warning(warning)
        if not compatibility.is_valid:
            missing.extend(compatibility.errors)
            self.logger.error(
                "Invalid ligand preparation profile/engine combination: %s",
                "; ".join(compatibility.errors),
            )
            return False, missing

        if profile_requires_autodocktools(profile, self.config.selected_engines):
            resolved_script = self._resolve_autodocktools_prepare_ligand4()
            if resolved_script is None:
                script_hint = (self.config.autodocktools_prepare_ligand4 or "").strip()
                if script_hint:
                    missing.append(f"prepare_ligand4.py@{script_hint}")
                else:
                    missing.append("prepare_ligand4.py")
            else:
                self.config.autodocktools_prepare_ligand4 = str(resolved_script)
                self.logger.info("AutoDockTools script resolved: %s", resolved_script)
            if missing:
                self.logger.error("Ligand preparation profile '%s' requires AutoDockTools but dependencies are missing: %s", profile, missing)
                return False, missing

        self.logger.info("All required dependencies found")
        return True, []
    
    def prepare_ligands_from_pdb(self, pdb_file: str, output_dir: str) -> Optional[str]:
        """
        Extract ligands from PDB file and prepare them for AutoDock
        
        Args:
            pdb_file: Path to PDB file
            output_dir: Output directory for prepared ligands
            
        Returns:
            Path to prepared ligand file or None if failed
        """
        try:
            from core_pipeline import MolecularDockingPipeline
            
            # Initialize pipeline
            pipeline = MolecularDockingPipeline(output_dir=output_dir)
            
            # Enumerate HETATMs
            hetatm_details, unique_hetatms = pipeline.enumerate_hetatms(pdb_file)
            
            if not unique_hetatms:
                self.logger.warning(f"No ligands found in {pdb_file}")
                return None
                
            # Get the first ligand (or let user choose)
            ligand_name = unique_hetatms[0]
            ligand_info = next(({"chain": chain, "res_id": resid} for name, chain, resid, _ in hetatm_details if name == ligand_name), None)
            if ligand_info is None:
                raise ValueError("Selected ligand instance was not found")
            
            # Extract PDB ID from filename (e.g., "1ABC.pdb" -> "1ABC")
            pdb_id = Path(pdb_file).stem.upper()
            if len(pdb_id) == 4 and pdb_id.isalnum():
                pdb_id = pdb_id
            else:
                pdb_id = None  # Fallback if can't extract
            
            # Extract ligand
            ligand_pdb = pipeline.save_hetatm_as_pdb(
                pdb_file, 
                ligand_name, 
                ligand_info['chain'], 
                ligand_info['res_id'],
                pdb_id=pdb_id
            )
            
            if ligand_pdb:
                base_name = Path(ligand_pdb).stem
                pdbqt_file = Path(output_dir) / f"{base_name}.pdbqt"

                # Meeko expects explicit hydrogens and 3D coordinates, so
                # normalize the extracted ligand before PDBQT conversion.
                prepare_ligand_for_vina_family(Path(ligand_pdb), pdbqt_file)
                self.logger.info(f"Prepared ligand: {pdbqt_file}")
                return str(pdbqt_file)
                
        except Exception as e:
            self.logger.error(f"Failed to prepare ligand from PDB: {e}")
            return None
    
    def run_enhanced_preparation(self, config_path: Optional[str] = None) -> bool:
        """
        Run the enhanced AutoDock preparation script
        
        Args:
            config_path: Path to configuration file
            
        Returns:
            True if successful, False otherwise
        """
        if not config_path:
            config_path = self.create_config_file()
            
        if not self.enhanced_script.exists():
            self.logger.error(f"Enhanced script not found: {self.enhanced_script}")
            return False
            
        try:
            # Make script executable
            os.chmod(self.enhanced_script, 0o755)

            env = dict(os.environ)
            profile = normalize_ligand_preparation_profile(
                self.config.ligand_preparation_profile or self.config.ligand_preparation_backend
            )
            env["PDBWIZARD_LIGAND_PREP_BACKEND"] = profile
            env["PDBWIZARD_LIGAND_PREP_PROFILE"] = profile
            env["PDBWIZARD_LIGAND_PREP_PH"] = str(self.config.ph)
            if self.config.selected_engines:
                env["PDBWIZARD_SELECTED_ENGINES"] = ",".join(
                    str(engine).strip().lower()
                    for engine in self.config.selected_engines
                    if str(engine).strip()
                )
            if self.config.autodocktools_prepare_ligand4:
                env["AUTODOCKTOOLS_PREPARE_LIGAND4"] = str(Path(self.config.autodocktools_prepare_ligand4).expanduser())
            if self.config.autodocktools_prepare_receptor4:
                env["AUTODOCKTOOLS_PREPARE_RECEPTOR4"] = str(Path(self.config.autodocktools_prepare_receptor4).expanduser())
            if self.config.autodocktools_python:
                env["AUTODOCKTOOLS_PYTHON"] = str(Path(self.config.autodocktools_python).expanduser())

            # Run the enhanced script and stream logs in real-time so long
            # preparations do not appear stuck in interactive mode.
            command = [str(self.enhanced_script), config_path]
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                env=env,
                bufsize=1,
            )
            streamed_lines: List[str] = []
            if process.stdout is not None:
                for line in process.stdout:
                    clean_line = line.rstrip()
                    streamed_lines.append(clean_line)
                    if clean_line:
                        self.logger.info(clean_line)
            return_code = process.wait()

            if return_code == 0:
                invalid_ligands = []
                ligands_output = Path(self.config.ligands_output).expanduser().resolve()
                if ligands_output.exists():
                    for ligand_file in sorted(ligands_output.glob("*.pdbqt")):
                        issues = validate_prepared_ligand_pdbqt(ligand_file)
                        if issues:
                            invalid_ligands.append(
                                {
                                    "ligand": ligand_file.name,
                                    "issues": [issue.to_dict() for issue in issues],
                                }
                            )
                if invalid_ligands:
                    report_path = ligands_output / "ligand_preparation_validation_report.json"
                    with open(report_path, "w", encoding="utf-8") as handle:
                        json.dump({"ligands_output": str(ligands_output), "invalid_ligands": invalid_ligands}, handle, indent=2)
                    self.logger.error(
                        "Prepared ligand validation failed. Invalid AutoDock-ready outputs were generated: %s",
                        ", ".join(item["ligand"] for item in invalid_ligands),
                    )
                    self.logger.error("Validation report saved to: %s", report_path)
                    return False
                self.logger.info("AutoDock preparation completed successfully")
                self.logger.info("Ligands output directory: %s", ligands_output)
                self.logger.info("Receptors output directory: %s", Path(self.config.receptors_output).expanduser().resolve())
                return True
            else:
                self.logger.error("Preparation failed with exit code %s", return_code)
                if streamed_lines:
                    self.logger.error("Last preparation log line: %s", streamed_lines[-1])
                return False
                
        except Exception as e:
            self.logger.error(f"Failed to run preparation script: {e}")
            return False
    
    def analyze_preparation_results(
        self,
        ligands_output_dir: Optional[str] = None,
        receptors_output_dir: Optional[str] = None,
    ) -> Dict:
        """
        Analyze the results of AutoDock preparation
        
        Args:
            ligands_output_dir: Directory containing prepared ligand PDBQT files
            receptors_output_dir: Directory containing prepared receptor PDBQT files
            
        Returns:
            Dictionary with analysis results
        """
        ligands_path = Path(ligands_output_dir or self.config.ligands_output).expanduser().resolve()
        receptors_path = Path(receptors_output_dir or self.config.receptors_output).expanduser().resolve()
        
        results = {
            "ligands": {
                "count": 0,
                "files": [],
                "total_size_mb": 0
            },
            "receptors": {
                "count": 0,
                "files": [],
                "total_size_mb": 0
            },
            "paths": {
                "ligands_output": str(ligands_path),
                "receptors_output": str(receptors_path),
            }
        }
        
        # Analyze ligands from ligand output directory
        ligand_files = list(ligands_path.glob("*.pdbqt")) if ligands_path.exists() else []
        results["ligands"]["count"] = len(ligand_files)
        results["ligands"]["files"] = [str(f) for f in ligand_files]
        results["ligands"]["total_size_mb"] = sum(f.stat().st_size for f in ligand_files) / (1024 * 1024)
        
        # Analyze receptors from receptor output directory
        receptor_files = list(receptors_path.glob("*.pdbqt")) if receptors_path.exists() else []
        results["receptors"]["count"] = len(receptor_files)
        results["receptors"]["files"] = [str(f) for f in receptor_files]
        results["receptors"]["total_size_mb"] = sum(f.stat().st_size for f in receptor_files) / (1024 * 1024)
        
        return results
    
    def generate_preparation_report(self, results: Dict, output_file: str = "preparation_report.json"):
        """
        Generate a comprehensive report of the preparation process
        
        Args:
            results: Results from analyze_preparation_results
            output_file: Path to save the report
        """
        report = {
            "preparation_summary": {
                "timestamp": str(Path().cwd()),
                "config": {
                    "force_field": self.config.force_field,
                    "ph": self.config.ph,
                    "ligand_preparation_profile": normalize_ligand_preparation_profile(
                        self.config.ligand_preparation_profile or self.config.ligand_preparation_backend
                    ),
                },
                "results": results
            },
            "recommendations": self._generate_recommendations(results)
        }
        
        with open(output_file, 'w') as f:
            json.dump(report, f, indent=2)
            
        self.logger.info(f"Preparation report saved: {output_file}")
    
    def _generate_recommendations(self, results: Dict) -> List[str]:
        """Generate recommendations based on preparation results"""
        recommendations = []
        
        if results["ligands"]["count"] == 0:
            recommendations.append("No ligands were prepared. Check input directory and file formats.")
        
        if results["receptors"]["count"] == 0:
            recommendations.append("No receptors were prepared. Check input directory and file formats.")
        
        if results["ligands"]["count"] > 0 and results["receptors"]["count"] > 0:
            recommendations.append("Ready for AutoDock Vina docking. Use the prepared PDBQT files.")
        
        return recommendations

def main():
    """Main function for command-line usage"""
    import argparse
    import tempfile
    
    parser = argparse.ArgumentParser(description="Enhanced AutoDock Preparation Pipeline")
    parser.add_argument("--config", help="Configuration file path")
    parser.add_argument("--ligands-input", help="Input directory or file for ligands")
    parser.add_argument("--receptors-input", help="Input directory or file for receptors")
    parser.add_argument("--ligands-output", help="Output directory for prepared ligands")
    parser.add_argument("--receptors-output", help="Output directory for prepared receptors")
    parser.add_argument("--force-field", default="AMBER", help="Force field for PDB2PQR")
    parser.add_argument("--ph", type=float, default=7.4, help="pH for protonation")
    parser.add_argument(
        "--ligand-profile",
        default="engine_aware_full",
        help=(
            "Ligand preparation profile: "
            "openbabel_only, meeko_only, autodocktools_only, openbabel_meeko, "
            "openbabel_meeko_autodock, openbabel_autodocktools, engine_aware_full"
        ),
    )
    parser.add_argument("--selected-engines", default="", help="Comma-separated engines for engine-aware profile")
    parser.add_argument("--create-config", action="store_true", help="Create configuration file and exit")
    
    args = parser.parse_args()

    def handle_input_path(input_path):
        if input_path and Path(input_path).is_file():
            temp_dir = tempfile.mkdtemp()
            shutil.copy(input_path, temp_dir)
            return temp_dir
        return input_path

    ligands_input = handle_input_path(args.ligands_input)
    receptors_input = handle_input_path(args.receptors_input)

    # Create configuration
    config = PreparationConfig(
        ligands_input=ligands_input or "./ligands_raw",
        receptors_input=receptors_input or "./receptors_raw",
        ligands_output=args.ligands_output or "./ligands_prep", 
        receptors_output=args.receptors_output or "./receptors_prep",
        force_field=args.force_field,
        ph=args.ph,
        ligand_preparation_profile=normalize_ligand_preparation_profile(args.ligand_profile),
        ligand_preparation_backend=normalize_ligand_preparation_profile(args.ligand_profile),
        selected_engines=[token.strip().lower() for token in str(args.selected_engines or "").split(",") if token.strip()],
    )
    
    # Initialize pipeline
    pipeline = AutoDockPreparationPipeline(config)
    
    if args.create_config:
        config_path = pipeline.create_config_file()
        print(f"Configuration file created: {config_path}")
        return
    
    # Check dependencies
    deps_ok, missing = pipeline.check_dependencies()
    if not deps_ok:
        print(f"Missing dependencies: {missing}")
        return 1
    
    # Run preparation
    success = pipeline.run_enhanced_preparation(args.config)
    
    if success:
        # Analyze results
        results = pipeline.analyze_preparation_results(
            ligands_output_dir=config.ligands_output,
            receptors_output_dir=config.receptors_output,
        )
        pipeline.generate_preparation_report(results)
        
        print("Preparation completed successfully!")
        print(f"Ligands prepared: {results['ligands']['count']}")
        print(f"Receptors prepared: {results['receptors']['count']}")
        
        return 0
    else:
        print("Preparation failed. Check logs for details.")
        return 1

if __name__ == "__main__":
    exit(main())
