#!/usr/bin/env python3
"""
PDB Prepare Wizard - Core Pipeline Module
========================================

Consolidated core pipeline for preparing PDB files for molecular docking studies.
This module provides functionality to:
- Download PDB files from RCSB PDB database
- Extract and analyze ligands (HETATMs)
- Clean PDB files by removing unwanted residues
- Analyze pocket properties and druggability
- Generate comprehensive reports

Author: Molecular Docking Pipeline
Version: 2.1.0
"""

import os
import re
import sys
import subprocess
import urllib.request
import copy
import shutil
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Any
import warnings
warnings.filterwarnings("ignore")

# Add Excel support
try:
    import openpyxl
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment
    EXCEL_AVAILABLE = True
except ImportError:
    EXCEL_AVAILABLE = False

try:
    from Bio.PDB import PDBList, PDBParser, MMCIFParser, Select, PDBIO
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain
except ImportError as e:
    print(f"❌ Error importing Biopython: {e}")
    print("Please install Biopython: pip install biopython")
    sys.exit(1)

class MolecularDockingPipeline:
    """
    A comprehensive molecular docking pipeline for PDB file preparation and analysis.
    
    This class provides methods to:
    - Download and process PDB files
    - Extract ligands and analyze binding sites
    - Clean and prepare structures for docking
    - Analyze pocket properties and druggability
    """
    
    def __init__(self, output_dir: str = "pipeline_output"):
        """
        Initialize the molecular docking pipeline.
        
        Args:
            output_dir (str): Directory to store all output files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        print(f"✓ Pipeline initialized. Output directory: {self.output_dir}")
        
        # Validate required packages
        self._validate_dependencies()
        
    def _validate_dependencies(self):
        """Validate that all required dependencies are available."""
        try:
            import numpy as np
            import pandas as pd
            from Bio.PDB import PDBParser
            print("✓ All core dependencies available")
        except ImportError as e:
            print(f"❌ Missing dependency: {e}")
            sys.exit(1)
            
        # Check for optional PLIP
        try:
            import plip
            print("✓ PLIP available for advanced interaction analysis")
            self.plip_available = True
        except ImportError:
            print("⚠️  PLIP not available - will use distance-based analysis")
            self.plip_available = False

    @staticmethod
    def _cleaning_counts(path: str) -> Dict[str, int]:
        counts = {"ATOM": 0, "HETATM": 0, "TER": 0}
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    rec = line[:6].strip()
                    if rec in counts:
                        counts[rec] += 1
        except Exception:
            return counts
        return counts

    @staticmethod
    def _find_pdb_tools() -> Optional[Dict[str, List[str]]]:
        required = ("pdb_selchain", "pdb_delresname", "pdb_tidy")
        resolved: Dict[str, List[str]] = {}
        all_on_path = True
        for tool_name in required:
            tool_path = shutil.which(tool_name)
            if not tool_path:
                all_on_path = False
                break
            resolved[tool_name] = [tool_path]
        if all_on_path:
            return resolved

        # Fallback to a bundled checkout under .workflow/tools/pdb-tools
        repo_root = Path(__file__).resolve().parent
        bundle_root = repo_root / ".workflow" / "tools" / "pdb-tools" / "pdbtools"
        if bundle_root.is_dir():
            bundled: Dict[str, List[str]] = {}
            for tool_name in required:
                script_path = bundle_root / f"{tool_name}.py"
                if not script_path.is_file():
                    return None
                bundled[tool_name] = [sys.executable, str(script_path)]
            return bundled

        return None

    @staticmethod
    def _run_pdb_tool(tool_cmd: List[str], args: List[str], input_text: str) -> str:
        cmd = [*tool_cmd, *args]
        completed = subprocess.run(
            cmd,
            input=input_text,
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0:
            stderr = (completed.stderr or "").strip()
            raise RuntimeError(
                f"{Path(tool_cmd[-1]).name} failed with code {completed.returncode}: {stderr or 'no stderr output'}"
            )
        return completed.stdout

    def _clean_pdb_with_pdb_tools(
        self,
        pdb_file: str,
        output_filename: str,
        to_remove_list: List[str],
        keep_chain_id: Optional[str],
        keep_chain_ids: Optional[List[str]] = None,
        tool_paths: Optional[Dict[str, List[str]]] = None,
    ) -> bool:
        """
        Try stream-based cleaning using pdb-tools.
        Returns True on success, False when pdb-tools are unavailable.
        """
        tool_paths = tool_paths or self._find_pdb_tools()
        if not tool_paths:
            return False

        with open(pdb_file, "r", encoding="utf-8", errors="replace") as handle:
            pdb_text = handle.read()

        keep_chains: List[str] = []
        if keep_chain_ids:
            keep_chains.extend([str(chain).strip() for chain in keep_chain_ids if str(chain).strip()])
        keep_chain = str(keep_chain_id or "").strip()
        if keep_chain and keep_chain not in keep_chains:
            keep_chains.append(keep_chain)
        if keep_chains:
            pdb_text = self._run_pdb_tool(tool_paths["pdb_selchain"], [f"-{','.join(keep_chains)}"], pdb_text)

        removable = sorted(
            {
                str(resname).strip().upper()
                for resname in to_remove_list
                if str(resname).strip() and len(str(resname).strip()) <= 3
            }
        )
        if removable:
            pdb_text = self._run_pdb_tool(
                tool_paths["pdb_delresname"],
                [f"-{','.join(removable)}"],
                pdb_text,
            )

        pdb_text = self._run_pdb_tool(tool_paths["pdb_tidy"], [], pdb_text)

        with open(output_filename, "w", encoding="utf-8", errors="replace") as handle:
            handle.write(pdb_text)

        return True

    @staticmethod
    def _load_structure(structure_id: str, structure_file: str):
        structure_path = Path(structure_file)
        suffix = structure_path.suffix.lower()
        if suffix in {".cif", ".mmcif"}:
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure(structure_id, str(structure_path))
        else:
            from io import StringIO
            from docking.preparation.structure_contract import select_pdb_lines
            parser = PDBParser(QUIET=True)
            lines = structure_path.read_text(encoding="utf-8").splitlines()
            if any(line.startswith(("ATOM  ", "HETATM")) for line in lines):
                lines, _ = select_pdb_lines(lines)
            structure = parser.get_structure(structure_id, StringIO("\n".join(lines) + "\n"))
        if len(list(structure)) > 1:
            raise ValueError("Multiple models: select one explicit model before preparation")
        return structure

    @staticmethod
    def _normalize_structure_for_pdbio(structure):
        normalized = copy.deepcopy(structure)
        valid_ids = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")

        for model in normalized:
            used_ids = set()
            chain_map: Dict[str, str] = {}
            for chain in model:
                original_id = str(chain.get_id()).strip() or "A"
                chain.xtra["original_chain_id"] = original_id
                candidate = original_id[0]
                if len(original_id) == 1 and candidate not in used_ids:
                    chain.id = candidate
                    chain.xtra["normalized_chain_id"] = candidate
                    used_ids.add(candidate)
                    chain_map[original_id] = candidate
                    continue

                if original_id in chain_map:
                    chain.id = chain_map[original_id]
                    chain.xtra["normalized_chain_id"] = chain_map[original_id]
                    continue

                for fallback_id in valid_ids:
                    if fallback_id not in used_ids:
                        chain.id = fallback_id
                        chain.xtra["normalized_chain_id"] = fallback_id
                        used_ids.add(fallback_id)
                        chain_map[original_id] = fallback_id
                        break

        return normalized
    
    def fetch_pdb(self, pdb_id: str) -> str:
        """
        Download PDB file from RCSB PDB database.
        
        Args:
            pdb_id (str): PDB identifier (e.g., '1ABC')
            
        Returns:
            str: Path to downloaded PDB file
        """
        print(f"🔄 Fetching PDB {pdb_id}...")
        try:
            pdbl = PDBList()
            filename = pdbl.retrieve_pdb_file(
                pdb_id.lower(), 
                pdir=str(self.output_dir), 
                file_format='pdb'
            )
            new_filename = self.output_dir / f"{pdb_id.upper()}.pdb"
            if filename:
                os.rename(filename, new_filename)
            else:
                raise FileNotFoundError(f"RCSB download returned no file for PDB {pdb_id.upper()}")
            print(f"✓ Downloaded: {new_filename}")
            return str(new_filename)
        except Exception as e:
            print(f"⚠️  Primary Biopython download failed for PDB {pdb_id}: {e}")
            direct_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            new_filename = self.output_dir / f"{pdb_id.upper()}.pdb"
            try:
                urllib.request.urlretrieve(direct_url, new_filename)
                print(f"✓ Downloaded via direct RCSB fallback: {new_filename}")
                return str(new_filename)
            except Exception as pdb_fallback_error:
                cif_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
                cif_filename = self.output_dir / f"{pdb_id.upper()}.cif"
                try:
                    urllib.request.urlretrieve(cif_url, cif_filename)
                    subprocess.run(
                        ["obabel", str(cif_filename), "-O", str(new_filename)],
                        check=True,
                        capture_output=True,
                        text=True,
                    )
                    print(f"✓ Downloaded via mmCIF fallback and converted to PDB: {new_filename}")
                    return str(new_filename)
                except Exception as cif_fallback_error:
                    print(f"❌ Failed to download PDB {pdb_id}: {pdb_fallback_error}; mmCIF fallback also failed: {cif_fallback_error}")
                    raise

    def enumerate_hetatms(self, pdb_file: str) -> Tuple[List[Tuple], List[str]]:
        """
        List all HETATM residues in the PDB file for user selection.
        
        Args:
            pdb_file (str): Path to PDB file
            
        Returns:
            Tuple[List[Tuple], List[str]]: (detailed hetatm info, unique hetatm types)
        """
        print(f"🔄 Enumerating HETATMs in {pdb_file}...")
        try:
            structure = self._load_structure('protein', pdb_file)
            
            hetatm_details = []
            hetatm_counts = {}
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.get_id()[0] != ' ':  # HETATM check
                            resname = residue.get_resname()
                            chain_id = chain.get_id()
                            res_id = residue.get_id()[1]
                            
                            hetatm_details.append((resname, chain_id, res_id, residue))
                            hetatm_counts[resname] = hetatm_counts.get(resname, 0) + 1
            
            unique_hetatms = list(hetatm_counts.keys())
            
            if not unique_hetatms:
                print("⚠️  No HETATMs found in this structure")
                return [], []
            
            print(f"✓ Found {len(hetatm_details)} HETATM instances of {len(unique_hetatms)} types:")
            for resname, count in hetatm_counts.items():
                print(f"   {resname}: {count} instance(s)")
            
            return hetatm_details, unique_hetatms
        except Exception as e:
            print(f"❌ Error enumerating HETATMs: {e}")
            return [], []

    def list_protein_chains(self, pdb_file: str) -> List[str]:
        """
        Return sorted chain IDs that contain ATOM records from the source structure.
        These are source PDB/mmCIF chain IDs (macromolecular chains), not inferred from
        ligand interactions.
        """
        chains: set[str] = set()
        try:
            structure = self._load_structure("protein", pdb_file)
            for model in structure:
                for chain in model:
                    has_atom = False
                    for residue in chain:
                        if residue.get_id()[0] == " ":
                            has_atom = True
                            break
                    if has_atom:
                        chain_id = str(chain.get_id()).strip()
                        if chain_id:
                            chains.add(chain_id)
            return sorted(chains)
        except Exception:
            return []

    def parse_remark350_chain_groups(self, pdb_file: str) -> List[Dict[str, Any]]:
        """
        Parse PDB REMARK 350 assembly chain-group hints.

        Returns a list like:
        [{"biomolecule": "1", "biological_unit": "HEXAMERIC", "chains": ["B","D",...]}]
        """
        groups: List[Dict[str, Any]] = []
        current_biomolecule: Optional[str] = None
        current_biological_unit: str = ""
        try:
            lines = Path(pdb_file).read_text(encoding="utf-8", errors="ignore").splitlines()
        except Exception:
            return groups

        idx = 0
        while idx < len(lines):
            line = lines[idx]
            if not line.startswith("REMARK 350"):
                idx += 1
                continue

            payload = line[10:].strip()
            if payload.startswith("BIOMOLECULE:"):
                current_biomolecule = payload.split(":", 1)[1].strip()
                idx += 1
                continue

            if payload.startswith("AUTHOR DETERMINED BIOLOGICAL UNIT:"):
                current_biological_unit = payload.split(":", 1)[1].strip()
                idx += 1
                continue

            if payload.startswith("APPLY THE FOLLOWING TO CHAINS:"):
                chain_text = payload.split(":", 1)[1].strip()
                look_ahead = idx + 1
                while look_ahead < len(lines):
                    nxt = lines[look_ahead]
                    if not nxt.startswith("REMARK 350"):
                        break
                    nxt_payload = nxt[10:].strip()
                    if nxt_payload.startswith("AND CHAINS:"):
                        chain_text += ", " + nxt_payload.split(":", 1)[1].strip()
                        look_ahead += 1
                        continue
                    break

                parsed_chains: List[str] = []
                for token in re.split(r"[,\s]+", chain_text):
                    chain = str(token).strip().strip(",;")
                    if chain and chain not in parsed_chains:
                        parsed_chains.append(chain)

                if parsed_chains:
                    groups.append(
                        {
                            "biomolecule": current_biomolecule or str(len(groups) + 1),
                            "biological_unit": current_biological_unit,
                            "chains": parsed_chains,
                        }
                    )
                idx = look_ahead
                continue

            idx += 1

        return groups

    def save_hetatm_as_pdb(self, pdb_file: str, selected_hetatm: str,
                          chain_id: str, res_id: int,
                          pdb_id: Optional[str] = None,
                          output_filename: Optional[str] = None,
                          insertion_code: Optional[str] = None) -> str:
        """Extract one unambiguous ligand instance without inventing chemistry."""
        from docking.preparation.structure_contract import write_selection, residue_key
        import json
        source = Path(pdb_file)
        if source.suffix.lower() not in {".pdb", ".ent"}:
            raise ValueError("Use an authoritative instance SDF for mmCIF ligand extraction; PDB conversion loses chemical identity.")
        lines = source.read_text(encoding="utf-8").splitlines()
        candidates = {residue_key(line) for line in lines if line.startswith("HETATM") and line[17:20].strip() == selected_hetatm and line[21:22].strip() == str(chain_id).strip() and line[22:26].strip() == str(res_id)}
        if insertion_code is not None:
            candidates = {key for key in candidates if key[2] == insertion_code.strip()}
        if len(candidates) != 1:
            raise ValueError("Ligand instance is missing or ambiguous; specify chain, residue number and insertion code")
        key = next(iter(candidates))
        destination = Path(output_filename) if output_filename else self.output_dir / f"{pdb_id + '_' if pdb_id else ''}ligand_{selected_hetatm}_{chain_id}_{res_id}{key[2]}.pdb"
        metadata = write_selection(source, destination, residue=key)
        metadata.update({"reference_source": "cocrystal" if pdb_id else "unverified", "reference_pdb_id": pdb_id or "", "selected_instance": list(key), "reference_pose_file": str(destination.resolve())})
        destination.with_suffix(destination.suffix + ".preparation.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        return str(destination)

    def clean_pdb(self, pdb_file: str, to_remove_list: List[str],
                  pdb_id: Optional[str] = None, output_filename: Optional[str] = None,
                  keep_chain_id: Optional[str] = None, keep_chain_ids: Optional[List[str]] = None,
                  remove_instances: Optional[List[Tuple[str, str, str, str]]] = None) -> str:
        """Build a receptor from explicit residue instances, preserving atom metadata."""
        from docking.preparation.structure_contract import write_selection, residue_key
        import json
        source = Path(pdb_file)
        if source.suffix.lower() not in {".pdb", ".ent"}:
            raise ValueError("Explicit mmCIF assembly/model selection is required before PDB receptor export")
        lines = source.read_text(encoding="utf-8").splitlines()
        chains = set(keep_chain_ids or [])
        if keep_chain_id:
            chains.add(keep_chain_id)
        remove = {tuple(map(str, key)) for key in (remove_instances or [])}
        for line in lines:
            if line.startswith(("ATOM  ", "HETATM")):
                key = residue_key(line)
                if key[3] in to_remove_list or (chains and key[0] not in chains):
                    remove.add(key)
        destination = Path(output_filename) if output_filename else self.output_dir / f"{pdb_id + '_' if pdb_id else ''}cleaned.pdb"
        metadata = write_selection(source, destination, remove_instances=remove)
        metadata.update({"removed_instances": sorted(remove), "receptor_frame_id": metadata["reference_frame_id"], "source_accession": pdb_id or "", "assembly_status": "source_asymmetric_unit_not_verified"})
        destination.with_suffix(destination.suffix + ".preparation.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        return str(destination)

    def distance_based_interaction_detection(self, pdb_file: str, ligand_name: str,
                                             chain_id: str, res_id: int, cutoff: float = 5.0) -> List[np.ndarray]:
        report = extract_residue_level_coordinates(pdb_file, ligand_name, chain_id, res_id, cutoff)
        return report["all_coords"]

    def extract_active_site_coords(self, cleaned_pdb: str, ligand_name: str,
                                  chain_id: str, res_id: int,
                                  method: str = "distance") -> Tuple[np.ndarray, int]:
        # Site geometry must not depend on whether PLIP is installed.
        report = extract_residue_level_coordinates(cleaned_pdb, ligand_name, chain_id, res_id)
        return report["overall_center"], report["num_interacting_atoms"]

    def analyze_pocket_properties(self, cleaned_pdb: str, center_coords: np.ndarray,
                                  ligand_pdb: Optional[str] = None) -> Dict[str, Any]:
        """Report descriptive contacts; unsupported physical properties are unevaluated."""
        center = np.asarray(center_coords, dtype=float)
        if center.shape != (3,) or not np.isfinite(center).all():
            raise ValueError("Pocket center must contain three finite coordinates")
        structure = self._load_structure("protein", cleaned_pdb)
        nearby = set()
        hydrophobic = set()
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[0] != " ":
                        continue
                    for atom in residue.get_atoms():
                        xyz = np.asarray(atom.get_coord(), dtype=float)
                        if not np.isfinite(xyz).all():
                            raise ValueError("Nonfinite receptor coordinates")
                        if np.linalg.norm(xyz - center) <= 8.0:
                            key = (model.id, chain.id, residue.id)
                            nearby.add(key)
                            if residue.get_resname() in {"PHE", "TRP", "TYR", "LEU", "ILE", "VAL", "ALA", "MET"}:
                                hydrophobic.add(key)
                            break
        return {"center_x": float(center[0]), "center_y": float(center[1]), "center_z": float(center[2]),
                "nearby_residue_count": len(nearby), "nearby_hydrophobic_residues": len(hydrophobic),
                "descriptive_radius_A": 8.0, "pocket_volume_A3": None, "electrostatic_score": None,
                "druggability_score": None, "druggability_interpretation": "Not evaluated",
                "pocket_properties_status": "not_evaluated_no_validated_physical_model"}

    def generate_summary_report(self, results: Dict[str, Any], pdb_id: str) -> str:
        """
        Generate a comprehensive summary report.
        
        Args:
            results (Dict[str, Any]): Analysis results
            pdb_id (str): PDB identifier
            
        Returns:
            str: Path to generated report file
        """
        print("📋 Generating summary report...")
        
        try:
            report_filename = self.output_dir / f"{pdb_id}_pipeline_results.csv"
            
            # Prepare data for CSV
            report_data = []
            for key, value in results.items():
                report_data.append([key, value])
            
            # Save as CSV
            df = pd.DataFrame(report_data, columns=['Property', 'Value'])
            df.to_csv(report_filename, index=False)
            
            print(f"✓ Summary report saved as: {report_filename}")
            return str(report_filename)
        except Exception as e:
            print(f"❌ Error generating summary report: {e}")
            raise

def parse_plip_text_report(report_file: str, ligand_name: str, chain_id: str, res_id: int) -> Optional[Dict]:
    """
    Parse PLIP text report to extract interaction data reliably.
    
    Args:
        report_file: Path to PLIP text report
        ligand_name: Name of the ligand (HETATM)
        chain_id: Chain identifier
        res_id: Residue number
    
    Returns:
        Dictionary containing parsed interaction data
    """
    try:
        with open(report_file, 'r') as f:
            content = f.read()
        
        # Find the section for our ligand
        ligand_section = f"{ligand_name}:{chain_id}:{res_id}"
        if ligand_section not in content:
            return None
        
        # Extract the section for this ligand
        start_idx = content.find(ligand_section)
        if start_idx == -1:
            return None
        
        # Find the end of this ligand's section using regex
        import re
        
        # Look for the next ligand section (format: LIGAND_NAME:CHAIN:RESID)
        pattern = rf'\n[A-Z0-9]+:{chain_id}:\d+'
        matches = list(re.finditer(pattern, content[start_idx + len(ligand_section):]))
        
        if matches:
            # Use the first match as the end
            end_idx = start_idx + len(ligand_section) + matches[0].start()
        else:
            # No next ligand found, use end of file
            end_idx = len(content)
        
        ligand_content = content[start_idx:end_idx]
        
        # Parse interactions - comprehensive list of all PLIP interaction types
        interactions = {
            'hydrophobic': [],
            'hydrogen_bonds': [],
            'halogen_bonds': [],
            'pi_stacking': [],
            'salt_bridges': [],
            'water_bridges': [],
            'metal_complexes': [],
            'pi_cation': []
        }
        
        # Define interaction type configurations
        interaction_configs = {
            '**Hydrophobic Interactions**': {
                'key': 'hydrophobic',
                'min_columns': 8,
                'distance_col': 6
            },
            '**Hydrogen Bonds**': {
                'key': 'hydrogen_bonds', 
                'min_columns': 10,
                'distance_col': 7  # DIST_H-A
            },
            '**Halogen Bonds**': {
                'key': 'halogen_bonds',
                'min_columns': 8,
                'distance_col': 7  # DIST
            },
            '**Water Bridges**': {
                'key': 'water_bridges',
                'min_columns': 10,
                'distance_col': 6  # DIST_A-W
            },
            '**pi-Stacking**': {
                'key': 'pi_stacking',
                'min_columns': 10,
                'distance_col': 7  # CENTDIST
            },
            '**Salt Bridges**': {
                'key': 'salt_bridges',
                'min_columns': 8,
                'distance_col': 6  # Estimated
            },
            '**Metal Complexes**': {
                'key': 'metal_complexes',
                'min_columns': 8,
                'distance_col': 6  # Estimated
            },
            '**pi-Cation Interactions**': {
                'key': 'pi_cation',
                'min_columns': 8,
                'distance_col': 6  # Estimated
            }
        }
        
        # Unified parsing for all interaction types
        for section_name, config in interaction_configs.items():
            if section_name in ligand_content:
                section_content = ligand_content.split(section_name)[1]
                if "**" in section_content:
                    section_content = section_content.split("**")[0]
                
                lines = section_content.split('\n')
                for line in lines:
                    if '|' in line and line.strip().startswith('|') and not line.strip().startswith('|---'):
                        parts = [p.strip() for p in line.split('|') if p.strip()]
                        if len(parts) >= config['min_columns']:
                            try:
                                resnr = int(parts[0])
                                restype = parts[1]
                                reschain = parts[2]
                                if reschain == chain_id:
                                    distance = float(parts[config['distance_col']]) if len(parts) > config['distance_col'] else 0.0
                                    interactions[config['key']].append({
                                        'resnr': resnr,
                                        'restype': restype,
                                        'reschain': reschain,
                                        'distance': distance
                                    })
                            except (ValueError, IndexError):
                                continue
        
        # Check for any unknown interaction types
        all_sections = []
        for line in ligand_content.split('\n'):
            if line.strip().startswith('**') and line.strip().endswith('**'):
                all_sections.append(line.strip())
        
        known_sections = set(interaction_configs.keys())
        found_sections = set(all_sections)
        unknown_sections = found_sections - known_sections
        
        if unknown_sections:
            print(f"⚠️ Found unknown interaction types: {unknown_sections}")
            print(f"   Please update interaction_configs to handle these types")
        
        return interactions
        
    except Exception as e:
        print(f"⚠️ Error parsing PLIP report: {e}")
        return None


def extract_residue_level_coordinates(pdb_file: str, ligand_name: str,
                                     chain_id: str, res_id: int, cutoff: float = 5.0) -> Optional[Dict]:
    """Define a reproducible ligand envelope and minimum-distance residue contacts.

    This is a geometric contact report, not an interaction classifier. A ligand
    instance and one model must be unambiguous; atom order cannot change the box.
    """
    structure = MolecularDockingPipeline._load_structure("protein", pdb_file)
    models = list(structure)
    if len(models) != 1:
        raise ValueError("Select one receptor model explicitly before site extraction")
    ligands = [(chain, residue) for chain in models[0] for residue in chain
               if str(chain.id).strip() == str(chain_id).strip() and residue.get_resname() == ligand_name
               and residue.get_id()[1] == res_id and residue.get_id()[0] != " "]
    if len(ligands) != 1:
        raise ValueError("Reference ligand instance is missing or ambiguous")
    def heavy_coords(residue):
        coords = [np.asarray(atom.get_coord(), dtype=float) for atom in residue.get_atoms()
                  if str(getattr(atom, "element", "")).strip().upper() not in {"H", "D"}]
        if coords and not np.isfinite(np.asarray(coords)).all():
            raise ValueError("Nonfinite structure coordinates")
        return coords
    ligand_coords = np.asarray(heavy_coords(ligands[0][1]))
    if not len(ligand_coords):
        raise ValueError("Reference ligand has no heavy atoms")
    center = (ligand_coords.min(axis=0) + ligand_coords.max(axis=0)) / 2
    sizes = ligand_coords.max(axis=0) - ligand_coords.min(axis=0) + 10.0
    residue_coords = {}
    for chain in models[0]:
        for residue in chain:
            if residue.get_id()[0] != " ":
                continue
            contacting = [xyz for xyz in heavy_coords(residue) if np.min(np.linalg.norm(ligand_coords - xyz, axis=1)) <= cutoff]
            if contacting:
                key = f"{chain.id}:{residue.id[1]}{residue.id[2].strip()}:{residue.get_resname()}"
                residue_coords[key] = contacting
    all_coords = [xyz for key in sorted(residue_coords) for xyz in residue_coords[key]]
    return {"overall_center": center, "ligand_center": ligand_coords.mean(axis=0),
            "residue_averages": {key: np.mean(coords, axis=0) for key, coords in residue_coords.items()},
            "all_coords": all_coords, "num_interacting_residues": len(residue_coords),
            "num_interacting_atoms": len(all_coords), "plip_enhanced": False,
            "interaction_method": "minimum_heavy_atom_distance", "contact_cutoff_A": cutoff,
            "size_x": float(sizes[0]), "size_y": float(sizes[1]), "size_z": float(sizes[2])}
