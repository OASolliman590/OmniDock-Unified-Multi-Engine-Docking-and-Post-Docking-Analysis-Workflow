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
        else:
            parser = PDBParser(QUIET=True)
        return parser.get_structure(structure_id, str(structure_path))

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
                          output_filename: Optional[str] = None) -> Optional[str]:
        """
        Save a specific HETATM residue as a separate PDB file.
        
        Args:
            pdb_file (str): Path to source PDB file
            selected_hetatm (str): HETATM residue name
            chain_id (str): Chain identifier
            res_id (int): Residue ID
            pdb_id (str, optional): PDB ID to include in filename
            output_filename (str, optional): Output filename
            
        Returns:
            Optional[str]: Path to saved ligand PDB file, or None if failed
        """
        print(f"🔄 Saving HETATM {selected_hetatm}_{chain_id}_{res_id} as separate PDB...")
        
        try:
            structure = self._load_structure('protein', pdb_file)
            
            # Create new structure containing only the selected HETATM
            new_structure = Structure('ligand')
            new_model = Model(0)
            output_chain_id = chain_id[0] if len(str(chain_id)) > 1 else chain_id
            new_chain = Chain(output_chain_id)
            
            # Find and copy the specific HETATM residue
            found_residue = None
            for model in structure:
                for chain in model:
                    if chain.get_id() == chain_id:
                        for residue in chain:
                            if (residue.get_resname() == selected_hetatm and 
                                residue.get_id()[1] == res_id and
                                residue.get_id()[0] != ' '):
                                found_residue = residue
                                break
            
            if found_residue:
                new_chain.add(found_residue.copy())
                new_model.add(new_chain)
                new_structure.add(new_model)
                
                # Save as PDB
                if output_filename is None:
                    if pdb_id:
                        output_filename = self.output_dir / f"{pdb_id}_ligand_{selected_hetatm}_{chain_id}_{res_id}.pdb"
                    else:
                        output_filename = self.output_dir / f"ligand_{selected_hetatm}_{chain_id}_{res_id}.pdb"
                
                io = PDBIO()
                io.set_structure(self._normalize_structure_for_pdbio(new_structure))
                io.save(str(output_filename))
                from docking.preparation.ligand_quality import sanitize_ligand_pdb_file

                sanitize_ligand_pdb_file(Path(output_filename), Path(output_filename))
                print(f"✓ Ligand saved as: {output_filename}")
                return str(output_filename)
            else:
                print(f"❌ Could not find HETATM {selected_hetatm}_{chain_id}_{res_id}")
                return None
                
        except Exception as e:
            print(f"❌ Error saving HETATM: {e}")
            return None

    def clean_pdb(
        self,
        pdb_file: str,
        to_remove_list: List[str],
        pdb_id: Optional[str] = None,
        output_filename: Optional[str] = None,
        keep_chain_id: Optional[str] = None,
        keep_chain_ids: Optional[List[str]] = None,
    ) -> str:
        """
        Clean PDB file by removing specified residues.
        
        Args:
            pdb_file (str): Path to input PDB file
            to_remove_list (List[str]): List of residue names to remove
            pdb_id (str, optional): PDB ID to include in filename
            output_filename (str, optional): Output filename
            
        Returns:
            str: Path to cleaned PDB file
        """
        print(f"🔄 Cleaning PDB file...")
        
        try:
            keep_chains: List[str] = []
            if keep_chain_ids:
                keep_chains.extend([str(chain).strip() for chain in keep_chain_ids if str(chain).strip()])
            keep_chain = str(keep_chain_id or "").strip()
            if keep_chain and keep_chain not in keep_chains:
                keep_chains.append(keep_chain)
            before_counts = self._cleaning_counts(str(pdb_file))

            if output_filename is None:
                if pdb_id:
                    output_filename = self.output_dir / f"{pdb_id}_cleaned.pdb"
                else:
                    output_filename = self.output_dir / "cleaned.pdb"
            output_filename = str(output_filename)

            cleaned_with_pdb_tools = False
            input_suffix = Path(pdb_file).suffix.lower()
            detected_pdb_tools: Optional[Dict[str, List[str]]] = None
            if input_suffix == ".pdb":
                detected_pdb_tools = self._find_pdb_tools()
                try:
                    if detected_pdb_tools:
                        cleaned_with_pdb_tools = self._clean_pdb_with_pdb_tools(
                            pdb_file=pdb_file,
                            output_filename=output_filename,
                            to_remove_list=to_remove_list,
                            keep_chain_id=keep_chain_id,
                            keep_chain_ids=keep_chains,
                            tool_paths=detected_pdb_tools,
                        )
                        if cleaned_with_pdb_tools:
                            print("✓ Cleaning backend: pdb-tools")
                    else:
                        print("ℹ️  pdb-tools not found in PATH; using Biopython backend")
                except Exception as pdb_tools_error:
                    print(f"⚠️  pdb-tools cleaning failed, falling back to Biopython: {pdb_tools_error}")

            if not cleaned_with_pdb_tools:
                structure = self._load_structure('protein', pdb_file)
                keep_chain_set = set(keep_chains)

                # Create a selector to keep only desired residues (and optional chain filter)
                class ResidueSelector(Select):
                    def accept_chain(self, chain):
                        if not keep_chain_set:
                            return True
                        chain_id = str(chain.get_id()).strip()
                        original_chain_id = str(chain.xtra.get("original_chain_id", "")).strip()
                        return chain_id in keep_chain_set or original_chain_id in keep_chain_set

                    def accept_residue(self, residue):
                        resname = residue.get_resname()
                        return resname not in to_remove_list

                io = PDBIO()
                io.set_structure(self._normalize_structure_for_pdbio(structure))
                io.save(str(output_filename), ResidueSelector())
                print("✓ Cleaning backend: Biopython")

            after_counts = self._cleaning_counts(str(output_filename))
            if keep_chains:
                print(f"✓ Chain filter applied: kept chain(s) {', '.join(keep_chains)}")
            print(
                "✓ Cleaning summary: "
                f"ATOM {before_counts['ATOM']} -> {after_counts['ATOM']}, "
                f"HETATM {before_counts['HETATM']} -> {after_counts['HETATM']}, "
                f"TER {before_counts['TER']} -> {after_counts['TER']}"
            )
            
            print(f"✓ Cleaned PDB saved as: {output_filename}")
            return str(output_filename)
        except Exception as e:
            print(f"❌ Error cleaning PDB: {e}")
            raise

    def distance_based_interaction_detection(self, pdb_file: str, ligand_name: str, 
                                          chain_id: str, res_id: int, 
                                          cutoff: float = 5.0) -> List[np.ndarray]:
        """
        Detect interacting residues using distance-based approach.
        
        Args:
            pdb_file (str): Path to PDB file
            ligand_name (str): Ligand residue name
            chain_id (str): Chain identifier
            res_id (int): Residue ID
            cutoff (float): Distance cutoff in Angstroms
            
        Returns:
            List[np.ndarray]: List of coordinates of interacting atoms
        """
        try:
            structure = self._load_structure('protein', pdb_file)
            
            # Find ligand atoms
            ligand_atoms = []
            for model in structure:
                for chain in model:
                    if chain.get_id() == chain_id:
                        for residue in chain:
                            if (residue.get_resname() == ligand_name and 
                                residue.get_id()[1] == res_id and
                                residue.get_id()[0] != ' '):
                                ligand_atoms = [atom for atom in residue.get_atoms()]
                                break
            
            if not ligand_atoms:
                raise ValueError(f"Ligand {ligand_name} not found in chain {chain_id}")
            
            # Find interacting atoms
            all_coords = []
            interacting_residues = set()
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.get_id()[0] == ' ':  # Protein residue
                            for res_atom in residue.get_atoms():
                                for lig_atom in ligand_atoms:
                                    distance = res_atom - lig_atom  # BioPython distance calculation
                                    if distance <= cutoff:
                                        all_coords.append(res_atom.get_coord())
                                        interacting_residues.add(residue.get_resname() + str(residue.get_id()[1]))
                                        break
            
            print(f"✓ Found {len(interacting_residues)} interacting residues")
            return all_coords
            
        except Exception as e:
            raise ValueError(f"Distance-based interaction detection failed: {e}")

    def extract_active_site_coords(self, cleaned_pdb: str, ligand_name: str, 
                                 chain_id: str, res_id: int, 
                                 method: str = 'distance') -> Tuple[np.ndarray, int]:
        """
        Extract active site coordinates using specified method.
        
        Args:
            cleaned_pdb (str): Path to cleaned PDB file
            ligand_name (str): Ligand residue name
            chain_id (str): Chain identifier
            res_id (int): Residue ID
            method (str): Method to use ('plip' or 'distance')
            
        Returns:
            Tuple[np.ndarray, int]: Active site center coordinates and number of atoms
        """
        print(f"🔄 Extracting active site coordinates using {method} method...")
        
        coords = []
        
        if method == 'plip' and self.plip_available:
            try:
                from plip.structure.preparation import PDBComplex
                from plip.exchange.report import BindingSiteReport

                my_mol = PDBComplex()
                my_mol.load_pdb(cleaned_pdb)
                my_mol.analyze()

                interactions = my_mol.interaction_sets
                if not interactions:
                    print("⚠️  No interactions found with PLIP, falling back to distance method")
                    method = 'distance'
                else:
                    for site in interactions.values():
                        report = BindingSiteReport(site)
                        all_interacting_residues = report.bs_res_interacting

                        for residue in all_interacting_residues:
                            for atom in residue.atoms:
                                coords.append(atom.get_coord())

                    print(f"✓ PLIP found {len(coords)} interacting atoms")
            except Exception as e:
                print(f"⚠️  PLIP failed: {e}, using distance method")
                method = 'distance'
        
        if method == 'distance' or len(coords) == 0:
            coords = self.distance_based_interaction_detection(
                cleaned_pdb, ligand_name, chain_id, res_id
            )
        
        if not coords:
            raise ValueError("No coordinates extracted")
        
        # Calculate average XYZ
        avg_xyz = np.mean(coords, axis=0)
        print(f"✓ Active site center: X={avg_xyz[0]:.2f}, Y={avg_xyz[1]:.2f}, Z={avg_xyz[2]:.2f}")
        
        return avg_xyz, len(coords)

    def analyze_pocket_properties(self, cleaned_pdb: str, center_coords: np.ndarray, 
                                ligand_pdb: Optional[str] = None) -> Dict[str, Any]:
        """
        Comprehensive pocket analysis using multiple approaches.
        
        Args:
            cleaned_pdb (str): Path to cleaned PDB file
            center_coords (np.ndarray): Pocket center coordinates
            ligand_pdb (str, optional): Path to ligand PDB file
            
        Returns:
            Dict[str, Any]: Dictionary containing pocket analysis results
        """
        print("🔄 Starting comprehensive pocket analysis...")
        
        results = {
            'center_x': center_coords[0],
            'center_y': center_coords[1], 
            'center_z': center_coords[2]
        }
        
        try:
            structure = self._load_structure('protein', cleaned_pdb)
            
            # Pocket Size and Shape Analysis
            print("📊 Analyzing pocket size and shape...")
            try:
                # Geometric estimation based on interaction sphere
                results['pocket_volume_A3'] = 4/3 * np.pi * (5.0**3)  # Assume 5Å interaction sphere
            except Exception as e:
                print(f"⚠️  Pocket volume analysis failed: {e}")
                results['pocket_volume_A3'] = 'N/A'
            
            # Electrostatic Potential Analysis
            print("⚡ Analyzing electrostatic potential...")
            try:
                charged_residues = {'ARG': 1, 'LYS': 1, 'HIS': 0.5, 'ASP': -1, 'GLU': -1}
                electrostatic_score = 0
                nearby_charges = 0
                
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if residue.get_id()[0] == ' ':  # Protein residue
                                resname = residue.get_resname()
                                if resname in charged_residues:
                                    # Check if residue is near the pocket center
                                    try:
                                        ca_atom = residue['CA']
                                        dist = np.linalg.norm(ca_atom.get_coord() - center_coords)
                                        if dist <= 10.0:  # Within 10Å of pocket center
                                            electrostatic_score += charged_residues[resname]
                                            nearby_charges += 1
                                    except:
                                        continue
                
                results['electrostatic_score'] = electrostatic_score
                results['nearby_charged_residues'] = nearby_charges
                print(f"✓ Electrostatic analysis: {electrostatic_score:.2f} (from {nearby_charges} charged residues)")
            except Exception as e:
                print(f"⚠️  Electrostatic analysis failed: {e}")
                results['electrostatic_score'] = 'N/A'
                results['nearby_charged_residues'] = 0
            
            # Hydrophobic Character Analysis
            print("🧪 Analyzing hydrophobic character...")
            try:
                hydrophobic_residues = {'PHE', 'TRP', 'TYR', 'LEU', 'ILE', 'VAL', 'ALA', 'MET'}
                hydrophobic_score = 0
                nearby_hydrophobic = 0
                
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if residue.get_id()[0] == ' ':  # Protein residue
                                resname = residue.get_resname()
                                if resname in hydrophobic_residues:
                                    try:
                                        ca_atom = residue['CA']
                                        dist = np.linalg.norm(ca_atom.get_coord() - center_coords)
                                        if dist <= 8.0:  # Within 8Å of pocket center
                                            hydrophobic_score += 1
                                            nearby_hydrophobic += 1
                                    except:
                                        continue
                
                results['hydrophobic_score'] = hydrophobic_score
                results['nearby_hydrophobic_residues'] = nearby_hydrophobic
                print(f"✓ Hydrophobic analysis: {hydrophobic_score} hydrophobic residues nearby")
            except Exception as e:
                print(f"⚠️  Hydrophobic analysis failed: {e}")
                results['hydrophobic_score'] = 'N/A'
                results['nearby_hydrophobic_residues'] = 0
            
            # Druggability Scoring
            print("💊 Calculating druggability score...")
            try:
                # Simple druggability scoring based on multiple factors
                druggability_factors = []
                
                # Factor 1: Pocket volume (normalized)
                if isinstance(results['pocket_volume_A3'], (int, float)):
                    vol_score = min(1.0, results['pocket_volume_A3'] / 1000.0)
                    druggability_factors.append(vol_score)
                
                # Factor 2: Hydrophobic character (normalized)
                if isinstance(results['hydrophobic_score'], (int, float)):
                    hydro_score = min(1.0, results['hydrophobic_score'] / 10.0)
                    druggability_factors.append(hydro_score)
                
                # Factor 3: Electrostatic balance (absolute value, inverted and normalized)
                if isinstance(results['electrostatic_score'], (int, float)):
                    elec_score = max(0.0, 1.0 - abs(results['electrostatic_score']) / 5.0)
                    druggability_factors.append(elec_score)
                
                if druggability_factors:
                    druggability_score = np.mean(druggability_factors)
                    results['druggability_score'] = round(druggability_score, 3)
                    
                    # Interpretation
                    if druggability_score >= 0.7:
                        interpretation = "Excellent"
                    elif druggability_score >= 0.5:
                        interpretation = "Good" 
                    elif druggability_score >= 0.3:
                        interpretation = "Moderate"
                    else:
                        interpretation = "Poor"
                    
                    results['druggability_interpretation'] = interpretation
                    print(f"✓ Druggability score: {druggability_score:.3f} ({interpretation})")
                else:
                    results['druggability_score'] = 'N/A'
                    results['druggability_interpretation'] = 'Unknown'
            except Exception as e:
                print(f"⚠️  Druggability scoring failed: {e}")
                results['druggability_score'] = 'N/A'
                results['druggability_interpretation'] = 'Unknown'
            
            print("✓ Pocket analysis completed successfully")
            return results
            
        except Exception as e:
            print(f"❌ Pocket analysis failed: {e}")
            raise

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
    """
    Extract binding site coordinates at residue level using PLIP's text report.
    This function provides reliable analysis using PLIP's text output format.
    
    Args:
        pdb_file (str): Path to PDB file
        ligand_name (str): Ligand residue name
        chain_id (str): Chain identifier
        res_id (int): Residue ID
        cutoff (float): Distance cutoff in Angstroms (fallback method)
        
    Returns:
        Optional[Dict]: Dictionary with coordinate analysis or None if failed
    """
    print(f"🔄 Extracting residue-level binding site coordinates for {ligand_name}...")
    
    try:
        structure_path = Path(pdb_file)
        structure = MolecularDockingPipeline._load_structure('protein', pdb_file)

        # Try PLIP text report analysis first
        if structure_path.suffix.lower() in {".cif", ".mmcif"}:
            print("⚠️  Skipping PLIP text analysis for mmCIF input; using distance-based fallback.")
        else:
            print(f"🔄 Using PLIP text report analysis...")
        
        # Create temporary directory for PLIP analysis
        import tempfile
        import subprocess
        import os
        
        if structure_path.suffix.lower() not in {".cif", ".mmcif"}:
            with tempfile.TemporaryDirectory() as temp_dir:
                # Copy PDB file to temp directory
                temp_pdb = os.path.join(temp_dir, os.path.basename(pdb_file))
                import shutil
                shutil.copy2(pdb_file, temp_pdb)
                
                # Run PLIP to generate text report
                try:
                    result = subprocess.run(['plip', '-f', temp_pdb, '-t'], 
                                          capture_output=True, text=True, cwd=temp_dir)
                    if result.returncode == 0:
                        # Find the report file
                        report_files = [f for f in os.listdir(temp_dir) if f.endswith('report.txt')]
                        if report_files:
                            report_file = os.path.join(temp_dir, report_files[0])
                            
                            # Parse the report
                            interactions = parse_plip_text_report(report_file, ligand_name, chain_id, res_id)
                            
                            if interactions:
                                interacting_residues = set()
                                all_coords = []
                                
                                # Collect coordinates from all interaction types
                                for interaction_type, interaction_list in interactions.items():
                                    for interaction in interaction_list:
                                        reskey = f"{interaction['restype']}_{interaction['resnr']}"
                                        interacting_residues.add(reskey)
                                        
                                        # Find the residue in the structure
                                        for model in structure:
                                            for chain in model:
                                                if chain.get_id() == chain_id:
                                                    for residue in chain:
                                                        if (residue.get_resname() == interaction['restype'] and 
                                                            residue.get_id()[1] == interaction['resnr'] and
                                                            residue.get_id()[0] == ' '):
                                                            for atom in residue.get_atoms():
                                                                all_coords.append(atom.get_coord())
                                                            break
                                
                                if all_coords:
                                    # Calculate overall center
                                    overall_center = np.mean(all_coords, axis=0)
                                    
                                    # Calculate residue averages
                                    residue_averages = {}
                                    for reskey in interacting_residues:
                                        res_coords = []
                                        restype, resnr = reskey.split('_')
                                        resnr = int(resnr)
                                        
                                        for model in structure:
                                            for chain in model:
                                                if chain.get_id() == chain_id:
                                                    for residue in chain:
                                                        if (residue.get_resname() == restype and 
                                                            residue.get_id()[1] == resnr and
                                                            residue.get_id()[0] == ' '):
                                                            for atom in residue.get_atoms():
                                                                res_coords.append(atom.get_coord())
                                                            break
                                        
                                        if res_coords:
                                            residue_averages[reskey] = np.mean(res_coords, axis=0)
                                    
                                    print(f"✓ PLIP found {len(interacting_residues)} interacting residues")
                                    print(f"✓ Total interacting atoms: {len(all_coords)}")
                                    print(f"✓ Binding site center: X={overall_center[0]:.2f}, Y={overall_center[1]:.2f}, Z={overall_center[2]:.2f}")
                                    
                                    return {
                                        'overall_center': overall_center,
                                        'residue_averages': residue_averages,
                                        'all_coords': all_coords,
                                        'ligand_center': overall_center,
                                        'num_interacting_residues': len(interacting_residues),
                                        'num_interacting_atoms': len(all_coords),
                                        'plip_enhanced': True,
                                        'interaction_types': {
                                            'hydrophobic': len(interactions['hydrophobic']),
                                            'hydrogen_bonds': len(interactions['hydrogen_bonds']),
                                            'pi_stacking': len(interactions['pi_stacking']),
                                            'salt_bridges': len(interactions['salt_bridges']),
                                            'halogen_bonds': len(interactions['halogen_bonds']),
                                            'water_bridges': len(interactions['water_bridges']),
                                            'metal_complexes': len(interactions['metal_complexes']),
                                            'pi_cation': 0
                                        },
                                        'detailed_interactions': {
                                            'hydrophobic_residues': sorted(set([f"{i['restype']}_{i['resnr']}" for i in interactions['hydrophobic']])),
                                            'hydrogen_bond_residues': sorted(set([f"{i['restype']}_{i['resnr']}" for i in interactions['hydrogen_bonds']])),
                                            'halogen_bond_residues': sorted(set([f"{i['restype']}_{i['resnr']}" for i in interactions['halogen_bonds']])),
                                            'pi_stacking_residues': sorted(set([f"{i['restype']}_{i['resnr']}" for i in interactions['pi_stacking']])),
                                            'salt_bridge_residues': sorted(set([f"{i['restype']}_{i['resnr']}" for i in interactions['salt_bridges']]))
                                        }
                                    }
                                else:
                                    print(f"⚠️ No interacting atoms found for {ligand_name}:{chain_id}:{res_id}")
                            else:
                                print(f"⚠️ No interactions found in PLIP report for {ligand_name}:{chain_id}:{res_id}")
                        else:
                            print(f"⚠️ PLIP report file not found")
                    else:
                        print(f"⚠️ PLIP command failed: {result.stderr}")
                except FileNotFoundError:
                    print(f"⚠️ PLIP command not found, using fallback method")
                except Exception as e:
                    print(f"⚠️ PLIP analysis failed: {e}, using fallback method")
        
        # Fallback to distance-based method
        print(f"🔄 Using distance-based fallback method...")
        # Find the ligand residue
        ligand_residue = None
        for model in structure:
            for chain in model:
                if chain.get_id() == chain_id:
                    for residue in chain:
                        if (residue.get_resname() == ligand_name and 
                            residue.get_id()[1] == res_id and
                            residue.get_id()[0] != ' '):
                            ligand_residue = residue
                            break
        
        if not ligand_residue:
            raise ValueError(f"Ligand {ligand_name} not found in chain {chain_id}")
        
        # Get ligand atom coordinates
        ligand_coords = [atom.get_coord() for atom in ligand_residue.get_atoms()]
        ligand_center = np.mean(ligand_coords, axis=0)
        
        # Find interacting residues within cutoff distance
        residue_coords = {}
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[0] == ' ':  # Protein residue
                        resname = residue.get_resname()
                        resnum = residue.get_id()[1]
                        reskey = f"{resname}_{resnum}"
                        
                        # Check if any atom in this residue is within cutoff
                        for atom in residue.get_atoms():
                            dist = np.linalg.norm(atom.get_coord() - ligand_center)
                            if dist <= cutoff:
                                if reskey not in residue_coords:
                                    residue_coords[reskey] = []
                                residue_coords[reskey].append(atom.get_coord())
                                break
        
        # Calculate average coordinates for each interacting residue
        residue_averages = {}
        all_coords = []
        
        for reskey, coords in residue_coords.items():
            avg_coord = np.mean(coords, axis=0)
            residue_averages[reskey] = avg_coord
            all_coords.extend(coords)
        
        # Calculate overall average binding site center
        overall_center = np.mean(all_coords, axis=0) if all_coords else ligand_center
        
        print(f"✓ Found {len(residue_averages)} interacting residues")
        print(f"✓ Total interacting atoms: {len(all_coords)}")
        print(f"✓ Binding site center: X={overall_center[0]:.2f}, Y={overall_center[1]:.2f}, Z={overall_center[2]:.2f}")
        
        return {
            'overall_center': overall_center,
            'residue_averages': residue_averages,
            'all_coords': all_coords,
            'ligand_center': ligand_center,
            'num_interacting_residues': len(residue_averages),
            'num_interacting_atoms': len(all_coords),
            'plip_enhanced': False
        }
        
    except Exception as e:
        print(f"❌ Error in residue-level coordinate extraction: {e}")
        return None
