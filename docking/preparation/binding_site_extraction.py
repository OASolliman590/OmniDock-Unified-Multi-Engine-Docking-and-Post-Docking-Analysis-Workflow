"""Binding-site extraction for PDB preparation (Spec 038 R2).

Moved verbatim from core_pipeline.py; behaviour is unchanged. core_pipeline keeps
delegating methods and re-exports so existing import paths still resolve.
"""

import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

from docking.preparation.binding_site_center import resolve_binding_site_center


def centre_method_from_provenance(provenance) -> str:
    """Centre method recorded by the binding-site resolver (Spec 036 R6)."""
    if isinstance(provenance, dict):
        return str(provenance.get("center_method") or "")
    return ""


def load_protein_structure(structure_id: str, structure_file: str):
    structure_path = Path(structure_file)
    suffix = structure_path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(structure_id, str(structure_path))


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
        structure = load_protein_structure('protein', pdb_file)
        center_provenance = resolve_binding_site_center(
            structure_file=structure_path,
            accession=structure_path.stem,
            ligand_name=ligand_name,
            chain_id=chain_id,
            residue_number=res_id,
        )
        if center_provenance.get('status') != 'completed':
            raise ValueError(
                f"Binding-site center unavailable: {center_provenance.get('status')}:"
                f"{center_provenance.get('reason')}"
            )
        binding_site_center = np.asarray(center_provenance['binding_site_center'], dtype=float)

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
                                    # Contact centroid is descriptive; the approved box center is the ligand centroid.
                                    contact_residue_centroid = np.mean(all_coords, axis=0)
                                    
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
                                    print(f"✓ Binding site center: X={binding_site_center[0]:.2f}, Y={binding_site_center[1]:.2f}, Z={binding_site_center[2]:.2f}")
                                    
                                    return {
                                        'binding_site_center': binding_site_center,
                                        'overall_center': binding_site_center,
                                        'residue_averages': residue_averages,
                                        'all_coords': all_coords,
                                        'ligand_center': binding_site_center,
                                        'ligand_centroid': binding_site_center,
                                        'contact_residue_centroid': contact_residue_centroid,
                                        'center_provenance': center_provenance,
                                        'binding_site_center_method': centre_method_from_provenance(center_provenance),
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
        ligand_center = binding_site_center
        
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
        
        contact_residue_centroid = np.mean(all_coords, axis=0) if all_coords else None
        
        print(f"✓ Found {len(residue_averages)} interacting residues")
        print(f"✓ Total interacting atoms: {len(all_coords)}")
        print(f"✓ Binding site center: X={binding_site_center[0]:.2f}, Y={binding_site_center[1]:.2f}, Z={binding_site_center[2]:.2f}")
        
        return {
            'binding_site_center': binding_site_center,
            'overall_center': binding_site_center,
            'residue_averages': residue_averages,
            'all_coords': all_coords,
            'ligand_center': ligand_center,
            'ligand_centroid': ligand_center,
            'contact_residue_centroid': contact_residue_centroid,
            'center_provenance': center_provenance,
            'binding_site_center_method': centre_method_from_provenance(center_provenance),
            'num_interacting_residues': len(residue_averages),
            'num_interacting_atoms': len(all_coords),
            'plip_enhanced': False
        }
        
    except Exception as e:
        print(f"❌ Error in residue-level coordinate extraction: {e}")
        return None


def distance_based_interaction_detection(pdb_file: str, ligand_name: str,
                                         chain_id: str, res_id: int,
                                         cutoff: float = 5.0, *,
                                         load_structure=load_protein_structure) -> List[np.ndarray]:
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
        structure = load_structure('protein', pdb_file)
        
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


def extract_active_site_coords(cleaned_pdb: str, ligand_name: str,
                              chain_id: str, res_id: int,
                              method: str = 'distance', *,
                              plip_available: bool = False,
                              distance_detector=distance_based_interaction_detection) -> Tuple[np.ndarray, int]:
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
    
    if method == 'plip' and plip_available:
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
        coords = distance_detector(
            cleaned_pdb, ligand_name, chain_id, res_id
        )
    
    if not coords:
        raise ValueError("No coordinates extracted")
    
    center_result = resolve_binding_site_center(
        structure_file=Path(cleaned_pdb),
        accession=Path(cleaned_pdb).stem,
        ligand_name=ligand_name,
        chain_id=chain_id,
        residue_number=res_id,
    )
    if center_result.get('status') != 'completed':
        raise ValueError(
            f"Binding-site center unavailable: {center_result.get('status')}:"
            f"{center_result.get('reason')}"
        )
    avg_xyz = np.asarray(center_result['binding_site_center'], dtype=float)
    print(f"✓ Active site center: X={avg_xyz[0]:.2f}, Y={avg_xyz[1]:.2f}, Z={avg_xyz[2]:.2f}")
    
    return avg_xyz, len(coords)
