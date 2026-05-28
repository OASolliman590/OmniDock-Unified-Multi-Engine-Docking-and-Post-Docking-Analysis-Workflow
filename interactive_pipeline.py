#!/usr/bin/env python3
"""
PDB Prepare Wizard - Interactive Pipeline
========================================

Interactive version of the PDB preparation pipeline with user-friendly prompts
and multi-PDB analysis capabilities.

Author: Molecular Docking Pipeline
Version: 2.1.0
"""

import os
import re
import sys
from pathlib import Path
from typing import Optional

# Import the core pipeline
from core_pipeline import MolecularDockingPipeline, EXCEL_AVAILABLE

if EXCEL_AVAILABLE:
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment

def get_user_choice(prompt: str, valid_choices: list, default: Optional[str] = None) -> str:
    """Get user choice with validation."""
    while True:
        if default:
            user_input = input(f"{prompt} (default: {default}): ").strip()
            if not user_input:
                return default
        else:
            user_input = input(f"{prompt}: ").strip()
        
        if user_input.lower() in [choice.lower() for choice in valid_choices]:
            return user_input.lower()
        
        print(f"❌ Invalid choice. Please select from: {', '.join(valid_choices)}")


def _parse_index_selection(raw: str, max_index: int):
    """
    Parse interactive index input supporting:
    - single index: 2
    - comma-separated: 1,3,4
    - ranges: 1-4
    - keyword: all
    Returns a de-duplicated, ordered list of zero-based indices.
    """
    value = str(raw or "").strip().lower()
    if not value:
        return []
    if value in {"all", "*"}:
        return list(range(max_index))

    selected = []
    seen = set()
    tokens = [token.strip() for token in re.split(r"[,\s]+", value) if token.strip()]

    for token in tokens:
        if "-" in token:
            parts = [part.strip() for part in token.split("-", 1)]
            if len(parts) != 2 or not parts[0].isdigit() or not parts[1].isdigit():
                return []
            start = int(parts[0])
            end = int(parts[1])
            if start <= 0 or end <= 0:
                return []
            step_range = range(min(start, end), max(start, end) + 1)
            for one_based in step_range:
                idx = one_based - 1
                if idx < 0 or idx >= max_index:
                    return []
                if idx not in seen:
                    selected.append(idx)
                    seen.add(idx)
            continue

        if not token.isdigit():
            return []
        idx = int(token) - 1
        if idx < 0 or idx >= max_index:
            return []
        if idx not in seen:
            selected.append(idx)
            seen.add(idx)

    return selected

def select_hetatm_interactive(hetatm_details, unique_hetatms):
    """Interactive HETATM selection with grouped display."""
    if not unique_hetatms:
        print("❌ No HETATMs found to select from")
        return None, None, None
    
    print("\n📋 Available HETATM types:")
    hetatm_counts = {}
    for resname, _, _, _ in hetatm_details:
        hetatm_counts[resname] = hetatm_counts.get(resname, 0) + 1
    
    for i, (resname, count) in enumerate(hetatm_counts.items(), 1):
        print(f"   {i}. {resname} ({count} instance(s))")
    
    # Get user selection by type
    while True:
        try:
            choice = input(f"\nSelect HETATM type (1-{len(hetatm_counts)}): ").strip()
            choice_num = int(choice) - 1
            if 0 <= choice_num < len(hetatm_counts):
                selected_type = list(hetatm_counts.keys())[choice_num]
                break
            else:
                print(f"❌ Please enter a number between 1 and {len(hetatm_counts)}")
        except ValueError:
            print("❌ Please enter a valid number")
    
    # Find instances of selected type
    instances = [(chain_id, res_id, residue) for resname, chain_id, res_id, residue 
                in hetatm_details if resname == selected_type]
    
    if len(instances) == 1:
        chain_id, res_id, residue = instances[0]
        print(f"✓ Selected: {selected_type} in PDB chain {chain_id}, residue {res_id}")
        return selected_type, chain_id, res_id
    else:
        print(f"\n📋 Multiple instances of {selected_type} found:")
        for i, (chain_id, res_id, _) in enumerate(instances, 1):
            print(f"   {i}. PDB chain {chain_id}, Residue {res_id}")
        
        while True:
            choice = input(
                f"Select instance (1-{len(instances)}; supports comma/range/all): "
            ).strip()
            selected_indices = _parse_index_selection(choice, len(instances))
            if not selected_indices:
                print(
                    f"❌ Please enter a valid selection between 1 and {len(instances)} "
                    "(examples: 1, 2, 1,3, 1-4, all)"
                )
                continue

            if len(selected_indices) > 1:
                summary = ", ".join(
                    f"{instances[idx][0]}:{instances[idx][1]}" for idx in selected_indices
                )
                print(
                    f"ℹ️  Multiple instances selected ({summary}). "
                    "This workflow analyzes one representative ligand instance, so the first selected instance will be used."
                )

            selected_idx = selected_indices[0]
            chain_id, res_id, _ = instances[selected_idx]
            print(f"✓ Selected: {selected_type} in PDB chain {chain_id}, residue {res_id}")
            return selected_type, chain_id, res_id

def select_receptor_chains_interactive(protein_chains, default_chains=None):
    """Choose one or more receptor chains from source PDB ATOM-record chains."""
    if not protein_chains:
        return []

    print("\n🧬 Protein chains detected in source PDB:")
    for i, chain in enumerate(protein_chains, 1):
        marker = " (default)" if default_chains and chain in default_chains else ""
        print(f"   {i}. Chain {chain}{marker}")

    default_indices = []
    if default_chains:
        for chain in default_chains:
            if chain in protein_chains:
                default_indices.append(str(protein_chains.index(chain) + 1))
    default_text = ",".join(default_indices) if default_indices else None

    while True:
        choice = input(
            f"Select receptor chain(s) to keep during cleaning (comma-separated 1-{len(protein_chains)})"
            + (f" [default {default_text}]" if default_text else "")
            + ": "
        ).strip()
        if not choice and default_text:
            return [protein_chains[int(idx) - 1] for idx in default_indices]
        try:
            parts = [part.strip() for part in choice.split(",") if part.strip()]
            indices = [int(part) - 1 for part in parts]
            if indices and all(0 <= idx < len(protein_chains) for idx in indices):
                selected = []
                seen = set()
                for idx in indices:
                    chain = protein_chains[idx]
                    if chain not in seen:
                        selected.append(chain)
                        seen.add(chain)
                return selected
        except ValueError:
            pass
        print(f"❌ Please enter a number between 1 and {len(protein_chains)}")


def choose_chain_restriction_mode(
    selected_chain: str,
    ligand_bound_chains,
    assembly_peer_chains=None,
):
    """
    Let users quickly choose how aggressively to reduce receptor chains.
    """
    assembly_peer_chains = assembly_peer_chains or []
    has_ligand_group = len(ligand_bound_chains) > 1
    has_assembly_group = len(assembly_peer_chains) > 1

    if not has_ligand_group and not has_assembly_group:
        return "manual"

    print("\n🔀 Symmetry-aware chain options:")
    print(f"   1. Keep only selected ligand chain ({selected_chain}) [fastest]")
    if has_ligand_group:
        print(f"   2. Keep all chains carrying this ligand type ({', '.join(ligand_bound_chains)})")
    if has_assembly_group:
        print(f"   3. Keep all chains in selected chain assembly ({', '.join(assembly_peer_chains)})")
    print("   4. Manually choose chain(s)")

    while True:
        choice = input("Choose mode (1-4, default 1): ").strip()
        if not choice or choice == "1":
            return "single"
        if choice == "2" and has_ligand_group:
            return "ligand_bound"
        if choice == "3" and has_assembly_group:
            return "assembly"
        if choice == "4":
            return "manual"
        print("❌ Please enter a valid option from the listed modes")

def get_cleaning_choice(unique_hetatms):
    """Interactive cleaning choice with smart options."""
    print("\n📝 Cleaning Options:")
    print("   1. Remove specific residues (enter comma-separated list)")
    print("   2. Remove ALL HETATMs (enter 'ALL')")
    print("   3. Remove common residues only (enter 'COMMON')")
    print("   4. Skip cleaning (enter 'SKIP')")
    print("   5. Keep full protein unchanged (enter 'UNCHANGED')")

    common_residues = ['HOH', 'NA', 'CL', 'SO4', 'CA', 'MG', 'ZN', 'FE', 'CU', 'MN']

    while True:
        choice = input("\nEnter your choice (1-5, comma-separated list, ALL, COMMON, SKIP, or UNCHANGED): ").strip().upper()
        if not choice:
            print("❌ Please enter a valid cleaning choice.")
            continue

        if choice in {'1', 'REMOVE'}:
            residues_input = input("Enter comma-separated residue names to remove: ").strip()
            return [r.strip().upper() for r in residues_input.split(',') if r.strip()]
        if choice in {'2', 'ALL'}:
            return unique_hetatms
        if choice in {'3', 'COMMON'}:
            return [r for r in common_residues if r in unique_hetatms]
        if choice in {'4', 'SKIP'}:
            return []
        if choice in {'5', 'UNCHANGED'}:
            return ["__PRESERVE_FULL_RECEPTOR__"]

        residues = [r.strip().upper() for r in choice.split(',') if r.strip()]
        if residues:
            return residues
        print("❌ Please enter 1-4, ALL, COMMON, SKIP, or a comma-separated residue list.")

def show_removal_summary(pdb_file, to_remove_list):
    """Show what will be removed before confirmation."""
    if not to_remove_list:
        print("✓ No residues will be removed")
        return True
    
    try:
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        
        removal_counts = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname()
                    if resname in to_remove_list:
                        removal_counts[resname] = removal_counts.get(resname, 0) + 1
        
        if not removal_counts:
            print("✓ No matching residues found to remove")
            return True
        
        print(f"\n📊 Removal Summary:")
        total_removed = sum(removal_counts.values())
        print(f"   Total residues to remove: {total_removed}")
        for resname, count in removal_counts.items():
            print(f"   - {resname}: {count} residues")
        
        confirm = get_user_choice("\nProceed with removal?", ['y', 'n'], 'y')
        return confirm == 'y'
        
    except Exception as e:
        print(f"⚠️  Could not generate removal summary: {e}")
        confirm = get_user_choice("Proceed with removal anyway?", ['y', 'n'], 'n')
        return confirm == 'y'

def run_single_pdb_analysis(pipeline, pdb_id, excel_workbook=None):
    """Run analysis for a single PDB."""
    print(f"\n🚀 Starting analysis for PDB {pdb_id}")
    print("=" * 50)
    
    try:
        # Step 1: Download PDB
        pdb_file = pipeline.fetch_pdb(pdb_id)
        
        # Step 2: Enumerate HETATMs
        hetatm_details, unique_hetatms = pipeline.enumerate_hetatms(pdb_file)
        if not unique_hetatms:
            print("⚠️  No HETATMs found. Proceeding with protein-only analysis.")
            selected_hetatm = None
            chain_id = None
            res_id = None
        else:
            # Step 3: Select ligand
            selected_hetatm, chain_id, res_id = select_hetatm_interactive(hetatm_details, unique_hetatms)
            
            if selected_hetatm:
                # Save ligand as separate PDB
                ligand_pdb = pipeline.save_hetatm_as_pdb(
                    pdb_file, selected_hetatm, chain_id, res_id, pdb_id=pdb_id
                )
        
        # Step 4: Extract active site coordinates BEFORE cleaning (if ligand selected)
        # This allows us to get binding site info even if we remove the ligand later
        results = {'pdb_id': pdb_id}
        coords = None
        num_atoms = 0
        
        if selected_hetatm:
            try:
                print("\n🔄 Extracting active site coordinates from original structure...")
                coords, num_atoms = pipeline.extract_active_site_coords(
                    pdb_file, selected_hetatm, chain_id, res_id  # Use original PDB with ligand
                )
                results.update({
                    'selected_ligand': f"{selected_hetatm}_{chain_id}_{res_id}",
                    'active_site_center_x': coords[0],
                    'active_site_center_y': coords[1],
                    'active_site_center_z': coords[2],
                    'interacting_atoms_count': num_atoms
                })
                print(f"✓ Active site coordinates extracted: X={coords[0]:.2f}, Y={coords[1]:.2f}, Z={coords[2]:.2f}")
            except Exception as e:
                print(f"⚠️  Active site analysis failed: {e}")
                print("Proceeding without active site coordinates.")
        
        # Step 5: Clean PDB (remove ligand and other unwanted residues)
        to_remove_list = get_cleaning_choice(unique_hetatms)
        
        preserve_full_receptor = "__PRESERVE_FULL_RECEPTOR__" in to_remove_list
        if preserve_full_receptor:
            to_remove_list = []

        if selected_hetatm and selected_hetatm in to_remove_list:
            print(f"\n⚠️  NOTE: You selected {selected_hetatm} as your ligand.")
            print("   Active site coordinates have been extracted from the original structure.")
            print("   Removing the ligand will create an apo (ligand-free) receptor for docking.")
            keep_ligand = get_user_choice("Remove ligand from cleaned structure? (y=apo receptor, n=keep ligand)", ['y', 'n'], 'y')
            if keep_ligand == 'n':
                to_remove_list.remove(selected_hetatm)

        keep_chain_id = None
        keep_chain_ids = None
        protein_chains = pipeline.list_protein_chains(pdb_file)
        ligand_bound_chains = []
        assembly_groups = pipeline.parse_remark350_chain_groups(pdb_file)
        assembly_peer_chains = []
        if selected_hetatm:
            ligand_bound_chains = sorted(
                {
                    str(chain).strip()
                    for resname, chain, _, _ in hetatm_details
                    if resname == selected_hetatm and str(chain).strip() in protein_chains
                }
            )
            for group in assembly_groups:
                group_chains = [chain for chain in group.get("chains", []) if chain in protein_chains]
                if chain_id in group_chains:
                    assembly_peer_chains = group_chains
                    break
        if preserve_full_receptor:
            print("✓ Full receptor unchanged mode enabled: cleaning step will be skipped")
        elif selected_hetatm and chain_id:
            if assembly_groups:
                print("\n🧩 Symmetry / biological assembly hints from PDB REMARK 350:")
                for group in assembly_groups:
                    biomolecule = group.get("biomolecule", "?")
                    biological_unit = str(group.get("biological_unit", "") or "").strip()
                    chains_in_model = [chain for chain in group.get("chains", []) if chain in protein_chains]
                    if not chains_in_model:
                        continue
                    suffix = f" [{biological_unit}]" if biological_unit else ""
                    print(f"   Biomolecule {biomolecule}{suffix}: {', '.join(chains_in_model)}")
                if assembly_peer_chains:
                    print(
                        f"   Selected chain {chain_id} belongs to assembly chain set: "
                        f"{', '.join(assembly_peer_chains)}"
                    )
            keep_selected_chain = get_user_choice(
                f"Keep only selected receptor chain(s) in cleaned structure? (y/n)",
                ['y', 'n'],
                'n',
            )
            if keep_selected_chain == 'y':
                mode = choose_chain_restriction_mode(
                    chain_id,
                    ligand_bound_chains,
                    assembly_peer_chains=assembly_peer_chains,
                )
                if mode == "single":
                    keep_chain_ids = [chain_id] if chain_id else []
                elif mode == "ligand_bound":
                    keep_chain_ids = ligand_bound_chains[:]
                elif mode == "assembly":
                    keep_chain_ids = assembly_peer_chains[:]
                else:
                    keep_chain_ids = select_receptor_chains_interactive(
                        protein_chains,
                        default_chains=[chain_id] if chain_id in protein_chains else None,
                    )
                if not keep_chain_ids and chain_id:
                    keep_chain_ids = [chain_id]
                if keep_chain_ids:
                    keep_chain_id = keep_chain_ids[0] if len(keep_chain_ids) == 1 else None
                    print(f"✓ Chain-restricted cleaning enabled for receptor chain(s): {', '.join(keep_chain_ids)}")
            else:
                print("✓ Full-receptor cleaning enabled (all chains kept, selected residues removed only)")
        
        if preserve_full_receptor:
            cleaned_pdb = pdb_file
        elif show_removal_summary(pdb_file, to_remove_list):
            cleaned_pdb = pipeline.clean_pdb(
                pdb_file,
                to_remove_list,
                pdb_id=pdb_id,
                keep_chain_id=keep_chain_id,
                keep_chain_ids=keep_chain_ids,
            )
        else:
            print("Cleaning cancelled. Using original structure.")
            cleaned_pdb = pdb_file
        
        # Step 6: Analyze pocket properties using cleaned structure (if coordinates were extracted)
        if selected_hetatm and coords is not None:
            try:
                # Use the coordinates extracted from original structure, but analyze cleaned structure
                pocket_results = pipeline.analyze_pocket_properties(cleaned_pdb, coords)
                results.update(pocket_results)
            except Exception as e:
                print(f"⚠️  Pocket analysis failed: {e}")
                print("Proceeding without pocket analysis.")
        
        # Step 7: Generate reports
        report_file = pipeline.generate_summary_report(results, pdb_id)
        
        # Add to Excel workbook if provided
        if excel_workbook and EXCEL_AVAILABLE:
            add_to_excel_workbook(excel_workbook, pdb_id, results)
        
        print(f"\n✅ Analysis completed successfully for PDB {pdb_id}!")
        return True
        
    except Exception as e:
        print(f"\n❌ Analysis failed for PDB {pdb_id}: {e}")
        return False

def add_to_excel_workbook(workbook, pdb_id, results):
    """Upsert results into the Summary sheet by PDB ID."""
    try:
        # Create or get worksheet
        if f"Summary" not in workbook.sheetnames:
            ws = workbook.create_sheet("Summary")
            # Add headers
            headers = ["PDB_ID", "Property", "Value"]
            ws.append(headers)
            
            # Style headers
            for cell in ws[1]:
                    cell.font = Font(bold=True)
                    cell.fill = PatternFill(start_color="CCCCCC", end_color="CCCCCC", fill_type="solid")
                    cell.alignment = Alignment(horizontal="center")
        else:
            ws = workbook["Summary"]

        normalized_pdb = str(pdb_id).strip().upper()

        # Replace prior entries for this PDB so reruns overwrite stale values.
        for row_idx in range(ws.max_row, 1, -1):
            cell_value = ws.cell(row=row_idx, column=1).value
            existing_pdb = str(cell_value).strip().upper() if cell_value is not None else ""
            if existing_pdb == normalized_pdb:
                ws.delete_rows(row_idx, 1)
        
        # Add results
        for key, value in results.items():
            ws.append([normalized_pdb, key, value])
        
    except Exception as e:
        print(f"⚠️  Failed to add results to Excel: {e}")

def run_interactive_pipeline(output_dir: Optional[str] = None):
    """Run the interactive pipeline.
    
    Args:
        output_dir: Optional output directory path. If None, user will be prompted.
    """
    print("🔬 PDB Prepare Wizard - Interactive Pipeline")
    print("=" * 50)
    
    # Initialize pipeline
    if output_dir is None:
        output_dir = input("Enter output directory (default: pipeline_output): ").strip()
        if not output_dir:
            output_dir = "pipeline_output"
    else:
        print(f"✓ Using output directory: {output_dir}")
    
    pipeline = MolecularDockingPipeline(output_dir)
    
    # Check if user wants multi-PDB analysis
    multi_pdb = get_user_choice("Do you want to analyze multiple PDBs?", ['y', 'n'], 'n')
    
    excel_workbook = None
    if multi_pdb == 'y' and EXCEL_AVAILABLE:
        excel_workbook = Workbook()
        # Remove default sheet
        if "Sheet" in excel_workbook.sheetnames:
            excel_workbook.remove(excel_workbook["Sheet"])
        print("✓ Excel workbook created for multi-PDB analysis")
    
    pdb_count = 0
    while True:
        pdb_count += 1
        
        if multi_pdb == 'y':
            pdb_id = input(f"\nEnter PDB ID #{pdb_count} (e.g., 1ABC) or 'quit' to exit: ").strip().upper()
        else:
            pdb_id = input("\nEnter PDB ID (e.g., 1ABC) or 'quit' to exit: ").strip()

        if pdb_id.lower() == 'quit':
            break

        pdb_id = pdb_id.upper()
        
        if not pdb_id:
            print("❌ Please enter a valid PDB ID")
            pdb_count -= 1
            continue
        
        # Validate PDB ID format
        if len(pdb_id) != 4 or not pdb_id.isalnum():
            print("❌ PDB ID should be 4 characters (letters and numbers)")
            pdb_count -= 1
            continue
        
        # Run analysis
        success = run_single_pdb_analysis(pipeline, pdb_id, excel_workbook)
        
        if multi_pdb == 'y':
            if success:
                continue_analysis = get_user_choice("Do you want to analyze another PDB?", ['y', 'n'], 'n')
                if continue_analysis == 'n':
                    break
            else:
                continue_anyway = get_user_choice("Continue with next PDB despite this failure?", ['y', 'n'], 'y')
                if continue_anyway == 'n':
                    break
        else:
            break
    
    # Save Excel file for multi-PDB analysis
    if excel_workbook and EXCEL_AVAILABLE and multi_pdb == 'y':
        excel_path = Path(output_dir) / "multi_pdb_analysis.xlsx"
        excel_workbook.save(excel_path)
        print(f"\n📊 All results saved to: {excel_path}")
    
    if multi_pdb == 'y':
        print(f"\n🎉 Analysis completed! Processed {pdb_count - 1} PDB structure(s).")
    
    print("\n✓ Pipeline execution completed!")

if __name__ == "__main__":
    try:
        run_interactive_pipeline()
    except KeyboardInterrupt:
        print("\n\n⚠️  Analysis interrupted by user")
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        sys.exit(1)
