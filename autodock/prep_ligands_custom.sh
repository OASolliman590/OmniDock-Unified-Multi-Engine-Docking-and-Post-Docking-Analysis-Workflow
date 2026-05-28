#!/usr/bin/env bash
set -euo pipefail

# Configurable directories (override with env vars).
# Example:
#   RAW_LIG=/path/to/raw_ligands OUT_LIG=/path/to/prepared_ligands ./prep_ligands_custom.sh
RAW_LIG="${RAW_LIG:-./1-Raw_Ligand}"
OUT_LIG="${OUT_LIG:-./2-Prepared_Ligands}"
mkdir -p "$OUT_LIG"

shopt -s nullglob

echo "🧪  Ligand preparation -----------------------------------"

# Process .mol2 and .sdf files (direct conversion)
for mol in "$RAW_LIG"/*.{mol2,sdf} 2>/dev/null; do
  [[ ! -f "$mol" ]] && continue
  base=${mol##*/}; base=${base%.*}
  out="$OUT_LIG/${base}.pdbqt"
  [[ -f $out ]] && { echo "skip $base"; continue; }
  echo "Processing: $base"
  mk_prepare_ligand.py -i "$mol" -o "$out"
done

# Process .mol files (direct conversion)
for mol in "$RAW_LIG"/*.mol 2>/dev/null; do
  [[ ! -f "$mol" ]] && continue
  base=${mol##*/}; base=${base%.*}
  out="$OUT_LIG/${base}.pdbqt"
  [[ -f $out ]] && { echo "skip $base"; continue; }
  echo "Processing: $base"
  mk_prepare_ligand.py -i "$mol" -o "$out"
done

# Process .pdb ligand files (convert to SDF first, then to PDBQT)
for pdb in "$RAW_LIG"/*.pdb 2>/dev/null; do
  [[ ! -f "$pdb" ]] && continue
  base=${pdb##*/}; base=${base%.*}
  out="$OUT_LIG/${base}.pdbqt"
  [[ -f $out ]] && { echo "skip $base"; continue; }
  echo "Processing: $base (PDB → SDF → PDBQT)"
  temp_sdf="$OUT_LIG/${base}_temp.sdf"
  # Convert PDB to SDF using obabel
  if obabel "$pdb" -O "$temp_sdf" -h 2>/dev/null; then
    # Convert SDF to PDBQT using Meeko
    if mk_prepare_ligand.py -i "$temp_sdf" -o "$out" 2>/dev/null; then
      rm -f "$temp_sdf"
    else
      rm -f "$temp_sdf" "$out"
      echo "⚠️  Failed to prepare: $base"
    fi
  else
    echo "⚠️  Failed to convert PDB to SDF: $base"
  fi
done

echo -e "\n✅  Finished."
echo "Ligands prepared: $(ls -1q $OUT_LIG/*.pdbqt 2>/dev/null | wc -l)"
