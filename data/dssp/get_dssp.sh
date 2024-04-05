#!/bin/bash

#use pydssp(https://github.com/ShintaroMinami/PyDSSP) to get secondary structure of Protein structurce

pdb_dir_path="/scratch/user/yxliu/dataset/landscape/ProteinGym/AF2_structures/"
# Declare an associative array

declare -A pdb_files



# Iterate over all PDB files in the directory

for pdb_file in "$pdb_dir_path"/*.pdb; do

  # Extract the ID from the file name

  pydssp "$pdb_dir_path"/*.pdb -o dssp.output
  id=$(basename "$pdb_file" .pdb)



  # Store the file path in the array, using the ID as the key

  pdb_files["$id"]="$pdb_file"

done



# Print all IDs and their corresponding file paths

for id in "${!pdb_files[@]}"; do

  echo "ID: $id"

  echo "File path: ${pdb_files[$id]}"

done

