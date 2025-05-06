#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=32g
#SBATCH --partition=gpus
#SBATCH -t 0-01:00
#SBATCH -J mpnn_identity_positions
#SBATCH -c 2

# Define base directories
# Update these paths to match your actual directory structure
PROTEIN_MPNN_DIR="/home/cfneira1/ProteinMPNN"  # Update this to the actual path of your ProteinMPNN directory
folder_with_pdbs="${PROTEIN_MPNN_DIR}/inputs/PDB_RfaH"
base_output_dir="${PROTEIN_MPNN_DIR}/Seqs/identity_positions"

# Create base output directory if it doesn't exist
if [ ! -d $base_output_dir ]; then
    mkdir -p $base_output_dir
fi

# Path to the conservation summary file
conservation_file="${PROTEIN_MPNN_DIR}/preanalysis/identity_conservation_aa_results_1/individual_region_analysis/region_115_162_summary.txt"

# Define paths for intermediate files
path_for_parsed_chains="$base_output_dir/parsed_pdbs.jsonl"
path_for_assigned_chains="$base_output_dir/assigned_pdbs.jsonl"
chains_to_design="A"

# Create a temporary Python script to extract positions by identity
cat > extract_positions.py << 'EOF'
#!/usr/bin/env python3
import re
import sys
import json

def extract_conserved_positions(file_content):
    """Extract conserved positions at each identity percentage"""
    results = {}
    
    # Find all identity percentage sections
    pattern = r"At (\d+)% identity:\n- .*?\n((?:- .*?\n)+)"
    matches = re.finditer(pattern, file_content, re.MULTILINE)
    
    for match in matches:
        percentage = int(match.group(1))
        section = match.group(2)
        
        positions = []
        for line in section.split('\n'):
            if line.startswith("- Conserved positions"):
                pos_matches = re.finditer(r'(\d+):[A-Z](?:->)?[A-Z]?', line)
                for pos_match in pos_matches:
                    position = int(pos_match.group(1))
                    positions.append(position)
        
        # Remove duplicates and sort
        positions = sorted(list(set(positions)))
        results[percentage] = positions
    
    # If regular extraction fails, try alternative method
    if not results:
        percentages = re.findall(r'At (\d+)% identity:', file_content)
        percentages = [int(p) for p in percentages]
        
        for percentage in percentages:
            section_pattern = rf'At {percentage}% identity:(.*?)(?=At \d+% identity:|$)'
            section_match = re.search(section_pattern, file_content, re.DOTALL)
            
            if section_match:
                section = section_match.group(1)
                positions = re.findall(r'(\d+):[A-Z]', section)
                positions = sorted(list(set(map(int, positions))))
                results[percentage] = positions
    
    return results

def create_complete_position_lists(conserved_positions):
    """Create comprehensive position lists by adding positions 1-114"""
    complete_positions = {}
    positions_1_to_114 = list(range(1, 115))
    
    for percentage, positions in conserved_positions.items():
        combined_positions = sorted(positions_1_to_114 + positions)
        complete_positions[percentage] = combined_positions
    
    return complete_positions

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_json = sys.argv[2]
    
    with open(input_file, 'r') as f:
        file_content = f.read()
    
    conserved_positions = extract_conserved_positions(file_content)
    complete_positions = create_complete_position_lists(conserved_positions)
    
    # Write the results to a JSON file
    with open(output_json, 'w') as f:
        json.dump(complete_positions, f, indent=2)
    
    # Also print a summary
    print("Identity percentages and number of positions:")
    for percentage, positions in sorted(complete_positions.items()):
        print(f"{percentage}%: {len(positions)} positions")
EOF

# Make the extraction script executable
chmod +x extract_positions.py

# Run the extraction script to get positions by identity
positions_json="$base_output_dir/positions_by_identity.json"
echo "Extracting conserved positions from $conservation_file..."
python extract_positions.py "$conservation_file" "$positions_json"

# Define paths for intermediate files
path_for_parsed_chains="$base_output_dir/parsed_pdbs.jsonl"
path_for_assigned_chains="$base_output_dir/assigned_pdbs.jsonl"
chains_to_design="A"

# Preprocessing PDB files (only need to do this once)
echo "Preprocessing PDB files..."
python ${PROTEIN_MPNN_DIR}/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

python ${PROTEIN_MPNN_DIR}/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

# Define identity thresholds to use
identity_thresholds=(57 62 67 72 77 82 87 92 97)

# Read the positions from the JSON file
positions_data=$(cat "$positions_json")

# Main loop to run design for each identity threshold
for threshold in "${identity_thresholds[@]}"; do
    # Create output directory for this threshold
    output_dir="$base_output_dir/identity_${threshold}"
    mkdir -p "$output_dir"
    
    # Create a temporary Python script to generate the fixed_positions.jsonl for this threshold
    cat > generate_fixed_positions.py << EOF
#!/usr/bin/env python3
import json
import sys

# Read the positions data from the JSON file
positions_data = $positions_data

# Get positions for the current threshold
threshold = $threshold
if str(threshold) in positions_data:
    positions = positions_data[str(threshold)]
    positions_str = ' '.join(map(str, positions))
    print(positions_str)
else:
    print("Error: Threshold $threshold not found in positions data")
    sys.exit(1)
EOF

    # Make the script executable
    chmod +x generate_fixed_positions.py
    
    # Get the fixed positions for this threshold
    fixed_positions=$(python generate_fixed_positions.py)
    
    if [[ $fixed_positions == Error* ]]; then
        echo "$fixed_positions"
        continue
    fi
    
    # Path for fixed positions for this threshold
    path_for_fixed_positions="$output_dir/fixed_positions.jsonl"
    
    # Define tied positions (the last 13 positions 150-162)
    tied_positions="150 151 152 153 154 155 156 157 158 159 160 161 162"
    path_for_tied_positions="$output_dir/tied_positions.jsonl"
    
    echo "Creating fixed positions file for identity threshold: $threshold"
    python ${PROTEIN_MPNN_DIR}/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"
    
    echo "Creating tied positions file"
    python ${PROTEIN_MPNN_DIR}/helper_scripts/make_tied_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_tied_positions --chain_list "$chains_to_design" --position_list "$tied_positions"
    
    # Run ProteinMPNN
    echo "--------------------------------------"
    echo "Running design with identity threshold: $threshold%"
    echo "Using $(echo $fixed_positions | wc -w) fixed positions"
    echo "Output will be saved to: $output_dir"
    
    python ${PROTEIN_MPNN_DIR}/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --tied_positions_jsonl $path_for_tied_positions \
        --out_folder $output_dir \
        --num_seq_per_target 10 \
        --sampling_temp "0.1 0.2 0.3" \
        --seed 37 \
        --batch_size 1
    
    echo "--------------------------------------"
    
    # Clean up temporary scripts
    rm generate_fixed_positions.py
done

# Clean up the extraction script
rm extract_positions.py

echo "All runs completed. Results are in $base_output_dir directory."