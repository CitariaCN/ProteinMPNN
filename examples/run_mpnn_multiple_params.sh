#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=32g
#SBATCH --partition=gpus
#SBATCH -t 0-01:00
#SBATCH -J test_tied_RfaH
#SBATCH -c 2

# Define base directories
folder_with_pdbs="../inputs/PDB_RfaH"
base_output_dir="../outputs/test_tied_RfaH"

# Create base output directory if it doesn't exist
if [ ! -d $base_output_dir ]; then
    mkdir -p $base_output_dir
fi

# Define paths for processed files - these will be shared across all runs
path_for_parsed_chains=$base_output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$base_output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$base_output_dir"/fixed_pdbs.jsonl"
path_for_tied_positions=$base_output_dir"/tied_pdbs.jsonl"

# Define input parameters
chains_to_design="A"
fixed_positions="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150"
tied_positions="150 151 152 153 154 155 156 157 158 159 160 161 162,150 151 152 153 154 155 156 157 158 159 160 161 162"

# Define PSSM thresholds to test
pssm_thresholds=(-0.5 -1.0 -1.5 -2.0 -2.5 -3.0)
pssm_jsonl_files=(
    "../inputs/pssm/rfah_pssm1.jsonl"
    "../inputs/pssm/rfah_pssm2.jsonl"
)
pssm_multi_values=(0.5 1.0 1.5)

# Process PDB files (only need to do this once for all parameter combinations)
echo "Preprocessing PDB files..."
python ../helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

python ../helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

python ../helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"

python ../helper_scripts/make_tied_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_tied_positions --chain_list "$chains_to_design" --position_list "$tied_positions"

# Run protein design with different parameter combinations
for pssm_jsonl in "${pssm_jsonl_files[@]}"; do
    # Extract filename without path for naming the output directory
    pssm_filename=$(basename "$pssm_jsonl" .jsonl)
    
    for pssm_multi in "${pssm_multi_values[@]}"; do
        for threshold in "${pssm_thresholds[@]}"; do
            # Create a unique output directory for each parameter combination
            output_dir="${base_output_dir}/${pssm_filename}_multi${pssm_multi}_thresh${threshold}"
            mkdir -p $output_dir
            
            echo "--------------------------------------"
            echo "Running design with the following parameters:"
            echo "PSSM file: $pssm_jsonl"
            echo "PSSM multiplier: $pssm_multi"
            echo "PSSM threshold: $threshold"
            echo "Output will be saved to: $output_dir"
            
            python ../protein_mpnn_run.py \
                --jsonl_path $path_for_parsed_chains \
                --chain_id_jsonl $path_for_assigned_chains \
                --fixed_positions_jsonl $path_for_fixed_positions \
                --tied_positions_jsonl $path_for_tied_positions \
                --pssm_jsonl $pssm_jsonl \
                --pssm_multi $pssm_multi \
                --pssm_threshold $threshold \
                --out_folder $output_dir \
                --num_seq_per_target 10 \
                --sampling_temp "0.1 0.2 0.3" \
                --seed 37 \
                --batch_size 1
            
            echo "Completed run with PSSM file: $pssm_jsonl, multiplier: $pssm_multi, threshold: $threshold"
            echo "--------------------------------------"
        done
    done
done

echo "All runs completed. Results are in $base_output_dir directory."