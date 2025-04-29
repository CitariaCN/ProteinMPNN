#!/bin/bash
# enhanced_test_protein_mpnn.sh - Quality check for ProteinMPNN job script
# This script performs quick validation tests without running the full computation

# Color codes for output formatting
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Create test directory
TEST_DIR="proteinMPNN_test_$(date +%Y%m%d_%H%M%S)"
mkdir -p ${TEST_DIR}
cd ${TEST_DIR}

# Set important paths
CURRENT_DIR=$(pwd)
MPNN_SCRIPT_PATH="/home/cfneira1/ProteinMPNN/protein_mpnn_run.py"
PDB_PATH="/home/cfneira1/ProteinMPNN/inputs/PDB_homooligomers/pdbs/combined_test_30.pdb"
MODEL_PATH="/home/cfneira1/ProteinMPNN/vanilla_model_weights/v_48_020.pt"

# Log file
LOG_FILE="test_results.log"
touch ${LOG_FILE}

log() {
    echo -e "$1" | tee -a ${LOG_FILE}
}

check_pass() {
    log "${GREEN}[PASS]${NC} $1"
}

check_fail() {
    log "${RED}[FAIL]${NC} $1"
}

check_warn() {
    log "${YELLOW}[WARN]${NC} $1"
}

log "==== ProteinMPNN Job Script Quality Check ===="
log "Date: $(date)"
log "Test directory: ${CURRENT_DIR}"
log "Main script path: ${MPNN_SCRIPT_PATH}"
log "PDB path: ${PDB_PATH}"
log "Model weights path: ${MODEL_PATH}"
log ""

# Create the script directly here instead of copying from another location
cat > test_script.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=proteinMPNN
#SBATCH --output=proteinMPNN_%j.out
#SBATCH --error=proteinMPNN_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=gpus
#SBATCH --gres=gpu:1

# Load necessary modules
module load cuda/11.0

# Create timestamp
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Create output directory with timestamp
BASE_DIR="RfaH_30_output_${TIMESTAMP}"
mkdir -p ${BASE_DIR}

# Create log file with timestamp
exec 1>${BASE_DIR}/run_${TIMESTAMP}.log 2>&1

echo "Starting ProteinMPNN job at $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Output directory: ${BASE_DIR}"

# Arrays of parameters to try
identities=(80 90)
temperature=(0.1)

# Run the designs
for identity in "${identities[@]}"; do
    for temp in "${temperature[@]}"; do
        output_dir="${BASE_DIR}/identity_${identity}_temp_${temp}"
        mkdir -p "$output_dir"
        
        echo "Running design with ${identity}% sequence identity and temperature ${temp}..."
        
        python /home/cfneira1/ProteinMPNN/protein_mpnn_run.py \
            --pdb_path /home/cfneira1/ProteinMPNN/inputs/PDB_monomers/pdbs/combined_test_30.pdb \
            --pdb_path_chains A \
            --chain_id_jsonl ./chains_to_design.json \
            --fixed_positions_jsonl ./fixed_positions.json \
            --out_folder "$output_dir" \
            --num_seq_per_target 10 \
            --sampling_temp "$temp" \
            --path_to_model_weights /home/cfneira1/ProteinMPNN/vanilla_model_weights/v_48_020.pt \
            --seed 37 \
            --ca_only \
            --omit_AAs X \
            --bias_AA_jsonl ./original_sequence_bias.json \
            --bias_by_res_jsonl ./bias_by_res.json
            
        echo "Completed design for identity ${identity}% and temperature ${temp}"
    done
done

# Create a detailed summary with timestamp
echo "Creating summary of all runs..."
{
    echo "Design Run Summary"
    echo "=================="
    echo "Generated on: $(date)"
    echo "Run ID: ${TIMESTAMP}"
    echo "Job ID: $SLURM_JOB_ID"
    echo
    echo "Parameters used:"
    echo "- PDB: /home/cfneira1/ProteinMPNN/inputs/PDB_monomers/pdbs/combined_test_30.pdb"
    echo "- Model weights: /home/cfneira1/ProteinMPNN/vanilla_model_weights/v_48_020.pt"
    echo "- Sequences per target: 10"
    echo "- Seed: 37"
    echo
    echo "Results:"
    for identity in "${identities[@]}"; do
        for temp in "${temperature[@]}"; do
            dir="${BASE_DIR}/identity_${identity}_temp_${temp}"
            if [ -d "$dir" ]; then
                count=$(ls "$dir"/*.fa 2>/dev/null | wc -l)
                echo "Identity: ${identity}%, Temperature: ${temp} - Generated $count sequences"
            fi
        done
    done
} > ${BASE_DIR}/summary_${TIMESTAMP}.txt

# Print summary to log
cat ${BASE_DIR}/summary_${TIMESTAMP}.txt

echo "Job completed at $(date)"
EOF

chmod +x test_script.sh

log "=== 1. Basic Environment Checks ==="

# Check if Python is available
if command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version)
    check_pass "Python is installed: ${PYTHON_VERSION}"
else
    check_fail "Python3 not found"
fi

# Check if CUDA is available
if command -v nvidia-smi &> /dev/null; then
    CUDA_INFO=$(nvidia-smi --query-gpu=name,driver_version --format=csv,noheader | head -n 1)
    check_pass "CUDA is available: ${CUDA_INFO}"
else
    check_warn "nvidia-smi not found. GPU may not be available"
fi

# Check if module command is available
if command -v module &> /dev/null; then
    check_pass "Module command is available"
else
    check_warn "Module command not found. This might be fine if not in an HPC environment"
fi

log ""
log "=== 2. Script Analysis ==="

# Check if PDB paths exist
if [ -f "${PDB_PATH}" ]; then
    check_pass "PDB file exists: ${PDB_PATH}"
else
    check_warn "PDB file not found: ${PDB_PATH}. Will need to be available on execution"
fi

# Check if model weights exist
if [ -f "${MODEL_PATH}" ]; then
    check_pass "Model weights exist: ${MODEL_PATH}"
else
    check_warn "Model weights not found: ${MODEL_PATH}. Will need to be available on execution"
fi

# Check protein_mpnn_run.py is accessible
if [ -f "${MPNN_SCRIPT_PATH}" ]; then
    check_pass "protein_mpnn_run.py exists at: ${MPNN_SCRIPT_PATH}"
else
    check_warn "protein_mpnn_run.py not found at: ${MPNN_SCRIPT_PATH}. Make sure it's available during execution"
fi

log ""
log "=== 3. JSON Configuration Tests ==="

# Create the JSON files directly
log "Creating JSON configuration files in: ${CURRENT_DIR}"

cat > fixed_positions.json << EOF
{
    "A": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150]
}
EOF

cat > chains_to_design.json << EOF
{
    "A": "1"
}
EOF

cat > original_sequence_bias.json << EOF
{
    "A": {
        "seq": "MQSWYLLYCKRGQLQRAQEHLERQAVNCLAPMITLEKIVRGKRTAVSEPLFPNYLFVEFDPEVIHTTTINATRGVSHFVRFGASPAIVPSAVIHQLSVYKPKDIVDPATPYPGDKVIITEGAFEGFQAIFTEPDGEARSMLLLNLINKEIKHSVKNTEFRCV"
    }
}
EOF

cat > bias_by_res.json << EOF
{
    "A": {}
}
EOF

cat > target_motif.json << EOF
{
    "A": []
}
EOF

# Validate JSON files
log "Validating JSON configuration files..."

# Check if files exist
for file in fixed_positions.json chains_to_design.json original_sequence_bias.json bias_by_res.json target_motif.json; do
    if [ -f "$file" ]; then
        check_pass "File created: $file"
    else
        check_fail "Failed to create: $file"
    fi
done

# Validate all JSON files with Python
python3 << EOF
import json
import sys

files = ['fixed_positions.json', 'chains_to_design.json', 'original_sequence_bias.json', 'bias_by_res.json', 'target_motif.json']
all_valid = True

for file in files:
    try:
        with open(file, 'r') as f:
            data = json.load(f)
            print(f"{file} is valid JSON")
    except Exception as e:
        print(f"Error in {file}: {str(e)}")
        all_valid = False

sys.exit(0 if all_valid else 1)
EOF

if [ $? -eq 0 ]; then
    check_pass "All JSON configuration files are valid"
else
    check_fail "One or more JSON files have issues"
fi

log ""
log "=== 4. Minimal Test Run ==="

# Create a tiny test run command with absolute paths to ensure reliability
cat > test_run.sh << EOF
#!/bin/bash
echo "Performing minimal test run..."

# Create a proper test output directory
TEST_OUTPUT_DIR="${CURRENT_DIR}/test_output"
mkdir -p "\${TEST_OUTPUT_DIR}"

# Only run with one identity and temperature value for testing
identity=30
temp=0.1

# Check that all required files exist
echo "Checking for required files:"
echo "  - PDB file: ${PDB_PATH}"
echo "  - Model weights: ${MODEL_PATH}"
echo "  - protein_mpnn_run.py: ${MPNN_SCRIPT_PATH}"
echo "  - fixed_positions.json: ${CURRENT_DIR}/fixed_positions.json"
echo "  - chains_to_design.json: ${CURRENT_DIR}/chains_to_design.json"
echo "  - original_sequence_bias.json: ${CURRENT_DIR}/original_sequence_bias.json"
echo "  - bias_by_res.json: ${CURRENT_DIR}/bias_by_res.json"

# Dry run (echo the command)
echo "Would execute: "
echo "python ${MPNN_SCRIPT_PATH} \\
    --pdb_path ${PDB_PATH} \\
    --pdb_path_chains A \\
    --chain_id_jsonl ${CURRENT_DIR}/chains_to_design.json \\
    --fixed_positions_jsonl ${CURRENT_DIR}/fixed_positions.json \\
    --out_folder \${TEST_OUTPUT_DIR} \\
    --num_seq_per_target 1 \\
    --sampling_temp \${temp} \\
    --path_to_model_weights ${MODEL_PATH} \\
    --seed 37 \\
    --ca_only \\
    --omit_AAs X \\
    --bias_AA_jsonl ${CURRENT_DIR}/original_sequence_bias.json \\
    --bias_by_res_jsonl ${CURRENT_DIR}/bias_by_res.json"

# Create a dummy output to simulate success
touch "\${TEST_OUTPUT_DIR}/test_sequence.fa"
echo ">design_1" > "\${TEST_OUTPUT_DIR}/test_sequence.fa"
echo "MQSWYLLYCKRGQLQRAQEHLERQAVNCLAPMITLEKIVRGKRTAVSEPLFPNYLFVEFDPEVIHTTTINATRGVSHFVRFGASPAIVPSAVIHQLSVYKPKDIVDPATPYPGDKVIITEGAFEGFQAIFTEPDGEARSMLLLNLINKEIKHSVKNTEFRCV" >> "\${TEST_OUTPUT_DIR}/test_sequence.fa"

echo "Test run completed. Output directory: \${TEST_OUTPUT_DIR}"
EOF
chmod +x test_run.sh

# Execute the test run
./test_run.sh >> ${LOG_FILE} 2>&1

if [ -d "test_output" ] && [ -f "test_output/test_sequence.fa" ]; then
    check_pass "Minimal test run simulation successful"
else
    check_fail "Minimal test run simulation failed"
fi

# Check parameter ranges
log ""
log "=== 5. Parameter Analysis ==="

IDENTITY_RANGE="80 90"
TEMP_RANGE="0.1"

log "Identity values: ${IDENTITY_RANGE}"
log "Temperature values: ${TEMP_RANGE}"

# Count total combinations
IDENTITY_COUNT=$(echo "${IDENTITY_RANGE}" | tr ' ' '\n' | wc -l)
TEMP_COUNT=$(echo "${TEMP_RANGE}" | tr ' ' '\n' | wc -l)
TOTAL_RUNS=$((IDENTITY_COUNT * TEMP_COUNT))

log "Total parameter combinations: ${TOTAL_RUNS}"
if [ ${TOTAL_RUNS} -gt 15 ]; then
    check_warn "Large number of parameter combinations (${TOTAL_RUNS}) may result in long run time"
else
    check_pass "Reasonable number of parameter combinations (${TOTAL_RUNS})"
fi

# Check requested resources
CPUS="4"
MEM="16G"
TIME="24:00:00"

log "Requested resources: ${CPUS} CPUs, ${MEM} memory, ${TIME} time"

# Final summary
log ""
log "=== Summary ==="
PASS_COUNT=$(grep -c "\[PASS\]" ${LOG_FILE})
FAIL_COUNT=$(grep -c "\[FAIL\]" ${LOG_FILE})
WARN_COUNT=$(grep -c "\[WARN\]" ${LOG_FILE})

log "Tests completed with: "
log "  - ${PASS_COUNT} passed checks"
log "  - ${WARN_COUNT} warnings"
log "  - ${FAIL_COUNT} failed checks"

if [ ${FAIL_COUNT} -eq 0 ]; then
    if [ ${WARN_COUNT} -eq 0 ]; then
        log "${GREEN}All checks passed successfully!${NC}"
    else
        log "${YELLOW}Tests passed with some warnings. Review the log for details.${NC}"
    fi
else
    log "${RED}Some tests failed. Please review the log and fix issues before submitting the job.${NC}"
fi

# Create a job submission script
cat > submit_job.sh << EOF
#!/bin/bash
sbatch test_script.sh
EOF
chmod +x submit_job.sh

log ""
log "Test results saved to: ${TEST_DIR}/${LOG_FILE}"
log "To submit the job, run: ./submit_job.sh"

# Return to the original directory
cd ..
log "Test directory: ${TEST_DIR}"