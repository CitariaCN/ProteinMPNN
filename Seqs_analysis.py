import re
import pandas as pd
from collections import defaultdict

def parse_mpnn_output(file_path):
    # Dictionary to store sequences by temperature
    sequences_by_temp = defaultdict(list)
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Get the native sequence first (the reference sequence)
    native_match = re.search(r'>combined_test_30.*?\n([A-Z]+)', content)
    native_sequence = native_match.group(1) if native_match else None
    
    # Find all sample sequences
    sample_pattern = r'>T=(\d+\.\d+), sample=(\d+), score=([0-9.]+), global_score=([0-9.]+), seq_recovery=([0-9.]+)\n([A-Z]+)'
    samples = re.finditer(sample_pattern, content)
    
    for sample in samples:
        temp = sample.group(1)
        sample_num = sample.group(2)
        score = sample.group(3)
        global_score = sample.group(4)
        seq_recovery = sample.group(5)
        sequence = sample.group(6)
        
        sequences_by_temp[temp].append({
            'sample': int(sample_num),
            'score': float(score),
            'global_score': float(global_score),
            'seq_recovery': float(seq_recovery),
            'sequence': sequence
        })
    
    # Convert to DataFrame
    dfs = {}
    for temp, seqs in sequences_by_temp.items():
        dfs[temp] = pd.DataFrame(seqs).sort_values('sample')
    
    return native_sequence, dfs

# Usage
file_path = "combined_test_30_pssm_threshold_-0.5"
native_sequence, sequences_by_temp = parse_mpnn_output(file_path)

print(f"Native sequence: {native_sequence}\n")

# Print sequences by temperature
for temp, df in sequences_by_temp.items():
    print(f"Temperature T={temp}:")
    for i, row in df.iterrows():
        print(f"  Sample {row['sample']}: {row['sequence']}")
    print()

# Optionally, save to CSV files
for temp, df in sequences_by_temp.items():
    df.to_csv(f"sequences_T{temp}.csv", index=False)