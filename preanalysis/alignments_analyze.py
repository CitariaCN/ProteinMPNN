import os
import glob
import numpy as np
from Bio import AlignIO
import matplotlib.pyplot as plt
import pandas as pd

# Define paths
alignment_dir = "/home/cfneira1/ProteinMPNN/preanalysis/MSA_per_identity"  # Directory with your alignment files
output_dir = "/home/cfneira1/ProteinMPNN/preanalysis/conservation_results"        # Where to save results
os.makedirs(output_dir, exist_ok=True)

# Get all alignment files (assuming FASTA format)
alignment_files = glob.glob(os.path.join(alignment_dir, "*.fasta"))
print(f"Found {len(alignment_files)} alignment files to process")

# Define amino acid property groups for method 3
hydrophobic = "AILMFWYV"
polar = "NQST"
positive = "RHK"
negative = "DE"
special = "CGP"

# METHOD 1: Shannon Entropy Conservation
def calculate_shannon_conservation(alignment):
    conservation_scores = []
    for i in range(alignment.get_alignment_length()):
        column = [seq[i] for seq in alignment]
        
        # Count amino acid frequencies (skipping gaps)
        aa_count = {}
        total = 0
        for aa in column:
            if aa != '-':
                aa_count[aa] = aa_count.get(aa, 0) + 1
                total += 1
        
        if total == 0:
            conservation_scores.append(0)
            continue
            
        # Calculate Shannon entropy
        entropy = 0
        for count in aa_count.values():
            p = count / total
            entropy -= p * np.log2(p)
        
        # Normalize to conservation score (1 = fully conserved)
        max_entropy = np.log2(min(20, total))
        conservation = 1 - (entropy / max_entropy) if max_entropy > 0 else 1.0
        conservation_scores.append(conservation)
    
    return conservation_scores

# METHOD 2: PSSM-like Position-Specific Conservation
def calculate_pssm_conservation(alignment):
    # Background frequencies of amino acids (approx. natural frequencies)
    bg_freqs = {
        'A': 0.074, 'R': 0.052, 'N': 0.045, 'D': 0.054, 'C': 0.025, 
        'Q': 0.034, 'E': 0.054, 'G': 0.074, 'H': 0.026, 'I': 0.068, 
        'L': 0.099, 'K': 0.058, 'M': 0.025, 'F': 0.047, 'P': 0.039, 
        'S': 0.057, 'T': 0.051, 'W': 0.013, 'Y': 0.032, 'V': 0.073
    }
    
    conservation_scores = []
    for i in range(alignment.get_alignment_length()):
        column = [seq[i] for seq in alignment if seq[i] != '-']
        if not column:
            conservation_scores.append(0)
            continue
            
        # Count amino acids in this column
        aa_counts = {}
        for aa in column:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        
        # Calculate position-specific scores
        total = len(column)
        max_score = 0
        for aa, count in aa_counts.items():
            if aa in bg_freqs:
                # Calculate log-odds score comparing observed to background
                freq = count / total
                odds = freq / bg_freqs[aa]
                score = np.log2(odds) if odds > 0 else -5
                max_score = max(max_score, score)
        
        # Normalize score to 0-1 range (heuristic)
        norm_score = min(1.0, max(0, (max_score + 5) / 10))
        conservation_scores.append(norm_score)
    
    return conservation_scores

# METHOD 3: Property-Based Conservation
def calculate_property_conservation(alignment):
    conservation_scores = []
    property_types = []
    
    for i in range(alignment.get_alignment_length()):
        column = [seq[i] for seq in alignment if seq[i] != '-']
        if not column:
            conservation_scores.append(0)
            property_types.append('none')
            continue
        
        # Count amino acids by property
        property_counts = {
            'hydrophobic': sum(1 for aa in column if aa in hydrophobic),
            'polar': sum(1 for aa in column if aa in polar),
            'positive': sum(1 for aa in column if aa in positive),
            'negative': sum(1 for aa in column if aa in negative),
            'special': sum(1 for aa in column if aa in special)
        }
        
        # Find dominant property and its conservation score
        dominant_property = max(property_counts.items(), key=lambda x: x[1])
        property_conservation = dominant_property[1] / len(column)
        
        conservation_scores.append(property_conservation)
        property_types.append(dominant_property[0])
    
    return conservation_scores, property_types

# Process each alignment file
for aln_file in alignment_files:
    alignment_name = os.path.basename(aln_file).split('.')[0]
    print(f"Processing alignment: {alignment_name}")
    
    try:
        # Load the alignment
        alignment = AlignIO.read(aln_file, "fasta")
        print(f"  {len(alignment)} sequences, {alignment.get_alignment_length()} positions")
        
        # Create result directory for this alignment
        aln_result_dir = os.path.join(output_dir, alignment_name)
        os.makedirs(aln_result_dir, exist_ok=True)
        
        # Apply all three methods
        # Method 1: Shannon Entropy
        shannon_scores = calculate_shannon_conservation(alignment)
        
        # Method 2: PSSM-like scoring
        pssm_scores = calculate_pssm_conservation(alignment)
        
        # Method 3: Property-based conservation
        property_scores, property_types = calculate_property_conservation(alignment)
        
        # Combine results
        results = pd.DataFrame({
            'Position': list(range(1, alignment.get_alignment_length() + 1)),
            'Shannon': shannon_scores,
            'PSSM': pssm_scores,
            'Property': property_scores,
            'PropertyType': property_types
        })
        
        # Save results to CSV
        results.to_csv(os.path.join(aln_result_dir, 'conservation_scores.csv'), index=False)
        
        # Identify conserved positions with each method (threshold = 0.8)
        highly_conserved = {}
        highly_conserved['Shannon'] = results[results['Shannon'] >= 0.8]['Position'].tolist()
        highly_conserved['PSSM'] = results[results['PSSM'] >= 0.8]['Position'].tolist()
        highly_conserved['Property'] = results[results['Property'] >= 0.8]['Position'].tolist()
        
        # Find positions identified by all methods
        consensus_positions = set(highly_conserved['Shannon']) & set(highly_conserved['PSSM']) & set(highly_conserved['Property'])
        
        # Save conserved positions to file
        with open(os.path.join(aln_result_dir, 'conserved_positions.txt'), 'w') as f:
            f.write(f"Alignment: {alignment_name}\n")
            f.write(f"Total positions: {alignment.get_alignment_length()}\n\n")
            
            f.write("Shannon Entropy Method:\n")
            f.write(f"Conserved positions (>= 0.8): {len(highly_conserved['Shannon'])}\n")
            f.write(', '.join(map(str, highly_conserved['Shannon'])) + "\n\n")
            
            f.write("PSSM Method:\n")
            f.write(f"Conserved positions (>= 0.8): {len(highly_conserved['PSSM'])}\n")
            f.write(', '.join(map(str, highly_conserved['PSSM'])) + "\n\n")
            
            f.write("Property Method:\n")
            f.write(f"Conserved positions (>= 0.8): {len(highly_conserved['Property'])}\n")
            f.write(', '.join(map(str, highly_conserved['Property'])) + "\n\n")
            
            f.write("Consensus (all methods):\n")
            f.write(f"Conserved positions: {len(consensus_positions)}\n")
            f.write(', '.join(map(str, sorted(consensus_positions))) + "\n")
        
        # Create conservation profile plot
        plt.figure(figsize=(12, 8))
        plt.plot(results['Position'], results['Shannon'], label='Shannon Entropy')
        plt.plot(results['Position'], results['PSSM'], label='PSSM')
        plt.plot(results['Position'], results['Property'], label='Property')
        plt.axhline(y=0.8, color='r', linestyle='--', label='Conservation Threshold (0.8)')
        plt.xlabel('Position')
        plt.ylabel('Conservation Score')
        plt.title(f'Conservation Profile - {alignment_name}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(aln_result_dir, 'conservation_profile.png'), dpi=300)
        plt.close()
        
        print(f"  Analysis complete. Results saved to {aln_result_dir}")
        
    except Exception as e:
        print(f"  Error processing {aln_file}: {e}")

print("All alignments processed!")