#!/usr/bin/env python3
"""
Conservation Analysis Across Identity-Filtered MSA Files

This script analyzes conservation patterns across multiple MSA files that have been
filtered by sequence identity percentage. It identifies how conservation patterns
change at different identity thresholds.
"""

import os
import re
import glob
import numpy as np
import pandas as pd
from Bio import AlignIO
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

# Configuration
alignment_dir = "/home/cfneira1/ProteinMPNN/preanalysis/MSA_per_identity"  # Directory with filtered MSA files
output_dir = "identity_conservation_aa_results_1/"                # Where to save results
os.makedirs(output_dir, exist_ok=True)

# Regular expression to extract identity percentage from filenames
# Adjust this based on your file naming convention
# Example: protein_family_70pct.fasta (for 70% identity)
identity_pattern = re.compile(r'(\d+)-(\d+)alignment\.fasta')
# Example: protein_family_(\d+)-(\d+)pct\.fasta (for 70-90% identity)

# Define amino acid property groups
hydrophobic = "AILMFWYV"
polar = "NQST"
positive = "RHK"
negative = "DE"
special = "CGP"

# Calculate conservation using Shannon entropy
def calculate_conservation(alignment):
    """Calculate conservation for each position in an alignment"""
    conservation_scores = []
    residues_at_positions = []
    
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
            residues_at_positions.append('')
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
        
        # Identify most common residue
        most_common = max(aa_count.items(), key=lambda x: x[1]) if aa_count else ('X', 0)
        residues_at_positions.append(most_common[0])
    
    return conservation_scores, residues_at_positions

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

# Calculate property-based conservation
def calculate_property_conservation(alignment):
    """Calculate property conservation for each position in an alignment"""
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

def analyze_conservation_differences(conserved_positions_by_identity, conserved_residues_by_identity, 
                                     identity_percentages, alignment_lengths, output_dir):
    """Analyze differences in conserved positions between identity ranges"""
    
    # Create output directory
    diff_dir = os.path.join(output_dir, "conservation_differences")
    os.makedirs(diff_dir, exist_ok=True)
    
    # Create a comprehensive table of all conserved positions across all identity ranges
    all_positions = set()
    for pct in identity_percentages:
        all_positions.update(conserved_positions_by_identity[pct])
    
    # Create DataFrame to track conservation status at each identity level
    conservation_table = []
    
    for pos in sorted(all_positions):
        row = {'Position': pos}
        
        # Add conservation status and residue for each identity level
        for pct in sorted(identity_percentages):
            is_conserved = pos in conserved_positions_by_identity[pct]
            residue = conserved_residues_by_identity[pct].get(pos, '-')
            
            row[f"{pct}%_Conserved"] = is_conserved
            row[f"{pct}%_Residue"] = residue
        
        conservation_table.append(row)
    
    # Convert to DataFrame and save
    conservation_df = pd.DataFrame(conservation_table)
    conservation_df.to_csv(os.path.join(diff_dir, 'all_conserved_positions.csv'), index=False)
    
    return conservation_df

def find_changing_positions(conservation_df, identity_percentages, diff_dir):
    """Find positions where conservation changes between identity ranges"""
    
    # Identify positions that are conserved at some levels but not others
    changing_positions = []
    
    for _, row in conservation_df.iterrows():
        pos = row['Position']
        conservation_status = [row[f"{pct}%_Conserved"] for pct in sorted(identity_percentages)]
        
        # Check if conservation status changes (not all True or all False)
        if not all(conservation_status) and any(conservation_status):
            # Get the residues at each level
            residues = [row[f"{pct}%_Residue"] for pct in sorted(identity_percentages)]
            
            # Check if the residue changes among the conserved positions
            conserved_residues = [res for is_cons, res in zip(conservation_status, residues) if is_cons]
            residue_changes = len(set(conserved_residues)) > 1
            
            changing_positions.append({
                'Position': pos,
                'ConservationChanges': True,
                'ResidueChanges': residue_changes,
                **{f"{pct}%_Conserved": row[f"{pct}%_Conserved"] for pct in sorted(identity_percentages)},
                **{f"{pct}%_Residue": row[f"{pct}%_Residue"] for pct in sorted(identity_percentages)}
            })
    
    # Save the results
    if changing_positions:
        changing_df = pd.DataFrame(changing_positions)
        changing_df.to_csv(os.path.join(diff_dir, 'changing_conservation.csv'), index=False)
        
        # Create visualization
        plt.figure(figsize=(14, len(changing_positions) * 0.3 + 2))
        
        # Create a binary heatmap of conservation
        conservation_matrix = np.zeros((len(changing_positions), len(identity_percentages)))
        for i, pos_data in enumerate(changing_positions):
            for j, pct in enumerate(sorted(identity_percentages)):
                conservation_matrix[i, j] = 1 if pos_data[f"{pct}%_Conserved"] else 0
        
        # Plot heatmap
        ax = plt.gca()
        im = ax.imshow(conservation_matrix, cmap='viridis', aspect='auto')
        
        # Add text annotations for residues
        for i in range(len(changing_positions)):
            for j, pct in enumerate(sorted(identity_percentages)):
                if changing_positions[i][f"{pct}%_Conserved"]:
                    text = changing_positions[i][f"{pct}%_Residue"]
                    ax.text(j, i, text, ha="center", va="center", color="white", fontweight="bold")
        
        plt.colorbar(im, label='Conserved (1) vs Not Conserved (0)')
        plt.yticks(range(len(changing_positions)), [f"Pos {d['Position']}" for d in changing_positions])
        plt.xticks(range(len(identity_percentages)), [f"{pct}%" for pct in sorted(identity_percentages)])
        plt.xlabel('Identity Range')
        plt.ylabel('Position')
        plt.title('Positions with Changing Conservation Patterns')
        plt.tight_layout()
        plt.savefig(os.path.join(diff_dir, 'changing_conservation.png'), dpi=300)
        plt.close()
        
        return changing_df
    
    return None

def categorize_conservation_changes(changing_df, identity_percentages, diff_dir):
    """Categorize positions based on their conservation change patterns"""
    
    if changing_df is None or len(changing_df) == 0:
        return
    
    # Create categories
    categories = {
        'lost_at_low_identity': [], # Conserved at high identity but lost at lower levels
        'gained_at_high_identity': [], # Not conserved at low identity but gained at higher levels
        'fluctuating': [], # Conservation pattern fluctuates
        'residue_switch': [] # Position stays conserved but the residue changes
    }
    
    sorted_percentages = sorted(identity_percentages)
    
    for _, row in changing_df.iterrows():
        pos = row['Position']
        
        # Get conservation pattern (lowest to highest identity)
        pattern = [row[f"{pct}%_Conserved"] for pct in sorted_percentages]
        residues = [row[f"{pct}%_Residue"] for pct in sorted_percentages]
        
        # Check for patterns
        if pattern[0] and not pattern[-1]:
            # Lost conservation at higher identity - unusual but possible
            categories['lost_at_low_identity'].append(pos)
        elif not pattern[0] and pattern[-1]:
            # Gained conservation at higher identity
            categories['gained_at_high_identity'].append(pos)
        elif pattern.count(True) > 1 and pattern.count(False) > 0:
            # Check if residue changes when conserved
            conserved_residues = [res for is_cons, res in zip(pattern, residues) if is_cons]
            unique_residues = set(conserved_residues)
            
            if len(unique_residues) > 1:
                # The residue changes among conserved positions
                categories['residue_switch'].append(pos)
            else:
                # Conservation fluctuates
                categories['fluctuating'].append(pos)
    
    # Save categorization results
    with open(os.path.join(diff_dir, 'conservation_categories.txt'), 'w') as f:
        f.write("Conservation Pattern Categories\n\n")
        
        f.write("1. Positions that GAIN conservation at higher identity levels:\n")
        f.write(', '.join(map(str, categories['gained_at_high_identity'])) + "\n\n")
        
        f.write("2. Positions that LOSE conservation at higher identity levels:\n")
        f.write(', '.join(map(str, categories['lost_at_low_identity'])) + "\n\n")
        
        f.write("3. Positions with FLUCTUATING conservation patterns:\n")
        f.write(', '.join(map(str, categories['fluctuating'])) + "\n\n")
        
        f.write("4. Positions where the conserved RESIDUE CHANGES between identity levels:\n")
        f.write(', '.join(map(str, categories['residue_switch'])) + "\n\n")
    
    # Create a simple bar chart showing number of positions in each category
    plt.figure(figsize=(10, 6))
    category_names = ['Gained at\nHigh Identity', 'Lost at\nHigh Identity', 'Fluctuating\nConservation', 'Residue\nSwitches']
    counts = [len(categories['gained_at_high_identity']), 
              len(categories['lost_at_low_identity']),
              len(categories['fluctuating']), 
              len(categories['residue_switch'])]
    
    plt.bar(category_names, counts)
    plt.ylabel('Number of Positions')
    plt.title('Types of Conservation Changes Across Identity Levels')
    plt.tight_layout()
    plt.savefig(os.path.join(diff_dir, 'conservation_categories.png'), dpi=300)
    plt.close()

def create_individual_sequence_visualizations(conserved_positions_by_identity, conserved_residues_by_identity, 
                                             identity_percentages, reference_sequence, region_start=115, region_end=162, 
                                             output_dir=None):
    """
    Create individual visualizations of conserved residues for a specific region
    for each identity range separately, including the reference sequence.
    
    Parameters:
    - reference_sequence: String containing the reference protein sequence
    - region_start: Starting position of the region (inclusive)
    - region_end: Ending position of the region (inclusive)
    """
    if not output_dir:
        output_dir = "identity_conservation_results"
    
    # Create directory for individual visualizations
    region_dir = os.path.join(output_dir, "individual_region_analysis")
    os.makedirs(region_dir, exist_ok=True)
    
    # Define the region positions
    region_positions = list(range(region_start, region_end + 1))
    region_length = len(region_positions)
    
    # Extract the relevant part of the reference sequence
    # Adjust indices for 0-based indexing in Python strings
    ref_seq_region = reference_sequence[region_start-1:region_end]
    
    # Create visualization of just the reference sequence
    plt.figure(figsize=(region_length * 0.4, 2))
    
    # Draw reference sequence boxes
    for i, pos in enumerate(region_positions):
        # Get the residue from reference sequence (adjusting for 0-based indexing)
        ref_idx = pos - region_start
        if 0 <= ref_idx < len(ref_seq_region):
            residue = ref_seq_region[ref_idx]
        else:
            residue = '?'  # In case position is out of bounds
            
        # Draw rectangle for this position
        rect = plt.Rectangle((i, 0), 0.9, 1, facecolor='skyblue', alpha=0.7, edgecolor='black')
        plt.gca().add_patch(rect)
        
        # Add residue letter
        plt.text(i + 0.45, 0.5, residue, ha='center', va='center', 
                 color='black', fontweight='bold', fontsize=12)
        
        # Add position number below
        plt.text(i + 0.45, -0.2, str(pos), ha='center', va='center', 
                 fontsize=8, rotation=90)
    
    # Set axis limits and remove ticks
    plt.gca().set_xlim(-0.1, len(region_positions))
    plt.gca().set_ylim(-0.5, 1.5)
    plt.gca().set_axis_off()
    
    # Add title
    plt.title(f'Reference Sequence (Region {region_start}-{region_end})')
    
    # Save the figure
    plt.tight_layout()
    plt.savefig(os.path.join(region_dir, f'reference_sequence_region_{region_start}_{region_end}.png'), dpi=300)
    plt.close()
    
    # Sort identity percentages
    sorted_identities = sorted(identity_percentages)
    
    # Create a combined figure showing reference seq and conservation for all identities
    plt.figure(figsize=(region_length * 0.4, len(sorted_identities) + 2))
    
    # First row: Reference sequence
    for i, pos in enumerate(region_positions):
        # Get the residue from reference sequence
        ref_idx = pos - region_start
        if 0 <= ref_idx < len(ref_seq_region):
            residue = ref_seq_region[ref_idx]
        else:
            residue = '?'
            
        # Draw rectangle for reference sequence
        rect = plt.Rectangle((i, len(sorted_identities)), 0.9, 1, facecolor='skyblue', 
                             alpha=0.7, edgecolor='black')
        plt.gca().add_patch(rect)
        
        # Add residue letter
        plt.text(i + 0.45, len(sorted_identities) + 0.5, residue, ha='center', va='center', 
                 color='black', fontweight='bold', fontsize=12)
    
    # Draw conservation for each identity level
    for idx, pct in enumerate(sorted_identities):
        row = len(sorted_identities) - 1 - idx  # Reverse order (highest identity at top)
        
        # Get conserved positions and residues for this identity
        conserved_pos = conserved_positions_by_identity[pct]
        conserved_res = conserved_residues_by_identity[pct]
        
        for i, pos in enumerate(region_positions):
            is_conserved = pos in conserved_pos
            
            if is_conserved:
                residue = conserved_res.get(pos, '?')
                color = 'green'
                text_color = 'white'
                alpha = 0.9
            else:
                residue = '●'  # Symbol for non-conserved
                color = 'lightgray'
                text_color = 'black'
                alpha = 0.5
                
            # Draw rectangle for this position
            rect = plt.Rectangle((i, row), 0.9, 1, facecolor=color, alpha=alpha, edgecolor='black')
            plt.gca().add_patch(rect)
            
            # Add residue letter
            plt.text(i + 0.45, row + 0.5, residue, ha='center', va='center', 
                     color=text_color, fontweight='bold', fontsize=12)
    
    # Add identity labels on the left
    for idx, pct in enumerate(sorted_identities):
        row = len(sorted_identities) - 1 - idx
        plt.text(-0.5, row + 0.5, f"{pct}%", ha='right', va='center', fontsize=10)
    
    # Add 'Reference' label
    plt.text(-0.5, len(sorted_identities) + 0.5, "Ref", ha='right', va='center', fontsize=10)
    
    # Add position numbers at the bottom
    for i, pos in enumerate(region_positions):
        plt.text(i + 0.45, -0.2, str(pos), ha='center', va='center', 
                 fontsize=8, rotation=90)
    
    # Set axis limits and remove ticks
    plt.gca().set_xlim(-1, len(region_positions))
    plt.gca().set_ylim(-0.5, len(sorted_identities) + 1.5)
    plt.gca().set_axis_off()
    
    # Add title
    plt.title(f'Conservation Compared to Reference Sequence (Region {region_start}-{region_end})')
    
    # Add legend
    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor='skyblue', alpha=0.7, edgecolor='black', label='Reference Sequence'),
        plt.Rectangle((0, 0), 1, 1, facecolor='green', alpha=0.9, edgecolor='black', label='Conserved Position'),
        plt.Rectangle((0, 0), 1, 1, facecolor='lightgray', alpha=0.5, edgecolor='black', label='Non-conserved Position')
    ]
    plt.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.1, 1.1))
    
    # Save the combined figure
    plt.tight_layout()
    plt.savefig(os.path.join(region_dir, f'combined_conservation_region_{region_start}_{region_end}.png'), dpi=300)
    plt.close()
    
    # Create individual visualizations for each identity level
    for pct in sorted_identities:
        # Get conserved positions and residues for this identity
        conserved_pos = conserved_positions_by_identity[pct]
        conserved_res = conserved_residues_by_identity[pct]
        
        # Create figure showing reference and this identity level
        plt.figure(figsize=(region_length * 0.4, 3))
        
        # Reference sequence (top row)
        for i, pos in enumerate(region_positions):
            # Get the residue from reference sequence
            ref_idx = pos - region_start
            if 0 <= ref_idx < len(ref_seq_region):
                residue = ref_seq_region[ref_idx]
            else:
                residue = '?'
                
            # Draw rectangle for reference sequence
            rect = plt.Rectangle((i, 1), 0.9, 1, facecolor='skyblue', alpha=0.7, edgecolor='black')
            plt.gca().add_patch(rect)
            
            # Add residue letter
            plt.text(i + 0.45, 1.5, residue, ha='center', va='center', 
                     color='black', fontweight='bold', fontsize=12)
        
        # Conservation at this identity level (bottom row)
        for i, pos in enumerate(region_positions):
            is_conserved = pos in conserved_pos
            
            if is_conserved:
                residue = conserved_res.get(pos, '?')
                color = 'green'
                text_color = 'white'
                alpha = 0.9
            else:
                residue = '●'  # Symbol for non-conserved
                color = 'lightgray'
                text_color = 'black'
                alpha = 0.5
                
            # Draw rectangle for this position
            rect = plt.Rectangle((i, 0), 0.9, 1, facecolor=color, alpha=alpha, edgecolor='black')
            plt.gca().add_patch(rect)
            
            # Add residue letter
            plt.text(i + 0.45, 0.5, residue, ha='center', va='center', 
                     color=text_color, fontweight='bold', fontsize=12)
        
        # Add labels for rows
        plt.text(-0.5, 1.5, "Ref", ha='right', va='center', fontsize=10)
        plt.text(-0.5, 0.5, f"{pct}%", ha='right', va='center', fontsize=10)
        
        # Add position numbers at the bottom
        for i, pos in enumerate(region_positions):
            plt.text(i + 0.45, -0.2, str(pos), ha='center', va='center', 
                     fontsize=8, rotation=90)
        
        # Set axis limits and remove ticks
        plt.gca().set_xlim(-1, len(region_positions))
        plt.gca().set_ylim(-0.5, 2.5)
        plt.gca().set_axis_off()
        
        # Add title
        plt.title(f'Reference vs {pct}% Identity (Region {region_start}-{region_end})')
        
        # Add legend
        legend_elements = [
            plt.Rectangle((0, 0), 1, 1, facecolor='skyblue', alpha=0.7, edgecolor='black', label='Reference Sequence'),
            plt.Rectangle((0, 0), 1, 1, facecolor='green', alpha=0.9, edgecolor='black', label='Conserved Position'),
            plt.Rectangle((0, 0), 1, 1, facecolor='lightgray', alpha=0.5, edgecolor='black', label='Non-conserved Position')
        ]
        plt.legend(handles=legend_elements, loc='upper right')
        
        # Save the figure
        plt.tight_layout()
        plt.savefig(os.path.join(region_dir, f'reference_vs_{pct}_region_{region_start}_{region_end}.png'), dpi=300)
        plt.close()
    
    # Create summary text file including reference sequence
    with open(os.path.join(region_dir, f'region_{region_start}_{region_end}_summary.txt'), 'w') as f:
        f.write(f"Conservation Summary for Region {region_start}-{region_end}\n\n")
        
        f.write("Reference sequence for this region:\n")
        for i, pos in enumerate(region_positions):
            ref_idx = pos - region_start
            if 0 <= ref_idx < len(ref_seq_region):
                residue = ref_seq_region[ref_idx]
                f.write(f"{pos}:{residue} ")
            if (i + 1) % 10 == 0:
                f.write("\n")
        f.write("\n\n")
        
        for pct in sorted_identities:
            conserved_pos = [p for p in conserved_positions_by_identity[pct] if region_start <= p <= region_end]
            conserved_count = len(conserved_pos)
            conservation_pct = (conserved_count / region_length) * 100
            
            f.write(f"At {pct}% identity:\n")
            f.write(f"- {conserved_count} out of {region_length} positions conserved ({conservation_pct:.1f}%)\n")
            
            if conserved_pos:
                # Check which conserved residues match the reference
                matching = []
                different = []
                
                for pos in sorted(conserved_pos):
                    cons_residue = conserved_residues_by_identity[pct].get(pos, '?')
                    ref_idx = pos - region_start
                    
                    if 0 <= ref_idx < len(ref_seq_region):
                        ref_residue = ref_seq_region[ref_idx]
                        
                        if cons_residue == ref_residue:
                            matching.append(f"{pos}:{cons_residue}")
                        else:
                            different.append(f"{pos}:{ref_residue}->{cons_residue}")
                
                if matching:
                    f.write(f"- Conserved positions matching reference: {', '.join(matching)}\n")
                if different:
                    f.write(f"- Conserved positions differing from reference: {', '.join(different)}\n")
            f.write("\n")
    
    print(f"Individual region analysis with reference sequence for positions {region_start}-{region_end} complete.")


# Get all alignment files
alignment_files = glob.glob(os.path.join(alignment_dir, "*.fasta"))
print(f"Found {len(alignment_files)} alignment files to process")

# Extract identity percentages and sort files
identity_alignments = []
for aln_file in alignment_files:
    filename = os.path.basename(aln_file)
    match = identity_pattern.search(filename)
    if match:
        min_identity = int(match.group(1))
        max_identity = int(match.group(2))
        identity_pct = (min_identity + max_identity) // 2
        identity_alignments.append((identity_pct, min_identity, max_identity, aln_file))
    else:
        print(f"Warning: Could not extract identity percentage from {filename}")

# Sort by identity percentage
identity_alignments.sort()
print(f"Identified {len(identity_alignments)} files with identity percentages")

# Process each alignment file
conservation_by_identity = {}
property_by_identity = {}
conserved_positions_by_identity = {}
conserved_residues_by_identity = {}
alignment_lengths = {}

for identity_pct, min_identity, max_identity,aln_file in identity_alignments:
    filename = os.path.basename(aln_file)
    print(f"Processing {filename} ({min_identity}-{max_identity}% identity)")
    
    try:
        # Load the alignment
        alignment = AlignIO.read(aln_file, "fasta")
        num_sequences = len(alignment)
        alignment_length = alignment.get_alignment_length()
        alignment_lengths[identity_pct] = alignment_length
        
        print(f"  {num_sequences} sequences, {alignment_length} positions")
        
        # Calculate conservation
        conservation_scores, conserved_residues = calculate_conservation(alignment)
        # Calculate PSSM-like conservation
        pssm_scores = calculate_pssm_conservation(alignment)
        # Calculate property-based conservation
        property_scores, property_types = calculate_property_conservation(alignment)

         # Combine results
        results = pd.DataFrame({
            'Position': list(range(1, alignment.get_alignment_length() + 1)),
            'Shannon': conservation_scores,
            'PSSM': pssm_scores,
            'Property': property_scores,
            'PropertyType': property_types
        })
        
        # Store results
        conservation_by_identity[identity_pct] = conservation_scores
        property_by_identity[identity_pct] = property_scores
        
        # Identify highly conserved positions (threshold = 0.8)
        conserved_positions = [i+1 for i, score in enumerate(conservation_scores) if score >= 0.8]
        conserved_positions_by_identity[identity_pct] = conserved_positions
        
        # Store conserved residues
        residue_dict = {pos: conserved_residues[pos-1] for pos in conserved_positions}
        conserved_residues_by_identity[identity_pct] = residue_dict
        
        # Create result directory for this identity percentage
        identity_dir = os.path.join(output_dir, f"{identity_pct}%_identity")
        if not os.path.exists(identity_dir):
            print(f"  Creating directory: {identity_dir}")
        os.makedirs(identity_dir, exist_ok=True)
        
        # Save conservation scores
        results = pd.DataFrame({
            'Position': list(range(1, alignment_length + 1)),
            'Conservation': conservation_scores,
            'PropertyConservation': property_scores,
            'PropertyType': property_types,
            'ConservedResidue': conserved_residues
        })
        results.to_csv(os.path.join(identity_dir, 'conservation_scores.csv'), index=False)
        
        # Save conserved positions
        with open(os.path.join(identity_dir, 'conserved_positions.txt'), 'w') as f:
            f.write(f"Identity: {identity_pct}%\n")
            f.write(f"Total positions: {alignment_length}\n")
            f.write(f"Conserved positions (>= 0.8): {len(conserved_positions)}\n\n")
            f.write("Positions: " + ', '.join(map(str, conserved_positions)) + "\n\n")
            f.write("Conserved residues at these positions:\n")
            for pos in conserved_positions:
                f.write(f"Position {pos}: {conserved_residues[pos-1]}\n")
        
        # Create conservation profile plot
        plt.figure(figsize=(12, 6))
        plt.plot(range(1, alignment_length + 1), conservation_scores, label='Sequence Conservation')
        plt.plot(range(1, alignment_length + 1), property_scores, label='Property Conservation')
        plt.axhline(y=0.8, color='r', linestyle='--', label='Conservation Threshold (0.8)')
        plt.xlabel('Position')
        plt.ylabel('Conservation Score')
        plt.title(f'Conservation Profile - {identity_pct}% Identity')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(identity_dir, 'conservation_profile.png'), dpi=300)
        plt.close()
        
        print(f"  Analysis complete. Results saved to {identity_dir}")
        
    except Exception as e:
        print(f"  Error processing {aln_file}: {e}")

# Cross-identity analysis
print("\nPerforming cross-identity analysis...")

# Create comparative analysis directory
comp_dir = os.path.join(output_dir, "comparative_analysis")
os.makedirs(comp_dir, exist_ok=True)

# 1. Compare number of conserved positions across identity levels
identity_percentages = sorted(conservation_by_identity.keys())
num_conserved = [len(conserved_positions_by_identity[pct]) for pct in identity_percentages]
conserved_percentage = [(num / alignment_lengths[pct]) * 100 for pct, num in zip(identity_percentages, num_conserved)]

# Create conservation stats dataframe
conservation_stats = pd.DataFrame({
    'Identity': identity_percentages,
    'Alignment_Length': [alignment_lengths[pct] for pct in identity_percentages],
    'Conserved_Positions': num_conserved,
    'Conservation_Percentage': conserved_percentage
})
conservation_stats.to_csv(os.path.join(comp_dir, 'conservation_stats.csv'), index=False)

# Plot conservation percentage across identity levels
plt.figure(figsize=(10, 6))
plt.plot(identity_percentages, conserved_percentage, marker='o')
plt.xlabel('Sequence Identity Percentage')
plt.ylabel('Percentage of Conserved Positions')
plt.title('Conservation Level vs. Sequence Identity')
plt.grid(True, alpha=0.3)
plt.savefig(os.path.join(comp_dir, 'conservation_vs_identity.png'), dpi=300)
plt.close()

# 2. Visualize conservation profiles across identity levels
# Normalize profiles to same length if needed
max_length = max(alignment_lengths.values())
normalized_profiles = {}

for pct in identity_percentages:
    profile = conservation_by_identity[pct]
    length = alignment_lengths[pct]
    
    if length == max_length:
        normalized_profiles[pct] = profile
    else:
        # Normalize to max length
        x_original = np.linspace(0, 1, length)
        x_new = np.linspace(0, 1, max_length)
        normalized_profiles[pct] = np.interp(x_new, x_original, profile)

# Create heatmap of conservation profiles
profile_matrix = np.zeros((len(identity_percentages), max_length))
for i, pct in enumerate(identity_percentages):
    profile_matrix[i, :] = normalized_profiles[pct]

plt.figure(figsize=(12, 8))
im = plt.imshow(profile_matrix, aspect='auto', cmap='viridis')
plt.colorbar(im, label='Conservation Score')
plt.yticks(range(len(identity_percentages)), [f"{pct}%" for pct in identity_percentages])
plt.xlabel('Alignment Position (Normalized)')
plt.ylabel('Sequence Identity')
plt.title('Conservation Profiles Across Identity Levels')
plt.savefig(os.path.join(comp_dir, 'conservation_profile_heatmap.png'), dpi=300)
plt.close()

# 3. Track specific positions across identity levels
# Find common positions that are conserved at the highest identity level
if identity_percentages:
    highest_identity = max(identity_percentages)
    highest_conserved = set(conserved_positions_by_identity[highest_identity])
    
    # Track how these positions behave at lower identity levels
    position_tracking = []
    
    for pos in sorted(highest_conserved):
        pos_data = {'Position': pos}
        
        # Get residue at highest identity
        highest_residue = conserved_residues_by_identity[highest_identity].get(pos, 'X')
        pos_data['HighestResidue'] = highest_residue
        
        # Track conservation and residue at each identity level
        for pct in identity_percentages:
            # Normalize position if alignment lengths differ
            if alignment_lengths[pct] != alignment_lengths[highest_identity]:
                norm_pos = (pos - 1) / alignment_lengths[highest_identity]
                adj_pos = int(norm_pos * alignment_lengths[pct]) + 1
            else:
                adj_pos = pos
            
            # Check if position exists and is in range
            if 1 <= adj_pos <= alignment_lengths[pct]:
                adj_idx = adj_pos - 1  # Convert to 0-based index
                is_conserved = adj_pos in conserved_positions_by_identity[pct]
                conservation = conservation_by_identity[pct][adj_idx]
                residue = conserved_residues_by_identity[pct].get(adj_pos, 'X')
                
                pos_data[f"{pct}%_Conserved"] = is_conserved
                pos_data[f"{pct}%_Score"] = conservation
                pos_data[f"{pct}%_Residue"] = residue
            else:
                pos_data[f"{pct}%_Conserved"] = False
                pos_data[f"{pct}%_Score"] = 0
                pos_data[f"{pct}%_Residue"] = 'X'
        
        position_tracking.append(pos_data)
    
    # Create position tracking dataframe
    if position_tracking:
        tracking_df = pd.DataFrame(position_tracking)
        tracking_df.to_csv(os.path.join(comp_dir, 'position_tracking.csv'), index=False)
        
        # Create visualization of position conservation across identity levels
        plt.figure(figsize=(12, 8))
        
        # Prepare data for heatmap
        pos_matrix = np.zeros((len(highest_conserved), len(identity_percentages)))
        for i, pos_data in enumerate(position_tracking):
            for j, pct in enumerate(identity_percentages):
                pos_matrix[i, j] = pos_data[f"{pct}%_Score"]
        
        # Create heatmap
        im = plt.imshow(pos_matrix, aspect='auto', cmap='viridis')
        plt.colorbar(im, label='Conservation Score')
        plt.yticks(range(len(position_tracking)), [f"{d['Position']} ({d['HighestResidue']})" for d in position_tracking])
        plt.xticks(range(len(identity_percentages)), [f"{pct}%" for pct in identity_percentages])
        plt.xlabel('Sequence Identity')
        plt.ylabel('Position (Residue)')
        plt.title('Conservation of Positions Across Identity Levels')
        plt.tight_layout()
        plt.savefig(os.path.join(comp_dir, 'position_conservation_heatmap.png'), dpi=300)
        plt.close()

# 4. Find positions that become conserved at specific identity thresholds
conservation_emergence = defaultdict(list)

# Go through identity levels (lowest to highest)
for i, pct in enumerate(sorted(identity_percentages)):
    conserved = set(conserved_positions_by_identity[pct])
    
    # For first identity level, all conserved positions "emerge" at this level
    if i == 0:
        conservation_emergence[pct] = list(conserved)
    else:
        # Find positions that weren't conserved at lower identity levels
        prev_pct = sorted(identity_percentages)[i-1]
        prev_conserved = set(conserved_positions_by_identity[prev_pct])
        
        # Positions conserved at this level but not previous
        new_conserved = conserved - prev_conserved
        conservation_emergence[pct] = list(new_conserved)

# Create summary of emergence data
emergence_data = []
for pct, positions in sorted(conservation_emergence.items()):
    if positions:
        for pos in sorted(positions):
            # Get normalized position (for comparison across alignments)
            norm_pos = pos / alignment_lengths[pct]
            
            # Get conserved residue at this position
            residue = conserved_residues_by_identity[pct].get(pos, 'X')
            
            emergence_data.append({
                'Identity': pct,
                'Position': pos,
                'NormalizedPosition': norm_pos,
                'Residue': residue
            })

if emergence_data:
    emergence_df = pd.DataFrame(emergence_data)
    emergence_df.to_csv(os.path.join(comp_dir, 'conservation_emergence.csv'), index=False)
    
    # Create visualization of conservation emergence
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=emergence_df, x='Identity', y='Position', hue='Residue', s=100)
    plt.xlabel('Sequence Identity Threshold (%)')
    plt.ylabel('Position')
    plt.title('Emergence of Conserved Positions at Different Identity Thresholds')
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(comp_dir, 'conservation_emergence.png'), dpi=300)
    plt.close()

print(f"Cross-identity analysis complete. Results saved to {comp_dir}")

# 5. Analyze differences in conservation patterns
print("Analyzing differences in conservation between identity ranges...")

# Analyze conservation differences
conservation_df = analyze_conservation_differences(
    conserved_positions_by_identity, 
    conserved_residues_by_identity,
    identity_percentages, 
    alignment_lengths, 
    output_dir
)

# Find positions with changing conservation patterns
changing_df = find_changing_positions(conservation_df, identity_percentages, 
                                     os.path.join(output_dir, "conservation_differences"))

# Categorize conservation changes
categorize_conservation_changes(changing_df, identity_percentages, 
                               os.path.join(output_dir, "conservation_differences"))

# Add this after your cross-identity analysis
print(f"\nCreating individual visualizations for region (aa 115-162) with reference sequence...")

# Get the reference sequence - adjust this based on how your data is stored
# Option 1: If you have the full sequence as a string variable
reference_sequence = "MQSWYLLYCKRGQLQRAQEHLERQAVNCLAPMITLEKIVRGKRTAVSEPLFPNYLFVEFDPEVIHTTTINATRGVSHFVRFGASPAIVPSAVIHQLSVYKPKDIVDPATPYPGDKVIITEGAFEGFQAIFTEPDGEARSMLLLNLINKEIKHSVKNTEFRKL"  # Replace with your sequence

# Option 2: If your reference sequence is in your first alignment file
# with open(alignment_files[0], 'r') as f:
#     lines = f.readlines()
#     if lines[0].startswith('>'):
#         reference_sequence = lines[1].strip()

create_individual_sequence_visualizations(
    conserved_positions_by_identity,
    conserved_residues_by_identity,
    identity_percentages,
    reference_sequence,
    region_start=115,
    region_end=162,
    output_dir=output_dir
)


print("Conservation difference analysis complete!")