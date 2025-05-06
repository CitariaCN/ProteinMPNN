#!/usr/bin/env python3
"""
Position-wise Conservation Analysis Script

This script analyzes a single MSA file to calculate the conservation percentage
of each position in the reference sequence.
"""

import os
import numpy as np
import pandas as pd
from Bio import AlignIO
import matplotlib.pyplot as plt
import seaborn as sns

# Configuration
alignment_file = "/home/cfneira1/ProteinMPNN/preanalysis/MSA_per_identity/clustalo-100-55alignment.fasta"  # Update with your alignment file path
output_dir = "conservation_results/"      # Where to save results
os.makedirs(output_dir, exist_ok=True)

# Position numbering configuration
# The script will label positions from start_position to (start_position + length of alignment - 1)
# For example, if you have 48 positions and start_position=115, positions will be labeled 115-162

# Define amino acid property groups
hydrophobic = "AILMFWYV"
polar = "NQST"
positive = "RHK"
negative = "DE"
special = "CGP"

def calculate_conservation(alignment):
    """Calculate conservation for each position in an alignment"""
    conservation_scores = []
    residues_at_positions = []
    conservation_percentages = []
    
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
            conservation_percentages.append(0)
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
        
        # Identify most common residue and its percentage
        most_common = max(aa_count.items(), key=lambda x: x[1]) if aa_count else ('X', 0)
        residues_at_positions.append(most_common[0])
        
        # Calculate conservation percentage for the most common residue
        conservation_percentage = (most_common[1] / total) * 100 if total > 0 else 0
        conservation_percentages.append(conservation_percentage)
    
    return conservation_scores, residues_at_positions, conservation_percentages

def calculate_property_conservation(alignment):
    """Calculate property conservation for each position in an alignment"""
    conservation_scores = []
    property_types = []
    property_percentages = []
    
    for i in range(alignment.get_alignment_length()):
        column = [seq[i] for seq in alignment if seq[i] != '-']
        if not column:
            conservation_scores.append(0)
            property_types.append('none')
            property_percentages.append(0)
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
        property_percentage = (dominant_property[1] / len(column)) * 100
        
        conservation_scores.append(property_conservation)
        property_types.append(dominant_property[0])
        property_percentages.append(property_percentage)
    
    return conservation_scores, property_types, property_percentages

def visualize_conservation_reference(reference_sequence, conservation_percentages, 
                                    conserved_residues, property_percentages, 
                                    property_types, output_dir, start_position=115):
    """
    Visualize conservation data mapped to reference sequence
    
    Parameters:
    -----------
    start_position : int
        The position number to assign to the first residue in the sequence (default: 115)
    """
    # Create a dataframe with all the data
    data = pd.DataFrame({
        'Position': list(range(start_position, start_position + len(reference_sequence))),
        'Residue': list(reference_sequence),
        'ConservationScore': conservation_percentages[:len(reference_sequence)],
        'ConservedResidue': conserved_residues[:len(reference_sequence)],
        'PropertyScore': property_percentages[:len(reference_sequence)],
        'PropertyType': property_types[:len(reference_sequence)]
    })
    
    # Save to CSV
    data.to_csv(os.path.join(output_dir, 'position_conservation.csv'), index=False)
    
    # Create a visualization of conservation by position
    plt.figure(figsize=(15, 8))
    
    # Conservation percentage plot
    ax = sns.barplot(x='Position', y='ConservationScore', data=data, alpha=0.7)
    
    # Add residue labels
    for i, row in data.iterrows():
        pos_idx = row['Position'] - start_position  # Convert to 0-based index for plotting
        plt.text(i, 5, row['Residue'], ha='center', fontweight='bold')
        
        # Add conserved residue if different from reference
        if row['Residue'] != row['ConservedResidue'] and row['ConservationScore'] > 50:
            plt.text(i, row['ConservationScore'] + 3, row['ConservedResidue'], 
                    ha='center', color='red', fontweight='bold')
    
    # Set x-ticks to show actual position numbers
    plt.xticks(range(len(data)), data['Position'], rotation=90 if len(data) > 30 else 0)
    
    plt.xlabel('Position in Reference Sequence')
    plt.ylabel('Conservation Percentage (%)')
    plt.title(f'Position-wise Conservation in Reference Sequence (Positions {start_position}-{start_position + len(reference_sequence) - 1})')
    
    # Create custom legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='blue', alpha=0.7, label='Conservation %'),
        Patch(facecolor='white', label='Black: Reference Residue'),
        Patch(facecolor='white', label='Red: Most Conserved Residue (if different)')
    ]
    plt.legend(handles=legend_elements, loc='upper right')
    
    # Save and close
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'position_conservation.png'), dpi=300)
    plt.close()
    
    # Create second plot for property conservation
    plt.figure(figsize=(15, 8))
    
    # Create a categorical color map for property types
    property_colors = {
        'hydrophobic': 'gold',
        'polar': 'skyblue',
        'positive': 'tomato',
        'negative': 'mediumseagreen',
        'special': 'mediumpurple',
        'none': 'lightgray'
    }
    
    # Map property types to colors
    bar_colors = [property_colors[prop] for prop in data['PropertyType']]
    
    # Plot with property coloring
    x_pos = range(len(data))
    bars = plt.bar(x_pos, data['PropertyScore'], color=bar_colors, alpha=0.7)
    
    # Add reference residue labels
    for i, row in data.iterrows():
        plt.text(i, 5, row['Residue'], ha='center', fontweight='bold')
    
    # Set x-ticks to show actual position numbers
    plt.xticks(x_pos, data['Position'], rotation=90 if len(data) > 30 else 0)
    
    plt.xlabel('Position in Reference Sequence')
    plt.ylabel('Property Conservation Percentage (%)')
    plt.title(f'Amino Acid Property Conservation by Position (Positions {start_position}-{start_position + len(reference_sequence) - 1})')
    
    # Create property type legend
    legend_elements = [Patch(facecolor=color, alpha=0.7, label=prop) 
                       for prop, color in property_colors.items() if prop != 'none']
    plt.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'property_conservation.png'), dpi=300)
    plt.close()
    
    # Create a heatmap of position vs conservation
    plt.figure(figsize=(20, 6))
    
    # Prepare data for heatmap
    heatmap_data = data[['Position', 'Residue', 'ConservationScore']].copy()
    # Reshape for the heatmap
    heatmap_data = heatmap_data.pivot_table(index=['Residue'], columns='Position', values='ConservationScore')
    
    # Plot heatmap
    sns.heatmap(heatmap_data, cmap='viridis', annot=True, fmt='.1f', cbar_kws={'label': 'Conservation %'})
    plt.title(f'Conservation Percentage by Position and Residue Type (Positions {start_position}-{start_position + len(reference_sequence) - 1})')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'conservation_heatmap.png'), dpi=300)
    plt.close()
    
    return data

def validate_and_fix_alignment(alignment_file):
    """
    Validates the alignment and offers options to fix issues with uneven sequence lengths.
    Returns a valid alignment object or None if the alignment cannot be fixed.
    """
    try:
        # Try to read the alignment and check if all sequences have the same length
        alignment = AlignIO.read(alignment_file, "fasta")
        return alignment
    except ValueError as e:
        if "Sequences must all be the same length" in str(e):
            print("\nERROR: Sequences in the alignment have different lengths.")
            print("This is likely because your file is not a proper multiple sequence alignment.")
            print("\nPossible solutions:")
            
            # Option 1: Try to read as individual sequences and manually check
            from Bio import SeqIO
            try:
                sequences = list(SeqIO.parse(alignment_file, "fasta"))
                if not sequences:
                    print("Could not parse any sequences from the file.")
                    return None
                
                print(f"\nFound {len(sequences)} sequences with varying lengths:")
                lengths = {}
                for seq in sequences[:10]:  # Show first 10 for brevity
                    length = len(seq.seq)
                    lengths[length] = lengths.get(length, 0) + 1
                    if len(lengths) <= 5:  # Only show details for first few different lengths
                        print(f"  {seq.id}: {length} nucleotides/amino acids")
                
                print("\nLength distribution:")
                for length, count in sorted(lengths.items()):
                    print(f"  Length {length}: {count} sequences")
                
                print("\nYour file contains unaligned sequences. You need to:")
                print("1. Perform a multiple sequence alignment with tools like MUSCLE, MAFFT, or Clustal Omega")
                print("2. Or use an online service like EBI's MUSCLE (https://www.ebi.ac.uk/Tools/msa/muscle/)")
                print("3. Update the script with the path to your aligned file")
                
                # Try to identify if there might be something else wrong with the file
                if len(sequences) == 1:
                    print("\nNOTE: Your file contains only one sequence. Did you mean to provide an alignment with multiple sequences?")
                
                return None
                
            except Exception as inner_e:
                print(f"Failed to parse the file as individual sequences: {inner_e}")
                return None
        else:
            # Some other ValueError
            print(f"Error reading alignment: {e}")
            return None
    except Exception as e:
        print(f"Unexpected error reading alignment file: {e}")
        return None

def main():
    """Main execution function"""
    print(f"Processing alignment file: {alignment_file}")
    
    # Define start position (set to 115 as requested)
    start_position = 115
    
    try:
        # Validate and get the alignment
        alignment = validate_and_fix_alignment(alignment_file)
        if alignment is None:
            print("Could not proceed with analysis due to alignment issues.")
            return
        
        # Get reference sequence (first sequence in alignment)
        reference_sequence = str(alignment[0].seq).replace('-', '')
        
        # Number of sequences in alignment
        num_sequences = len(alignment)
        print(f"Alignment contains {num_sequences} sequences of length {alignment.get_alignment_length()}")
        print(f"Reference sequence length: {len(reference_sequence)}")
        print(f"Position numbering will start at {start_position} (ending at {start_position + len(reference_sequence) - 1})")
        
        # Calculate conservation scores
        conservation_scores, conserved_residues, conservation_percentages = calculate_conservation(alignment)
        
        # Calculate property conservation
        property_scores, property_types, property_percentages = calculate_property_conservation(alignment)
        
        # Visualize results
        results = visualize_conservation_reference(
            reference_sequence, 
            conservation_percentages, 
            conserved_residues,
            property_percentages,
            property_types,
            output_dir,
            start_position=start_position
        )
        
        # Identify highly conserved positions (e.g., > 80%)
        high_conservation = results[results['ConservationScore'] > 80].sort_values('Position')
        
        # Save highly conserved positions
        high_conservation.to_csv(os.path.join(output_dir, 'high_conservation_positions.csv'), index=False)
        
        # Print summary
        print(f"\nAnalysis complete!")
        print(f"Found {len(high_conservation)} highly conserved positions (>80% conservation)")
        
        # Print a sample of highly conserved positions
        if not high_conservation.empty:
            sample_size = min(10, len(high_conservation))
            print(f"\nSample of highly conserved positions:")
            print(high_conservation[['Position', 'Residue', 'ConservationScore', 'ConservedResidue']].head(sample_size))
        
        print(f"\nResults saved to {output_dir}")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()