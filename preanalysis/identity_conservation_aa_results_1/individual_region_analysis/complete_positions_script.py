#!/usr/bin/env python3
"""
Complete Position Analysis Script

This script extracts amino acid positions for each identity percentage from 
a conservation summary file and combines them with positions 1-114.
It generates comprehensive lists for each identity threshold.

Usage: 
    python3 complete_position_script.py input_file.txt

Output:
    Prints lists of positions for each identity percentage.
    Creates CSV files with the results.
"""

import re
import sys
import csv

def extract_conserved_positions(file_content):
    """
    Extract conserved positions at each identity percentage from the conservation summary.
    
    Args:
        file_content (str): Content of the conservation summary file
        
    Returns:
        dict: Dictionary with identity percentages as keys and lists of conserved positions as values
    """
    # Dictionary to store the results
    results = {}
    
    # Find all identity percentage sections
    pattern = r"At (\d+)% identity:\n- .*?\n((?:- .*?\n)+)"
    matches = re.finditer(pattern, file_content, re.MULTILINE)
    
    for match in matches:
        percentage = int(match.group(1))
        section = match.group(2)
        
        # Find all conserved positions mentioned in the section
        positions = []
        
        # Process each line that starts with "- Conserved positions"
        for line in section.split('\n'):
            if line.startswith("- Conserved positions"):
                # Extract position numbers from formats like "115:K->G, 117:I->R, ..."
                # or "162:L" (for positions matching reference)
                pos_matches = re.finditer(r'(\d+):[A-Z](?:->)?[A-Z]?', line)
                for pos_match in pos_matches:
                    position = int(pos_match.group(1))
                    positions.append(position)
        
        # Remove duplicates and sort
        positions = sorted(list(set(positions)))
        
        # Store in results dictionary
        results[percentage] = positions
    
    return results

def create_complete_position_lists(conserved_positions):
    """
    Create comprehensive position lists by adding positions 1-114 to 
    each identity percentage's conserved positions.
    
    Args:
        conserved_positions (dict): Dictionary with identity percentages as keys 
                                    and lists of conserved positions as values
        
    Returns:
        dict: Dictionary with identity percentages as keys and comprehensive 
              position lists as values
    """
    complete_positions = {}
    
    # Create a list of positions from 1 to 114
    positions_1_to_114 = list(range(1, 115))
    
    # For each identity percentage, add positions 1-114 to the conserved positions
    for percentage, positions in conserved_positions.items():
        # Combine the positions and sort
        combined_positions = sorted(positions_1_to_114 + positions)
        # Store in the complete_positions dictionary
        complete_positions[percentage] = combined_positions
    
    return complete_positions

def main():
    """Main function to process the input file and output results."""
    if len(sys.argv) < 2:
        print("Please provide an input file.")
        print("Usage: python3 complete_position_script.py input_file.txt")
        return
    
    input_file = sys.argv[1]
    
    try:
        with open(input_file, 'r') as f:
            file_content = f.read()
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        return
    except Exception as e:
        print(f"Error reading file: {e}")
        return
    
    # Extract the conserved positions (115-162)
    conserved_positions = extract_conserved_positions(file_content)
    
    if not conserved_positions:
        print("No conserved positions were found. The script might need adjustment for this file format.")
        
        # Alternative approach for difficult formats
        print("\nAttempting alternative extraction method...")
        
        # Extract identity percentages
        percentages = re.findall(r'At (\d+)% identity:', file_content)
        percentages = [int(p) for p in percentages]
        
        # Manual extraction for each percentage
        for percentage in percentages:
            # Find the section for this percentage
            section_pattern = rf'At {percentage}% identity:(.*?)(?=At \d+% identity:|$)'
            section_match = re.search(section_pattern, file_content, re.DOTALL)
            
            if section_match:
                section = section_match.group(1)
                
                # Extract all position numbers using a simpler pattern
                positions = re.findall(r'(\d+):[A-Z]', section)
                positions = sorted(list(set(map(int, positions))))
                
                # Store the results
                conserved_positions[percentage] = positions
    
    # Create complete position lists by adding positions 1-114
    complete_positions = create_complete_position_lists(conserved_positions)
    
    # Output the results
    print("\nComplete position lists for each identity percentage:")
    for percentage in sorted(complete_positions.keys()):
        positions = complete_positions[percentage]
        print(f"\n{percentage}% identity: {len(positions)} positions")
        # Print first few and last few positions as an example
        first_5 = ', '.join(map(str, positions[:5]))
        last_5 = ', '.join(map(str, positions[-5:]))
        print(f"First 5: {first_5} ... Last 5: {last_5}")
        print(f"Total: {len(positions)} positions")
    
    # Write results to CSV
    output_file = "complete_positions.csv"
    try:
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Identity_Percentage", "Positions"])
            for percentage in sorted(complete_positions.keys()):
                writer.writerow([percentage, ','.join(map(str, complete_positions[percentage]))])
        print(f"\nResults written to {output_file}")
    except Exception as e:
        print(f"Error writing to CSV: {e}")
        
    # Also write a more detailed CSV with each position on a separate row
    detailed_output = "complete_positions_detailed.csv"
    try:
        with open(detailed_output, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Identity_Percentage", "Position"])
            
            for percentage in sorted(complete_positions.keys()):
                for position in complete_positions[percentage]:
                    writer.writerow([percentage, position])
                    
        print(f"Detailed results written to {detailed_output}")
    except Exception as e:
        print(f"Error writing to detailed CSV: {e}")

if __name__ == "__main__":
    main()