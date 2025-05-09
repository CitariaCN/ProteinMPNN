import os
import numpy as np
from Bio import pairwise2
import datetime
import matplotlib.pyplot as plt

class MPNNResultsChecker:
    def __init__(self, reference_sequence, output_dir, region_start=None, region_end=None):
        self.reference_sequence = reference_sequence
        self.output_dir = output_dir
        self.region_start = region_start  # 1-based index
        self.region_end = region_end      # 1-based index
        
        # Create a meaningful name for this analysis
        if region_start and region_end:
            self.analysis_name = f"region_{region_start}_{region_end}"
        else:
            self.analysis_name = "full_sequence"
            
        # Create specific output directory for this analysis
        self.specific_output_dir = os.path.join(output_dir, self.analysis_name)
        if not os.path.exists(self.specific_output_dir):
            os.makedirs(self.specific_output_dir)

    def check_results(self, mpnn_results, target_identities):
        """Check the results of MPNN runs for the specified target identities."""
        print(f"\nChecking MPNN results for {self.analysis_name}...")
        
        # Verify the identity percentages
        all_results = self.verify_identity_percentages(mpnn_results, target_identities)
        
        # Check if all targets were processed
        for target in target_identities:
            if target not in all_results:
                print(f"  No results found for {target}% identity target.")
        
        # Generate summary report
        if all_results:
            self.generate_summary_report(all_results, target_identities)
        
        print(f"MPNN results check completed for {self.analysis_name}.")
        return all_results
    
    def verify_identity_percentages(self, mpnn_results, target_identities):
        """Verify the actual identity percentages of the generated sequences"""
        print(f"Verifying sequence identities for {self.analysis_name}...")
        
        # If region is specified, prepare the reference subsequence
        if self.region_start is not None and self.region_end is not None:
            # Convert from 1-based to 0-based indexing
            ref_start = self.region_start - 1
            ref_end = self.region_end
            region_reference = self.reference_sequence[ref_start:ref_end]
            print(f"Analyzing region from position {self.region_start} to {self.region_end}")
            print(f"Reference region: {region_reference}")
        else:
            region_reference = self.reference_sequence
            ref_start = 0
            ref_end = len(self.reference_sequence)
        
        all_results = {}
        
        for target in target_identities:
            output_dir = mpnn_results.get(target)
            if not output_dir or not os.path.exists(output_dir):
                print(f"  No results directory found for {target}% identity target.")
                continue
            
            # Check for all possible threshold files
            distance_thresholds = [30, 60, 90]
            all_sequences = []
            
            for threshold in distance_thresholds:
                fasta_file = os.path.join(output_dir, f"combined_test_{threshold}.fa")
                
                if not os.path.exists(fasta_file):
                    print(f"  No FASTA file found for {target}% identity target at threshold {threshold}.")
                    continue
                
                print(f"  Processing threshold file combined_test_{threshold}.fa for {target}% identity target...")
                
                try:
                    sequences_in_file = []  # Initialize here to avoid the UnboundLocalError
                    with open(fasta_file, 'r') as f:
                        current_seq = ""
                        for line in f:
                            if line.startswith(">"):
                                if current_seq:  # If we've collected a sequence, add it
                                    sequences_in_file.append(current_seq)
                                    current_seq = ""
                                continue
                            current_seq += line.strip()
                        
                        # Add the last sequence if there is one
                        if current_seq:
                            sequences_in_file.append(current_seq)
                        
                    print(f"    Found {len(sequences_in_file)} sequences in threshold file {threshold}")
                    all_sequences.extend(sequences_in_file)
                
                except Exception as e:
                    print(f"  Error processing threshold file {threshold} for {target}% identity: {e}")
            
            if not all_sequences:
                print(f"  No sequences found for {target}% identity target in any threshold file.")
                continue
            
            print(f"  Found a total of {len(all_sequences)} sequences for {target}% identity target across all thresholds.")
            
            # Calculate identities for the specific region or whole sequence
            identities = []
            valid_sequences = []
            
            for seq in all_sequences:
                try:
                    if self.region_start is not None and self.region_end is not None:
                        # Extract the same region from the generated sequence if it's long enough
                        if len(seq) >= self.region_end:
                            seq_region = seq[ref_start:ref_end]
                            alignment = pairwise2.align.globalxx(region_reference, seq_region, one_alignment_only=True)[0]
                            matches = sum(a == b for a, b in zip(alignment.seqA, alignment.seqB))
                            identity = (matches / len(alignment.seqA)) * 100
                            identities.append(identity)
                            valid_sequences.append(seq)
                        else:
                            print(f"    Skipping sequence that is too short: {len(seq)} < {self.region_end}")
                    else:
                        # Process the full sequence
                        alignment = pairwise2.align.globalxx(self.reference_sequence, seq, one_alignment_only=True)[0]
                        matches = sum(a == b for a, b in zip(alignment.seqA, alignment.seqB))
                        identity = (matches / len(alignment.seqA)) * 100
                        identities.append(identity)
                        valid_sequences.append(seq)
                except Exception as e:
                    print(f"    Error processing sequence: {e}")
            
            # Only continue if we have identities to analyze
            if not identities:
                print(f"  No valid sequences found for {target}% identity target.")
                continue
                
            all_results[target] = {
                "sequences": valid_sequences,
                "identities": identities,
                "mean_identity": np.mean(identities),
                "std_identity": np.std(identities),
                "min_identity": min(identities),
                "max_identity": max(identities)
            }
            
            print(f"  Mean identity: {np.mean(identities):.2f}% (Target: {target}%)")
            print(f"  Standard deviation: {np.std(identities):.2f}%")
            print(f"  Min identity: {min(identities):.2f}%")
            print(f"  Max identity: {max(identities):.2f}%")
            
            # Save sequences with their actual identities
            fasta_file = os.path.join(self.specific_output_dir, f"sequences_{target}percent_target.fasta")
            with open(fasta_file, "w") as f:
                for i, (seq, identity) in enumerate(zip(valid_sequences, identities)):
                    f.write(f">Seq_{target}percent_target_{i+1}_actual_{identity:.2f}\n")
                    f.write(f"{seq}\n")
        
        # Plot results if we have data
        if all_results:
            self._plot_identity_distributions(all_results, target_identities)
        
        return all_results
    
    def _plot_identity_distributions(self, all_results, target_identities):
        """Plot the identity distributions for the given results."""
        plt.figure(figsize=(10, 6))
        
        for target in sorted(target_identities):
            if target in all_results:
                result = all_results[target]
                plt.hist(result["identities"], bins=20, alpha=0.5, label=f"{target}% identity")
        
        title_suffix = f" for {self.analysis_name}" if self.analysis_name != "full_sequence" else ""
        plt.title(f"Identity Distributions{title_suffix}")
        plt.xlabel("Identity (%)")
        plt.ylabel("Frequency")
        plt.legend()
        plt.grid(True)
        
        plot_file = os.path.join(self.specific_output_dir, "identity_distributions.png")
        plt.savefig(plot_file)
        print(f"  Identity distributions saved to {plot_file}")
        
        plt.close()
        
        # Create an additional plot comparing mean identities across targets
        plt.figure(figsize=(10, 6))
        
        targets = []
        means = []
        stds = []
        
        for target in sorted(target_identities):
            if target in all_results:
                targets.append(target)
                means.append(all_results[target]["mean_identity"])
                stds.append(all_results[target]["std_identity"])
        
        plt.errorbar(targets, means, yerr=stds, fmt='o-', capsize=5)
        plt.title(f"Mean Identity by Target{title_suffix}")
        plt.xlabel("Target Identity (%)")
        plt.ylabel("Actual Mean Identity (%)")
        plt.grid(True)
        
        # Add a diagonal line representing perfect identity match
        min_target = min(targets) if targets else 0
        max_target = max(targets) if targets else 100
        plt.plot([min_target, max_target], [min_target, max_target], 'k--', alpha=0.5, label="Ideal")
        plt.legend()
        
        plot_file = os.path.join(self.specific_output_dir, "mean_identity_by_target.png")
        plt.savefig(plot_file)
        print(f"  Mean identity comparison saved to {plot_file}")
        
        plt.close()

    def generate_summary_report(self, all_results, target_identities):
        """Generate a summary report of the analysis results."""
        print("Generating summary report...")
        
        report_file = os.path.join(self.specific_output_dir, "analysis_summary.txt")
        
        with open(report_file, "w") as f:
            f.write("MPNN SEQUENCE ANALYSIS SUMMARY\n")
            f.write("============================\n\n")
            
            f.write(f"Date of analysis: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            
            if self.region_start is not None and self.region_end is not None:
                f.write(f"Analyzing region: {self.region_start}-{self.region_end}\n")
                region_len = self.region_end - self.region_start + 1
                f.write(f"Region length: {region_len} amino acids\n")
                f.write(f"Full reference sequence length: {len(self.reference_sequence)} amino acids\n\n")
            else:
                f.write(f"Analyzing full sequence\n")
                f.write(f"Reference sequence length: {len(self.reference_sequence)} amino acids\n\n")
            
            f.write("RESULTS BY TARGET IDENTITY\n")
            f.write("-------------------------\n\n")
            
            # Summary table header
            f.write(f"{'Target':>10} | {'Mean':>10} | {'Std Dev':>10} | {'Min':>10} | {'Max':>10} | {'Count':>10}\n")
            f.write(f"{'-'*10:>10} | {'-'*10:>10} | {'-'*10:>10} | {'-'*10:>10} | {'-'*10:>10} | {'-'*10:>10}\n")
            
            # Add row for each target identity
            for target in sorted(target_identities):
                if target in all_results:
                    result = all_results[target]
                    identities = result["identities"]
                    
                    f.write(f"{target:>10}% | {result['mean_identity']:>10.2f}% | {result['std_identity']:>10.2f}% | ")
                    f.write(f"{min(identities):>10.2f}% | {max(identities):>10.2f}% | {len(identities):>10}\n")
                else:
                    f.write(f"{target:>10}% | {'N/A':>10} | {'N/A':>10} | {'N/A':>10} | {'N/A':>10} | {'N/A':>10}\n")
            
            f.write("\n\nDETAILED ANALYSIS\n")
            f.write("----------------\n\n")
            
            for target in sorted(target_identities):
                if target in all_results:
                    result = all_results[target]
                    identities = result["identities"]
                    
                    f.write(f"Target {target}% Identity:\n")
                    f.write(f"  - Mean actual identity: {result['mean_identity']:.2f}%\n")
                    f.write(f"  - Standard deviation: {result['std_identity']:.2f}%\n")
                    f.write(f"  - Minimum identity: {min(identities):.2f}%\n")
                    f.write(f"  - Maximum identity: {max(identities):.2f}%\n")
                    f.write(f"  - Number of sequences: {len(identities)}\n")
                    
                    # Add histogram representation as ASCII art
                    bins = np.linspace(min(identities), max(identities), 20)
                    hist, _ = np.histogram(identities, bins=bins)
                    max_count = max(hist) if hist.size > 0 else 0
                    
                    if max_count > 0:  # Guard against division by zero
                        f.write("\n  Identity distribution (ASCII histogram):\n")
                        f.write(f"    {min(identities):.1f}% {'-'*30} {max(identities):.1f}%\n")
                        
                        for count in hist:
                            bar_len = int((count / max_count) * 30)
                            f.write(f"    {'|' * bar_len}\n")
                    
                    f.write("\n  Output files:\n")
                    f.write(f"    - Analyzed sequences: {os.path.join(self.specific_output_dir, f'sequences_{target}percent_target.fasta')}\n")
                    
                    f.write("\n")
                else:
                    f.write(f"Target {target}% Identity: No results found\n\n")
            
            f.write("\nANALYSIS NOTES\n")
            f.write("-------------\n\n")
            f.write("The identity percentage is calculated as the number of matching residues\n")
            f.write("divided by the length of the alignment between the reference sequence and\n")
            f.write("the generated sequence, multiplied by 100.\n\n")
            
            f.write("Visualizations saved:\n")
            f.write(f"1. Identity distributions: {os.path.join(self.specific_output_dir, 'identity_distributions.png')}\n")
            f.write(f"2. Mean identity by target: {os.path.join(self.specific_output_dir, 'mean_identity_by_target.png')}\n")
        
        print(f"Summary report saved to {report_file}")
        return report_file


# Execution script
def run_analysis():
    # Define your reference sequence
    reference_sequence = "MQSWYLLYCKRGQLQRAQEHLERQAVNCLAPMITLEKIVRGKRTAVSEPLFPNYLFVEFDPEVIHTTTINATRGVSHFVRFGASPAIVPSAVIHQLSVYKPKDIVDPATPYPGDKVIITEGAFEGFQAIFTEPDGEARSMLLLNLINKEIKHSVKNTEFRKL"
    
    # Set up output directory
    output_dir = "/home/cfneira1/ProteinMPNN/Seqs/identity_positions_analysis"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Define the MPNN results (paths to directories containing the FASTA files)
    mpnn_results = {
        57: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/identity_57/seqs/",
        62: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/identity_62/seqs/",
        67: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/identity_67/seqs/",
        72: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/identity_72/seqs/",
        77: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/identity_77/seqs/",
        82: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/identity_82/seqs/",
        87: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/identity_87/seqs/",
        92: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/identity_92/seqs/",
        97: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/identity_97/seqs/",
    }
    
    # Define target identities to check
    target_identities = [57, 62, 67, 72, 77, 82, 87, 92, 97]
    
    # Run analysis for full sequence
    full_checker = MPNNResultsChecker(
        reference_sequence=reference_sequence,
        output_dir=output_dir
    )
    full_results = full_checker.check_results(mpnn_results, target_identities)
    
    # Run analysis for region 115-162
    region_checker = MPNNResultsChecker(
        reference_sequence=reference_sequence,
        output_dir=output_dir,
        region_start=115,
        region_end=162
    )
    region_results = region_checker.check_results(mpnn_results, target_identities)
    
    # Print overall summary
    print("\n\n===== ANALYSIS COMPLETE =====")
    print("\nSummary of Results for Full Sequence:")
    for target, data in full_results.items():
        print(f"Target {target}% Identity:")
        print(f"  - Mean actual identity: {data['mean_identity']:.2f}%")
        print(f"  - Standard deviation: {data['std_identity']:.2f}%")
        print(f"  - Number of sequences: {len(data['sequences'])}")
    
    print("\nSummary of Results for Region 115-162:")
    for target, data in region_results.items():
        print(f"Target {target}% Identity:")
        print(f"  - Mean actual identity: {data['mean_identity']:.2f}%")
        print(f"  - Standard deviation: {data['std_identity']:.2f}%")
        print(f"  - Number of sequences: {len(data['sequences'])}")
    
    # Create a comparison plot between full sequence and region identities
    plt.figure(figsize=(12, 8))
    
    targets_full = []
    means_full = []
    
    targets_region = []
    means_region = []
    
    for target in sorted(target_identities):
        if target in full_results:
            targets_full.append(target)
            means_full.append(full_results[target]["mean_identity"])
        
        if target in region_results:
            targets_region.append(target)
            means_region.append(region_results[target]["mean_identity"])
    
    plt.plot(targets_full, means_full, 'o-', label="Full Sequence", color="blue")
    plt.plot(targets_region, means_region, 's-', label="Region 115-162", color="red")
    
    # Add a diagonal line representing perfect identity match
    all_targets = sorted(set(targets_full + targets_region))
    if all_targets:
        min_target = min(all_targets)
        max_target = max(all_targets)
        plt.plot([min_target, max_target], [min_target, max_target], 'k--', alpha=0.5, label="Ideal")
    
    plt.title("Mean Identity Comparison: Full Sequence vs Region 115-162")
    plt.xlabel("Target Identity (%)")
    plt.ylabel("Actual Mean Identity (%)")
    plt.legend()
    plt.grid(True)
    
    comparison_file = os.path.join(output_dir, "full_vs_region_comparison.png")
    plt.savefig(comparison_file)
    print(f"\nComparison plot saved to {comparison_file}")
    
    plt.close()

if __name__ == "__main__":
    run_analysis()