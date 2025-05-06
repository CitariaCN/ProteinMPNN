import os
import numpy as np
from Bio import Align
from Bio import pairwise2
from Bio import SeqIO
from Bio.SeqIO import parse
from Bio.SeqIO import write
from Bio.SeqIO import FastaIO
import datetime
from Bio.Seq import Seq

class MPNNResultsChecker:
    def __init__(self, reference_sequence, output_dir):
        self.reference_sequence = reference_sequence
        self.output_dir = output_dir

    def check_results(self, mpnn_results, target_identities):
        """Check the results of MPNN runs for the specified target identities."""
        print("Checking MPNN results...")
        
        # Verify the identity percentages
        all_results = self.verify_identity_percentages(mpnn_results, target_identities)
        
        # Check if all targets were processed
        for target in target_identities:
            if target not in all_results:
                print(f"  No results found for {target}% identity target.")
        
        # Generate summary report
        if all_results:
            self.generate_summary_report(all_results, target_identities)
        
        print("MPNN results check completed.")
        return all_results
    
    def verify_identity_percentages(self, mpnn_results, target_identities):
        """Verify the actual identity percentages of the generated sequences"""
        print("Verifying sequence identities...")
        
        all_results = {}
        
        for target in target_identities:
            output_dir = mpnn_results.get(target)
            if not output_dir or not os.path.exists(output_dir):
                print(f"  No results directory found for {target}% identity target.")
                continue
                
            distance_threshold = [30,60,90]
            all_sequences = [] 

            for threshold in distance_threshold:
                fasta_file = os.path.join(output_dir, f"sequences_{target}percent_target_{threshold}.fasta")
                if not os.path.exists(fasta_file):
                    print(f"  No FASTA file found for {target}% identity target at threshold {threshold}.")
                    continue
                
                print(f"  Processing threshold {threshold} for {target}% identity target...")
                try:
                    with open(fasta_file, 'r') as f:
                        sequences_in_file = []
                        for line in f:
                            if line.startswith(">"):
                                continue  # Skip header lines in FASTA format
                    
                            sequences_in_file.append(line.strip())  # Add sequence lines to the list
                        print(f"  Found {len(sequences_in_file)} sequences for {target}% identity target at threshold {threshold}.")
                
                except Exception as e:
                    print(f"  Error processing sequences for {target}% identity: {e}")
                    continue
            
            if not sequences_in_file:
                print(f"  No sequences found for {target}% identity target.")
                continue
            
            print(f"  Found {len(sequences_in_file)} sequences for {target}% identity target.")

            #Calculate identities
            identities = []
            for seq in all_sequences:
                    alignment = pairwise2.align.globalxx(self.reference_sequence, seq, one_alignment_only=True)[0]
                    matches = sum(a == b for a, b in zip(alignment.seqA, alignment.seqB))
                    identity = (matches / len(alignment.seqA)) * 100
                    identities.append(identity)
                
            all_results[target] = {
                    "sequences": all_sequences,
                    "identities": identities,
                    "mean_identity": np.mean(identities),
                    "std_identity": np.std(identities)
                }
            
            print(f"  Mean identity for {target}% target: {all_results[target]['mean_identity']:.2f}%")
            print(f"  Standard deviation for {target}% target: {all_results[target]['std_identity']:.2f}%")

            fasta_file = os.path.join(self.output_dir, f"sequences_{target}percent_target.fasta")
            with open(fasta_file, 'w') as f:
                for i, (seq, identity) in enumerate(zip(all_sequences, identities)):
                    f.write(f">Sequence_{i+1} Identity: {identity:.2f}%\n")
                    f.write(f"{seq}\n")
        
        if all_results:
            # Plot the identity distributions
            self._plot_identity_distributions(all_results, target_identities)
        
        return all_results
            

    
    def _plot_identity_distributions(self, all_results, target_identities):
        """Plot the identity distributions for the given results."""
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        
        for target, result in all_results.items():
            plt.hist(result["identities"], bins=20, alpha=0.5, label=f"{target}% identity")
        
        plt.title("Identity Distributions")
        plt.xlabel("Identity (%)")
        plt.ylabel("Frequency")
        plt.legend()
        plt.grid()
        
        plot_file = os.path.join(self.output_dir, "identity_distributions.png")
        plt.savefig(plot_file)
        print(f"  Identity distributions saved to {plot_file}")
        
        plt.close()

    def generate_summary_report(self, all_results, target_identities):
        """Generate a summary report of the analysis results."""
        print("Generating summary report...")
    
        report_file = os.path.join(self.output_dir, "analysis_summary.txt")
    
        with open(report_file, "w") as f:
            f.write("MPNN SEQUENCE ANALYSIS SUMMARY\n")
            f.write("============================\n\n")
        
            f.write(f"Date of analysis: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
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
                    max_count = max(hist)
                
                    f.write("\n  Identity distribution (ASCII histogram):\n")
                    f.write(f"    {min(identities):.1f}% {'-'*30} {max(identities):.1f}%\n")
                
                    for count in hist:
                        bar_len = int((count / max_count) * 30)
                        f.write(f"    {'|' * bar_len}\n")
                
                    f.write("\n  Output files:\n")
                    f.write(f"    - Analyzed sequences: {os.path.join(self.output_dir, f'sequences_{target}percent_target.fasta')}\n")
                
                    f.write("\n")
                else:
                    f.write(f"Target {target}% Identity: No results found\n\n")
        
            f.write("\nANALYSIS NOTES\n")
            f.write("-------------\n\n")
            f.write("The identity percentage is calculated as the number of matching residues\n")
            f.write("divided by the length of the alignment between the reference sequence and\n")
            f.write("the generated sequence, multiplied by 100.\n\n")
        
            f.write("A visualization of the identity distributions has been saved as:\n")
            f.write(f"{os.path.join(self.output_dir, 'identity_distributions.png')}\n")
    
        print(f"Summary report saved to {report_file}")
        return report_file