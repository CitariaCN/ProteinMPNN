def verify_identity_percentages(mpnn_results, reference_sequence, target_identities):
    import json
    from Bio import pairwise2
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    
    all_results = {}
    
    for target in target_identities:
        output_dir = mpnn_results[target]
        jsonl_file = os.path.join(output_dir, "seqs.jsonl")  # Adjust path if needed
        
        sequences = []
        with open(jsonl_file, 'r') as f:
            for line in f:
                data = json.loads(line)
                for seq in data["sampled_sequences"]:
                    sequences.append(seq)
        
        # Calculate identities
        identities = []
        for seq in sequences:
            alignment = pairwise2.align.globalxx(reference_sequence, seq, one_alignment_only=True)[0]
            matches = sum(a == b for a, b in zip(alignment.seqA, alignment.seqB))
            identity = (matches / len(alignment.seqA)) * 100
            identities.append(identity)
        
        all_results[target] = {
            "sequences": sequences,
            "identities": identities,
            "mean_identity": np.mean(identities),
            "std_identity": np.std(identities)
        }
        
        # Save sequences with their actual identities
        with open(f"sequences_{target}percent_target.fasta", "w") as f:
            for i, (seq, identity) in enumerate(zip(sequences, identities)):
                f.write(f">Seq_{target}percent_target_{i+1}_actual_{identity:.2f}\n")
                f.write(f"{seq}\n")
    
    # Plot results
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for target in target_identities:
        identities = all_results[target]["identities"]
        ax.hist(identities, alpha=0.5, label=f"Target: {target}%")
    
    ax.set_xlabel("Sequence Identity (%)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of Sequence Identities")
    ax.legend()
    
    plt.savefig("identity_distributions.png")
    plt.close()
    
    return all_results

# Usage
identity_results = verify_identity_percentages(mpnn_results, reference_sequence, 
                                              [90, 85, 80, 75, 70, 65, 60, 55])