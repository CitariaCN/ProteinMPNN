import os
import datetime
from check_identity_results import MPNNResultsChecker  # Assuming you saved the class in check_identity_results.py

# Define your reference sequence
reference_sequence = "MQSWYLLYCKRGQLQRAQEHLERQAVNCLAPMITLEKIVRGKRTAVSEPLFPNYLFVEFDPEVIHTTTINATRGVSHFVRFGASPAIVPSAVIHQLSVYKPKDIVDPATPYPGDKVIITEGAFEGFQAIFTEPDGEARSMLLLNLINKEIKHSVKNTEFRKL"

# Set up output directory
output_dir = "/home/cfneira1/ProteinMPNN/Seqs/identity_positions"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Initialize the checker
checker = MPNNResultsChecker(reference_sequence, output_dir)

# Define the MPNN results (paths to directories containing seqs.fa files)
mpnn_results = {
    57: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/seqs/identity_57",
    62: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/seqs/identity_62",
    67: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/seqs/identity_67",
    72: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/seqs/identity_72",
    77: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/seqs/identity_77",
    82: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/seqs/identity_82",
    87: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/seqs/identity_87",
    92: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/seqs/identity_92",
    97: "/home/cfneira1/ProteinMPNN/Seqs/identity_positions/seqs/identity_97",
}

# Define target identities to check
target_identities = [57, 62, 67, 72, 77, 82, 87, 92, 97]

# Run the checker
results = checker.check_results(mpnn_results, target_identities)

# Print a summary
print("\nSummary of Results:")
for target, data in results.items():
    print(f"Target {target}% Identity:")
    print(f"  - Mean actual identity: {data['mean_identity']:.2f}%")
    print(f"  - Standard deviation: {data['std_identity']:.2f}%")
    print(f"  - Number of sequences: {len(data['sequences'])}")