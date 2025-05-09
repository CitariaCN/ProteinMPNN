import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def analyze_sequence_identities(csv_file):
    # Read the CSV file
    df = pd.read_csv(csv_file)
    
    # Convert percent identity to numeric values
    df['Per. ident'] = pd.to_numeric(df['Per. ident'], errors='coerce')
    
    # Count the total number of sequences
    total_sequences = len(df)
    
    # Create a value counts series of the exact percent identities
    identity_counts = df['Per. ident'].value_counts().sort_index(ascending=False)
    
    # Calculate the percentage of sequences with each identity value
    identity_percentages = (identity_counts / total_sequences) * 100
    
    # Calculate cumulative percentages (what percentage of sequences have at least X% identity)
    cumulative_percentages = identity_counts.sort_index().cumsum() / total_sequences * 100
    cumulative_percentages = cumulative_percentages.sort_index(ascending=False)
    
    # Create a detailed histogram of percent identities
    plt.figure(figsize=(12, 6))
    sns.histplot(df['Per. ident'], bins=50, kde=True)
    plt.title('Distribution of Percent Identities')
    plt.xlabel('Percent Identity')
    plt.ylabel('Count')
    plt.savefig('identity_distribution_detailed.png')
    
    # Create a cumulative distribution plot
    plt.figure(figsize=(12, 6))
    plt.plot(cumulative_percentages.index, cumulative_percentages.values)
    plt.title('Cumulative Distribution of Percent Identities')
    plt.xlabel('Percent Identity')
    plt.ylabel('Percentage of Sequences with ≥ X% Identity')
    plt.grid(True)
    plt.savefig('cumulative_distribution.png')
    
    return df, identity_counts, identity_percentages, cumulative_percentages

def main():
    csv_file = '/home/cfneira1/ProteinMPNN/1W6HPFAE013-Alignment-Descriptions.csv'  # Replace with your CSV file path
    df, identity_counts, identity_percentages, cumulative_percentages = analyze_sequence_identities(csv_file)
    
    # Print results
    print(f"Total sequences analyzed: {len(df)}")
    
    print("\nExact distribution of percent identities:")
    print("=======================================")
    print("Percent Identity | Count | Percentage of Total")
    print("-----------------------------------------------")
    for identity, count in identity_counts.items():
        percentage = identity_percentages[identity]
        print(f"{identity:.2f}% | {count} | {percentage:.2f}%")
    
    print("\nCumulative distribution (percentage of sequences with at least X% identity):")
    print("======================================================================")
    print("Percent Identity | Percentage of Sequences")
    print("----------------------------------------")
    for identity, percentage in cumulative_percentages.items():
        print(f"≥ {identity:.2f}% | {percentage:.2f}%")
    
    # Additional analysis: Calculate statistics by taxonomic group
    if 'Scientific Name' in df.columns:
        print("\nPercent identity statistics by organism:")
        print("=====================================")
        species_stats = df.groupby('Scientific Name')['Per. ident'].agg(['count', 'mean', 'min', 'max'])
        species_stats = species_stats.sort_values('count', ascending=False)
        
        for species, stats in species_stats.iterrows():
            print(f"{species}: {stats['count']} sequences, Mean: {stats['mean']:.2f}%, Range: {stats['min']:.2f}%-{stats['max']:.2f}%")
    
if __name__ == "__main__":
    main()