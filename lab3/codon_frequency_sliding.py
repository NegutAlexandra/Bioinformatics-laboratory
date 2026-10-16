import collections

CODON_TABLE = {
    'AUG': 'Methionine (M)', 'UUU': 'Phenylalanine (F)', 'UUC': 'Phenylalanine (F)',
    'UUA': 'Leucine (L)', 'UUG': 'Leucine (L)', 'UCU': 'Serine (S)',
    'UCC': 'Serine (S)', 'UCA': 'Serine (S)', 'UCG': 'Serine (S)',
    'UAU': 'Tyrosine (Y)', 'UAC': 'Tyrosine (Y)', 'UAA': 'STOP',
    'UAG': 'STOP', 'UGA': 'STOP', 'UGU': 'Cysteine (C)',
    'UGC': 'Cysteine (C)', 'UGG': 'Tryptophan (W)', 'CUU': 'Leucine (L)',
    'CUC': 'Leucine (L)', 'CUA': 'Leucine (L)', 'CUG': 'Leucine (L)',
    'CCU': 'Proline (P)', 'CCC': 'Proline (P)', 'CCA': 'Proline (P)',
    'CCG': 'Proline (P)', 'CAU': 'Histidine (H)', 'CAC': 'Histidine (H)',
    'CAA': 'Glutamine (Q)', 'CAG': 'Glutamine (Q)', 'CGU': 'Arginine (R)',
    'CGC': 'Arginine (R)', 'CGA': 'Arginine (R)', 'CGG': 'Arginine (R)',
    'AUU': 'Isoleucine (I)', 'AUC': 'Isoleucine (I)', 'AUA': 'Isoleucine (I)',
    'ACU': 'Threonine (T)', 'ACC': 'Threonine (T)', 'ACA': 'Threonine (T)',
    'ACG': 'Threonine (T)', 'AAU': 'Asparagine (N)', 'AAC': 'Asparagine (N)',
    'AAA': 'Lysine (K)', 'AAG': 'Lysine (K)', 'AGU': 'Serine (S)',
    'AGC': 'Serine (S)', 'AGA': 'Arginine (R)', 'AGG': 'Arginine (R)',
    'GUU': 'Valine (V)', 'GUC': 'Valine (V)', 'GUA': 'Valine (V)',
    'GUG': 'Valine (V)', 'GCU': 'Alanine (A)', 'GCC': 'Alanine (A)',
    'GCA': 'Alanine (A)', 'GCG': 'Alanine (A)', 'GAU': 'Aspartic acid (D)',
    'GAC': 'Aspartic acid (D)', 'GAA': 'Glutamic acid (E)', 'GAG': 'Glutamic acid (E)',
    'GGU': 'Glycine (G)', 'GGC': 'Glycine (G)', 'GGA': 'Glycine (G)',
    'GGG': 'Glycine (G)'
}

SAMPLE_VIRAL_GENOME = """
>sample_viral_genome_sequence
ATGGCCATTGTAATGGGCCGCTGAAAGGGTCCCGGGAAATTTGGCCATGAATAGCGATCG
ATCGATCGATCGATCGCATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
ATGGCCATTGTAATGGGCCGCTGAAAGGGTCCCGGGAAATTTGGCCATGAATAGCGATCG
GAAAGAAGAAGAAGAAGAAATGGCCATTGTAATGGGCCGCTGAAAGGGTCCCGGGAAATT
"""

def get_clean_sequence(fasta_data: str) -> str:
    """Parses FASTA formatted string data and returns a clean sequence."""
    sequence = []
    for line in fasta_data.strip().split('\n'):
        line = line.strip()
        if not line.startswith('>'):
            sequence.append(line.upper())
    return "".join(sequence)

def analyze_overlapping_codons(dna_sequence: str, source_label: str = "Embedded Sequence"):
    rna_sequence = dna_sequence.replace('T', 'U')
    
    codon_counts = collections.Counter()
    
    window_size = 3
    for i in range(len(rna_sequence) - window_size + 1):
        codon = rna_sequence[i:i + window_size]
        if all(base in 'ACGU' for base in codon):
            codon_counts[codon] += 1

    total_windows = sum(codon_counts.values())
    top_10 = codon_counts.most_common(10)

    print("=" * 68)
    print(f"SLIDING WINDOW CODON ANALYSIS ({source_label})")
    print(f"Total 3-bp windows analyzed: {total_windows}")
    print("=" * 68)
    print(f"{'Rank':<5} | {'Codon (RNA)':<11} | {'Frequency':<10} | {'Percentage':<10} | {'Amino Acid'}")
    print("-" * 68)

    for rank, (codon, count) in enumerate(top_10, 1):
        percentage = (count / total_windows) * 100 if total_windows > 0 else 0
        amino_acid = CODON_TABLE.get(codon, 'Unknown')
        print(f"{rank:<5} | {codon:<11} | {count:<10} | {percentage:>6.2f}%    | {amino_acid}")
    print("=" * 68)

if __name__ == "__main__":
    clean_dna = get_clean_sequence(SAMPLE_VIRAL_GENOME)
    
    analyze_overlapping_codons(clean_dna, source_label="Sample Genome")