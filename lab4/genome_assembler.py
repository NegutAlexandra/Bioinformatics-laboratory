import random
import os
import urllib.request

def get_sample_genome(accession: str = "NC_000913.3", length_limit: int = 1500) -> str:
    """Descarcă un genom de pe NCBI și returnează o secvență de test."""
    fasta_file = "genome_test.fasta"
    if not os.path.exists(fasta_file):
        print(f"Se descarcă genomul {accession} de pe NCBI...")
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
        urllib.request.urlretrieve(url, fasta_file)

    sequence_chunks = []
    with open(fasta_file, "r") as f:
        for line in f:
            if not line.startswith(">"):
                sequence_chunks.append(line.strip().upper())
    
    full_genome = "".join(sequence_chunks)
    return full_genome[:length_limit]


def generate_overlapping_reads(genome: str, min_len: int = 100, max_len: int = 150, min_overlap: int = 10):
    """Sparge genomul în subsecvențe cu suprapunere garantată."""
    reads = []
    pos = 0
    genome_len = len(genome)

    while pos < genome_len:
        read_len = random.randint(min_len, max_len)
        end_pos = min(pos + read_len, genome_len)
        
        read = genome[pos:end_pos]
        reads.append(read)

        if end_pos == genome_len:
            break

        max_step = max(1, (end_pos - pos) - min_overlap)
        step = random.randint(1, max_step)
        pos += step

    print(f"-> Genomul de {genome_len} bp a fost spart în {len(reads)} subsecvențe.")
    return reads


def find_overlap(seq1: str, seq2: str, min_overlap: int = 10) -> int:
    """Calculează lungimea suprapunerii unde sufixul seq1 se potrivește cu prefixul seq2."""
    max_len = 0
    start = 0
    while True:
        start = seq1.find(seq2[:min_overlap], start)
        if start == -1:
            break
        if seq2.startswith(seq1[start:]):
            overlap_len = len(seq1) - start
            if overlap_len > max_len:
                max_len = overlap_len
            break
        start += 1
    return max_len


def print_overlap_alignment(seq1: str, seq2: str, overlap_len: int, step_num: int):
    """Afișează grafic modul în care cele două secvențe se suprapun."""
    offset = len(seq1) - overlap_len
    indent = " " * offset
    match_line = " " * offset + "|" * overlap_len
    
    print(f"\n--- PASUL DE UNIRE #{step_num} (Suprapunere: {overlap_len} bp) ---")
    print(f"Read 1:  {seq1}")
    print(f"Potriv:  {match_line}")
    print(f"Read 2:  {indent}{seq2}")
    print(f"Rezultat:{seq1 + seq2[overlap_len:]}")


def assemble_genome_with_visualization(reads: list, min_overlap: int = 10, max_visualizations: int = 3) -> str:
    """Reconstruiește genomul și afișează primele alinieri de suprapunere."""
    reads_pool = reads.copy()
    step = 0

    while True:
        max_ov = 0
        best_pair = (-1, -1)

        n = len(reads_pool)
        for i in range(n):
            for j in range(n):
                if i != j:
                    ov = find_overlap(reads_pool[i], reads_pool[j], min_overlap)
                    if ov > max_ov:
                        max_ov = ov
                        best_pair = (i, j)

        if max_ov < min_overlap:
            break

        i, j = best_pair
        step += 1

        if step <= max_visualizations:
            print_overlap_alignment(reads_pool[i], reads_pool[j], max_ov, step)

        merged_seq = reads_pool[i] + reads_pool[j][max_ov:]
        
        if i > j:
            reads_pool.pop(i)
            reads_pool.pop(j)
        else:
            reads_pool.pop(j)
            reads_pool.pop(i)
            
        reads_pool.append(merged_seq)

    reads_pool.sort(key=len, reverse=True)
    return reads_pool[0]


if __name__ == "__main__":
    a, b = 100, 150
    min_overlap = 10

    original_genome = get_sample_genome(length_limit=1000)
    print(f"Lungime genom original: {len(original_genome)} bp")

    reads = generate_overlapping_reads(original_genome, min_len=a, max_len=b, min_overlap=min_overlap)
    random.shuffle(reads)

    print("\nSe începe asamblarea genomului...")
    reconstructed_genome = assemble_genome_with_visualization(reads, min_overlap=min_overlap, max_visualizations=3)

    print("\n================ REZULTAT FINAL ================")
    print(f"Lungime genom reconstruit: {len(reconstructed_genome)} bp")
    if reconstructed_genome == original_genome:
        print(" SUCCESS: Genomul a fost RECONSTRUIT PERFECT!")
    else:
        print(" Partial/Mismatch: Genomul s-a asamblat diferit (pot apărea regiuni repetitive).")