import os
import random
import urllib.request
from collections import defaultdict, Counter

NCBI_ACCESSION = "NC_000913.3"
GENOME_FASTA = "genome_reference.fasta"
SAMPLES_FASTA = "sampled_reads_20k.fasta"
OUTPUT_RECONSTRUCTED = "reconstructed_genome.fasta"

NUM_SAMPLES = 20000
MIN_READ_LEN = 100
MAX_READ_LEN = 150
K_MER_SIZE = 31  


def generate_samples_if_missing():
    """Descarcă genomul și generează cele 20.000 de eșantioane dacă nu există deja."""
    if os.path.exists(SAMPLES_FASTA):
        print(f"-> Găsit fișierul existent '{SAMPLES_FASTA}'. Se trece la asamblare.")
        return

    print(f"-> Fișierul '{SAMPLES_FASTA}' nu a fost găsit. Se generează eșantioanele...")

    if not os.path.exists(GENOME_FASTA):
        print(f"Downloading genom {NCBI_ACCESSION} de pe NCBI...")
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={NCBI_ACCESSION}&rettype=fasta&retmode=text"
        urllib.request.urlretrieve(url, GENOME_FASTA)

    seq_chunks = []
    with open(GENOME_FASTA, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq_chunks.append(line.strip().upper())
    
    full_genome = "".join(seq_chunks)
    genome_len = len(full_genome)
    print(f"Genom încărcat cu succes (Lungime totală: {genome_len:,} bp).")

    print(f"Se extrag {NUM_SAMPLES:,} eșantioane aleatoare...")
    samples = []
    for i in range(1, NUM_SAMPLES + 1):
        sample_len = random.randint(MIN_READ_LEN, MAX_READ_LEN)
        max_start = genome_len - sample_len
        start_idx = random.randint(0, max_start)
        samples.append((i, start_idx, sample_len, full_genome[start_idx : start_idx + sample_len]))

    # 4. Salvează eșantioanele în fișier FASTA
    with open(SAMPLES_FASTA, "w") as out:
        for sample_id, start_pos, length, seq in samples:
            out.write(f">sample_{sample_id}|start={start_pos}|len={length}\n{seq}\n")

    print(f" Generat cu succes fișierul '{SAMPLES_FASTA}'.\n")


def load_reads_from_fasta(filepath: str) -> list:
    """Încarcă secvențele din fișierul FASTA."""
    reads = []
    with open(filepath, "r") as file:
        current_seq = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    reads.append("".join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line.upper())
        if current_seq:
            reads.append("".join(current_seq))

    print(f"S-au încărcat {len(reads):,} secvențe (reads).")
    return reads


def build_de_bruijn_graph(reads: list, k: int):
    """Construiește graful De Bruijn pe baza k-merilor."""
    edges = defaultdict(Counter)
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)

    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i : i + k]
            left_node = kmer[:-1]
            right_node = kmer[1:]

            edges[left_node][right_node] += 1
            out_degree[left_node] += 1
            in_degree[right_node] += 1

    print(f"Graf De Bruijn construit: {len(edges):,} noduri unice (k-1-meri).")
    return edges, in_degree, out_degree


def assemble_contigs(edges, in_degree, out_degree, k: int) -> list:
    """Traversează căile din graf pentru a reconstrui contigurile genomului."""
    visited_edges = set()
    contigs = []

    start_nodes = [node for node in edges if in_degree[node] < out_degree[node]]
    if not start_nodes:
        start_nodes = list(edges.keys())

    for start in start_nodes:
        curr = start
        while True:
            next_candidates = [
                (target, count)
                for target, count in edges[curr].items()
                if (curr, target) not in visited_edges
            ]

            if not next_candidates:
                break

            next_candidates.sort(key=lambda x: x[1], reverse=True)
            next_node = next_candidates[0][0]

            visited_edges.add((curr, next_node))
            contig = curr
            curr_step = next_node

            while True:
                contig += curr_step[-1]
                candidates = [
                    (t, c)
                    for t, c in edges[curr_step].items()
                    if (curr_step, t) not in visited_edges
                ]

                if not candidates or len(candidates) > 1 or in_degree[curr_step] > 1:
                    break

                candidates.sort(key=lambda x: x[1], reverse=True)
                next_step = candidates[0][0]
                visited_edges.add((curr_step, next_step))
                curr_step = next_step

            if len(contig) >= k:
                contigs.append(contig)

    contigs.sort(key=len, reverse=True)
    return contigs


if __name__ == "__main__":
    print("=== DE NOVO GENOME ASSEMBLER (20,000 SAMPLES) ===")

    generate_samples_if_missing()

    reads = load_reads_from_fasta(SAMPLES_FASTA)

    edges, in_deg, out_deg = build_de_bruijn_graph(reads, k=K_MER_SIZE)

    print("Se asamblează genomul de novo...")
    contigs = assemble_contigs(edges, in_deg, out_deg, k=K_MER_SIZE)

    if contigs:
        primary_contig = contigs[0]
        print("\n=== REZULTAT FINAL ===")
        print(f"Număr total contiguri generate: {len(contigs)}")
        print(f"Lungimea celui mai lung contig reconstruit: {len(primary_contig):,} bp")
        print(f"Primele 100 nucleotide: {primary_contig[:100]}...")

        with open(OUTPUT_RECONSTRUCTED, "w") as f:
            f.write(f">reconstructed_contig_1|len={len(primary_contig)}\n{primary_contig}\n")
        print(f"Rezultat salvat în '{OUTPUT_RECONSTRUCTED}'.")
    else:
        print("Nu s-au putut genera contiguri. Încearcă să micșorezi K_MER_SIZE (ex: k=21).")