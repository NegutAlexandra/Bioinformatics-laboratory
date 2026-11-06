import urllib.request
from collections import defaultdict

def fetch_ncbi_genome(accession="NC_000908.2", output_file="genome.fasta") -> str:
    """
    Descarcă un genom complet de pe NCBI în format FASTA.
    """
    print(f"[+] Descărcăm genomul de pe NCBI (Accession: {accession})...")
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    
    try:
        urllib.request.urlretrieve(url, output_file)
        print(f"[+] Genom salvat cu succes în '{output_file}'!")
    except Exception as e:
        print(f"[-] Eroare la descărcare: {e}")
        return ""
    
    return read_fasta(output_file)

def read_fasta(file_path: str) -> str:
    """Citește fișierul FASTA și returnează secvența curățată."""
    sequence = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('>'):
                sequence.append(line.upper())
    return "".join(sequence)

COMPLEMENT = str.maketrans('ATCGN', 'TAGCN')

def get_reverse_complement(seq: str) -> str:
    """Calculează inversul complementat al unei secvențe ADN."""
    return seq.translate(COMPLEMENT)[::-1]

def scan_unknown_transposons(dna_seq: str, min_ir_len=10, max_ir_len=15, min_transposon_len=100, max_transposon_len=3000):
    """
    Detectează transposoni necunoscuți folosind două ferestre glisante:
    - S1 parcurge genomul cu pas 1.
    - S2 caută complementul invers (Inverted Repeat) la o distanță de CEL PUȚIN 100 de baze față de S1.
    """
    seq_len = len(dna_seq)
    detected_transposons = []
    
    print(f"\n[+] Începem scanarea genomului (Lungime totală: {seq_len} bp)...")
    print(f"[+] Dimensiune Inverted Repeat: {min_ir_len}-{max_ir_len} bp")
    print(f"[+] Distanță minimă între S1 și S2: {min_transposon_len} bp (Maxim: {max_transposon_len} bp)\n")

    kmer_positions = defaultdict(list)
    for i in range(seq_len - min_ir_len + 1):
        kmer = dna_seq[i:i + min_ir_len]
        if 'N' not in kmer:
            kmer_positions[kmer].append(i)

    scanned_count = 0
    for s1_pos in range(0, seq_len - min_transposon_len):
        scanned_count += 1
        if scanned_count % 100000 == 0:
            print(f"   --> Procesat {scanned_count} / {seq_len} bp...")

        s1_pattern = dna_seq[s1_pos : s1_pos + min_ir_len]
        if 'N' in s1_pattern:
            continue
            
        rev_comp_pattern = get_reverse_complement(s1_pattern)
        max_search_pos = min(s1_pos + max_transposon_len, seq_len - min_ir_len)
        
        if rev_comp_pattern in kmer_positions:
            for s2_pos in kmer_positions[rev_comp_pattern]:
                if s1_pos + min_transposon_len <= s2_pos <= max_search_pos:
                    
                    ir_len = min_ir_len
                    while (s1_pos + ir_len < s2_pos) and (s2_pos + ir_len < seq_len) and (ir_len < max_ir_len):
                        if dna_seq[s1_pos + ir_len] == get_reverse_complement(dna_seq[s2_pos + ir_len - 1]):
                            ir_len += 1
                        else:
                            break
                    
                    left_ir = dna_seq[s1_pos : s1_pos + ir_len]
                    right_ir = dna_seq[s2_pos : s2_pos + ir_len]
                    
                    detected_transposons.append({
                        'start': s1_pos,
                        'end': s2_pos + ir_len,
                        'length': (s2_pos + ir_len) - s1_pos,
                        'left_ir': left_ir,
                        'right_ir': right_ir
                    })

    unique_transposons = []
    for t in detected_transposons:
        if not any(abs(t['start'] - u['start']) < 10 and abs(t['end'] - u['end']) < 10 for u in unique_transposons):
            unique_transposons.append(t)

    return unique_transposons

if __name__ == "__main__":
    fasta_file = "genome_nc_000908.fasta"
    dna_genome = fetch_ncbi_genome(accession="NC_000908.2", output_file=fasta_file)
    
    if dna_genome:
        sample_region = dna_genome[:50000]
        
        results = scan_unknown_transposons(
            sample_region, 
            min_ir_len=10, 
            max_ir_len=15, 
            min_transposon_len=100,
            max_transposon_len=2000
        )
        
        print("=" * 80)
        print(f"CĂUTARE DE NOVO FINALIZATĂ — S-AU GĂSIT {len(results)} POTENȚIALI TRANSPOZONI (MIN 100 BP)")
        print("=" * 80)
        print(f"{'No.':<4} | {'Start':<7} | {'End':<7} | {'Length':<8} | {'Left IR (S1)':<15} | {'Right IR (S2)'}")
        print("-" * 80)
        
        for idx, res in enumerate(results[:15], 1):
            print(f"{idx:<4} | {res['start']:<7} | {res['end']:<7} | {res['length']:<8} | {res['left_ir']:<15} | {res['right_ir']}")
        print("=" * 80)