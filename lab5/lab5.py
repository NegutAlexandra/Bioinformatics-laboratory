import random
import matplotlib.pyplot as plt
import re
from collections import defaultdict

NCBI_DNA_SEQUENCE = (
    "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCT"
    "GCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCAACTCAG"
    "GACTTGTTCTTACCTTTCTTTTCCAATGTTACTTGGTTCCATGCTATACATGTCTCTGGGACCAATGGTACTAAGAGG"
    "TTTGATAACCCTGTCCTACCATTTAATGATGGTGTTTATTTTGCTTCCACTGAGAAGTCTAACATAATAAGAGGCTGG"
    "ATCTTTGGTACTACTTTAGATTCGAAGACCCAGTCCCTACTTATTGTTAATAACGCTACTAATGTTGTTATTAAAGTC"
    "TGTGAATTTCAATTTTGTAATGATCCATTTTTGGGTGTTTATTACCACAAAAACAACAAAAGTTGGATGGAAAGTGAG"
    "TTCAGAGTTTATTCTAGTGCGAATAATTGCACTTTTGAATATGTCTCTCAGCCTTTTCTTATGGACCTTGAAGGAAAA"
    "CAGGGTAATTTCAAAAATCTTAGGGAATTTGTGTTTAAGAATATTGATGGTTATTTCAAAATATATTCTAAGCACACC"
    "CCTATTAATCTTGTTCGCGATCTCCCTCAGGGTTTTCGGCTTTAGAACCATTGGTAGATTTGCCAATAGGTATTAAC"
    "ATCACTAGGTTTCAAACTCTACTTGCCCTACATAGAAGTTATTTGACTCCTGGTGATTCTTCTTCAGGTTGGACAGCT"
    "GGTGCTGCAGCTTATTATGTGGGTTATCTTCAACCTAGGACTTTTCTATTAAAATATAATGAAAATGGAACCATTACA"
    "GATGCTGTAGACTGTGCACTTGACCCTCTCTCAGAAACAAAGTGTACGTTGAAATCCTTCACTGTAGAAAAAGGAATC"
    "TATCAAACTTCTAACTTTAGAGTCCAACCAACAGAATCTATTGTTAGATTTCCTAATATTACAAACTTGTGCCCTTTT"
    "GGTGAAGTTTTTAACGCCACCAGATTTGCATCTGTTTATGCTTGGAACAGGAAGAGAATCAGCAACTGTGTTGCTGAC"
    "TATTCTGTCCTATATAATTCCACATCATTTTCCACTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTC"
)

def simulate_gel_electrophoresis(dna_seq: str, num_samples: int = 10):
    """
    Simulates DNA gel electrophoresis:
    - Takes random DNA fragments between 100 and 3000 bp.
    - Simulates fragment migration: Smaller fragments move faster/further.
    """
    seq_len = len(dna_seq)
    samples = []

    for _ in range(num_samples):
        frag_len = random.randint(100, min(3000, seq_len))
        start_idx = random.randint(0, max(0, seq_len - frag_len))
        fragment = dna_seq[start_idx : start_idx + frag_len]
        samples.append(fragment)

    lengths = [len(s) for s in samples]

    print("=" * 60)
    print("EXERCISE 1: GEL ELECTROPHORESIS SIMULATION")
    print("=" * 60)
    print(f"Generated {num_samples} DNA fragments:")
    for i, length in enumerate(lengths, 1):
        print(f"  Lane {i}: {length} bp")

    fig, ax = plt.subplots(figsize=(6, 8))
    ax.set_facecolor("black")

    num_lanes = len(lengths)
    ax.set_xlim(0, num_lanes + 1)
    ax.set_ylim(0, 3200)
    ax.invert_yaxis()

    for lane in range(1, num_lanes + 1):
        ax.plot([lane - 0.25, lane + 0.25], [50, 50], color="gray", lw=4)

    for lane, length in enumerate(lengths, 1):
        ax.plot([lane - 0.2, lane + 0.2], [length, length], color="cyan", lw=3, alpha=0.9)

    ax.set_title("Simulated Agarose Gel Electrophoresis", fontsize=14, color="black")
    ax.set_xlabel("Sample Lanes", fontsize=12)
    ax.set_ylabel("Fragment Size (Base Pairs / Distance)", fontsize=12)
    ax.set_xticks(range(1, num_lanes + 1))
    
    plt.tight_layout()
    plt.savefig("gel_electrophoresis_simulation.png", dpi=300)
    print("\n[+] Plot saved as 'gel_electrophoresis_simulation.png'")
    plt.show()


def detect_dna_repeats(dna_seq: str, min_repeat_len=6, max_repeat_len=10, min_unit_len=3, max_unit_len=5):
    """
    Detects repetitions of total length between 6 and 10 bases,
    where the repeated unit (motif) length is between 3 and 5 bases.
    """
    print("\n" + "=" * 60)
    print("EXERCISE 2: DNA REPEAT DETECTION")
    print("=" * 60)

    found_repeats = []

    for unit_len in range(min_unit_len, max_unit_len + 1):
        pattern = f"(([A-Z]{{{unit_len}}})\\2+)"
        
        for match in re.finditer(pattern, dna_seq):
            full_repeat = match.group(1)
            repeat_unit = match.group(2)
            total_len = len(full_repeat)

            if min_repeat_len <= total_len <= max_repeat_len:
                start_pos = match.start()
                end_pos = match.end()
                num_copies = total_len // unit_len
                found_repeats.append({
                    "unit": repeat_unit,
                    "copies": num_copies,
                    "full_repeat": full_repeat,
                    "total_len": total_len,
                    "start": start_pos,
                    "end": end_pos
                })

    if found_repeats:
        print(f"{'Motif':<8} | {'Repeat Seq':<12} | {'Copies':<8} | {'Total Len':<10} | {'Position (Start-End)'}")
        print("-" * 60)
        for r in found_repeats:
            print(f"{r['unit']:<8} | {r['full_repeat']:<12} | {r['copies']:<8} | {r['total_len']:<10} | {r['start']}-{r['end']}")
    else:
        print("No repeats satisfying the exact criteria (total length 6-10 bp, unit length 3-5 bp) were found.")

if __name__ == "__main__":
    simulate_gel_electrophoresis(NCBI_DNA_SEQUENCE, num_samples=10)
    
    detect_dna_repeats(NCBI_DNA_SEQUENCE)
    