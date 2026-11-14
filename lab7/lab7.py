import urllib.request
import re
import matplotlib.pyplot as plt


RESTRICTION_ENZYMES = {
    "EcoRI": {
        "site": "GAATTC",
        "cut_offset": 1,
        "description": "5'-G^AATTC-3'"
    },
    "BamHI": {
        "site": "GGATCC",
        "cut_offset": 1,
        "description": "5'-G^GATCC-3'"
    },
    "HindIII": {
        "site": "AAGCTT",
        "cut_offset": 1,
        "description": "5'-A^AGCTT-3'"
    },
    "TaqI": {
        "site": "TCGA",
        "cut_offset": 1,
        "description": "5'-T^CGA-3'"
    },
    "HaeIII": {
        "site": "GGCC",
        "cut_offset": 2,
        "description": "5'-GG^CC-3'"
    }
}

FALLBACK_BACTERIAL_DNA = (
    "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCT"
    "GAATTCACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCAACTCAG"
    "GACTTGTTCTTACCTTTCTTTTCCAATGTTACTTGGTTCCATGCTATACATGTCTCTGGGACCAATGGTACTAAGAGG"
    "TTTGATAACCCTGTCCTACCATTTAATGATGGTGTTTATTTTGCTTCCACTGAGAAGTCTAACATAATAAGAGGCTGG"
    "GGATCCGGTACTACTTTAGATTCGAAGACCCAGTCCCTACTTATTGTTAATAACGCTACTAATGTTGTTATTAAAGTC"
    "TGTGAATTTCAATTTTGTAATGATCCATTTTTGGGTGTTTATTACCACAAAAACAACAAAAGTTGGATGGAAAGTGAG"
    "TTCAGAGTTTATTCTAGTGCGAATAATTGCACTTTTGAATATGTCTCTCAGCCTTTTCTTATGGACCTTGAAGGAAAA"
    "AAGCTTATTTCAAAAATCTTAGGGAATTTGTGTTTAAGAATATTGATGGTTATTTCAAAATATATTCTAAGCACACC"
    "CCTATTAATCTTGTTCGCGATCTCCCTCAGGGTTTTCGGCTTTAGAACCATTGGTAGATTTGCCAATAGGTATTAAC"
    "ATCACTAGGTTTCAAACTCTACTTGCCCTACATAGAAGTTATTTGACTCCTGGTGATTCTTCTTCAGGTTGGACAGCT"
    "GGCCGCAGCTTATTATGTGGGTTATCTTCAACCTAGGACTTTTCTATTAAAATATAATGAAAATGGAACCATTACA"
    "GATGCTGTAGACTGTGCACTTGACCCTCTCTCAGAAACAAAGTGTACGTTGAAATCCTTCACTGTAGAAAAAGGAATC"
    "TATCAAACTTCTAACTTTAGAGTCCAACCAACAGAATCTATTGTTAGATTTCCTAATATTACAAACTTGTGCCCTTTT"
    "TCGATTTTTAACGCCACCAGATTTGCATCTGTTTATGCTTGGAACAGGAAGAGAATCAGCAACTGTGTTGCTGAC"
    "TATTCTGTCCTATATAATTCCACATCATTTTCCACTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTC"
    "GAATTCGGATCCAAGCTTTCGAGGCCATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATC"
) * 3

def get_bacterial_dna(accession="NC_000908.2", length=4000) -> str:
    """
    Descarcă un fragment dintr-un genom bacterian de pe NCBI (Mycoplasma genitalium).
    """
    print(f"[+] Încercăm descărcarea genomului bacterian de pe NCBI ({accession})...")
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    
    try:
        req = urllib.request.urlopen(url, timeout=8)
        lines = req.read().decode('utf-8').splitlines()
        seq = "".join(line.strip().upper() for line in lines if not line.startswith('>'))
        print(f"[+] Secvență descărcată cu succes! Extragem primele {length} bp pentru analiză.")
        return seq[:length]
    except Exception as e:
        print(f"[!] Nu s-a putut accesa NCBI ({e}). Folosim secvența bacteriană locală ({len(FALLBACK_BACTERIAL_DNA)} bp).")
        return FALLBACK_BACTERIAL_DNA

def digest_dna(dna_seq: str, enzymes: dict):
    """
    Găsește toate pozițiile de tăiere (cleavage sites) și calculează lungimea fragmentelor rezultate.
    """
    seq_len = len(dna_seq)
    results = {}
    all_cuts = set()

    print("\n" + "=" * 75)
    print(f"REZULTATELE DIGESTIEI ENZIMATICE (Lungime ADN: {seq_len} bp)")
    print("=" * 75)

    for name, info in enzymes.items():
        site = info["site"]
        offset = info["cut_offset"]
        
        cut_positions = []
        for match in re.finditer(f"(?={site})", dna_seq):
            cut_pos = match.start() + offset
            cut_positions.append(cut_pos)
            all_cuts.add(cut_pos)
            
        boundaries = [0] + sorted(cut_positions) + [seq_len]
        fragments = [boundaries[i+1] - boundaries[i] for i in range(len(boundaries) - 1)]
        
        results[name] = {
            "num_cleavages": len(cut_positions),
            "cut_positions": cut_positions,
            "fragments": sorted(fragments, reverse=True)
        }
        
        print(f"\nEnzima: {name} ({info['description']}) | Secvență recunoscută: 5'-{site}-3'")
        print(f"  -> Număr de tăieturi (cleavages): {len(cut_positions)}")
        print(f"  -> Pozițiile tăieturilor (bp):    {cut_positions if cut_positions else 'Nicio tăietură'}")
        print(f"  -> Lungimea fragmentelor (bp):    {results[name]['fragments']}")

    all_boundaries = [0] + sorted(list(all_cuts)) + [seq_len]
    combined_fragments = sorted([all_boundaries[i+1] - all_boundaries[i] for i in range(len(all_boundaries) - 1)], reverse=True)
    results["All 5 Combined"] = {
        "num_cleavages": len(all_cuts),
        "cut_positions": sorted(list(all_cuts)),
        "fragments": combined_fragments
    }
    
    print("\n" + "-" * 75)
    print(f"DIGESTIE COMBINATĂ (Toate cele 5 enzime simultan):")
    print(f"  -> Număr total de tăieturi:       {len(all_cuts)}")
    print(f"  -> Lungimea fragmentelor (bp):    {combined_fragments}")
    print("=" * 75)
    
    return results, seq_len

def simulate_electrophoresis_gel(digest_results: dict, total_dna_len: int):
    """
    Simulează migrarea fragmentelor ADN pe un gel de agaroză.
    Fragmentele mici migrează mai rapid (spre partea de jos a gelului),
    iar fragmentele mari rămân aproape de godeuri (partea de sus).
    """
    lanes = ["DNA Ladder"] + list(digest_results.keys())
    num_lanes = len(lanes)
    
    ladder_bands = [4000, 3000, 2000, 1500, 1000, 750, 500, 250, 100]

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_facecolor("#121212") 
    ax.set_yscale("log")
    ax.set_ylim(50, max(total_dna_len + 1000, 5000))
    ax.invert_yaxis()

    ax.set_xlim(0, num_lanes + 1)
    ax.set_xticks(range(1, num_lanes + 1))
    ax.set_xticklabels(lanes, rotation=20, ha="right", fontsize=10, fontweight="bold")
    ax.set_ylabel("Fragment Size (bp - Log Scale)", fontsize=11, fontweight="bold")
    ax.set_title("Agarose Gel Electrophoresis Simulation - Restriction Digest", fontsize=13, fontweight="bold", pad=15)

    top_y = ax.get_ylim()[0]
    for lane_idx in range(1, num_lanes + 1):
        ax.plot([lane_idx - 0.3, lane_idx + 0.3], [top_y * 1.05, top_y * 1.05], color="#555555", lw=5)

    for band in ladder_bands:
        ax.plot([1 - 0.25, 1 + 0.25], [band, band], color="#ffaa00", lw=2.5, alpha=0.85)
        ax.text(0.65, band, f"{band}", color="#ffaa00", fontsize=8, va="center", ha="right")

    for idx, (enzyme_name, data) in enumerate(digest_results.items(), start=2):
        fragments = data["fragments"]
        for frag in fragments:
            if frag >= 50:
                intensity = min(1.0, 0.5 + (frag / total_dna_len) * 0.5)
                ax.plot([idx - 0.25, idx + 0.25], [frag, frag], color="#00ffff", lw=2.8, alpha=intensity)

    ax.grid(axis="y", linestyle="--", alpha=0.2, color="white")
    
    plt.tight_layout()
    output_img = "restriction_digest_gel.png"
    plt.savefig(output_img, dpi=300)
    print(f"\n[+] Simularea gelului de electroforeză a fost salvată ca '{output_img}'.")
    plt.show()

if __name__ == "__main__":
    bacterial_dna = get_bacterial_dna(accession="NC_000908.2", length=4000)
    
    digest_data, dna_length = digest_dna(bacterial_dna, RESTRICTION_ENZYMES)
    
    simulate_electrophoresis_gel(digest_data, dna_length)