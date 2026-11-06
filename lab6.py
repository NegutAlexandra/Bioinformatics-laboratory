import random

# Complementul bazelor azotate pentru Inverted Repeats
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def get_reverse_complement(seq: str) -> str:
    """Returnează inversul complementat al unei secvențe ADN."""
    return "".join(COMPLEMENT[base] for base in reversed(seq))

def generate_random_dna(length: int) -> str:
    """Generează o secvență ADN aleatorie."""
    return "".join(random.choice(['A', 'C', 'G', 'T']) for _ in range(length))
\
def create_artificial_dna_with_transposons():
    base_dna = list(generate_random_dna(300))
    
    ir_motif_1 = "ACGTAC"
    ir_motif_2 = "GTTACGA"
    ir_motif_3 = "TCGATC"
    
    dr_1 = "ATGCA"  
    t1_content = dr_1 + ir_motif_1 + "GATCGATCGATC" + get_reverse_complement(ir_motif_1) + dr_1
    
    dr_2 = "GCTAC"
    t2_content = dr_2 + ir_motif_2 + "AAAAATTTTTCCCCC" + get_reverse_complement(ir_motif_2) + dr_2

    dna_seq = "".join(base_dna[:40]) + t1_content + "".join(base_dna[40:150]) + t2_content + "".join(base_dna[150:])
    
    dr_3 = "CGTA"
    t3_content = dr_3 + ir_motif_3 + "CCGGTTGGCC" + get_reverse_complement(ir_motif_3) + dr_3
    
    insert_pos = 80
    dna_seq = dna_seq[:insert_pos] + t3_content + dna_seq[insert_pos:]
    
    return dna_seq

def detect_transposable_elements(dna_seq: str, min_ir_len=5, max_ir_len=8, min_dr_len=3, max_dr_len=6):
    """
    Detectează transposonii căutând structura: Direct Repeat (DR) + Inverted Repeat (IR)
    """
    detected_elements = []
    seq_len = len(dna_seq)
    
    for i in range(seq_len - 15):
        for ir_len in range(min_ir_len, max_ir_len + 1):
            left_ir = dna_seq[i:i + ir_len]
            expected_right_ir = get_reverse_complement(left_ir)
            
            for j in range(i + ir_len + 10, min(i + ir_len + 100, seq_len - ir_len)):
                right_ir = dna_seq[j:j + ir_len]
                
                if right_ir == expected_right_ir:
                    for dr_len in range(min_dr_len, max_dr_len + 1):
                        if i >= dr_len and (j + ir_len + dr_len) <= seq_len:
                            left_dr = dna_seq[i - dr_len : i]
                            right_dr = dna_seq[j + ir_len : j + ir_len + dr_len]
                            
                            if left_dr == right_dr:
                                start_pos = i - dr_len
                                end_pos = j + ir_len + dr_len
                                
                                detected_elements.append({
                                    "start": start_pos,
                                    "end": end_pos,
                                    "length": end_pos - start_pos,
                                    "dr": left_dr,
                                    "ir_left": left_ir,
                                    "ir_right": right_ir
                                })
                                break
                                
    unique_elements = []
    for elem in detected_elements:
        if not any(abs(elem['start'] - u['start']) < 3 and abs(elem['end'] - u['end']) < 3 for u in unique_elements):
            unique_elements.append(elem)
            
    return unique_elements

if __name__ == "__main__":
    dna_sequence = create_artificial_dna_with_transposons()
    
    print("=" * 70)
    print(f"LUNGIME SECVENȚĂ GENERATĂ: {len(dna_sequence)} bp")
    print("=" * 70)
    print(f"Secvență: {dna_sequence}\n")
    
    elements = detect_transposable_elements(dna_sequence)
    
    print("=" * 70)
    print(f"ELEMENTE TRANSPOZABILE DETECTATE ({len(elements)} găsite):")
    print("=" * 70)
    print(f"{'No.':<4} | {'Start':<6} | {'End':<6} | {'Lungime':<8} | {'Direct Repeat':<14} | {'Inverted Repeat'}")
    print("-" * 70)
    
    for idx, elem in enumerate(elements, 1):
        print(f"{idx:<4} | {elem['start']:<6} | {elem['end']:<6} | {elem['length']:<8} | {elem['dr']:<14} | {elem['ir_left']}...{elem['ir_right']}")
    print("=" * 70)