# Standard RNA Genetic Code Table
CODON_TABLE = {
    'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP',
    'UGU': 'C', 'UGC': 'C', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L',
    'CUA': 'L', 'CUG': 'L', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AUU': 'I',
    'AUC': 'I', 'AUA': 'I', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T',
    'ACG': 'T', 'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V',
    'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'GCU': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E',
    'GAG': 'E', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

raw_input = input("Enter DNA/RNA sequence: ")
clean_input = ''.join([c.upper() for c in raw_input if c.isalpha()])

rna_sequence = clean_input.replace('T', 'U')

start_index = rna_sequence.find('AUG')

if start_index == -1:
    print("\nError: No START codon (ATG / AUG) was found in the sequence!")
else:
    coding_region = rna_sequence[start_index:]
    protein = []
    stop_found = False
    used_nucleotides_count = 0
    
    for i in range(0, len(coding_region) - 2, 3):
        codon = coding_region[i:i+3]
        amino_acid = CODON_TABLE.get(codon, '?')
        
        if amino_acid == 'STOP':
            stop_found = True
            used_nucleotides_count = i + 3
            break
            
        protein.append(amino_acid)
        used_nucleotides_count = i + 3

    clean_cds = coding_region[:used_nucleotides_count]
 
    print(f"Cleaned sequence:     {clean_input}")
    print(f"Extracted CDS region: {clean_cds}")
    print(f"Amino acid sequence:  {''.join(protein)}")
    
    if stop_found:
        print("Status: STOP codon successfully identified!")
    else:
        print("Note: No STOP codon was found. Translation reached the end of the sequence.")