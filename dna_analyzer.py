import os
import matplotlib.pyplot as plt
from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO
import streamlit as st

def calculate_gc_content(dna_sequence):
    """Calculate the GC content of a DNA sequence"""
    gc_count = dna_sequence.count('G') + dna_sequence.count('C')
    return round((gc_count / len(dna_sequence)) * 100, 2)

def count_nucleotide_frequencies(dna_sequence):
    """Count the occurences of each nucleotide in the sequence"""
    return Counter(dna_sequence)

def transcribe_dna(dna_sequence):
    """Transcribe DNA to mRNA by replacing Tymine (T) with Uracil (U)"""
    return dna_sequence.replace('T', 'U')

def plot_nucleotide_frequencies(frequencies):
    """Plot pie chart of nucleotide composition"""
    plt.pie(frequencies.values(), labels=frequencies.keys(), autopct='%1.1f%%', startangle=140)
    plt.title("Nucleotide Composition")
    plt.show()

def read_fasta(file_path):
    """Read a DNA sequence from a FASTA file"""
    with open(file_path, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            return str(record.seq)
        
def find_orfs(dna_sequence):
    """Identify Open Reading Frames (ORFs) in the given DNA sequence"""
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []

    for i in range(len(dna_sequence) - 2):
        codon = dna_sequence[i:i+3]
        if codon == start_codon:
            for j in range(i, len(dna_sequence) - 2, 3):
                stop_codon = dna_sequence[j:j+3]
                if stop_codon in stop_codons:
                    orfs.append(dna_sequence[i:j+3])
                    break
    return orfs

def reverse_complement(dna_sequence):
    """Generate the reverse complement of a DNA sequence"""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complement[base] for base in reversed(dna_sequence))

def translate_dna(dna_sequence):
    """Translate DNA sequence into an amino acid sequence."""
    return str(Seq(dna_sequence).translate())


def main():
    choice = input("Do you want to enter a sequence manually (M) or use a FASTA file (F)? ").upper()

    if choice == "F":
        file_path = input("enter the path to the FASTA file: ").strip()
        if os.path.exists(file_path):
            dna_sequence = read_fasta(file_path)
            print(f"Loaded sequence: {dna_sequence}")
        else:
            print("Error: File not found.")
            return
    elif choice == "M":
        dna_sequence = input("Enter a DNA sequence: ").upper()
    else:
        print("Invalid option. Exiting.")
        return


    print(f"\nGC Content: {calculate_gc_content(dna_sequence)}%")

    frequencies = count_nucleotide_frequencies(dna_sequence)
    print("\nNucleotide Frequencies: ", frequencies)

    print("\nTranscribed mRNA:", transcribe_dna(dna_sequence))


    orfs = find_orfs(dna_sequence)
    print("Open Reading Frames (ORFs): ")
    for idx, orf in enumerate(orfs, 1):
        print(f"ORF {idx}: {orf}")

    print("Reverse Complement: ", reverse_complement(dna_sequence))

    print("Amino Acid Translation: ", translate_dna(dna_sequence))

    plot_nucleotide_frequencies(frequencies)



if __name__ == "__main__":
    main()