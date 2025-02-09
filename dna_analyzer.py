import matplotlib.pyplot as plt
from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO
import os

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
    """Plot a bar graph of nucleotide frequencies"""
    plt.bar(frequencies.keys(), frequencies.values(), color=['blue', 'red', 'green', 'yellow'])
    plt.xlabel('Nucleotides')
    plt.ylabel('Frequency')
    plt.title('Nucleotide Frequency in DNA Sequence')
    plt.show()

def read_fasta(file_path):
    """Read a DNA sequence from a FASTA file"""
    with open(file_path, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            return str(record.seq)

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

    plot_nucleotide_frequencies(frequencies)



if __name__ == "__main__":
    main()