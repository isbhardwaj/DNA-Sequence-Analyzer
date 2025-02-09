import matplotlib.pyplot as plt
from collections import Counter
from Bio.Seq import Seq

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

def main():
    dna_sequence = input("Enter a DNA sequence: ").upper()

    print(f"\nGC Content: {calculate_gc_content(dna_sequence)}%")

    frequencies = count_nucleotide_frequencies(dna_sequence)
    print("\nNucleotide Frequencies: ", frequencies)

    print("\nTranscribed mRNA:", transcribe_dna(dna_sequence))

    plot_nucleotide_frequencies(frequencies)



if __name__ == "__main__":
    main()