import matplotlib.pyplot as plt
from collections import Counter
from Bio.Seq import Seq

def calculate_gc_content(dna_sequence):
    """Calculate the GC content of a DNA sequence"""
    gc_count = dna_sequence.count('G') + dna_sequence.count('C')
    return round((gc_count / len(dna_sequence)) * 100, 2)

def main():
    dna_sequence = input("Enter a DNA sequence: ").upper()

    print(f"GC Content: {calculate_gc_content(dna_sequence)}%")


if __name__ == "__main__":
    main()