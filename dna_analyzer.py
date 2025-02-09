import os
import unittest
import streamlit as st
import matplotlib.pyplot as plt
from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO

def calculate_gc_content(dna_sequence):
    """Calculate the GC content of a DNA sequence"""
    cleaned_sequence = dna_sequence.strip().replace("\n", "").replace(" ", "")
    gc_count = cleaned_sequence.count('G') + cleaned_sequence.count('C')
    return round((gc_count / len(dna_sequence)) * 100, 2) if cleaned_sequence else 0

def count_nucleotide_frequencies(dna_sequence):
    """Count the occurences of each nucleotide in the sequence"""
    cleaned_sequence = dna_sequence.strip().replace("\n", "").replace(" ", "")
    return Counter(cleaned_sequence)

def transcribe_dna(dna_sequence):
    """Transcribe DNA to mRNA by replacing Tymine (T) with Uracil (U)"""
    cleaned_sequence = dna_sequence.strip().replace("\n", "").replace(" ", "")
    return cleaned_sequence.replace('T', 'U')

def plot_nucleotide_frequencies(frequencies):
    """Plot a bar graph of nucleotide frequencies."""
    fig, ax = plt.subplots()
    ax.bar(frequencies.keys(), frequencies.values(), color=['blue', 'red', 'green', 'yellow'])
    ax.set_xlabel('Nucleotides')
    ax.set_ylabel('Frequency')
    ax.set_title('Nucleotide Frequency in DNA Sequence')
    st.pyplot(fig)

def plot_nucleotide_pie(frequencies):
    """Plot pie chart of nucleotide composition"""
    fig, ax = plt.subplots()
    ax.pie(frequencies.values(), labels=frequencies.keys(), autopct='%1.1f%%', startangle=140)
    ax.set_title("Nucleotide Composition")
    st.pyplot(fig)

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
    cleaned_sequence = dna_sequence.strip().replace("\n", "").replace(" ", "")

    for i in range(len(cleaned_sequence) - 2):
        codon = cleaned_sequence[i:i+3]
        if codon == start_codon:
            for j in range(i, len(cleaned_sequence) - 2, 3):
                stop_codon = cleaned_sequence[j:j+3]
                if stop_codon in stop_codons:
                    orfs.append(cleaned_sequence[i:j+3])
                    break
    return orfs

def reverse_complement(dna_sequence):
    """Generate the reverse complement of a DNA sequence"""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    cleaned_sequence = dna_sequence.strip().replace("\n", "").replace(" ", "")
    return "".join(complement[base] for base in reversed(cleaned_sequence))

def translate_dna(dna_sequence):
    """Translate DNA sequence into an amino acid sequence."""
    cleaned_sequence = dna_sequence.strip().replace("\n", "").replace(" ", "")
    return str(Seq(cleaned_sequence).translate())


# class TestDataAnalyzer(unittest.TestCase):
#     def tes_gc_content(self):
#         self.assertEqual(calculate_gc_content("GCGCGC"), 100.0)
#         self.assertEqual(calculate_gc_content("ATATAT"), 0.0)

#     def test_transcription(self):
#         self.assertEqual(transcribe_dna("ATGC"), "AUGC")

#     def test_reverse_complement(self):
#         self.assertEqual(reverse_complement("ATGC"), "GCAT")

#     def test_translation(self):
#         self.assertEqual(translate_dna("ATGCGT"), "MR")


def main():
    st.title("DNA Sequence Analyzer")
    st.write("Upload a FASTA file or enter a DNA sequence below.")

    uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta"])
    dna_sequence = ""

    if uploaded_file is not None:
        dna_sequence = read_fasta(uploaded_file)
        st.write(f"Loaded sequence: {dna_sequence[:50]}")
    else:
        dna_sequence = st.text_area("Or enter a DNA sequence: ").strip().upper()

    if dna_sequence:
        st.write(f"GC Content: {calculate_gc_content(dna_sequence)}%")
        frequencies = count_nucleotide_frequencies(dna_sequence)
        st.write("Nucleotide Frequencies: ", frequencies)
        st.write("Transcribed mRNA:", transcribe_dna(dna_sequence))


        orfs = find_orfs(dna_sequence)
        st.write("Open Reading Frames (ORFs): ")
        for idx, orf in enumerate(orfs, 1):
            st.write(f"ORF {idx}: {orf}")

        dna_sequence = dna_sequence.strip().replace("\n", "").replace(" ", "")

        st.write("Reverse Complement: ", reverse_complement(dna_sequence))

        st.write("Amino Acid Translation: ", translate_dna(dna_sequence))

        plot_nucleotide_frequencies(frequencies)
        plot_nucleotide_pie(frequencies)



if __name__ == "__main__":
    main()