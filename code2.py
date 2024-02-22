from Bio import SeqIO
from Bio.Seq import Seq
import random

# Function to introduce random mutation for one nucleotide
def introduce_mutation(seq):
    length = len(seq)
    mutation_index = random.randint(0, length - 1)
    nucleotides = ['A', 'C', 'G', 'T']
    mutation_nucleotide = random.choice([n for n in nucleotides if n != seq[mutation_index]])
    mutated_seq = seq[:mutation_index] + mutation_nucleotide + seq[mutation_index + 1:]
    return mutation_index, seq[mutation_index], mutation_nucleotide, mutated_seq

# Read the FASTA file
file_path = "extracted_subsequences.fasta"
sequences = list(SeqIO.parse(file_path, "fasta"))

# Iterate over each sequence
for seq_record in sequences:
    # Introduce mutation
    mutation_index, original_nucleotide, mutation_nucleotide, mutated_seq = introduce_mutation(str(seq_record.seq))
    seq_record.seq = Seq(mutated_seq)
    
    # Calculate mutation percentage
    seq_length = len(seq_record.seq)
    mutation_percentage = round(1 / seq_length * 100, 2)
    
    # Print mutation details
    print(f"Mutation introduced in sequence {seq_record.id} at position {mutation_index + 1}: {original_nucleotide} -> {mutation_nucleotide}")
    print(f"Mutation percentage for sequence {seq_record.id}: {mutation_percentage}%")

# Write the mutated sequences to a new FASTA file
output_file = "mutated_sequences.fasta"
with open(output_file, "w") as out_handle:
    SeqIO.write(sequences, out_handle, "fasta")

print(f"Mutated sequences written to {output_file}")

# write back the mutated sequences into original file
with open("influenza.fna.1", "a") as out_handle:
    SeqIO.write(sequences, out_handle, "fasta")
print("Mutated sequences written to original multifasta")
