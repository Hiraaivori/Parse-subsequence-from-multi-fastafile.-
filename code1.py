# Read the file and fetch sequences

from Bio import SeqIO
file = open("influenza.fna")
scaffolds = list(SeqIO.parse(file, 'fasta'))
file.close()
print(f"Found {len(scaffolds)} scaffold sequences")

# extract one-by-one scaffold sequences from scaffold list
print("\nThe first scaffold\n")
first_scaffold = scaffolds[0] # as Python counts from zero index
print(first_scaffold.id)
print(repr(first_scaffold.seq)) # only print the sequence without header
print(len(first_scaffold))

print("\nThe second scaffold\n")
second_scaffold = scaffolds[1]
print(second_scaffold.id)
print(repr(second_scaffold.seq))
print(len(second_scaffold))

print("\nThe Last record\n")
last_scaffold = scaffolds[-1]
print(last_scaffold.id)
print(repr(last_scaffold.seq))
print(len(last_scaffold))

# first_scaffold, second_scaffold and last_scaffold are the extracted subsequences
# write the subsequences into a new FASTA (.fasta) file

new_fasta = "extracted_subsequences.fasta"
with open(new_fasta, 'w') as new_file:
    SeqIO.write([first_scaffold, second_scaffold, last_scaffold], new_file, 'fasta')
print("Extracted subsequences written to new FASTA file:", new_fasta)

# Remove the extracted subsequences from the original multifasta file
scaffolds_to_remove = [first_scaffold, second_scaffold, last_scaffold]
remaining_scaffolds = [scaffold for scaffold in scaffolds if scaffold.id not in [s.id for s in scaffolds_to_remove]] # this command comapres the SeqIO objects'
                                                                                                                     # IDs with extracted SeqIO objct ID (s.id) 
                                                                                                                     # scaffold.id -> orginal file scaffolds
                                                                                                                     # s.id -> extracted sequence id
# Write the remaining sequences back to the original file
file_path = "influenza.fna" # at first only file=open(filename) was used, but python expects 'str' type in file_path object 
with open(file_path, 'w') as original_file:
    SeqIO.write(remaining_scaffolds, original_file, 'fasta')
print("Extracted sequences removed from the original multifasta")

# Reverse complement of each sequences
first_reverse = first_scaffold.reverse_complement()
second_reverse = second_scaffold.reverse_complement()
last_reverse = last_scaffold.reverse_complement()
print(first_reverse, second_reverse, last_reverse)
