# Parse-subsequence-from-multi-fastafile.-
Extraction of Subsequence from multi-fasta file
Extract a subsequence from a given primary multi-fasta file (genome fasta) and a bed file. The extracted subsequence must be removed from the primary sequence (multifasta). You can create your own BED file. It should contain coordinates for at least 3 scaffolds.

Reverse complement the extracted sequences.

Introduce random 1 nucleotide change in the above extracted sequence(s); this becomes the processed subsequence(s).

Insert the processed subsequence(s) back to the primary sequence(multifasta).
