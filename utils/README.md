## align
This program generates global alignment of DNA sequences and report
the Hamming distance(s) between each pair of sequences.

Input mode 1: stdin. In this mode, each line represents a sequence
and the alignment will be performed for each pair.
Input mode 2: an FASTA file. In this mode, all the sequences with
corresponding names will be loaded; then each pair will be aligned.

Output: a symmetric matrix of the Hamming distances. If stdin is used,
then sequence names will be assigned as numbers.
