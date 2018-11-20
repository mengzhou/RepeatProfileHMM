# Toy test example
This directory includes a toy example for demonstrating the programs of
`RepeatProfileHMM`. Below is a list of files in this example:

*  consensus.fa : the sequences of repeats inserted into the genome.

* `chrF.fa`: a fake chromosome with 100 simulated repeat insertions.

* `true_locations.bed`: the insertion sites of simulated copies. All copies
were inserted to the positive strand. The name of each repeat indicates
the type and number of mutations simulated. `S` means substitution; `I` means
insertion; `D` means deletion; `T` means 5'-end truncation.

* `repeat_sequences.fa`: the sequences of all inserted repeats, simulated from
the consensus sequence by adding truncation, indel and substitution.

* `train.aln`: the multiple sequence alignment of the first 50 copies simulated.

* `model.params`: the model built from the above alignment using

    `../bin/construct -o model.params train.aln`

Turn on the -d flag for construct to see which columns are chosen.

* `repeat_sequences.hits`: identification results using

    `../bin/scan -c chrF.fa -o repeat_sequences.hits model.params`

Each region is found by posterior decoding and a posterior z-score
is calculated. Regions with z-score greater than 3 are labeled by "COPY".
