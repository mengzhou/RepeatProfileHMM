# RepeatProfileHMM

A profile-HMM package designed for repeat identification. Please see
INSTALL for installation instructions.

# Usage
## `construct`
A program to construct profile-HMM based on the input multiple sequence
alignment result.

### Input
The input must be in [ClustalW format][1]. Both `clustalw2` and `muscle`
can generate this format of multiple alignment.

### Parameters
1. `-l`: set lambda value
By deault, the MAP algorithm is used to determine which columns in the multiple
alignment should be included in the model. This parameter is the penalizing
factor for the MAP algorithm. The value of this parameter indicates the prior
probability in logarithm for the case that all columns of the multiple alignment
should be included for the model construction. Therefore, in most cases, this
parameter should be set to a non-positive value, and the bigger this value is,
the larger the model length will be. 

2. `-u`: turn on heuristic mode
If this option is turned on, the algorithm for choosing columns in the multiple
alignment will be switched to picking the ones with more than 50% non-gap.

3. `-m`: use guide sequence for model construction
If this option is turned on, instead of using MAP or heuristic algorithms, this
program will look for a sequence named as "consensus" in the multiple alignment,
and use it for determining match or indel states for all other sequences.

4. `-n`: do not train truncation probabilities
If the input multiple alignment does not include sequences that are known to
have 5'-end sequence, you can turn on this option to keep the truncation
probabilities of the profile-HMM unchanged.

5. `-v: show verbose information

6. `-d`: show debug information
The input multiple alignment will be visualized together with columns chosen by
both MAP and heuristic algorithms labeled specifically.

### Output
The output is a text file which includes: model length, transition probabilities
and emission probabilities. All probabilities are in log scale.

## `scan`
A program to scan the genome using the model file generated by `construct`.

### Input
The output of `construct`, an profile-HMM file.

### Parameters
1. `-c`: reference genome
Specify the path to a directory of FASTA files of the reference genome, or a
signle FASTA file with multiple sequences.

2. `-o`: output file name
The name of output file. If not set, result will be printed to the screen.

3. `-n`: do not use log-odds ratio as score
If this option is turned on, the scoring scheme during posterior decoding will
be switched to using emission probabilities only from the match states.

4. `-z`: threshold for z-score filtering for repeat identification.
After posterior decoding, all consecutive match states will be combined as
potential repeats, and the score of each individual region is its corresponding
posterior likelihood. All scores will be z-transformed and regions with z-scores
greater than the threshold will be reported as repeats. If this value is set
to 0, the filtering functionality will be disabled.

5. `-b`: use state bits instead of CIGAR string in the output
By default, a CIGAR string will be generated for each identified repeat based
on the posterior decoding result, to represent the indel states of the repeat.
If this option is turned on, a string composed of 0 or 1 will be reported, with
1 indicating a present match state and 0 indicating an absent one. The length of
this string is equal to the length of profile-HMM. This state bits string is
useful for determining which match states are present in the identified repeat.

6. `-p`: the number of processes
This option is only available when openmp is installed. Setting a number greater
than 1 can parallelize the computation.

7. `-d`: print debug information
The debug information includes the detailed posterior decoding results, e.g.
state type and index for each position.

### Output
A BED file including identified repeats, and the corresponding z-score.

## `fisher`
A program to calculate the Fisher score vector given a profile-HMM.

### Input
The output of `construct`, and profile-HMM file.

### Parameters
1. `-c`: input sequences
Specify the path to a FASTA files of multiple sequences for the computation.

2. `-o`:
Specifty the path to an output file. If this option is not set, the output will
be printed to the screen.

### Output
For each input sequence, a Fisher score vector will be reported. The dimensions
are separated by tab, and the values are in log scale. The total number of
dimensions will be 4 times the model length, including all emission parameters.

In addition, the expected state count values are also computed and appended to
the Fisher score vector. This includes (model length + 1) numbers, including
all match states and the background state, I_0. These values are not in log
scale.

# Reference
[1]: http://web.mit.edu/meme_v4.11.4/share/doc/clustalw-format.html
