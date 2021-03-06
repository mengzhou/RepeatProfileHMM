# Pipelines of genome-wide monomer identification
The scripts used for identifying monomers in the genome are included in this
directory.

# Dependencies
Here is a list of required tools for this pipeline:

* [TandemRepeatsFinder](https://tandem.bu.edu/trf/trf.html)
* [MUSCLE](https://www.drive5.com/muscle/)
* [HMMER](http://hmmer.org/)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)

# Scripts in this directory
## config.sh
The file for setting up paths of dependencies. Please change according to
your environment before running other pipeline scripts.

## model_construct.sh
The script for constructing profile-HMM based on a set of L1Md promoters.

## iteration_rphmm.sh
The script for training model parameters using sequences identified as
monomers.

# Instructions
1. Make sure all the above dependencies have been installed.

2. Change the variables in `config.sh` accordingly to your environment.
Note that some variables correspond to *executable programs themselves*,
but others correspond to *directories* which includes multiple programs.
Also you need to properly set the path in both `model_construct.sh` and
`iteration_rphmm.sh`.

4. If you used `nhmmer` for the initial detection of monomers, you need
to reformat the output of `nhmmer` to BED format. For example, below is
the command which was used for the monomer detection in mm10 genome:
```
nhmmer --dna --dfamtblout <nhmmer_out> --max --cut_ga <monomer.hmm> mm10.fa
```
Here the output format is set to the Dfam tabulated format. To convert this
format to BED, you can use the following `awk` command:
```
awk 'BEGIN{OFS="\t"}$1!~/#/{if($9=="+"){start=$10-1; end=$11;}
  else{start=$11-1; end=$10} if($4 >= 60)print $1, start, end, $3, end-start,
  $9, $4, $5}' <nhmmer_out> | sort -k1,1 -k2,2n -k3,3n > <nhmmer_bed>
```
Replace `<nhmmer_out>` and `<nhmmer_bed>` to corresponding file names.
Note here we applied a filter of minimum bit socre of 60, to remove potential
noise.

3. Run `model_construct.sh` to get an initial profile-HMM. It takes 3
parameters:
  * The result of `nhmmer` scan, which includes individual monomers.
  * The length of extension for each individual promoter. Recommended
    value is 500.
  * The file which includes chromosome sizes. It can be obtained via
    `fetchChromSizes`. You can find this program in the UCSC kent tools:
    http://hgdownload.soe.ucsc.edu/admin/exe/

This script will estimate the optimal model length using information
obtained from the detected monomers. The output files include a profile-HMM
file with `.params` extension, a `consensus.fa` file with the consensus
sequence of the model, and a FASTA file with the sequences of individual
extended promoters, named as `[monomer_name].ext[extension_length].fa`.

4. Run `iteration_rphmm.sh` to train model parameters through an iterative
process. This scripts takes 3 parameters:
  * The `consensus.fa` generated by `model_construct.sh`.
  * The `nhmmer` scan output which is also the input of `model_construct.sh`.
  * The FASTA file of extended promoter sequences generated before.

This script will train the emission probabilities of the initial profile-HMM.
It will generate a number of intermediate files. When the training process
terminates, a finalized profile-HMM file with name
`[monomer_name].ext[extension_length].final.params` will be generated. This
file includes parameters of the finalized model.
