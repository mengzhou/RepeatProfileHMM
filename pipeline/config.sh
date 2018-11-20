#!/bin/bash
# Set number of threads for parallelization
NTHREAD=16
# Set path to dependencies here
# 1. path to Tandem Repeats Finder executable binary
TRF=/home/rcf-40/mengzhou/panfs/repeats/cyclic/sim/trf407b.linux64
# 2. path to MUSCLE executable binary
MUSCLE=/home/rcf-40/mengzhou/panfs/tools/muscle3.8.31_i86linux64
# 3. path to HMMER *directory* of binaries. E.g.: /home/user/HMMER/binaries
HMMER=/home/rcf-40/mengzhou/panfs/tools/hmmer-3.1b1/binaries
# 4. path to bedtools *directory* of binaries
BEDTOOLS=/home/rcf-40/mengzhou/bin/bedtools/bin

# 5. path to RepeatProfileHMM *root* directory
RPHMM=/home/rcf-40/mengzhou/panfs/repeats/RepeatProfileHMM

# set path to reference genome, must be a single FASTA file
#REF=/home/rcf-40/mengzhou/panfs/repeats/mm10/RepeatMasker/L1Base/mm10.fa
REF=/home/rcf-40/mengzhou/panfs/genome/mm10/chr19.fa

# validation of paths. no need to change things below
REV_COMP=$RPHMM/utils/revcomp_fa.py
GET_MEDIAN=$RPHMM/utils/get_median.awk
NUMBERING=$RPHMM/utils/numbering.awk
for i in $TRF $MUSCLE $HMMER $NHMMER $BEDTOOLS $RPHMM $REF $REV_COMP
do
  if [ ! -e $i ]
  then
    echo "Path $i does not exist! Please check config.sh"
    exit
  fi
done
