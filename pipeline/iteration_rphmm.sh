#!/bin/bash
# TO-DO: deambiguity_GT.bed needs to be changed to something more general
NHMMER=~mengzhou/panfs/repeats/mm10/RepeatMasker/L1Base/scan/model_construction/hmmer/nhmmer
MUSCLE=~mengzhou/panfs/tools/muscle3.8.31_i86linux64
ESL=~mengzhou/panfs/tools/hmmer-3.1b1/binaries/esl-reformat
BUILD=~mengzhou/panfs/tools/hmmer-3.1b1/binaries/hmmbuild
RPHMM=/home/rcf-40/mengzhou/panfs/repeats/RepeatProfileHMM/bin
NUMBER=/home/rcf-40/mengzhou/panfs/repeats/RepeatProfileHMM/pipeline/get_number.sh
NTHREAD=16

#GENOME=/home/rcf-40/mengzhou/panfs/repeats/ms1509/ms1509.fa
GENOME=/home/rcf-40/mengzhou/panfs/repeats/mm10/RepeatMasker/L1Base/mm10.fa
RAND=500
MAX_ITR=9
JACCARD=0

function float_lt () {
  echo $1"<"$2 | bc -l
}
function prep_seq() {
  # $1: consensus.fa
  # $2: nhmmer_scan.bed of previous iteration
  # $3: name of the output, usually the name of current iteration
  # $4: number of sequences to sample
  # $5: length of the consensus
  echo "[EXTRACTING SEQ]"
  REG=${3}.bed
  #echo $1, $2, $REG, $4, $5
  TEMP=$(mktemp tmp.$$.XXXXXXXX)
  # the score is length normalized -log(E-value)
  cat $2 | awk -v LEN=$5 -v EXT=15 'BEGIN{OFS="\t"}{\
      score=-log($8)/$5; if($5<=LEN+EXT && $5>=LEN-EXT)\
      {print $1, $2, $3, $4"|"$1":"$2"-"$3"|"$6, score, $6}}' | sort -k5,5gr > $TEMP
  TOTAL=$(cat $TEMP | wc -l)
  echo $TOTAL
  # keep the first 1/3 scoring hits for random sampling
  cat $TEMP | sort -k5,5gr | awk -v t=$TOTAL 'BEGIN{OFS="\t"; srand()}NR<=t/1.5\
      {print rand(), $0}' | sort -k1,1g | cut -f 2- | awk -v s=$4 'NR<=s' | \
    sort -k1,1 -k2,2n -k3,3n > $REG
  rm $TEMP
  bedtools getfasta -fi $GENOME -s -name -bed $REG -fo ${3}.fa
  cat $1 >> ${3}.fa
}

function muscle() {
  echo "[ALIGNING]"
  $MUSCLE -in $1 -clw -quiet > ${1%.*}.aln
  $ESL stockholm ${1%.*}.aln | sed 's/consensus/#=GC RF/' > ${1%.*}.sto
}

function build_nhmmer() {
  $BUILD --hand --dna --plaplace --wnone --enone ${1%.*}.hmm $1
}

function build_rphmm() {
  $RPHMM/construct -n -m -o ${1%.*}.params $1
}

function scan_nhmmer() {
  # $1: profile-HMM
  # $2: reference genome FASTA
  # $3: output BED name
  # $4: (optional) filtering region for removing ambiguity in F type monomer detection
  $NHMMER -o /dev/null --cpu 16 --dna --dfamtblout ${3%.*}.out --max --cut_ga $1 $2
  awk 'BEGIN{OFS="\t"}$1!~/#/{if($9=="+"){start=$10-1; end=$11;}\
      else{start=$11-1; end=$10} len=end-start; score=-log($5)/len;\
      if(score >= .2)print $1, start, end, $3, end-start, $9, $4, $5}' ${3%.*}.out | \
    sort -k1,1 -k2,2n -k3,3n > $3
  if [ $4 ]
  then
    mv $3 ${3}.all
    intersectBed -a ${3}.all -b $4 -v > $3
  fi
}

function scan_rphmm() {
  # $1: profile-HMM
  # $2: reference genome FA
  # $3: output BED
  $RPHMM/scan -v -p $NTHREAD -c $2 $1 | \
    sort -k1,1 -k2,2n -k3,3n > ${3%.*}.out
  awk 'BEGIN{OFS="\t"}{split($1, f, "|"); split(f[2], ff, "[:-]");\
      chr=ff[1]; start=ff[2]; end=ff[3];\
      if(length(f)==2){strand=$6; print chr, start+$2, start+$3, $1, $3-$2, strand, $7}\
      else{strand=f[3];if(strand=="+"){print chr, start+$2, start+$3, $1, $3-$2, strand, $7}\
      else{print chr, end-$3, end-$2, $1, $3-$2, strand, $7}}}' ${3%.*}.out | \
    sort -k1,1 -k2,2n -k3,3n > $3
}

function number_monomers() {
  # $1: scan result .bed
  NAME=${1%%.*}
  cp ${1%.*}.out ${NAME}.final.out
  $NUMBER ${NAME}.final.out
}

if [ $# -lt 3 ]
then
  echo "Usage: $0 <CONSENSUS.FA> <HMMER_OUTPUT_BED> <SEQ_OF_ROI>"
  exit
else
  export -f scan_rphmm
  LEN=$(tail -n +2 $1 | awk '{print length($1)}')

  PREV_OUT=$2
  REF=$3
  ITR=1
  CONVERGE=$(float_lt $JACCARD 0.99)
  while [ \( $ITR -le $MAX_ITR \) -a $CONVERGE -eq 1 ]
  do
    echo "[ITERATION $ITR]"
    NAME=${PREV_OUT%%.*}.${ITR}.20bp.rand${RAND}
    CURR_OUT=${NAME}.scan.bed
    prep_seq $1 $PREV_OUT $NAME $RAND $LEN
    muscle ${NAME}.fa

    #nhmmer
    #build_nhmmer ${NAME}.sto
    #echo "[SCANNING]"
    #scan_nhmmer ${NAME}.hmm $GENOME $CURR_OUT

    #rphmm
    echo "[BUILDING]"
    build_rphmm ${NAME}.aln
    echo "[SCANNING]"
    scan_rphmm ${NAME}.params $REF $CURR_OUT

    INTER=$(intersectBed -a $PREV_OUT -b $CURR_OUT | sort -k1,1 -k2,2n -k3,3n | mergeBed\
      | awk 'BEGIN{t=0}{t+=$3-$2}END{print t}')
    UNION=$(cat $PREV_OUT $CURR_OUT | cut -f 1-3 | sort -k1,1 -k2,2n -k3,3n | \
      mergeBed | awk 'BEGIN{t=0}{t+=$3-$2}END{print t}')
    JACCARD=$(echo "$INTER / $UNION" | bc -l)
    echo "[JACCARD=$JACCARD]"
    CONVERGE=$(float_lt $JACCARD 0.95)
    let ITR=ITR+1
    PREV_OUT=$CURR_OUT
    echo ""
  done
  number_monomers $CURR_OUT
fi
