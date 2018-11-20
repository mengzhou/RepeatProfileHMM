#!/bin/bash
# Set the below to RepeatProfileHMM *root* directory
RPHMM=/home/rcf-40/mengzhou/panfs/repeats/RepeatProfileHMM/
source ${RPHMM}/pipeline/config.sh
RAND=500
MAX_ITR=9

function float_lt () {
  echo $1"<"$2 | bc -l
}
function prep_seq() {
  # $1: consensus.fa
  # $2: nhmmer_scan.bed of previous iteration
  # $3: name of the output, usually the name of current iteration
  # $4: number of sequences to sample
  # $5: length of the consensus
  # $6: reference genome
  echo "[EXTRACTING SEQ]"
  REG=${3}.bed
  #echo $1, $2, $REG, $4, $5
  TEMP=$(mktemp tmp.$$.XXXXXXXX)
  # the score is length normalized -log(E-value)
  cat $2 | awk -v LEN=$5 -v EXT=20 'BEGIN{OFS="\t"}{\
      score=-log($8)/$5; if($5<=LEN+EXT && $5>=LEN-EXT)\
      {print $1, $2, $3, $4"|"$1":"$2"-"$3"|"$6, score, $6}}' | sort -k5,5gr > $TEMP
  TOTAL=$(cat $TEMP | wc -l)
  # keep the first 1/3 scoring hits for random sampling
  cat $TEMP | sort -k5,5gr | awk -v t=$TOTAL 'BEGIN{OFS="\t"; srand()}NR<=t/1.5\
      {print rand(), $0}' | sort -k1,1g | cut -f 2- | awk -v s=$4 'NR<=s' | \
    sort -k1,1 -k2,2n -k3,3n > $REG
  rm $TEMP
  $BEDTOOLS/bedtools getfasta -fi $6 -s -name -bed $REG -fo ${3}.fa
  cat $1 >> ${3}.fa
}

function muscle() {
  echo "[ALIGNING]"
  $MUSCLE -in $1 -clw -quiet > ${1%.*}.aln
  $HMMER/esl-reformat stockholm ${1%.*}.aln | sed 's/consensus/#=GC RF/' > ${1%.*}.sto
}

function build_nhmmer() {
  $HMMER/hmmbuild --hand --dna --plaplace --wnone --enone ${1%.*}.hmm $1
}

function build_rphmm() {
  $RPHMM/bin/construct -n -m -o ${1%.*}.params $1
}

function scan_nhmmer() {
  # $1: profile-HMM
  # $2: reference genome FASTA
  # $3: output BED name
  # $4: (optional) filtering region for removing ambiguity in F type monomer detection
  $HMMER/nhmmer -o /dev/null --cpu $NTHREAD --dna --dfamtblout ${3%.*}.out --max --cut_ga $1 $2
  awk 'BEGIN{OFS="\t"}$1!~/#/{if($9=="+"){start=$10-1; end=$11;}\
      else{start=$11-1; end=$10} len=end-start; score=-log($5)/len;\
      if(score >= .2)print $1, start, end, $3, end-start, $9, $4, $5}' ${3%.*}.out | \
    sort -k1,1 -k2,2n -k3,3n > $3
  if [ $4 ]
  then
    mv $3 ${3}.all
    $BEDTOOLS/intersectBed -a ${3}.all -b $4 -v > $3
  fi
}

function scan_rphmm() {
  # $1: profile-HMM
  # $2: input FA
  # $3: output BED
  $RPHMM/bin/scan -v -p $NTHREAD -c $2 $1 | awk '$6=="+"' | \
    sort -k1,1 -k2,2n -k3,3n > ${3%.*}.out
  awk 'BEGIN{OFS="\t"}{split($1, f, "|"); split(f[2], ff, "[:-]");\
      chr=ff[1]; start=ff[2]; end=ff[3];\
      if(length(f)==2){strand=$6; print chr, start+$2, start+$3, $1, $3-$2, \
        strand, $7}else{strand=f[3];if(strand=="+"){print chr, start+$2, \
        start+$3, $1, $3-$2, strand, $7}else{print chr, end-$3, end-$2, $1, \
        $3-$2, strand, $7}}}' ${3%.*}.out | sort -k1,1 -k2,2n -k3,3n > $3
}

function number_monomers() {
  # $1: scan result
  cp ${1%.*}.out ${1%%.*}.final.out
  IN=${1%%.*}.final.out
  cat $IN | awk 'BEGIN{OFS="\t"}{split($1, f, "|"); split(f[2], ff, "[:-]");\
      chr=ff[1]; start=ff[2]; end=ff[3];\
      if(length(f)==2){strand=$6; print chr, start+$2, start+$3, $1, $3-$2, \
        strand, $7}else{strand=f[3];if(strand=="+"){print chr, start+$2, \
        start+$3, $1, $3-$2, strand, $7}else{print chr, end-$3, end-$2, $1, \
        $3-$2, strand, $7}}}' | sort -k1,1 -k2,2n -k3,3n -k4,4 > ${IN%.*}.data
  cat ${IN%.*}.data | awk '$6=="+"' | sort -k4,4 -k1,1 -k2,2n -k3,3n | \
    awk -v chr="%" -f $NUMBERING > ${IN%.*}.pos.5p
  cat ${IN%.*}.data | awk '$6=="-"' | sort -k4,4 -k1,1 -k3,3nr -k2,2nr | \
    awk -v chr="%" -f $NUMBERING > ${IN%.*}.neg.5p
  cat ${IN%.*}.pos.5p | sort -k4,4 -k1,1 -k3,3nr -k2,2nr | \
    awk -v chr="#" -f $NUMBERING > ${IN%.*}.pos.3p
  cat ${IN%.*}.neg.5p | sort -k4,4 -k1,1 -k2,2n -k3,3n | awk '$6=="-"' | \
    awk -v chr="#" -f $NUMBERING > ${IN%.*}.neg.3p
  cat ${IN%.*}.???.3p | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4"|"$8, \
    $5, $6, $7}' | sort -k1,1 -k2,2n -k3,3n > ${IN%.*}.bed
  
  rm ${IN%.*}.data ${IN%.*}.???.?p
}

if [ $# -lt 3 ]
then
  echo "Usage: $0 <consensus.fa> <nhmmer output in BED format> <FASTA of extended promoter sequecnes>"
  exit
else
  export -f scan_rphmm
  LEN=$(tail -n +2 $1 | awk '{print length($1)}')

  PREV_OUT=$2
  SEQ=$3
  ITR=1
  JACCARD=0
  CONVERGE=$(float_lt $JACCARD 0.99)
  while [ \( $ITR -le $MAX_ITR \) -a $CONVERGE -eq 1 ]
  do
    echo "[ITERATION $ITR]"
    NAME=${PREV_OUT%%.*}.${ITR}.20bp.rand${RAND}
    CURR_OUT=${NAME}.scan.bed
    prep_seq $1 $PREV_OUT $NAME $RAND $LEN $REF
    muscle ${NAME}.fa

    echo "[BUILDING]"
    build_rphmm ${NAME}.aln
    echo "[SCANNING]"
    scan_rphmm ${NAME}.params $SEQ $CURR_OUT

    INTER=$($BEDTOOLS/intersectBed -a $PREV_OUT -b $CURR_OUT | sort -k1,1 \
      -k2,2n -k3,3n | $BEDTOOLS/mergeBed| awk 'BEGIN{t=0}{t+=$3-$2}END{print t}')
    UNION=$(cat $PREV_OUT $CURR_OUT | cut -f 1-3 | sort -k1,1 -k2,2n -k3,3n | \
      $BEDTOOLS/mergeBed | awk 'BEGIN{t=0}{t+=$3-$2}END{print t}')
    JACCARD=$(echo "$INTER / $UNION" | bc -l)
    echo "[JACCARD=$JACCARD]"
    CONVERGE=$(float_lt $JACCARD 0.95)
    let ITR=ITR+1
    PREV_OUT=$CURR_OUT
    echo ""
  done
  number_monomers $CURR_OUT
  cp ${NAME}.params ${NAME%%.*}.final.params
fi
