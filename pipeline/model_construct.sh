#!/bin/bash
# Set the below to RepeatProfileHMM *root* directory
RPHMM=/home/rcf-40/mengzhou/panfs/repeats/RepeatProfileHMM/
source ${RPHMM}/pipeline/config.sh

function muscle_align() {
  $MUSCLE -in $1 -quiet -clw > ${1%.*}.aln
}

function trf_by_revcomp() {
  INPUT_FA=$1
  EXPECTED_START=$2
  EXTENSION=$3
  NAME=${INPUT_FA%.fa}
  TRF_PARS="2 5 5 80 10 50 500 -h -ngs"
  MONOMER_MIN_LEN=180
  MONOMER_MAX_LEN=220
  MONOMER_MIN_CN=2.0
  MONOMER_MIN_START=$(expr $EXPECTED_START - $EXTENSION)
  MONOMER_MAX_START=$(expr $EXPECTED_START + $EXTENSION)
  INPUT_REV=${NAME}.rev.fa
  $REV_COMP $INPUT_FA > $INPUT_REV
  $TRF $INPUT_REV $TRF_PARS > ${NAME}.rev.monomers.list
  cat ${NAME}.rev.monomers.list | awk 'BEGIN{OFS="\t"; name="";}\
    {if($1~/^@/){name=substr($1, 2);}else{print name, $1, $5, $4, $14}}' | \
    sort -k3,3n -k1,1 > ${NAME}.rev.monomers.data
  awk "BEGIN{flag=0;l=0}{if(substr(\$0,1,1)==\"@\")\
      {l=\$0; flag=1}else if(flag==1 && \$5>=$MONOMER_MIN_LEN && \
      \$5<=$MONOMER_MAX_LEN && \$4>=$MONOMER_MIN_CN\
      && \$1>=$MONOMER_MIN_START && \$1<=$MONOMER_MAX_START){print l; \
      print \$14; flag=0}}" ${NAME}.rev.monomers.list | sed 's/^@/>/g' | \
    awk 'BEGIN{srand(); OFS="\t"}{if($1~/^>/){name=$1}else{print rand(), name, $1}}' | \
    sort -k1,1g | cut -f 2- | awk 'NR<=500{print $1; print $2}' > ${NAME}.rev.monomers.rand500.fa
  $REV_COMP ${NAME}.rev.monomers.rand500.fa > ${NAME}.monomers.fa
  rm ${NAME}.rev.monomers.rand500.fa $INPUT_REV
}

function get_init_promoter() {
  # get potential promoter regions by merging all monomers and extend by
  # PROMOTER_EXT bp. The promoter regions are stranded for now to allow TRF.
  echo "[GENERATING EXTENDED PROMOTERS]"
  IN=$1
  NAME=${IN%.bed}
  NHMMER_EXT=$2
  PROMOTER_EXT=$3
  CHR_SIZES=$4
  REF=$5
  INIT_PROMOTER=${FAMILY}.ext${PROMOTER_EXT}
  cat $IN | $BEDTOOLS/mergeBed -d $NHMMER_EXT -s | sort -k1,1 -k2,2n -k3,3n | \
    $BEDTOOLS/bedtools slop -b $PROMOTER_EXT -g $CHR_SIZES -s | \
    awk -v family=$NAME \
    'BEGIN{OFS="\t"}{print $1, $2, $3, family"|"$1":"$2"-"$3"|"$4, $3-$2, $4;}'\
    > ${INIT_PROMOTER}.bed
  $BEDTOOLS/bedtools getfasta -s -name -fi $REF -bed ${INIT_PROMOTER}.bed -fo \
    ${INIT_PROMOTER}.fa
}

function get_init_monomer() {
  # get monomers using TRF, and construct the initial model
  echo "[CONSTRUCTING INITIAL MODEL]"
  IN=$1
  trf_by_revcomp $IN $2 50
  muscle_align ${1%.*}.monomers.fa
  $RPHMM/bin/construct -o ${1%.*}.initial.params ${1%.*}.monomers.aln
}

function scan_rpmm() {
  IN_FA=$1
  MODEL=$2
  OUT_PREFIX=$3
  echo "[SCANNING]"
  $RPHMM/bin/scan -v -p $NTHREAD -c $IN_FA -o ${OUT_PREFIX}.scan $MODEL
  cat ${OUT_PREFIX}.scan | awk 'BEGIN{OFS="\t"}{split($1, f, "|"); \
      split(f[2], ff, "[:-]"); strand=f[3]; if(strand=="+"){start=ff[2]+$2; \
      end=ff[2]+$3}else{start=ff[3]-$3; end=ff[3]-$2;} \
      print ff[1], start, end, $1, end-start, strand, $NF }' | \
    sort -k1,1 -k2,2n -k3,3n -k6,6 > ${OUT_PREFIX}.bed
}

function get_top_n_median() {
  NUM=$1
  IN_BED=$2
  CORE_LEN=$(cut -f 5 $IN_BED | sort | uniq -c | \
    awk 'BEGIN{OFS="\t"}{print $2, $1}' | sort -k2nr | head -n $NUM | \
    cut -f 1 | sort -k1n | awk -f $GET_MEDIAN)
}

function sample_by_len() {
  IN_BED=$1
  LEN=$2
  EXT=$3
  NUM=$4
  cat $IN_BED | awk -v l=$LEN -v ext=$EXT \
    'BEGIN{OFS="\t";srand()}{if($5>=l-ext && $5<=l+ext){print rand(), $0}}' | \
    sort -k1,1g | head -n $NUM | cut -f 2- | sort -k1,1 -k2,2n -k3,3n > \
    ${IN_BED%.*}.rand${NUM}.bed
  $BEDTOOLS/bedtools getfasta -s -name -fi $REF -bed ${IN_BED%.*}.rand${NUM}.bed \
    -fo ${IN_BED%.*}.rand${NUM}.fa
}

if [ $# -lt 3 ]
then
  echo "Usage: $0 <nhmmer output in BED format> <extension length for each promoter> <file for chr sizes>"
  exit
fi
if [ ! -e $1 ]
then
  echo "Input file [$1] not found!"
  exit
fi
if [ ! -e $3 ]
then
  echo "File for chr sizes: [$3] does not exist!"
  exit
fi

FAMILY=${1%.*}
get_init_promoter $1 20 $2 $3 $REF
get_init_monomer ${INIT_PROMOTER}.fa $2
# first scan using the initial model to determine monomer bounds
scan_rpmm ${INIT_PROMOTER}.fa ${INIT_PROMOTER}.initial.params ${INIT_PROMOTER}.split
# estimate the model length by getting the median of top 8 abundant lengths
get_top_n_median 8 ${INIT_PROMOTER}.split.bed
# Randomly sample monomers that are close to the core length to build
# a refined model
echo "Model lenths is estimated to be ${CORE_LEN}."
sample_by_len ${INIT_PROMOTER}.split.bed $CORE_LEN 20 500
muscle_align ${INIT_PROMOTER}.split.rand${NUM}.fa

$RPHMM/bin/construct -n -u -o ${INIT_PROMOTER}.split.rand${NUM}.params \
  ${INIT_PROMOTER}.split.rand${NUM}.aln
tail -n +2 ${INIT_PROMOTER}.split.rand${NUM}.params | head -n 1 | \
  awk '{print ">consensus"; print $2}' > consensus.fa
