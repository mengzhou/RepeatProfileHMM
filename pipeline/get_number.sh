#!/bin/bash
# Input format:
# <Sequence name with coordinates> <start> <end> <COPY> <score> <strand> <state bits>
# assuming input not stranded; therefore everything is on the forward strand, and
# the strand of the orignial interval needs to be inferred from the name
if [ $# -lt 1 ]
then
  echo "Usage: $0 <scan output>"
  exit 0
fi

IN=$1
NUMBERING=/home/rcf-40/mengzhou/panfs/repeats/RepeatProfileHMM/pipeline/numbering.awk
cat $IN | awk 'BEGIN{OFS="\t"}{split($1, f, "|"); split(f[2], ff, "[:-]");\
    chr=ff[1]; start=ff[2]; end=ff[3];\
    if(length(f)==2){strand=$6; print chr, start+$2, start+$3, $1, $3-$2, strand, $7}\
    else{strand=f[3];if(strand=="+"){print chr, start+$2, start+$3, $1, $3-$2, strand, $7}\
    else{print chr, end-$3, end-$2, $1, $3-$2, strand, $7}}}' |\
  sort -k1,1 -k2,2n -k3,3n -k4,4 > ${IN%.*}.data
cat ${IN%.*}.data | awk '$6=="+"' | sort -k4,4 -k1,1 -k2,2n -k3,3n | awk -v chr="%" -f $NUMBERING > ${IN%.*}.pos.5p
cat ${IN%.*}.data | awk '$6=="-"' | sort -k4,4 -k1,1 -k3,3nr -k2,2nr | awk -v chr="%" -f $NUMBERING > ${IN%.*}.neg.5p
cat ${IN%.*}.pos.5p | sort -k4,4 -k1,1 -k3,3nr -k2,2nr | awk -v chr="#" -f $NUMBERING > ${IN%.*}.pos.3p
cat ${IN%.*}.neg.5p | sort -k4,4 -k1,1 -k2,2n -k3,3n | awk '$6=="-"' | awk -v chr="#" -f $NUMBERING > ${IN%.*}.neg.3p
cat ${IN%.*}.???.3p | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4"|"$8, $5, $6, $7}' | sort -k1,1 -k2,2n -k3,3n > ${IN%.*}.bed

rm ${IN%.*}.data ${IN%.*}.???.?p
