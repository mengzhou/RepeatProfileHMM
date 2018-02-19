#!/bin/awk
BEGIN {
  OFS="\t";
  counter=0;
  prev="";
}
{
  if ($4==prev) ++counter;
  else counter=1;
  print $1, $2, $3, $4, $5, $6, $7, $8""chr""counter;
  prev=$4;
}
