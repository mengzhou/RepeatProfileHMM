#!/usr/bin/env python
"""Convert all sequences of a FASTA file to their
reverse complements.
"""
import sys

def revc( seq ):
  comp = {"A":"T","C":"G","G":"C","T":"A","a":"t","c":"g","g":"c","t":"a",\
      "N":"N","n":"n"}
  return "".join([comp[nt] for nt in seq[::-1]])

def fold(string, width = 50):
  return [string[i:i+width]+"\n" for i in range(0,len(string),width)]

def main():
  if len(sys.argv) < 2:
    sys.stdout.write("Usage: %s <input_fa>\n"%sys.argv[0])
    sys.exit(0)

  inf_seq = open(sys.argv[1], 'r')

  seq = ""
  for l in inf_seq:
    if l.startswith(">"):
      if not seq == "":
        for s in fold(revc(seq)):
          sys.stdout.write(s);
      sys.stdout.write(l)
      seq = ""
    else:
      seq += l.strip()

  for s in fold(revc(seq)):
    sys.stdout.write(s);

if __name__ == "__main__":
  main()
