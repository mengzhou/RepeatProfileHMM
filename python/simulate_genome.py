#!/usr/bin/env python
"""Simulate insertions with a genome and repeat sequences.
Can simulate some mutations while inserting.
"""
import sys
import random
from optparse import OptionParser

from utils import *
from mutation_simulation import *

def generate_insertion_sites( genome_length, num ):
  """Generate insertion sites by random sampling.
  """
  return sorted(random.sample(xrange(genome_length), num))

def parse_str_list(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

def parse_int_list(option, opt, value, parser):
  setattr(parser.values, option.dest, [int(i) for i in value.split(',')])

def parse_float_list(option, opt, value, parser):
  setattr(parser.values, option.dest, [float(i) for i in value.split(',')])

def opt_validation(parser, opt):
  if not opt.genome or not opt.repeat or not opt.out_genome \
    or not opt.annotation or not opt.out_repeats:
    parser.print_help()
    sys.exit(0)

  family_count = len(opt.repeat)
  for par in ["number", "trunc_par", "sub_par", "indel_par", "age"]:
    if len(getattr(opt,par)) < family_count:
      setattr(opt, par, getattr(opt, par) * family_count)

  return opt

def main():
  usage = "Usage: %prog"
  parser = OptionParser(usage=usage)
  parser.add_option("-g", "--genome", action="store", type="string", \
    dest="genome", help="Reference genome, a FASTA file.")
  parser.add_option("-r", "--repeat", action="callback", type="string",
    dest="repeat", help="Repeat consensus, a FASTA file or a comma " +\
    "separated list.", callback=parse_str_list)
  parser.add_option("-o", "--output-genome", action="store", type="string",
    dest="out_genome", help="Output path for the new genome sequence.")
  parser.add_option("-a", "--annotation", action="store", type="string",
    dest="annotation", help="Output path for annotation of inserted repeats.")
  parser.add_option("-s", "--output-repeats", action="store", type="string",
    dest="out_repeats", help="Output path for sequence of inserted repeats.")
  parser.add_option("-n", "--num", action="callback", type="string",
    default=[50], dest="number", help="Number of insertions. Default: 50.",\
    callback=parse_int_list)
  parser.add_option("-t", "--truncation", action="callback", type="string",
    dest="trunc_par", help="Parameter for truncation probability. "+\
    "Defalut: 0.35", default=[0.35], callback=parse_float_list)
  parser.add_option("-m", "--substitution", action="callback", type="float",
    dest="sub_par", help="Parameter for substitution probability. "+\
    "Defalut: 0.02", default=[0.02], callback=parse_float_list)
  parser.add_option("-i", "--indel", action="callback", type="float",
    dest="indel_par", help="Parameter for indel probability. "+\
    "Defalut: 0.0008", default=[0.0008], callback=parse_float_list)
  parser.add_option("-e", "--age", action="callback", type="string",
    dest="age", help="Parameter for family age (in Mya). "+\
    "Defalut: 10.0", default=[10.0], callback=parse_float_list)
  (opt, args) = parser.parse_args(sys.argv)

  opt = opt_validation(parser, opt)

  genome_name, genome_seq = load_fasta(opt.genome)
  genome_name = opt.out_genome.rstrip(".fa")
  consensus_name = []
  consensus_seq = []
  families = []
  mutators = [mutations(opt.sub_par[i], opt.indel_par[i], opt.trunc_par[i]) \
      for i in range(len(opt.repeat))]

  for idx, infile in enumerate(opt.repeat):
    name, seq = load_fasta(infile)
    consensus_name.append(name)
    consensus_seq.append(seq)
    new_family = repeat_family(name, seq, opt.age[idx], opt.number[idx], \
        mutators[idx])
    new_family.generate_copies(True, False)
    families.append(new_family)

  insertion_sites = generate_insertion_sites(len(genome_seq), sum(opt.number))
  # get a merged family
  all_copy_seq = []
  all_copy_name = []
  all_copy_mut_sites = []
  for i in families:
    for j in i.copies:
      r = i.copies[j]
      all_copy_seq.append(r.actual_seq)
      all_copy_name.append(r.name+"_"+r.label)
      all_copy_mut_sites.append(len(r.subs)+len(r.indels))

  shuffle_order = range(len(all_copy_seq))
  r = random.shuffle(shuffle_order)

  new_genome = [genome_seq[:insertion_sites[0]]]
  for i in xrange(len(shuffle_order)-1):
    new_genome.append(all_copy_seq[shuffle_order[i]])
    new_genome.append(genome_seq[insertion_sites[i]:insertion_sites[i+1]])
  new_genome.append(all_copy_seq[shuffle_order[-1]])
  new_genome.append(genome_seq[insertion_sites[-1]:])
  new_genome_seq = text_wrap("".join(new_genome))

  # output the new genome
  outg = open(opt.out_genome, 'w')
  outg.write(">" + genome_name + "\n")
  outg.write(new_genome_seq)
  outg.close()

  # output the annotation of repeat insertion
  annotation_fh = open(opt.annotation, 'w')
  extra_length = [len(all_copy_seq[shuffle_order[i]]) \
    for i in xrange(len(shuffle_order))]
  for i in xrange(len(shuffle_order)):
    start = insertion_sites[i] + sum(extra_length[:i])
    end = start + len(all_copy_seq[shuffle_order[i]])
    annotation_fh.write("\t".join((genome_name, str(start), str(end), \
      all_copy_name[shuffle_order[i]], \
      str(all_copy_mut_sites[shuffle_order[i]]), "+")) + "\n")

  annotation_fh.close()

  # output the full length sequences in FASTA for model training,
  # or output the inserted sequences
  outr = open(opt.out_repeats, 'w')
  for i in xrange(len(all_copy_seq)):
    outr.write(">" + all_copy_name[i] + "\n")
    outr.write(text_wrap(all_copy_seq[i]))
  outr.close()

if __name__ == "__main__":
  main()
