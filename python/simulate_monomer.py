#!/usr/bin/env python
"""Simulate tandem repeat sequences (monomers) with 5' truncation.
"""
import sys
import numpy
from optparse import OptionParser

from utils import *
from mutation_simulation import *

def rand_padding(length=100):
  """Generate random sequence at a given length.
  """
  alphabet = ["a", "g", "c", "t"]
  cdf = (0.33, 0.54, 0.74, 1)
  seq = ""
  for x in xrange(length):
    p = numpy.random.random()
    nt = next(i[0] for i in enumerate(cdf) if i[1] > p)
    seq = seq + alphabet[nt]

  return seq

def truncate_geo(p, lower, seq):
  """Truncate the input sequence from 5' end based on a geometric
  distribution.
  """
  length = numpy.random.geometric(p)
  trunc_length = lower + \
    length%(len(seq)-lower) + 1

  return seq[-trunc_length:]

def truncate_gamma(seq, shape=12.374916, scale=6.749954):
  """Truncate the input sequence from 5' end based on a gamma
  distribution fitted from L1Md_A1.
  """
  trunc_length = int(numpy.round(\
    numpy.random.gamma(shape=shape, scale=scale)))
  return seq[-trunc_length:]

def opt_validation(parser, opt):
  if not opt.consensus or not opt.joint_output or not opt.sep_output \
      or not opt.mut_info:
    parser.print_help()
    sys.exit(0)

  if opt.consensus.endswith(".fa"):
    opt.con_name, opt.con_seq = load_fasta(opt.consensus)
  else:
    opt.con_name = "NULL"
    opt.con_seq = opt.consensus

  return opt

def main():
  usage = "Usage: %prog"
  parser = OptionParser(usage=usage)
  parser.add_option("-c", "--consensus", action="store", type="string", \
    dest="consensus", help="The consensus sequence of the monomer." +\
    " Can be a FASTA file or a string.")
  parser.add_option("-j", "--joint-output", action="store", type="string",
    dest="joint_output", help="Output file for the simulated promoter sequence.")
  parser.add_option("-m", "--monomer-output", action="store", type="string",
    dest="sep_output", help="Output file for the simulated monomer sequence, "+\
        "with monomers listed separatedly.")
  parser.add_option("-u", "--mut-info", action="store", type="string",
    dest="mut_info", help="Output file for the simulated mutation information.")
  parser.add_option("-n", "--num", action="store", type="int",
    default=100, dest="number", help="Number of simulated monomers." +\
    " Default: 100")
  parser.add_option("-s", "--substitution", action="store", type="float",
    dest="sub_par", help="Parameter for substitution probability. "+\
    "Defalut: 0.02", default=0.02)
  parser.add_option("-i", "--indel", action="store", type="float",
    dest="indel_par", help="Parameter for indel probability. "+\
    "Defalut: 0.0008", default=0.0008)
  parser.add_option("-t", "--truncation", action="store", type="float",
    dest="trunc_par", help="Parameter for average truncation proportion. "+\
    "Defalut: 0.35", default=0.35)
  parser.add_option("-f", "--frac", action="store", type="float",\
    dest="frac", help="Parameter for fraction of repeats to apply familywise"+\
    " mutation. Default: 1.0 (all)", default=1.0)
  parser.add_option("-e", "--age", action="store", type="float",\
    dest="age", help="Parameter for family age (in Mya). "+\
    "Defalut: 10.0", default=10.0)
  (opt, args) = parser.parse_args(sys.argv)

  opt = opt_validation(parser, opt)

  joint_outf = open(opt.joint_output, 'w')
  separate_outf = open(opt.sep_output, 'w')
  mut_outf = open(opt.mut_info, 'w')

  unit_counts = []
  for i in xrange(opt.number):
    # the number of monomer units per monomer sequence seems to follow
    # a gamma distribution. The following parameters are computed by
    # fitting a gamma distribution for the number of monomer units of
    # L1Md_A3
    unit_counts.append(int(round(\
      numpy.random.gamma(shape=9.4967591, scale=0.2188508)))+3)

  max_count = max(unit_counts)
  con_len = len(opt.con_seq)
  repeat_families = []
  mutator = mutations(opt.sub_par, opt.indel_par, opt.trunc_par)
  for i in range(max_count):
    # it seems the monomers of L1Md_A are grouped by their order starting
    # from the 3' end: the "right-prime-most" monomers share the same
    # mutation patterns, and the second-to-rightmost ones, etc. Therefore
    # the number of the simulated families is the maximum unit count
    subs = mutator.substitution(opt.age, opt.con_seq)
    indels = mutator.indel(opt.age, con_len)
    tandems = mutator.tandem_repeat(con_len)
    rfamily = repeat_family("M%d"%i, opt.con_seq, opt.age,\
        unit_counts.count(i), mutator, subs, indels, tandems)
    repeat_families.append(rfamily)

  for repeat_idx,count in enumerate(unit_counts):
    for monomer_idx in range(count):
      subs = mutator.substitution(opt.age, opt.con_seq)
      indels = mutator.indel(opt.age, con_len)
      tandems = {}
      if monomer_idx == count - 1:
        trunc_len = mutator.trunc_len_geo(con_len)
      else:
        trunc_len = 0
      monomer_name = "%s_%d%%%d#%d"%(opt.con_name, repeat_idx, 
          count - monomer_idx, monomer_idx+1)
      monomer_unit = repeat_copy(monomer_name, opt.con_seq, repeat_idx,
          subs, indels, tandems, trunc_len)
      repeat_families[monomer_idx].copies[repeat_idx] = monomer_unit

  for family in repeat_families:
    family.familywise_mutate(frac=opt.frac)
    family.fasta(separate_outf)

  header = "\t".join(("#Seq_name", "Family_ID", "Family_label", \
    "Copy_ID", "Copy_label")) + "\n"
  mut_outf.write(header)
  for repeat_idx,count in enumerate(unit_counts):
    repeat_name = "%s_%d"%(opt.con_name, repeat_idx)
    repeat_seq = ""
    for i in range(count):
      fam_label, copy_label = \
          repeat_families[count-i-1].copies[repeat_idx].label.split("|")
      mut_outf.write("\t".join((repeat_name, str(count-i-1), \
          fam_label.split("_")[-1], str(repeat_idx), \
          copy_label.split("_")[-1])) + "\n")
      # i=0 is the 3'-most monomer, therefore the last one;
      # then the last repeat family is the 5'-most one.
      repeat_seq = repeat_seq + \
          repeat_families[count-i-1].copies[repeat_idx].actual_seq
      # 5' truncation here
    joint_outf.write(">%s\n"%repeat_name)
    joint_outf.write(text_wrap(rand_padding() + repeat_seq + rand_padding()))

  joint_outf.close()
  separate_outf.close()
  mut_outf.close()

if __name__ == "__main__":
  main()
