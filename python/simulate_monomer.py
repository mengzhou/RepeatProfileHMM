#!/usr/bin/env python
"""Simulate tandem repeat sequences (monomers) with 5' truncation.
"""
import sys
import random, numpy
from optparse import OptionParser
from simulate_genome import repeat_family
from simulate_genome import load_fasta
from simulate_genome import text_wrap

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
  if not opt.consensus or not opt.joint_output or not opt.sep_output:
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
  parser.add_option("-o", "--output", action="store", type="string",
    dest="joint_output", help="Output file for the simulated promoter sequence.")
  parser.add_option("-s", "--sep-output", action="store", type="string",
    dest="sep_output", help="Output file for the simulated monomer sequence, "+\
        "with monomers listed separatedly.")
  parser.add_option("-n", "--num", action="store", type="int",
    default=100, dest="number", help="Number of simulated monomers." +\
    " Default: 100")
  parser.add_option("-t", "--truncation", action="store", type="float",
    dest="trunc_par", help="Parameter for truncation probability. "+\
    "Defalut: 0.005", default=0.005)
  parser.add_option("-i", "--indel", action="store", type="float",
    dest="indel_par", help="Parameter for indel probability. "+\
    "Defalut: 0.0005", default=0.0005)
  parser.add_option("-l", "--lower", action="store", type="int",
    default=10, dest="trunc_lower", \
    help="Lower bound of length left after truncation. Default: 10.")
  parser.add_option("-e", "--age", action="store", type="float",
    dest="age", help="Parameter for family age (in Mya). "+\
    "Defalut: 10.0", default=10.0)
  (opt, args) = parser.parse_args(sys.argv)

  opt = opt_validation(parser, opt)

  monomer = repeat_family(opt.con_seq, opt.con_name,\
    opt.trunc_par, opt.trunc_lower)

  joint_outf = open(opt.joint_output, 'w')
  separate_outf = open(opt.sep_output, 'w')

  for i in xrange(opt.number):
    # the number of monomer units per monomer sequence seems to follow
    # a gamma distribution. The following parameters are computed by
    # fitting a gamma distribution for the number of monomer units of
    # L1Md_A3
    unit_count = int(round(\
      numpy.random.gamma(shape=9.4967591, scale=0.2188508)))
    monomer.generate_full_length(unit_count, opt.age, opt.indel_par)
    #monomer.copies_fl[0]["seq"] = truncate_geo(opt.trunc_par, opt.trunc_lower, \
    #  monomer.copies_fl[0]["seq"])
    monomer.copies_fl[0]["seq"] = truncate_gamma(monomer.copies_fl[0]["seq"])

    for j in xrange(unit_count):
      separate_outf.write(">%s_%d#%d|%s\n"%(opt.con_name, i+1, j+1, \
        monomer.copies_fl[j]["name"]))
      separate_outf.write(text_wrap(monomer.copies_fl[j]["seq"]))

    monomer_seq = "".join([unit["seq"] for unit in monomer.copies_fl])
    joint_outf.write(">%s_%d#%d|%s\n"%(opt.con_name, i+1, unit_count, \
      "|".join([unit["name"] for unit in monomer.copies_fl])))
    joint_outf.write(text_wrap(monomer_seq))

if __name__ == "__main__":
  main()
