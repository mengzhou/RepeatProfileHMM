#!/usr/bin/env python
"""Simulate hypothetical expansion models of mouse L1 monomer arrays.
The first model assumes a single template; the second one assumes
multiple template.
"""
import sys,copy
from optparse import OptionParser

from utils import *
from mutation_simulation import *

def opt_validation(parser, opt):
  if not opt.consensus or not opt.output_single or not opt.output_multi:
    parser.print_help()
    sys.exit(0)

  if opt.consensus.endswith(".fa"):
    opt.con_name, opt.con_seq = load_fasta(opt.consensus)
  else:
    opt.con_name = "NULL"
    opt.con_seq = opt.consensus

  opt.outfh_single = open(opt.output_single, 'w')
  opt.outfh_multi = open(opt.output_multi, 'w')

  return opt

def main():
  usage = "Usage: %prog"
  parser = OptionParser(usage=usage)
  parser.add_option("-c", "--consensus", action="store", type="string", \
    dest="consensus", help="The consensus sequence of the monomer." +\
    " Can be a FASTA file or a string.")
  parser.add_option("-1", "--output-1", action="store", type="string",
    dest="output_single", help="Output file for the single template model.")
  parser.add_option("-2", "--output-2", action="store", type="string",
    dest="output_multi", help="Output file for the multiple template model.")
  parser.add_option("-n", "--num", action="store", type="int",
    default=100, dest="number", help="Number of simulated promoters." +\
    " Default: 100")
  parser.add_option("-s", "--substitution", action="store", type="float",
    dest="sub_par", help="Parameter for substitution probability. "+\
    "Defalut: 0.02", default=0.02)
  parser.add_option("-i", "--indel", action="store", type="float",
    dest="indel_par", help="Parameter for indel probability. "+\
    "Defalut: 0.0008", default=0.0008)
  parser.add_option("-e", "--age", action="store", type="float",\
    dest="age", help="Parameter for separation per new insertion (in Mya). "+\
    "Defalut: 0.02", default=0.02)
  (opt, args) = parser.parse_args(sys.argv)

  opt = opt_validation(parser, opt)

  mutator = mutations(opt.sub_par, opt.indel_par, 0)

  # initial copy: 10 monomers with age = 10
  init_promoter = repeat_family("R0", opt.con_seq, 0.0, \
    10, mutator, {}, {}, {})
  init_promoter.generate_copies(True, True)

  clade = [init_promoter]
  for i in xrange(opt.number-1):
    for promoter in clade:
      promoter.aging(opt.age)
    template = clade[-1].copies[0].actual_seq
    new_promoter = repeat_family("R%d"%(i+1), template, 0.0, 10, \
      mutator, {}, {}, {})
    new_promoter.generate_copies(True, True)
    clade.append(new_promoter)

  for promoter in clade:
    promoter.fasta(opt.outfh_single)

  init_promoter = repeat_family("R0", opt.con_seq, 0.0, \
    10, mutator, {}, {}, {})
  init_promoter.generate_copies(True, True)
  clade = [init_promoter]
  for i in xrange(opt.number-1):
    for promoter in clade:
      promoter.aging(opt.age)
    new_promoter = copy.deepcopy(clade[-1])
    new_promoter.name = "R%d"%(i+1)
    new_promoter._get_label()
    for idx in xrange(len(new_promoter.copies)):
      monomer = new_promoter.copies[idx]
      monomer.name = new_promoter._get_copy_name(idx)
    clade.append(new_promoter)

  for promoter in clade:
    promoter.fasta(opt.outfh_multi)

if __name__ == "__main__":
  main()
