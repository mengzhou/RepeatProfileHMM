#!/usr/bin/env python
"""Simulate a clade of tandem repeat sequences (monomers) with 5' truncation,
using a uni/multi-mastercopy model.
"""
import sys
import numpy
from optparse import OptionParser

from utils import *
from mutation_simulation import *

def uni_ancestor(clade, mutator, active_list, generation, \
    age_increment, descendant_count):
  """Sample from the current family as the ancestor, and generate descendants
  assuming the only ancestor gave rise to all descendants.
  """
  clade_size = len(clade)
  family_idx = 0
  for active_idx in active_list:
    active_family = clade[active_idx]
    ancestor_id = numpy.random.choice(active_family.copies.keys())
    ancestor = active_family.copies[ancestor_id]
    ancestor.name = ancestor.name.split("|")[0] \
        + "|master_%s"%active_family.name
    new_cons = ancestor.actual_seq
    new_name = "F%d.%d"%(generation, family_idx)
    new_family = repeat_family(new_name, new_cons, age_increment, \
        descendant_count, mutator, {}, {}, {})
    new_family.generate_copies(True, True)
    clade.append(new_family)
    family_idx += 1

  for idx in xrange(clade_size-len(active_list)):
    old_family = clade[idx]
    for copy_id in old_family.copies.keys():
      cons = old_family.copies[copy_id].init_seq
      subs = mutator.substitution(age_increment, cons)
      indels = mutator.indel(age_increment, len(cons))
      tandems = {}
      trunc_len = 0
      old_family.copies[copy_id].additional_mutation(subs, indels,\
          tandems, trunc_len)

  return clade

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
  parser.add_option("-n", "--num", action="store", type="int",
    default=5, dest="number", help="Number of simulated hierarchical groups." +\
    " Default: 5")
  parser.add_option("-e", "--age-step", action="store", type="float",\
    dest="step", help="Parameter for age difference between families (in Mya). "+\
    "Defalut: 2.0", default=2.0)
  (opt, args) = parser.parse_args(sys.argv)

  opt = opt_validation(parser, opt)

  #joint_outf = open(opt.joint_output, 'w')
  separate_outf = open(opt.sep_output, 'w')
  #mut_outf = open(opt.mut_info, 'w')

  con_len = len(opt.con_seq)
  size = 50
  mutator = mutations(opt.sub_par, opt.indel_par, opt.trunc_par)
  init_family = repeat_family("F0", opt.con_seq, opt.step,\
        size, mutator, {}, {}, {})
  init_family.generate_copies(True, True)

  clade = [init_family]
  curr_active = [0]
  for generation in xrange(1, opt.number):
    # only the latest family is active, i.e. subject to sampling of ancenstor
    clade = uni_ancestor(clade, mutator, curr_active, generation, opt.step, size)
    curr_active = [generation-1, generation]

  for family in clade:
    family.fasta(separate_outf)

  separate_outf.close()

if __name__ == "__main__":
  main()
