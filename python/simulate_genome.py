#!/usr/bin/env python
"""Simulate insertions with a genome and repeat sequences.
Can simulate some mutations while inserting.
"""
import sys
import random, numpy
from optparse import OptionParser

class repeat_family:
  def __init__( self, consensus, name = "Copy", trunc_par = 0.01, \
      trunc_lower = 20 ):
    self.consensus = consensus
    self.truncation_par = trunc_par
    self.truncation_lower = trunc_lower
    self.name = name
    
  def generate( self, size, age = 10 ):
    """Generate the sequences of the whole family.
    """
    self.copies = [{"name":"", "trunc_len":0, "mut_sites":[], "seq":""} \
      for i in xrange(size)]

    for i in xrange(size):
      trunc_seq, trunc_len = \
        self._apply_truncated(self.truncation_par, self.truncation_lower)
      self.copies[i]["trunc_len"] = len(self.consensus) - trunc_len
      self.copies[i]["mut_sites"], self.copies[i]["seq"] = \
        self._apply_spontaneous_mutation( trunc_seq, age)
      indel_count, self.copies[i]["seq"] = \
        self._apply_simple_indel(self.copies[i]["seq"],age)
      self.copies[i]["name"] =  "%s_%d_%s"%(self.name, i+1, "T%d_M%d_I%d"%\
        (self.copies[i]["trunc_len"], len(self.copies[i]["mut_sites"]), \
        indel_count))

  def generate_full_length( self, size, age = 10 ):
    """Generate sequences of full length copies (minimal truncation).
    """
    self.copies_fl = [{"name":"", "trunc_len":0, "mut_sites":[], "seq":""} \
      for i in xrange(size)]

    for i in xrange(size):
      trunc_len = 0
      trunc_seq = self.consensus
      trunc_seq, trunc_len = \
        self._apply_truncated(self.truncation_par, int(len(self.consensus)*0.95))
      self.copies_fl[i]["trunc_len"] = len(self.consensus) - trunc_len
      self.copies_fl[i]["mut_sites"], self.copies_fl[i]["seq"] = \
        self._apply_spontaneous_mutation( trunc_seq, age)
      indel_count, self.copies_fl[i]["seq"] = \
        self._apply_simple_indel(self.copies_fl[i]["seq"],age)
      self.copies_fl[i]["name"] =  "%s_%d_%s"%(self.name, i+1, "T%d_M%d_I%d"%\
        (self.copies_fl[i]["trunc_len"], len(self.copies_fl[i]["mut_sites"]), \
        indel_count))
      
  def _apply_truncated( self, p, lower):
    """Generate one repeat sequence. This simulates the process of
    L1 reverse-transcription, i.e. truncation.
    """
    # assuming the truncation follows a truncated geometric distribution
    length = numpy.random.geometric(p)
    # trunc_length is the length of remaining sequence after truncation
    trunc_length = lower + \
      length%(len(self.consensus)-lower) + 1
    truncated = self.consensus[-trunc_length:]

    return truncated, trunc_length

  def _apply_spontaneous_mutation( self, seq, age ):
    """Simulate spontaneous mutation.
    """
    p = age*0.01
    site_count = numpy.random.binomial(len(seq), p)
    sites = sorted(list(numpy.random.choice(len(seq), site_count)))
  
    mutated = [seq[i] for i in xrange(len(seq))]
    for i in sites:
      mutated[i] = self.__point_mutation(mutated[i])
  
    return sites, "".join(mutated)

  def __point_mutation( self, nt ):
    alphabet = ["a", "c", "g", "t"]
    nt = nt.lower()
    alphabet.remove(nt)
    
    return random.choice(alphabet)

  def _apply_simple_indel( self, seq, age ):
    """Simulate simple short indel mutations.
    """
    # assume the number of indel instances follows a binomial distribution
    # also assume the length of an indel follows a Poisson disttribution
    p_count = age*0.0004
    indel_length = 3
    # favoring insertion or deletion
    ins_vs_del = 0.4
    site_count = numpy.random.binomial(len(seq), p_count)
    if site_count == 0:
      return site_count, seq
    else:
      # simulate all mutations
      sites = sorted(list(numpy.random.choice(len(seq), site_count)))
      indels = {}.fromkeys(sites)
      for i in sites:
        if random.random() < ins_vs_del:
          # insertion
          indels[i] = [0, self.__poisson_seq(indel_length)]
        else:
          # deletion
          indels[i] = [1, numpy.random.poisson(indel_length)]

      # apply these mutations
      mutated = seq[:sites[0]]
      for i in xrange(len(sites)-1):
        if indels[sites[i]][0] == 0:
          mutated += indels[sites[i]][1]
          mutated += seq[sites[i]:sites[i+1]]
        else:
          mutated += seq[sites[i]+indels[sites[i]][1]:sites[i+1]]
      i = len(sites) - 1
      if indels[sites[i]][0] == 0:
        mutated += indels[sites[i]][1]
        mutated += seq[sites[i]:]
      else:
        mutated += seq[sites[i]+indels[sites[i]][1]:]

      return site_count, mutated

  def __poisson_seq( self, p = 3 ):
    length =  numpy.random.poisson(p)
    alphabet = ["a", "c", "g", "t"]
    seq = "".join([alphabet[i] for i in numpy.random.randint(0,3,length)])
    return seq

def load_fasta( inf ):
  """Load fasta into a string.
  """
  infh = open(inf, 'r')
  l = infh.readline()
  name = l[1:].strip()
  seq = ""
  for l in infh:
    seq += l.strip()
  infh.close()

  return name, seq

def generate_insertion_sites( genome_length, num ):
  """Generate insertion sites by random sampling.
  """
  return sorted(random.sample(xrange(genome_length), num))

def text_wrap( text, width = 50 ):
  if len(text) < width:
    return text + "\n"

  bin = len(text) / width
  new_text = ""
  for i in xrange(bin):
    new_text += text[i*width:(i+1)*width] + "\n"
  new_text += text[(i+1)*width:]
  if new_text[-1] != "\n":
    new_text += "\n"

  return new_text

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
  for par in ["number", "trunc_par", "trunc_lower", "age"]:
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
    default=[10], dest="number", help="Number of insertions. Default: 10.",\
    callback=parse_int_list)
  parser.add_option("-t", "--truncation", action="callback", type="string",
    dest="trunc_par", help="Parameter for truncation probability. "+\
    "Defalut: 0.005", default=[0.005], callback=parse_float_list)
  parser.add_option("-l", "--lower", action="callback", type="string",
    default=[10], dest="trunc_lower", \
    help="Lower bound of length left after truncation. Default: 10.",\
    callback=parse_int_list)
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

  for idx, infile in enumerate(opt.repeat):
    name, seq = load_fasta(infile)
    consensus_name.append(name)
    consensus_seq.append(seq)
    new_family = repeat_family(seq, name, opt.trunc_par[idx], opt.trunc_lower[idx])
    new_family.generate(opt.number[idx], opt.age[idx])
    # generate full length copies for model training
    #new_family.generate_full_length(100, opt.age[idx])
    families.append(new_family)

  insertion_sites = generate_insertion_sites(len(genome_seq), sum(opt.number))
  # get a merged family
  all_copy_seq = []
  all_copy_name = []
  all_copy_mut_sites = []
  for i in families:
    for j in i.copies:
      all_copy_seq.append(j["seq"])
      all_copy_name.append(j["name"])
      all_copy_mut_sites.append(j["mut_sites"])

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
      str(len(all_copy_mut_sites[shuffle_order[i]])), "+")) + "\n")

  annotation_fh.close()

  # output the full length sequences in FASTA for model training,
  # or output the inserted sequences
  outr = open(opt.out_repeats, 'w')
  #for i in new_family.copies_fl:
  for i in xrange(len(all_copy_seq)):
    outr.write(">" + all_copy_name[i] + "\n")
    outr.write(text_wrap(all_copy_seq[i]))
  outr.close()

if __name__ == "__main__":
  main()
