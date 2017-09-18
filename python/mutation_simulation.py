"""Simulation related classes.
"""
import sys
import numpy

from utils import *

class mutations:
  def __init__(self, sub_p = 0.01, indel_p = 0.0005, \
      trunc_p = 0.30, seed = None):
    self.sub_p = sub_p
    self.indel_p = indel_p
    self.trunc_p = trunc_p
    if seed:
      numpy.random.seed(seed)

  def substitution(self, age, seq, fixed_nt=None):
    """p is the probability of mutation per site.
    """
    subs = {}
    site_count = numpy.random.binomial(len(seq), age*self.sub_p)
    if site_count > 0:
      sites = list(numpy.random.choice(len(seq), site_count, replace=False))
      for i in sites:
        if not fixed_nt:
          subs[i] = self._point_mutation(seq[i])
        else:
          subs[i] = fixed_nt

    return subs

  def indel(self, age, seq_length):
    # the following constants are set based on results from
    # Zhang Z and Gerstein M, 2003.
    # mean 4.2, median 2 in human for both ins and del
    # the ratio of del:ins in human is around 3:1
    I_VS_D = 0.25

    # assume the number of indel instances follows a binomial distribution
    site_count = numpy.random.binomial(seq_length, age*self.indel_p)
    indels = {}
    # assume the length of an indel follows a log-normal disttribution,
    # truncated in [1, MAX_INDEL]. The mean is close to 4 in this setup, but
    # the distribution does not look similar.
    MAX_INDEL = 10
    MU = 1.2
    SIG = 0.3
    if site_count > 0:
      sites = list(numpy.random.choice(seq_length, site_count, replace=False))
      for i in sites:
        length = int(numpy.round(numpy.random.lognormal(MU, SIG)))
        while length < 1 or length > MAX_INDEL:
          length = int(numpy.round(numpy.random.lognormal(MU, SIG)))
        if numpy.random.random() > I_VS_D:
          # deletion
          indels[i] = [1, length]
        else:
          # insertion
          indels[i] = [0, self._random_seq(length)]

    return indels

  def tandem_repeat(self, seq_length, SITE_MIN = 0, SITE_MAX = 3):
    SITE_ALPHA = 2
    SITE_BETA = 7
    site_count = int(numpy.round(\
      numpy.random.beta(SITE_ALPHA, SITE_BETA)*(SITE_MAX - SITE_MIN))+SITE_MIN)
    tandems = {}
    if site_count > 0:
      sites = list(numpy.random.choice(seq_length, site_count, replace=False))
      LEN_MIN = 3
      LEN_MAX = 8
      LEN_ALPHA = 3
      LEN_BETA = 4
      for i in sites:
        length = int(numpy.round(\
            numpy.random.beta(LEN_ALPHA, LEN_BETA)*(LEN_MAX - LEN_MIN)+LEN_MIN))
        copy_number = numpy.random.randint(2, 4)
        tandems[i] = [length, copy_number]

    return tandems

  def trunc_len_beta(self, seq_length):
    ALPHA = 3
    BETA = (1-self.trunc_p)/self.trunc_p*ALPHA
    MIN = 0
    MAX = int(seq_length*0.4)
    trunc_len = int(\
      numpy.round(numpy.random.beta(ALPHA, BETA)*(MAX-MIN)+MIN))

    return trunc_len

  def trunc_len_geo(self, seq_length, lower=1):
    length = numpy.random.geometric(self.trunc_p)
    trunc_len = min(seq_length - length % seq_length + 1, seq_length - lower)

    return trunc_len

  def _point_mutation(self, nt):
    alphabet = ["a", "c", "g", "t"]
    nt = nt.lower()
    alphabet.remove(nt)
    
    return alphabet[numpy.random.choice(len(alphabet))]

  def _random_seq(self, length):
    alphabet = ["a", "c", "g", "t"]
    seq = "".join([alphabet[i] for i in numpy.random.randint(0,4,length)])
    return seq

class repeat_copy:
  def __init__(self, name, seq, rep_id, subs = {}, indels = {}, tandems = {}, \
      trunc_len = 0):
    self.name = name
    self.init_seq = seq
    self.actual_seq = self.init_seq
    self.rep_id = rep_id
    self.subs = subs
    self.indels = indels
    self.tandems = tandems
    self.trunc_len = trunc_len
    self._get_label()
    self._self_mutate()

  def __repr__(self):
    return "<repeat_copy> %s %s %s"%(self.name, self.init_seq, self.label)

  def __str__(self):
    return "\t".join([self.name, self.label, self.init_seq, self.actual_seq])

  def fasta(self, out = sys.stdout):
    out.write(">%s|%s\n"%(self.name, self.label))
    out.write(text_wrap(self.actual_seq))

  def add_mutation(self, subs, indels, tandems):
    merged_subs = merge_dict(subs, self.subs)
    merged_indels = merge_dict(indels, self.indels)
    merged_tandems = merge_dict(tandems, self.tandems)
    self.subs = merged_subs
    self.indels = merged_indels
    self.tandems = merged_tandems

    self._self_mutate()
    self._get_label()

  def external_mutation(self, subs, indels, tandems, trunc_len):
    self.actual_seq = self.init_seq
    merged_subs = merge_dict(subs, self.subs)
    merged_indels = merge_dict(indels, self.indels)
    merged_tandems = merge_dict(tandems, self.tandems)
    self._apply_indel(merged_indels)
    self._apply_tandem_repeat(merged_tandems)
    if trunc_len == 0:
      self._apply_truncation(self.trunc_len)
    else:
      self._apply_truncation(trunc_len)
    self._apply_substitution(merged_subs)

  def _self_mutate(self):
    self.actual_seq = self.init_seq
    self._apply_indel(self.indels)
    self._apply_tandem_repeat(self.tandems)
    self._apply_truncation(self.trunc_len)
    self._apply_substitution(self.subs)

  def _get_label(self):
    num_sub = len(self.subs)
    count = lambda x:sum([i[0] == x for i in self.indels.values()])
    num_insertion = count(0)
    num_deletion = count(1)
    num_tandems = len(self.tandems)

    self.label = "S%dI%dD%dR%dT%d"%(\
        num_sub, num_insertion, num_deletion, num_tandems, self.trunc_len)

  def _apply_substitution(self, subs):
    """subs = {location:str for seq}
    """
    if len(subs) == 0:
      return
    sites = sorted([i for i in subs.keys() if i < len(self.actual_seq)])
    mutated = [self.actual_seq[i] for i in xrange(len(self.actual_seq))]
    for i in sites:
      mutated[i] = subs[i]
    self.actual_seq = "".join(mutated)

  def _apply_indel(self, indels):
    """indels = {location:[0/1 for insertion/deletion, str for seq]}
    """
    if len(indels) == 0:
      return
    sites = sorted(indels.keys())
    mutated = self.actual_seq[:sites[0]]
    for i in xrange(len(sites)-1):
      if indels[sites[i]][0] == 0:
        # insertion
        mutated += indels[sites[i]][1]
        mutated += self.actual_seq[sites[i]:sites[i+1]]
      else:
        # deletion
        mutated += self.actual_seq[sites[i]+indels[sites[i]][1]:sites[i+1]]
    if indels[sites[-1]][0] == 0:
      mutated += indels[sites[-1]][1]
      mutated += self.actual_seq[sites[-1]:]
    else:
      mutated += self.actual_seq[sites[-1]+indels[sites[-1]][1]:]

    self.actual_seq = mutated

  def _apply_tandem_repeat(self, tandems):
    """tandems = {location:[int for length, int for copy number]}
    """
    if len(tandems) == 0:
      return
    sites = sorted(tandems.keys())
    mutated = self.actual_seq[:sites[0]]
    for i in xrange(len(sites)-1):
      unit = self.actual_seq[sites[i]:sites[i]+tandems[sites[i]][0]]
      mutated += unit * tandems[sites[i]][1]
      mutated += self.actual_seq[sites[i]:sites[i+1]]
    unit = self.actual_seq[sites[-1]:sites[-1]+tandems[sites[-1]][0]]
    mutated += unit * tandems[sites[-1]][1] \
      + self.actual_seq[sites[-1]+tandems[sites[-1]][0]:]

    self.actual_seq = mutated

  def _apply_truncation(self, trunc_len):
    trunc_len = int(1.0 * trunc_len / len(self.init_seq) * len(self.actual_seq))
    self.actual_seq = self.actual_seq[trunc_len:]

class repeat_family:
  def __init__(self, name, consensus, age, num,\
      mutator, subs = {}, indels = {}, tandems = {}):
    self.name = name
    self.consensus = consensus
    self.age = age
    self.num = num
    self.mutator = mutator
    self.copies = {}
    self.subs = subs
    self.indels = indels
    self.tandems = tandems
    self._get_label()

  def __repr__(self):
    return "<repeat_family:%s age:%d size:%d>"%(self.name, self.age, self.num)

  def __str__(self):
    return "\n".join([i.__str__() for i in self.copies.values()])

  def fasta(self, out = sys.stdout):
    for i in self.copies.values():
      i.fasta(out)
      #out.write(">%s\n"%(i.name))
      #out.write(text_wrap(i.actual_seq))

  def generate_copies(self, no_tandem = False, no_trunc = False):
    for i in xrange(self.num):
      subs = self.mutator.substitution(self.age, self.consensus)
      indels = self.mutator.indel(self.age, len(self.consensus))
      if no_tandem:
        tandems = {}
      else:
        tandems = self.mutator.tandem_repeat(len(self.consensus))
      if no_trunc:
        trunc_len = 0
      else:
        trunc_len = self.mutator.trunc_len_beta(len(self.consensus))
      rep_id = i
      new = repeat_copy(self._get_copy_name(i), self.consensus,\
          rep_id, subs, indels, tandems, trunc_len)
      self.copies[rep_id] = new

  def familywise_mutate(self, frac = 1.0):
    for i in self.copies.values():
      if numpy.random.rand() <= frac:
        i.external_mutation(self.subs, self.indels, self.tandems, 0)
        i.label = self.label + "|" + i.label
      else:
        i.label = "NULL|" + i.label

  def aging(self, increment):
    """Make the whole family older by adding random mutations to individual
    copies according to the age increment.
    """
    for copy_id in self.copies.keys():
      repeat = self.copies[copy_id]
      new_subs = self.mutator.substitution(increment, repeat.actual_seq)
      new_indels = self.mutator.indel(increment, len(repeat.actual_seq))
      new_tandems = {}
      repeat.add_mutation(new_subs, new_indels, new_tandems)

    self.age += increment

  def get_most_likely_mc(self, count = 1, max_trunc_bp = 10):
    """Sort all copies based on the truncation length, and sample from those
    that were not heavily truncated, i.e. the most likely master copies.
    Criteria:
      1. sorted descending by total amount of mutations from the consensus
      2. trunc length smaller than max_trunc_bp
    Note that the returned list can be empty
    """
    pool = sorted(self.copies.items(), reverse=True, \
        key=lambda x:len(x[1].subs) + len(x[1].indels) + len(x[1].tandems))
    candidates = [repeat[0] for repeat in pool \
        if repeat[1].trunc_len <= max_trunc_bp]
    if len(candidates) >= count:
      return numpy.random.choice(candidates, count, replace=False)
    elif len(candidates) > 0:
      return candidates
    else:
      return []

  def _get_label(self):
    num_sub = len(self.subs)
    count = lambda x:sum([i[0] == x for i in self.indels.values()])
    num_insertion = count(0)
    num_deletion = count(1)
    num_tandems = len(self.tandems)

    self.label = "S%dI%dD%dR%d"%(\
        num_sub, num_insertion, num_deletion, num_tandems)

  def _get_copy_name(self, idx):
    name = self.name
    width = int(numpy.log10(self.num))+1
    idx_format = "%0"+str(width)+"d"
    #return "%s_%s|%s"%(name, idx_format%idx, self.label)
    return "%s_%s"%(name, idx_format%idx)
