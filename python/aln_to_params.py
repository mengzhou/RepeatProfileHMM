#!/usr/bin/env python
from optparse import OptionParser
import sys
from math import log

LOG_ZERO = -1000.0

def index_m(model_len, idx):
  # M_0 ~ M_L
  return idx

def index_i(model_len, idx):
  # I_0 ~ I_L-1
  return model_len + idx + 1

def index_d(model_len, idx):
  # D_1 ~ D_L
  return model_len * 2 + idx

def col_summary(data, col):
  """Generate a summary for the given column about how many
  non "-" characters there are. used for determining match or insert state
  """
  match = 0
  for name in data.keys():
    if not data[name][col] == '-':
      match += 1

  return match > len(data.keys()) / 2.0

def skip_truncation(col_states, ref_states = None):
  # ref is usually the consensus sequence. because I want to start from the
  # first homologous state (M) for each sequence, all those truncations and
  # the proceeding insertions should be skipped
  start = 0
  end = len(col_states) - 1
  while col_states[start] == 1 or (ref_states and ref_states[start] == 1):
    start += 1
  while col_states[end] == 1 or (ref_states and ref_states[end] == 1):
    end -= 1
  return start, end

def state_to_model_index(model_len, state, col):
  if state == 'M':
    return col
  elif state == 'I':
    return col + model_len + 1
  elif state == 'D':
    return col + model_len * 2
  else:
    return -1

def emission_index(nt):
  alphabet = {'A':0, 'C':1, 'G':2, 'T':3}
  return alphabet[nt]

def pseudo_count_and_normalize(counts, only_zero = False):
  PSEUDO_COUNT = min(10.0, max(sum(counts)*0.1, 0.1))
  if only_zero:
    fractions = [i if i > 0 else i + PSEUDO_COUNT for i in counts]
  else:
    fractions = [i + PSEUDO_COUNT for i in counts]

  return [i/sum(fractions) for i in fractions]

def three_prime_truncation_and_normalize(fractions, model_len, position):
  # 3' truncation probability: negatively proportional to the distance to
  # 3' end of the repeat
  assert position > 0 and position < model_len
  WEIGHT = 0.1 * position / model_len + 0.05
  new_fractions = [i * (1.0 - WEIGHT) for i in fractions]
  new_fractions.append(WEIGHT)
  assert abs(sum(new_fractions) - 1) < 1e-4

  return new_fractions

def get_consensus(data, consensus_states):
  alphabet = {'A':0, 'C':1, 'G':2, 'T':3}
  translate = {0:'A', 1:'C', 2:'G', 3:'T'}
  consensus = []
  for i, state in enumerate(consensus_states):
    if state == 1:
      counts = [0] * 4
      for name in data.keys():
        if not data[name][i] == '-':
          counts[alphabet[data[name][i]]] += 1
      consensus.append(translate[counts.index(max(counts))])

  return "".join(consensus)
  
def main():
  global LOG_ZERO
  parser = OptionParser()
  parser.add_option("-a", "--align", action="store", type="string",
                    dest="alnfile", help="CLUSTAL/MUSCLE .aln file",
                    metavar="<FILE>")
  parser.add_option("-o", "--out", action="store", type="string",
                    dest="outfile", help="output HMM parameter file",
                    metavar="<FILE>")
  (opt, args) = parser.parse_args(sys.argv)

  if not opt.alnfile or not opt.outfile:
    parser.print_help()
    sys.exit(1)

  # 1. load sequences
  fhd = open(opt.alnfile,'r')
  data = {}
  for l in fhd:
    f = l.split()
    if not len(f) == 2:
      continue
    if l.find("*") > 0:
      continue
    
    name = f[0]
    seq = f[1]
    if not data.has_key(name):
      data[name] = [seq]
    else:
      data[name].append(seq)

  fhd.close()
  for name in data.keys():
    data[name] = "".join(data[name])

  # 2. initialize parameters
  if "consensus" not in data.keys():
    # use MLE to decide the consensus states for each column
    consensus_states = [0 if col_summary(data, i) else 1 \
      for i in xrange(len(data[data.keys()[0]]))]
    model_len = consensus_states.count(0)
  else:
    model_len = len(data["consensus"]) - data["consensus"].count('-')
    # 0 for match, 1 for insertion
    consensus_states = [1 if data["consensus"][i] == '-' else 0 \
      for i in xrange(len(data["consensus"]))]

  total_size = model_len * 3 + 2;
  transition = [[0] * total_size for i in xrange(total_size)]
  emission = [[0] * 4 for i in xrange(total_size)]

  # 3. find the column indices (i for the M_i) through the consensus
  ref_start, ref_end = skip_truncation(consensus_states)
  col = 0
  col_idx = []
  for i in xrange(ref_start, ref_end+1):
    if consensus_states[i] == 0:
      col += 1
    col_idx.append(col)

  # 4. skip flanking truncations. this is done using the skip_truncation
  # function. basically what it does is locating the leftmost and rightmost
  # matching state of the sequence compared to the consensus
  for name in data.keys():
    if name == "consensus":
      continue
    current_state = []
    for i in xrange(len(data[name])):
      if consensus_states[i] == 0 and not data[name][i] == '-':
        current_state.append(0)
      else:
        current_state.append(1)

    start, end = skip_truncation(current_state, consensus_states)
    pool = []
    for i in xrange(start, end):
      if not (consensus_states[i] == 1 and data[name][i] == '-'):
        pool.append(i)
    # 5. now get the counts of transition and emission
    for i in xrange(len(pool)-1):
      state_idx = col_idx[pool[i] - ref_start]
      if consensus_states[pool[i]] == 0:
        if current_state[pool[i]] == 0:
          state = 'M'
        else:
          state = 'D'
      else:
        state = 'I'
      idx_i = state_to_model_index(model_len, state, state_idx)
      if not state == 'D':
        emission[idx_i][emission_index(data[name][pool[i]])] += 1

      j = i + 1
      next_state_idx = col_idx[pool[j] - ref_start]
      if consensus_states[pool[j]] == 0:
        if current_state[pool[j]] == 0:
          state = 'M'
        else:
          state = 'D'
      else:
        if data[name][pool[j]] == '-':
          continue
        else:
          state = 'I'
      idx_j = state_to_model_index(model_len, state, next_state_idx)
      transition[idx_i][idx_j] += 1
      if pool[j] == end and not state == 'D':
        emission[idx_j][emission_index(data[name][pool[j]])] += 1

  # 6. add pseudo counts, and convert counts to probabilities and normalize
  # to make sure the sum is 1
  # 6.1 transitions
  for i in xrange(1, model_len-1):
    # M_1~M_L-2
    transition[index_m(model_len,i)][index_i(model_len,i)],\
      transition[index_m(model_len,i)][index_m(model_len,i+1)],\
      transition[index_m(model_len,i)][index_d(model_len,i+1)],\
      transition[index_m(model_len,i)][index_d(model_len,model_len)]= \
      three_prime_truncation_and_normalize(\
        pseudo_count_and_normalize((\
        transition[index_m(model_len,i)][index_i(model_len,i)],\
        transition[index_m(model_len,i)][index_m(model_len,i+1)],\
        transition[index_m(model_len,i)][index_d(model_len,i+1)])), model_len, i)

    # I_1~I_L-2
    transition[index_i(model_len,i)][index_i(model_len,i)],\
      transition[index_i(model_len,i)][index_m(model_len,i+1)],\
      transition[index_i(model_len,i)][index_d(model_len,i+1)] = \
      pseudo_count_and_normalize((\
        transition[index_i(model_len,i)][index_i(model_len,i)],\
        transition[index_i(model_len,i)][index_m(model_len,i+1)],\
        transition[index_i(model_len,i)][index_d(model_len,i+1)]))

    if i > 1:
      # D_2~D_L-2
      transition[index_d(model_len,i)][index_i(model_len,i)],\
        transition[index_d(model_len,i)][index_m(model_len,i+1)],\
        transition[index_d(model_len,i)][index_d(model_len,i+1)] = \
        pseudo_count_and_normalize((\
          transition[index_d(model_len,i)][index_i(model_len,i)],\
          transition[index_d(model_len,i)][index_m(model_len,i+1)],\
          transition[index_d(model_len,i)][index_d(model_len,i+1)]))

  # M_L-1
  transition[index_m(model_len,model_len-1)][index_i(model_len,model_len-1)],\
    transition[index_m(model_len,model_len-1)][index_m(model_len,model_len)],\
    transition[index_m(model_len,model_len-1)][index_d(model_len,model_len)] = \
    three_prime_truncation_and_normalize(\
      pseudo_count_and_normalize((\
      transition[index_m(model_len,model_len-1)][index_i(model_len,model_len-1)],\
      transition[index_m(model_len,model_len-1)][index_m(model_len,model_len)])),\
      model_len, i)
  # I_L-1
  transition[index_i(model_len,model_len-1)][index_i(model_len,model_len-1)],\
    transition[index_i(model_len,model_len-1)][index_m(model_len,model_len)] = \
    pseudo_count_and_normalize((\
    transition[index_i(model_len,model_len-1)][index_i(model_len,model_len-1)],\
    transition[index_i(model_len,model_len-1)][index_m(model_len,model_len)]))
  # D_L-1
  transition[index_d(model_len,model_len-1)][index_i(model_len,model_len-1)],\
    transition[index_d(model_len,model_len-1)][index_m(model_len,model_len)] = \
    pseudo_count_and_normalize((\
    transition[index_d(model_len,model_len-1)][index_i(model_len,model_len-1)],\
    transition[index_d(model_len,model_len-1)][index_m(model_len,model_len)]))

  # M_L
  transition[index_m(model_len,model_len)][index_d(model_len,model_len)] = 1.0

  # M_0, D_1, I0, and D_L
  transition[index_m(model_len,0)][index_d(model_len,1)], \
    transition[index_m(model_len,0)][index_i(model_len,0)] = \
    (0.1, 0.9)
  transition[index_d(model_len,1)]\
    [index_m(model_len,1):index_m(model_len,model_len)+1] = \
    [0.9] + [0.1/(model_len-1) for i in xrange(model_len-1)]
    #pseudo_count_and_normalize(\
    #  [min(0.01*i+0.02, 0.05) for i in xrange(model_len-1)], True)
    # just some crazy weight function
  transition[index_i(model_len,0)][index_i(model_len,0)], \
    transition[index_i(model_len,0)][index_d(model_len,1)], \
    transition[index_i(model_len,0)][total_size-1] = (0.7, 0.2, 0.1)
  transition[index_d(model_len,model_len)][index_i(model_len,0)], \
    transition[index_d(model_len,model_len)][total_size-1] = (0.9, 0.1)

  # 6.2 emissions
  # human chrM
  emission[index_i(model_len,0)] = [0.308551,0.313318,0.131555,0.246575]
  for i in xrange(1, model_len+1):
    # M_1~M_L
    emission[index_m(model_len,i)] = \
      pseudo_count_and_normalize(emission[index_m(model_len,i)])
    if i < model_len:
      # I_1~I_L-1
      emission[index_i(model_len,i)] = \
        pseudo_count_and_normalize(emission[index_i(model_len,i)])

  # 7. log transformation
  for i in xrange(len(transition)):
    for j in xrange(len(transition[i])):
      if transition[i][j] - 0 < 1e-6:
        transition[i][j] = LOG_ZERO
      else:
        transition[i][j] = log(transition[i][j])
  for i in xrange(len(emission)):
    for j in xrange(len(emission[i])):
      if emission[i][j] - 0 < 1e-6:
        emission[i][j] = LOG_ZERO
      else:
        emission[i][j] = log(emission[i][j])

  # 8. output
  outfh = open(opt.outfile, 'w')
  outfh.write("# Model length: %d\n"%model_len)
  consensus_seq = get_consensus(data, consensus_states)
  consensus_seq_wrapped = "\n#".join(\
    [consensus_seq[i:i+50] for i in xrange(0,len(consensus_seq), 50)])
  outfh.write("#%s\n"%consensus_seq_wrapped)
  outfh.write("#%s\n"%",".join([str(i) for i in consensus_states]))
  for row in xrange(len(transition)):
    for i in transition[row]:
      outfh.write("%f\t"%i)
    outfh.write("\n")
  outfh.write("//\n")

  for row in xrange(len(emission)):
    for i in emission[row]:
      outfh.write("%.4f\t"%i)
    outfh.write("\n")

  outfh.close()

if __name__ == "__main__":
  main()
