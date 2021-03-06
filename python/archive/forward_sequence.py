from profileHMM import *
from optparse import OptionParser
import sys
import time

def loadSeq(file):
  """Load sequences from txt file. File format:
  seq_name	sequence(ACGT)
  Returns a dictionary of emit symbol list (0-3) for each seq.
  """
  fh = open(file,'r')
  alphabet = {"A":0, "C":1, "G":2, "T":3}
  pool = {}
  for l in fh:
    f = l.split()
    name = f[0]
    seq = f[1]
    symbol = [alphabet[seq[i]] for i in range(len(seq))]
    pool[name] = symbol

  fh.close()
  return pool

def loadALN(file):
  """Load alignment from HMMER .aln file. Will look for 
  the sequence with name consensus to mark states.
  """
  alphabet = {"A":0, "C":1, "G":2, "T":3, "-":4, "$":5}
  fhd = open(file,'r')
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

  # find consensus alignment, and set the length of model to 
  # its length with only matching
  length = len(data['consensus']) - data['consensus'].count('-')
  # N is for linear model
  N = length*3 + 3

  # mark consensus states: 0=M, 1=I, 2=D
  marked_states = []
  for i in range(len(data['consensus'])):
    if data['consensus'][i] == '-':
      marked_states.append(1)
    else:
      marked_states.append(0)

  # get state and emit in each position
  states = {}.fromkeys(data.keys())
  emits = {}.fromkeys(data.keys())
  for name in data.keys():
    seq = data[name]
    for i in range(len(seq)):
      if states[name] == None:
        states[name] = []
        emits[name] = []
      if marked_states[i] == 1:
        if seq[i] == '-':
          continue
        else:
          states[name].append(1)
      else:
        if seq[i] == '-':
          states[name].append(2)
        else:
          states[name].append(0)
      emits[name].append(alphabet[seq[i]])

    # convert states to state path, and add $ to start and end of
    # emission sequence; add B, E to state path
    states[name] = [0] + [StateIndToTransProbInd(states[name],i) \
        for i in range(len(states[name]))] + [N-1]
    emits[name] = [5] + emits[name] + [5]

  return states, emits

def StateIndToTransProbInd(state, ind):
  """Convert state index of a state sequence to its corresponding
  index in transistion matrix.
  """
  counter = 0
  # get the order of give state index, e.g. M2 or D3 or ...
  for i in range(ind+1):
    if not state[i] == 1:
      counter += 1

  if counter == 0:
    # given state is I0
    return 1
  else:
    return 2+(counter-1)*3+state[ind]

def LoadBackground(file):
  """File format (value is [0-1]):
  A	frac_a
  C	frac_c
  G	frac_g
  T	frac_t
  """
  fh = open(file, 'r')
  d = {}.fromkeys(["A","C","G","T"])
  for l in fh:
    f = l.split()
    d[f[0]] = float(f[1])

  fh.close()
  return [d[key] for key in sorted(d.keys())]

def main():
  parser = OptionParser()
  parser.add_option("-a", "--align", action="store", type="string",
                    dest="AlnFile", help="ALN file", metavar="<FILE>")
  parser.add_option("--trans", action="store", type="string",
                    dest="TransFile", help="output transistion matrix File", metavar="<FILE>")
  parser.add_option("--emit", action="store", type="string",
                    dest="EmitFile", help="output emission matrix File", metavar="<FILE>")
  parser.add_option("-l", "--likelihood", action="store", type="string",
                    dest="LhFile", help="output file. Can be stdout for screen output", metavar="<FILE>")
  parser.add_option("-i", "--input", action="store", type="string",
                    dest="InFile", help="input sequence file.", metavar="<FILE>")
  parser.add_option("-b", "--background", action="store", type="string",
                    dest="BackgroundFile", help="Background fraction file", metavar="<FILE>")
  parser.add_option("-s", "--start", type="string", default="-",
      dest="start", help="Specify the number of alignment to start from. Use - for the first one. Default: -") 
  parser.add_option("-e", "--end", type="string", default="-",
      dest="end", help="Specify the number of alignment to stop computing. Use - for the last one. Default: -") 
  parser.add_option("--local-alignment", action="store_true",
                    dest="LocalAlignment", help="Convert model to local alignment model.") 
  (opt, args) = parser.parse_args(sys.argv)

  if not opt.AlnFile or not opt.LhFile or not opt.InFile:
    sys.stderr.write("ERROR: must specify ALN file and output file!\n")
    parser.print_help()
    sys.exit(1)
  if opt.LocalAlignment:
    sys.stderr.write("Not available now.\n")
    sys.exit(1)

  hmm = ProfileHMM(4,opt.AlnFile)
  if opt.TransFile and opt.EmitFile:
    hmm.OutputTransProb(opt.TransFile)
    hmm.OutputEmitProb(opt.EmitFile)
  if opt.LocalAlignment:
    hmm.SimpleLinearModelToLocalAlignment()
  if opt.BackgroundFile:
    frac = LoadBackground(opt.BackgroundFile)
    hmm.SetNullModel(frac)

  seq = loadSeq(opt.InFile)
  if opt.LhFile == 'stdout':
    outfh = sys.stdout
  else:
    outfh = open(opt.LhFile, 'w')
  if opt.end == '-':
    end = len(seq.keys())
  else:
    end = min(int(opt.end), len(seq.keys()))
  if opt.start == '-':
    start = 0
  else:
    start = max(int(opt.start), 0)

  pool = sorted(seq.keys())
  for name in pool[start:end]:
    log_lh_no_del = hmm.ForwardLogLikelihoodNew(seq[name], True)
    outfh.write("%s\t%f\n"%(name, log_lh_no_del))

  if not outfh == sys.stdout:
    outfh.close()

if __name__ == '__main__':
  main()
