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

def merge_dict(first, second):
  new = first
  for i in second.keys():
    if not i in first.keys():
      new[i] = second[i]

  return new
