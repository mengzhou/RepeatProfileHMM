import numpy as np
import sys

class ProfileHMM(object):
  """Simulate model based on given length.
  Type 1 for linear. Option is (len);
  type 2 for cyclic. Option is (cyclicLength, linearLength);
  type 3 for file input. Option is (file handler for transProb, fhd for emitProb);
  type 4 for HMMER .aln file input. Option is file handler for aln file.
  """
  def __init__(self, type, option, seed = None):
    # ACGT-$
    self.K = 6
    self.type = type
    if type == 1:
      self.linearLength = option[0]
      (self.transProb, self.emitProb) = \
          self.SimpleLinearModel(self.linearLength, seed = seed)
    elif type == 2:
      self.cyclicLength = option[0]
      self.linearLength = option[1]
      (self.transProb, self.emitProb) = \
          self.SimpleCyclicModel(self.cyclicLength, self.linearLength, seed)
    elif type == 3:
      self.__loadTransProbFromFile(option[0])
      self.__loadEmitProbFromFile(option[1])
    elif type == 4:
      self.__loadALN(option)

    self.N = self.transProb.shape[0]
    self.initialProb = self.transProb[0,:]
    if self.type == 2:
      self.transProbDecoding = self._NonEmittingStateReductionSimpleCyclic()
    else:
      self.transProbDecoding = self.transProb
    self.logInitialProb = np.log(self.initialProb)
    self.logTransProb = np.log(self.transProb)
    self.logEmitProb = np.log(self.emitProb)

    self.__initiateStates()

  def SimpleCyclicModel(self, cyclicLength, linearLength, seed = None):
    """Generate a simple cyclic profile HMM model by setting the 1st linear model
    as cyclic.
    """
    np.random.seed(seed)
    (seed1, seed2) = np.random.randint(0,65535,2)
    cycleRandomPar = (1,2)
    # len*(M, I, D)+2*(B, I0, E)
    totalLength = (self.cyclicLength + self.linearLength) * 3 + 2 * 3

    (cyclicTransProb, cyclicEmitProb) = \
        self.SimpleLinearModel(cyclicLength, (1,5,1), (1,8,8,1), seed1)
    (linearTransProb, linearEmitProb) = \
        self.SimpleLinearModel(linearLength, (1,5,1), (8,1,1,8), seed2)

    ####################Transistion matrix####################
    transProb = np.zeros((totalLength, totalLength))
    cyclicStateSize = cyclicTransProb.shape[0]
    transProb[:cyclicStateSize, :cyclicStateSize] = cyclicTransProb
    transProb[cyclicStateSize:, cyclicStateSize:] = linearTransProb
    # update E state of cyclic part, adding the path to linear B
    (cycle1, cycle2) = np.random.dirichlet(cycleRandomPar)
    transProb[cyclicStateSize-1, 0] = cycle1
    transProb[cyclicStateSize-1, cyclicStateSize-1] = 0
    transProb[cyclicStateSize-1, cyclicStateSize] = cycle2
    ####################Transistion matrix####################

    ####################Emission matrix####################
    emitProb = np.zeros((totalLength, self.K))
    emitProb[:cyclicStateSize, :] = cyclicEmitProb
    emitProb[cyclicStateSize:, :] = linearEmitProb
    ####################Emission matrix####################

    return (transProb, emitProb)


  def SimpleLinearModel(self, length, transRandomPar = (1,3,1), \
      emitRandomPar = (2,1,1,2), seed = None):
    """Generate simple random transition and emission matrices.
    """
    np.random.seed(seed)
    # Dirichlet parameters for transition simulation
    #transRandomPar = (1,3,1)

    ####################Transistion matrix####################
    # length is the number of M states. So the total size is B + I0 + L*3 + E.
    transProb = np.zeros((length*3+3, length*3+3))

    # transition prob of B to I0, M1, D1
    prob = np.random.dirichlet(transRandomPar)
    transProb[0,1] = prob[0]
    transProb[0,2] = prob[1]
    transProb[0,4] = prob[2]
    # transition prob of I0 to I0, M1, D1
    prob = np.random.dirichlet(transRandomPar)
    transProb[1,1] = prob[0]
    transProb[1,2] = prob[1]
    transProb[1,4] = prob[2]
    # transition prob of M_i to I_i, M_i+1, D_i+1. i=1 to L-1
    for i in xrange(0, length-1):
      prob = np.random.dirichlet(transRandomPar)
      transProb[i*3+2,i*3+3] = prob[0]
      transProb[i*3+2,(i+1)*3+2] = prob[1]
      transProb[i*3+2,(i+1)*3+4] = prob[2]
    # transition prob of I_i to I_i, M_i+1, D_i+1. i=1 to L-1
    for i in xrange(0, length-1):
      prob = np.random.dirichlet(transRandomPar)
      transProb[i*3+3,i*3+3] = prob[0]
      transProb[i*3+3,(i+1)*3+2] = prob[1]
      transProb[i*3+3,(i+1)*3+4] = prob[2]
    # transition prob of D_i to I_i, M_i+1, D_i+1. i=1 to L-1
    for i in xrange(0, length-1):
      prob = np.random.dirichlet(transRandomPar)
      transProb[i*3+4,i*3+3] = prob[0]
      transProb[i*3+4,(i+1)*3+2] = prob[1]
      transProb[i*3+4,(i+1)*3+4] = prob[2]
    # transistion prob of M_L to I_L, E
    prob = np.random.dirichlet(transRandomPar[:2])
    transProb[(length-1)*3+2,(length-1)*3+3] = prob[0]
    transProb[(length-1)*3+2,(length-1)*3+5] = prob[1]
    # transistion prob of I_L to I_L, E
    prob = np.random.dirichlet(transRandomPar[:2])
    transProb[(length-1)*3+3,(length-1)*3+3] = prob[0]
    transProb[(length-1)*3+3,(length-1)*3+5] = prob[1]
    # transistion prob of D_L to I_L, E
    prob = np.random.dirichlet(transRandomPar[:2])
    transProb[(length-1)*3+4,(length-1)*3+3] = prob[0]
    transProb[(length-1)*3+4,(length-1)*3+5] = prob[1]
    # prob of E to E = 0
    transProb[(length-1)*3+5,(length-1)*3+5] = 0
    ####################Transistion matrix####################

    ####################Emission matrix####################
    # Observation alphabet size K = 6 for ACGT-$
    emitProb = np.zeros((length*3+3, self.K))
    #emitRandomPar = (2,1,1,2)
    # B and E emits $
    emitProb[0, self.K-1] = 1
    emitProb[emitProb.shape[0]-1, self.K-1] = 1
    # I_0
    emitProb[1, :self.K-2] = np.random.dirichlet(emitRandomPar)
    for i in xrange(0,length):
      # M_i, I_i, D_i
      emitProb[i*3+2, :self.K-2] = np.random.dirichlet(emitRandomPar)
      emitProb[i*3+3, :self.K-2] = np.random.dirichlet(emitRandomPar)
      emitProb[i*3+4, self.K-2] = 1
    ####################Emission matrix####################

    return (transProb, emitProb)

  def EmitSequence(self):
    """Emit a sequence according to trans/emitProb matrix. Also output hidden
    states.
    Format: (obsSequence, hiddenStates)
    """
    obsSeq = [5]
    stateSeq = [0]
    state = 0
    while state < self.transProb.shape[1]-1:
      state = self.__choice(self.transProb[state,:])
      if abs(self.emitProb[state,:].sum() - 1) < 1e-5:
        stateSeq.append(state)
        obs = self.__choice(self.emitProb[state,:])
        obsSeq.append(obs)

    return (obsSeq,stateSeq)

  def TranslateEmitToDNA(self, emitSeq):
    """Translate emisstion states according to alphabet (nucleotides).
    """
    alphabet = ["A","C","G","T","-","$"]
    return "".join([alphabet[i] for i in emitSeq])
  
  def TranslateStatesSimpleLinearModel(self, transSeq):
    """Translate hidden states according to alphabet (M, I, D).
    """
    alphabet = ["M","I","D"]
    output = []
    for i in xrange(len(transSeq)):
      j = transSeq[i]
      if j == self.transProb.shape[0] - 1:
        output.append("E")
      elif j == 0:
        output.append("B")
      elif j == 1:
        output.append("I0")
      else:
        state = alphabet[(j-2)%3]
        number = (j-2)/3 + 1
        output.append(state + str(number))
  
    return " ".join(output)
  
  def TranslateStatesLocalAlignment(self, transSeq):
    """Translate hidden states according to alphabet (M, I, D).
    """
    alphabet = ["M","I","D"]
    output = []
    for i,j in enumerate(transSeq):
      if j == self.N - 1:
        output.append("E_A")
      elif j == 0:
        output.append("B_A")
      elif j == 1:
        output.append("G0")
      elif j == 2:
        output.append("N0")
      elif j == 3:
        output.append("B")
      elif j == 4:
        output.append("I0")
      elif j == self.N - 2:
        output.append("G1")
      elif j == self.N - 3:
        output.append("N1")
      elif j == self.N - 4:
        output.append("E")
      else:
        state = alphabet[(j-5)%3]
        number = (j-5)/3 + 1
        output.append(state + str(number))
  
    return " ".join(output)

  def TranslateStatesSimpleCyclicModel(self, transSeq):
    """Translate hidden states according to alphabet:
    (M, I, D) for linear;
    (N, U, S) for cyclic.
    """
    cyclicAlphabet = ["N","U","S"]
    linearAlphabet = ["M","I","D"]
    output = []
    for i in xrange(len(transSeq)):
      j = transSeq[i]
      if j == self.transProb.shape[0] - 1:
        output.append("E")
      elif j == 0:
        output.append("B")
      elif j == 1:
        output.append("U0")
      elif j == self.cyclicLength*3+4:
        output.append("I0")
      elif j >= 2 and j <= self.cyclicLength*3+1:
        state = cyclicAlphabet[(j-2)%3]
        number = (j-2)/3 + 1
        output.append(state + str(number))
      else:
        state = linearAlphabet[(j+1)%3]
        number = (j-5-self.cyclicLength*3)/3 + 1
        output.append(state + str(number))
  
    return " ".join(output)

  def SetNullModel(self, background = [.25,.25,.25,.25]):
    """Set null model nucleotide fraction.
    """
    log_bg = np.log(background)
    for i in self.NonDeletionEmitStates:
      self.logEmitProb[i,:4] = self.__replaceNaN(self.logEmitProb[i,:4] -\
          log_bg)

    self.emitProb = np.exp(self.logEmitProb)
    
  def BaumWelchTraining(self, data, numsteps = 20, verbose = False):
    """Train model using Baum-Welch algorithm.
    Data is expected to be in one-hot encoding. 
    """
    #learn on a single sequence
    if type(data) is not type([]):
      numclasses, numpoints = data.shape
      assert numclasses == self.K
      logXi = \
          np.zeros((numpoints-1,self.N,self.N),dtype=float)
      logGamma = np.zeros((numpoints, self.N), dtype = float)

      lastlogprob = -np.inf
      for iteration in range(numsteps):
        if verbose:
          sys.stderr.write("\rIteration: %d/%d"%(iteration+1, numsteps))
          sys.stderr.flush()
        #E-step:
        logAlpha, logBeta = self.__logForwardBackward(data)
        #compute xi
        for t in range(numpoints-1):
          for i in range(self.N):
            for j in range(self.N):
              logXi[t, i, j] = logAlpha[i, t]+\
                  self.logTransProb[i, j]+\
                  self.logEmitProb[j,np.where(data[:,t+1]==1)[0][0]]+\
                  logBeta[j, t+1]
                  #self.__replaceNaN(self.logEmitProb[j,:]*data[:,t+1],0).sum()+\
          logXi[t, :, :] -= self.__logSumExp(logXi[t, :, :].flatten())
        logXi = self.__replaceNaN(logXi)

        #compute gamma
        for t in range(numpoints):
          for i in range(self.N):
            logGamma[t, i] = logAlpha[i,t] + logBeta[i,t]
          logGamma[t,:] -= self.__logSumExp(logGamma[t,:].flatten())
          logGamma = self.__replaceNaN(logGamma)

        #check convergence
        logprob = self.__logSumExp(logAlpha[:, -1])
        if abs(logprob-lastlogprob)<=10**-6:
          if verbose:
            sys.stderr.write("\n")
          break
        lastlogprob = logprob

        #M-step:
        self.logInitialProb[:] = logGamma[0, :]
        for i in range(self.N):
          for j in range(self.N):
            numerator = self.__logSumExp(logXi[:,i,j].flatten())
            denominator = self.__logSumExp(logGamma[:,i])
            self.logTransProb[i,j] = numerator - denominator
        self.logTransProb = self.__replaceNaN(self.logTransProb)

        for i in range(self.N):
          for j in range(self.K):
            emit = np.where(data[j,:]==1)
            if len(emit[0]) == 0:
              numerator = -np.inf
            else:
              numerator = self.__logSumExp(logGamma[emit[0],i])
            denominator = self.__logSumExp(logGamma[:,i])
            #print numerator, denominator
            self.logEmitProb[i,j] = self.__replaceNaN(numerator - denominator)

        #manually set special emission states
        #self.logEmitProb[self.StateDelete,:] = \
        #    np.array([-np.inf,-np.inf,-np.inf,-np.inf,0,-np.inf])
        #self.logEmitProb[[0,self.N-1],:] = \
        #    np.array([-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,0])

      if verbose:
        sys.stderr.write("\n")

    #got a list -- learn on multiple sequences
    else:
      numseqs = len(data)
      lastaverageLogprob = -np.inf
      logAlphas = [None] * numseqs
      logBetas = [None] * numseqs
      logXis = [None] * numseqs
      logGammas = [None] * numseqs
      logprobs = [None] * numseqs
      dataarray = np.concatenate(data, 1)
      for iteration in xrange(numsteps):
        if verbose:
          sys.stderr.write("\rIteration: %d/%d"%(iteration+1, numsteps))
          sys.stderr.flush()
        #E-step:
        for seqindex, d in enumerate(data):
          numclasses, numpoints = d.shape
          logAlphas[seqindex], logBetas[seqindex] = self.__logForwardBackward(d)
          #compute xi
          assert numclasses == self.K
          logXis[seqindex] = np.zeros(
              (numpoints-1,self.N,self.N), dtype=float)
          for t in range(numpoints-1):
            for i in range(self.N):
              for j in range(self.N):
                logXis[seqindex][t, i, j] = \
                    logAlphas[seqindex][i, t]+\
                    self.logTransProb[i, j]+\
                    self.logEmitProb[j,np.where(d[:,t+1]==1)[0][0]]+\
                    logBetas[seqindex][j, t+1]
            logXis[seqindex][t, :, :] -= self.__logSumExp(
                logXis[seqindex][t, :, :].flatten())
          logXis[seqindex] = self.__replaceNaN(logXis[seqindex])

          #compute gamma
          logGammas[seqindex] = np.zeros((numpoints, self.N), dtype = float)
          for t in range(numpoints):
            for i in range(self.N):
              logGammas[seqindex][t, i] = logAlphas[seqindex][i,t] +\
                  logBetas[seqindex][i,t]
            logGammas[seqindex][t,:] -=\
                self.__logSumExp(logGammas[seqindex][t,:].flatten())
          logGammas[seqindex] = self.__replaceNaN(logGammas[seqindex])

          logprobs[seqindex] = self.__logSumExp(logAlphas[seqindex][:, -1])

        #check convergence
        averageLogprob = np.mean(logprobs)
        if abs(averageLogprob-lastaverageLogprob)<=10**-6:
          if verbose:
            sys.stderr.write("\n")
          break
        lastaverageLogprob = averageLogprob

        #M-step:
        self.logInitialProb = logGammas[0][0,:]
        logXisArray = np.concatenate(logXis, 0)
        logGammasArray = \
            np.concatenate(map(lambda x: x[:-1,:],logGammas), 0)
        self.logTransProb = self.__replaceNaN(self.__logSumExp(logXisArray, 0) \
            - self.__logSumExp(logGammasArray, 0)[:, np.newaxis])
        logGammasArray = np.concatenate(logGammas, 0)
        G = np.exp(self.__replaceNaN(logGammasArray-\
            self.__logSumExp(logGammasArray,0)[np.newaxis,:]))

        #manually set special emission states
        for k in range(self.N):
          if k in self.StateDelete:
            self.logEmitProb[k,:] = \
              np.array([-np.inf,-np.inf,-np.inf,-np.inf,0,-np.inf])
          elif k == 0 or k == self.N - 1:
            self.logEmitProb[k,:] = \
              np.array([-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,0])
          else:
            self.logEmitProb[k,:] = \
                self.__replaceNaN(np.log(np.sum(G[:, k][np.newaxis,:]*dataarray, 1))-\
                np.log(np.sum(G[:,k])))

      if verbose:
        sys.stderr.write("\n")

  def ViterbiDecoding(self, obs):
    """The input obs includes B and E as an emitting state, which emits a $.
    This means this method can only deal with the case that all states are emitting.
    """
    #initialProb = self.transProbDecoding[0,:,None]
    initialProb = np.zeros((self.N,1))
    initialProb[0] = 1
    trellis = np.log(np.zeros((self.transProbDecoding.shape[0], len(obs))))
    backpt = np.ones((self.transProbDecoding.shape[0], len(obs)), 'int32') * -1
    logInitialProb = np.log(initialProb)
    logTransProb = np.log(self.transProbDecoding)
    logEmitProb = self.logEmitProb
    
    # initialization
    trellis[:, 0] = np.squeeze(logInitialProb + self.__logEmit(obs[0]))
    
    # decoding
    for t in xrange(1, len(obs)):
      trellis[:, t] = (self.__logDot(trellis[:, t-1, None],self.__logEmit(obs[t]).T) + \
          logTransProb).max(0)
      backpt[:, t] = (np.tile(trellis[:, t-1, None], [1, self.N]) + \
          logTransProb).argmax(0)

    # back track
    tokens = [trellis[:, -1].argmax()]
    for i in xrange(len(obs)-1, 0, -1):
      tokens.append(backpt[tokens[-1], i])
  
    return tokens[::-1]

  def ViterbiDecodingNew(self, obs, no_del = False):
    """Viterbi decoding method which can deal with silent states. The input obs
    must conform with this: no -, no $.
    """
    trellis = np.log(np.zeros((self.transProb.shape[0], len(obs)+1)))
    backpt = np.ones((self.transProb.shape[0], len(obs)+1), 'int32') * -1
    if no_del:
      EmitStates = self.NonDeletionEmitStates
      SilentStates = self.WithDeletionSilentStates
    else:
      EmitStates = self.EmitStates
      SilentStates = self.SilentStates
    
    # initialization
    trellis[0, 0] = 0
    
    # decoding
    for t in xrange(1, len(obs)+1):
      for i in EmitStates:
        trellis[i, t] = (trellis[:, t-1] + self.logTransProb[:,i]).max() + \
            self.logEmitProb[i,obs[t-1]]
      for ind, state in enumerate(SilentStates):
        pool = EmitStates + SilentStates[:ind]
        trellis[state, t] = (trellis[pool, t] + \
            self.logTransProb[pool, state]).max()

      for i in EmitStates:
        backpt[i, t] = (trellis[:, t-1] + self.logTransProb[:, i]).argmax()
      for ind, state in enumerate(SilentStates):
        pool = EmitStates + SilentStates[:ind]
        temp = trellis[:, t] + self.logTransProb[:, state]
        for j in xrange(self.N):
          if j not in pool:
            temp[j] = -np.inf
        backpt[state, t] = temp.argmax()

    tokens = [self.N-1, (trellis[:, -1] + self.logTransProb[:, -1]).argmax()]
    i = len(obs)
    while i > 0:
      tokens.append(backpt[tokens[-1], i])
      if tokens[-1] == 0:
        break
      if tokens[-2] in EmitStates:
        i -= 1

    if tokens[-1] > 0:
      tokens.append(0)
    return tokens[::-1]

  def LogLikelihood(self, states, seqs):
    """Compute log likelihood from given state path and emission sequence
    based on current model.
    """
    # state path includes B and E, so emission sequence includes $
    assert len(states) == len(seqs)
    log_lh = 0
    for i in range(len(states)-1):
      log_lh += self.logTransProb[states[i]][states[i+1]]
      log_lh += self.logEmitProb[states[i]][seqs[i]]

    return log_lh

  def OutputTransProb(self, outf = sys.stdout):
    """Print transistion matrix to file/stdout.
    """
    if not outf == sys.stdout:
      outf = open(outf,'w')

    data = np.exp(self.logTransProb)
    for i in range(self.N):
      outline = "\t".join([str(j) for j in data[i,:]]) + "\n"
      outf.write(outline)

    if not outf == sys.stdout:
      outf.close()

  def OutputEmitProb(self, outf = sys.stdout):
    """Print emission matrix to file/stdout.
    """
    if not outf == sys.stdout:
      outf = open(outf,'w')

    data = np.exp(self.logEmitProb)
    for i in range(self.N):
      outline = "\t".join([str(j) for j in data[i,:]]) + "\n"
      outf.write(outline)

    if not outf == sys.stdout:
      outf.close()
      
  def GetModelLength(self):
    length = 0
    for i in range(self.emitProb.shape[0]):
      if self.emitProb[i, :4].sum() > 0:
        length += 1

    return (length - 1)/2

  def ForwardLogLikelihood(self, seq):
    """Compute the total likelihood of given sequence using forward
    algorithm. Input sequence is in alphabet index form.
    """
    # N: number of states; K: size of alphabet; L: length of obeservation sequence
    numpoints = len(seq)
    logAlpha = np.zeros((self.N, numpoints), dtype=float)

    for k in xrange(self.N):
      #logAlpha[k, 1] = self.logInitialProb[k] +\
      #    self.logEmitProb[k, seq[0]]
      logAlpha[k, 0] = -np.inf
    logAlpha[0, 0] = 0.0

    for t in range(numpoints-1):
      for k in range(self.N):
        logAlpha[k, t+1] = self.__logSumExp(self.logTransProb[:, k]+\
            logAlpha[:, t])+\
            self.logEmitProb[k,seq[t+1]]

    log_lh = self.__logSumExp(logAlpha[:,numpoints-1])

    return log_lh

  def ForwardLogLikelihoodNew(self, seq, no_del = False):
    """Compute the total likelihood of given sequence using forward
    algorithm. Input sequence is in alphabet index form.
    If no_del is True, input must have no - nor $.
    """
    # N: number of states; K: size of alphabet; L: length of obeservation sequence
    numpoints = len(seq)
    if no_del:
      EmitStates = self.NonDeletionEmitStates
      SilentStates = self.WithDeletionSilentStates
    else:
      EmitStates = self.EmitStates
      SilentStates = self.SilentStates
    logAlpha = np.zeros((self.N, numpoints+1), dtype=float)
    if not np.isinf(self.logEmitProb[1,4]):
      # model is local alignment
      isLocal = True
    else:
      isLocal = False

    for k in xrange(self.N):
      logAlpha[k, 0] = -np.inf
    if isLocal:
      logAlpha[:, 0] = self.logTransProb[0,:]
    else:
      logAlpha[0, 0] = 0.0

    ## if there are more than 1 leading silent states: recalculate column 1
    #if EmitStates[0] > SilentStates[1]:
    #  temp = np.log(np.zeros((self.N, 1), dtype = float))
    #  for ind, state in enumerate(SilentStates):
    #    temp[state] = self.__logSumExp(self.logTransProb[:, state]+\
    #        logAlpha[:, 0])
    #  logAlpha[:,0,None] = temp

    for t in range(numpoints):
      for k in EmitStates:
        logAlpha[k, t+1] = self.__logSumExp(self.logTransProb[:, k]+\
            logAlpha[:, t])+\
            self.logEmitProb[k,seq[t]]
      for ind, state in enumerate(SilentStates):
        pool = EmitStates + SilentStates[:ind]
        logAlpha[state, t+1] = self.__logSumExp(self.logTransProb[pool, state]+\
            logAlpha[pool, t+1])

    log_lh = logAlpha[-1,-1]
    #if True:
    #  for i in range(self.N):
    #    print "%d\t"%i, logAlpha[i,:]
    #  print "#############"
    #  for i in range(self.N):
    #    print "%d\t"%i, self.logTransProb[i,:]

    return log_lh

  def SimpleLinearModelToLocalAlignment(self):
    """Add flanking deletion states to enable local alignment modeling.
    This local alignment can be set to only allow flanking gaps, or also
    allow partial alignment of profile HMM by enabling transition to
    all matching states.
    """
    # probability of starting a flanking gap
    GapStart = 0.6
    # probability of extening the flanking gap
    GapExt = 0.9
    # probability of alignment starts from the beginning of the model
    AlignStart = 0.5
    # probability of alignment ends at the end of the model
    AlignEnd = AlignStart

    # changing emission matrix
    self.emitProb = np.concatenate((np.zeros((3, self.K)), self.emitProb))
    self.emitProb = np.concatenate((self.emitProb, np.zeros((3, self.K))))
    self.emitProb[1, 4] = 1.0
    self.emitProb[-2, 4] = 1.0
    self.logEmitProb = np.log(self.emitProb)

    # add 3 states before the beginning and after ending respectively
    augment = np.concatenate((np.zeros((3, self.N)), np.eye(self.N)))
    augment = np.concatenate((augment, np.zeros((3,self.N))))
    self.transProb = np.dot(np.dot(augment, self.transProb), augment.T)

    # State 0, the very beginning. Set gap start probability
    self.transProb[0, 1] = GapStart
    self.transProb[0, 2] = 1 - GapStart
    # State 1, the looping state for flanking gap
    self.transProb[1, 1] = GapExt
    self.transProb[1, 2] = 1 - GapExt
    # State 2, the silent state for connecting gap and match states of HMM.
    # there are in total Length + 1 connections from this state to B, M_i
    # note this state does not connect E, because I don't want 0 match
    self.transProb[2, 3] = AlignStart
    for i in range(self.linearLength):
      self.transProb[2, i*3 + 5] = (1 - AlignStart) / (self.linearLength + 0)
    #self.transProb[2, self.linearLength*3 + 5] = (1 - AlignStart) /\
    #    (self.linearLength + 1)
    # State -3, the silent state for connecting match states and right gap.
    self.transProb[3, :] = self.__RescaleProbabilities(self.transProb[3, :],\
        (1-AlignEnd)/(self.linearLength+0) )
    self.transProb[3, -3] = (1 - AlignEnd) / (self.linearLength + 1)
    for i in range(self.linearLength):
      self.transProb[i*3 + 5, :] = self.__RescaleProbabilities(\
          self.transProb[i*3 + 5, :],\
          (1-AlignEnd)/(self.linearLength+0) )
      self.transProb[i*3 + 5, -3] = (1 - AlignEnd) / (self.linearLength + 0)
    self.transProb[self.linearLength*3 + 5, -3] = 1
    # State -2, the right flanking gap state
    self.transProb[-3, -2] = GapStart
    self.transProb[-3, -1] = 1 - GapStart
    self.transProb[-2, -2] = GapExt
    self.transProb[-2, -1] = 1 - GapExt

    self.logTransProb = np.log(self.transProb)
    self.N = self.transProb.shape[0]
    self.__initiateStates()

  def _NonEmittingStateReductionSimpleCyclic(self):
    """Modify transition probability matrix to remove non-emitting states
    because they will not be needed in decoding.
    """
    transProb = self.transProb
    # there are two states to be removed: E state of the cyclic part and 
    # B state of the linear part.
    E_idx = self.cyclicLength * 3 + 2
    B_idx = self.cyclicLength * 3 + 3

    # change transition probability of last M, I, D of cyclic part
    addAmount = self.transProb[E_idx,0]*self.transProb[0,:] + \
        self.transProb[E_idx,B_idx]*self.transProb[B_idx,:]
    for i in range(3):
      stateIndex = E_idx - i - 1
      transProb[stateIndex,:] = transProb[stateIndex,:] +\
          addAmount * transProb[stateIndex,E_idx]
      transProb[stateIndex,E_idx] = 0

    return transProb

  def __choice(self, prob):
    cutoff = np.cumsum(prob)
    index = cutoff.searchsorted(np.random.uniform())
    return index

  def __replaceNaN(self, x, rep = -np.inf):
    """Replace NaN in x by rep. This is needed in logarithm
    computations.
    """
    if not type(x) == type(np.array([])):
      if np.isnan(x):
        x = rep
    else:
      x[np.isnan(x)] = rep

    return x

  def __logSumExp(self, x, dim=-1):
    """Compute log(sum(exp(x))) in a numerically stable way.
    Use second argument to specify along which dimensions the logsumexp
    shall be computed. If -1 (which is the default), logsumexp is 
    computed along the last dimension. 
    """
    if len(x.shape) < 2:
      #xmax = np.nanmax(x)
      #return xmax + np.log(np.sum(np.exp(x-xmax)))
      if np.isnan(x.max()) or np.isinf(x.max()):
        return -np.inf
      else:
        xmax = x.max()
        return xmax + np.log(np.sum(np.exp(x-xmax)))
    else:
      if dim != -1:
        x = x.transpose(range(dim) + range(dim+1, len(x.shape)) + [dim])
      lastdim = len(x.shape)-1
      #xmax = x.max(lastdim)
      xmax = np.nanmax(x,lastdim)
      return xmax + np.log(np.sum(
        np.exp(self.__replaceNaN(x-xmax[...,np.newaxis])),lastdim))

  def __logForwardBackward(self, data):
    """Compute log-alpha and log-beta tables with Forward-backward
    algorithm.
    Data is expected to be in one-hot encoding. 
    """
    # N: number of states; K: size of alphabet; L: length of obeservation sequence
    numclasses, numpoints = data.shape
    logAlpha = np.zeros((self.N, numpoints), dtype=float)
    logBeta = np.zeros((self.N, numpoints), dtype=float)

    for k in xrange(self.N):
      #logAlpha[k, 1] = self.logInitialProb[k] +\
      #    self.logEmitProb[k,np.where(data[:,0]==1)[0][0]]
      logAlpha[k, 0] = -np.inf
      logBeta[k, -1] = 0.0
    logAlpha[0, 0] = 0.0

    for t in range(numpoints-1):
      for k in range(self.N):
        logAlpha[k, t+1] = self.__logSumExp(self.logTransProb[:, k]+\
            logAlpha[:, t])+\
            self.logEmitProb[k,np.where(data[:,t+1]==1)[0][0]]
            #self.__replaceNaN(self.logEmitProb[k,:]*data[:,t+1],0).sum()
    for t in range(numpoints-2, -1, -1):
      for k in range(self.N):
        logBeta[k, t] = self.__logSumExp(logBeta[:,t+1]+\
            self.logTransProb[k, :]+\
            #self.__replaceNaN(self.logEmitProb[k,:]*data[:,t+1],0).sum())
            self.logEmitProb[:,np.where(data[:,t+1]==1)[0][0]])

    return logAlpha, logBeta

  def __logEmit(self, obs):
    """Return the probability of emitting given observation.
    Dimension: (N,1).
    """
    return self.logEmitProb[:, obs, None]

  def __logDot(self, a, b):
    """Compute dot product of a and b with logarithm values.
    Will use __logSumExp function for summation.
    """
    assert a.shape[1] == b.shape[0]
    sum = np.zeros((a.shape[0],b.shape[1]))

    for i in xrange(a.shape[0]):
      for j in xrange(b.shape[1]):
        sum[i,j] = a[i,:] + b[:,j]

    return sum

  def __initiateStates(self):
    """Initiate a list for state indices in transition probability
    matrix by state types.
    """
    #if self.type == 1:
    #  # linear model
    #  self.StateBegin = [0]
    #  self.StateEnd = [len(self.transProb) - 1]
    #  self.StateMatch = [i*3+2 for i in range(self.linearLength)]
    #  self.StateInsert = [1] + [i*3+3 for i in range(self.linearLength)]
    #  self.StateDelete = [i*3+4 for i in range(self.linearLength)]
    #elif self.type == 2:
    #  # 1 cyclic model + 1 linear model:
    #  self.StateBegin = [0, self.cyclicLength*3+3]
    #  self.StateEnd = [self.cyclicLength*3+2, len(self.transProb) - 1]
    #  self.StateMatch = [i*3+2 for i in range(self.cyclicLength)] +\
    #      [i*3+self.cyclicLength*3+5 for i in range(self.linearLength)]
    #  self.StateInsert = [1] + [i*3+3 for i in range(self.cyclicLength)] +\
    #      [self.cyclicLength*3+4] + \
    #      [i*3+self.cyclicLength*3+6 for i in range(self.linearLength)]
    #  self.StateDelete = [i*3+4 for i in range(self.cyclicLength)] +\
    #      [i*3+self.cyclicLength*3+7 for i in range(self.linearLength)]

    self.SilentStates = []
    self.EmitStates = []
    self.NonDeletionEmitStates = []
    self.WithDeletionSilentStates = []
    for i in range(self.N):
      if self.emitProb[i,:4].sum() > 1e-8:
        self.NonDeletionEmitStates.append(i)
      else:
        self.WithDeletionSilentStates.append(i)
      if self.emitProb[i,:5].sum() > 1e-8:
        self.EmitStates.append(i)
      else:
        self.SilentStates.append(i)

      #print i
      #print 'wd:', abs(self.emitProb[i,4:].sum())
      #if abs(self.emitProb[i,4:].sum()) > 1e-8:
      #  print '!wd'
      #  self.WithDeletionSilentStates.append(i)
      #print 's:', abs(self.emitProb[i,:].sum()-1) ,abs(self.emitProb[i,5:]) 
      #if abs(self.emitProb[i,:].sum()-1) > 1e-8 or \
      #    abs(self.emitProb[i,5:]) > 1e-8:
      #  print '!s'
      #  self.SilentStates.append(i)
      #print 'e:', abs(self.emitProb[i,0:5].sum()) 
      #if abs(self.emitProb[i,0:5].sum()) > 1e-8:
      #  print '!e'
      #  self.EmitStates.append(i)
      #  print 'nd:', abs(self.emitProb[i,4]) 
      #  if abs(self.emitProb[i,4]) < 1e-8 :
      #    print '!nd'
      #    self.NonDeletionEmitStates.append(i)

  def __loadTransProbFromFile(self, file):
    """Load transition probability matrix from file handler.
    Values must NOT be in logarithm form.
    """
    fhd = open(file,'r')
    row = 0
    l = fhd.readline()
    while l:
      if l.startswith("#"):
        pass
      else:
        f = l.split()
        self.N = len(f)
        break
    fhd.seek(0)
    self.transProb = np.zeros((self.N,self.N))
    for l in fhd:
      #skip comment lines
      if l.startswith("#"):
        continue
      f = l.split()
      assert len(f) == self.N
      if row >= self.N:
        sys.stderr.write("Transition matrix is not a sqaure matrix!\n")
        sys.exit(1)
      for col in range(self.N):
        self.transProb[row,col] = float(f[col])
      row += 1

    fhd.close()

  def __loadEmitProbFromFile(self, file):
    """Load emission probability matrix from file handler.
    Values must NOT be in logarithm form.
    """
    fhd = open(file,'r')
    self.emitProb = np.zeros((self.N, self.K))
    row = 0
    for l in fhd:
      if l.startswith("#"):
        continue
      else:
        f = l.split()
        assert len(f) == self.K
        if row >= self.N:
          sys.stderr.write("Emission matrix has wrong dimension!\n")
          sys.exit(1)
        for col in range(self.K):
          self.emitProb[row, col] = float(f[col])
        row += 1

    fhd.close()

  def __loadALN(self, file):
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
    self.transProb = np.zeros((N,N))
    self.emitProb = np.zeros((N,self.K))
    self.emitProb[0,-1] = 1
    self.emitProb[-1,-1] = 1

    # mark consensus states: 0=M, 1=I, 2=D
    marked_states = []
    for i in range(len(data['consensus'])):
      if data['consensus'][i] == '-':
        marked_states.append(1)
      else:
        marked_states.append(0)

    # get counts of states and emits in each position
    states = {}.fromkeys(data.keys())
    emits = {}.fromkeys(data.keys())
    for name in data.keys():
      if name == 'consensus':
        continue
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
        emits[name].append(seq[i])

    #for name in data.keys():
    #  if not name == 'consensus':
    #    print name, self.TranslateStatesSimpleLinearModel(\
    #      [self.__StateIndToTransProbInd(states[name], i) for i in range(len(states[name]))]\
    #      )
    #    print data[name]
    #    print data['consensus']
    #sys.exit(0)

    # pseudo count
    PSEUDO = min(int(0.05*len(data.keys())), 20)
    # B, I0
    self.transProb[0,[1,2,4]] += PSEUDO
    self.transProb[1,[1,2,4]] += PSEUDO
    # M_i, I_i, D_i
    for i in range(length-1):
      self.transProb[i*3+2, [i*3+3,i*3+5,i*3+7]] += PSEUDO
      self.transProb[i*3+3, [i*3+3,i*3+5,i*3+7]] += PSEUDO
      self.transProb[i*3+4, [i*3+3,i*3+5,i*3+7]] += PSEUDO
    # M_L to I_L, E
    self.transProb[(length-1)*3+2,\
        [(length-1)*3+3,(length-1)*3+5]] += PSEUDO
    # I_L to I_L, E
    self.transProb[(length-1)*3+3,\
        [(length-1)*3+3,(length-1)*3+5]] += PSEUDO
    # D_L to I_L, E
    self.transProb[(length-1)*3+4,\
        [(length-1)*3+3,(length-1)*3+5]] += PSEUDO
    # Emission matrix
    # M, I
    self.emitProb[[2+i*3 for i in range(length)], :-2] += PSEUDO
    self.emitProb[[1]+[3+i*3 for i in range(length)], :-2] += PSEUDO
    # D
    self.emitProb[[4+i*3 for i in range(length)], self.K-2] += PSEUDO

    # calculate transition and emission counts
    for name in states.keys():
      if name == 'consensus':
        continue
      prev_ind = 0
      for i in range(len(states[name])):
        ind = self.__StateIndToTransProbInd(states[name], i)
        self.emitProb[ind][alphabet[emits[name][i]]] += 1
        self.transProb[prev_ind][ind] += 1
        prev_ind = ind
      self.transProb[prev_ind][-1] += 1

    # convert counts to probabilities
    self.transProb = self.__replaceNaN(self.transProb /\
        self.transProb.sum(1)[:,None],0)
    self.emitProb = self.__replaceNaN(self.emitProb /\
        self.emitProb.sum(1)[:,None],0)
    self.linearLength = self.GetModelLength()

  def __StateIndToTransProbInd(self, state, ind):
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

  def __RescaleProbabilities(self, entry, delta):
    """Rescale probabilites of given entry to make room for delta,
    maiking sure the sum of the entry is still zero.
    """
    for ind, val in enumerate(entry):
      entry[ind] = val*(1-delta)

    return entry


def FastaToSequence(file):
  """Load a fasta file and return a list of tuples:
  (sequence name, sequence).
  """
  fhd = open(file, 'r')
  data = {}
  for l in fhd:
    if l.startswith(">"):
      name = l.strip()[1:]
      if not data.has_key(name):
        data[name] = []
    else:
      data[name].append(l.strip())

  fhd.close()
  seqs = [(name, "$"+"".join(data[name])+"$") for name in data.keys()]

  return seqs

def NucleotideToOneHot(seq):
  """Covert a sequence of nucleotide (ACGT-) to one-hot encoded
  numpy array.
  N(row) = alphabet size; N(column) = sequence length
  """
  alphabet = {"A":0, "C":1, "G":2, "T":3, "-":4, "$":5}
  data = np.zeros((len(alphabet.keys()), len(seq)))

  for i in xrange(len(seq)):
    data[alphabet[seq[i]], i] = 1

  return data

def main():
  #hmm = ProfileHMM(1, [10])
  #hmm = ProfileHMM(2, [5,7])
  #hmm = ProfileHMM(2, [1, 2])
  #print hmm.transProb
  #print hmm.emitProb
  #print hmm.transProb.shape
  #print hmm.emitProb.shape
  #print emitSeq
  #print hmm.TranslateEmitToDNA(emitSeq[0])
  #print hmm.TranslateStatesSimpleLinearModel(emitSeq[1])
  #print hmm.TranslateStatesSimpleCyclicModel(emitSeq[1])

  # Viterbi decoding
  #hmm = ProfileHMM(1, [3], int(sys.argv[1]))
  #emitSeq = hmm.EmitSequence()
  #print hmm.TranslateEmitToDNA(emitSeq[0])
  #print hmm.TranslateStatesSimpleLinearModel(emitSeq[1])
  ##print hmm.TranslateStatesSimpleLinearModel(hmm.ViterbiDecoding(emitSeq[0]))
  #print hmm.ForwardLogLikelihood(emitSeq[0])
  #emitBases = []
  #for i in emitSeq[0]:
  #  if i < 4:
  #    emitBases.append(i)
  ##print hmm.TranslateStatesSimpleLinearModel(hmm.ViterbiDecodingNew(emitBases, True))
  #print hmm.ForwardLogLikelihoodNew(emitSeq[0][1:-1])
  ##print hmm.ForwardLogLikelihoodNew(emitSeq[0])
  #print hmm.ForwardLogLikelihoodNew(emitBases, True)

  # Baum-Welch training
  # one seq
  #emitSeq = hmm.EmitSequence()
  #hmmb = ProfileHMM(1,[2])
  #print hmm.TranslateEmitToDNA(emitSeq[0])
  #print hmm.TranslateStatesSimpleLinearModel(emitSeq[1])
  #print hmm.TranslateStatesSimpleLinearModel(hmm.ViterbiDecoding(emitSeq[0]))
  #hmmb.BaumWelchTraining(NucleotideToOneHot(hmm.TranslateEmitToDNA(emitSeq[0])))
  #print hmmb.TranslateStatesSimpleLinearModel(hmmb.ViterbiDecoding(emitSeq[0]))
  #print hmmb.logTransProb
  # groups
  #obs_list = [hmm.EmitSequence() for i in range(10)]
  #hmmb = ProfileHMM(2,[2,3])
  ##hmmb = ProfileHMM(3,["trans.txt", "emit.txt"])
  #hmmb.BaumWelchTraining([NucleotideToOneHot(hmm.TranslateEmitToDNA(i[0])) \
  #    for i in obs_list], 5, True)
  #for i in obs_list:
  #  print hmm.TranslateEmitToDNA(i[0])
  #  #print hmm.TranslateStatesSimpleLinearModel(i[1])
  #  #print hmmb.TranslateStatesSimpleLinearModel(hmmb.ViterbiDecoding(i[0]))
  #  print hmm.TranslateStatesSimpleCyclicModel(i[1])
  #  print hmmb.TranslateStatesSimpleCyclicModel(hmmb.ViterbiDecoding(i[0]))
  #  print ""

  # local alignment
  hmm_short = ProfileHMM(1, [3])
  #test = [4,4,4,4,3,4,2,3,1,4,4,4]
  test = [3,4,2,3,1]
  test_base = [3]
  hmm_a = ProfileHMM(1, [2], 1)
  hmm_b = ProfileHMM(1, [5], 1)
  hmm_b.SimpleLinearModelToLocalAlignment()

  emit = hmm_short.EmitSequence()
  emitBases = []
  for i in emit[0]:
    if i < 4:
      emitBases.append(i)
  #print hmm_short.TranslateEmitToDNA(emit[0])
  #print hmm_short.TranslateStatesSimpleLinearModel(emit[1])
  #print hmm_a.TranslateStatesSimpleLinearModel(hmm_a.ViterbiDecoding(emit[0]))
  #print hmm_a.TranslateStatesSimpleLinearModel(hmm_a.ViterbiDecodingNew(emit[0][1:-1]))
  #print hmm_a.ForwardLogLikelihoodNew(emit[0][1:-1])
  print hmm_b.ForwardLogLikelihoodNew(test)
  hmm_b.SetNullModel([0.7,0.1,0.1,0.1])
  print hmm_b.ForwardLogLikelihoodNew(test)
  #print hmm_b.ForwardLogLikelihoodNew(test_base, True)
  print hmm_a.ForwardLogLikelihood(test)
  #print hmm_b.TranslateStatesLocalAlignment(hmm_b.ViterbiDecodingNew(emit[0][1:-1]))
  #print hmm_b.TranslateStatesLocalAlignment(hmm_b.ViterbiDecodingNew([4,4,1,3,4,4,4,4]))


  #hmmb.OutputTransProb()
  #hmmb.OutputEmitProb()
  #print np.exp(hmmb.logInitialProb)
  #print np.exp(hmmb.logTransProb)
  #print np.exp(hmmb.logEmitProb)

  #hmm = ProfileHMM(4, sys.argv[1])
  #hmm.OutputTransProb()
  #hmm.OutputEmitProb()

  #hmm = ProfileHMM(1, [10])
  #for i in range(100):
  #  emit = hmm.EmitSequence()
  #  print hmm.TranslateEmitToDNA(emit[0])
  

if __name__ == '__main__':
  np.seterr("ignore")
  main()
