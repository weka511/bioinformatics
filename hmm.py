#!/usr/bin/env python

#    Copyright (C) 2019-2023 Simon Crase
#
#    This is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This software is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>

'''
   Hidden Markov Models
'''
from unittest import TestCase, main, skip

import numpy as np
from numpy.testing import assert_array_almost_equal

class StateSet:
    '''
    Manage states: S, I0, M1, D1, I1, M2, D2, I2, ...E
    '''

    def __init__(self,m):
        '''
        Initialize state set

        Parameters:
            m       Number of states
        '''
        self.m = m

    def get_pair(self,index):
        '''
        Convert index number to State, expresses as pair, e.g. 2->('M',1)
        '''
        if index   == 0:        return ('S',None)
        if index   == self.m-1: return ('E',None)
        if index%3 == 0:        return ('D', index//3)
        if index%3 == 1:        return ('I', index//3)
        if index%3 == 2:        return ('M', index//3+1)


    def get_index(self,state):
        '''
        Convert State to index number, e.g. ('M',1)->2
        '''
        function,i = state
        if function == 'S': return 0
        if function == 'E': return self.m-1
        if function == 'D': return 3*i
        if function == 'I': return 3*i+1
        if function == 'M': return 3*i-1


    def get_key(self,op,seq):
        '''
        Format state for printing, e.g. ('M',1)->'M1'
        '''
        return f'{op}-{seq}'

    def split_key(self,key):
        '''
        Used to read a state name, e.g. 'M1'->('M',1)
        '''
        parts =key.split('-')
        if parts[0] in ['M','D','I']:
            return parts[0],int(parts[1])
        else:
            return parts[0],parts[1]

    @classmethod
    def create_path(cls,s,mask):
        '''
        Parse a substring into a sequence of Match and Delete states

        Parameters:
              s     String from which subset is to be extracted
              mask  Indicates which elements of are to be included in subset

        Returns:
           A set of States, starting with ('S',) and ending with ('E'). The other elements are matches or
           deletes, depending on whether corresponding character is a non-space character, each accomanied by the
           corresponding character
        '''
        product = [('S',None,None)]
        index   = 0
        assert len(s)==len(mask)
        for i in range(len(s)):
            match mask[i],is_space(s[i]):
                case (True,False):
                    index += 1
                    product.append(('M',index,s[i]))
                case (True,True):
                    index += 1
                    product.append(('D',index,s[i]))
                case (False,True):  # Nothing needs to happen skipping spaces
                    pass
                case (False,False):
                    product.append(('I',index,s[i]))

        product.append(('E',None,None))
        return product

    def get_successors(self,state):
        block = (state+1)//3
        for i in range(3):
            successor = 3*block+i+1
            if successor<self.m:
                yield(successor)

    @classmethod
    def format(cls,s):
        code= s
        if code in ['S','E']:
            return code
        else:
            return f'{code}'


def get_indices(S,Alphabet='AB'):
    '''
    Used to convert a string into integers representing the position of each character in the alphabet.
    It can also be used to convert a string of states into state indices.

    Parameters:
       S         String to be converrted
       Alphabet

    Returns:
       List of indices
    '''
    IndexTable = {}
    for i in range(len(Alphabet)):
        IndexTable[Alphabet[i]]=i
    return [IndexTable[s] for s in S]



def is_space(ch):
    '''
    Check whether character is a space (represented by hyphen)
    '''
    return ch=='-'

def normalize_rows(A):
    '''
    Normalize a matrix so each row sums to 1, unless all elements of row are 0
    '''
    m,_        = A.shape
    row_totals = A.sum(axis=1)

    for i in range(m):
        if row_totals[i]==0:
            row_totals[i]=1

    return A/row_totals.reshape(m,1)

def float2str(x,
              precision=3,
              p0=0,
              p1=1):
    '''
     float2str

     Format a floating point number in the manner that the
     grader wants for matrices when dealing with HMM

     Parameters:
         x         Value to be displayed
         precision Number of digits (after decimal point)
    '''
    format_str = f'{{:.{precision}f}}'
    if x==0:
        format_str = f'{{:.{p0}f}}'
    if x==1:
        format_str = f'{{:.{p1}f}}'
    format3= format_str.format(x)
    while len(format3)>3 and format3[-1]=='0':
        format3 = format3[:-1]
    return format3




def formatEmission(Emission,States,Alphabet,precision=2):
    '''
    formatEmission
    '''
    m,n = Emission.shape
    yield '\t'.join(Alphabet)
    for i in range(m):
        row = []
        for j in range(n):
            row.append(float2str(Emission[i,j],precision))
        yield StateSet.format(States[i]) + '\t' + '\t'.join(row)


def formatTransition(Transition,States,precision=2):
    '''
    formatTransition
    '''
    m,n = Transition.shape
    yield  '\t'.join(StateSet.format(state) for state in States)
    for i in range(m):
        row = []
        for j in range(n):
            row.append(float2str(Transition[i,j],precision))
        yield StateSet.format(States[i]) + '\t' + '\t'.join(row)


def ProbabilityHiddenPath(path,States,Transition):
    '''
    ProbabilityHiddenPath

    Solves: BA10A Compute the Probability of a Hidden Path

    Given: A hidden path followed by the states States and transition matrix Transition of an HMM.

    Return: The probability of this path. You may assume that initial probabilities are equal.
    '''
    def logP_Transition(i,path_indices):
        return logTransition[(path_indices[i-1],path_indices[i])]
    path_indices  = get_indices(path)
    logTransition = np.log(Transition)
    _,n           = logTransition.shape
    return np.exp(sum([logP_Transition(i,path_indices) for i in range(1,len(path_indices))]))/n

def ProbabilityOutcomeGivenHiddenPath(string,alphabet,path,states,Emission):
    '''
    ProbabilityOutcomeGivenHiddenPath

    Solves: BA10B Compute the Probability of an Outcome Given a Hidden Path

    Inputs: string
            alphabet          from which x was constructed
            path              hidden path
            states
            Emission

    A string x, followed by the alphabet from which x was constructed, followed by a hidden path,
    followed by the states  and emission matrix Emission of an HMM.

    Return: The conditional probability  that string x will be emitted by the HMM given the hidden path
    '''

    logEmission = np.log(Emission)
    return np.exp(sum([logEmission[a,b] for a,b in zip(get_indices(path,Alphabet=states),
                                                       get_indices(string,Alphabet=alphabet))]))


def Viterbi(x,Alphabet,States,Transition,Emission):
    '''
    Viterbi

    Solve:    BA10C Implement the Viterbi Algorithm

    Input:
        x           A string
        Alphabet    alphabet from which s was constructed,
        States      states of HMM
        Transition  transition matrix ,
        Emission    emission matrix

    Return:
        A path that maximizes the (unconditional) probability Pr(x, p) over all possible paths p
    '''

    def create_weights(s_source = 1):
        '''
        create_weights

        Calculate array of weights by position in string and state

        '''
        def get_max_weight(state,i):
            '''
            Calculate weight-- see textbook "The Viterbi algorithm

            Parameters:
                state     A state
                i         An index used to iterate through string x

            Returns:
                The maximum weight of all transitions from a possible state at position i-1 to i,
                where the weight is the product of the transion probability and emission probability
                of the character at position i.
            '''
            return np.multiply(s[i-1,:],Transition[:,state]* Emission[state,x_indices[i]]).max()

        m         = len(x)
        n         = len(States)
        x_indices = get_indices(x,Alphabet)
        s         = np.zeros((m,n))         # Product weights - see textbook "The Viterbi algorithm"
        s[0,:]    = s_source  * Emission[:,x_indices[0]] / n
        for i in range(1,len(x_indices)):
            s[i,:] = np.array([get_max_weight(state,i) for state in range(n)])

        return s

    def find_most_likely_path(s):
        '''
        find_most_likely_path

        Find most likely path through state space by backtracking.

        We start with last position, build backwards, and then reverse.

        '''
        m,_    = s.shape
        m     -= 1
        index_most_likely_state  = np.argmax(s[m])
        path   = [States[index_most_likely_state]]
        for i in range(m-1,-1,-1):
            index_most_likely_state = np.argmax([s[i,j] * Transition[j,index_most_likely_state]  for j in range(len(States))])
            path.append(States[index_most_likely_state])
        return path[::-1]

    return ''.join(find_most_likely_path(create_weights()))



def Likelihood(xs,Alphabet,States,Transition,Emission):
    '''
    Likelihood

      BA10D 	Compute the Probability of a String Emitted by an HMM

      Given: A string x, followed by the alphabet from which x was constructed, followed by the states,
            transition matrix, and emission matrix  of an HMM

      Return: The probability Pr(x) that the HMM emits x.
    '''
    def create_weights(xs, s_source = 1):
        def product_weight(k,x):
            return sum([s[-1][j] * Transition[j,k] * Emission[k,x] for j in range(n)])

        n = len(States)
        s = []  # first index is position, 2nd state

        s.append([s_source * (1/len(States)) * Emission[k,xs[0]] for k in range(n)])
        for x in xs[1:]:
            s.append([product_weight(k,x) for k in range(n)])
        return s

    return sum(create_weights(get_indices(xs,Alphabet=Alphabet))[-1])

def ConstructProfileHMM(theta,Alphabet,Alignment,sigma=0):
    '''
    ConstructProfileHMM

    BA10E Construct a Profile HMM and
    BA10F Construct a Profile HMM with Pseudocounts

    Parameters:
        theta      Threshold: exclude columns from an alignment if the fraction of spaces is greater than or equal to theta
        Alphabet   Symbols used for Alignment
        Alignment  A list of strings representing aligned proteins
        sigma      Minimum probability assigned to legal transitions

    Returns:
        Transition    transition matrix
        Emission      emission matrix
        States        A list of states, e.g. [
                      ('S', None), ('I', 0), ('M', 1), ('D', 1), ('I', 1), ..., ('E', None)], so row 0 of each
                      matrix refers to ('S', None), etc.

   '''

    def create_mask():
        '''
        Construct a mask to exclude columns from an alignment
        if the fraction of spaces is greater than or equal to theta
        '''
        def get_count(i):
            '''
            Returns:
                Number of spaces in column i
            '''
            return sum([is_space(c) for s in Alignment for c in s[i]])

        fractions = [get_count(i)/len(Alignment) for i in range(len(Alignment[0]))]
        return [fraction<theta for fraction in fractions]


    def create_transition():
        '''
        Create matrix of transition probabilities
        '''
        def create_census():
            States = {}
            for path in Paths:
                for (op1, seq1, _),(op2, seq2, _) in zip(path[:-1],path[1:]):
                    if state_set.get_key(op1,seq1)  not in States:
                        States[state_set.get_key(op1,seq1)] = []
                    States[state_set.get_key(op1,seq1)].append((op2, seq2))
            return States

        product = np.zeros((m,m))
        for i in range(m):
            for j in state_set.get_successors(i):
                product[i,j] = sigma

        States  = create_census()
        for key1,successors in States.items():
            state1,index1 = state_set.split_key(key1)
            counts      = {}
            for succ,seq in successors:
                if state_set.get_key(succ,seq) not in counts:
                    counts[state_set.get_key(succ,seq)] = 0
                counts[state_set.get_key(succ,seq)]   += 1
            fractions    = {key:count/len(successors) for key,count in counts.items()}
            state_index1 = state_set.get_index((state1,index1))
            for key2,fraction in fractions.items():
                state2,index2 = state_set.split_key(key2)
                state_index2 = state_set.get_index((state2,index2))
                product[state_index1,state_index2] = max(fraction,sigma)

        return normalize_rows(product)


    def create_emission():
        '''
        Create matrix of emission probabilities
        '''
        def create_census():
            States = {}
            for path in Paths:
                for op, seq, ch in path:
                    if state_set.get_key(op,seq)  not in States:
                        States[state_set.get_key(op,seq)] = []
                    States[state_set.get_key(op,seq)].append(ch)
            return States

        product = np.zeros((m,n))

        States  = create_census()
        for key,chars in States.items():
            state,index = state_set.split_key(key)
            if state in ['M','I']:
                counts = {ch:0 for ch in Alphabet}
                for ch in chars:
                    counts[ch] += 1
                fractions = {ch:max(count/len(chars),sigma) for ch,count in counts.items()}
                state_index = state_set.get_index((state,index))
                for j in range(n):
                    product[state_index,j] = fractions[Alphabet[j]]
        return normalize_rows(product)

    mask          = create_mask()
    Paths         = [StateSet.create_path(s,mask) for s in Alignment]
    m             = 3 * (sum([1 for m in mask if m]) + 1)
    state_set     = StateSet(m)
    n             = len(Alphabet)
    return create_transition(), create_emission(), [state_set.get_pair(i) for i in range(m)]

def AlignSequence(self,string,Alphabet,Emission):
    '''
    AlignSequenceWithProfileHMM

    BA10G   Sequence Alignment with Profile HMM Problem

    Given: A string Text, a multiple alignment Alignment, a threshold θ, and a pseudocount σ.

    Return: An optimal hidden path emitting Text in HMM(Alignment,θ,σ).
    '''
    raise Exception('not implemented')

def ViterbiLearning(s, Alphabet,  States, Transition, Emission,N=3):
    '''
        BA10I Implement Viterbi Learning

        Parameters:
        s           A String
        Alphabet    Characters from which string constructed
        States      States for HMM
        Transition  Transition probabilities
        Emission    Probabilities of symbols being emitted in each state

    '''

    for k in range(N):
        Path                   = Viterbi(s, Alphabet,  States, Transition, Emission)
        Transition1, Emission1 = EstimateParameters(s, Alphabet, Path, States)
        Transition = Transition1.copy()
        Emission   = Emission1.copy()
    return Transition, Emission

def SoftDecode(s, Alphabet,  States, Transition, Emission):
    '''
    BA10I Soft Decoding Problem

    Parameters:
        s                A String
        Alphabet         Characters from which string constructed
        States           States for HMM
        Transition       Transition probabilities
        Emission         Probabilities of symbols being emitted in each state

    Return:
        A Result object comprising the probability that the HMM was in state k at step i
        (for each state k and each step i) and the forward and backward matrices
    '''
    class Result:
        def __init__(self,probabilities,forward,backward):
            self.probabilities = probabilities/probabilities.sum(axis=1).reshape(-1,1)
            self.forward       = forward
            self.backward      = backward

    def get_forward():
        message       = np.full((m,n),np.nan)
        message[0][:] = Emission[:, x[0]]
        for i in range(1, m):
            for k in range(n):
                message[i][k] = np.dot(message[i-1][:],Transition[:, k])*Emission[k, x[i]]

        return message

    def get_backward():
        message        = np.full((m,n),np.nan)
        message[-1][:] = 1
        for i in range(m-2, -1, -1):
            for k in range(n):
                message[i][k] = (message[i+1,:]*Transition[k, :]*Emission[:, x[i+1]]).sum()

        return message

    m             = len(s)
    n             = len(States)
    x             = get_indices(s,Alphabet)
    forward       = get_forward()
    backward      = get_backward()
    probabilities = np.multiply(forward,backward)
    return Result(probabilities,forward,backward)

def get_edge_responsibilities(forward,backward, Transition, Emission, xindices):
    m,_ = Transition.shape
    edge_responsibilities = np.full((len(forward)-1,m,m),np.nan)
    for i in range(0,len(forward)-1):
        for l in range(m):
            for k in range(m):
                weight = Transition[l,k] * Emission[k,xindices[i+1]]
                edge_responsibilities[i,l,k] = forward[i,l] * weight * backward[i+1,k]

        edge_responsibilities[i,:,:] /= edge_responsibilities[i,:,:].sum()
    return edge_responsibilities

def BaumWelch(s, Alphabet,  States, Transition, Emission, N = 5):
    '''  BA10K 	Implement Baum-Welch Learning'''
    m,n = Emission.shape
    for k in range(N):
        decoded                           = SoftDecode(s,Alphabet, States,Transition,Emission)
        responsibilities,forward,backward = decoded.probabilities,decoded.forward,decoded.backward
        xindices                          = get_indices(s,Alphabet=Alphabet)
        edge_responsibilities             = get_edge_responsibilities(forward,backward, Transition, Emission, xindices)
        Transition1                       = edge_responsibilities.sum(axis=0)
        Transition1                      /= Transition1.sum(axis=1).reshape(-1,1)
        Transition                        = Transition1.copy()
        Emission1                         = np.zeros_like(Emission)
        for i in range(n):   # symbols
            for k in range(m): #states
                for j in range(len(xindices)):
                    if xindices[j] == i:
                        Emission1[k,i] += responsibilities[j,k]

        Emission1                      /= Emission1.sum(axis=1).reshape(-1,1)
        Emission                        = Emission1.copy()
    return Transition, Emission


def EstimateParameters(s,Alphabet,path,States):
    '''
    BA10H 	Estimate the Parameters of an HMM

    Parameters:
        s            A string
        Alphabet     Characters from which string has been built
        path         Path through HMM that generates string
        States       States used by HMM

    Returns:
        Estimated transition and emission probabilities

    '''
    def create_Transitions(path_indices,n):
        Transitions = np.zeros((len(States),len(States)))
        for i in range(1,n):
            Transitions[path_indices[i-1],path_indices[i]]+= 1
        if Transitions[-1,:].sum()==0:   # This is a hack, to prevent lines of all zeros
            Transitions[-1,:] = 1

        return normalize_rows(Transitions)

    def create_Emissions(path_indices,n):
        str_indices = get_indices(s,Alphabet)
        Emissions   = np.zeros((len(States),len(Alphabet)))
        for i in range(n):
            Emissions[path_indices[i],str_indices[i]]+= 1
        if Emissions[-1,:].sum()==0:                        # This is a hack, to prevent lines of all zeros
            Emissions[-1,:] = 1

        return normalize_rows(Emissions)

    n = len(path)
    assert n==len(s),f'String {s} and Path {path} should have the same length'
    path_indices = get_indices(path,States)
    return create_Transitions(path_indices,n),create_Emissions(path_indices,n)



if __name__=='__main__':
    class Test_10_HMM(TestCase):

        def test_ba10a(self):
            '''BA10A Compute the Probability of a Hidden Path'''

            self.assertAlmostEqual(5.01732865318,
                                        1e19*ProbabilityHiddenPath('AABBBAABABAAAABBBBAABBABABBBAABBAAAABABAABBABABBAB',
                                             'AB',
                                             np.array([[ 0.194, 0.806],
                                                       [ 0.273, 0.727]])),
                                        places=5)

        def test_ba10b(self):
            '''BA10B Compute the Probability of an Outcome Given a Hidden Path'''
            self.assertAlmostEqual(1.93157070893,
                                   1e28*ProbabilityOutcomeGivenHiddenPath('xxyzyxzzxzxyxyyzxxzzxxyyxxyxyzzxxyzyzxzxxyxyyzxxzx',
                                                                          'xyz',
                                                                          'BBBAAABABABBBBBBAAAAAABAAAABABABBBBBABAABABABABBBB',
                                                                          'AB',
                                                                          np.array([[0.612, 0.314, 0.074 ],
                                                                                    [0.346, 0.317, 0.336]])),
                                   places=5)

        def test_ba10c(self):
            '''BA10C Implement the Viterbi Algorithm'''
            Transition = np.array([
                [ 0.641,    0.359],
                [ 0.729,   0.271]])

            Emission = np.array([
                [0.117,  0.691, 0.192],
                [0.097,  0.42,  0.483]])

            self.assertEqual('AAABBAAAAA',
                             Viterbi('xyxzzxyxyy','xyz','AB',Transition,Emission))

        def test_ba10d(self):
            '''BA10D 	Compute the Probability of a String Emitted by an HMM'''
            self.assertAlmostEqual(1.1005510319694847,
                              1e6*Likelihood('xzyyzzyzyy',
                                         'xyz',
                                         'AB',
                                         np.array([[0.303, 0.697],
                                                   [0.831, 0.169]]),
                                         np.array([[0.533, 0.065, 0.402],
                                                   [0.342, 0.334, 0.324]])))

        def test_ba10e1(self):
            '''BA10E Construct a Profile HMM'''
            Transition, Emission,StateNames = ConstructProfileHMM(0.35,
                                                                  'ACDEF',
                                                                  ['ACDEFACADF',
                                                                   'AFDA---CCF',
                                                                   'A--EFD-FDC',
                                                                   'ACAEF--A-C',
                                                                   'ADDEFAAADF'])
            m1,m2 = Transition.shape
            m,n   = Emission.shape
            self.assertEqual(5,n)
            self.assertEqual(27,m)
            self.assertEqual(m1,m)
            self.assertEqual(m2,m)
            self.assertEqual(1, Transition[0,2])
            self.assertEqual(0.8, Transition[2,5])
            self.assertEqual(0.2, Transition[2,6])
            self.assertEqual(1, Transition[5,8])
            self.assertEqual(1, Transition[6,9])
            self.assertEqual(1, Transition[8,11])
            self.assertEqual(1, Transition[9,11])
            self.assertEqual(0.8, Transition[11,14])
            self.assertEqual(0.2, Transition[11,15])
            self.assertEqual(0.75, Transition[14,16])
            self.assertEqual(0.25, Transition[14,17])
            self.assertEqual(1, Transition[15,17])
            self.assertEqual(0.4, Transition[16,16])
            self.assertEqual(0.6, Transition[16,17])
            self.assertEqual(0.8, Transition[17,20])
            self.assertEqual(0.2, Transition[17,21])
            self.assertEqual(1, Transition[20,23])
            self.assertEqual(1, Transition[21,23])
            self.assertEqual(1, Transition[23,26])

        def test_ba10e2(self): #BA10E Construct a Profile HMM
            Transition, Emission,StateNames = ConstructProfileHMM(0.289,
                                                                  'ABCDE',
                                                                  ['EBA',
                                                                   'EBD',
                                                                   'EB-',
                                                                   'EED',
                                                                   'EBD',
                                                                   'EBE',
                                                                   'E-D',
                                                                   'EBD'])
            m1,m2 = Transition.shape
            m,n   = Emission.shape
            self.assertEqual(5,n)
            self.assertEqual(12,m)
            self.assertEqual(m1,m)
            self.assertEqual(m2,m)
            self.assertEqual(1, Transition[0,2])
            self.assertEqual(0.875, Transition[2,5])
            self.assertEqual(0.125, Transition[2,6])
            self.assertAlmostEqual(0.857, Transition[5,8], places=3)
            self.assertAlmostEqual(0.143, Transition[5,9],places=3)
            self.assertEqual(1, Transition[6,8])
            self.assertEqual(1, Transition[8,11])
            self.assertEqual(1, Transition[9,11])
            self.assertEqual(1, Emission[2,4])
            self.assertAlmostEqual(0.857, Emission[5,1], places=3)
            self.assertAlmostEqual(0.143, Emission[5,4], places=3)
            self.assertAlmostEqual(0.143, Emission[8,0], places=3)
            self.assertAlmostEqual( 0.714, Emission[8,3], places=3)
            self.assertAlmostEqual(0.143, Emission[8,4], places=3)

        def test_ba10e3(self): #BA10E Construct a Profile HMM
            Transition, Emission,StateNames = ConstructProfileHMM(0.252,
                                                                  'ABCDE',
                                                                  ['DCDABACED',
                                                                   'DCCA--CA-',
                                                                   'DCDAB-CA-',
                                                                   'BCDA---A-',
                                                                   'BC-ABE-AE'])
            m1,m2 = Transition.shape
            m,n   = Emission.shape
            self.assertEqual(5,n)
            self.assertEqual(18,m)
            self.assertEqual(m1,m)
            self.assertEqual(m2,m)

        def test_ba10f1(self):
            '''BA10F Construct a Profile HMM with Pseudocounts'''
            Transition, Emission,StateNames = ConstructProfileHMM(0.358,
                                                             'ABCD',
                                                             ['ADA',
                                                              'ADA',
                                                              'AAA',
                                                              'ADC',
                                                              '-DA',
                                                              'D-A'],
                                                             sigma=0.01)

            m1,m2 = Transition.shape
            m,n   = Emission.shape
            self.assertEqual(4,n)
            self.assertEqual(12,m)
            self.assertEqual(m1,m)
            self.assertEqual(m2,m)
            self.assertAlmostEqual(0, Transition[0,0],places=3)
            self.assertAlmostEqual(0.01, Transition[0,1],places=2)
            self.assertAlmostEqual( 0.819 , Transition[0,2],delta=0.01)
            self.assertAlmostEqual(  0.172, Transition[0,3],delta=0.01)
            self.assertAlmostEqual(0.333, Transition[7,7],places=2)
            self.assertAlmostEqual(0.333, Transition[7,8],places=2)
            self.assertAlmostEqual(0.333, Transition[7,9],places=2)
            self.assertAlmostEqual(0.01, Emission[2,1],places=3)

        @skip('TBP')
        def test_ba10g(self):
            '''
            BA10G   Sequence Alignment with Profile HMM Problem
            '''
            Transition, Emission,StateNames = ConstructProfileHMM(0.4,
                                                                  'ABCDEF',
                                                                  ['ACDEFACADF',
                                                                   'AFDA---CCF',
                                                                   'A--EFD-FDC',
                                                                   'ACAEF--A-C',
                                                                   'ADDEFAAADF'],
                                                                  sigma=0.01)
            Path                            = AlignSequence('AEFDFDC','ABCDEF',Emission)


        def test_ba10h(self):
            Transitions,Emissions = EstimateParameters('yzzzyxzxxx','xyz','BBABABABAB','ABC')
            self.assertEqual(1, Transitions[0,1])
            self.assertEqual(0.8, Transitions[1,0])
            self.assertEqual(0.2, Transitions[1,1])
            self.assertAlmostEqual(1/3, Transitions[2,0])
            self.assertAlmostEqual(1/3, Transitions[2,1])
            self.assertAlmostEqual(1/3, Transitions[2,2])
            self.assertEqual(0.25, Emissions[0,0])
            self.assertEqual(0.25, Emissions[0,1])
            self.assertEqual(0.5, Emissions[0,2])
            self.assertEqual(0.5, Emissions[1,0])
            self.assertAlmostEqual(1/6, Emissions[1,1])
            self.assertAlmostEqual(1/3, Emissions[1,2])
            self.assertAlmostEqual(1/3, Emissions[2,0])
            self.assertAlmostEqual(1/3, Emissions[2,1])
            self.assertAlmostEqual(1/3, Emissions[2,2])


        def test_ba10i_sample(self):
            Transition, Emission = ViterbiLearning('xxxzyzzxxzxyzxzxyxxzyzyzyyyyzzxxxzzxzyzzzxyxzzzxyzzxxxxzzzxyyxzzzzzyzzzxxzzxxxyxyzzyxzxxxyxzyxxyzyxz',
                                                   'xyz',
                                                   'AB',
                                                   np.array([ [0.582,   0.418],
                                                              [0.272,   0.728]]),
                                                   np.array([[0.129,   0.35,    0.52],
                                                             [0.422,   0.151,   0.426]]))
            assert_array_almost_equal(np.array([ [0.875,   0.125],
                                                 [0.011,   0.989]]),
                                      Transition,
                                      decimal = 2)
            assert_array_almost_equal(np.array([[0.0, 0.75 ,   0.25],
                                                [0.402,   0.174 ,  0.424]]),
                                      Emission,
                                      decimal = 2)


        def test_ba10i_sample(self):
            Transition, Emission = ViterbiLearning('xxxzyzzxxzxyzxzxyxxzyzyzyyyyzzxxxzzxzyzzzxyxzzzxyzzxxxxzzzxyyxzzzzzyzzzxxzzxxxyxyzzyxzxxxyxzyxxyzyxz',
                                                   'xyz',
                                                   'AB',
                                                   np.array([ [0.582,   0.418],
                                                              [0.272,   0.728]]),
                                                   np.array([[0.129,   0.35,    0.52],
                                                             [0.422,   0.151,   0.426]]))
            assert_array_almost_equal(np.array([ [0.875,   0.125],
                                                 [0.011,   0.989]]),
                                      Transition,
                                      decimal = 2)
            assert_array_almost_equal(np.array([[0.0, 0.75 ,   0.25],
                                                [0.402,   0.174 ,  0.424]]),
                                      Emission,
                                      decimal = 2)

        def test_ba10i_extra(self):
            Transition, Emission = ViterbiLearning('zzzyyxyzzzxyxzxxzzyxxzzzzzzyzxyzxxxzxxxzyzzxzxzzzxyzyyyxxxxzxyyyyyxzzzyxyzzxxzxxzxyxyxyzxzxzxzyxyzzz',
                                                   'xyz',
                                                   'AB',
                                                   np.array([ [0.436,	0.564],
                                                              [0.953,	0.047]]),
                                                   np.array([[0.367,	0.248,	0.385],
                                                             [0.401,	0.361,	0.238]]))
            assert_array_almost_equal(np.array([ [0,   1],
                                                 [1,   0]]),
                                      Transition,
                                      decimal = 2)
            assert_array_almost_equal(np.array([[0.36	,0.14,	0.5 ],
                                                [0.34,	0.34,	0.32]]),
                                      Emission,
                                      decimal = 2)

        def test_ba10j_extra(self):
            probabilities = SoftDecode('zyxxxxyxzz',
                                       'xyz',
                                       'AB',
                                       np.array([[0.911,   0.089],
                                                 [0.228,   0.772]]),
                                       np.array([[0.356,   0.191,   0.453 ],
                                                 [0.04,    0.467,   0.493]]))
            expected      = np.array([[0.5438,  0.4562 ],
                                      [0.6492,  0.3508 ],
                                      [0.9647,  0.0353 ],
                                      [0.9936,  0.0064 ],
                                      [0.9957,  0.0043 ],
                                      [0.9891,  0.0109 ],
                                      [0.9154,  0.0846 ],
                                      [0.964,   0.036  ],
                                      [0.8737,  0.1263 ],
                                      [0.8167,  0.1833 ]
             ])
            assert_array_almost_equal(expected,probabilities,decimal=4)

        def test_ba10j_extra(self):
            probabilities = SoftDecode('xyyzxzyxyy',
                                       'xyz',
                                       'ABCD',
                                       np.array([[0.401,	0.009,	0.195,	0.396		],
                                                 [0.375,	0.237,	0.269,	0.119	],
                                                 [0.283,	0.25,	0.259,	0.207	],
                                                 [0.108,	0.529,	0.107,	0.256	]]),
                                       np.array([[0.414,	0.335,	0.251	 ],
                                                 [0.233,	0.172,	0.596	],
                                                 [0.284,	0.355,	0.361	],
                                                 [0.028,	0.638,	0.334]])).probabilities
            expected      = np.array([[0.5003,	0.2114,	0.2662,	0.0220],
                                      [0.3648,	0.053,	0.1909,	0.3913	],
                                      [0.1511,	0.1251,	0.1553,	0.5685],
                                      [0.1297,	0.5359,	0.1542,	0.1802	],
                                      [0.4414,	0.2628,	0.2673,	0.0285],
                                      [0.3031,	0.2213,	0.2339,	0.2417	],
                                      [0.2789,	0.1536,	0.2139,	0.3537],
                                      [0.5088,	0.269,	0.1975,	0.0247	],
                                      [0.3695,	0.0578,	0.1978,	0.3748],
                                      [0.2231,	0.1356,	0.1658,	0.4755]
             ])

            assert_array_almost_equal(expected,probabilities,decimal=3)


        def test_ba10k_responsibility(self):
            '''
            This test has been written to help me understand the responsibility matrices introduced on page 224
            '''
            def get_weight(i,l,k,x, Path, Transition,Emission):
                return Transition[Path[i-1],Path[i]] * Emission[Path[i],x[i]]  #page 192
            Transition               = np.array([[9/10, 1/10],   # Figure 10.5
                                                 [1/10, 9/10]])
            Emission                 = np.array([[1/2, 1/2],      # Figure 10.5
                                                 [3/4,1/4]])
            x                        = 'THTHHHTHTTH'   # Figure 10.26
            Alphabet                 = 'HT'    # Figure 10.5
            States                   = 'FB'     # Figure 10.5
            xindices                 = get_indices(x,Alphabet=Alphabet)
            decoded                  = SoftDecode(x,Alphabet, States,Transition,Emission)
            responsibilities,forward,backward = decoded.probabilities,decoded.forward,decoded.backward
            edge_responsibilities =  get_edge_responsibilities(forward,backward, Transition, Emission, xindices)
            assert_array_almost_equal(np.array([[0.636, 0.364],     # Figure 10.26
                                                [0.593, 0.407],
                                                [0.600, 0.400],
                                                [0.533, 0.467],
                                                [0.515, 0.485],
                                                [0.544, 0.456],
                                                [0.627, 0.373],
                                                [0.633, 0.367],
                                                [0.692, 0.308],
                                                [0.686, 0.314],
                                                [0.609, 0.391]
                                                ]),
                                      responsibilities,
                                      decimal = 3)

            assert_array_almost_equal([[0.562,0.074],[0.031,0.333]],edge_responsibilities[0,:,:],decimal=3)
            assert_array_almost_equal([[0.523,0.022],[0.104,0.351]],edge_responsibilities[5,:,:],decimal=3)
            assert_array_almost_equal([[0.588,0.098],[0.022,0.293]],edge_responsibilities[-1,:,:],decimal=3)


        def test_ba10k_sample(self):
            Transition, Emission = BaumWelch('xzyyzyzyxy',
                                             'xyz',
                                             'AB',
                                             np.array([[0.019 ,  0.981],
                                                       [ 0.668,   0.332 ]]),
                                             np.array([[ 0.175,   0.003,   0.821  ],
                                                       [ 0.196,   0.512 ,  0.293]]))

            assert_array_almost_equal(np.array([[0.000, 1.000],
                                                [0.786,   0.214]]),
                                      Transition,
                                      decimal = 3)

            assert_array_almost_equal(np.array([[0.242,   0.000,   0.758],
                                                [0.172,   0.828,  0.000]]),
                                      Emission,
                                      decimal = 3)

        def test_ba10k_extra(self):
            Transition, Emission = BaumWelch('yzzxzzzyxyxzyxzzyyzzxxzyzyyyyyyxzyxxyzzzyzxyxxxxyxzzzzyzxyxxzyyyyxyzyxzzyzyxyxzzxyxxxzxyxyyxzxxyzxzz',
                                             'xyz',
                                             'ABCD',
                                             np.array([[0.056,	0.443,	0.334,	0.167	],
                                                       [0.826,	0.052,	0.011,	0.111	 ],
                                                       [0.022,	0.163,	0.811,	0.004	],
                                                       [0.141,	0.696,	0.055,	0.108	]]),
                                             np.array([[0.340,	0.084,	0.576	  ],
                                                       [ 0.196,	0.592,	0.212	],
                                                       [0.344,	0.609,	0.047	],
                                                       [0.320,	0.487,	0.193]]),
                                             N=100)

            assert_array_almost_equal(np.array([[0.000,	0.384,	0.186,	0.430],
                                                [0.989,	0.000,	0.000,	0.011],
                                                [0.000,	0.000,	0.596,	0.404	],
                                                [0.403,	0.596,	0.000,	0.001	]]),
                                      Transition,
                                      decimal = 3)

            assert_array_almost_equal(np.array([[0.447,	0.000,	0.553],
                                                [0.032,	0.648,	0.321	],
                                                [0.000,	1.000,	0.000	],
                                                [0.714,	0.000,	0.286]]),
                                      Emission,
                                      decimal = 3)


    main()
