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
from unittest import TestCase, main

import numpy as np


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
        code,index = s
        if code in ['S','E']:
            return code
        else:
            return f'{code}{index}'


def get_indices(S,Alphabet='AB'):
    '''
    Used to convert a string into integers representing the position of each character in the alphabet.

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

    Given: A string x, followed by the alphabet from which x was constructed, followed by a hidden path, followed by the states  and emission matrix Emission of an HMM.

    Return: The conditional probability  that string x will be emitted by the HMM given the hidden path
    '''

    logEmission = np.log(Emission)
    return np.exp(sum([logEmission[a,b] for a,b in zip(get_indices(path,Alphabet=states),
                                                       get_indices(string,Alphabet=alphabet))]))


def Viterbi(xs,alphabet,States,Transition,Emission):
    '''
    Viterbi

    Solve:    BA10C Implement the Viterbi Algorithm

    Input: A string x, followed by the alphabet from which x was constructed,
        followed by the states States, transition matrix Transition,
        and emission matrix Emission of an HMM.

    Return: A path that maximizes the (unconditional) probability Pr(x, p) over all possible paths p
    '''

    def create_weights(s_source = 1):
        '''
        calculateproduct_weights

        Calculate array of weights by position in string and state

        '''
        def get_weight(state,xs,i):
            return np.multiply(s[i-1,:],Transition[:,state]* Emission[state,xs[i]]).max()

        m         = len(xs)
        n         = len(States)
        x_indices = get_indices(xs,alphabet)
        s         = np.zeros((m,n))
        s[0,:]    = s_source  * Emission[:,x_indices[0]] / n
        for i in range(1,len(x_indices)):
            s[i,:] = np.array([get_weight(state,x_indices,i) for state in range(n)])

        return s

    def find_most_likely_path(s):
        '''
        backtrack

        Find most likely path through state space by backtracking

        '''
        m,_    = s.shape
        m     -= 1
        state = np.argmax(s[m])
        path  = [States[state]]
        for i in range(m-1,-1,-1):
            ps    = [s[i,j] * Transition[j,state]  for j in range(len(States))]
            state = np.argmax(ps)
            path.append(States[state])
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
        Alphabet   Symbols sed for Alignment
        Alignment  A list of strings representing aligned proteins
        sigma      Minimum probability assigned to legal transitions
   '''

    def create_mask():
        '''
        Construct a mask to exclude columns from an alignment
        if the fraction of spaces is greater than or equal to theta
        '''
        def get_count(i):
            '''
            Returns:
                Number of spacen in column i
            '''
            return sum([is_space(c) for s in Alignment for c in s[i]])

        fractions = [get_count(i)/len(Alignment) for i in range(len(Alignment[0]))]
        return [fraction<theta for fraction in fractions]


    def create_transition(m,Paths):
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


    def create_emission(m,n,Paths):
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
    return create_transition(m,Paths), create_emission(m,n,Paths), [state_set.get_pair(i) for i in range(m)]


def SoftDecode(s, Alphabet,  States, Transition, Emission):
    result = np.zeros((len(s),len(States)))
    return result


class ViterbiGraph:

    def __init__(self, Transition, States):
        def format(state):
            a,b = state
            return f'{a}' if b==None else  f'{a}{b}'
        self.Transition = Transition
        self.m,_        = Transition.shape
        self.k          = self.m//3 - 1
        self.state_index = {format(state):i for i,state in enumerate(States)}
        self.n           = None

    def __getitem__(self,index):
        '''
        Implement mapping from rows and columns to states that is depicted in Figure 10.21
        '''
        i,j = index
        if self.n != None and j>self.n:
            return 'E'
        if i%3==0:
            return 'S' if j==0 else f'I{i//3}'
        if i%3==1: return f'M{i//3+1}'
        if i%3==2: return f'D{i//3+1}'

    def get_successors(self,i,j):
        '''
        Implement mapping from rows and columns to successors that is depicted in Figure 10.21
        '''
        def should_include(node):
            i_1,j_1 = node
            if self.n != None and j_1>self.n:
                return False
            return True
        state = self[i,j]
        result = None
        if self.n != None and j>self.n:
            pass
        else:
            match state[0]:
                case 'S':
                    result = [(2,0), # I
                              (1,0), # M
                              (1,1)] #D
                case 'D':
                    if i<3*self.k -1:
                        result = [(i+3,j), #D
                                  (i+1,j+1), # I
                                  (i+2,j+1)] # M
                    else:
                        result = [(0,j+1)] # I
                case 'M':
                    if i<3*self.k -2:
                        result = [(i+3,j+1), #M
                                  (i+2,j+1), # I
                                  (i+4,j)] # D
                    else:
                        result = [(i+2,j+1)] # I

                case 'I':
                    if i<3*self.k:
                        result =[(i+2,j), #D
                                (i,j+1), # I
                                (i+1,j+1)] # M
                    else:
                        result = [(0,j+1)] # I

                case 'E':
                    pass

        return [node for node in result if should_include(node)]

    def set_text_length(self,n):
        '''Used to establish last column of table'''
        self.n = n

    def AlignSequence(self,string,Alphabet,Emission):
        '''
        AlignSequenceWithProfileHMM

        BA10G   Sequence Alignment with Profile HMM Problem

        Given: A string Text, a multiple alignment Alignment, a threshold θ, and a pseudocount σ.

        Return: An optimal hidden path emitting Text in HMM(Alignment,θ,σ).
        '''
        def get_weight(state,xs,i):
            return np.multiply(s[i-1,:],self.Transition[:,state]* Emission[state,xs[i]]).max()
        self.set_text_length(len(string))
        start_state = self[0,0]
        start_index = self.state_index[start_state]
        Paths       = [j for j in range(self.m) if self.Transition[start_index,j]!=0]
        z=0



def EstimateParameters(s,Alphabet,path,States):
    '''BA10H 	Estimate the Parameters of an HMM'''
    def create_Transitions(path_indices,n):
        Transitions = np.zeros((len(States),len(States)))
        for i in range(1,n):
            Transitions[path_indices[i-1],path_indices[i]]+= 1
        Transitions[-1,:] = 1

        return normalize_rows(Transitions)

    def create_Emissions(path_indices,n):
        str_indices = get_indices(s,Alphabet)
        Emissions   = np.zeros((len(States),len(Alphabet)))
        for i in range(n):
            Emissions[path_indices[i],str_indices[i]]+= 1
        Emissions[-1,:] = 1

        return normalize_rows(Emissions)

    n = len(path)
    assert n==len(s)
    path_indices=get_indices(path,States)
    return create_Transitions(path_indices,n),create_Emissions(path_indices,n)



if __name__=='__main__':
    class Test_10_HMM(TestCase):

        def test_ba10a(self):
            '''BA10A Compute the Probability of a Hidden Path'''

            self.assertAlmostEqual(5.01732865318,
                                        1e19*ProbabilityHiddenPath('AABBBAABABAAAABBBBAABBABABBBAABBAAAABABAABBABABBAB',
                                             'AB',
                                             np.array([[ 0.194,   0.806],[ 0.273,   0.727]])),
                                        places=5)

        def test_ba10b(self):
            '''BA10B Compute the Probability of an Outcome Given a Hidden Path'''
            self.assertAlmostEqual(1.93157070893,
                                   1e28*ProbabilityOutcomeGivenHiddenPath('xxyzyxzzxzxyxyyzxxzzxxyyxxyxyzzxxyzyzxzxxyxyyzxxzx',
                                                                          'xyz',
                                                                          'BBBAAABABABBBBBBAAAAAABAAAABABABBBBBABAABABABABBBB',
                                                                          'AB',
                                                                          np.array([[0.612,   0.314,   0.074 ],
                                                                                   [0.346,   0.317 ,  0.336]])),
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
                                                                  ['EBA', 'EBD', 'EB-', 'EED', 'EBD', 'EBE', 'E-D','EBD'])
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

        def test_ba10j(self):
            probabilities = SoftDecode('zyxxxxyxzz',
                                       'xyz',
                                       'AB',
                                       np.array([[0.911,   0.089],
                                                 [0.228,   0.772]]),
                                       np.array([[0.356,   0.191,   0.453 ],
                                                 [0.04, 0.467, 0.493]]))
            expected      = np.array([
                 [0.5438,  0.4562 ],
                 [0.6492,  0.3508 ],
                 [0.9647,  0.0353 ],
                 [0.9936,  0.0064 ],
                 [0.9957,  0.0043 ],
                 [0.9891,  0.0109 ],
                 [0.9154,  0.0846 ],
                 [0.964,   0.036 ],
                 [0.8737,  0.1263 ],
                 [0.8167,  0.1833  ]
             ])
            m1,n1 = expected.shape
            m2,n2 = probabilities.shape
            self.assertEqual(m1,m2)
            self.assertEqual(n1,n2)
            for i in range(m1):
                self.assertEqual(expected[i,0],probabilities[i,0])
                self.assertEqual(expected[i,1],probabilities[i,1])

    class Test_10_Viterbi(TestCase):
        '''
        BA10G   Sequence Alignment with Profile HMM Problem
        test the ViterbiGraph class
        '''
        def setUp(self):
            Transition, Emission,StateNames = ConstructProfileHMM(0.35,
                                                                  'ACDEF',
                                                                  ['ACDEFACADF',
                                                                   'AFDA---CCF',
                                                                   'A--EFD-FDC',
                                                                   'ACAEF--A-C',
                                                                   'ADDEFAAADF'])
            self.viterbi = ViterbiGraph(Transition,StateNames)

        def test_I(self):
            '''From Figure 10.21'''
            self.assertEqual('I0', self.viterbi[0,1])
            successors = self.viterbi.get_successors(0,1)
            self.assertEqual(3,len(successors))
            self.assertEqual((2,1),successors[0])
            self.assertEqual((0,2),successors[1])
            self.assertEqual((1,2),successors[2])
            self.assertEqual('I8', self.viterbi[8*3,1])
            successors = self.viterbi.get_successors(8*3,1)
            self.assertEqual(1,len(successors))
            self.assertEqual((0,2),successors[0])

        def test_S(self):
            '''From Figure 10.21'''
            self.assertEqual('S',  self.viterbi[0,0])
            successors = self.viterbi.get_successors(0,0)
            self.assertEqual(3,len(successors))
            self.assertEqual((2,0),successors[0])
            self.assertEqual((1,0),successors[1])
            self.assertEqual((1,1),successors[2])

        def test_D(self):
            '''From Figure 10.21'''
            self.assertEqual('D1', self.viterbi[2,0])

            successors = self.viterbi.get_successors(2,0)
            self.assertEqual(3,len(successors))
            self.assertEqual((5,0),successors[0])
            self.assertEqual((3,1),successors[1])
            self.assertEqual((4,1),successors[2])
            self.assertEqual('D1', self.viterbi[2,1])
            self.assertEqual('D8', self.viterbi[2+7*3,1])
            successors = self.viterbi.get_successors(2+7*3,1)
            self.assertEqual(1,len(successors))
            self.assertEqual((0,2),successors[0])
            self.assertEqual('I0', self.viterbi[0,2])

        def test_M(self):
            '''From Figure 10.21'''
            self.assertEqual('M1', self.viterbi[1,1])
            successors = self.viterbi.get_successors(1,1)
            self.assertEqual(3,len(successors))
            self.assertEqual((4,2),successors[0]) #M
            self.assertEqual((3,2),successors[1]) #I
            self.assertEqual((5,1),successors[2]) #D
            self.assertEqual('M8', self.viterbi[1+7*3,0])
            successors = self.viterbi.get_successors(1+7*3,0)
            self.assertEqual(1,len(successors))
            self.assertEqual((3+7*3,1),successors[0])
            self.assertEqual('I8', self.viterbi[3+7*3,1])


        def test_Last_Column(self):
            '''From Figure 10.21'''
            self.viterbi.set_text_length(7)
            self.assertEqual('I0', self.viterbi[0,7])
            successors = self.viterbi.get_successors(0,7)
            self.assertEqual(1,len(successors))
            self.assertEqual((2,7),successors[0])
            self.assertEqual('M1', self.viterbi[1,7])
            successors = self.viterbi.get_successors(1,7)
            self.assertEqual(1,len(successors))
            self.assertEqual((5,7),successors[0])
            self.assertEqual('D2', self.viterbi[5,7])
            self.assertEqual('D1', self.viterbi[2,7])
            successors = self.viterbi.get_successors(2,7)
            self.assertEqual(1,len(successors))
            self.assertEqual((5,7),successors[0])
            self.assertEqual('D2', self.viterbi[5,7])

            self.assertEqual('I8', self.viterbi[8*3,1])
            successors = self.viterbi.get_successors(8*3,1)
            self.assertEqual(1,len(successors))
            self.assertEqual((0,2),successors[0])
            self.assertEqual('E', self.viterbi[0,8])

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
            viterbi                         = ViterbiGraph(Transition,StateNames)
            Path                            = viterbi.AlignSequence('AEFDFDC','ABCDEF',Emission)
            # self.assertEqual(9,len(Path))

    main()
