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

def get_indices(S,Alphabet='AB'):
    '''
    Used to convert a string into integers representing the position of each charcter in the alphabet.

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

def get_reduced(s,mask):
    return [s[i] for i in range(len(s)) if mask[i]]

def is_space(ch):
    '''
    Check whether chatacter is a space (represented by hyphen)
    '''
    return ch=='-'

class StateSet:
    '''
    Manage states for, S, Mi, Ii, Di, E
    '''
    def __init__(self,m):
        self.m = m

    def get_pair(self,index):
        if index   == 0:   return ('S',None)
        if index == self.m-1: return ('E',None)
        if index%3 == 0:   return ('D', index//3)
        if index%3 == 1:   return ('I', index//3)
        if index%3 == 2:   return ('M', index//3+1)

    def get_index(self,state):
        function,i = state
        if function == 'S': return 0
        if function == 'E': return self.m-1
        if function == 'D': return 3*i
        if function == 'I': return 3*i+1
        if function == 'M': return 3*i-1


    def get_key(self,op,seq):
        return f'{op}-{seq}'

    def split_key(self,key):
        parts =key.split('-')
        if parts[0] in ['M','D','I']:
            return parts[0],int(parts[1])
        else:
            return parts[0],parts[1]

    @classmethod
    def create_path(cls,s,mask):
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

def ConstructProfileHMM(theta,Alphabet,Alignment,sigma=0):
    '''
    ConstructProfileHMM

    BA10E Construct a Profile HMM and
    BA10F Construct a Profile HMM with Pseudocounts

    Parameters:
        theta      Threshold.
        Alphabet
        Alignment
        sigma      Minimum probability assigned to legal transitions
   '''

    def create_mask():
        '''
        Construct a mask to exclude columns from an alignment
        if the freaction of spaces exceeds theta
        '''
        def get_count(i):
            return sum([is_space(c) for s in Alignment for c in s[i]])

        fractions = [get_count(i)/len(Alignment) for i in range(len(Alignment[0]))]
        return [fractions[i]<theta for i in range(len(Alignment[0]))]



    def normalize_rows(A):
        '''
        Normalize a matrix so each row sums to 1, unless all elements of row are 0
        '''
        row_totals = A.sum(axis=1)

        for i in range(m):
            if row_totals[i]==0:
                row_totals[i]=1

        return A/row_totals.reshape(m,1)

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

def float2str(x,
              precision=3,
              p0=0,
              p1=1):
    '''
     float2str

     Format a floatinmg point number in the manner that the
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

def format_state(s):
    code,index = s
    if code in ['S','E']:
        return code
    else:
        return f'{code}{index}'


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
        yield format_state(States[i]) + '\t' + '\t'.join(row)


def formatTransition(Transition,States,precision=2):
    '''
    formatTransition
    '''
    m,n = Transition.shape
    yield  '\t'.join(format_state(state) for state in States)
    for i in range(m):
        row = []
        for j in range(n):
            row.append(float2str(Transition[i,j],precision))
        yield format_state(States[i]) + '\t' + '\t'.join(row)


def AlignSequenceWithProfileHMM(theta,Alphabet,Alignment,sigma=0):
    '''
    AlignSequenceWithProfileHMM

    BA10G   Sequence Alignment with Profile HMM Problem

    Given: A string Text, a multiple alignment Alignment, a threshold θ, and a pseudocount σ.

    Return: An optimal hidden path emitting Text in HMM(Alignment,θ,σ).
    '''
    pass


def EstimateParameters(s,Alphabet,path,States):
    '''BA10H 	Estimate the Parameters of an HMM'''
    def create_Transitions():
        Transitions = np.zeros((len(States),len(States))) #{(i,j):0 for i in States for j in States}
        for i in range(1,n):
            Transitions[path[i-1],path[i]]+= 1
        Sums = {i: sum(Transitions[(i,j)] for j in States) for i in States}
        for i in States:
            for j in States:
                if Sums[i]>0:
                    Transitions[i,j]/= Sums[i]
                else:
                    Transitions[i,j] = 1/len(States)

        return Transitions

    def create_Emissions():
        Emissions   = {(ch,state): 0 for state in States for ch in Alphabet}
        for i in range(n):
            Emissions[(s[i],path[i])]+= 1
        Sums = {j:sum(Emissions[(ch,j)] for ch in Alphabet) for j in States }
        for ch in Alphabet:
            for j in States:
                if Sums[j]>0:
                    Emissions[(ch,j)]/= Sums[j]
                else:
                    Emissions[(ch,j)] = 1/len(States)
        return Emissions

    n = len(path)
    assert n==len(s)
    path_indices=get_indices(path,States)
    return create_Transitions(),create_Emissions()



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

        # def test_ba10g(self):
            # '''
            # BA10G   Sequence Alignment with Profile HMM Problem
            # '''
            # Path = AlignSequenceWithProfileHMM(0.4,
                                                # 'ABCDEF',
                                                # ['ACDEFACADF',
                                                 # 'AFDA---CCF',
                                                 # 'A--EFD-FDC',
                                                 # 'ACAEF--A-C',
                                                 # 'ADDEFAAADF'],
                                                # sigma=0.01)
            # self.assertEqual(9,len(Path))

        # def test_ba10h(self):
            # Transitions,Emissions = EstimateParameters('yzzzyxzxxx','xyz','BBABABABAB','ABC')
            # x=0



    main()
