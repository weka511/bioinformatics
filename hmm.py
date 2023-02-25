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

def ConstructProfileHMM(theta,Alphabet,Alignment,sigma=0):
    '''
    ConstructProfileHMM

    BA10E Construct a Profile HMM

    Parameters:
        theta      Threshold. This isn't the same as the theta in the textbook
                   See David Eccles and fanta's comments - http://rosalind.info/problems/ba10e/questions/
                   Is this true?
        Alphabet
        Alignment
   '''

    def is_space(c):
        return c=='-'

    def create_mask():
        def get_count(i):
            return sum([is_space(c) for s in Alignment for c in s[i]])

        fractions = [get_count(i)/len(Alignment) for i in range(len(Alignment[0]))]
        return [fractions[i]<theta for i in range(len(Alignment[0]))]


    def get_state_name(i):
        if i   == 0: return ('S',None)
        if i   == m-1: return ('E',None)
        if i%3 == 0: return ('D', i//3)
        if i%3 == 1: return ('I', i//3)
        if i%3 == 2: return ('M', i//3+1)

    def get_state_index(state):
        function,i = state
        if function == 'S': return 0
        if function == 'E': return m-1
        if function == 'D': return 3*i
        if function == 'I': return 3*i+1
        if function == 'M': return 3*i-1

    def get_key(op,seq):
        return f'{op}-{seq}'

    def split_key(key):
        parts =key.split('-')
        if parts[0] in ['M','D','I']:
            return parts[0],int(parts[1])
        else:
            return parts[0],parts[1]

    def create_path(s,mask):
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

    def verify_constraints(matrix):
        m,_ = matrix.shape
        for i in range(m):
            row_total = matrix[i,:].sum()
            assert row_total==0 or row_total==1,f'Row {i} total = {row_total}'

    def create_transition(m,Paths):
        def create_census():
            States = {}
            for path in Paths:
                for (op1, seq1, _),(op2, seq2, _) in zip(path[:-1],path[1:]):
                    if get_key(op1,seq1)  not in States:
                        States[get_key(op1,seq1)] = []
                    States[get_key(op1,seq1)].append((op2, seq2))
            return States

        product = np.zeros((m,m))
        States  = create_census()
        for key1,successors in States.items():
            state1,index1 = split_key(key1)
            counts      = {}
            for succ,seq in successors:
                if get_key(succ,seq) not in counts:
                    counts[get_key(succ,seq)] = 0
                counts[get_key(succ,seq)]   += 1
            fractions    = {key:count/len(successors) for key,count in counts.items()}
            state_index1 = get_state_index((state1,index1))
            for key2,fraction in fractions.items():
                state2,index2 = split_key(key2)
                state_index2 = get_state_index((state2,index2))
                product[state_index1,state_index2] = fraction

        verify_constraints(product)
        return product

    def create_emission(m,n,Paths):
        def create_census():
            States = {}
            for path in Paths:
                for op, seq, ch in path:
                    if get_key(op,seq)  not in States:
                        States[get_key(op,seq)] = []
                    States[get_key(op,seq)].append(ch)
            return States

        product = np.zeros((m,n))
        States  = create_census()
        for key,chars in States.items():
            state,index = split_key(key)
            if state in ['M','I']:
                counts = {ch:0 for ch in Alphabet}
                for ch in chars:
                    counts[ch] += 1
                fractions = {ch:count/len(chars) for ch,count in counts.items()}
                state_index = get_state_index((state,index))
                for j in range(n):
                    product[state_index,j] = fractions[Alphabet[j]]

        verify_constraints(product)
        return product

    mask          = create_mask()
    Paths         = [create_path(s,mask) for s in Alignment]
    m             = 3 * (sum([1 for m in mask if m]) + 1)
    n             = len(Alphabet)
    return create_transition(m,Paths), create_emission(m,n,Paths), [get_state_name(i) for i in range(m)]

# ConstructProfileHMM
#
# ConstructProfileHMM
#
# Parameters:
#    theta      Threshold. This isn't the same as the theta in the textbook
#               See David Eccles and fanta's comments - http://rosalind.info/problems/ba10e/questions/
#    Alphabet
#    Alignment

def ConstructProfileHMM666(theta,Alphabet,Alignment,sigma=0):
    #   CountChars
    #
    #   Used to count alphabetical characters in specified column of alignment
    #
    #   Parameters:
    #       m     Number of sequences
    #       j     The column to be counted
    #       K     Number of symbols in Alphabet
    #
    #  Returns:
    #       Number of symbols from alphabet in column i

    def CountChars(m,j,K):
        Counts = [0]*K
        for i in range(m):
            if Alignment[i][j] in Alphabet:
                Counts[Alphabet.index(Alignment[i][j])]+=1
        return Counts

    def create_states(Conserved):
        Product        = [[] for column in Conserved if column]
        Product.append([])
        return Product

    def create_state_indices(Conserved):
        def get_symbol(State):
            return (State,index)

        index   = 0
        Product = [get_symbol('S'),get_symbol('I')]
        for column in Conserved:
            if not column: continue
            index += 1
            Product.append(get_symbol('M'))
            Product.append(get_symbol('D'))
            Product.append(get_symbol('I'))
        index += 1
        Product.append(get_symbol('E'))
        return Product

    def create_state_counts(StateIndices):
        Product  = {}
        for s,i in StateIndices:
            for t,j in StateIndices:
                if s=='S':
                    if t=='I' and j==0:
                        Product[(s,i),(t,j)] = 0
                    elif t in ['M','D'] and j==1:
                        Product[(s,i),(t,j)] = 0
                elif s=='I':
                    if i==0:
                        if t=='I' and j==0:
                            Product[(s,i),(t,j)] = 0
                        elif t in ['M','D'] and j==1:
                            Product[(s,i),(t,j)] = 0
                    else:
                        if (t=='I' and j==i) or (t in ['M','D'] and j==i+1):
                            Product[(s,i),(t,j)] = 0
                        if t=='E' and j==i+1:
                            Product[(s,i),(t,j)] = 0
                elif s in ['M','D']:
                    if (t=='I' and j==i) or (t in ['M','D'] and j==i+1):
                        Product[(s,i),(t,j)] = 0
                    if t=='E' and j==i+1:
                        Product[(s,i),(t,j)] = 0
                else:
                    assert(s=='E')

        return Product

    def create_emission_counts(StateIndices,Alphabet):
        Product  = {}
        for index in StateIndices:
            for ch in Alphabet:
                Product[(index,ch)] = 0
        return Product

    def create_state_frequencies(StateCounts,StateIndices):
        Totals  = {i:0 for i in StateIndices}

        for key,count in StateCounts.items():
            ((s,i),_) = key
            Totals[(s,i)] += count

        Product = {}
        for key,count in StateCounts.items():
            ((s,i),_) = key
            Product[key] = count/Totals[(s,i)] if Totals[(s,i)]>0 else 0
        return Product

    def create_emission_frequencies(EmissionCounts, StateIndices):
        Totals  = {i:0 for i in StateIndices}
        for key,count in EmissionCounts.items():
            ((s,i),_) = key
            Totals[(s,i)] += count
        Product = {}
        for key,count in EmissionCounts.items():
            ((s,i),_) = key
            Product[key] = count/Totals[(s,i)] if Totals[(s,i)]>0 else 0

        return Product

    #   Useful constants - lengths of arrays

    K              = len(Alphabet)      # Number of symbols in alphabet
    m              = len(Alignment)     # Number of strings in alignment
    n              = len(Alignment[0])  # Number of symbols in each alignment
    for Sequence in Alignment[1:]:      # All sequences should be the same length
        assert(n == len(Sequence))

    #   construct profile - Number of symbols from alphabet in each column

    Counts         = [CountChars(m,j,K) for j in range(n)]

    #   Indicate whether or not symbols in column are over threshold
    #   If theta is maximum proportion of deleted symbols, 1-theta is
    #   maximum number of conserved.

    Conserved      = [sum(Count) > (1-theta)*K for Count in Counts]

    column_count   = sum(1 for column in Conserved if column)

    Merges         = create_states(Conserved)
    Inserts        = create_states(Conserved)
    Deletes        = create_states(Conserved)
    StateIndices   = create_state_indices(Conserved)
    StateCounts    = create_state_counts(StateIndices)
    EmissionCounts = create_emission_counts(StateIndices,Alphabet)

    Runs        = []
    for Sequence in Alignment:
        previous =  'S'
        States   = [previous]
        j = 0
        for i in range(n):
            ch = Sequence[i]
            if Conserved[i]:
                j+=1
                if ch in Alphabet:
                    Merges[j].append((previous,ch))
                    previous =  'M'
                    States.append(previous)
                elif ch == '-':
                    Deletes[j].append((previous,ch))
                    previous =  'D'
                    States.append(previous)
                else:
                    raise RosalindException(f'Invalid {ch}')
            else:
                if ch in Alphabet:
                    Inserts[j].append((previous,ch))
                    previous =  'I'
                    States.append(previous)
                elif ch == '-':
                    pass
                else:
                    raise RosalindException(f'Invalid {ch}')
        States.append('E')
        Runs.append(States)
        print (f'Counting {Sequence} {"".join(States)}')
        index     = 0
        previous  = None
        seq_index = 0
        for state in States:
            while seq_index < len(Sequence)-1 and Sequence[seq_index]=='-':
                seq_index += 1
            if state in ['M','D','E']:
                index = index+1
            if previous != None:
                StateCounts[previous,(state,index)] += 1
                if state in ['M','I']:
                    EmissionCounts[(state,index),Sequence[seq_index]] += 1
                    seq_index += 1
            previous = (state, index)

    return (StateIndices,
            create_state_frequencies(StateCounts, StateIndices),
            create_emission_frequencies(EmissionCounts, StateIndices))



# float2str
#
# Format a floatinmg point number for display
#
# Parameters:
#     x         Value to be displayed
#     precision Number of digits (after decimal point)

def float2str(x,precision=2,p0=0,p1=1):
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

# formatEmission

def formatEmission(Emission,States,Alphabet,precision=2):
    yield '\t'.join(Alphabet)
    for state in States:
        row = []
        for symbol in Alphabet:
            probability = Emission[(state,symbol)] if (state,symbol) in Emission else 0
            row.append(float2str(probability,precision))
        yield format_state(state) + '\t' + '\t'.join(row)


# formatTransition



def formatTransition(Transition,States,precision=2):
    yield  '\t'.join(format_state(state) for state in States)
    for state1 in States:
        row = []
        for state2 in States:
            probability = Transition[(state1,state2)] if (state1,state2) in Transition else 0
            row.append(float2str(probability,precision))
        yield format_state(state1) + '\t' + '\t'.join(row)

#  BA10H 	Estimate the Parameters of an HMM

def EstimateParameters(s,Alphabet,path,States):
    def create_Transitions():
        Transitions = {(i,j):0 for i in States for j in States}
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
    return create_Transitions(),create_Emissions()

if __name__=='__main__':
    class Test_10_HMM(TestCase):

        def test_ba10e1(self): #BA10E Construct a Profile HMM
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

    main()
