#    Copyright (C) 2019-2020 Greenweaves Software Limited
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
#
#    Hidden Markov Models

from numpy import argmax
from math  import isqrt

# BA10A Compute the Probability of a Hidden Path

# ProbabilityHiddenPath
#
# Given: A hidden path followed by the states States and transition matrix Transition of an HMM.
#
# Return: The probability of this path. You may assume that initial probabilities are equal.

def ProbabilityHiddenPath(path,Transition):
    n      = isqrt(len(Transition.keys()))  # Number of states: Transition matrix is n by n
    result = 1/n                            # Initial probability for each state 
    for i in range(1,len(path)):
        result *= Transition[(path[i-1],path[i])]
    return result

# BA10B Compute the Probability of an Outcome Given a Hidden Path 

# ProbabilityOutcomeGivenHiddenPath

# Given: A string x, followed by the alphabet from which x was constructed, followed by a hidden path, followed by the states  and emission matrix Emission of an HMM.

# Return: The conditional probability  that string x will be emitted by the HMM given the hidden path

def ProbabilityOutcomeGivenHiddenPath(string,path,Emission):
    result = 1
    for i in range(0,len(path)):
        result *= Emission[(path[i],string[i])]
    return result 

# BA10C Implement the Viterbi Algorithm 

# Viterbi
#
# Input: A string x, followed by the alphabet from which x was constructed, 
#        followed by the states States, transition matrix Transition,
#        and emission matrix Emission of an HMM.
#
# Return: A path that maximizes the (unconditional) probability Pr(x, p) over all possible paths p

def Viterbi(xs,alphabet,States,Transition,Emission):
    
    # calculateproduct_weights
    #
    # Calculate array of weights by position in string and state
    
    def calculateproduct_weights(s_source = 1):
        def product_weight(k,x):
            return max([s[-1][l] * Transition[(States[l],k)] * Emission[(k,x)] for l in range(len(States))])
        
        s = []  # first index is position, 2nd state
    
        s.append([s_source * (1/len(States)) * Emission[(k,xs[0])] for k in States])
        for x in xs[1:]:
            s.append([product_weight(k,x) for k in States])
        return s
    
    # backtrack
    # 
    # Find most likely path through state space by backtracking
    
    def backtrack(s):
        n     = len(s) - 1
        state = argmax(s[n])
        path  = [States[state]]
        while True:
            ps = [s[n-1][l] * Transition[(States[l],States[state])]  for l in range(len(States))]
            state = argmax(ps)
            path.append(States[state])
            n-=1
            if n<=0: return path[::-1]
    
    return ''.join(backtrack(calculateproduct_weights()))

#  BA10D 	Compute the Probability of a String Emitted by an HMM 
#
# Likelihood
#
# Given: A string x, followed by the alphabet from which x was constructed, followed by the states,
# transition matrix, and emission matrix  of an HMM 
#
# Return: The probability Pr(x) that the HMM emits x. 

def Likelihood(xs,Alphabet,States,Transition,Emission):
    def calculateproduct_weights(s_source = 1):
        def product_weight(k,x):
            return sum([s[-1][l] * Transition[(States[l],k)] * Emission[(k,x)] for l in range(len(States))])
        
        s = []  # first index is position, 2nd state
    
        s.append([s_source * (1/len(States)) * Emission[(k,xs[0])] for k in States])
        for x in xs[1:]:
            s.append([product_weight(k,x) for k in States])
        return s 
    return sum(calculateproduct_weights()[-1])


# ConstructProfileHMM
#
# ConstructProfileHMM
#
# Parameters:
#    theta      Threshold. This isn't the same as the theta in the textbook
#               See David Eccles and fanta's comments - http://rosalind.info/problems/ba10e/questions/
#    Alphabet
#    Alignment

def ConstructProfileHMM(theta,Alphabet,Alignment,sigma=0):
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