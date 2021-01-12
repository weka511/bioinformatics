#   Copyright (C) 2020 Greenweaves Software Limited

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

#  BA10E Construct a Profile HMM

import abc
import argparse
import os
import time
from   helpers import read_strings
from   rosalind import RosalindException

# ConstructProfileHMM
#
# ConstructProfileHMM
#
# Parameters:
#    theta      Threshold. This isn't the same as the theta in the textbook
#               See David Eccles and fanta's comments - http://rosalind.info/problems/ba10e/questions/
#    Alphabet
#    Alignment

def ConstructProfileHMM(theta,Alphabet,Alignment):
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
    yield '\t' + '\t'.join(Alphabet)
    for state in States:
        row = []    
        for symbol in Alphabet:
            probability = Emission[(state,symbol)] if (state,symbol) in Emission else 0  
            row.append(float2str(probability,precision))
        yield format_state(state) + '\t' + '\t'.join(row)            
 

# formatTransition


    
def formatTransition(Transition,States,precision=2):
    yield '\t' + '\t'.join(format_state(state) for state in States)
    for state1 in States:
        row = []
        for state2 in States:
            probability = Transition[(state1,state2)] if (state1,state2) in Transition else 0
            row.append(float2str(probability,precision))
        yield format_state(state1) + '\t' + '\t'.join(row)
        
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA10E Construct a Profile HMM')
    parser.add_argument('--sample',    default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',     default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind',  default=False, action='store_true', help='process Rosalind dataset')
    parser.add_argument('--text',      default=False, action='store_true', help='process dataset from textbook')
    parser.add_argument('--precision', default=3,                          help='Controls display of probabilities')

    args = parser.parse_args()
    if args.sample:
        States,Transition,Emission = ConstructProfileHMM(0.289,
                                                         ['A',   'B',   'C',   'D',   'E'],
                                                         ['EBA', 'EBD', 'EB-', 'EED', 'EBD', 'EBE','E-D','EBD'])
        
        for row in formatTransition(Transition,States,precision=args.precision):
            print (row)
        print ('--------')
        for row in formatEmission(Emission,States,['A',   'B',   'C',   'D',   'E'],precision=args.precision):
            print (row)
            
    if args.text:
        Alphabet = ['A',   'C',   'D',   'E', 'F']
        States,Transition,Emission = ConstructProfileHMM(0.35,
                                                         Alphabet,
                                                         ['ACDEFACADF',
                                                          'AFDA---CCF',
                                                          'A--EFD-FDC',
                                                          'ACAEF--A-C',
                                                          'ADDEFAAADF'])
        
        for row in formatTransition(Transition,States,precision=args.precision):
            print (row)
        print ('--------')
        for row in formatEmission(Emission,States,Alphabet,precision=args.precision):
            print (row)
            
            
    if args.extra:
        Input,Expected             = read_strings(f'data/ProfileHMM.txt',init=0)
        States,Transition,Emission = ConstructProfileHMM(float(Input[0]),
                                                         Input[2].split(),
                                                         Input[4:-1]) 
        for row in formatTransition(Transition,States,precision=args.precision):
            print (row)
        print ('--------')
        for row in formatEmission(Emission,States,Input[2].split(),precision=args.precision):
            print (row)
               
    if args.rosalind:
        Input                      = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        States,Transition,Emission = ConstructProfileHMM(float(Input[0]),
                                                         Input[2].split(),
                                                         Input[4:-1]) 

        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for row in formatTransition(Transition,States,precision=args.precision):
                print (row)
                f.write(f'{row}\n')
            print ('--------')
            f.write('--------\n')
            for row in formatEmission(Emission,States,Input[2].split(),precision=args.precision):
                print (row)            
                f.write(f'{row}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
