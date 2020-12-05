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

#  BA10E 	Construct a Profile HMM

import argparse
import os
import time
from   helpers import read_strings

# ConstructProfileHMM

def ConstructProfileHMM(theta,Alphabet,Alignment):
    def CountChars(m,j,K):
        RowCounts = [0]*K
        for i in range(m):
            if Alignment[i][j] in Alphabet:
                RowCounts[Alphabet.index(Alignment[i][j])]+=1
        return RowCounts
    
    def CreateStates(m):
        Product = [('S',0),('I',0)]
        max_pos = -1
        for i in range(m):
            Product.append(('M',i+1))
            Product.append(('D',i+1))
            Product.append(('I',i+1))
            max_pos = i+1
        Product.append(('E',None))
        return Product,max_pos
    
    def StateGenerator(States):
        for l in range(L): # iterate through states
            operation,pos = States[l]
            yield operation,pos
   
    def get_index_of_state(operation,pos):
        for i in range(len(States)):
            if (operation,pos)==States[i]:
                return i
            
    n              = len(Alignment[0])
    m              = len(Alignment)
    K              = len(Alphabet)
    Counts         = [CountChars(m,j,K) for j in range(n)]
    Conserved      = [sum(Count) > (1-theta)*K for Count in Counts]
    LinkBack       = [j for j in range(n) if Conserved[j]]
    Profile        = [[c/sum(Counts[i]) for c in Counts[i]] for i in range(len(Counts)) if  Conserved[i]]
    States,max_pos = CreateStates(len(Profile))
    L              = len(States)
    Transition     = [] # row--state HMM is in, column = P(move to)
    Emission       = [] # row--state HMM is in, column = P emit this symbol
    seq            = 0
    previous_del   = False
    for operation,pos in StateGenerator(States):
        if operation == 'S':
            Transition.append([1 if op=='M' and pos2==1 else 0 for op,pos2 in StateGenerator(States)] ) # FIXME - I0
            Emission.append([0]*K)   
        elif operation == 'M':
            next_row = [0]*L
            if pos<len(LinkBack):   #FIXME
                emission_count = sum(Counts[LinkBack[pos]])
                if emission_count>0:
                    next_row[seq+3] = emission_count/K
                skip_count     = K - emission_count
                if skip_count>0:
                    next_row[seq+4] = skip_count/K
                    previous_del   = True
            Transition.append(next_row) 
            Emission.append(Profile[pos-1])
        elif operation == 'D':
            next_row = [0]*L
            if previous_del:
                next_row[seq+2] = 1
                previous_del   = False
            Transition.append(next_row)
            Emission.append([0]*K)
        elif operation == 'I':
            print (pos,seq,Conserved[seq])
            Transition.append([0]*L)       # TODO
            Emission.append([0]*K)       # TODO        
        else:
            Transition.append([0]*L)
            Emission.append([0]*K)
        seq += 1
    return States,Transition,Emission

def formatState(state):
    operation,pos = state
    if operation in ['S','E']:
        return operation
    else:
        return f'{operation}{pos}'
    
def formatEmission(Emission,States,Alphabet):  
    yield '\t' + '\t'.join(Alphabet)
    for row,state in zip(Emission,States):
        yield formatState(state) + '\t' + '\t'.join(f'{r:.2f}' for r in row)

def formatTransition(Transition,States):
    yield '\t' + '\t'.join(formatState(s) for s in States)
    for row,state in zip(Transition,States):
        yield formatState(state) + '\t' + '\t'.join(f'{r:.2f}' for r in row)
        
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA10E 	Construct a Profile HMM')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        States,Transition,Emission = ConstructProfileHMM(0.289,
                                                         ['A',   'B',   'C',   'D',   'E'],
                                                         ['EBA', 'EBD', 'EB-', 'EED', 'EBD', 'EBE','E-D','EBD'])
        
        for row in formatTransition(Transition,States):
            print (row)
        print ('--------')
        for row in formatEmission(Emission,States,['A',   'B',   'C',   'D',   'E']):
            print (row)
               
    if args.extra:
        Input,Expected             = read_strings(f'data/ProfileHMM.txt',init=0)
        States,Transition,Emission = ConstructProfileHMM(float(Input[0]),
                                                         Input[2].split(),
                                                         Input[4:-1]) 
        for row in formatTransition(Transition,States):
            print (row)
        print ('--------')
        for row in formatEmission(Emission,States,Input[2].split()):
            print (row)
               
    if args.rosalind:
        Input                      = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        States,Transition,Emission = ConstructProfileHMM(float(Input[0]),
                                                         Input[2].split(),
                                                         Input[4:-1]) 

        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for row in formatTransition(Transition,States):
                print (row)
                f.write(f'{row}\n')
                print ('--------')
                f.write('--------\n')
            for row in formatEmission(Emission,States,Input[2].split()):
                print (row)            
                f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
