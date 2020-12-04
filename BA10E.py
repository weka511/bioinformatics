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

def ConstructProfileHMM(theta,Alphabet,Alignment):
    def CountChars(m,j,K):
        RowCounts = [0]*K
        for i in range(m):
            if Alignment[i][j] in Alphabet:
                RowCounts[Alphabet.index(Alignment[i][j])]+=1
        return RowCounts
    
    def CreateStates(m):
        Product = [('S',0),('I',0)]
        for i in range(m):
            Product.append(('M',i+1))
            Product.append(('D',i+1))
            Product.append(('I',i+1))
        Product.append(('E',m+1))
        return Product
    
    n       = len(Alignment[0])
    m       = len(Alignment)
    K       = len(Alphabet)
    Counts  = [CountChars(m,j,K) for j in range(n)]
    Profile = [[c/sum(Count) for c in Count] for Count in Counts if  sum(Count) > (1-theta)*K]
    States  = CreateStates(len(Profile))
    x=0
    
def old():    
    def CreateStates(m):
        Product = [('S',0),('I',0)]
        for i in range(m):
            Product.append(('M',i+1))
            Product.append(('D',i+1))
            Product.append(('I',i+1))
        Product.append(('E',m+1))
        return Product

    def get_index(operation,pos):
        for i in range(len(States)):
            if (operation,pos)==States[i]:
                return i
    
    def normalize(Row):
        total = sum(Row)
        if total>0:
            for i in range(len(Row)):
                Row[i]/=total
            
    def CreateCounts(States,n=0,L=0,m=0,K=0):
        StateCounts = [] # row--state HMM is in, column = P(move to)
        AlphaCounts = [] # row--state HMM is in, column = P emit this symbol
        for i in range(L): # iterate through states
            StateCounts.append([0]*L)
            AlphaCounts.append([0]*K)
            operation,pos = States[i]
            for j in range(n):
                operation,pos = States[i]
                if operation == 'S':
                    ch = Alignment[j][pos]
                    k  = get_index('M',pos+1) if ch in Alphabet else get_index('I',pos)
                    StateCounts[-1][k]+=1
                elif operation =='M':
                    ch = Alignment[j][pos-1]
                    if ch in Alphabet:
                        k  = Alphabet.index(ch)
                        AlphaCounts[-1][k]+=1
                    x=0
                elif operation =='D':
                    pass
                elif operation =='I':
                    pass
                else:      #E
                    pass
            normalize(StateCounts[-1])
            normalize(AlphaCounts[-1])
            x=0
 
        return StateCounts,AlphaCounts
   
    m                       = len(Alignment[0])
    n                       = len(Alignment)
    K                       = len(Alphabet)
    States                  = CreateStates(m)
    L                       = len(States)
    StateCounts,AlphaCounts = CreateCounts(States,n=n,L=L,m=m,K=K)
    return States,StateCounts,AlphaCounts

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
