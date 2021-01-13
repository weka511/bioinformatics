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

#  BA10H 	Estimate the Parameters of an HMM 

import argparse
import os
import time
from   helpers import read_strings
from   hmm     import float2str

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

# formatEmission

def formatEmission(Emission,States,Alphabet,precision=2):  
    yield '\t'.join(Alphabet)
    for state in States:
        row = []    
        for symbol in Alphabet:
            probability = Emission[(symbol,state)]  
            row.append(float2str(probability,precision))
        yield state + '\t' + '\t'.join(row)            
 

# formatTransition


    
def formatTransition(Transition,States,precision=2):
    yield  '\t'.join(state for state in States)
    for state1 in States:
        row = []
        for state2 in States:
            probability = Transition[(state1,state2)] if (state1,state2) in Transition else 0
            row.append(float2str(probability,precision))
        yield state1 + '\t' + '\t'.join(row)

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA10H 	Estimate the Parameters of an HMM ')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    parser.add_argument('--precision', default=3,                          help='Controls display of probabilities')
    args = parser.parse_args()
    if args.sample:
        Transitions,Emissions = EstimateParameters('yzzzyxzxxx',
                           ['x',  'y',   'z'],
                           'BBABABABAB',
                           ['A', 'B', 'C'])
        for row in formatTransition(Transitions,['A', 'B', 'C'],precision=args.precision):
            print (row)
        for row in formatEmission(Emissions,['A', 'B', 'C'], ['x',  'y',   'z'],precision=args.precision):
            print (row)
    

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = None
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
