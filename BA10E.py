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
    # State
    #
    # This class and its children represent the states of the HMM
    
    class State:
        INSERT  = 0
        MATCH   = 1 
        DELETE  = 2
        END     = 3
        N_TYPES = 4
            
        def __init__(self,index=None):
            self.index              = index
            self.emissions         = {}
            self.next_state_counts = [0]*State.N_TYPES
            
        def record_transition(self,conserved,ch):
            next_state = None
            if conserved:
                if ch in Alphabet:
                    next_state,offset =  self.match(ch)
                else:
                    next_state,offset =  self.delete()
            else:
                if ch in Alphabet:
                    next_state,offset =  self.insert(ch)
                else:
                    next_state,offset = State.MATCH,1
            self.next_state_counts[next_state] += 1
            return next_state,self.index+offset
        
        def match(self,ch):
            return State.MATCH,1
        
        def insert(self,ch):
            self.match(ch)
            return State.INSERT,0
        
        def delete(self):
            return State.DELETE,1
        
        def record_emission(self,ch,Alphabet):
            if ch not in Alphabet: return
            if ch in self.emissions:
                self.emissions[ch] += 1
            else:
                self.emissions[ch] = 1 
                  
        
        def get_normalized_state_counts(self):
            total = sum(self.next_state_counts)
            return [c/total for c in self.next_state_counts] if total>0 else self.next_state_counts
        
        def get_emission(self,Alphabet):
            K           = len(Alphabet)
            Frequencies = [0]*K
            count_chars = sum([count for _,count in self.emissions.items()])
            if count_chars >0:
                for i in range(K):
                    if Alphabet[i] in self.emissions:
                        Frequencies[i] = self.emissions[Alphabet[i]]/count_chars
                #assert (sum(Frequencies)==1.0)
            return Frequencies
        
    class Start(State):
        def __init__(self):
            super().__init__(0)
            
        def __str__(self):
            return 'S'
        
        def get_transition(self,L,pos):
            Result                 = [0]*L
            transition_frequencies = self.get_normalized_state_counts()
            for i in range(3):
                Result[1+i]        = transition_frequencies[i]
            return Result
        
    class Match(State):
        def __init__(self,index):
            super().__init__(index)
            
        def __str__(self):
            return f'M{self.index}'
        
        def get_transition(self,L,pos):
            Result                 = [0]*L
            transition_frequencies = self.get_normalized_state_counts()
            for i in range(3):
                if pos+i+2<L:
                    Result[pos+i+2]    = transition_frequencies[i]
            return Result        
        
    class Insert(State):
        def __init__(self,index):
            super().__init__(index)
            
        def __str__(self):
            return f'I{self.index}' 
        
        def get_transition(self,L,pos):
            Result                 = [0]*L
            transition_frequencies = self.get_normalized_state_counts()
            for i in range(3):
                if pos+i<L:
                    Result[pos+i]        = transition_frequencies[i]
            return Result         
        
    class Delete(State):
        def __init__(self,index):
            super().__init__(index)
            
        def __str__(self):
            return f'D{self.index}'        
        
        def get_transition(self,L,pos):
            Result                 = [0]*L
            transition_frequencies = self.get_normalized_state_counts()
            for i in range(3):
                if pos+i+1<L:
                    Result[pos+i+1] = transition_frequencies[i]
            return Result 
        
    class End(State):
        def __init__(self):
            super().__init__(0)
            
        def __str__(self):
            return 'E'   
        
        def get_transition(self,L,pos):
            Result                 = [0.0]*L
            return Result 
        
#   CountChars

    def CountChars(m,j,K):
        RowCounts = [0]*K
        for i in range(m):
            if Alignment[i][j] in Alphabet:
                RowCounts[Alphabet.index(Alignment[i][j])]+=1
        return RowCounts 
    
# create_states
#
#   construct list of states

    def create_states(Conserved):
        index  = 0
        Product = [Start(),Insert(index)]
        for conserved in Conserved:
            if conserved:
                index += 1
                Product.append(Match(index))
                Product.append(Delete(index))
                Product.append(Insert(index))
        Product.append(End())
        return Product

#   get_state_index

    def get_state_index(state_type,index):
        if State.MATCH == state_type:
            return 3*index -1
        elif State.INSERT == state_type:
            return 3*index + 1 
        else:  # STATE.DELETE
            return 3*index

#   Accumulate statistics

    def accumulate_statistics(States,Alignment,n):
        for Sequence in Alignment:
            state_index = 0   # Index in states array
            for str_index in range(n):
                # next_state_type is State.MATCH(1), State.INSERT(0), State.DELETE (2), or State.END (3)
                # index tells is whether we are dealing with I0, I1/M1/D1, etc
                next_state_type,index = States[state_index].record_transition(Conserved[str_index],Sequence[str_index])
                state_index           = get_state_index(next_state_type,index)
                States[state_index].record_emission(Sequence[str_index],Alphabet)
                
    n              = len(Alignment[0])
    m              = len(Alignment)
    K              = len(Alphabet)
    
#   construct profile
    Counts         = [CountChars(m,j,K) for j in range(n)]
    Conserved      = [sum(Count) > (1-theta)*K for Count in Counts]
    
#   construct list of states

    States = create_states(Conserved)
    L      = len(States)
    
    accumulate_statistics(States,Alignment,n)    
 
    Transition = []
    Emission   = []

    for i in range(L):
        Transition.append(States[i].get_transition(L,i))
        Emission.append(States[i].get_emission(Alphabet))
        
    return States,Transition,Emission
    
def formatEmission(Emission,States,Alphabet):  
    yield '\t' + '\t'.join(Alphabet)
    for row,state in zip(Emission,States):
        yield str(state) + '\t' + '\t'.join(f'{r:.3f}' for r in row)

def formatTransition(Transition,States):
    yield '\t' + '\t'.join(str(s) for s in States)
    for row,state in zip(Transition,States):
        yield str(state) + '\t' + '\t'.join(f'{r:.3f}' for r in row)
        
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
                f.write(f'{row}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
