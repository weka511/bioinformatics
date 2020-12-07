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
            self.next_state_counts[next_state] += 1
            return next_state,self.index+offset
        
        def match(self,ch):
            return State.MATCH,1
        
        def insert(self,ch):
            self.match(ch)
            return State.INSERT,0
        
        def delete(self):
            return State.DELETE,1
        
        def record_emission(self,ch):
            if ch in self.emissions:
                self.emissions[ch] += 1
            else:
                self.emissions[ch] = 1 
                
        def display(self):
            print (self.emissions)
            print (self.next_state_counts)
            print ('----------')        
        
    class Start(State):
        def __init__(self):
            super().__init__(0)
            
        def __str__(self):
            return 'S'
        
    class Match(State):
        def __init__(self,index):
            super().__init__(index)
            
        def __str__(self):
            return f'M{self.index}'
        
    class Insert(State):
        def __init__(self,index):
            super().__init__(index)
            
        def __str__(self):
            return f'I{self.index}' 
        
    class Delete(State):
        def __init__(self,index):
            super().__init__(index)
            
        def __str__(self):
            return f'D{self.index}'        
            
    class End(State):
        def __init__(self):
            super().__init__()
            
        def __str__(self):
            return 'E'        
            
    def CountChars(m,j,K):
        RowCounts = [0]*K
        for i in range(m):
            if Alignment[i][j] in Alphabet:
                RowCounts[Alphabet.index(Alignment[i][j])]+=1
        return RowCounts  
    
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
 
    def get_state_index(state_type,index):
        if State.MATCH == state_type:
            return 3*index -1
        elif State.INSERT == state_type:
            return 3*index + 1 
        else:  # STATE.DELETE
            return 3*index

#   Accumulate statistics

    def accumulate_statistics(States,Alignment,n):
        for Sequence in Alignment[0:1]:
            state_index = 0   # Index in states array
            for str_index in range(n):
                # next_state_type is State.MATCH(1), State.INSERT(0), State.DELETE (2), or State.END (3)
                # index tells is whether we are dealing with I0, I1/M1/D1, etc
                next_state_type,index = States[state_index].record_transition(Conserved[str_index],Sequence[str_index])
                state_index           = get_state_index(next_state_type,index)
                States[state_index].record_emission(Sequence[str_index])
                
    n              = len(Alignment[0])
    m              = len(Alignment)
    K              = len(Alphabet)
    
#   construct profile
    Counts         = [CountChars(m,j,K) for j in range(n)]
    Conserved      = [sum(Count) > (1-theta)*K for Count in Counts]
    
#   construct list of states

    States = create_states(Conserved)
    
    accumulate_statistics(States,Alignment,n)

    for state in States:
        print (state)
        print (state.display())
    
 
    Transition = []
    Emission   = []


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
        #for row in formatTransition(Transition,States):
            #print (row)
        #print ('--------')
        #for row in formatEmission(Emission,States,Input[2].split()):
            #print (row)
               
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


    #def CountChars(m,j,K):
        #RowCounts = [0]*K
        #for i in range(m):
            #if Alignment[i][j] in Alphabet:
                #RowCounts[Alphabet.index(Alignment[i][j])]+=1
        #return RowCounts
    
    #def CreateStates(m):
        #Product = [('S',0),('I',0)]
        #max_pos = -1
        #for i in range(m):
            #Product.append(('M',i+1))
            #Product.append(('D',i+1))
            #Product.append(('I',i+1))
            #max_pos = i+1
        #Product.append(('E',None))
        #return Product,max_pos
    
    #def StateGenerator(States):
        #for l in range(L): # iterate through states
            #operation,pos = States[l]
            #yield operation,pos
   
    #def get_index_of_state(operation,pos):
        #for i in range(len(States)):
            #if (operation,pos)==States[i]:
                #return i
  
    #def create_actions(Conserved,Incomplete):
        #Product = ['S']
        #for conserved,incomplete in zip(Conserved,Incomplete):
            #if conserved:
                #Product.append('M')
                #if incomplete:
                    #Product.append('D')
            #else:
                #Product.append('I')
        #Product.append('E')
        #return Product
    
        
    #n              = len(Alignment[0])
    #m              = len(Alignment)
    #K              = len(Alphabet)
    #Counts         = [CountChars(m,j,K) for j in range(n)]
    #Conserved      = [sum(Count) > (1-theta)*K for Count in Counts]
    #Incomplete     = [sum(Count) < K for Count in Counts]
    #Actions        = create_actions(Conserved,Incomplete)
    
    #LinkBack       = [j for j in range(n) if Conserved[j]]
    #Profile        = [[c/sum(Counts[i]) for c in Counts[i]] for i in range(len(Counts)) if  Conserved[i]]
    #States,max_pos = CreateStates(len(Profile))
    #L              = len(States)
    #Transition     = [] # row--state HMM is in, column = P(move to)
    #Emission       = [] # row--state HMM is in, column = P emit this symbol
    #action         = Actions.pop(0)
    #next_action    = Actions[0]
    #state_seq      = -1
    #row_seq        = -1
    #for operation,pos in StateGenerator(States):
        #Transition.append([0]*L)
        #Emission.append([0]*K)
        #if operation!=action: continue     
        #if operation == 'S':
            #state_seq +=1
            #row_seq   += 1
            #if next_action=='I':
                #Transition[-1][state_seq+2] =sum(Counts[row_seq])/K
                #Transition[-1][state_seq+1] =1.0-Transition[-1][state_seq+1]
            #else:
                #Transition[-1][state_seq+2] =1.0
        #elif operation == 'M':
            #state_seq +=1
            #normalizerEmission = sum(Counts[row_seq])
            #for k in range(K):
                #Emission[-1][k] = Counts[row_seq][k]/normalizerEmission
            #row_seq   += 1
            #if next_action=='I':
                #pass
            #else:
                #Transition[-1][state_seq+3] = sum(Counts[row_seq])/len(Counts[row_seq])
                #Transition[-1][state_seq+4] = (len(Counts[row_seq])-sum(Counts[row_seq]))/len(Counts[row_seq])
        #elif operation == 'D':
            #state_seq +=1
            #row_seq   += 1
            #if next_action=='I':
                #pass
            #else:
                #Transition[-1][state_seq+3] = sum(Counts[row_seq])/len(Counts[row_seq])
                #Transition[-1][state_seq+4] = (len(Counts[row_seq])-sum(Counts[row_seq]))/len(Counts[row_seq])            
        #elif operation == 'I':
            #state_seq +=1
            #row_seq   += 1            
        #else:
            #state_seq +=1
            #row_seq   += 1
        #action = Actions.pop(0)
        #next_action = Actions[0]
        
    #seq            = 0
    #previous_del   = False
    #for operation,pos in StateGenerator(States):
        #if operation == 'S':
            #Transition.append([1 if op=='M' and pos2==1 else 0 for op,pos2 in StateGenerator(States)] ) # FIXME - I0
            #Emission.append([0]*K)   
        #elif operation == 'M':
            #next_row = [0]*L
            #if pos<len(LinkBack):   #FIXME
                #emission_count = sum(Counts[LinkBack[pos]])
                #if emission_count>0:
                    #next_row[seq+3] = emission_count/K
                #skip_count     = K - emission_count
                #if skip_count>0:
                    #next_row[seq+4] = skip_count/K
                    #previous_del   = True
            #Transition.append(next_row) 
            #Emission.append(Profile[pos-1])
        #elif operation == 'D':
            #next_row = [0]*L
            #if previous_del:
                #next_row[seq+2] = 1
                #previous_del   = False
            #Transition.append(next_row)
            #Emission.append([0]*K)
        #elif operation == 'I':
            #print (pos,seq,Conserved[seq])
            #Transition.append([0]*L)       # TODO
            #Emission.append([0]*K)       # TODO        
        #else:
            #Transition.append([0]*L)
            #Emission.append([0]*K)
        #seq += 1