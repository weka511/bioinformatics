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

def ConstructProfileHMM(theta,Alphabet,Alignment,trace=True):
    # State
    #
    # This class and its children represent the states of the HMM
    
    class State:
        INSERT   = 0
        MATCH    = INSERT + 1 
        DELETE   = MATCH  + 1
        N_STATES = DELETE + 1        # Number of states

        def __init__(self,index=None):
            self.index             = index
            self.emissions         = {}
            self.next_state_counts = [0]*State.N_STATES

#       record_transition
#
#       Used to record the transition resulting from each character
#       Parameters:
#            ch         The character
#            conserved  Whether or not this is a conserved position
        @abc.abstractmethod
        def record_transition(self,ch,conserved):
            pass 
        
        def get_match_plus_offset(self,ch):
            return State.MATCH,1
        
        def get_insert_plus_offset(self,ch):
            return State.INSERT,0
        
        def get_delete_plus_offset(self):
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
        
        def record_transition(self,ch,conserved):
            next_state = None
            if conserved:
                if ch in Alphabet:
                    next_state,offset = self.get_match_plus_offset(ch)
                else:
                    next_state,offset = self.get_delete_plus_offset()
            else:
                if ch in Alphabet:
                    next_state,offset = self.get_insert_plus_offset(ch)
                else:
                    next_state,offset = self.get_match_plus_offset(ch)
            self.next_state_counts[next_state] += 1
            return next_state,self.index+offset
        
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
        
        def record_transition(self,ch,conserved):
            next_state = None
            if ch=='$':
                next_state,offset = State.MATCH,1
                self.next_state_counts[next_state] += 1
                return next_state,self.index+offset
            if conserved:
                if ch in Alphabet:
                    next_state,offset = self.get_match_plus_offset(ch)
                else:
                    next_state,offset = self.get_delete_plus_offset()
            else:
                if ch in Alphabet:
                    next_state,offset = self.get_insert_plus_offset(ch)
                else:
                    next_state,offset = State.MATCH,0
            if next_state!=State.MATCH or offset!=0:
                self.next_state_counts[next_state] += 1
            return next_state,self.index+offset        
        
    class Insert(State):
        def __init__(self,index):
            super().__init__(index)
            
        def __str__(self):
            return f'I{self.index}' 
 
        def record_transition(self,ch,conserved):
            if ch=='$':
                next_state,offset = State.MATCH,1
                self.next_state_counts[next_state] += 1
                return next_state,self.index+offset            
            next_state = None
            if ch=='-':
                next_state = State.INSERT
                offset     = 1 if conserved else 0         
            else:
                if conserved:
                    next_state,offset = self.get_match_plus_offset(ch)
                else:
                    next_state,offset = self.get_insert_plus_offset(ch)
            if ch!='-':#next_state!=State.INSERT or offset!=0:
                self.next_state_counts[next_state] += 1
            return next_state,self.index+offset 
        
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
        
        def record_transition(self,ch,conserved):
            if ch=='$':
                next_state,offset = State.MATCH,1
                self.next_state_counts[next_state] += 1
                return next_state,self.index+offset            
            next_state = None
            if ch=='-':
                next_state = State.DELETE
                offset     = 1 if conserved else 0                   
            else:
                if conserved:
                    next_state,offset = self.get_match_plus_offset(ch)
                else:
                    next_state,offset = self.get_insert_plus_offset(ch)
            if next_state!=State.DELETE or offset!=0:
                self.next_state_counts[next_state] += 1
            return next_state,self.index+offset        
            
    class End(State):
        def __init__(self):
            super().__init__(0)
            
        def __str__(self):
            return 'E'   
        
        def get_transition(self,L,pos):
            Result                 = [0.0]*L
            return Result 
        
        def record_transition(self,ch,conserved):
            raise RosalindException('Should never get here')
        
#   Tracer
#
#   Used to map flow through states

    class Tracer:
        def __init__(self,trace,m,theta):
            self.StateTrace = [[(-1,0,0,'^')] for _ in range(m)] if trace else None
            self.trace      = trace
            self.m          = m
            if trace:
                print (f'Theta = {theta}. There are {len(States)} States')
                
        def trace_state(self,next_state_type,index,state_index,ch,i):
            if self.trace:
                self.StateTrace[i].append((next_state_type,index,state_index,ch)) 
        
        def display(self):
            if not trace: return
            for i in range(self.m):
                print ('-'.join(self.format_trace(rec) for rec in self.StateTrace[i]))
                    
        def format_trace(self,trace_record):
            name = ['S', 'I', 'M', 'D', 'E'][trace_record[0]+1]
            return f'{name}{trace_record[1]}({trace_record[2]},{trace_record[3]})' 
        
        def trace_boxen(self,States):
            if not trace: return
            for state in States:
                print (state, ' '.join(str(c) for c in state.next_state_counts))
                
        def trace_exception(self,e,i,Sequence,str_index,State):
            print (f'Exception {e} row {i} position {str_index}')
            print (Sequence)
            print (State)
                
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
    
# create_states
#
# Construct list of states.
# Parameters:
#    Conserved         List indicating whether or not count of symbols in each position exceeds threshold.
#
# Returns: list of states [S, I0, M1,D1,I1, ..., E], with one MDI group for each conserved position

    def create_states(Conserved):
        conserved_index  = 0
        Product = [Start(),Insert(conserved_index)]
        for conserved in Conserved:
            if conserved:
                conserved_index += 1
                Product.append(Match(conserved_index))
                Product.append(Delete(conserved_index))
                Product.append(Insert(conserved_index))
                   
        Product.append(End())
        return Product

#   get_state_index
#
#   Used to locate states in State array.
#
#   Parameters:
#       state_type           M, D or I
#       conserved_index      Indicates which group MDI belongs to 
#
#   Returns:   index in States array

    def get_state_index(state_type,conserved_index):
        def get_offset():
            if State.MATCH == state_type:
                return -1
            elif State.INSERT == state_type:
                return  +1 
            elif State.DELETE == state_type:
                return 0
            else:
                raise RosalindException(f'Could not get state index {state_type} {conserved_index}')
        return   3*conserved_index  + get_offset()

    
#   Accumulate statistics

    def accumulate_statistics(States,Alignment,n,trace=False):
    
        tracer = Tracer(trace,m,theta)
        
        for i in range(m):
            Sequence    = Alignment[i]
            state_index = 0   # Index in states array
            for str_index in range(n):
                # next_state_type is State.MATCH(1), State.INSERT(0), or State.DELETE (2)
                # index tells is whether we are dealing with I0, I1/M1/D1, etc
                try:
                    next_state_type,index = States[state_index].record_transition(Sequence[str_index],Conserved[str_index])
                except RosalindException as e:
                    tracer.trace_exception(e,i,Sequence,str_index,States[state_index])
                state_index           = min(get_state_index(next_state_type,index),len(States)-1) # FIXME
                States[state_index].record_emission(Sequence[str_index],Alphabet)
                tracer.trace_state(next_state_type,index,state_index,Sequence[str_index],i)
            try:
                States[state_index].record_transition('$',None)
            except RosalindException as e:
                tracer.trace_exception(e,i,Sequence,str_index,States[state_index])
        tracer.display()
        tracer.trace_boxen(States)
 

#   Useful constants - lengths of arrays

    K              = len(Alphabet)
    m              = len(Alignment)
    n              = len(Alignment[0])   
    for Sequence in Alignment[1:]:      # All sequences should be the same length
        assert(n == len(Sequence))
    
#   construct profile
    Counts         = [CountChars(m,j,K) for j in range(n)]
    
#   Indicate whether or not symbols in column are over threshold
    Conserved      = [sum(Count) > (1-theta)*K for Count in Counts]
    
#   construct list of states

    States = create_states(Conserved)

    L      = len(States)
    
    accumulate_statistics(States,Alignment,n,trace=trace)    
 
    Transition = []
    Emission   = []

    for i in range(L):
        Transition.append(States[i].get_transition(L,i))
        Emission.append(States[i].get_emission(Alphabet))
        
    return States,Transition,Emission

# float2str
#
# Format a floatinmg point number for display
#
# Parameters:
#     x         Value to be displayed
#     procision Number of digits (after decimal point)

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

def formatEmission(Emission,States,Alphabet,precision=2):  
    yield '\t' + '\t'.join(Alphabet)
    for row,state in zip(Emission,States):
        yield str(state) + '\t' + '\t'.join(float2str(r,precision) for r in row)

def formatTransition(Transition,States,precision=2):
    yield '\t' + '\t'.join(str(s) for s in States)
    for row,state in zip(Transition,States):
        yield str(state) + '\t' + '\t'.join(float2str(r,precision) for r in row)
        
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA10E Construct a Profile HMM')
    parser.add_argument('--sample',    default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',     default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind',  default=False, action='store_true', help='process Rosalind dataset')
    parser.add_argument('--text',      default=False, action='store_true', help='process dataset from textbook')
    parser.add_argument('--precision', default=3,                          help='Controls display of probabilities')
    parser.add_argument('--trace',     default=False, action='store_true', help='Trace progression through states')
    args = parser.parse_args()
    if args.sample:
        States,Transition,Emission = ConstructProfileHMM(0.289,
                                                         ['A',   'B',   'C',   'D',   'E'],
                                                         ['EBA', 'EBD', 'EB-', 'EED', 'EBD', 'EBE','E-D','EBD'],
                                                         trace=args.trace)
        
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
                                                          'ADDEFAAADF'],
                                                         trace=args.trace)
        
        for row in formatTransition(Transition,States,precision=args.precision):
            print (row)
        print ('--------')
        for row in formatEmission(Emission,States,Alphabet,precision=args.precision):
            print (row)
            
            
    if args.extra:
        Input,Expected             = read_strings(f'data/ProfileHMM.txt',init=0)
        States,Transition,Emission = ConstructProfileHMM(float(Input[0]),
                                                         Input[2].split(),
                                                         Input[4:-1],
                                                         trace=args.trace) 
        for row in formatTransition(Transition,States,precision=args.precision):
            print (row)
        print ('--------')
        for row in formatEmission(Emission,States,Input[2].split(),precision=args.precision):
            print (row)
               
    if args.rosalind:
        Input                      = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        States,Transition,Emission = ConstructProfileHMM(float(Input[0]),
                                                         Input[2].split(),
                                                         Input[4:-1],
                                                         trace=args.trace) 

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
