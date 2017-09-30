# Copyright (C) 2017 Greenweaves Software Pty Ltd

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

# BA7F Implement SmallParsimony 

import operator,random

from rosalind import LabelledTree,hamm

def SmallParsimony(T,alphabet='ATGC'):

    def SmallParsimonyC(Character):
        def get_ripe():
            for v in T.get_nodes():
                if not processed[v] and v in T.edges:
                    for e,_ in T.edges[v]:
                        if e>v: continue
                        if not processed[e]: break
      
                    return v
            return None
        
        def calculate_s(symbol,v):
            '''
             Calculate score if node v is set to a specified symbol
                Parameters:
                    symbol The symbol, e.g. 'A', not the index in alphabet
                    v      The node
             '''
            def delta(i):
                '''
                Complement of Kronecker delta
                '''
                return 0 if symbol==alphabet[i] else 1 
            def get_min(e):
                return min(s[e][i]+delta(i) for i in range(len(alphabet)))
            
            return sum([get_min(e) for e,_ in T.edges[v]])
        
        def update_assignments(v,s):
            if not v in assignments.labels:
                assignments.labels[v]=''
            index=0
            min_s=float('inf')
            for i in range(len(s)):
                if s[i]<min_s:
                    min_s=s[i]
                    index=i
            assignments.set_label(v,assignments.labels[v]+alphabet[index])
            return alphabet[index]
        
        def backtrack(v, current_assignment):
            for v_next,_ in T.edges[v]:
                if T.is_leaf(v_next): continue
                if not v_next in assignments.labels:
                    assignments.labels[v_next]=''
                min_score=min([s[v_next][i] for i in range(len(alphabet))])
                indices=[i for i in range(len(alphabet)) if s[v_next][i]==min_score ]
                matched=False
                for i in indices:
                    if alphabet[i]==current_assignment:
                       
                        matched=True
                        assignments.set_label(v_next,assignments.labels[v_next]+current_assignment)
                        backtrack(v_next,current_assignment)
                if not matched:
                    # Black magic alert: I am not clear why the introduction of random numbers
                    # helps here. Maybe it stops the tree being biased towatds the first strings
                    # in the alphabet.
                    next_assignment=alphabet[indices[random.randrange(0,(len(indices)))]]               
                    assignments.set_label(v_next,assignments.labels[v_next]+next_assignment)
                    backtrack(v_next,next_assignment)        
            
        processed={}
        s={}
        for v in T.get_nodes():      
            if T.is_leaf(v):
                processed[v]=True
                s[v]=[0 if symbol==Character[v] else float('inf') for symbol in alphabet]
                #print(v,s[v])
            else:
                processed[v]=False
                 
        v = get_ripe()
        while not v == None:
            processed[v]=True
            s[v]=[calculate_s(symbol,v) for symbol in alphabet ]
            v_last=v
            v = get_ripe()

        backtrack(v_last,update_assignments(v_last,s[v_last]))
        return min([s[v_last][c] for c in range(len(alphabet))])
       
    assignments=LabelledTree(T.N)
    assignments.initialize_from(T)
    
    return sum([SmallParsimonyC([v[i] for l,v in T.labels.items()]) for i in range(len(T.labels[0]))]),assignments



if __name__=='__main__':
       
    #N=4
    #T=LabelledTree.parse(4,
                  #['4->CAAATCCC',
                   #'4->ATTGCGAC',
                   #'5->CTGCGCTG',
                   #'5->ATGGACGA',
                   #'6->4',
                   #'6->5'],
                  #bidirectional=False
                  #)
    #T.print()
    N=-1
    rest=[]
    with open('c:/Users/Weka/Downloads/rosalind_ba7f(5).txt') as f:
        for line in f:
            if N==-1:
                N=int(line.strip())
            else:
                rest.append(line.strip())
    T=LabelledTree.parse(N,rest,bidirectional=False)            
    
    score,assignments=SmallParsimony(T)
    print (score)
    assignments.nodes.sort()
    #for l,v in assignments.labels.items():
        #print (l,v)
    for node in assignments.nodes:
        if node in assignments.edges:
            for edge in assignments.edges[node]:
                end,weight=edge
                if node in assignments.labels and end in assignments.labels:
                    print ('{0}->{1}:{2}'.format(assignments.labels[node],
                                                 assignments.labels[end],
                                                 hamm(assignments.labels[node],assignments.labels[end])))    
 
  
