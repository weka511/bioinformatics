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

from rosalind import LabelledTree

def SmallParsimony(T):
    def SmallParsimonyC(Character,alphabet='ATGC'):
        def get_ripe():
            for v in T.get_nodes():
                if not processed[v] and v in T.edges:
                    for e,w in T.edges[v]:
                        if e>v: continue
                        if not processed[e]: break
      
                    return v
            return None
        
        def calculate_s(symbol,v,alphabet):
            def get_min(e):
                minimum=float('inf')
                for i in range(len(alphabet)):
                    delta=0 if symbol==alphabet[i] else 1
                    candidate=s[e][i]+delta
                    if candidate<minimum:
                        minimum=candidate
                return minimum
            return sum([get_min(e) for e,_ in T.edges[v]])
        
        print (Character)
        processed={}
        s={}
        for v in T.get_nodes():      
            if T.is_leaf(v):
                processed[v]=True
                values={}
                s[v]=[0 if symbol==Character[v] else float('inf') for symbol in alphabet]
                print(v,s[v])
            else:
                processed[v]=False
                 
        v = get_ripe()
        while not v == None:
            processed[v]=True
            s[v]=[calculate_s(symbol,v,alphabet) for symbol in alphabet ]
            print (s[v])
            v = get_ripe()
        return sum([min(s[v]) for v in T.get_nodes()])
    
    m=len(T.labels[0])
    score=0
    for i in range(m):
        score+=SmallParsimonyC([v[i] for l,v in T.labels.items()])
    return score



if __name__=='__main__':
    N=4
    T=LabelledTree.parse(4,
                  ['4->CAAATCCC',
                   '4->ATTGCGAC',
                   '5->CTGCGCTG',
                   '5->ATGGACGA',
                   '6->4',
                   '6->5']  )
    T.print()
    
    s=SmallParsimony(T)
    print (s)
  
