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

def SmallParsimony(T, Character,alphabet='ATGC'):
    def get_ripe():
        for v in T.get_nodes():
            if not processed(v):
                for e,w in T.edges[v]:
                    if e>v: continue
                    if not processed(e): break
                return v
        return None
    
    processed={}
    s={}
    for v in T.get_nodes():      
        if T.is_leaf(v):
            processed[v]=True
            values={}
            for k in alphabet:
                values[k]=0 if k==Character else float('inf')
                s[v]=values
            else:
                processed[v]=false
    v = get_ripe()
    while not v == Node:
        processed[v]=True
        #do stuff
        v = get_ripe()
    

if __name__=='__main__':
    N=4
    T=LabelledTree.parse(4,
                  ['4->CAAATCCC',
                   '4->ATTGCGAC',
                   '5->CTGCGCTG',
                   '5->ATGGACGA',
                   '6->4',
                   '6->5']  )
    SmallParsimony(T, '',alphabet='ATGC')
  
