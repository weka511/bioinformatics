#    Copyright (C) 2019 Greenweaves Software Limited
#
#    This is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This software is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

# HDAG Hamiltonian Path in DAG 
#

from helpers import create_adjacency
from align import topological_order

def hdag(graph):  
    m,_,adj = create_adjacency(graph,back=False,self=False)
    t = topological_order(adj)
    m,_,adj = create_adjacency(graph,back=False,self=False)
    for a,b in zip(t[:-1],t[1:]):
        if not b in adj[a]:
            return (-1,[])
    return (1,t)
        

if __name__=='__main__':
    from helpers import create_strings
       
    gs = []
    g  = []
    for line in create_strings(ext=1):
        if len(line)==0:
            if len(g)>0:
                gs.append(g)
            g=[]
            continue
        numbers = [int(s) for s in line.split(' ')]
        if len(numbers)==1:
            continue
        elif len(numbers)==2:
            
            g.append((numbers[0],numbers[1]))

    if len(g)>0:
        gs.append(g)           
 
    
    for g in gs:
        result,path=hdag(g)
        if result==-1:
            print (result)
        else:
            print (' '.join(str(x) for x in [result]+path) )
                
 