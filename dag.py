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

#    dag Testing Acyclicity

from helpers import format_list, create_adjacency, parse_graphs

# dag  Testing Acyclicity
#
# Input: a graph
#
# Output: 1 if graph acyclic, otherwise -1

def dag(graph):
    
    # explore: find cyclic path
    #
    # Input: node
    #        path
    #
    # Output: True iff there is a cyclic path starting at node, and
    #         whole 1st elements coincide with path
    
    def explore(node,path=[]):
        explored.add(node)
        linked = adjacency[node]
        for succ in linked:
            if succ == node: continue
            if succ in path: return True
            if explore(succ,path+[node]): return True
            
        return False

    m,_,adjacency = create_adjacency(graph,False)
    explored   = set()
    
    for a,_ in adjacency.items():
        if not a in explored:
            if explore(a): return -1
            
    return 1

if __name__=='__main__':  
    with open(r'C:\Users\Simon\Downloads\rosalind_dag(1).txt') as f:
        print (format_list([dag(g) for g in parse_graphs(f)]))    