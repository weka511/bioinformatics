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

# BIP Testing Bipartiteness

from helpers import format_list, create_adjacency, parse_graphs

# bip
#
# Input: a graph in edgelist form
#
# Output:  1 if praph is bipartite, otherwise -1
#
# A graph is bipartite if and only if it does not contain an odd cycle
# https://en.wikipedia.org/wiki/Bipartite_graph#Characterization

def bip(graph):

    # explore
    #
    # Inputs: node
    #         edges
    #
    # Outputs: True iff there is an odd cycle
    
    def explore(node,edges):
        for linked in adjacency[node]:
            if linked==node: continue
            
            if  (node,linked) in edges or (linked,node) in edges: continue
           
            if len(edges)>0 and linked==edges[0][0]:   # we have a cycle
                explored.add(node)
                edges     = edges + [(node,linked)]
                if len(edges)%2==1:                    # odd cycle
                    print (edges)
                    return True
            else:
                if linked in explored:                 # we've already seen this link
                    explored.add(node)
                    continue
                if explore(linked,edges+[(node,linked)]): # add this link and continue
                    return True
        explored.add(node)        
        return False
   
    m,_,adjacency = create_adjacency(graph,back=True)

    explored = set()
   
    for node in range(1,m+1):
        if  len(adjacency[node])>1:
            if explore(node,[]): return -1
 
    return 1




if __name__=='__main__':
    import sys
    simple = [[[3, 3],
               [1, 2],
               [3, 2],
               [3, 1]],
        
              [[4, 3],
               [1, 4],
               [3, 1],
               [1, 2]]]
    print (format_list([bip(g) for g in simple]))
    sys.setrecursionlimit(2000)
    with open(r'C:\Users\Simon\Downloads\rosalind_bip(10).txt') as f:
        print (format_list([bip(g) for g in parse_graphs(f)]))
 