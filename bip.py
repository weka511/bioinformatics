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
from collections import deque

# bip
#
# Input: a graph in edgelist form
#
# Output:  1 if graph is bipartite, otherwise -1
#

def bip(graph):
    red  = set()
    blue = set()
    
    # colour
    #
    # Attempt to assign this node, and all reachable nodes, to one colour or t'other
    def colour(node,isBlue):
        if isBlue:
            if node in red: return False
            if node in blue: return True
            blue.add(node)
        else:
            if node in blue: return False
            if node in red: return True
            red.add(node)            
        for link in adjacency[node]:
            if not colour(link,not isBlue): return False
        return True
            
    _,_,adjacency = create_adjacency(graph,back=True,self=False)
    for k in [k for k,v in adjacency.items() if len(v)==0]:
        adjacency.pop(k)
    #for k,v in adjacency.items():
        #print (k,v)

    for node in adjacency.keys():
        if node in red or node in blue: continue
        red.add(node)
        for link in adjacency[node]:
            coloured = colour(link,True)
            if not coloured:
                return -1
            
    return 1  # assume bipartite, i.e. no odd cycles

if __name__=='__main__':

    simple = [[[3, 3],
               [1, 2],
               [3, 2],
               [3, 1]],
        
              [[4, 3],
               [1, 4],
               [3, 1],
               [1, 2]]]
    

    #print (format_list([bip(g) for g in simple]))
    
    with open(r'C:\Users\Simon\Downloads\rosalind_bip(3).txt') as f:
        #gs = parse_graphs(f)
        #bip(gs[0])
        print (format_list([bip(g) for g in parse_graphs(f)]))
 