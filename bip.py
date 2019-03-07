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
# Strategy: try to partition graph into two colours

def bip(graph):
    red  = set()
    blue = set()
    
    # colour
    #
    # Attempt to assign this node, and all reachable nodes, to one colour or t'other
    #
    # Inputs:  node    The node we are assigning
    #          isBlue  Indicates whether we are trying Red or Blue
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
    
    # Purge isolated nodes
    
    for k in [k for k,v in adjacency.items() if len(v)==0]:
        adjacency.pop(k)

    # Try to colour nodes
    
    for node in adjacency.keys():
        if node in red or node in blue: continue
        red.add(node)
        for link in adjacency[node]:
            coloured = colour(link,True)
            if not coloured:
                return -1
            
    return 1  # assume bipartite unless we fail

if __name__=='__main__':

    
    with open(r'C:\Users\Simon\Downloads\rosalind_bip(3).txt') as f:
        print (format_list([bip(g) for g in parse_graphs(f)]))
 