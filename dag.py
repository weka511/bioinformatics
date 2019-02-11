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

from helpers import format_list, create_adjacency
from bip import parse_graphs

def explore(node,adjacency,explored,path=[]):
    explored.add(node)
    linked = adjacency[node]
    for succ in linked:
        if succ == node: continue
        if succ in path: return True
        if explore(succ,adjacency,explored,path+[node]): return True
    return False

def dag(graph):
    m,_,adjacency = create_adjacency(graph,False)
    explored   = set()
    
    for a,_ in adjacency.items():
        if not a in explored:
            if explore(a,adjacency,explored): return -1
            
    return 1

if __name__=='__main__':
    simple = [
        [[2, 1],
         [1, 2]],
        
        [[4, 4],
         [4, 1],
         [1, 2],
         [2, 3],
         [3, 1]],
        
        [[4, 3],
         [4, 3],
         [3, 2],
         [2, 1]]]
    print (format_list([dag(g) for g in simple]))
    
    with open(r'C:\Users\Simon\Downloads\rosalind_dag.txt') as f:
        print (format_list([dag(g) for g in parse_graphs(f)]))    