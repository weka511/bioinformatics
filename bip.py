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

# A graph is bipartite if and only if it does not contain an odd cycle
# https://en.wikipedia.org/wiki/Bipartite_graph#Characterization

import cc

def bip(graph):
    print ('bip')
    def explore(node,adjacency,path,explored):
        #if node in explored: return (False,[])
        print ('Explore',node,path)
        if node in path:
            if len(path)>=2 and len(path)%2 == 1 and node==path[0]:
                print ('Odd cycle',path+[node])
                return (True,path+[node])
        else:
            explored.add(node)
            for next_node in adjacency[node]:
                (found_odd,odd_cycle) = explore(next_node,adjacency,path+[node],explored)
                if found_odd: return (found_odd,odd_cycle)
 
        return (False,[])
        
    m,_,adjacency = cc.create_adjacency(graph)
    cycles = {}
    explored = set()
    for node in range(1,m+1):
        found_odd,odd_cycle = explore(node,adjacency,[],explored)
        if found_odd: return -1
    return 1
    
if __name__=='__main__':
    graphs = ([
            [[3, 3],
             [1, 2],
             [3, 2],
             [3, 1]],
        
            [[4, 3],
             [1, 4],
             [3, 1],
             [1, 2]]])
    print ([bip(g) for g in graphs])