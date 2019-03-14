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

# CTE Shortest Cycle Through a Given Edge
#
# Returns: For each graph, output the length of a shortest cycle going
#           through the first specified edge if there is a cycle and "-1" otherwise.

from dij import dij

def cte(edges):
    def relabel(edges):
        def subst(x):
            if   x==1: return b
            elif x==b: return 1
            else:
                return x
        def transform(x,y,w):
            return (subst(x),subst(y),w)


        return [(n,k-1)] + [transform (x,y,w) for (x,y,w) in edges[2:]]
    
    n,k   = edges[0]
    a,b,w = edges[1]        
    new_edges = relabel(edges)
    #print (new_edges)
    path_length = dij(new_edges)[a-1]
    return path_length + w if path_length>-1 else -1

if __name__=='__main__':
    from helpers import create_strings

    #k = 2
    #gs = [
          ##[(4, 5),
           ##(2, 4, 2),
           ##(3, 2, 1),
           ##(1, 4, 3),
           ##(2, 1, 10),
           ##(1, 3, 4)],

          #[(4, 5),
           #(3, 2, 1),
           #(2, 4, 2),
           #(4, 1, 3),
           #(2, 1, 10),
           #(1, 3, 4)]]
    #print (' '.join([str(cte(g)) for g in gs]))   
    gs = []
    g  = []
    for line in create_strings(ext=1):
        numbers = [int(s) for s in line.split(' ')]
        if len(numbers)==1:
            pass
        elif len(numbers)==2:
            if len(g)>0:
                gs.append(g)
            g = [(numbers[0],numbers[1])]
        else:
            g.append((numbers[0],numbers[1],numbers[2]))
    if len(g)>0:
        gs.append(g)           
 
        
    print (' '.join([str(cte(g)) for g in gs])) 
  