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

def cte(edges):
    # relabel
    #
    # Exploit Dijkstra's algorithm by relabelling nodes so that first edge starts at 1
    def relabel(edges):
        def subst(x):
            if   x==1: return b
            elif x==b: return 1
            else:
                return x
        def transform(x,y,w):
            return (subst(x),subst(y),w)


        return [(n,k-1)] + [transform (x,y,w) for (x,y,w) in edges[2:]]
    
    n,k         = edges[0]
    a,b,w       = edges[1]        
    new_edges   = relabel(edges)
    path_length = dij(new_edges)[a-1]
    return path_length + w if path_length>-1 else -1

def ddeg(n,M,A):
    lookup=deg(n,m,A)
    sums=[0 for a in range(n)]
    for (a,b) in A:
        sums[a-1]+=lookup[b-1]
        sums[b-1]+=lookup[a-1]    
    return sums

def deg(n,m,A):
    degrees=[0 for a in range(n)]
    for (a,b) in A:
        print (a,b)
        degrees[a-1]+=1
        degrees[b-1]+=1
    return degrees


# DIJ  Dijkstra's Algorithm: compute single-source shortest distances 
#                            in a directed graph with positive edge weights.
def dij(g):
    n,_             = g[0]
    open_edges      = [edge for edge in g[1:]]
    D               = [-1 for i in range(n+1)]
    D[1]            = 0
    for i in range(n):
        for a,b,w in open_edges:
            if D[a]>-1:
                proposed_distance = D[a]+w
                if D[b]==-1 or D[b]>proposed_distance:
                    D[b]=proposed_distance
                
    return D[1:]
