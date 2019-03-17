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

from helpers import create_adjacency
from align import topological_order
from collections import deque

# bf
#
# Bellman-Ford Algorithm
#
# Given: A simple directed graph with integer edge weights and fewer than 1,000 vertices in the edge list format.
#
# Return: An array D[1..n] where D[i] is the length of a shortest path from the vertex 1 to the vertex i (D[1]=0).
#         If i is not reachable from 1 set D[i] to x

def bf(edges,s=1):

    n,_         = edges[0]
    dist        = [float('inf') for i in range(1,n+1)]
    previous    = [None for i in range(1,n+1)]

    dist[s-1]=0
    for i in range(n-1):
        for u,v,w in edges[1:]:
            dist[v-1] = min(dist[v-1],dist[u-1]+w)

    return [d if d<float('inf') else 'x' for d in dist]

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

# cc
#
# Connected Components
#
# Input: A simple graph with fewere than 1,000 vertices in the edge list format.
#
# Return: The number of connected components in the graph.

def cc(graph): 
    def explore(a,adjacency,explored,component):
        explored.add(a)
        component.append(a)
        for b in adjacency[a]:
            if not b in explored:
                explore(b,adjacency,explored,component)
        return sorted(list(set(component)))    
    m,_,adjacency = create_adjacency(graph)
  
    explored   = set()

    components = {}  # The connected components, with one element as key
    for a,_ in adjacency.items():
        if not a in explored:
            component = explore(a,adjacency,explored,[])
            for c in component:
                components[c]=component
                    
    count = 0
    uniques = []
    duplicates = []
    for k,v in components.items():
        if k in uniques:
            duplicates.append(k)
        else:
            for vv in v:
                uniques.append(vv)
            count+=1
    for d in duplicates:
        del components[d]

#   a few sanity checks

    nodes = sorted(list(set([v for k,vs in components.items() for v in vs])))
    assert len(nodes) ==m, '{0} not {1}'.format(len(nodes), m) 
    v0=0
    for v in nodes:
        assert v == v0+1
        v0=v
        
    return count,components

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

# ddeg
#
# Double Degree Array
#
# Input : A simple graph with fewer than 1,000 vertices in the edge list format.

# Return: An array D[1..n] where D[i] is the sum of the degrees of i's neighbors.

def ddeg(n,M,A):
    lookup=deg(n,m,A)
    sums=[0 for a in range(n)]
    for (a,b) in A:
        sums[a-1]+=lookup[b-1]
        sums[b-1]+=lookup[a-1]    
    return sums

# deg
#
# Degree array
#
# Input: A simple graph with fewer than 1,000 vertices in the edge list format.
#
# Return: An array D[1..n] where D[i] is the degree of vertex i.

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

# hdag
#
# Hamiltonian Path in DAG 
#
# Key idea: we can always form a topological sort of a DAG. Once a graphs has been
#           sorted topologically, the only way it can fail to be Hamiltonian is if
#           there exists a pair of nodes in the sorted data that is not an edge.

def hdag(graph):
    def clone(adj):
        copy = {}
        for k,v in adj.items():
            copy[k]=v
        return copy
    _,_,adj = create_adjacency(graph,back=False,self=False)
    ordered = topological_order(clone(adj)) # use clone because topological_order destroys its parameter
    for a,b in zip(ordered[:-1],ordered[1:]):
        if not b in adj[a]:
            return (-1,[])
    return (1,ordered)

# sq
#
# Square in a Graph

def sq(g):
    def explore(node,path=[]):
        if len(path)>4: return -1
 
        for next_node in adj[node]:
            if next_node==path[0] and len(path)==4:
                #print (path+[next_node])
                return 1
            if next_node in path: continue
            if explore(next_node,path+[next_node])==1: return 1
        return -1
    
        
    adj = {}
    for a,b in g[1:]:
        if not a in adj: adj[a]=[]
        adj[a].append(b)
        if not b in adj: adj[b]=[]
        adj[b].append(a)
    
    for node in adj.keys():
        if explore(node,path=[node])==1:
            return 1
        
    return -1

# dfs
#
# Depth First Search
#
# Inputs:    adj
#            n
#            sequence 
#            previsit 
#            postvisit  
#            preexplore
def dfs(adj          = None,
        n            = None,
        sequence     = None,
        previsit     = lambda v:None,
        postvisit    = lambda v:None,
        preexplore   = lambda v:None,
        list_visited = True):

    def explore(v):
        visited[v] = True
        previsit(v)

        for u in adj[v]:
            if not visited[u]:
                explore(u)

        postvisit(v)   

    visited = [False for v in range(n+1)]  # Zero element won't be used, but it does simplify indexing

    for v in sequence:
        if not visited[v]:
            preexplore(v)
            explore(v)
    return [i for i in range(1,len(visited)) if visited[i]==list_visited]

def create_adj(edges,reverse=False):
    n,_=edges[0]
    adj = {}
    for i in range(1,n+1):
        adj[i]=[]
    for (a,b) in edges[1:]:
        if reverse:
            (a,b)=(b,a)
        adj[a].append(b)
    return adj


#  scc
#
# Strongly Connected Component
#
# Input: A simple directed graph with fewer than 1000 vertices in the edge list format.
#
# Return: The number of strongly connected components in the graph.

def scc(edges):

    def counter():
        count = 0
        while True:
            yield count
            count+=1
            
    def decreasing(post):
        pairs = sorted(zip(post,range(1,len(post)+1)),reverse=True)
        for a,b in pairs:
            yield b

    def incr_pre(v):
        pre[v]     = next(clock)

    def incr_post(v):
        post[v]     = next(clock)

    def incr_ccnum(v):
        ccnum[v-1]  = next(cc)

    n,_     = edges[0]
    clock   = counter()
    pre     = [-1 for v in range(n+1)]   # Zero element won't be used, but it does simplify indexing
    post    = [-1 for v in range(n+1)]   # Zero element won't be used, but it does simplify indexing 
    adj_R      = create_adj(edges,reverse=True)
    dfs(adj_R,n,
      sequence = range(1,n+1),
      previsit = incr_pre,
      postvisit = incr_post)

    cc       = counter()
    ccnum    = [-1 for v in range(n+1)]   # Zero element won't be used, but it does simplify indexing
    adj         = create_adj(edges)
    dfs(adj,n,
      sequence    = decreasing(post[1:]),
      preexplore  = incr_ccnum)

    return [cc+1 for cc in ccnum if cc>-1],adj_R,adj 

#    sc
#
#    Semi-Connected Graph
#    A directed graph is semi-connected if for all pairs of vertices i,j there is either a path from i to j or a path from j to i
# Input: a simple directed graphs with at most 1 thousand
#               vertices each in the edge list format.
#
# Return: 1 if the graph is semi-connected and -1  otherwise.

def sc(edges):
    n,_ = edges[0]
    ccnum,adj,adjr = scc(edges)
    pairs = {}
    for i in range(len(ccnum)):
        for j in range(i+1,len(ccnum)):
            pairs[(ccnum[i],ccnum[j])] = False
    for component in ccnum:
        visited = dfs(adj,n,sequence=(i for i in [component]),list_visited=True)
        visitedr = dfs(adjr,n,sequence=(i for i in [component]),list_visited=True)
        for node in visited:
            a,b = min(node,component),max(node,component),
            pairs[(a,b)]=True
        for node in visitedr:
            a,b = min(node,component),max(node,component),
            pairs[(a,b)]=True      
    for k,v in pairs.items():
        if not v: return -1
    return 1