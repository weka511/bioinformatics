#!/usr/bin/env python
#    Copyright (C) 2019-2023 Greenweaves Software Limited
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
#    along with this program.  If not, see <http://www.gnu.org/licenses/>

''' Functions to solve graphs problems from http://rosalind.info/problems/topics/graphs/'''

from functools   import reduce
from collections import deque
from math        import isinf
from unittest    import TestCase, main

# from   tarjan      import tarjan
from   deprecated import deprecated
from helpers     import create_adjacency
from align       import create_topological_order
from rosalind    import grph

def bf(edges,
       s=1):    # Set to zero for nwc - see Ayaan Hossain's comment http://rosalind.info/problems/nwc/questions/
                # For anyone having a problem solving this, observe that in the problem statement,
                #unlike the Bellman-Ford Problem it has not been mentioned that you have to take
                #  vertex 1 as the source node. If you do take vertex 1 or any other vertex for that matter,
                # as the source node and if there is no out-going edge from that vertex or
                # if the negative-weight cycle is unreachable from that vertex, then there
                # will be no way to possibly detect the negative-weight cycle.
                # In this case we have to come up with a strategy to ensure that our source node
                # is a vertex from which there are out-going edges and that the negative-weight
                # cycle is reachable from our source node. The simplest way to do this is to add
                # a dummy "source_node" to our graph and add an edge from this "source_node"
                # to every other vertex in the graph with cost/weight/length 0.
                # Then if we begin our Bellman-Ford Algorithm on this "source_node", surely enough
                # we will be able to detect any negative-weight cycle in the graph if it is present.
    '''
    bf

    Bellman-Ford Algorithm

    Given: A simple directed graph with integer edge weights and fewer than 1,000 vertices in the edge list format.

    Return: neg    -1 if there is a negative cycle, else +1
            dists   distances from origin to each node
            pre     predecessor of each node
    '''
    n,_         = edges[0]

    if s==0:
        for i in range(1,n):
            edges.append([0,i,0])

    dist        = [float('inf') for i in range(n+1)]
    predecessor = [None for i in range(n+1)]

    dist[s]  = 0 # distance of source from itself is zero

    # See https://stackoverflow.com/questions/28857918/bellman-ford-algorithm-explanation
    # Bellman--Ford has two relevant invariants that hold for all vertices u.
    # There exists a path from the source to u of length dist[u] (unless dist[u] is INT_MAX).
    # After i iterations of the outer loop, for all paths from the source to u with
    # i or fewer edges, the length of that path is no less than dist[u].
    #
    # After V-1 iterations, the second invariant implies that no simple path
    # from the source to u is shorter than dist[u]. The first hence implies
    # that the path that we found is shortest.

    for _ in range(n):
        for u,v,w in edges[1:]:
            if dist[u] + w     < dist[v]:
                dist[v]        = dist[u] + w
                predecessor[v] = u

    return (max(1 if dist[u]+w<dist[v] else -1 for u,v,w in edges[1:]),
            dist[1:],
            predecessor[1:])


def bip(graph):
    '''
    bip

    Input: a graph in edgelist form

    Output:  1 if graph is bipartite, otherwise -1

    Strategy: try to partition graph into two colours
    '''
    red  = set()
    blue = set()

    def colour(node,isBlue):
        '''
        colour

         Attempt to assign this node, and all reachable nodes, to one colour or t'other

         Inputs:  node    The node we are assigning
                  isBlue  Indicates whether we are trying Red or Blue
        '''
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
# Input: A simple graph with fewer than 1,000 vertices in the edge list format.
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



def dij(g):
    '''
     DIJ  Dijkstra's Algorithm: compute single-source shortest distances
    in a directed graph with positive edge weights.

    I have followed the version described in Wikipedia -- https://en.wikipedia.org/wiki/Dijkstra's_algorithm
    1. Mark all nodes unvisited. Create a set of all the unvisited nodes called the unvisited set.
    2. Assign to every node a tentative distance value: set it to zero for our initial node and to
       infinity for all other nodes. Set the initial node as current.
    3. For the current node, consider all of its unvisited neighbours and calculate their
        tentative distances through the current node. Compare the newly calculated tentative
        distance to the current assigned value and assign the smaller one. For example, if the
        current node A is marked with a distance of 6, and the edge connecting it with a neighbour
        B has length 2, then the distance to B through A will be 6 + 2 = 8. If B was previously
        marked with a distance greater than 8 then change it to 8. Otherwise, the current value will be kept.
    4. When we are done considering all of the unvisited neighbours of the current node,
       mark the current node as visited and remove it from the unvisited set. A visited
       node will never be checked again.
    5. If the destination node has been marked visited (when planning a route between
       two specific nodes) or if the smallest tentative distance among the nodes in the
       unvisited set is infinity (when planning a complete traversal; occurs when there
       is no connection between the initial node and remaining unvisited nodes), then stop. The algorithm has finished.
    6. Otherwise, select the unvisited node that is marked with the smallest tentative distance,
       set it as the new "current node", and go back to step 3.
    '''
    n,_             = g[0]                     # Number of nodes
    open_edges      = [edge for edge in g[1:]] # Edges that haven't been used-FIXME (issue #112(
    D               = [float('inf')]*(n+1)     # Distance from 1 to each node
    D[1]            = 0

    for _ in range(n):
        for a,b,w in open_edges:
            if not isinf(D[a]):
                D[b] = min(D[b],D[a] + w)

    return [d if not isinf(d) else -1 for d in D[1:]]

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
    ordered = get_topological_order(adj)
    for a,b in zip(ordered[:-1],ordered[1:]):
        if not b in adj[a]:
            return (-1,[])
    return (1,ordered)

# sq
#
# Square in a Graph
#
# Given: A positive integer k<=20 and k simple undirected graphs with n<=400 vertices in the edge list format.
#
# Return: For each graph, output 1 if it contains a simple cycle (that is, a cycle which doesn�t intersect
# itself) of length 4 and -1 otherwise.

def sq(g):
    def explore(node,path=[]):
        if len(path)>4: return -1

        for next_node in adj[node]:
            if next_node==path[0] and len(path)==4:  return 1
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
        if not a in adj:
            adj[a] = []
        adj[a].append(b)
    return adj


#  scc
#
# Strongly Connected Component
#
# Input: A simple directed graph with fewer than 1000 vertices in the edge list format.
#
# Return: A tuple ( ccs,adj,adj_R), where
#         nscc  The number of strongly connected components in the graph.
#         adj   Adjacency list from forward links
#         adj_r Adjacency list from reverse links

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

    # Find sink nodes by looking for source nodes of the reversed graph

    clock      = counter()
    pre        = [-1 for v in range(n+1)]   # Zero element won't be used, but it does simplify indexing
    post       = [-1 for v in range(n+1)]   # Zero element won't be used, but it does simplify indexing
    adj_R      = create_adj(edges,reverse=True)
    dfs(adj_R,n,
      sequence  = range(1,n+1),
      previsit  = incr_pre,
      postvisit = incr_post)

    # construct components

    cc       = counter()
    ccnum    = [-1 for v in range(n+1)]   # Zero element won't be used, but it does simplify indexing
    adj      = create_adj(edges)
    dfs(adj,n,
      sequence    = decreasing(post[1:]),
      preexplore  = incr_ccnum)

    return len([cc+1 for cc in ccnum if cc>-1]),adj,adj_R

#    sc
#
#    Semi-Connected Graph
#
#    A directed graph is semi-connected if for all pairs of vertices i,j there is either a path from i to j or a path from j to i
#
# Input: a simple directed graphs with at most 1 thousand
#               vertices in the edge list format.
#
# Return: 1 if the graph is semi-connected and -1  otherwise.

def sc(edges):
    n,_            = edges[0]
    ccnum,adj,adjr = scc(edges)  # Graph is semi connected iff the DAG of its strongly connected components is

    pairs          = {}   # Build a structure of pairs of strongly connected components
    for i in range(len(ccnum)):
        for j in range(i+1,len(ccnum)):
            pairs[(ccnum[i],ccnum[j])] = False

    for component in ccnum:
        for node in dfs(adj,n,sequence=(i for i in [component])):  # find links in one direction
            a,b = min(node,component),max(node,component),
            pairs[(a,b)]=True
        for node in dfs(adjr,n,sequence=(i for i in [component])): # now check t'other
            a,b = min(node,component),max(node,component),
            pairs[(a,b)]=True

    for k,v in pairs.items():  # Is any pair unlinked?
        if not v: return -1

    return 1

# gs
#
# General Sink
#
# Input: a simple directed graphs with at most 1,000 vertices and 2,000 edges in the edge list format.
#
# Return: a vertex from which all other vertices are reachable (if such a vertex exists) and -1 otherwise.
#
#
def gs(edges):
    n,_         = edges[0]
    ccnum,adj,_ = scc(edges) # We need test only one vertex from each simply connected component

    for component in ccnum:
        if len(dfs(adj,n,sequence=(i for i in [component])))==n:
            return component

    return -1

# SDAG Shortest Paths in DAG

def sdag(m,adjacency,weights):
    t    = topological_order(adjacency.copy())
    D    = [None]*(m+1)
    D[1] = 0
    for i  in t:
        if D[i] == None: continue
        for j in adjacency[i]:
            if D[j]==None:
                D[j] = D[i] +weights[(i,j)]
            else:
                trial = D[i] +weights[(i,j)]
                if trial < D[j]:
                    D[j] = trial
    return D[1:]

# two_sat
#
# 2SAT 2-Satisfiability
#
# Given: A positive integer k<=20 and k 2SAT formulas represented as follows.
# The first line gives the number of variables n and the number of clauses m,
# each of the following m lines gives a clause of length 2 by specifying two different literals.
# Return: For each formula, output 0 if it cannot be satisfied or 1 followed by a satisfying assignment otherwise.
#
# See description in Wikipedia -- https://en.wikipedia.org/wiki/2-satisfiability#Strongly_connected_components
#
# Aspvall, Bengt; Plass, Michael F.; Tarjan, Robert E. (1979),
# "A linear-time algorithm for testing the truth of certain quantified boolean formulas"
#  Information Processing Letters, 8 (3): 121�123, doi:10.1016/0020-0190(79)90002-4.

def two_sat(problem):

    # create_assignment
    #
    # Assign values to variable, given the condensation of the implication graph,
    # a smaller graph that has one vertex for each strongly connected component,
    # and an edge from component i to component j whenever the implication graph
    # contains an edge uv such that u belongs to component i and v belongs to component j.
    # The condensation is automatically a directed acyclic graph and,
    # like the implication graph from which it was formed, it is skew-symmetric.
    def create_assignment(condensation):
        assignment = []
        for component in condensation:
            for v in component:
                if not v in assignment and not -v in assignment:
                    assignment.append(v)
        return sorted(assignment,key=lambda x: abs(x))

    n,m,clauses = problem
    edges       = [(-a,b) for a,b in clauses] + [(-b,a) for a,b in clauses]
    scc         = tarjan(create_adj([[n,len(edges)]] + edges))
    for component in scc:
        for i in range(len(component)):
            for j in range(i+1,len(component)):
                if component[i]==-component[j]:
                    return 0,[]

    return 1,create_assignment(scc)





#  SUFF Creating a Suffix Tree

def createTrie(string):

    class Node:

        seq = 0

        def __init__(self):
            self.seq       = Node.seq
            self.edges     = {}
            self.positions = {}
            self.label     = None
            Node.seq       += 1

        def bfs(self,
                path          = [],
                propagatePath = lambda path,symbol: path,
                visitLeaf     = lambda label,path:None,
                visitInternal = lambda node,symbol,position,path: None):
            if len(self.edges)==0:
                visitLeaf(self.label,propagatePath(path,''))
            else:
                for symbol,node in self.edges.items():
                    visitInternal(node,symbol,self.positions[symbol],path)
                    node.bfs(propagatePath(path,symbol),
                             propagatePath,
                             visitInternal=visitInternal,
                             visitLeaf=visitLeaf)



    Trie = Node()
    for i in range(len(string)):
        currentNode = Trie
        for j in range(i,len(string)):
            currentSymbol = string[j]
            if currentSymbol in currentNode.edges:
                currentNode = currentNode.edges[currentSymbol]
            else:
                newNode                              = Node()
                currentNode.edges[currentSymbol]     = newNode
                currentNode.positions[currentSymbol] = j
                currentNode                          = newNode
        if len(currentNode.edges)==0:
            currentNode.label = i
    return Trie

def convertToTree(Trie):
    Starts = []
    def identifyBranchPoints(node,symbol,position,prefix):
        if len(node.edges)>1:
            Starts.append(node)
            #print (symbol,position,node.seq,len(node.edges))

    Trie.bfs(visitInternal=identifyBranchPoints)

    for node in Starts:
        branches  = {}
        jumps     = {}
        positions = {}
        for symbol,nextNode in node.edges.items():
            candidateBranch = [symbol]
            while len(list(nextNode.edges.items()))==1:
                nextSymbol    = list(nextNode.edges.keys())[0]
                nextNode      = nextNode.edges[nextSymbol]
                candidateBranch.append(nextSymbol)
            if len(candidateBranch)>1:
                branches[symbol]            = ''.join(candidateBranch)
                jumps[branches[symbol]]     = nextNode
                positions[branches[symbol]] = node.positions[symbol]
        for symbol,symbols in branches.items():
            node.edges[symbols]     = jumps[symbols]
            node.positions[symbols] = positions[symbols]
            del node.edges[symbol]
            del node.positions[symbol]

def suff(string):
    Trie = createTrie(string)
    convertToTree(Trie)
    Result = []
    Trie.bfs(path        = '',
            propagatePath = lambda path,symbol: path,
            visitLeaf     = lambda label,prefix:None,
            visitInternal = lambda node,symbol,position,prefix: Result.append(symbol))

    return Result



#    bfs

#   Perform Breadth First search

def bfs(graph,                           # Graph to be searched
        start,                          # Starting node
        isGoal    = lambda v   : False, # Used to terminate search
                                        # before visiting all nodes
        visit     = lambda v   : None,  # Action to be performed when node processed
        pre_visit = lambda v,w : None,  # Action to be performed when node discovered
        selector  = 0                   # set to -1 for depth first instead of breadth first
        ):
    Q               = [start]
    discovered      = {start}
    while len(Q)>0:
        v = Q.pop(selector)
        visit(v)
        if isGoal(v): return v
        for w in graph[v]:
            if w not in discovered:
                Q.append(w)
                discovered.add(w)
                pre_visit(v,w)

# ShortestDistances
#
# use breadth-first search to compute single-source shortest distances in an unweighted directed graph.

def ShortestDistances(graph):

    def update_distance(v,w):
        distances[w-1] = distances[v-1]+1

    def start(v):
        if v>1: return
        distances[0] = 0

    max_node,graph    = graph
    distances        = [-1]*max_node
    bfs(graph,1,visit=start,pre_visit=update_distance)
    return distances


### Deprecated code ###

#TRIE  Pattern matching

# trie
#
# Given: A list of at most 100 DNA strings of length at most 100 bp, none of which is a prefix of another.
#         one_based Indicates whether numnering of nodes should start at 1 or zero
#
# Return: The adjacency list corresponding to the trie T
#         for these patterns, in the following format.
#         If T has n nodes, first label the root with 1 and then label the remaining nodes
#         with the integers 2 through n in any order you like.
#         Each edge of the adjacency list of T will be encoded by a triple
#         containing the integer representing the edge's parent node, followed by the integer
#         representing the edge's child node, and finally the symbol labeling the edge.
@deprecated(reason="Use snp.create_tree instead")
def trie(strings,one_based=True):

    def find_string_in_adjacency_list(string, adjacency_list):
        index=0
        parent=0
        path=[]
        for cc in string:
            matched=False
            while index<len(adjacency_list) and not matched:
                a,b,c=adjacency_list[index]
                if a==parent and c==cc:
                    matched=True
                    path.append(index)
                    parent=b
                else:
                    index+=1
        return (parent,path)

    def create_suffix(parent,path,string,b):
        result=[]
        a = parent
        for i in range(len(path),len(string)):
            b+=1
            result.append((a,b,string[i]))
            a=b

        return result

    def merge_string_with_adjacency_list(adjacency_list,string):
        parent,path=find_string_in_adjacency_list(string, adjacency_list)
        return adjacency_list+create_suffix(parent,path,string,len(adjacency_list))

    incr = 1 if one_based else 0

    def increment_indices(adjacency_list):
        return [(a+incr,b+incr,c) for (a,b,c) in adjacency_list]

    return increment_indices(reduce(merge_string_with_adjacency_list,strings,[]))

if __name__=='__main__':
    class Test_graphs(TestCase):

        def test_bf(self):
            edges = [(9, 13),
                     (1, 2, 10),
                     (3, 2, 1),
                     (3, 4, 1),
                     (4, 5, 3),
                     (5, 6, -1),
                     (7, 6, -1),
                     (8, 7, 1),
                     (1, 8, 8),
                     (7, 2, -4),
                     (2, 6, 2),
                     (6, 3, -2),
                     (9, 5 ,-10),
                     (9, 4, 7)]
            _,dists,_ = bf(edges)

            self.assertEqual([0, 5, 5, 6, 9, 7, 9, 8],dists[:-1])
            self.assertEqual(float('inf'),dists[-1])

        def test_bip(self):
            self.assertEqual(-1,bip([(3, 3),
                                     (1, 2),
                                     (3, 2),
                                     (3, 1)]))

            self.assertEqual(1,bip([(4, 3),
                                    (1, 4),
                                    (3, 1),
                                    (1, 2)]))

        def test_cc(self):
            count,_ = cc([(12, 13),
                          ( 1, 2),
                          ( 1, 5),
                          ( 5, 9),
                          (5, 10),
                          (9, 10),
                          (3, 4),
                          (3, 7),
                          (3, 8),
                          (4, 8),
                          (7, 11),
                          (8, 11),
                          (11, 12),
                          (8, 12)])
            self.assertEqual(3,count)

        def test_dij(self):
            self.assertEqual([0, 3, 2, 5, 6, -1],
                             dij([[6 ,10],
                                  [1, 2, 4],
                                  [1, 3, 2],
                                  [2, 3, 3],
                                  [6, 3, 2],
                                  [3, 5, 5],
                                  [5, 4, 1],
                                  [3, 2, 1],
                                  [2, 4, 2],
                                  [2, 5, 3]]))

        def test_nwc(self):
            n1,_,_ =bf([[4, 5],
                       [1, 4, 4],
                       [4, 2, 3],
                       [2, 3, 1],
                       [3, 1, 6],
                       [2, 1, -7]])
            self.assertEqual(-1,n1)
            n2,_,_ = bf([[3, 4],
                         [1, 2, -8],
                         [2, 3, 20],
                         [3, 1, -1],
                         [3, 2, -30]])
            self.assertEqual(1,n2)

    main()
