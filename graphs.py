#!/usr/bin/env python

#    Copyright (C) 2019-2024 Greenweaves Software Limited
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

from functools import reduce
from collections import deque
from math import isinf
from unittest import TestCase, main, skip
import numpy as np
from numpy.testing import assert_array_equal
from helpers import create_adjacency
from align import create_topological_order
from tarjan import tarjan

def bf(edges,s=1):
    '''
    Bellman-Ford Algorithm, and
    NWC Negative Weight Cycle

    Check whether a given graph contains a cycle of negative weight.
    Given: A positive integer k, and k simple directed graphs with integer edge weights from -1000 to 1000
    Return: For each graph, output 1 if it contains a negative weight cycle and -1 otherwise.

    Parameters:
        edges    A simple directed graph with integer edge weights and fewer than 1,000 vertices in the edge list format.
        s        Set to zero for nwc

    Return: neg    -1 if there is a negative cycle, else +1
            dists   distances from origin to each node
            pre     predecessor of each node
    '''
    n,_         = edges[0]

    # see Ayaan Hossain's comment http://rosalind.info/problems/nwc/questions/
    # ..., observe that in the problem statement, unlike the Bellman-Ford Problem,
    # it has not been mentioned that you have to take vertex 1 as the source node.
    # If you do take vertex 1 or any other vertex for that matter,
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
    if s == 0:
        for i in range(1,n):
            edges.append([0,i,0])

    dist = np.full((n+1),float('inf'))
    predecessor = np.full((n+1),-1,dtype=int)
    dist[s]  = 0 # distance of source from itself is zero

    # See https://stackoverflow.com/questions/28857918/bellman-ford-algorithm-explanation
    # Bellman--Ford has two relevant invariants that hold for all vertices u.
    # There exists a path from the source to u of length dist[u] (unless dist[u] is INT_MAX).
    # After i iterations of the outer loop, for all paths from the source to u with
    # i or fewer edges, the length of that path is no less than dist[u].

    # After V-1 iterations, the second invariant implies that no simple path
    # from the source to u is shorter than dist[u]. The first hence implies
    # that the path that we found is shortest.

    for _ in range(n):
        for u,v,w in edges[1:]:
            if dist[u] + w < dist[v]:
                dist[v] = dist[u] + w
                predecessor[v] = u

    return (max(1 if dist[u] + w < dist[v] else -1 for u,v,w in edges[1:]), dist[1:], predecessor[1:])


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


def cc(graph):
    '''
    Connected Components

    Input: A simple graph with fewer than 1,000 vertices in the edge list format.

    Return: The number of connected components in the graph.
    '''
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

def cte(edges):
    '''
    CTE Shortest Cycle Through a Given Edge

    Returns: For each graph, output the length of a shortest cycle going
           through the first specified edge if there is a cycle and "-1" otherwise.
    '''

    def relabel(edges):
        '''
         Exploit Dijkstra's algorithm by relabelling nodes so that first edge starts at 1
        '''
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


def dag(graph):
    '''
    dag  Testing Acyclicity

    Input: a graph

    Output: 1 if graph acyclic, otherwise -1
    '''

    def explore(node,path=[]):
        '''
        explore: find cyclic path

        Input: node
               path

        Output: True iff there is a cyclic path starting at node, and
             whole 1st elements coincide with path
        '''
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

def ddeg(n,m,A):
    '''
    Double Degree Array

    Input : A simple graph with fewer than 1,000 vertices in the edge list format.

    Return: An array D[1..n] where D[i] is the sum of the degrees of i's neighbors.
    '''
    lookup = deg(n,m,A)
    sums = np.zeros((n),dtype=int)
    for (a,b) in A:
        sums[a-1] += lookup[b-1]
        sums[b-1] += lookup[a-1]
    return sums

def deg(n,m,A):
    '''
     Degree array

     Input: A simple graph with fewer than 1,000 vertices in the edge list format.

     Return: An array D[1..n] where D[i] is the degree of vertex i.
    '''
    degrees = np.zeros((n), dtype=int)
    for (a,b) in A:
        degrees[a-1] += 1
        degrees[b-1] += 1
    return degrees



def dij(g):
    '''
     DIJ  Dijkstra's Algorithm: compute single-source shortest distances in a directed graph with positive edge weights.

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
    n,_ = g[0]                     # Number of nodes
    unvisited = list(range(n+1))
    D  = np.full((n+1),float('inf'))     # Distance from 1 to each node
    D[1] = 0
    Adj = np.full((n+1,n+1),float('inf'))  # Weighted adjacency matrix
    for i,j,distance in g[1:]:
        Adj[i,j] = distance

    while True:
        if len(unvisited) == 0: break
        unvisited_distances = [D[i] for i in unvisited]
        i_current = np.argmin(unvisited_distances)
        if np.isinf(unvisited_distances[i_current]): break
        current_node = unvisited[i_current]
        d_current = D[current_node]
        unvisited.remove(current_node)
        for j in unvisited:
            D[j] = min(D[j],d_current + Adj[current_node,j])

    return [d if not isinf(d) else -1 for d in D[1:]]


def hdag(graph):
    '''
    Hamiltonian Path in DAG

    Key idea: we can always form a topological sort of a DAG. Once a graphs has been
              sorted topologically, the only way it can fail to be Hamiltonian is if
              there exists a pair of nodes in the sorted data that is not an edge.
    '''

    _,_,adj = create_adjacency(graph,back=False,self=False)
    ordered = create_topological_order(adj)
    for a,b in zip(ordered[:-1],ordered[1:]):
        if not b in adj[a]:
            return (-1,[])
    return (1,ordered)


def sq(g):
    '''
    Square in a Graph

    Given: A positive integer k<=20 and k simple undirected graphs with n<=400 vertices in the edge list format.

    Return: For each graph, output 1 if it contains a simple cycle (that is, a cycle which doesn�t intersect
    itself) of length 4 and -1 otherwise.
    '''
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


def dfs(adj = None,
        n = None,
        sequence = [],
        previsit = lambda v:None,
        postvisit = lambda v:None,
        preexplore = lambda v:None):
    '''
    Depth First Search

    Parameters:
        adj        Adjacency list for graph
        n          Number of nodes in graph
        sequence   Optional parameter to force some nodes to be visited first
        previsit   Function to be executed before visiting each node
        postvisit  Function to be executed after visiting each node
        preexplore Function to be executed before visiting each node if sequence specified
    '''
    def explore(v):
        '''
        Perform depth-first search recursively
        '''
        visited[v] = True
        previsit(v)

        for u in adj[v]:
            if not visited[u]:
                explore(u)

        postvisit(v)

    visited = np.full((n+1),False,dtype=bool)  # Keep track of which nodes have been visited
                                               # Zero element won't be used, but it does simplify indexing

    # Visit any nodes that need to be given priority
    for v in sequence:
        if not visited[v]:
            preexplore(v)
            explore(v)

    # Visit any nodes that have been left over
    return [i for i in range(1,len(visited)) if visited[i]]

def create_adj(edges,reverse=False):
    '''
    Create adjacency list for graph

    Parameters:
        edges      Graph in edge list format
        reverse    Used to reverse all edges

    Returns:
       Dict containting adjacency list: a->[b,c,d...] where edges are (a,b), (a,c), (a,d),...
    '''
    n,_= edges[0]
    Product = {i:[] for i in range(1,n+1)}
    for (a,b) in edges[1:]:
        if reverse:
            (a,b) = (b,a)
        if not a in Product:
            Product[a] = []
        Product[a].append(b)
    return Product


def scc(edges):
    '''
    Strongly Connected Component

    Input: A simple directed graph with fewer than 1000 vertices in the edge list format.

    Return: A tuple ( ccs,adj,adj_R), where
        nscc  The number of strongly connected components in the graph.
        adj   Adjacency list from forward links
        adj_r Adjacency list from reverse links
    '''
    def VisitCounter():
        '''
        Generator used to keep track of the times when nodes are visited or departed from
        '''
        count = 0
        while True:
            yield count
            count += 1

    def decreasing(post):
        pairs = sorted(zip(post,range(1,len(post)+1)),reverse=True)
        for a,b in pairs:
            yield b

    def incr_pre(v):
        '''
        Used during depth-first traversal of reversed adjacency list
        to record time when node visited
        '''
        pre[v] = next(clock)

    def incr_post(v):
        '''
        Used during depth-first traversal of reversed adjacency list
        to record time when node departed from
        '''
        post[v] = next(clock)

    def incr_ccnum(v):
        '''
        Used during depth-first traversal of adjacency list
        to record time when node visited
        '''
        ccnum[v-1] = next(clock)

    n,_ = edges[0]

    # Find sink nodes by looking for source nodes of the reversed graph

    clock = VisitCounter()
    pre = np.full((n+1),-1)   # Used to keep track of time of first visit to each node
                              # see https://rosalind.info/glossary/algo-depth-first-search/
                              # Zero element won't be used, but it does simplify indexing
    post = np.full((n+1),-1)  # Used to keep track of time of last departure from each node
                              # see https://rosalind.info/glossary/algo-depth-first-search/
                              # Zero element won't be used, but it does simplify indexing
    adj_R  = create_adj(edges,reverse=True)
    dfs(adj_R,n,
      sequence = range(1,n+1),
      previsit = incr_pre,
      postvisit = incr_post)

    # Run the undirected connected components algorithm on G, and during the depth-first search,
    # process the vertices in decreasing order of their post numbers from step 1.

    ccnum = np.full(n,-1)   # Used to record one node from each connected component
    adj = create_adj(edges)
    dfs(adj,n,
      sequence = decreasing(post[1:]),
      preexplore = incr_ccnum)

    return (ccnum > -1).sum(),adj,adj_R

def create_node2component(m,n,cc):
    '''Used with connect components to map each node onto the corresponding component'''
    Product = np.zeros((n+1), dtype=np.int64)
    for i in range(m):
        for j in cc[i]:
            Product[j] = i + 1
    return Product

def create_adj_for_cc(m,edges,cc_index):
    '''Used with connect components to create an adjacency matrix for components'''
    Product = np.eye(m, dtype=np.int64)
    for i,j in edges[1:]:
        Product[cc_index[i]-1,cc_index[j]-1] = 1
    return Product

def iterate_adj(m,A):
    '''
    Used with connect components to to iterate adjacency matrix
    to work out which components are accessible
    '''
    T = np.eye(m, dtype=np.int64)
    for _ in range(m):
        np.matmul(A,T,out=T)
    return T

def sc(edges):
    '''
     Semi-Connected Graph

    A directed graph is semi-connected if for all pairs of vertices i,j there is either a path from i to j or a path from j to i.
    We shall use the observation that graph is semi-connected iff graph of connected components is semi-connected

    Input: a simple directed graphs with at most 1 thousand
               vertices in the edge list format.

    Return: 1 if the graph is semi-connected and -1  otherwise.
    '''
    n,_ = edges[0]
    adj = create_adj(edges)
    # Create connected components
    cc = tarjan(adj)
    m = len(cc)
    T =  iterate_adj(m,
                     create_adj_for_cc(m,edges,
                                       create_node2component(m,n,cc)))

    for i in range(m):
        for j in range(m):
            if T[i,j] == 0 and T[j,i] == 0: return -1 # No path between i and j

    return 1


def gs(edges):
    '''
    General Sink

    Input: a simple directed graph with at most 1,000 vertices and 2,000 edges in the edge list format.

    Return: a vertex from which all other vertices can be reached (if such a vertex exists) and -1 otherwise.
    '''
    n,_ = edges[0]
    adj = create_adj(edges)
    cc = tarjan(adj)
    m = len(cc)
    T =  iterate_adj(m,
                     create_adj_for_cc(m,edges,
                                       create_node2component(m,n,cc)))
    T1 = np.clip(T,a_min=0,a_max=1)
    T2 = np.sum(T1,axis=1)
    index_gs = np.argwhere(T2 == m)  # Can reach any node from here
    return -1 if len(index_gs)==0 else index_gs[0] + 1


def sdag(m,adjacency,weights):
    '''
    DAG Shortest Paths in DAG
    '''
    t  = create_topological_order(adjacency.copy())
    D = [None]*(m+1)
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


def two_sat(problem):
    '''
    2SAT 2-Satisfiability

    Given: A positive integer k<=20 and k 2SAT formulas represented as follows.
    The first line gives the number of variables n and the number of clauses m,
    each of the following m lines gives a clause of length 2 by specifying two different literals.
    Return: For each formula, output 0 if it cannot be satisfied or 1 followed by a satisfying assignment otherwise.

    See description in Wikipedia -- https://en.wikipedia.org/wiki/2-satisfiability#Strongly_connected_components

    Aspvall, Bengt; Plass, Michael F.; Tarjan, Robert E. (1979),
    "A linear-time algorithm for testing the truth of certain quantified boolean formulas"
    Information Processing Letters, 8 (3): 121-123, doi:10.1016/0020-0190(79)90002-4.
    '''
    def create_assignment(condensation):
        '''
         Assign values to variable, given the condensation of the implication graph,
         a smaller graph that has one vertex for each strongly connected component,
         and an edge from component i to component j whenever the implication graph
         contains an edge uv such that u belongs to component i and v belongs to component j.
         The condensation is automatically a directed acyclic graph and,
         like the implication graph from which it was formed, it is skew-symmetric.
        '''
        assignment = []
        for component in condensation:
            for v in component:
                if not v in assignment and not -v in assignment:
                    assignment.append(v)
        return sorted(assignment,key=lambda x: abs(x))

    n,m,clauses = problem
    edges = [(-a,b) for a,b in clauses] + [(-b,a) for a,b in clauses]
    scc = tarjan(create_adj([[n,len(edges)]] + edges))
    for component in scc:
        for i in range(len(component)):
            for j in range(i+1,len(component)):
                if component[i]==-component[j]:
                    return 0,[]

    return 1,create_assignment(scc)

def createTrie(string):
    '''
    SUFF Creating a Suffix Tree
    '''
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
    '''Creating a Suffix Tree'''
    Trie = createTrie(string)
    convertToTree(Trie)
    Result = []
    Trie.bfs(path        = '',
            propagatePath = lambda path,symbol: path,
            visitLeaf     = lambda label,prefix:None,
            visitInternal = lambda node,symbol,position,prefix: Result.append(symbol))

    return Result



def bfs(graph,                           # Graph to be searched
        start,                          # Starting node
        isGoal    = lambda v   : False, # Used to terminate search
                                        # before visiting all nodes
        visit     = lambda v   : None,  # Action to be performed when node processed
        pre_visit = lambda v,w : None,  # Action to be performed when node discovered
        selector  = 0                   # set to -1 for depth first instead of breadth first
        ):
    '''
     Perform Breadth First search
    '''
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

def ShortestDistances(graph):
    '''
    Use breadth-first search to compute single-source shortest distances in an unweighted directed graph.
    '''
    def update_distance(v,w):
        distances[w-1] = distances[v-1]+1

    def start(v):
        if v > 1: return
        distances[0] = 0

    max_node,graph = graph
    distances = [-1]*max_node
    bfs(graph,1,visit=start,pre_visit=update_distance)
    return distances

if __name__=='__main__':
    class Test_graphs(TestCase):

        def test_2sat(self):
            '''2-Satisfiability'''
            status,Solution = two_sat((2, 4, [(1, 2), (-1, 2), (1, -2), (-1, -2)]) )
            self.assertEqual(0,status)
            status,Solution = two_sat((3, 4, [(1, 2), (2, 3), (-1, -2), (-2, -3)]))
            self.assertEqual(1,status)
            self.assertEqual([1, -2, 3],Solution)


        def test_bf(self):
            '''Bellman-Ford Algorithm,'''
            _,dists,_ = bf([(9, 13),
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
                            (9, 4, 7)])
            assert_array_equal([0, 5, 5, 6, 9, 7, 9, 8],dists[:-1])
            self.assertEqual(float('inf'),dists[-1])


        def test_bfs(self):
            ''' use breadth-first search to compute single-source shortest distances in an unweighted directed graph.'''
            def create_tree(links):
                max_node,_ = links[0]
                Product = {a:[] for a in range(1,max_node+1)}
                for (a,b) in links[1:]:
                    Product[a].append(b)
                return (max_node,Product)

            self.assertEqual([0, -1, 2, 1, 3, 2],
                             ShortestDistances(
                                 create_tree([(6, 6),
                                              (4, 6),
                                              (6, 5),
                                              (4, 3),
                                              (3, 5),
                                              (2, 1),
                                              (1, 4)])))

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

        def test_cte(self):
            ''' CTE Shortest Cycle Through a Given Edge'''
            self.assertEqual(-1,
                             cte([[4, 5],
                                  [2, 4, 2],
                                  [3, 2, 1],
                                  [1, 4, 3],
                                  [2, 1, 10],
                                  [1, 3, 4]]))
            self.assertEqual(10,cte([[4, 5],
                                [3, 2, 1],
                                [2, 4, 2],
                                [4, 1, 3],
                                [2, 1, 10],
                                [1, 3, 4]]))


        def test_dag(self):
            '''dag  Testing Acyclicity'''
            self.assertEqual(1,dag([[2, 1],
                                    [1, 2]]))
            self.assertEqual(-1,dag([[4, 4],
                                     [4, 1],
                                     [1, 2],
                                     [2, 3],
                                     [3, 1]]))
            self.assertEqual(1,dag([[4, 3],
                                    [4, 3],
                                    [3, 2],
                                    [2, 1]]))

        def test_ddeg(self):
            ''' Double Degree Array'''
            assert_array_equal(
                np.array([3, 5, 5, 5, 0],dtype=int),
                ddeg(5, 4,
                     [[1, 2],
                     [2, 3],
                     [4, 3],
                     [2, 4]]))


        def test_deg(self):
            '''Degree array'''
            assert_array_equal(
                np.array([2, 4, 2, 2, 2, 2],dtype=int),
                deg(6, 7,
                    [[1, 2],
                     [2, 3],
                     [6, 3],
                     [5, 6],
                     [2, 5],
                     [2, 4],
                     [4, 1]]))


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

        def test_gs(self):
                ''' General Sink'''
                self.assertEqual(3,
                                 gs([[3,2],[3, 2],[2, 1]]))
                self.assertEqual(-1,
                                 gs([[3,2],[3, 2],[1, 2]]))

        def test_hdag(self):
            '''Hamiltonian Path in DAG'''
            self.assertEqual((1,[1,2,3]),
                             hdag([[3, 3],
                                   [1, 2],
                                   [2, 3],
                                   [1, 3]]))
            self.assertEqual((-1,[]),
                             hdag([[4, 3],
                                   [4, 3],
                                   [3, 2],
                                   [4, 1]]))

        def test_nwc(self):
            n1,_,_ =bf([[4, 5],
                       [1, 4, 4],
                       [4, 2, 3],
                       [2, 3, 1],
                       [3, 1, 6],
                       [2, 1, -7]],
                       s=0)
            self.assertEqual(-1,n1)
            n2,_,_ = bf([[3, 4],
                         [1, 2, -8],
                         [2, 3, 20],
                         [3, 1, -1],
                         [3, 2, -30]],
                        s=0)
            self.assertEqual(1,n2)


        def test_sc(self):
            '''Semi-Connected Graph'''
            self.assertEqual(1,sc( [(3, 2), (3 ,2), (2, 1)]))
            self.assertEqual(-1,sc( [(3, 2), (3 ,2), (1, 2)]))


        def test_scc(self):
            '''Strongly Connected Component '''
            nscc,adj,adj_r = scc([[6, 7],
                                  [4, 1],
                                  [1, 2],
                                  [2, 4],
                                  [5, 6],
                                  [3, 2],
                                  [5, 3],
                                  [3, 5]])
            self.assertEqual(3,nscc)

        def test_sdag(self):

            def create_adjacency(edges):
                m,n       = edges[0]

                product = {}
                weights = {}

                for a in range(1,m+1):
                    product[a]   = []

                for a,b,w in edges[1:]:
                    product[a].append(b)
                    weights[(a,b)] = w

                for a in product.keys():
                    product[a]=sorted(list(set(product[a])))

                return m,n,product,weights

            m,_,adjacency,weights = create_adjacency([[5, 6],
                [2, 3, 4],
                [4, 3, -2],
                [1, 4, 1],
                [1, 5, -3],
                [2, 4, -2],
                [5, 4, 1]
            ])
            self.assertEqual([0, None, -4, -2, -3],sdag(m,adjacency,weights))

        def test_sq(self):
            ''' Square in a Graph'''
            self.assertEqual(1,sq([(4, 5),
                                   (3, 4),
                                   (4, 2),
                                   (3, 2),
                                   (3, 1),
                                   (1, 2)]))
            self.assertEqual(-1,sq( [(4, 4),
                                     (1, 2),
                                     (3, 4),
                                     (2, 4),
                                     (4, 1)]))

        def test_suff(self):
            '''
            Creating a Suffix Tree

            I don't match the Sample test, but do on my first Rosalind sibmission,
            do I have used that as the basis for my test
            '''
            Edges = suff('TGCCTGAATGGGGGCTGCAGGATGGACGGCCGTCCATTGAAATGAATGGCCAGCTCTCGCCGGTGGAGGCTGTCGTCTAGTTACACGGCGTACTGCTTCAT'
                         'CAACTGGAAAGTCGTATTATGGAATGGATGTCAGCCAGTTCTCTGGTGAATCATTGAATTTAGACAAAGAGGCGTAGGCTCCGCGATAATGTTATGTTCCA'
                         'CATTACTTTTAATCGGTTGTCCCATCGAGGTCGGAAACTTGTGCTACAGCCACAGCCGCAATTCCGTAGACGGCCCCCATAACCTACCCTATGCCTAACCA'
                         'TCTGATTCCTCACGCACGCACCATGGTGGAAGAGAGAGTCCTGCTATGGGGGATTGGGTCGGGGTCGGTGTGGCCTGTGTTAGCTCCACAAGCAACGCATA'
                         'AAGTGCAGATACTTGGAGCTCAATTGGCCGCGCGCTGGCTAGTGATTTTCTATAGTGGTTTTAAAAGTCGGCCTGGTTTAGACAATTTGGCTCGCCGAAAT'
                         'TGCCGCTACCCGAGTGGCGCAAACGTCCATTCTCATCCCCAGGAATGAAACCGGGCGGACACTGCCATTCCCGCGAGAGATCCTCAGTGAAACTTTGTGGA'
                         'TAAGAATAAATTACCCTGCCACCGGTGATATTTACGAAGGCCTCGACATGCACAATGAGTCAAATGCACGGACCCGGTTACACGCTATAACATTCCGCTAC'
                         'TTCCTTTTGAATACGGTTTGTCGACTGGTTAATACGAGAAAGCCTGTCAAAGTAACCATGAGCGACCGTATGTTCTTTCCGCTGTAGCCCCATTGTTCCTT'
                         'ATTAG$')
            self.assertEqual(1327,len(Edges))
            self.assertIn('TACCCTATGCCTAACCATCTGATTCCTCACGCACGCACCATGGTGGAAGAGAGAGTCCTGCTATGGGGGATTGGGTCGGGGTCGGTGTGGCCTGTGTTAG'
                          'CTCCACAAGCAACGCATAAAGTGCAGATACTTGGAGCTCAATTGGCCGCGCGCTGGCTAGTGATTTTCTATAGTGGTTTTAAAAGTCGGCCTGGTTTAGA'
                          'CAATTTGGCTCGCCGAAATTGCCGCTACCCGAGTGGCGCAAACGTCCATTCTCATCCCCAGGAATGAAACCGGGCGGACACTGCCATTCCCGCGAGAGAT'
                          'CCTCAGTGAAACTTTGTGGATAAGAATAAATTACCCTGCCACCGGTGATATTTACGAAGGCCTCGACATGCACAATGAGTCAAATGCACGGACCCGGTTA'
                          'CACGCTATAACATTCCGCTACTTCCTTTTGAATACGGTTTGTCGACTGGTTAATACGAGAAAGCCTGTCAAAGTAACCATGAGCGACCGTATGTTCTTTC'
                          'CGCTGTAGCCCCATTGTTCCTTATTAG$',
                          Edges)
            self.assertIn('TTCCCGCGAGAGATCCTCAGTGAAACTTTGTGGATAAGAATAAATTACCCTGCCACCGGTGATATTTACGAAGGCCTCGACATGCACAATGAGTCAAATGC'
                          'ACGGACCCGGTTACACGCTATAACATTCCGCTACTTCCTTTTGAATACGGTTTGTCGACTGGTTAATACGAGAAAGCCTGTCAAAGTAACCATGAGCGACC'
                          'GTATGTTCTTTCCGCTGTAGCCCCATTGTTCCTTATTAG$',Edges)
            self.assertIn('A',Edges)

        def test_ts(self):
            '''Topological Sorting'''
            _,_,adj = create_adjacency( [[4, 5],
                                         [1, 2],
                                         [3, 1],
                                         [3, 2],
                                         [4, 3],
                                         [4, 2]],
                                        back=False)
            for k,v in adj.items():
                if k in v:
                    v.remove(k)

            self.assertEqual([4, 3, 1, 2],create_topological_order(adj))



    main()
