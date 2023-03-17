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
#    along with the program.  If not, see <http://www.gnu.org/licenses/>

'''Common code for alignment problems'''

from sys        import float_info
from unittest   import TestCase, main, skip
from deprecated import deprecated
import numpy as np

from reference_tables import createSimpleDNASubst,  BLOSUM62, PAM250
from helpers          import sign


def get_number_of_coins(money,coins):
    '''
    BA5A 	Find the Minimum Number of Coins Needed to Make Change

    Input: An integer money and an array Coins of positive integers.

    Return: The minimum number of coins with denominations Coins that changes money.

    http://rosalind.info/problems/ba5a/

    We will use Dynamic Programming, and solve the problem for each amount up to and including money
    '''
    def get_next_number_of_coins(m):
        '''
        Find a coin such that we can make change using it
        plus a previously computed value
        '''
        number_of_coins = 1 + money//min(coins)  # This is an an upper bound: we will try to reduce it
        for coin in coins:
            if m>=coin and number[m-coin]+1<number_of_coins:
                number_of_coins = number[m-coin] + 1
        return number_of_coins

    number = np.zeros((money+1),dtype=int)
    for m in range(1,money+1):             # solve for m
        number[m] = get_next_number_of_coins(m)
    return number[money]

def longest_manhattan_path(n,m,down,right):
    '''
    BA5B 	Find the Length of a Longest Path in a Manhattan-like Grid

    Input: Integers n and m, followed by an n*(m+1) matrix Down and an
           (n+1)*m matrix Right. The two matrices are separated by the "-" symbol.

    Return: The length of a longest path from source (0, 0) to sink (n, m)
        in the n*m rectangular grid whose edges are defined by the matrices
        Down and Right.

    http://rosalind.info/problems/ba5b/
    '''
    s = np.zeros((n+1,m+1))

    for i in range(1,n+1):
        s[i][0] = s[i-1,0] + down[i-1,0]

    for j in range(1,m+1):
        s[0][j] = s[0,j-1] + right[0,j-1]

    for i in range(1,n+1):
        for j in range(1,m+1):
            s[i,j] = max(s[i-1,j] + down[i-1,j],
                         s[i,j-1] + right[i,j-1])

    return s[n,m]



def get_longest_common_subsequence(v,w):
    '''
     BA5C 	Find a Longest Common Subsequence of Two Strings

     Input: Two strings.

     Return: A longest common subsequence of these strings.

     http://rosalind.info/problems/ba5a/
    '''
    m            = len(v)
    n            = len(w)
    s            = np.zeros((m+1,n+1))
    predecessors = -1 * np.ones_like(s)
    for i in range(1,m+1):
        for j in range(1,n+1):
            choices           = [s[i,j-1], s[i-1,j],s[i-1,j-1]]
            if v[i-1]==w[j-1]:
                choices[-1] += 1
            index             = np.argmax(choices)
            s[i,j]            = choices[index]
            predecessors[i,j] = index

    i       = m
    j       = n
    matches = []
    while i>0 and j > 0:
        if predecessors[i,j]==0:
            j -= 1
        elif predecessors[i,j]==1:
            i -= 1
        elif predecessors[i,j]==2:
            i -= 1
            j -= 1
            if v[i]==w[j]:
                matches.append(v[i])
        else:
            raise Exception(f'Predecessors[{i},{j}] referenced, but not defined')
    return ''.join(matches[-1::-1])

def get_longest_path(source,sink,graph):
    '''
    BA5D 	Find the Longest Path in a DAG

    Input:
        source An integer representing the source node of a graph
        sink   an integer representing the sink node of the graph
        graph  an edge-weighted graph of triples (a,b,w) indicating that a connects to be with weight w

    Return:
        The length of a longest path in the graph, followed by a longest path.


    http://rosalind.info/problems/ba5d/
    '''

    def create_adjacency_list():
        '''
        Convert graph (list of tuples) to an adjacency list, ignoring weights
        '''
        product  = {}
        for a,b,_ in graph:
            if not a in product:
                product[a] = []
            product[a].append(b)

        return product

    def create_nodes_adjacency_list(relevant):
        '''
        Convert graph (list of tuples) to an adjacency list, ignoring weights

        Parameters:
            relevant   Used to restrict the set of nodes that are to be extracted
        '''
        product = {}
        for a,b,_ in graph:
            if a in relevant and b in relevant:
                if b not in product:
                    product[b] = []
                product[b].append(a)
        return product

    def crop(ordered):
        '''
        We want only the data between source and sink
        '''
        pos1 = 0
        while ordered[pos1] != source: pos1+=1
        pos2 = len(ordered)-1
        while ordered[pos2] != sink: pos2 -= 1
        return ordered[pos1:pos2+1]

    def backtrack(s,predecessors):
        '''
        Used when to construct solution by backtracking through array of distances from source to each node
        '''
        path = [sink]
        while True:
            if path[-1] in predecessors:
                path.append(predecessors[path[-1]])
            else:
                return s[-1],path[-1::-1]

    nodes        = crop(create_topological_order(create_adjacency_list()))
    index          = {nodes[i]:i for i in range(len(nodes))}
    weights        = {(a,b):w for a,b,w in graph if a in nodes and b in nodes}
    backward       = create_nodes_adjacency_list(nodes)
    s              =  np.full((len(nodes)),-float_info.max)
    s[0]           = 0
    assert index[source]==0
    predecessors   = {}
    for i,b in enumerate(nodes):
        if b in backward:
            candidates            = [(a,s[index[a]] + weights[(a,b)]) for a in backward[b]]
            predecessors[b], s[i] = candidates[np.argmax([score for _,score in candidates])]

    return backtrack(s,predecessors)


def get_highest_scoring_alignment(v,w,
                                     weights = BLOSUM62(),
                                     sigma   = 5,
                                     local    = False):
    '''
    BA5E 	Find a Highest-Scoring Alignment of Two Strings
    BA5F 	Find a Highest-Scoring Local Alignment of Two Strings
    Find the highest-scoring alignment between two strings using a scoring matrix.

    Input: v         an amino acid string
           w         another amino acid string
           weights   The matrix to us for scoring matches
           sigma     Indel penalty
           local     If true perform local alignment, otherwise global

    Return: The maximum alignment score of these strings followed by an
            alignment achieving this maximum score. Use the BLOSUM62 scoring matrix
            and indel penalty σ = 5. (If multiple alignments achieving the maximum
            score exist, you may return any one.)
    '''

    m            = len(v)
    n            = len(w)
    s            = np.zeros((m+1,n+1),dtype=int)
    s[0,0]       = 0
    predecessors = -1 * np.ones_like(s,dtype=int)
    if not local:
        for i in range(1,m+1):
            s[i,0] = -i*sigma
            predecessors[i,0] = 1
        for j in range(1,n+1):
            s[0,j] = -j*sigma
            predecessors[0,j] = 0

    for i in range(1,m+1):
        for j in range(1,n+1):
            choices           = [s[i,j-1]   - sigma,
                                 s[i-1,j]   - sigma,
                                 s[i-1,j-1] + weights.get_score(v[i-1],w[j-1])]
            if local:
                choices.append(0)

            predecessors[i,j] = np.argmax(choices)
            s[i,j]            = choices[ predecessors[i,j]]

    i,j       = np.unravel_index(np.argmax(s), s.shape)
    s_max     = s[i,j]
    v1      = []
    w1      = []
    while True:
        if predecessors[i,j]==0:
            j -= 1
            v1.append('-')
            w1.append(w[j])
        elif predecessors[i,j]==1:
            i -= 1
            v1.append(v[i])
            w1.append('-')
        elif predecessors[i,j]==2:
            i -= 1
            j -= 1
            v1.append(v[i])
            w1.append(w[j])
        else:
            return s_max,''.join(v1[-1::-1]),''.join(w1[-1::-1])


def edit(s,t,
         indel_cost   = 1,
         replace_cost = lambda a,b: 1):
    '''
    EDIT 	Edit Distance http://rosalind.info/problems/edit/

    Parameters:
        s            A string
        t            A string
        indel_cost   Cost of insertion/deletion
        replace_cost Cost of replacing a symbol

    Returns: the edit distance between s and t

    '''
    m           = len(s)
    n           = len(t)
    matrix      = np.full((m+1,n+1),np.nan)
    matrix[0,:] = list(range(n+1))
    matrix[:,0] = list(range(m+1))

    for i in range(1,len(s)+1):
        for j in range(1,len(t)+1):
            matrix[i,j] = min(
                matrix[i-1,j]   + indel_cost,
                matrix[i,j-1]   + indel_cost,
                matrix[i-1,j-1] + (0 if s[i-1]==t[j-1] else replace_cost(s[i-1],t[j-1])))

    return matrix[m,n],matrix

def edta(s,t,
         indel_cost   = 1,
         replace_cost = lambda a,b: 1):
    '''
    EDTA Edit Distance Alignment http://rosalind.info/problems/edta/

    Given: Two protein strings s and t (with each string having length at most 1000 aa).

    Return: The edit distance dE(s,t) followed by two augmented strings s' and t' representing an optimal alignment of s and t.

    Most of the work is already done by EDIT; we merely need to backtrack through
    the matrix and to construct the alignment.
    '''
    def backtrack():
        m,n = matrix.shape
        m  -=1
        n  -=1
        s1 = []
        t1  = []
        while m>0 and n>0:
            moves  = [(m-1,n), (m,n-1), (m-1,n-1)]
            scores = [matrix[m-1,n]+indel_cost,
                      matrix[m,n-1]   + indel_cost,
                      matrix[m-1,n-1] + (0 if s[m-1]==t[n-1] else replace_cost(s[m-1],t[n-1]))]
            ss     = [s[m-1],'-',s[m-1]]
            ts     = ['-',t[n-1],t[n-1]]
            index  = np.argmin(scores)
            m,n    = moves[index]
            s1.append(ss[index])
            t1.append(ts[index])
        s1.reverse()
        t1.reverse()
        return ''.join(s1),''.join(t1)

    dist,matrix = edit(s,t,indel_cost,replace_cost)
    s1,t1       = backtrack()
    return (dist,s1,t1)

def  FindHighestScoringFittingAlignment(s,t,
                                        replace_score  = createSimpleDNASubst(),
                                        indel_cost     = 1                                        ):
    '''
    BA5H Find a Highest-Scoring Fitting Alignment of Two Strings
    '''
    def build_matrix():
        matrix = np.zeros((len(s)+1,len(t)+1))
        moves = {}
        for i in range(len(s)+1):
            for j in range(len(t)+1):
                if i==0 and j==0: pass
                elif i==0:
                    matrix[i,j] = 0
                    moves[(i,j)] = (0,0,0,-1)
                elif j==0:
                    matrix[i,j] = 0
                    moves[(i,j)]  =(0,0,-1,0)
                else:
                    scores       = [matrix[i-1,j]   - indel_cost,
                                    matrix[i,j-1]   - indel_cost,
                                    matrix[i-1,j-1] + replace_score[(s[i-1],t[j-1])]]
                    froms        = [(i-1, j,   -1,  0),
                                              (i,   j-1,  0, -1),
                                              (i-1, j-1, -1, -1)]
                    index        = np.argmax(scores)
                    matrix[i,j]  = scores[index]
                    moves[(i,j)] = froms[index]

        return matrix,moves

    def backtrack(matrix,moves):

        score = max([matrix[i,-1] for i in range(len(s)+1)])
        i     = -1
        j     = len(t)
        for k in range(len(s)-1,-1,-1):
            if matrix[k,-1]==score:
                i = k
                break
        s1    = []
        t1    = []
        while i>0 or j>0:
            i,j,di,dj = moves[(i,j)]
            if di==0:
                s1.append('-')
                t1.append(t[j])
            elif dj==0:
                s1.append(s[i])
                t1.append('-')
            else:
                s1.append(s[i])
                t1.append(t[j])

        return score,s1[:-1][::-1],t1[:-1][::-1]

    distances,moves = build_matrix()

    dist,s1,t1 = backtrack(distances,moves)

    return (dist,''.join(s1),''.join(t1))

#-------------------------- Everything above this line has been refactored

def alignUsingLinearSpace(v,w,
                          replace_score = BLOSUM62,
                          indel_cost    = 5):
    '''
    alignUsingLinearSpace

    Align Two Strings Using Linear Space

    Inputs:
        v
        w
            replace_score
            indel_cost
    '''

    pass

# BA5F 	Find a Highest-Scoring Local Alignment of Two Strings
#
# common code

def calculate_scores_for_alignment(s,string1, string2, weights,sigma,
                                   init_predecessors = None):
    predecessors={}
    for i in range(len(string1)+1):
        for j in range(len(string2)+1):
            predecessors[(i,j)]=init_predecessors
            if i>0:
                s_new = s[i-1,j]-sigma
                if s_new>s[i,j]:
                    s[i,j]             = s_new
                    predecessors[(i,j)] = (-1,0,i-1,j)
            if j>0:
                s_new = s[i,j-1]-sigma
                if s_new>s[i,j]:
                    s[i,j]             = s_new
                    predecessors[(i,j)] = (0,-1,i,j-1)
                if i>0:
                    s_new = s[i-1,j-1]+ weights.get_score(string1[i-1],string2[j-1])
                    if s_new>s[i,j]:
                        s[i,j]             = s_new
                        predecessors[(i,j)] = (-1,-1,i-1,j-1)
    return (s,predecessors)

def create_alignment(string1, string2,s_predecessors,
                     i_start = -1,
                     j_start = -1):
    s,predecessors = s_predecessors
    result1        = []
    result2        = []
    i              = len(string1) if i_start==-1 else i_start
    j              = len(string2) if j_start==-1 else j_start
    while i>0 or j>0:
        if predecessors[(i,j)]==None: break
        x,y,i,j=predecessors[(i,j)]
        if x==-1 and y==0:
            result1.append(string1[i])
            result2.append('-')
        elif x==0 and y==-1:
            result1.append('-')
            result2.append(string2[j])
        elif x==-1 and y==-1:
            result1.append(string1[i])
            result2.append(string2[j])

    return (s[len(string1)][len(string2)],\
            ''.join(result1[::-1]),       \
            ''.join(result2[::-1]))

def FindMiddleEdge(s,t,
                   replace_score = BLOSUM62(),
                   indel_cost    = 5):
    '''
    FindMiddleEdge

    BA5K Find a Middle Edge in an Alignment Graph in Linear Space

    Inputs:
        s               an amino acid string
        t               an amino acid string
        replace_score   scoring matrix
        indel_cost      linear indel penalty
        Returns: A middle edge in the alignment graph for s and t, in the form ((i, j) (k, l)),
                 where (i, j) connects to (k, l).
    '''

    def update(j,previous,current,s,t):
        '''
        update

        Calculate scores in current column using values from previous column

        parameters:
            j           Index of column that is being updated
            previous    Data from previous coulumn
            current     Data in current column
            s           First protein string.
            t           Second protein string.
        '''
        current[0] = - j * indel_cost
        for i in range(1,len(current)):
            scores     = [previous[i-1] + replace_score.get_score(s[i-1],t[j-1]),
                          current[i-1]  - indel_cost,
                          previous[i]   - indel_cost]

            best       = np.argmax(scores)
            current[i] = scores[best]

        return current,previous


    def explore(s,t,limit):   # explore
        '''
        Find lengths of all paths from source that end at specified column

        Parameters:

           s           First protein string. This is an explicit parameter
                       because FindMiddleEdge reverse the string during the second call
                       to calculate lengths from sink
           t           Second protein string.
           limit       Last column to be explored
        '''

        column_A = [-(i * indel_cost) for i in range(len(s)+1)]
        column_B = [0 for i in range(len(s)+1)]
        for j in range(1,limit+1):
            scores,previous = update(j,column_A,column_B,s,t) if j%2 ==1 else update(j,column_B,column_A,s,t)
        return scores,previous

    middle_column = len(t)//2
    from_source,_ = explore(s,t,middle_column)
    to_sink,_     = explore(s[::-1],t[::-1],len(t) - middle_column)
    length        = [a+b for (a,b) in zip(from_source,to_sink[::-1])]

    return ((np.argmax(length),   middle_column),
            (np.argmax(length)+1, middle_column+1))

@deprecated
def reverse(chars):
    '''
    reverse

    Input:    a string,
    Return:   new string with characters in reverse order
    '''

    return chars[-1::-1]


@deprecated
def get_highest_scoring_local_alignment(string1,string2,
                                    weights = PAM250(),
                                    sigma = 5):
    '''
    BA5F 	Find a Highest-Scoring Local Alignment of Two Strings

    Input: Two amino acid strings.

    Return: The maximum score of a local alignment of the strings, followed by
    a local alignment of these strings achieving the maximum score. Use the
    PAM250 scoring matrix and indel penalty  5. (If multiple local alignments
    achieving the maximum score exist, you may return any one.)
    '''
    def find_best_substring(s):
        s_right=0
        i_best=-1
        j_best=-1
        predecessor=(0,0,0,0)
        for i in range(len(string1)+1):
            for j in range(len(string2)+1):
                if s[i][j]>s_right:
                    i_best=i
                    j_best=j
                    s_right=s[i][j]
                    predecessor=(sign(len(string1)-i_best),\
                                 sign(len(string2)-j_best),\
                                 i_best,\
                                 j_best)
        if s_right>s[len(string1)][len(string2)]:
            s[len(string1)][len(string2)]=s_right
        else:
            i_best=len(string1)
            j_best=len(string2)
            predecessor=(0,0,0,0)
        return (i_best,j_best,s,predecessor)

    s,predecessors = calculate_scores_for_alignment(
                                            np.zeros((len(string1)+1,len(string2)+1)),
                                            string1,
                                            string2,
                                            weights,
                                            sigma,
                                            (0,0,0,0))

    i_best,j_best,s,predecessor=find_best_substring(s)
    if predecessor!=(0,0,0,0):
        predecessors[(len(string1),len(string2))]=predecessor

    return create_alignment(string1,string2,(s,predecessors),i_best,j_best)




def get_indel_cost(sigma,delta,i,j,di,dj,moves):
    return di,dj,sigma


def build_matrix(s,t,matrix,replace_score=createSimpleDNASubst(),indel_cost=1,get_indel_cost=get_indel_cost):
    def init_indel_costs():
        if isinstance(indel_cost,tuple):
            return indel_cost
        return indel_cost,None

    sigma,delta = init_indel_costs()
    moves = {}

    for i in range(len(s)+1):
        for j in range(len(t)+1):
            if i==0 and j==0: next
            elif i==0:
                di,dj,cost   = get_indel_cost(sigma,delta,i,j,0,-1,moves)
                matrix[i][j] = matrix[i+di][j+dj] - cost
                moves[(i,j)] = (i+di,j+dj,di,dj)
            elif j==0:
                di,dj,cost   = get_indel_cost(sigma,delta,i,j,-1,0,moves)
                matrix[i][j] = matrix[i+di][j+dj] - cost
                moves[(i,j)] = (i+di,j+dj,di,dj)
            else:
                di_down,dj_down,cost_down   = get_indel_cost(sigma,delta,i,j,-1,0,moves)
                di_left,dj_left,cost_left   = get_indel_cost(sigma,delta,i,j,0,-1,moves)
                scores                      = [matrix[i+di_down][j+dj_down]   - cost_down,
                                               matrix[i+di_left][j+dj_left]   - cost_left,
                                               matrix[i-1][j-1] + replace_score.get_score(s[i-1],t[j-1])]
                froms                       = [(i+di_down,j+dj_down,di_down,dj_down),
                                               (i+di_left,j+dj_left,di_left,dj_left),
                                               (i-1,j-1,-1,-1)]
                index                       = np.argmax(scores)
                matrix[i][j]                = scores[index]
                moves[(i,j)]                = froms[index]

    return matrix,moves

def backtrack(s,t,matrix,moves,showPath=False):
    i     = len(s)
    j     = len(t)
    score = matrix[i][j]
    s1    = []
    t1    = []
    if showPath:
        print ('Path')
        print (i,j)
    while i>0 and j>0:
        i,j,di,dj = moves[(i,j)]
        if di==0:
            s1.append('-')
            t1.append(t[j])
        elif dj==0:
            s1.append(s[i])
            t1.append('-')
        else:
            s1.append(s[i])
            t1.append(t[j])
        if showPath:
            print (i,j,di,dj,s1[-1],t1[-1] )
    return score,s1[::-1],t1[::-1]

def align(s,t,
          replace_score= createSimpleDNASubst(),
          indel_cost     = 1,
          build_matrix   = build_matrix,
          backtrack      = backtrack,
          get_indel_cost = get_indel_cost,
          showScores     = False,
          showPath       = False):
    distances = np.zeros((len(s)+1,len(t)+1))
    distances,moves = build_matrix(s,t,distances,
                                   replace_score  = replace_score,
                                   indel_cost     = indel_cost,
                                   get_indel_cost = get_indel_cost)
    if showScores:
        for row in distances:
            print (row)
        for k,v in moves.items():
            print (k,v)
    return backtrack(s,t,distances,moves,showPath=showPath)





def overlap_assignment(v,w,match_bonus=+1,mismatch_cost=2,indel_cost=2):
    def dynamic_programming(v,w):
        distances = np.zeros((len(v)+1,len(w)+1))
        path      = {}
        for i in range(1,len(v)+1):
            for j in range(1,len(w)+1):
                moves           = [(i-1,j),(i,j-1),(i-1,j-1)]
                scores          = [distances[i-1][j]   - indel_cost,
                                   distances[i][j-1]   - indel_cost,
                                   distances[i-1][j-1] + (match_bonus if v[i-1]==w[j-1] else -mismatch_cost)]
                index           = np.argmax(scores)
                distances[i][j] = scores[index]
                path[(i,j)]      = moves[index]

        i        = len(v)
        j        = np.argmax(distances[i])
        distance = distances[i][j]
        v1       = []
        w1       = []
        while i>0 and j>0:
            i1,j1 = path[(i,j)]
            v1.append(v[i1] if i1<i else '-')
            w1.append(w[j1] if j1<j else '-')
            i,j=i1,j1

        return distance,v1[::-1],w1[::-1]

    score,u1,v1=dynamic_programming([vv for vv in v],[ww for ww in w])
    return score,''.join(u1),''.join(v1)

def create_topological_order(graph):
    '''
    BA5N   Find a Topological Ordering of a DAG

    Input: The adjacency list of a graph (with nodes represented by integers).

    Return: A topological ordering of this graph.
    '''
    def get_number_incoming(node):
        '''
        Find number of incoming edges to specified node
        '''
        return sum([node in out for out in mygraph.values()])

    mygraph    = graph.copy() # To avoid destroying original graph - issue #105
    product    = []
    candidates = [node for node in mygraph.keys() if get_number_incoming(node)==0]

    while len(candidates) > 0:
        a = candidates.pop()
        product.append(a)
        if a in mygraph:
            bs = [b for b in mygraph[a]]
            del mygraph[a]
            for b in bs:
                if get_number_incoming(b)==0:
                    candidates.append(b)

    if len(mygraph)>0:
        raise Exception('Input graph is not a DAG')

    return product


def unwind_moves(moves,score,i,j):
    ss = []
    ts = []

    while i>0 and j > 0:
        i,j,s0,t0=moves[(i,j)]
        ss.append(s0)
        ts.append(t0)
    return score,ss[::-1],ts[::-1]

def san_kai(s,t, replace_score=BLOSUM62(),sigma=11,epsilon=1,backtrack=unwind_moves):

    def match(pair,replace_score=replace_score):
        def reverse(pair):
            a,b=pair
            return (b,a)
        return replace_score[pair] if pair in replace_score else replace_score[reverse(pair)]

    lower        = np.zeros((len(s)+1,len(t)+1))
    middle       = np.zeros((len(s)+1,len(t)+1))
    upper        = np.zeros((len(s)+1,len(t)+1))

    moves        = {}
    lower[0][0]  = -float('inf')
    middle[0][0] = 0
    upper[0][0]  = -float('inf')

    for i in range(1,len(s)+1):
        lower[i][0]  = - (sigma + epsilon *(i-1))
        middle[i][0] =  - (sigma + epsilon *(i-1)) #-float('inf')
        upper[i][0]  =  - (sigma + epsilon *(i-1))# -float('inf')
    for j in range(1,len(t)+1):
        lower[0][j]  =  - (sigma + epsilon *(j-1))#-float('inf')
        middle[0][j] =  - (sigma + epsilon *(j-1)) #-float('inf')
        upper[0][j]  = - (sigma + epsilon *(j-1))

    for i in range(1,len(s)+1):
        for j in range(1,len(t)+1):
            lower[i][j]  = max(lower[i-1][j] - epsilon,
                               middle[i-1][j] - sigma)

            upper[i][j]  =  max(upper[i][j-1] - epsilon,
                               middle[i][j-1] - sigma)

            choices      = [lower[i][j],
                            middle[i-1][j-1] + match((s[i-1],t[j-1])),
                            upper[i][j]]
            index        = np.argmax(choices)
            middle[i][j] = choices[index]
            moves[(i,j)] = [(i-1, j,   s[i-1], '-'),     # Comes from lower
                            (i-1, j-1, s[i-1], t[j-1]),  # Comes from middle
                            (i,   j-1, '-',    t[j-1]    # Comes from upper
                             )][index]

    return backtrack(moves,middle[len(s)][len(t)],len(s),len(t))





# FindHighestScoringMultipleSequenceAlignment
#
# BA5M Find a Highest-Scoring Multiple Sequence Alignment

def FindHighestScoringMultipleSequenceAlignment (u,
                                                 v,
                                                 w,
                                                 score=lambda x,y,z: 1 if x==y and y==z and x!='-' else 0):
    def build_matrix():
        s = [[[0 for i in range(len(w)+1)] for j in range(len(v)+1)] for k in range(len(u)+1) ]
        path = {}

        for i in range(1,len(u)+1):
            for j in range(1,len(v)+1):
                for k in range(1,len(w)+1):
                    scores     = [
                        s[i-1][j-1][k-1] + score(u[i-1], v[j-1], w[k-1]),
                        s[i-1][j][k]     + score(u[i-1], '-',    '-'),
                        s[i][j-1][k]     + score('-',    v[j-1], '-'),
                        s[i][j][k-1]     + score('-',    '-',    w[k-1]),
                        s[i][j-1][k-1]   + score('-',    v[j-1], w[k-1]),
                        s[i-1][j][k-1]   + score(u[i-1], '-',    w[k-1]),
                        s[i-1][j-1][k]   + score(u[i-1], v[j-1], '-'),

                    ]
                    moves = [
                        (-1, -1, -1),
                        (-1,  0,  0),
                        ( 0, -1,  0),
                        ( 0,  0, -1),
                        ( 0, -1, -1),
                        (-1,  0, -1),
                        (-1, -1,  0),

                    ]
                    index          = np.argmax(scores)
                    s[i][j][k]     = scores[index]
                    path[(i,j,k)] = moves[index]
        return s,path

    def backtrack(path):
        i  = len(u)
        j  = len(v)
        k  = len(w)
        u1 = []
        v1 = []
        w1 = []
        while i>0 and j>0 and k>0:
            di,dj,dk = path[(i,j,k)]
            i        += di
            j        += dj
            k        += dk
            if dj==0 and dk==0:
                u1.append(u[i])
                v1.append('-')
                w1.append('-')
            elif di==0 and dk==0:
                u1.append('-')
                v1.append(v[j])
                w1.append('-')
            elif di==0 and dj==0:
                u1.append('-')
                v1.append('-')
                w1.append(w[k])
            elif di==0:
                u1.append('-')
                v1.append(v[j])
                w1.append(w[k])
            elif dj==0:
                u1.append(u[i])
                v1.append('-')
                w1.append(w[k])
            elif dk==0:
                u1.append(u[i])
                v1.append(v[j])
                w1.append('-')
            else:
                u1.append(u[i])
                v1.append(v[j])
                w1.append(w[k])
        while i>0:
            i-=1
            u1.append(u[i])
            v1.append('-')
            w1.append('-')
        while j>0:
            j-=1
            u1.append('-')
            v1.append(v[j])
            w1.append('-')
        while k>0:
            k-=1
            u1.append('-')
            v1.append('-')
            w1.append(w[k])
        return u1,v1,w1

    s,path   = build_matrix()
    u1,v1,w1 = backtrack(path)

    return s[len(u)][len(v)][len(w)],''.join(u1[::-1]),''.join(v1[::-1]),''.join(w1[::-1])

# FindMultipleSequenceAlignment
#
#   mult Multiple Alignment
#
# A multiple alignment of a collection of three or more strings is formed by adding gap symbols
# to the strings to produce a collection of augmented strings all having the same length.
#
# A multiple alignment score is obtained by taking the sum of an alignment score
# over all possible pairs of augmented strings. The only difference in scoring the alignment
# of two strings is that two gap symbols may be aligned for a given pair (requiring us to
# specify a score for matched gap symbols).

# Given: A collection of  DNA strings of length at most 10 bp in FASTA format.

# Return: A multiple alignment of the strings having maximum score, where we score matched symbols 0
#         (including matched gap symbols) and all mismatched symbols -1 (thus incorporating a linear gap penalty of 1).
#
#  This function is a straightforward extension of the dynamic programming algorithm of Pevzner and Compeau, with
#  the rectangle replaced by a N-cube.

def FindMultipleSequenceAlignment(
          Strings,
          score = lambda ch:sum(0 if ch[i]==ch[j] else  -1 for i in range(len(ch)) for j in range(i+1,len(ch)))):
    # indices
    #
    # A generator for iterating through positions in hypercube. Iterate along last axis,
    # then wrap around to zero an incrment next axis, etc
    #
    # Parameters:
    #      ns Dimesnion of space (number of strings to be aligned)

    def indices(ns=[len(S)+1 for S in Strings]):
        N  = 1
        for n in ns:
            N *= n
        index_set = [0]*len(ns)
        for _ in range(N-1):
            for j in range(len(ns)-1,-1,-1):
                index_set[j] += 1
                if index_set[j]==ns[j]:
                    index_set[j] = 0
                else:
                    break
            yield tuple(index_set)

    # create_moves
    #
    # Create list of possible moves from current cells to its predecessors, for use in backtracking.
    # Used as we update current cell to record the best past from a predecessor. It will be used later
    # to reconstruct path for backtracking.
    #
    #   Parameters:
    #       m         Number of dimensions
    #       options   Predecssor is either one less along current axis, or the same
    #
    #   Returns:   List of combinations of options, except for all zeroes (which is filtered out).
    #              Each combination is a tuple, so it can be used as an array index

    def create_moves(m,options=[0,-1]):
        # create_raw_moves
        #
        # Create list of moves: each move is a list (not yet a tuple), and the trivial move,
        #  all zeroes, has not yet been filtered
        def create_raw_moves(m):
            return [[o] for o in options] if m==1 else [[o] + c for o in options for c in create_raw_moves(m-1)]
        return [move for move in tuple(create_raw_moves(m)) if len([x for x in move if x!=0])>0]

    # add
    #
    # Used to add a move to a tuple for backtracking
    def add(u,v):
        return tuple([a + b for (a,b) in zip(list(u),list(v))])

    # build_matrix
    #
    # Create hypercube of scores
    #
    # Returns: scores, and also a dict of moves for backtracking
    def build_matrix():

        # calculate_scores
        #
        # For each cell (other than the origin) calculate score using predecessors
        # Allow for case of cell being on the boundary (it doesn't have as many predecessors)
        def calculate_scores(i):

            # get_score
            #
            # Calculate score from a possible move
            #
            # Returns: score of predecssor, plus effect of move
            #          If predecessor outside boundary, return None instead

            def get_score(move):
                predecessor = add(i,move)
                if len([p for p in predecessor if p<0]):  # check to see whther we are on boundary
                    return None
                scorable = []               # Construct vector of characters for scoring
                for j in range(len(move)):
                    scorable.append(Strings[j][predecessor[j]]  if move[j]<0 else '-')
                return s[predecessor] + score(scorable)

            raw_scores = [(get_score(move),move) for move in available_moves]
            return [(score,move) for score,move in raw_scores if score != None] # filter out moves outside hypercube

        s               = np.zeros([len(S)+1 for S in Strings],dtype=int)
        path            = {}
        m               = len(Strings)
        available_moves = create_moves(m)

        for index_set in indices():
            scores_moves     = calculate_scores(index_set)
            scores           = [score for score,_ in scores_moves]
            moves            = [move  for _,move  in scores_moves]
            index_best_score = np.argmax(scores)
            s[index_set]     = scores[index_best_score]
            path[index_set]  = moves[index_best_score]
        return s,path

    # backtrack
    #
    # Backtrack through hypercube to isentify best path from origin
    #
    #    Parameters:
    #        history      Scores and path generated by build_matrix()

    def backtrack(history):
        # reverse
        #
        # Input:   List of characters
        # Returns: String which is the  characters in reverse

        def reverse(S):
            return ''.join([s for s in S[::-1]])

        # is_not_origin
        #
        # Used when we backtrack - verifies that we have not reached origin

        def is_not_origin(position):
            return len([p for p in position if p!=0]) >0

        s,path          = history
        position        = tuple([len(S) for S in Strings])  #Start with vertex opposite origin
        alignment_score = s[position]
        Alignments      = [[] for S in Strings]          # We will constuct alignment starting
                                                         #from last character in each position
        while (is_not_origin(position)):
            move = path[position]
            # Add next character to each string in alignment
            for j in range(len(move)):
                Alignments[j].append('-' if move[j]==0 else Strings[j][position[j]-1])
            position = add(position,move)

        return alignment_score,[reverse(s) for s in Alignments]

    return backtrack(build_matrix())

# sims

# Finding Mutated Motifs

def sims(s,t, match=1,mismatch=-1):
    def get_score(a,b, match=1,mismatch=-1):
        return match if a==b else mismatch

    def score_string(s,t):
        return sum(get_score(a,b) for (a,b) in zip(s,t))

    def dynamic_programming():
        scores      = np.zeros((len(t),len(s)),dtype=int)
        predecessor = {}
        for j in range(len(s)):
            scores[0]  = get_score(s[j],t[0],match=match,mismatch=mismatch)
        for i in range(1,len(t)):
            for j in range(i-1):
                scores[i,j]  = get_score(s[j],t[i],match=match,mismatch=mismatch)
            for j in range(i,len(s)):
                score_diag       = scores[i-1][j-1] + get_score(s[j],t[i],match=match,mismatch=mismatch)
                score_horizontal = scores[i][j-1]   + mismatch
                score_vertical   = scores[i-1][j]   + mismatch
                scores[i,j]      = max(score_diag,
                                       score_horizontal,
                                       score_vertical)

                predecessor[(i,j)] = (i-1, j-1) if scores[i,j] == score_diag       else \
                                     (i,   j-1) if scores[i,j] == score_horizontal else \
                                     (i-1, j)   if scores[i,j] == score_vertical   else \
                                     None
        return scores,predecessor

    def trace_back(scores,predecessor):
        i      = len(scores)-1
        j      = np.argmax(scores[i])

        print (f'Number of maxima={sum(1 for s in scores[i] if s>=scores[i,j])}')
        s_match = [s[j]]
        t_match = [t[i]]
        while (i,j) in predecessor:
            i_pre,j_pre = predecessor[(i,j)]
            if i==i_pre:
                assert j != j_pre
                t_match.append('-')
                s_match.append(s[j_pre])
            elif j==j_pre:
                t_match.append(t[i_pre])
                s_match.append('-')
            else:
                t_match.append(t[i_pre])
                s_match.append(s[j_pre])
            print (f'{(i_pre,j_pre)}=>{(i,j)} s={s_match[-1]} t={t_match[-1]}')
            i,j = i_pre,j_pre

        return (s_match,t_match)

    scores,predecessor = dynamic_programming()

    print (max(scores[-1]))
    print (scores[-1])

    s_match,t_match = trace_back(scores,predecessor)

    return (max(scores[-1]), ''.join(s_match[::-1]), ''.join(t_match[::-1]))


#  ITWV Finding Disjoint Motifs in a Gene

# itwv
#
# Given: A text DNA string s of length at most 10 kbp, followed by a collection of n DNA strings of
#        length at most 10 bp acting as patterns.
#
# Return: An nxn matrix M for which M[j,k]==1 if the jth and kth pattern strings can be interwoven into s and Mj,k=0 otherwise.



def itwv(s,patterns):

    # itwv2
    #
    # Check to see whether two specific strings can be interweaved
    #
    def itwv2(u,v):
        # score
        # Compare one character from s with one from u or v
        def score(a,b):
            return 1 if a==b else 0

        # match
        #
        # Try to match u and v with substrict of s starting at k.

        # We will adapt the regular Manhattan algorithm to use a 2D grid along u and v
        # to keep track of the score when we try to aling successive characters chosen
        # from u and v with s.
        def match(k):
            # update_score
            #
            # Update each position with match-so-far
            # We wil keep track of updated positions by storing the in Closed

            def update_score(i,j):
                if (i,j) in Closed: return
                s0 = s[k+i+j-1]
                u0 = u[i-1]
                v0 = v[j-1]
                S[i,j] = max(S[i-1,j] + score(s0,u0),
                             S[i,j-1] + score(s0,v0))
                Closed.add((i,j))

            S      = np.zeros((len(u)+1,len(v)+1),dtype=int)
            N      = max(len(u),len(v))+1
            Closed = set()
            Closed.add((0,0))
            for n in range(1,N+1):
                for i in range(0,min(n,len(u)+1)):
                    for j in range(0,min(n,len(v)+1)):
                        update_score(i,j)

            return np.amax(S)==len(u)+len(v)

        # Only try to match if 1st character matches
        for i in range(len(s)-len(u)-len(v)+1):
            if s[i]==u[0] or s[i]==v[0]:
                if match(i)>0:
                    return 1
        return 0

    return [[itwv2(u,v)  for v in patterns] for u in patterns]


#  CTEA Counting Optimal Alignments

# ctea
#
# Given: Two protein strings s and t in FASTA format, each of length at most 1000 aa.
#
# Return: The total number of optimal alignments of s and t with respect to edit alignment score, modulo 134,217,727.

def ctea(s,t,indel_cost=1,replace_cost=lambda a,b: 1, mod=134217727):

    # create_adjacency_list
    #
    # Backtrack through matrix and create DAG made from all possible paths
    # Inputs: matrix from computing edit distance
    # Returns:   Adjacency matrix for traversing from last posotion
    #            Leaves of graph

    def create_adjacency_list(matrix):
        Closed       = set()      # Positions in matrix that have been visited
        Open         = []         # Positions in matrix that need to be visited
        Leaves       = []         # Terminal positions
        Adj          = {}         # Adjacency list

        # explore
        # Used to build adjacancy list for graph moving back to origin of matrix
        # It processes first element from Open List, and adds its successors to Open if necessary
        def explore():
            mn = Open.pop(0)
            if mn in Closed: return
            Closed.add(mn)
            m,n = mn
            Adj[(m,n)] = []
            if m>0 and n>0:    # if not a leaf
                moves      = [(m-1,n),(m,n-1),(m-1,n-1)]
                scores     = [matrix[m-1][n]+indel_cost,
                              matrix[m][n-1]+indel_cost,
                              matrix[m-1][n-1] + (0 if s[m-1]==t[n-1] else replace_cost(s[m-1],t[n-1]))]
                index      = np.argmin(scores)
                lowest     = scores[index]
                # Find all possible nodes links to this node, i.e. all that have lowest score
                candidates = [i for i in range(len(scores)) if scores[i]==lowest]
                Adj[(m,n)] = [moves[i] for i in candidates]
                # We need to process any that we haven't already processed
                for i in candidates:
                    if i not in Closed:
                        Open.append((moves[i]))
            else:
                Leaves.append((m,n))

        Open.append((m,n))    # Start with last node from edit
        while len(Open)>0:
            explore()
        return Adj,Leaves

    # invert
    # Invert adcancy list - build list of backward links

    def invert(Adj):
        Inverse = {a:[] for links in Adj.values() for a in links}
        for source,destinations in Adj.items():
            for destination in destinations:
                Inverse[destination].append(source)
        return Inverse

    # count_paths

    # Used the algorithm of https://cs.stackexchange.com/questions/118799/counting-number-of-paths-between-two-vertices-in-a-dag
    # to count paths through DAG
    def count_paths(Adj,Inverse,start=None):
        C  = {ground:1}
        for node in create_topological_order(dict(Inverse)):
            if node == ground: continue
            C[node] = sum([C[x] for x in Adj[node]])

        return C[(m,n)]

    # Build matrix of distances

    _,matrix   = edit(s,t,indel_cost,replace_cost)
    m,n        = len(matrix)-1, len(matrix[0])-1
    Adj,Leaves = create_adjacency_list(matrix)

    ground     = (-1,-1) # Link all leaves to dummy node, so we can just do one pass through tree to count paths
    for leaf in Leaves:
        Adj[leaf] = [ground]

    return  count_paths(Adj,invert(Adj),start=ground) % mod

if __name__=='__main__':

    class Test_5_Alignment(TestCase):

        def test_simple(self):
            '''
            A simple test of the align function, based on https://rosalind.info/problems/glob/
            '''
            score,s1,s2 = align('PLEASANTLY','MEANLY',
                                replace_score = BLOSUM62(),
                                indel_cost    = 5)
            self.assertEqual(8,score)
            self.assertEqual('LEASANTLY',''.join(s1))
            self.assertEqual('MEA--N-LY',''.join(s2))


        def test_ba5a(self):
            '''BA5A 	Find the Minimum Number of Coins Needed to Make Change'''
            self.assertEqual(2, get_number_of_coins(40,[1,5,10,20,25,50]))
            self.assertEqual(338, get_number_of_coins(8074,[24,13,12,7,5,3,1]))

        def test_ba5b_sample(self):
            '''BA5B 	Find the Length of a Longest Path in a Manhattan-like Grid'''
            self.assertEqual(34,\
                             longest_manhattan_path(4,\
                                                    4,\
                                                    np.array([[1, 0, 2, 4, 3],\
                                                     [4, 6, 5, 2, 1],\
                                                     [4, 4, 5, 2, 1],\
                                                     [5, 6, 8, 5, 3]]),\
                                                    np.array([[3, 2, 4, 0],\
                                                     [3, 2, 4, 2],\
                                                     [0, 7, 3, 3],\
                                                     [3, 3, 0, 2],\
                                                     [1, 3, 2, 2]])))
        def test_ba5b_rosalind(self):
            '''BA5B 	Find the Length of a Longest Path in a Manhattan-like Grid'''
            self.assertEqual(84,\
                             longest_manhattan_path(17,\
                                                    9,\
                                                    np.array([[2,3,4,0,3,1,1,1,1,1],
                                                     [4,2,3,4,3,3,0,4,1,1],\
                                                     [4,4,0,1,4,3,2,0,2,2],\
                                                     [4,3,0,3,4,4,3,2,4,4],\
                                                     [0,1,0,1,3,0,3,0,3,4],\
                                                     [3,2,4,4,4,3,1,0,0,0],\
                                                     [3,4,3,1,2,3,0,0,4,0],\
                                                     [2,4,3,4,1,2,0,3,2,0],\
                                                     [1,4,4,1,4,4,3,1,1,4],\
                                                     [3,1,2,2,3,3,0,4,0,0],\
                                                     [0,2,1,4,1,3,1,3,1,0],\
                                                     [1,0,4,0,4,3,3,2,3,1],\
                                                     [2,0,0,4,3,4,0,3,0,0],\
                                                     [4,1,0,4,3,2,1,1,3,1],\
                                                     [2,4,4,3,3,4,0,0,4,3],\
                                                     [1,0,2,3,3,0,4,0,2,0],\
                                                     [3,1,0,3,2,3,2,2,1,4]]),\
                                                    np.array([[1,0,4,4,3,3,1,0,4],\
                                                     [0,2,0,3,3,0,1,2,1],\
                                                     [3,2,3,1,1,4,2,4,4],\
                                                     [1,3,4,4,2,1,1,1,4],\
                                                     [1,4,2,2,3,1,3,2,3],\
                                                     [0,3,1,0,1,0,4,1,4],\
                                                     [1,3,4,4,1,0,3,2,1],\
                                                     [2,3,1,2,3,2,2,2,3],\
                                                     [3,2,1,4,0,2,4,2,4],\
                                                     [4,0,2,0,1,3,1,4,4],\
                                                     [1,3,0,2,2,1,0,3,2],\
                                                     [1,4,0,4,4,1,2,4,2],\
                                                     [0,2,4,3,4,0,3,2,2],\
                                                     [2,3,4,4,0,4,3,0,4],\
                                                     [1,0,4,1,3,3,1,4,2],\
                                                     [4,3,4,3,2,3,2,2,0],\
                                                     [0,1,2,2,4,4,2,4,2],\
                                                     [2,3,1,4,4,3,4,0,3]])))


        def test_ba5c_sample(self):
            '''
            BA5C 	Find a Longest Common Subsequence of Two Strings

            My result doesn't match Rosalind, but it is the same length, and is a common subsequence,
            so is a longest common subsequence.
            '''
            self.assertEqual('AACTTG',
                             get_longest_common_subsequence(
                                 'AACCTTGG',
                                 'ACACTGTGA'))
        @skip('Issue  #117')
        def test_ba5c_rosalind(self):
            '''BA5C 	Find a Longest Common Subsequence of Two Strings'''
            self.assertEqual('AGAATCCATCGACAACTGGATAAGTAGAGGTTGATTGTAACGAGGGTCTCGAATAAGCCTGGCACCGCTACATCGGAGTTGCTGAACAAGATATAGGCCTAGGACAACGCGGGCTATCTGCTTAAATAGAGAAGTGTGATGGCTGGGTCTTGTAGGGAGACCCCCGTGACCAAGTGGTCCGGCATAGAAACCGCCGTGGCTTGAGGTCCGGCAGCGCTTATCCATCCCCATCGCAGAGACTCATTACAAGTTTGCCCCACGCAGGACTTGGTGGCCGGGATACAGGGTGGAACTCGCTCAAAATTTGGACGAAAATTGGTCCGCCAGTGCGACAGAACGGGAGCTCATTTCAAGTGGAGTGAGGAACGTGCACGTTTAAAATAATCCTGTGGTAGACCAATCCCGTGTGCAGGTAAAACTCCGTTAGTCATGTGACCTTTGATCGCCCGTCAATGGCCAACTCGTCTCGGTGTTCCATTGAAAGATCGGAATACCCTATTAGCGCCCATCTCTTACAGGAGAGCGCAGGCATTGAAAGGGAGTTGCAACTTCAGTGCACCCGCATTGA',
                             get_longest_common_subsequence(
                                 'CAAGAACAATTCGCCAATCGTAAACTACACTTTGGAGTCTAAGCTTCAGAGGCTTGGAATCTGCGCCCTTAACGAGGGTCGTCGAATACATTGCTCTGTTGCGACCGGCTATTCTAGTCGCAAGAGTCTTGGGGCTGTAACGGAATGATAATATGGGCCTAAGTGGGACGGAGAGAGCGGCGGGCTGAATGTCTGACCTTGAAATAGAAGAAGATCGATAGATGGCTGGGTCTCTGTAGGGAGACTCCCCTGTTGACGCGAAGTGGTCCGGCATACGACACCACCGCCTTTGTGGCTTCGGACCCCGGGTACCGGACAGCGGCGTTCATCATTCATCCACAAGCATACAAGAAAAGACAGGTGAATGGACGGGTGCATTACAAGTCTTTGCCGCCTACGCGTAGTCAGACTTGGGTGGGGCCGGCGGCCATAGCCAGGTGTGGCAAGTTCTACAAGCTCATATTAATTGTGGACAGACAAAATTGGTCCGACCAGTAGCGACAGAAGCCGGAATAGGAAGCTATCGATTATTTCGAAAGGCAACTGGACGTGAAGATGACACCAACCGTCGCACGTTAATAATAATAAGTCCTATGATTTGGATAGGACCATATTCTCCTTTTGTGCCTCACGCAGGTATAGGGATGAATCTCCTGGTTGAGCTCAGTAGATTAGACCTTGTTAAAGTCGACTACGGTCCCGATCGCAAGTGGAACCACATACCACCTCCGTCTCCCCGGTGATCTGCCATTGGAAAGATGCGAGAATGAGCGCCTAGGATGTTAAGGCGGACCCAGTAACTACATTACACCGCTGTTTCAGAGTCCGCATGGTCAATTCGAACCGGTAGGTCGAGTTTTGTCAACTGTACACCGTGCACTTTCGGGGCTTGTTTTAGGTCGATGTTGCCA',
                                 'AGAATCCTACTCCGACAGACCTGGATAAAGGGTGAGAGCACCGGTTGATTGTGAAGGCGACAGAGAAAGTCCCTCAGAATAAGCCATGGGCAATTCACGCCGTGACCCACTCGGTAGGTTGATCTAGAAACTAAGCACCCTACTAGCTAGCCTAGCGCTACCTCAACGCGGTGTCTACTCTGCGTTAACCGCCCCCGGCCTATCAGGAGAAGTGTCCGAGTGGCTGGGTTCACTTGTCGTAGCGGACGAACCCCAGAACAGAGGTGACCAAAGGTGGTCTTCGGGCACATGTAGGAAAACCGCGGAAAGCAGTGGCCGTTATTGATTTTGAGTGCTCGTTTGTTCTCAGCTGCTTGATGCCATACCCCCTATTCTGCAGAGACCTCAGTATACACAGGGTTATGCACCCACTGCAGGGACGTTGAAGTCGAGCCTGGTGAAAATAAACGATGCGAGATCCGGAACCCTCGCCGTGACACGAAATTTGCCGCATCGGAAATGTAGGCTTCTATGGTTTCACGCCAGTGCAAGAGTGCAGAACGGGCAGGCCTCCATCTTCATATGTGGTAGTGAGGTAATCGTTGCACGATTGTCAAAATAATCACTGGCTAGGTGTTACGACCAAATGACACCGTGTGCAGTGTAAACAGGCTCACCGTGTAAGTCATGTGCTACCTCTTGATTCGCTCCGTCAACTTCGGCCAACTCTGTCATACGGGTGTTTCCCTAGTACTGAATACTAGTTACTCGGTACATACCCTACTAAATATGTCGCCCACTCTGCCTCTGAGCAGGAGAGGCGTTACACCGGGACAGTTGGAAAGGGGAGTTCGCAAACATTGCTAGTGTCAGCCCCCGCAATCTGA'
                             ))

        def test_ba5d_sample(self):
            '''BA5D 	Find the Longest Path in a DAG'''
            n,path = get_longest_path(0,4,
                                [(0,1,7),
                                 (0,2,4),
                                 (2,3,2),
                                 (1,4,1),
                                 (3,4,3)
                                ])
            self.assertEqual(9,n)
            self.assertEqual([0,2,3,4],path)

        def test_ba5d_extra(self):
            '''BA5D 	Find the Longest Path in a DAG'''
            n,path = get_longest_path(0, 44,
                                      [[6, 26, 32],
                                       [10, 39, 30],
                                       [26, 28, 24],
                                       [3, 16, 19],
                                       [10, 35, 35],
                                       [10, 37, 19],
                                       [10, 31, 36],
                                       [10, 33, 32],
                                       [10, 32, 4],
                                       [15, 23, 0],
                                       [15, 21, 0],
                                       [22, 24, 0],
                                       [22, 27, 31],
                                       [1, 3, 36],
                                       [5, 43, 37],
                                       [8, 30, 23],
                                       [19, 34, 11],
                                       [12, 13, 38],
                                       [39, 40, 35],
                                       [12, 15, 29],
                                       [27, 29, 13],
                                       [1, 42, 31],
                                       [24, 25, 2],
                                       [1, 10, 4],
                                       [4, 30, 11],
                                       [13, 35, 17],
                                       [24, 28, 2],
                                       [23, 25, 37],
                                       [31, 43, 7],
                                       [31, 40, 17],
                                       [3, 28, 2],
                                       [5, 12, 39],
                                       [5, 11, 37],
                                       [3, 4, 4],
                                       [2, 31, 23],
                                       [14, 29, 13],
                                       [19, 27, 21],
                                       [27, 36, 20],
                                       [31, 33, 23],
                                       [30, 40, 27],
                                       [28, 42, 29],
                                       [21, 35, 33],
                                       [21, 37, 5],
                                       [20, 37, 24],
                                       [2, 9, 38],
                                       [0, 14, 19],
                                       [4, 20, 0],
                                       [1, 41, 8],
                                       [8, 14, 28],
                                       [19, 20, 13],
                                       [4, 43, 3],
                                       [14, 31, 25],
                                       [14, 30, 22],
                                       [13, 41, 19],
                                       [13, 40, 32],
                                       [14, 35, 10],
                                       [10, 11, 5],
                                       [14, 38, 23],
                                       [2, 23, 9],
                                       [2, 25, 1],
                                       [24, 40, 37],
                                       [12, 38, 38],
                                       [20, 23, 34],
                                       [20, 21, 29],
                                       [12, 30, 10],
                                       [12, 37, 12],
                                       [29, 44, 30],
                                       [33, 35, 15],
                                       [33, 37, 22],
                                       [0, 36, 8],
                                       [37, 38, 17],
                                       [10, 29, 13],
                                       [17, 44, 11],
                                       [6, 14, 5],
                                       [10, 22, 8],
                                       [22, 37, 19],
                                       [22, 34, 3],
                                       [32, 43, 4],
                                       [15, 36, 28],
                                       [11, 35, 20],
                                       [2, 16, 27],
                                       [7, 10, 22],
                                       [11, 31, 19],
                                       [16, 41, 24],
                                       [15, 30, 25],
                                       [32, 37, 29],
                                       [0, 27, 9],
                                       [0, 28, 7],
                                       [32, 38, 0],
                                       [12, 43, 5],
                                       [22, 35, 37],
                                       [24, 30, 7],
                                       [24, 32, 19],
                                       [24, 38, 38] ])
            self.assertEqual(62,n)
            self.assertEqual([0, 14, 29, 44],path)


        def test_ba5d_rosalind(self):
            n,path = get_longest_path(11,18,
                                [(21,33,29),
                                 (5,19,8),
                                 (11,17,25),
                                 (10,35,33),
                                 (14,27,19),
                                 (6,27,20),
                                 (6,20,32),
                                 (14,23,4),
                                 (7,8,29),
                                 (12,18,9),
                                 (20,27,16),
                                 (15,35,39),
                                 (1,27,2),
                                 (9,19,37),
                                 (12,14,20),
                                 (30,34,9),
                                 (2,9,5),
                                 (24,31,37),
                                 (8,18,33),
                                 (8,13,20),
                                 (17,23,20),
                                 (6,34,13),
                                 (9,29,20),
                                 (16,29,12),
                                 (1,22,8),
                                 (2,32,27),
                                 (13,26,25),
                                 (3,21,2),
                                 (13,20,22),
                                 (3,23,20),
                                 (8,9,13),
                                 (12,16,30),
                                 (3,29,22),
                                 (1,3,2),
                                 (19,35,7),
                                 (14,19,36),
                                 (10,21,38),
                                 (19,30,38),
                                 (23,32,3),
                                 (5,6,19),
                                 (9,33,2),
                                 (6,10,17),
                                 (21,24,8),
                                 (9,15,36),
                                 (10,27,17),
                                 (25,33,20),
                                 (9,11,1),
                                 (0,10,39),
                                 (6,31,39),
                                 (6,16,24),
                                 (0,6,34),
                                 (14,15,30),
                                 (1,7,0),
                                 (0,22,3),
                                 (4,30,3),
                                 (4,15,20),
                                 (0,21,18),
                                 (0,26,28),
                                 (22,27,20),
                                 (3,32,32),
                                 (3,33,10),
                                 (5,28,13),
                                 (17,32,33),
                                 (7,13,25),
                                 (17,18,37),
                                 (8,31,33),
                                 (7,18,20),
                                 (20,35,20)
                                 ])

            self.assertEqual(62,n)
            self.assertEqual([11,17,18],path)


        def test_ba5e_sample(self):
            ''' BA5E 	Find a Highest-Scoring Alignment of Two Strings'''
            score,s1,s2 = get_highest_scoring_alignment('PLEASANTLY','MEANLY')
            self.assertEqual(8,score)
            self.assertEqual('PLEASANTLY',s1)
            self.assertEqual('-MEA--N-LY',s2)

        def test_ba5e_rosalind(self):
            ''' BA5E 	Find a Highest-Scoring Alignment of Two Strings'''
            score,s1,s2 = get_highest_scoring_alignment('HSQLDPPHCYCVPPWSWLDYACRLNDCSCHCMGKFWKKDGTTPYPLFTFEGWKQVTFVCDCFKYITINAEPQSTRSMPKIGMGELPHFINDFHNGGRDKICRTHQICALYLVSGPQGKCWSGQDNTKVKSFHRFMFLTYEIMYSSPWVNGYFREFSACPSPEVQYSDQLNEPNFTNGLHIMIMREWKRWNVQDYLSMICDQCCGWCSCGKPVIDPAVHLYSVWCRCTGMDKEMLQTDRWAFMWHCCSFRNERTEHKFFWFGSDVAVRPRLKGWMNTACNRTKNMCHCATHWCNEMGPRIAVIGGFQECLDHWGHHGRLLEFGISRLGQAKSASCMICMKCGGVLGGKVFTKKVNCVFCGDPYCYLYAKCCAEPTNTGTPKLPAGRLNNVNESYYGGSEITDFIEALSCKTMNKCAFNDAWTKHIAESAMGDEVAKFYSRGSPFNINQERFCWPRHTNDDWDIHMFPINHDTYHRACGWTFNHQKVGIRWAAASSPVQETCTDIPHEFSVPASQYIMIRLRIHCKYPQYWPYGKDMKRGTAIRMVAQNKHCPDLADRQRSEYNVDLMSYHVPFEYRNQRGFATDDNSWKPNDLHIERHDDRQCMPHITTWIAQSPPNSCGTCNAMFWAPRRGYSSNKVWHNPDDHHNQDWCHQEGAMLAQQGELHDCIRSDWRGNPDPAGFPYDNQLSHSQHIRDDWLEEFEEMAYYTQFLQWHRGYIEDGIKQECIQDGPNENIAPREVLKIEHFQMWGETKEHWYHIYYFTEMISRDRVAVEYWKFPFKISTYMVADIGEFCKTHQYGDSCPHQDPLDIKFSQMLKDEKQKKGICKQWINVCHLPMDCIFYRNGTTQQFNCHWCTCNWMFRGNIVATTMWWEDNHNYDLWDLDDCAVQRLFVYVTTGHMVQHYDCFDIECSDDSEDELTRLCFSGLKSMRSEYTQHSCKPMWFYLCNDCCKRHRFCTMEAYHSVMQGHCMTMGTTRQVYELVRVPMSQVYEHWIAHVREHLLVATTIFCAETFQGKQYQRQHHQTAAFCLNHNEEDGSAFIENQYWYLCKNNIWCYYQAAEKFQEKVPELLNRARPGIKNDPCMEMPLFSMKIMHYNYMCCMACYRKRCRPTFAMFVDPCRIVLESEQVNWMFLFLCGLYNLRYTSVTPVNCWGGLYECRYDRHWPYIRKEEGRDHQMYYNQLLRTNTADGDEAHWCAFLEKVPLLEMYKNVKEAGPNSHKKVYDCAMSPQFCQEKQVAMTYQGNRMSQKCAETDYGAHLQGPMDNLGNLHEMARIFRLSNSIPPKSPTIMWWMDAARIINQRHKRSDKLRWGKCGGPIMAFINEIRFEYYPVNICKRIASTLVERYTFSFYKPYCVRWVTPQKTVVCEISCAGHGCCSPWHKHWTMCYSVCCLRQQWTVTCCKPLNAIYGTSPQSGIFHIAMVYNKVFQDQLMYYSKITHAEKTSTTAEMSEDYWPRIIFLKSCLDWLINDDMIELYKLYQMWRCTASQMKDDLVSDEIDDWVSLGRQMGHIDGKYATSCRIAAHMAVKECEACCHFCEIVPCKAHGQTNFPAGGLQCCYHRNMVQRHSMCFNQENHWDMANFIVQFVTEICTYGFMEVKMYWEGHKMHVKPAYQCIWSQWEMFCRANFGQWWKHGRQRGYRSNMKTSYPGMWFYAWVCAGPWWLLHLEFGRCNNCLSQSFCQLRNDWTPCGLCPAVTGDMDWVTEGAKYDRDTLQKPFMASFFEYIHNVIYKDFAQKPCPEVEQHERCFRKYARSKRCNFFGMWSLRISGCHWDENKIYIDIGEHHRHLSMAMSPTKKKIYMRIFLSPHARKSAMNVADAIYHANHGPMEWLRQSFFRHHGHDSKNQPCHDTEKLDQNSFFDDKNRSFNTWHGCRTTAAKFNWYMFAVHPQKLYVTAEHFFVGQIFLVMFSIMLFVYMAHVFACHPNASHSGLVPWAYFNMYFRLEKVDCQCAFYLFFDLACRIKHAQMKCPDQICTQPGQLTLLKYRAGMINDDQSCEFTHMWRNYTYCVTHTPEHNQYKRHNGESYCRTVTDFWRWSRHNQFGYMDDKWGMSNRVQYNDSHGAVRGCLVRMGYIYDVLCFYGFFNHKFACSWCSDFHHTWDTKIPLMFCDWEGFNEVSHRHVCIRSEVFTISYWHASLRAQIRGALLNRRTRIPPRIRYEMQGCCVNNMRWVPFSTLNTVMHNQYQEECNTPCCGDSYQSPAGQDCGIVLGAKEFRIMRINYTWELDRCQLMKGTMKLDVAHKYWKKMPMPHEVPMTRRCGFQFSSPRTMWLNIGLAGTGLYKLCHNCCPKIRTVEAAGWGCWLQMCDAKSFCKDTKWKMARFHNIRNRHIKYPGKCFMHNMLLWPFAFPQWYDANYYWFYPDTGPYIFPCGGLFYFDGFTAFCAQLNGYRMQTHAGLFFHHAKHRYCHSVSREHLWAMVGAGVEQPRMNGTRSYYKQSWFDWSENSCGYKQDQILLRNNGTAWRRPEHHVMYIGDAFYLENKMRYQNIHCPESQHATQIWAAYAFQPCSVSNSTSTSKRFTRCMYLSAYAITQTYYDYWKGALNNFNRGDAPQYNHPWEFSNWRTVNPWAMWDQTREGKCMFLAEACNWFVSGSHQDWFIPDLANTWGCSKEKKYMHNPKIVTYNEMLAVAFRQEQRYSIHSPCKTKGCIGHKCMNWKIIKWMTLCGDHCGMVQIYFKLCKPGVPINRCFTAILHEIPPLYTGGAVVKMNHYSRTGIWHDGPQYCYWKIDIHYFKFNAWNPIATWDWVRWLASSAWSHGHSKIQRVVMHREGWMVAWYIKWVFKGKDEQDNWEQPHVAMASEVAWAGLKGFCPKGTIPGSHPKILYKRTEKGPDWMYYRAMTQGGKGGVIQGSWGVEHDSHANFIKCVTLALSVTKWQHRRHEAKKQVAEPSMDMVYWDYEVWCPCSPIAPRSKPAIMSLNCWNDRHPGAIQHQLGKLFSTQGAGFWFVIHKCGWAADRRQIHSLNMYRFQPDLCIPKVRGANPLYMADAWWPCAACVKYFSRLTFLYRMNFLKWIHGAHNVIYCDHTWHNCMRASRTDMDAHARAKACDLMPHHVDGRRWPIDSHCPNGKSPAPAIYQAWDLECKQLAPKESWKWSIVRFHIASASEMTWHMWYCVLCAMTHHPDAVDGANLWDSPGEMNIMHWEFEEFPHEWARCVGPDKHWPKMFERGFPMEGNWKCDRKEWQKIERLFITRLCINRPSKNVSYQEYNDICDMKTSIRYKYDQPSQPERAGYLDAMYVYECGSDGLAHNSIEGEWPRVNQAFFARQSLYHNQWRPEMEKYAICLHKFKNWHTYKDQITRPQMVATNDSQFNDWNMMELQIEDVLAGQINDPMPSFKQVMPMPYCKNIDTKEVMYMETEPVYRFYDHHMVHIQNEYTAMDNCESVQNRKSWLVFRMFMIDKEGQKVQRGCYHPKRTIEHKNVDSFIMEVNSCEFGENHDATDTLKRDMPKKNDCGVHEDQVEQVSRHGLPYYQNMCKNHTPVWVNVMIWIAPCVIFMKTGFHSNIHPWFWVCYCSYTYQIGPRGMIWLPITHSCQYQTGQQQWFPRCSLAPEYTQAPEYLGYGPLFDAKMEFHFAWRCMMEIWGAPVKPCYITSNMKACCCIKLNGGYISSFISHVVMKDPLKMDMPQFEMEADNDDMDNNTKCETLSLFCCMEYTWQPCTTFAKLRKECCNGDGKTVSEHQVFWPFPNACQHSPYKTCVYAITLNMVTRLASMPISSANIWKGMMHIPSLDTLPWEAPGVIMAKWFPCKTWINWYIYCVATIGEPDPSQPECSVPSIDKMIKKMHERPAIQPHRTYTAHPTADVYRGMTWFADSPVHYLKTGMPKGHYYFVRTHSEGGSLESTMKNWANCYDDQSAVYCHACTKKTVYGDFNMYKKKHQEHHGYTRLGIRPPIYEHPAWMLMHYRKVMGQPEACLETDCAPIASQIMTCPAKTPYELWNRRRSYGMMPAHWCLFLTCEESAVKDRPDSHTKDMWAGYICKRGAMHPVKDSDRTNNKQISCMNCEWTYLRMCRETYSRQNNKGASQPMSGMREHICRCWMMHAHHAGYECPQTGMFLNFCNFYNCHDTVLMENIYHETVQNGCLLKHLF',
                                                               'HSQCHRMWVDPCQGKNVFHCYPWSWLGNDNDGIRESSCHCNGKGWKQVEFVCDLFSYITKAPATRSMPKIGWIRRYFNFGELPHFAYDFHNGGRDKKCRTHQITALYLVSGPQGKCWGNDKVKGRIGTAEFHRFMFLTFWVNGYFREFSKNPCPSFEVQYHDQLNEPNFTNGLHIKIMREWKVTWNVQDYLSMITCQCIRYDDQTWPQTCGWCSCGKPVIQDWMLQTDRWAFMWHCCSFREERTEHKFPDWWFGSDVAVTACNRQSNMWKCATSWCNEMPGERIAYIRGVMCISSECFLEYSRLGQWKVVNHNMKSASCMICMKCGGVLGGKVMEYYKKVNCVFCHVKDPQCYGMEPKLPAGNNVNESYYGGSEITSCKTMNKCAFRPLDDAWTKHIASSAMGDEVAWISYEHFGDYFPPARPTMSRGSNNNINQERFCDDWDIIHNCFESDLSINHDTYHRACGQKPGIRWAAASSTIYIMIRLRIHCKYPQYWPYGRMMIDMKRGTAIRMVAQDRQRSEYNMSYHYRNQRQDKTDDNSRKPNDLHIEDHDDRMPTNTPYCMPHDYYNKDFIKICIAMFEASSPEHKVWHTFLTMKNPWEDVPLSLDDHHNQAWCHQEGAMLAQQYELHDCIRSDWRGNEDPSHSQQNMVNLEEFREMAYYTQNLQWPRGYIEITIKQHCAPQEVLKIEYDIQTFFQMWGQVKHWYEMRIQVTMIYYFTEMIVGYWKFPKDDIECQAGKTPPSTYMVADGCTPGQYGDSCMHDVQFSQMKGIFMGPTMLAVKQWICWCHWPMDCIFRNGTTQQTPLCCHWCCNWMFRGFIVATTMGWEDNHNYDLWDLDDCAVQRLFVYVTTGHMVQHYDCFDIECSDDELTRLCFSGLKCKGCWMRSEYTQHCCKPMGFYLQYLSQTCVHSVMQGRVPMSHVRERLLVATTIMFCTWPTEKEASKQYQRQHWAQTAAFCTALYRECNENQVRYHCGMLMINNIWCYYYYEKFQEKVPELLARARIKNQGPSYRSPCMEMPLFSMKIMSQHGNDKRCRPTKAMCVDSCAIVLENEQVNWMFHCFKVGLYNLRYDSVTLQINCWGGLKERMRYIRKRNAPCLAKSSCLWSKKTYNQLLRTNTADGDELERSYGQQKSHKKVCMCTNRDCAGSAQDHQEKQVAMTYQWNRMKCAEINPDMYPKKYAAHLGPMDNLANLHEMARINSIPPKSPTPMVWKDARRIINQRAKYSDKLRWGKEGEPIMAFWCHRQRFEYYVNQCCNHLYVERYTFIWATNQYTVICEISCALKRIHWPWGQHKVCYLRQQWGVTCCKPLNAIYSAQSGIFHFSGFQWYAFQDQFMYYLKIEHAEKTSTSLYWVQVYVAEMSRHDYWPSYKTVHVGVNYDLINKDMIELYKLYRSFTDNMWRCCASQMKDDDGCLARQHIPWFIHPCIVGKYAAAMAVKEACKYHMDAMMKAEFWACSPINQYYKSNCCYHRNMVQRHSHCFNQENHWDMANFIVQFVTEIMTYLFMEYKMYWEGHWCGQQGTTFCIWSQWEMFCRANFLRDDQQWTWHLWNHGRQYGYRSNMWLIWKPAYRTSYPGMWQIAWDDMSHKCWMADGGVAGPWQLVIHLEKKFWFGRCNNCLQQSFCQCRLFTFLTNFMNAFNWWTPCGLCPATNGPWGCMDWVTEGAFRIDFMASFFEYLHNVIYKDFASEYQNENPIPEAEQHEYARSKRQHTDERKIYIDGFSCSKVTKKKIYMRIFKSAMNVAKAIYDPNMVWGLRGSFCRTHGHDSKNQPCDTEKLDQCGTITARAFLYHDFTWHGCRTTAAKFNWYAFSVDDVTAEWFFFGQIFRVMFSIMLFVYMAIASHSGLVPWAYFNMYFRCAFYLFDMEWKFDLACWIKHAQMKCCKTFMQICTQQLTLLKYRAGMSNPLRMGIFRMWRNYCHWEQSCVTHTPQVKSKLKSNNGESYCRTVTDNVVQWLSDDWRWSRHNHIWKWGMSDRVQYTSHGRMPYIYDVLCFYGFFECMMENGSEHKFLICSDFHHLMFCDWRNEVSHRHVYMLRGFFTISYAHASWRAQLSRRTRIPPRIRYDLHDIWRSRVPHNKINMRWFMTLNTVMHNQYQEECNTPCSAQHPAGVDCGIVLGADEFRIMRINYTWELDRCQGTMVAENYKVQKYWKKMPYPHLEYMTRRCGFQQSSSPRLAGTGLYKFPTEQCHNCLAGMVGYLQWLQMCDAKFCDDTKWKMARFHNGRNRHQKKMKYPGNCWYLLQHNAVTQWYDANYEHEDKNLCFYYDEGPYIFPCGGYCRFQSCVWYDGFYAFCKVKMMEYWCGMGATCLRMHPPNEDDLKHRYCHSTSREHFWIHMVNGTRQTNHCYKQPMRILGQNWFDSENSCGSKQDQILLWNAWWLHIYLQVPEVHVMDIGSKYNMDNVHCPESQEATQIMNKILQQYEAAYSKRKYFVTVCTSCMYLSAYAITQTQFVCYDYWKGACGNDNNFNRPIQRPEEFSNWATVNPWAMWDREGKCMFFCSGSHQDWFKMHNPKIVTYNKMLAVAFNRNCDMANAFGEEYWNQPKRRNQVDWPTKGCIGHKCMNWKYIKWGMVQIYFCSCKLDVFYLCDKHLHINRCFTAIPRSAHKIPVIYTGGAVNKMNHYSRTICYWKHDIHYFKFNAWWPIATSDWVRAMTLPPPALASQAHGHSKIQRVVMHRCGWMVAWYIKWVFKGKDEQDCIEQPAWAGLKGFCPKGTINVGMKNGSHPKILYKRTEKGPYPGWWHKIYWMYYRAMKRGGKGGVIQGPKMVGVEHGSHANFIKFQGYFCCVPLALSVMQDIKFSKKKLTLFPLTTQHRRHEAKKQVILYEQYIPREPSPMYEIAPRSLWYYVPAIMSVHWLVLPNFNCWNDRHPGAIQHQNGEDFSEAESGTQGAVFKFVIHKCGWAADRSLNMYRFQPDLWAMAACPLYMWDAWWPARRLTHLPRMNFLKWIHGTHITDKANYNQWCQHTKGWLWKTYNIQVSNWMRASRDAHKWRAKACDMPHSVDGRRWPIDSHCPNGHSPAVAIYQAWETIHPLILECKQLAPKESWKWSIVRFQCDHGPTINIASASEMTWTMWYCVLCAMTGALWDSPGEMNIMHGAFSDAEREGWLHPKQNVCVKHEMARCVGPDKHWPKMFERGAPMEGNWKCKRKEDQKIERLFITRLCINRPSKNVSDHKWFQEYNDICDMTTFIRYKYDQPMSHVTGYLDEMYVYEYGSDGLAHNSIEGEWPRVNQAFFARQSLYHNQWRPKERWPYFKYKDQPTRPQFAPAMEGHATRDSQFNDWNMMNLQIEDVLAGQINDALDVPSRVTWQHKQVMAMPYCKKEPDFAYWVMYRETECRYGDWHHMVHIQNEYTAMDCGSVQRYSWLVAAQSMIDTEGQKVQRGCYHPMKRTIEHKNVVNSDEFEFGENHDATSTLKRDMPKKVHEDQVEQVSRSGLPNYLNMYKNHTPWQKSQFWVNVMIWIAPCVIVFTQCYEPIIKTFFHSNIHPWFWVCYCSYTYQIGPGIKGMWWLPITHSCQSRSHIHSLQQQQWFPRCSLFPEYTQAPLGYGPLFDAKMEFHFAWRCWSAKVKPIYITSNMKCCIKLISHVVMGQDPLKMDLPVFEMEADFRPCQPPDDMDNNTKCKTLSLECCMEYTWCPCTTFAKLRKECCNGDGVFWPQLFRAYKTQVYNMVTRLMDHLVVYSSANIWKGMMQIPSNDNLWEAPGWAKWFPWTDWKKTWINWYIYCVAPHFPRGSIGEPDPSQPPSIDLMRPACQRRTYTAHYRGGTWFADHGDHYNNNPPVHYLKTGMVRYSTMKNWANCYDDQSAVTVYEFRHISKDFNMYKKKHQEQHGINLNKTRLGIRPPIYEHDAWMLMHKRKVPGSACKCWCTETCHEIGHCGWLVKPIASQIMKPVTTLEDTPYELWNRRRYGGYGWMPAHWCLFCEESQTWRDPVKRPSSHTKDMWAGYISPYVDSDRTNLLHCKQLSCMNCEWTYLRMCRHGATNKCNDQSCKVTYYDEATGPHREGEAQLTCRCCMMHAHHAGYENFCNRVCQYSYNVGHDTVLMENIYHFDITVQQMHPRHGCLLKHLF')
            self.assertEqual(9842,score)


        def test_ba5f_sample(self):
            '''BA5F 	Find a Highest-Scoring Local Alignment of Two Strings'''
            score,s1,s2 = get_highest_scoring_alignment('MEANLY','PENALTY',
                                                               weights = PAM250(),
                                                               local   = True)
            self.assertEqual(15,score)
            self.assertEqual('EANL-Y',s1)
            self.assertEqual('ENALTY',s2)

        def test_ba5f_rosalind(self):
            '''BA5F 	Find a Highest-Scoring Local Alignment of Two Strings'''
            score,s1,s2 = get_highest_scoring_alignment(\
                'AMTAFRYRQGNPRYVKHFAYEIRLSHIWLLTQMPWEFVMGIKMPEDVFQHWRVYSVCTAEPMRSDETYEQKPKPMAKWSGMTIMYQAGIIRQPPRGDRGVSDRNYSQCGKQNQAQLDNNPTWTKYEIEWRVQILPPGAGVFEGDNGQNQCLCPNWAWEQPCQWGALHSNEQYPNRIHLWAPMSKLHIKIEKSSYNRNAQFPNRCMYECEFPSYREQVDSCHYENVQIAFTIFSGAEQKRKFCSCHFWSNFIDQAVFSTGLIPWCYRRDDHSAFFMPNWNKQYKHPQLQFRVAGEGTQCRPFYTREMFTKVSAWRIAGRFAGPYERHHDAHLELWYQHHKVRTGQQLGIIWNNRDKTRNPCPFSAYYNKLPWWKINQNAFYNCLQNIAHSTHDETHEFNPVKCIDWLQGTMVPTECKKGFVHEKCECYRNPGPPLHDMYHQMEDIFGVRFDCLTGWKHLSDYNPCQERRNINDFYIFAYEIAPAVKNLVLSPQPLADATKKCAFNYTPLDQSPVVIACKWYIHQPICMLLIVLICAMDKYNAHMIVIRTTEGQQPMHACRMTEGPGMCMKEPLVTFTLPAQWQWPNHEFKYVYMYVLNYHLSQYTYTDEGHAGGQHYSFNVAVDVGMAWGHNRCYCQPACYSQQETQTRTIDYEKWQYMKHQAFKWGLWFCEQERHAWFKGQNRCEMFTAKMTRMGADSNLDQYKLMLAQNYEEQWEQPIMECGMSEIIEIDPPYRSELIFTFWPFCTYSPWQNLIKCRCNNVIEEMDQCVPLTFIGFGVKQAGGIQAWAFYKEEWTSTYYLMCQCMKSDKAQYPYEIILFWMQPMDTGEQEPPQQNMWIFLPHSWFFDWCCNAPWSEICSSRHDHGQCQDAFYPCELFTVFDDIFTAEPVVCSCFYDDPM',\
                'WQEKAVDGTVPSRHQYREKEDRQGNEIGKEFRRGPQVCEYSCNSHSCGWMPIFCIVCMSYVAFYCGLEYPMSRKTAKSQFIEWCDWFCFNHWTNWAPLSIVRTSVAFAVWGHCWYPCGGVCKTNRCKDDFCGRWRKALFAEGPRDWKCCKNDLQNWNPQYSQGTRNTKRMVATTNQTMIEWKQSHIFETWLFCHVIIEYNWSAFWMWMNRNEAFNSIIKSGYPKLLLTQYPLSQGSTPIVKPLIRRDQGKFWAWAQMWWFREPTNIPTADYCHSWWQSRADLQNDRDMGPEADASFYVEFWYWVRCAARTYGQQLGIIWNNRLKTRNPCPYSADGIQNKENYVFWWKNMCTKSHIAFYYCLQNVAHYTHDVTAEFNPVKCIDWLQGHMVLSSWFKYNTECKKLFVHEKCECYRMFCGVVEDIFGVRFHTGWKHLSTAKPVPHVCVYNPSVQERRNINDFYIFYEIAPAVKNLVLSAQPLHDYTKKCAFNYTPITITRIISTRNQIIWAHVVIACQFYSPHQMLLIELAMDKYCADMNVRRSTEGHQPMHACRSTFGPGMAAKEPLVTFTLVAFWQWPNHEFQYVYMYTEDKIIQIGPHLSNGCEMVEYCVDCYAKRPCYRAYSAEAQYWRMITEAEDYSYKTRNAIAATATVRGQYCHPFRWLGIVWMAHHDCFFANECGTICIPQMAEMRPPETTPYEIDIIFMMFWKEHMSTTILDVVGMYRPATFSHWHDAHHQCEPYLTPLMCQSKLVFDAAFTQVGVKGVWYHTEKLELMAGFNHMKFKKEEAQQSCFYWFQDCPDYDPPDAVRKTDEKHIRAHGEIWWLMRYYCMYHILHIASRHEWMHLRWDQACTNPGYELFEFIPWVLRRYVVYDKIRYNYSYRNSASMEFV',
                weights = PAM250(),
                local   = True)
            self.assertEqual(1062,score)
            self.maxDiff=None
            # self.assertEqual('YQAGIIRQPPRGD-RGVSDRNYSQCGKQ-NQ-AQLDNNPTWTKYEIEWRVQI-LPPGAGVFEGDNGQNQCLCPNW--A-W-EQPCQW----GALHS-NEQYPNRIHLWAPMSKLHIKIEKSSYN-RNAQ-FPNRCMYECE-FPSY-REQVDSCHYENVQIAF-TIFSGAEQKRKFCSCHFWSNFIDQAVFSTGLI-PWCYRRDDHSAFFMPNWNKQ--YKHPQLQFRVAGEGTQCRPFYTREMFTKVSAWRIAGRFAGPYERHHDAHLELWY-QHHKVRT-GQQLGIIWNNRDKTRNPCPFSAY-Y-NK--LP-WWK-I-NQ-N-AFYNCLQNIAHSTHDETHEFNPVKCIDWLQGTMV-P------TECKKGFVHEKCECYRNPGPPLHDMYHQMEDIFGVRFDCLTGWKHLS------D---YNPC-QERRNINDFYIFAYEIAPAVKNLVLSPQPLADATKKCAFNYTPLDQSPVVIACK---WYIHQPI-CMLL----IVLIC-AMDKYNAHMIVIRTTEGQQPMHACRMTEGPGMCMKEPLVTFTLPAQWQWPNHEFKYVYMYVLNYHLSQYTYTDEGHAGGQHYSFNVAVDVGMAWGHNRCYCQPACYSQQETQTRTIDYEKWQYMKHQAFKWGLWFCEQER-HA--WFKGQNRCEMFTAKMTRMGADSNLDQYKLMLAQNYEEQWEQPIMECGMSEIIEIDPPYRSELIFTFWPFCTYSPWQNLIKCRCNNVIEEMDQCVP-LTF-IGFGVKQAGGIQA-WAFYKE--EWTSTYYLMCQCMKSDKAQYPYEIILFWMQ--P-MDTGE--QEPPQQNMWIFLPHSWFFDWCCNAPWSEICSSRHD--H---GQ-CQDAFYPCELFTVF',s1)
            # self.assertEqual('Y-P-MSRKTAKSQFIEWCDW-F--CFNHWTNWAPLSIVRTSVAFAV-W-GHCWYPCG-GVCKTNRCKDD-FCGRWRKALFAEGPRDWKCCKNDLQNWNPQYSQGTR--NTK-RMVATTNQTMIEWKQSHIFETW-LF-CHVIIEYNWSAF-W-MWMNRNEAFNSIIKSGYPKLLL-T-QY-P-L-SQG--STPIVKPL-IRRD-QGKFW-A-WAQMWWFREPT-NIPTA-D-Y-CHSW--WQ--SR-ADLQ-NDRDMGP-EADASFYVEFWYWVRCAARTYGQQLGIIWNNRLKTRNPCPYSADGIQNKENYVFWWKNMCTKSHIAFYYCLQNVAHYTHDVTAEFNPVKCIDWLQGHMVLSSWFKYNTECKKLFVHEKCECYRM----FCGV---VEDIFGVRFH--TGWKHLSTAKPVPHVCVYNPSVQERRNINDFYIF-YEIAPAVKNLVLSAQPLHDYTKKCAFNYTPITITRIISTRNQIIW-AHVVIACQFYSPHQMLLIELAMDKYCADMNVRRSTEGHQPMHACRSTFGPGMAAKEPLVTFTLVAFWQWPNHEFQYVYMYTED-KIIQIG-PHLSN-GCEMVEYCVDC-YAK-RPCYRAYSAEAQYWRMITEAEDYSYKTRNAIAATATVRGQ-YCHPFRWLGIVWM-AHHDC-FFANECGTICI-PQMAEMRPPETTPYEI--DIIFMMF-WKE--HMSTTIL-DVVGMYRP-ATFSHWHDAHH-QCEPYLTPL-MCQSKLVFDAAFT--QVG-VKGVW-YHTEKLELMAGFNHM-K-FKKEEAQ---QSCFYWFQDCPDYDPPDAVRKTDEKHIRAHGEIWWLMRYYCMYHILHI-ASRHEWMHLRWDQACTNPGY--ELFE-F',s2)


        def test_ba5f_rosalind1(self):
            '''BA5F 	Find a Highest-Scoring Local Alignment of Two Strings'''
            score,s1,s2 = get_highest_scoring_alignment('KGWLILYGMTAEEIVLVMSVNFTDPVERCWDCEYHEHNLIWEGAKNNHNGFVKHNSWGGTGENQWCYANATMCVTSRCHLRKFLHEPWASNSRLEKKPMWHEPLAHQVAVLMWLMDGQCGWQSMVKKPMDAYTHCAASFWMVMYIIYQQCFRHYNICQVQCDPHVRVFCLWNSRHHIFNEDYWDRLRQSHQVQILPMMMQTHDCAKVQWIRRQFGNYNYEPSRRIKTTEGVCFKFPNWDWVATDGEHSEETPYFPHNNLMHDWADQFSAKRTQSHTDEQTGQFEEHCFFRPWCTWNWEKKGDASGWQLEVSSMYKGNGRGHANFLLDRDGDMSFCGKDGMHDYAQMFFTRYLGCRFMMIICQNMCDQRKDKHKHMQCPCGISHKCAKENHHCKMIAITIRPDAIEVAAYAGFPSCICWFWIPHLCWRSEIAWAVGESGDKIEQVMNRWWEIWASRAEPPMNSNASMQHISGYFQCGPTDSVCDAKPCILIYHSQTIPRTQYDTFNIYERGSNTKRQICHGGLTCSNSVLHVCREMWGMLPSVFNFPWWLWWFHQTIQEGMIARVWRLCFCYWCNPQKFYTHSRDSAAYMNHTAEWPYGGKRRQGGCSEFMDDEQMTDHGRLYIDHCMQQLEDGSISNYGFTFNEKYNCCISDDCYANCVAIEKMANVPWAHRTLTTLLIQHWECMNADYNTGSTHPFCEPADNPVHMDQCTNVENFDSRCSKVKCHHHESHLHDADKYAITELQAWHWWMLNVQEITDWWRMAAQANMELFYKCDPWGNWTKLALMVFLPSWNGRYRMCFCCFSARMCCTMVPEAICGPITVSAASQIHTRDLILCNTQQQLVVHAQSDQVKGLNKHMREFGAPHWWIRDRSPHHVPICSRWCKITNNHTVCDEWSYNGMPMEFLVNIRYLNLSDQYCHDRARHSKTSVIRCIDPGLWYAGSGYQKAGVDESFCVMNCDEGKSQYKTQWHMPMRNFIHPGEFWMEQEVFSYSCDLMIVYTGYPIAWAPVDFPTIHMQEGKIKQSKYLHMIIIAGLYYHCPWPFVISDYRFPLYCRINDCMRESHEFGEQYLGYKNVANKCWDFVMDSVSEKLYWTFVKAYKTWHRWCMMCMDGISFEDNHMDKWGKTIWVSKGCITQCMLKSKSTRTPRTTWADAIGDPHNYCEMNRERPSPVVGAPQIRCDSCMTKGTPDFIHAMNGPDTTQYYRRSRSLWDIYMAIECTYLQECHENTDWLCRQAPGRNGEWGPAVSVVYINHLGHKDNYNSLQHIVRHWWVPEIMEPIYKCPSVEWVQQWFRRHHDQIATDMNCAAPVGLAECQLGVEFQMMWNEWDPFDSNLMFCCGCQTEATPPAMLCEPVDLHNCIALCAMISMWNYSWYLRSYGKPNKFWILCSFTEIYVPYRAIWQVSDHQCWGRVYNYGFLDFCREWHRPKDTGGQFLSKSSGFMLVEFRKSACFMFGDFPNCLWTNGFMAWHKLVPCNMSHHNSIQIGPKAYRHKMFGFNYLQQFMFMWKMMRKDHANERNYKLGWCYWWITWGDMTRQAQVGYVIYELVNKVWFTMQRVTKMETDASRQCFDVERHNDPTVFRSQDAGDHVYLVRTRSFATSSYCGMAVPSHQDKFFQKDCDCMQRSKENSFQSGETDKYNGYDHRFRKGWDNEFPTSPGVIAHWWFNNNNGMTQLGRVVSWMFNFYVSRTPTNLTISAEFNEVCWQQYLGMIASYCHQLFVVMKVVTRVDFASQGFENEMPPQSVASYNVWQCNGMYWHDEKTQTSSRMFTMDGGNDQVSIIGVLRVEMEVQARWCWRSFFNIIMNSHGNHSAIDFPGTVPEMDNPRYHEKDALIYEDACLFPWYNMKNQFIIASHWQYWDMTWNRCLDFDVQHTNYDYQRFPQTNEYPAKEMCRKNTSHQAPRIETRQMYTHPYNKQCFHYTPSAFCAAWTSTVRYTCAKRETYVLPIVRDMLPGDMNYFHRYRFCIQDVFNWVCGWQIVPMMYVRYILENTSPPLWCFYIIEVKGWIMYNEEKFWRDANHSMTIAQYCCQCWEWQLEMNARKDQCAKLETIDWHQMESKYSMRRQHKKFFKWGYRMEGHVNERKSTAHFKKVANWLLEAQLVGNLHGNEQYMAQLKVLVMFFSNKFPCRVADQAAMRKVMCRAVHEFVLYKEQTVHRTVDIINEEGEGAYAYWVHEQGHWGLQPAFVWGFRQAYWIKYRQCDSMKHQWWDAPFHSPFNWCTYFWILWFTMFTVKYWDVCQCMYWDMHIRNIVVSQCLCWGGYRQVEYQSAEWAKAFKLQCVWKCLKWSRVALTCANTASIMCDEHSTRNHVANYSMQNHSAQSTPCIWAAFWAKQCMTFLPYVDFAEIIFTELAYHRGEDSCKGKMKAGFRHCWEVNRGNHGNEGCQGRMIPWLKHAVVLPYNRGWGYYCFKKARPVVAEADRGAVDIHACDWRCKSALTKIYRIKKQIGDYVQMVCAKFPWKVAVINVLEWPSQLDDYAFTIDSRNRNCESVHACEPHHRNQGWNEMDHNNHQGPWRENWWDLFKFDNHFGVKLWWQLPARNVMDKLVIFADIRDWFYKKVHYNIERQMGDCEIEMMHTKMDMFTSCQNIHYKHTMCMVYHIYIFMVHCGRENAPISNKMMPLISKENGFHRWGFLVNCYQFGYCGGIIMMNERGVMGRQMREMLMGFRLVGKPNCHNMCVLNMMCHNTKHMKME',
                                                        'SVNILVTWMQKWPHTTLMLRGLWSKSMDRTYTHDYMGRKRHIAQQLKSKAKHSAKPWVNMRISWMYQIYRGRKVYEMERRCQPGAMWVGNECGCMIEHNWEYFRLYAWVKYFSIRYTQCIFHNIKDPDMDGQHLLKRMDKVRVKHFDTYDNHLMATYMLLPCIHESIGCNTMSFDGEHKPAILGNAKEDTDKGASHHLYGKALDMDPFTWCVNDDQPMMFVKVSYWTYRRTQKYGCFRCSLGWGHDPCSTCPVRSNNITELDEFTNEKFGQCEWDLEGRCVKQECWIMCSWAVSQVQTMVFMWLDQLKIFFNAMMWYAEMVEANHMEFWPEEDETLTGMTHENEMECSQFGPSWVDGMDQSWGEGEGPFCHQNCGPLVMYFAWNHWSHRWVFIYRLTGNYITIMKCRKMDVVKMSCNTCGWYQRCFEPWFLAVCACVGFCTHYDVWTFKDIAYMRSRSYQIFMHVYCTGNGVEFAYMNKYAKDHWVEEEKMIQVIWPLMAMWFHFMYRIQFKYKAHEHLRFAECHHSATMAAKREHAQTNFRFKQVKGTDTENLTDWQKLARCPYIQQICVLYDEHLYPGNTYWYEWTCNRAARILKGTISYIMHKRWGEFLCETITAKDVKSEEWYYEGDVSYLECGSYTAWGCDTKTESSHRWRVIFWHYCFWPNTQFIIQNQGSPCTFDAYPSCYLMFFHACKWEPAGAARCQWPIPERFSHTERHRAHHENPWSHHPGKCNMGCRSLYTWMSIVPTLGVNNCWQSLPDAGVAYCMYGYAYAVGHIIQGQSMEHETINEHFRVCADSIVFNECCEAVNPAPDQCWQWCWGIMACLVIEMGDCHRCLNNKLNNHMANWAEEKCHMFNRQHPHVFNDTINTGQQGAPGWKKEFDVMTYAISHGLAGPIYRVYKTGWLTSAAARIHTRDLTLCNTQQQLVVHQSQLVQWLMGLNKHMRELASMHMEFGAPDKYKSESWPICSRWCKITNNHKPVCNSKRPACMHMEFLVNINLSDQYCHDRARHSSGYWKAIPQTYVDESFCVLECDEGKSQYKTQWHMGKCAHNIMRRFIHPGEFWMYWGQWVFSSSCDLIAWADVDFPTEHQMKIIAGLYYHCPHGFYLPGSEHGAYFEVISDYWFPLYNENDFINACMRESHEFGEQYATKCWDFVSDSVSCKLYWTFVKAYITWHRWECGMCKVMIHWNICMDGISFCDNWVSRYSNWLNGELCATCKTTQCMKGGTKSKTHLDFRWWADAIGDPHNWCEPNDERPSPVVGAPQIRCDSCMTTPDFIHAMNGPDTPDLQCYEVPQYCYSQRSRSLWDIYMAIECNYLQDNLCRQAPGMNGEWGPAVANCGVTAAHKDNYNSRQHIVRHFCFWWVPEIMKPIYKCPSIEWVQNWLRRLHIKDTACHDQPVENTLPAATDDACAAPVGLACGMCPHNNCQPGVHEQMMWREWDPFDSNLMFCIGCATPPAMLCEPVDLHNCIALRMAPEPNKFWIKCSFTEIYVPGKYRAIWQVSDMTWCMGSLIRQYNYFLMFCREWHSPKFIWVMYIEGTGGQFLERGDNKMLVNFSWDRCPLNCFMFGDFGFMAWHMDFSLVYNKCEHMYVFPMSHHNSMQIGPKAYRHKMFGFNYLQQFMDHLNEGRGCFSNYKLGWCYWWITWGAKCFIGHPMPRQQQVGYVIYELQNKSQFVMQRVTKMETDASRQCDPTVFRSVYLVEASYCGMGFCHHPGTVPSHQDDFFQKDCDNMQRSFQDGETDKYNGYDVYTQASEARMHESNKKPSMCCLRVEVHAAHWWFNNNNGMTQLGRVVSWMFMKWMNPMEQFYWAKRTLKHSRTPTNLEVCWWFCPFYSLPQAEYCHQWFVVRKHRFVTWYKVVTRVDFASQGFENRMPPQSVASYNVWQCNGMYWHDEKTWFAMTPKTLNDQYIGHARVELEVLMQSNSFSAHFYNKQASYGTLLPWCNNAGCVPYGSIWCHVLGIWTRYNPFDPGPHFIVGIYHHPNSPHHSCPQQWTVYTAPIFRIFMSDVFMKMVYSSRTSDTSTPVLRSRIHDDNSKHLGNCSMLTEIGKRAAWFWQAVLQQDNRSNFRCYRHLAWFRWKFQRVTIWWSRRFAHRQTKKGCFGTIHPSGQYNFMTTWERQGCEIFQVNRTTKHTVKDEIDPHHCVNNCHDHQTEPDNPGYHIHFEKRSVHRENQYVECGRGWQTQLTGVYEGDSTQCRHEGHTDLIGPYAFSSNNKNELKQCALRDIWKEGVDKYLRRQEDQEIRSYIHNPDQFTRISYHCNMSPAVCWDPNYEWSRQMTIDMHLGFRRWPLYNQKVDWSAALKVRSALTAETADSNMCVRRNIMVAPADMAPGAEPRQQNNFLYIVCTFIESGACSRNFNHQAQGDNNKHIQMMWMWVNSWEKGIIKPKHDNAAPMLGYSPSWVITPLAEHFPPRCNGACFCNDGWLLNAVNVCFDHTGQYWEKFPMHMIKHEHSSPEQEMGFTKYPFDRVDHSKAHKTGAKAWCSRKAQFALEYFYYERYGQDLDCWHQVLNPRIISVFITFTHLKQYWVLSSCVPMMTLTARPETSKIVWYTRALPIMPPWTLEMEPKDWEHAQVNNTHADKKLHQKSWSEHMEAYEGRDKNLKAQPYACQMPGHVYCMSGIFKDYECVPAGVGICAFFLAFQFQVFWLVIGLGRSVKIIWELLYIYHYKMVTIPNPFGYKWIEKFNIVFKEIVFHNTLCNENGTVLELVWKHAPIMNFSHDEQHCMWRRGMQQECDTGPDHCETIYGVYPRPTPMAIANYETHKQGTASMGVSKFQVNQHRVNNTVSFM',
                                                        weights = PAM250(),
                                                        local   = True)
            self.assertEqual(3272,score)
            # self.assertEqual('EANL-Y',s1)
            # self.assertEqual('ENALTY',s2)

        def test_ba5h_sample(self):
            ''' BA5H Find a Highest-Scoring Fitting Alignment of Two Strings '''
            score, s1,s2 = FindHighestScoringFittingAlignment('GTAGGCTTAAGGTTA','TAGATA')
            self.assertEqual(2,score)
            # self.assertEqual('TAGGCTTA',s1)
            # self.assertEqual('TAGA--TA',s2)

        def test_ba5k_sample(self):
            '''BA5K Find a middle edge in the alignment graph in linear space.'''
            self.assertEqual(((4, 3), (5, 4)),
                             FindMiddleEdge('PLEASANTLY',
                                            'MEASNLY'))

        def test_ba5k_rosalind(self):
            '''BA5K Find a middle edge in the alignment graph in linear space.'''
            self.assertEqual(((514,519), (515,520)),
                             FindMiddleEdge('AMEIATWSDTKMSPDTNMVYHEFDNCWALYNNTEFCKYDDKVFFRELNYHQISMWEMNHDIHEETCADDTDSNTQAIRVLCMKIGCWTEFMSSAQRYDFEPHNPESMQDFLREEYHLKWVAWCCSEECYFRQIRPLVALRCDIMYQCSRFIYYNYQNYMEFCKQQEDDIHFIKYNVYSSSLGHFQLRFPSCFFIMPCPVIVSDSVCTPDFHADSFHHPFVLYRGYYWEPFKVGTSRQQKMVGQKDDGCMPCIPCRLIIWLGMMKGWRPTLDMNAKMWHSAKTMAPWIFIADWAFHKYTKSSWMKTQEWHDCTGFAHGLCDALTEIGDHMVATSKSAPSHLITNGLFNRLPELPFTSHFRPDCCWFEKAYEMCISHFFVCTTYEKPSENGNYKNRTNRQITVSEDKHPDVMAEGVDRAYHREQDWCIGQGKQYKAISMDCVIMHWMPLGPTTNQLQPWYKFKNIAPVDFWDKHWCHKGQNSQWQPNRCSFPPPVIGHHINPACDFWESYSYYNFWNHDCYTDWSETMYTLVGETPDVWEFLYKFGIQWEKGGKEKTSSTNTFDFKDKYRMQHYDYDSPNSPTHFMREEKRYSQMDGLENHNYDDHHWKIRNRYYVHGLPWCWKADRPRYFTHSPICSMLEHVTYNNNPMGGSCQCLYETMCWAHFYRPEADFYDMQEQIHVNFWYLNMIDSYNQEHIILYGSGQNAWYCDGTWSIYKGAGAITLDNFMQSARSKQFGTTCQYMPTVRSQLIQASWWFPFASEHPSRYRDIPAECQCHPRKGDKMVWDPQHGAHIDWMNDSYYQCDHDYVALMDYSGFNNLHMFECKTCPPAIIEQIMLYIDYGDRYGDYETTPVRMGSKDMFCASSNDTKDHLRSAGYVEMVGWGTVLDVHQGCQEGWHSRGPNCIVPFDNTMPLVHWHFIHWGQQIRASMGQNFFVCIMLYDPDNPSYNFSHVDQHLVYHFGHVCEEPHGVIIIHEMLQHNSMHEIMYDMVCYYWASLDCDWVAVLGYECGGYGAEYQANKGYLEMDDNEVHDTV',
                                            'RCVKPFIAYGSDDLGHQDCHPEPTHDFWARFEETDHEMHAEFWTNCFDMFDPHGRNDLNCPTHSRGWLMSNPFCSDESETGFMKSHRWHGNQTMSWVKDEYRRSANCRPADGQFTRVASFIRSSQDFDYGARGPCKWTEFTFICDYKCTWVCVCTCMWSWWTCKSPEEHHEIQSKDTFFPQQKEQTHIHKSITRPALSHGIHWALCLQWFHWAEMGPCFNTTAVMNPSVYFMFTAAWQPQWNHFHKNKHVREDSHQADWETMFYRRQAVYTGCQWEIYISVWHFVEGFIMEMDRLMATSDDQNFVHDIEELRVGMNADHIQVDLRAHWMDVTVHPKGGSNFGPIPYWHSGLCGSANEESMKLADLYVKMELKDMAEVAAPGFTNDVLYLEPMCYMINYHIDHFLREQKKTERLGYMTDIVYLCMNQHHKQLIVWFVYFMGTQRNQVKNNKTNQHQRWYKFPNIADEVDFWNKHWVHKGQNSQWQTEMTRPPPVIVHHINPFWESYSYYNFWLFWPRRGRNHGCYTDWSETYMTLVGETPDVWEFLYKCKAEIAYNIQWSSTNTFDFKDKYRMQHYDYDSPNMPTHFMRCSMSRDNTNGFHCFEKPACTHHFWEFLKGVMYFGCWNLFSWEQLPGCMNEPTPVGWWFWPSNWYAHTNEHNWPIDSLPLGLNHCACHQEHYCTACKMQPGMMWHDEPKHYRVCQDRALQYAEICEPNNKYMMVPQCWQCHIVVREFNATAPNGYTGSPSRVSNFLIGWFMCAKKLMTRENACQHVHMCTQYDCENPQPFIYPRGSYTLPAHYEQRVTGIVFRVRFTVSREDCAWPYARKPMVYAGLINMYSDPEMVDPMTPDVDIHTWVENMLTQHGTLFGKHRNKMGRPKGYDFSQILSKMRNQPAACKWATGCEHCNKLMMCGRVCHLMNQGDIIAWACYVINARNLYVHDCYGYMHFGNYGCRDYRIDEMHPIFNGQKIKHWGRILDQMPFTANFVVRPLHDYREDQWKSPIRYQFPMDGMGDVDYNVAWTQDSRKIRHTCFRNIDFA'))

        # def test_ba5l_sample(self):
            # '''BA5L Align Two Strings Using Linear Space'''
            # score, s1,s2 = alignUsingLinearSpace('PLEASANTLY','MEANLY')
            # self.assertEqual(8,score)
            # self.assertEqual('PLEASANTLY',s1)
            # self.assertEqual('-MEA--N-LY',s2)

        def test_ba5m_sample(self):
            '''BA5M Find a Highest-Scoring Multiple Sequence Alignment'''
            s,u,v,w = FindHighestScoringMultipleSequenceAlignment('ATATCCG','TCCGA','ATGTACTG')
            self.assertEqual(3,s)
            self.assertEqual('ATATCC-G-',u)
            self.assertEqual('---TCC-GA',v)
            self.assertEqual('ATGTACTG-',w)

        def test_ba5n_sample(self):
            '''BA5N 	Find a Topological Ordering of a DAG'''
            self.assertEqual([5, 4, 1, 2, 3],
                             create_topological_order({
                                 1 : [2],
                                 2 : [3],
                                 4 : [2],
                                 5 : [3]
                             }))

        def test_ba5n_rosalind(self):
            self.assertEqual([5, 19, 11, 4, 22, 2, 6, 12, 0, 10, 1, 9, 7, 8, 20, 18, 14, 23, 21, 24, 3, 13, 25, 16, 15, 17],
                             create_topological_order({
                                 0 :[ 1,10,13,16,17,21],
                                 1 :[ 13,15,3,7,8,9],
                                 10 :[ 16,20,25],
                                 11 :[ 13,15,23],
                                 12 :[ 14,20,25],
                                 13 :[ 16,17,25],
                                 14 :[ 21,23],
                                 15 :[ 17],
                                 16 :[ 17],
                                 18 :[ 21],
                                 19 :[ 20,21],
                                 2 :[ 12,15,16,20,3,6,7,9],
                                 20 :[ 21],
                                 21 :[ 24],
                                 22 :[ 23],
                                 3 :[ 13],
                                 4 :[ 10,22,24,25],
                                 5 :[ 11,19,21,23,7],
                                 6 :[ 14,7],
                                 7 :[ 14,18,20,23,25,8],
                                 8 :[ 13],
                                 9 :[ 17,20,21]
                             }))

        def test_edit_sample(self):
            '''EDIT Edit Distance'''
            dist,_ = edit('PLEASANTLY','MEANLY')
            self.assertEqual(5,dist)

        def test_edit_rosalind(self):
            '''EDIT Edit Distance'''
            dist,_ = edit('IAEWFANIDRHKMCSFDSPNEVPNSWIWQWYRVEYFRMHWDYSSFACDYAPKTRTLQGDH'
                                      'RYMMANKEENEVTGVSVWNVEETYLCPFDLSRPHVMPKFTGEDNRLERMDRKYLPQCSAI'
                                      'GKDLCMVFEFNSDTSIREAIVDAAVFHRAGHGGNSNKEARHQTAHFNDTAWYELGNCEVI'
                                      'GRITNDGKTSKHRYMGVLCWECCLEYGPPSHHFEVATKIWPPTGQWHTIKSVTNPHPKSF'
                                      'RTDSCVNNRYNSMWTESGKVASTIAVLYKQYYIQVCKRDWNMMSRNHVTLTSLNKTMFAS'
                                      'FFWFRSCNTHSLHHHHMGECDWSEPMRPCWISMYPYCATHSQVNDTHIYVEHGARMKEST'
                                      'EDTQGKGVGGIAGDARRYNSQTLTITWFIKISKWKPVHGNPFRPPWVETHEYAYWHQKVI'
                                      'NGFVRRCGIKDHHKNCAYYITTEFVLNGLFKDANSPSYLDADYYDVQHGGGPFRVHTGTW'
                                      'MNAIITHFFLPSVPDCQAWLCAPKRNSLVTREVHCWNGRQYYPATLWHMVLIKGPYNAYN'
                                      'IGSMFCAWRCPEFNWSPACPCVYNNLWHWCCCCPQRLKMPFCWATPMPKMHSTFNDTTCE'
                                      'PCMRVGAWTEIITMNYCYGHDSFTGVFMIYMVRFVTQGRIPPRGCVPFHAIWFAQNRRIC'
                                      'DQMEDVSWQTMMHDYSKRCDPRQTTNNCVSNLEGQAQESKQFPHVQMPAEFTHHENGWMD'
                                      'CMYQKHGPIIGDPQTTTDVRWGREKCMEDIGSARQDWPFIHNDNRWKDPMSICSAHGIND'
                                      'ALAKYKGAAWISH',
                                      'IANQPPKSLYHWFINIPRQINGFTSFDFPNEFVYAHWEITDRHWYRVEYFRAHWDYSSFF'
                                      'CMYAPDYHYWWVTRELQGDHRYMMANKDEEVIRVPPSKVSVWNVEEGPFCLSRPHSHKGL'
                                      'MMPKFTGEDNLTTWWYVYAKYLPPHAIGKDLCMVTSSREAIVDVSVFDCAQWRAGHDGAL'
                                      'SNKEARWQTAHSCDLRNGLRAWYEWVVAHNCEVIGRITNDGKTWEHRYMGVLWKSECPQA'
                                      'GVVDVATKIWHTIKSVTPHPKSCVNNNSMWTGKVWCTIAVLYKQYYIQWNNMWADRVVTL'
                                      'HSLNKTMFASFFWFRDHSKHHHRMGEHDDWSEKIAEYMCWISMYPYCATHAQVNDTHIYV'
                                      'EHGARLVWGGKISTHEGDDTQGKGMEWHYILLCGGIMGDARRYNITDAGFIKIHGPMAKW'
                                      'KPVHHTHEYAYVHQKVINGWVRRCFPRISKDHHYYITTLSNTWDYYDVQHGGGTGTWMNA'
                                      'IQTHFFLPSVRHTWCVCVMYLVAIPCHLEFEVHCWNSATLWHMVQMGQRKPHDAKGPYNA'
                                      'YNIGSYFCPKCFWRCPEQNWSPACPPVYDNLWHWRCPFCWATPMDVMHSTTNDTTMWKTM'
                                      'MVVEPCTDESIPRVGSYSWMIQEFMIYMVRPPVGCVPGHAIWFAQGRRICDQMYIVGRDV'
                                      'SWQTMMHDYHKRCDPRCVSNLEVADFVYFEAEAQAQSQMYAENGWMDCMKHGPIIGDTQT'
                                      'TTDVRWGGKIIGDARQDWPFIHNDKRWKDPMSIHWNDALAKYFGAAWVCQEHSISL')
            self.assertEqual(376,dist)

        def test_edta_sample(self):
            '''EDTA Edit Distance Alignment'''
            dist,s,t= edta('PRETTY','PRTTEIN')
            self.assertEqual(4,dist)
            self.assertEqual('PRETTY--',s)
            self.assertEqual('PR-TTEIN',t)

        def test_edta_rosalind(self):
            '''EDTA Edit Distance Alignment'''
            dist,s,t= edta('MEIPNQCTYTMNNGERMLKMDNKYADNMGVRQEENKQHREDISRWDDWGQFTTCEGFTRW'
                           'FKLHHLNAHQISPCLLRMPTWLQCKMTFAQMDLAIDPPTCHTFPMFTMHSHCIFDPSFYA'
                           'ARRCHMNKVCQFRYDHQIHCSMMQSGTTWWYMYNGCETMSWDYLVEMAVGWWVYRNRHNG'
                           'MMSYVDHGDDDCWHRGHHSSLPWKLDAYSALYCVHGFEYCSAYFKVDAFNLTCERPDITN'
                           'TQCRCEKPTEGPAYFNYTKTMSKGPEEWEYNCPCVYYEEVEGLANNRESSIPCSKMSNRC'
                           'RMRCLGHKWFSILPEMHKYNAQSQVKTREKAGRSYNWGNQPLVKFDSLEPMWFHRQHCQS'
                           'RLETWSTDPPSQDWQFFHICNTQKDAINSTWAGGRIFWTENDSRRKVCCMENKNKNLFRR'
                           'GWISGIDITMHYHLTSSHGQFGICDPEVHIVPCHMLGLTTFMQERNNLGQCPMVAIFEHH'
                           'QLPPPCEGYVLQKVHPSWHMSMMKKNEFDRPWAPVFARDAEHMQPCEPMIAPHIKFWNYC'
                           'HKKSGMCLGGQGGIDFHGYMYQNYHLCLFCTYKAMKFLVQDTFEDLISGEKYPMEATMTE'
                           'TFAWVCMTVFTTMRTFIVFWGDPYRIHELYMFSCAAHGYLRTADRSCKWTNKPHLSYGCH'
                           'LSKTYCSCKMIDAPMVCRCLMWIDHPFVHYPHMVWCMNYCCYNTTCRMHEPTFWHLALVF'
                           'CPCLYWQQNLEINEIRQGSLGHVGQSEMLQMTIMKKYMQPCMMWDCGFHCCGTYWDAMGE'
                           'TSGTEGIVTTRTRVTPCYRWQHLKKEGHALANTVPKMDLYQMRHHGLMGPVLDRAPKPCH'
                           'AQTEVDLVNIVPYNWES',
                           'SEIPNQCNEFMYMNYTVRIIYHPNNGERMHEQSDNHHFTAPMDNYHGILLRQRMGVRQEE'
                           'NRDIDLQCIVNHSYLREDISRWDDWDRAQFTTCEGFTRWFRYPKHLSEQLHHLNADQISP'
                           'CLLCMGTWLQCKMTFARFTMIWQHWMDLAIPPYTCHTNHSHCIFDPSFYLARKVLQFCSM'
                           'MQSGTTWGQQLVEMAPGWWVYRNRHNGMMSLAYVLTDVDHYDQMQHDDDDCWHRGHHSSA'
                           'PWKLDAYSHGFEYFSAYRFVDAFNLTCERDDTNTQCREGPAYFNYTSQNNPLHDFEMSKG'
                           'PEIAMWQEVIPCSKMSNRKRWRCLGNKIKSGDVGGKTASVSYNWGNQPLVKFDSIEMTFH'
                           'KQHCQSRLETWQQKDAINSTWAGGFRFWTECWNKNAYEHTVPWFMLFRCLAGGHGWISGA'
                           'TMHYHLTSSHIQFGICDPEMHIKLCFMLSYESNNLGQCPMVAIPPPCEGHPSGNMSMMKK'
                           'NEFTRPWAPFEARYAEHMQPCIPGFCQMREYRKSGMCLGGQGGIDFHGLMYQNYHLCLFC'
                           'TGKFMKFLGQDTFEYPMEMTESFAKVSMEVFTTMRTFIVFWGDYYKLFYRITELYMFDEV'
                           'DTDNNWHGYLRTAQKWECHLSYGCHLSKTYCSPMVCIDPFVHYPNGACELHVWCMAYVHR'
                           'TCCCYNTAWRMHEFCPCLYWQQNLEINELGHQNEMLQMQPCMFYEFMGKLQHCCGPYGET'
                           'ATEFFGTCSTCAIVTTTDCYRWQHLKKEGFKPSALQAHNTVPKMDLYQMRHHGLMGPVLD'
                           'RAPKPCHAQKFVIRDAVWHLVNIVPYNWES')
            self.assertEqual(406,dist)
            self.assertEqual('MEIPNQCT---Y---TM------NNGERMLKM-DNKY--A--DN---------MGVRQEENK----Q----H---REDISRWDDWG--QFTTCEGFTRWF---K-L----HHLNAHQISPCLLRMPTWLQCKMTFA------Q--MDLAIDPP-TCHTFPMFTMHSHCIFDPSFYAARRCHMNKVCQFRYDHQIHCSMMQSGTTWWYMYNGCETMSWDYLVEMAVGWWVYRNRHNGMMS--YV--D--HGD----DD--CWHRGHHSSLPWKLDAYSALYCVHGFEYCSAY-FKVDAFNLTCERPDITNTQCRCEKPTEGPAYFNYTKTMSKGPEEWEYNCPCVYYEEVEGL--ANNRESSIPCSKMSNRCRMRCLGHKWFSILPEMHKYNAQSQVKTREKAGRSYNWGNQPLVKFDSLEPMWFHRQHCQSRLETWSTDPPSQDWQFFHICNTQKDAINSTWAGG-RIFWTENDSRRKVCCMENKNKNLFRR--G---WISGIDITMHYHLTSSHGQFGICDPEVHIVPCHMLGLTTFMQERNNLGQCPMVAIFEHHQLPPPCEGYVLQKVHPSWHMSMMKKNEFDRPWAPVF-ARDAEHMQPCEPMIAPHIKFWNYCHKKSGMCLGGQGGIDFHGYMYQNYHLCLFCTYKAMKFLVQDTFEDLISGEKYPMEATMTETFAWVCMTVFTTMRTFIVFWGDPY----RIHELYMFSCAA-----HGYLRTADRSCKWTNKPHLSYGCHLSKTYCSCKMIDAPMVCRCLMWIDHPFVHYP------HMVWCMNY----CC-YNTTCRMHEPTFWHLALVFCPCLYWQQNLEINEIRQGSLGHVGQSEMLQMTI-MKKYMQPCMMWDCGFHCCGTYWD-AMGETSGTEG---IVTTRTRVTPCYRWQHLKKEGH---AL-A-NTVPKMDLYQMRHHGLMGPVLDRAPKPCHAQTEV--D----LVNIVPYNWES',s)
            self.assertEqual('SEIPNQCNEFMYMNYTVRIIYHPNNGERMHEQSDNHHFTAPMDNYHGILLRQRMGVRQEENRDIDLQCIVNHSYLREDISRWDDWDRAQFTTCEGFTRWFRYPKHLSEQLHHLNADQISPCLLCMGTWLQCKMTFARFTMIWQHWMDLAI-PPYTCHTN-----HSHCIFDPSFYLAR-----KVLQF-------CSMMQSGTTW-----GQQ------LVEMAPGWWVYRNRHNGMMSLAYVLTDVDHYDQMQHDDDDCWHRGHHSSAPWKLDAYS-----HGFEYFSAYRF-VDAFNLTCERDD-TNTQCR-E----GPAYFNYT---SQN------N-PLHDFEMSKGPEIAMWQEV-IPCSKMSNRKRWRCLGNKIKSGDVGG-KT-A-S-V--------SYNWGNQPLVKFDSIE-MTFHKQHCQSRLETW------Q--Q--------KDAINSTWAGGFR-FWTECWNKNAYEHTVPWFM-LFRCLAGGHGWISGA--TMHYHLTSSHIQFGICDPEMHIKLCFMLSY-----ESNNLGQCPMVAI------PPPCEG------HPSGNMSMMKKNEFTRPWAP-FEARYAEHMQPCIPGFCQMRE---YR-K-SGMCLGGQGGIDFHGLMYQNYHLCLFCTGKFMKFLGQDTFE-------YPME--MTESFAKVSMEVFTTMRTFIVFWGDYYKLFYRITELYMFDEVDTDNNWHGYLRTAQ---KWEC--HLSYGCHLSKTYCS------PMVC-----ID-PFVHYPNGACELH-VWCMAYVHRTCCCYNTAWRMHE--F-------CPCLYWQQNLEINE-----LGH--QNEMLQMQPCMF-YEF--MGKLQ--HCCGPYGETAT-EFFGTCSTCAIVTT-TD---CYRWQHLKKEGFKPSALQAHNTVPKMDLYQMRHHGLMGPVLDRAPKPCHAQKFVIRDAVWHLVNIVPYNWES',t)
    main()
