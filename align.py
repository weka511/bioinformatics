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

# Common code for alignment problems

from numpy import argmax,argmin
from sys import float_info
from Bio.SubsMat.MatrixInfo import blosum62
from reference_tables import createSimpleDNASubst,createBLOSUM62,createPAM250
from helpers import zeroes,sign

# BA5A 	Find the Minimum Number of Coins Needed to Make Change 	
#
# Input: An integer money and an array Coins of positive integers.
#
# Return: The minimum number of coins with denominations Coins that changes money.
#
# http://rosalind.info/problems/ba5a/

def number_of_coins(money,coins):
    number = [0]                           # We will use Dynamic Programming, and solve 
                                           # the problem for each amount up to and including money
    for m in range(1,money+1):             # solve for m
        nn = float_info.max            # Number of coins: assume that we haven't solved
        for coin in coins:                 # Find a coin such that we can make change using it
                                           # plus a previoudly comuted value
            if m>=coin:
                if number[m-coin]+1<nn:
                    nn = number[m-coin]+1
        number.append(nn)
    return number[money]

# BA5B 	Find the Length of a Longest Path in a Manhattan-like Grid 
#
# Input: Integers n and m, followed by an n*(m+1) matrix Down and an
#        (n+1)*m matrix Right. The two matrices are separated by the "-" symbol.
#
# Return: The length of a longest path from source (0, 0) to sink (n, m)
#        in the n*m rectangular grid whose edges are defined by the matrices
#        Down and Right.
#
# http://rosalind.info/problems/ba5a/

def longest_manhattan_path(n,m,down,right):
    s=[]
    for i in range(n+1):
        s.append(zeroes(m+1))

    for i in range(1,n+1):
        s[i][0]=s[i-1][0]+down[i-1][0]
        
    for j in range(1,m+1):
        s[0][j]=s[0][j-1]+right[0][j-1]
        
    for i in range(1,n+1):    
        for j in range(1,m+1):
            s[i][j]=max(s[i-1][j]+down[i-1][j],s[i][j-1]+right[i][j-1])

    return s[n][m]

# BA5C 	Find a Longest Common Subsequence of Two Strings
#
# Input: Two strings.
#
# Return: A longest common subsequence of these strings.
#
# http://rosalind.info/problems/ba5a/

def longest_common_subsequence(string1,string2):

    # Calculate longest path through "map" defined by the two strings
    #
    # At each point we have a state, s, which is defined as follows
    # 1. Count of matches (==path lebth from (0,0) to here
    # 2. horizonal position of predecessor of this point
    # 3. vertical position of predecessor of this point
    # 4. If this point corresponds to a matching character (in both strings)
    #    this should be that chracter, otherwise an empty string
    
    def longest_path():
        
        # Used to update distamce at each node     
        def new_s(i,j,s):
            count_horizonal,_,_,_=s[i-1][j]
            count_vertical,_,_,_=s[i][j-1]
            count_diagonal,_,_,_=s[i-1][j-1]
            if string1[i-1]==string2[j-1]:
                count_diagonal+=1
            count=max(count_horizonal,count_vertical,count_diagonal)
            if count_diagonal==count:
                return (count,\
                        i-1,\
                        j-1,
                        string1[i-1] if string1[i-1]==string2[j-1] else '')
            elif count_vertical==count:
                return (count,i,j-1,'')
            else: #count_horizonal==count
                return (count,i-1,j,'')
            
        m1=len(string1)+1
        m2=len(string2)+1
        s=[]
        for i in range(m1):
            ss=[]
            for j in range(m2):
                ss.append((0,-1,-1,''))
            s.append(ss)
        for i in range(1,m1):    
            for j in range(1,m2):
                s[i][j]=new_s(i,j,s)
        return s
    
    # Build string from status
    
    def construct_string(s):
        i=len(string1)
        j=len(string2)
        result=[]
        while i>-1 and j>-1:
            _,i,j,chars=s[i][j]
            result.append(chars)
        return ''.join(result[::-1])
    
    return construct_string(longest_path())

# BA5D 	Find the Longest Path in a DAG  	
#
# Input: An integer representing the source node of a graph, followed by an integer
#        representing the sink node of the graph, followed by an edge-weighted graph. 
#        The graph is represented by a modified adjacency list in which the notation "0->1:7"
#        indicates that an edge connects node 0 to node 1 with weight 7.
#
# Return: The length of a longest path in the graph, followed by a longest path. 
#         (If multiple longest paths exist, you may return any one.)
#
# http://rosalind.info/problems/ba5d/

def longest_path(source,sink,graph):
    def initialize_s():
        s={}
        for a,b,_ in graph:
            s[a]=-float_info.max
            s[b]=-float_info.max
        s[source]=0
        return s
    
    def create_adjacency_list():
        adjacency_list={}
        for a,b,w in graph:
            if not a in adjacency_list:
                adjacency_list[a]=[]
            adjacency_list[a].append(b)
        return adjacency_list
   
    def create_weights():
        weights={}
        for a,b,w in graph:
            weights[(a,b)]=w
        return weights
        
    def calculate_distances(ordering):
        s=initialize_s()
        weights=create_weights()        
        predecessors={}
        for b in ordering:
            for a in ordering:
                if a==b:
                    break
                
                new_s=max(s[b],s[a]+(weights[(a,b)] if (a,b) in weights else 0))
                if new_s>s[b]:
                    s[b]=new_s
                    predecessors[b]=a
        return (s,predecessors)
   
    def create_path(predecessors):
        path=[sink]
        node=sink
        while node in predecessors:
            node=predecessors[node]
            path.append(node)
        return path
    
    s,predecessors=calculate_distances(topological_order(create_adjacency_list()))
    
    return (s[sink],create_path(predecessors)[::-1])

# BA5F 	Find a Highest-Scoring Local Alignment of Two Strings  
#
# common code

# BA5E 	Find a Highest-Scoring Alignment of Two Strings
# create_distance_matrix
def create_distance_matrix(nrows,ncolumns,initial_value=-float_info.max):
    s=[]
    for i in range(nrows):
        row=[]
        for j in range(ncolumns):
            row.append(initial_value)        
        s.append(row)
    s[0][0]=0
    return s

def calculate_scores_for_alignment(s,string1, string2, weights,sigma,init_predecessors=None):
    predecessors={}
    for i in range(len(string1)+1):
        for j in range(len(string2)+1):
            predecessors[(i,j)]=init_predecessors
            if i>0:
                s_new=s[i-1][j]-sigma
                if s_new>s[i][j]:
                    s[i][j]=s_new
                    predecessors[(i,j)]=(-1,0,i-1,j)
            if j>0:
                s_new=s[i][j-1]-sigma
                if s_new>s[i][j]:
                    s[i][j]=s_new
                    predecessors[(i,j)]=(0,-1,i,j-1)            
                if i>0:
                    s_new=s[i-1][j-1]+weights[(string1[i-1],string2[j-1])]
                    if s_new>s[i][j]:
                        s[i][j]=s_new
                        predecessors[(i,j)]=(-1,-1,i-1,j-1)
    return (s,predecessors)

def create_alignment(string1, string2,s_predecessors,i_start=-1,j_start=-1):
    s,predecessors = s_predecessors
    result1        = []
    result2        = []
    i              = len(string1) if i_start==-1 else i_start
    j              = len(string2) if j_start==-1 else j_start
    while i>0 or j>0:
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

# BA5E 	Find a Highest-Scoring Alignment of Two Strings
# Find the highest-scoring alignment between two strings using a scoring matrix.
#
# Input: Two amino acid strings.
#
# Return: The maximum alignment score of these strings followed by an
#         alignment achieving this maximum score. Use the BLOSUM62 scoring matrix
#         and indel penalty σ = 5. (If multiple alignments achieving the maximum 
#         score exist, you may return any one.)

def highest_scoring_global_alignment(string1,string2,weights=createBLOSUM62(),sigma=5):
    return create_alignment(string1,
                            string2,
                            calculate_scores_for_alignment(
                                create_distance_matrix(len(string1)+1,len(string2)+1),\
                                string1,\
                                string2,\
                                weights,\
                                sigma))


# BA5F 	Find a Highest-Scoring Local Alignment of Two Strings 
#
# Input: Two amino acid strings.
#
# Return: The maximum score of a local alignment of the strings, followed by
# a local alignment of these strings achieving the maximum score. Use the
# PAM250 scoring matrix and indel penalty  5. (If multiple local alignments
# achieving the maximum score exist, you may return any one.)

def highest_scoring_local_alignment(string1,string2,weights=createPAM250(),sigma=5):
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
   
    s,predecessors=calculate_scores_for_alignment(
        create_distance_matrix(len(string1)+1,len(string2)+1),
        string1,
        string2,
        weights,
        sigma,
        (0,0,0,0))
   
    i_best,j_best,s,predecessor=find_best_substring(s)
    if predecessor!=(0,0,0,0):
        predecessors[(len(string1),len(string2))]=predecessor
    
    return create_alignment(string1,string2,(s,predecessors),i_best,j_best)


# create_distance_matrix
def create_distance_matrix(nrows,ncolumns):
    distances = []
    for i in range(nrows):
        distances.append([0]*ncolumns)
 
    return distances

def get_indel_cost(sigma,delta,i,j,di,dj,moves):
    return di,dj,sigma

def score(pair,replace_score=blosum62):
    def reverse(pair):
        a,b=pair
        return (b,a)
    return replace_score[pair] if pair in replace_score else replace_score[reverse(pair)] 

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
                                               matrix[i-1][j-1]               + score((s[i-1],t[j-1]),
                                                                                      replace_score=replace_score)]
                froms                       = [(i+di_down,j+dj_down,di_down,dj_down),
                                               (i+di_left,j+dj_left,di_left,dj_left),
                                               (i-1,j-1,-1,-1)]
                index                       = argmax(scores)
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
    distances = create_distance_matrix(len(s)+1,len(t)+1)
    distances,moves = build_matrix(s,t,distances,replace_score=replace_score,indel_cost=indel_cost,get_indel_cost=get_indel_cost)
    if showScores:
        for row in distances:
            print (row)
        for k,v in moves.items():
            print (k,v)
    return backtrack(s,t,distances,moves,showPath=showPath)

# EDIT 	Edit Distance http://rosalind.info/problems/edit/

def edit(s,t,indel_cost=1,replace_cost=lambda a,b: 1,show_matrix=False):
    
    def dynamic_programming(s,t):
        matrix=[[0 for j in range(len(t)+1)] for i in range(len(s)+1)]
    
        for j in range(len(t)+1):
            matrix[0][j]=j
        for i in range(len(s)+1):
            matrix[i][0]=i

        for i in range(1,len(s)+1):
            for j in range(1,len(t)+1):
                matrix[i][j] = min(
                    matrix[i-1][j]   + indel_cost,
                    matrix[i][j-1]   + indel_cost,
                    matrix[i-1][j-1] + (0 if s[i-1]==t[j-1] else replace_cost(s[i-1],t[j-1])))
                    
        if show_matrix:
            for i in range(0,len(s)+1):
                ii = len(matrix)-i-1
                print (s[ii] if i>0 else '#',matrix[ii])
            print (' ',['#']+t)
        
        return matrix[len(s)][len(t)],matrix
            
    return dynamic_programming([s0 for s0 in s], [t0 for t0 in t])

# EDTA Edit Distance Alignment http://rosalind.info/problems/edta/

def edta(s,t,indel_cost=1,replace_cost=lambda a,b: 1):
    def extract(s,t,matrix):
        m  = len(matrix)-1
        n  = len(matrix[0])-1
        s1 = []
        t1 = []
        while m>0 and n>0:
            moves  = [(m-1,n),(m,n-1),(m-1,n-1)]
            scores = [matrix[m-1][n]+indel_cost,
                      matrix[m][n-1]+indel_cost,
                      matrix[m-1][n-1] + (0 if s[m-1]==t[n-1] else replace_cost(s[m-1],t[n-1]))]
            ss     = [s[m-1],'-',s[m-1]]
            ts     = ['-',t[n-1],t[n-1]]
            index  = argmin(scores)
            m,n    = moves[index]
            s1.append(ss[index])
            t1.append(ts[index])
        s1.reverse()
        t1.reverse()
        return ''.join(s1),''.join(t1)
    d,matrix = edit(s,t,indel_cost,replace_cost)
    s1,t1    = extract([s0 for s0 in s], [t0 for t0 in t],matrix)
    return (d,s1,t1)

def overlap_assignment(v,w,match_bonus=+1,mismatch_cost=2,indel_cost=2):
    def dynamic_programming(v,w):
        distances = create_distance_matrix(len(v)+1,len(w)+1)
        path      = {}
        for i in range(1,len(v)+1):
            for j in range(1,len(w)+1):
                moves           = [(i-1,j),(i,j-1),(i-1,j-1)]
                scores          = [distances[i-1][j]   - indel_cost,
                                   distances[i][j-1]   - indel_cost,
                                   distances[i-1][j-1] + (match_bonus if v[i-1]==w[j-1] else -mismatch_cost)]
                index           = argmax(scores)
                distances[i][j] = scores[index]
                path[(i,j)]      = moves[index]
        
        i        = len(v)
        j        = argmax(distances[i])
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
    

    
# BA5N 	Find a Topological Ordering of a DAG 
#
# Input: The adjacency list of a graph (with nodes represented by integers).
#
# Return: A topological ordering of this graph.

def topological_order(graph):
    def number_incoming(node):
        n=0
        for out in graph.values():
            if node in out:
                n+=1
        return n
    
    ordering=[]
    candidates=[node for node in graph.keys() if number_incoming(node)==0]
    while len(candidates)>0:
        a=candidates.pop()
        ordering.append(a)
        if a in graph:
            bs=[b for b in graph[a]]
            del graph[a]
            for b in bs:
                if number_incoming(b)==0:
                    candidates.append(b)
                    
    if len(graph)>0:
        raise RosalindException('Input graph is not a DAG')
    
    return ordering


def unwind_moves(moves,score,i,j):
    ss = []
    ts = []

    while i>0 and j > 0:
        i,j,s0,t0=moves[(i,j)]
        ss.append(s0)
        ts.append(t0)
    return score,ss[::-1],ts[::-1]
    
def san_kai(s,t, replace_score=blosum62,sigma=11,epsilon=1,backtrack=unwind_moves):
    
    def match(pair,replace_score=replace_score):
        def reverse(pair):
            a,b=pair
            return (b,a)
        return replace_score[pair] if pair in replace_score else replace_score[reverse(pair)]

    lower        = create_distance_matrix(len(s)+1,len(t)+1)
    middle       = create_distance_matrix(len(s)+1,len(t)+1)
    upper        = create_distance_matrix(len(s)+1,len(t)+1)

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
            index        = argmax(choices)
            middle[i][j] = choices[index]
            moves[(i,j)] = [(i-1, j,   s[i-1], '-'),     # Comes from lower
                            (i-1, j-1, s[i-1], t[j-1]),  # Comes from middle
                            (i,   j-1, '-',    t[j-1]    # Comes from upper
                             )][index]

    return backtrack(moves,middle[len(s)][len(t)],len(s),len(t))

def FindMiddleEdge(s,t,replace_score=blosum62,indel_cost=5):
 
    def update(j,col1,col2,s,t):
        col2[0] = - (j * indel_cost)
        for i in range(1,len(col2)):
            scores = [col1[i-1] + score((s[i-1],t[j-1]),replace_score=replace_score),
                      col2[i-1] - indel_cost,
                      col1[i]   - indel_cost
            ]
 
            best        = argmax(scores)
            col2[i]     = scores[best]
        return col2,col1
    
    def explore(s,t,middle):
        column_A = [-(i * indel_cost) for i in range(len(s)+1)]
        column_B = [0 for i in range(len(s)+1)]
        for j in range(1,middle+1):
            scores,previous = update(j,column_A,column_B,s,t) if j%2 ==1 else update(j,column_B,column_A,s,t)
        return scores,previous
    
    middle        = len(t)//2  
    from_source,_ = explore(s,t,middle)
    to_sink,pre   = explore(s[::-1],t[::-1],middle)
    length        = [a+b for (a,b) in zip(from_source,to_sink[::-1])]
    aux           = [0 for i in range(len(s)+1)]
    post,_        = update(middle+1,from_source,aux,s,t)
    next_steps    = [a+b for (a,b) in zip(post,pre[::-1])]
    
    return ((argmax(length), middle), ( argmax(next_steps), middle+1))
        
if __name__=='__main__':
    from Bio.SubsMat.MatrixInfo import blosum62
    import unittest

    class Test_5_Alignment(unittest.TestCase):
        
        def test_simple(self): # Simple test
            score,s1,s2=align('PLEASANTLY','MEANLY',replace_score=blosum62,indel_cost=5)
            self.assertEqual(8,score)
            self.assertEqual('LEASANTLY',''.join(s1))
            self.assertEqual('MEA--N-LY',''.join(s2))
            
        def test_ba5a(self): # BA5A 	Find the Minimum Number of Coins Needed to Make Change 	
            self.assertEqual(2,number_of_coins(40,[1,5,10,20,25,50]))
            self.assertEqual(338,number_of_coins(8074,[24,13,12,7,5,3,1]))
            
        def test_ba5b(self): # BA5B 	Find the Length of a Longest Path in a Manhattan-like Grid  	
            self.assertEqual(34,\
                             longest_manhattan_path(4,\
                                                    4,\
                                                    [[1, 0, 2, 4, 3],\
                                                     [4, 6, 5, 2, 1],\
                                                     [4, 4, 5, 2, 1],\
                                                     [5, 6, 8, 5, 3]],\
                                                    [[3, 2, 4, 0],\
                                                     [3, 2, 4, 2],\
                                                     [0, 7, 3, 3],\
                                                     [3, 3, 0, 2],\
                                                     [1, 3, 2, 2]]))
            self.assertEqual(84,\
                             longest_manhattan_path(17,\
                                                    9,\
                                                    [[2,3,4,0,3,1,1,1,1,1],
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
                                                     [3,1,0,3,2,3,2,2,1,4]],\
                                                    [[1,0,4,4,3,3,1,0,4],\
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
                                                     [2,3,1,4,4,3,4,0,3]]))            
                             
	
        def test_ba5c(self): # BA5C 	Find a Longest Common Subsequence of Two Strings  
            self.assertEqual('ACCTTG',
                             longest_common_subsequence('AACCTTGG',
                                                    'ACACTGTGA'))            
	
        def test_ba5d(self): # BA5D 	Find the Longest Path in a DAG  
            n,path=longest_path(0,4,
                                [
                                    (0,1,7),
                                    (0,2,4),
                                    (2,3,2),
                                    (1,4,1),
                                    (3,4,3)
                                ])
            self.assertEqual(9,n)
            self.assertEqual([0,2,3,4],path)
            n,path=longest_path(11,18,
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
        
  	
        def test_ba5e(self): # BA5E 	Find a Highest-Scoring Alignment of Two Strings
            score,s1,s2=highest_scoring_global_alignment('PLEASANTLY','MEANLY')
            self.assertEqual(8,score)
            self.assertEqual('PLEASANTLY',s1)
            self.assertEqual('-MEA--N-LY',s2)
            

        def test_ba5f(self): # BA5F 	Find a Highest-Scoring Local Alignment of Two Strings 
            score,s1,s2=highest_scoring_local_alignment('MEANLY','PENALTY')
            self.assertEqual(15,score)
            self.assertEqual('EANL-Y',s1)
            self.assertEqual('ENALTY',s2)
            score,s1,s2=highest_scoring_local_alignment(\
                'AMTAFRYRQGNPRYVKHFAYEIRLSHIWLLTQMPWEFVMGIKMPEDVFQHWRVYSVCTAEPMRSDETYEQKPKPMAKWSGMTIMYQAGIIRQPPRGDRGVSDRNYSQCGKQNQAQLDNNPTWTKYEIEWRVQILPPGAGVFEGDNGQNQCLCPNWAWEQPCQWGALHSNEQYPNRIHLWAPMSKLHIKIEKSSYNRNAQFPNRCMYECEFPSYREQVDSCHYENVQIAFTIFSGAEQKRKFCSCHFWSNFIDQAVFSTGLIPWCYRRDDHSAFFMPNWNKQYKHPQLQFRVAGEGTQCRPFYTREMFTKVSAWRIAGRFAGPYERHHDAHLELWYQHHKVRTGQQLGIIWNNRDKTRNPCPFSAYYNKLPWWKINQNAFYNCLQNIAHSTHDETHEFNPVKCIDWLQGTMVPTECKKGFVHEKCECYRNPGPPLHDMYHQMEDIFGVRFDCLTGWKHLSDYNPCQERRNINDFYIFAYEIAPAVKNLVLSPQPLADATKKCAFNYTPLDQSPVVIACKWYIHQPICMLLIVLICAMDKYNAHMIVIRTTEGQQPMHACRMTEGPGMCMKEPLVTFTLPAQWQWPNHEFKYVYMYVLNYHLSQYTYTDEGHAGGQHYSFNVAVDVGMAWGHNRCYCQPACYSQQETQTRTIDYEKWQYMKHQAFKWGLWFCEQERHAWFKGQNRCEMFTAKMTRMGADSNLDQYKLMLAQNYEEQWEQPIMECGMSEIIEIDPPYRSELIFTFWPFCTYSPWQNLIKCRCNNVIEEMDQCVPLTFIGFGVKQAGGIQAWAFYKEEWTSTYYLMCQCMKSDKAQYPYEIILFWMQPMDTGEQEPPQQNMWIFLPHSWFFDWCCNAPWSEICSSRHDHGQCQDAFYPCELFTVFDDIFTAEPVVCSCFYDDPM',\
                'WQEKAVDGTVPSRHQYREKEDRQGNEIGKEFRRGPQVCEYSCNSHSCGWMPIFCIVCMSYVAFYCGLEYPMSRKTAKSQFIEWCDWFCFNHWTNWAPLSIVRTSVAFAVWGHCWYPCGGVCKTNRCKDDFCGRWRKALFAEGPRDWKCCKNDLQNWNPQYSQGTRNTKRMVATTNQTMIEWKQSHIFETWLFCHVIIEYNWSAFWMWMNRNEAFNSIIKSGYPKLLLTQYPLSQGSTPIVKPLIRRDQGKFWAWAQMWWFREPTNIPTADYCHSWWQSRADLQNDRDMGPEADASFYVEFWYWVRCAARTYGQQLGIIWNNRLKTRNPCPYSADGIQNKENYVFWWKNMCTKSHIAFYYCLQNVAHYTHDVTAEFNPVKCIDWLQGHMVLSSWFKYNTECKKLFVHEKCECYRMFCGVVEDIFGVRFHTGWKHLSTAKPVPHVCVYNPSVQERRNINDFYIFYEIAPAVKNLVLSAQPLHDYTKKCAFNYTPITITRIISTRNQIIWAHVVIACQFYSPHQMLLIELAMDKYCADMNVRRSTEGHQPMHACRSTFGPGMAAKEPLVTFTLVAFWQWPNHEFQYVYMYTEDKIIQIGPHLSNGCEMVEYCVDCYAKRPCYRAYSAEAQYWRMITEAEDYSYKTRNAIAATATVRGQYCHPFRWLGIVWMAHHDCFFANECGTICIPQMAEMRPPETTPYEIDIIFMMFWKEHMSTTILDVVGMYRPATFSHWHDAHHQCEPYLTPLMCQSKLVFDAAFTQVGVKGVWYHTEKLELMAGFNHMKFKKEEAQQSCFYWFQDCPDYDPPDAVRKTDEKHIRAHGEIWWLMRYYCMYHILHIASRHEWMHLRWDQACTNPGYELFEFIPWVLRRYVVYDKIRYNYSYRNSASMEFV')
            self.assertEqual(1062,score)
            self.maxDiff=None
            self.assertEqual('YQAGIIRQPPRGD-RGVSDRNYSQCGKQ-NQ-AQLDNNPTWTKYEIEWRVQI-LPPGAGVFEGDNGQNQCLCPNW--A-W-EQPCQW----GALHS-NEQYPNRIHLWAPMSKLHIKIEKSSYN-RNAQ-FPNRCMYECE-FPSY-REQVDSCHYENVQIAF-TIFSGAEQKRKFCSCHFWSNFIDQAVFSTGLI-PWCYRRDDHSAFFMPNWNKQ--YKHPQLQFRVAGEGTQCRPFYTREMFTKVSAWRIAGRFAGPYERHHDAHLELWY-QHHKVRT-GQQLGIIWNNRDKTRNPCPFSAY-Y-NK--LP-WWK-I-NQ-N-AFYNCLQNIAHSTHDETHEFNPVKCIDWLQGTMV-P------TECKKGFVHEKCECYRNPGPPLHDMYHQMEDIFGVRFDCLTGWKHLS------D---YNPC-QERRNINDFYIFAYEIAPAVKNLVLSPQPLADATKKCAFNYTPLDQSPVVIACK---WYIHQPI-CMLL----IVLIC-AMDKYNAHMIVIRTTEGQQPMHACRMTEGPGMCMKEPLVTFTLPAQWQWPNHEFKYVYMYVLNYHLSQYTYTDEGHAGGQHYSFNVAVDVGMAWGHNRCYCQPACYSQQETQTRTIDYEKWQYMKHQAFKWGLWFCEQER-HA--WFKGQNRCEMFTAKMTRMGADSNLDQYKLMLAQNYEEQWEQPIMECGMSEIIEIDPPYRSELIFTFWPFCTYSPWQNLIKCRCNNVIEEMDQCVP-LTF-IGFGVKQAGGIQA-WAFYKE--EWTSTYYLMCQCMKSDKAQYPYEIILFWMQ--P-MDTGE--QEPPQQNMWIFLPHSWFFDWCCNAPWSEICSSRHD--H---GQ-CQDAFYPCELFTVF',s1)
            self.assertEqual('Y-P-MSRKTAKSQFIEWCDW-F--CFNHWTNWAPLSIVRTSVAFAV-W-GHCWYPCG-GVCKTNRCKDD-FCGRWRKALFAEGPRDWKCCKNDLQNWNPQYSQGTR--NTK-RMVATTNQTMIEWKQSHIFETW-LF-CHVIIEYNWSAF-W-MWMNRNEAFNSIIKSGYPKLLL-T-QY-P-L-SQG--STPIVKPL-IRRD-QGKFW-A-WAQMWWFREPT-NIPTA-D-Y-CHSW--WQ--SR-ADLQ-NDRDMGP-EADASFYVEFWYWVRCAARTYGQQLGIIWNNRLKTRNPCPYSADGIQNKENYVFWWKNMCTKSHIAFYYCLQNVAHYTHDVTAEFNPVKCIDWLQGHMVLSSWFKYNTECKKLFVHEKCECYRM----FCGV---VEDIFGVRFH--TGWKHLSTAKPVPHVCVYNPSVQERRNINDFYIF-YEIAPAVKNLVLSAQPLHDYTKKCAFNYTPITITRIISTRNQIIW-AHVVIACQFYSPHQMLLIELAMDKYCADMNVRRSTEGHQPMHACRSTFGPGMAAKEPLVTFTLVAFWQWPNHEFQYVYMYTED-KIIQIG-PHLSN-GCEMVEYCVDC-YAK-RPCYRAYSAEAQYWRMITEAEDYSYKTRNAIAATATVRGQ-YCHPFRWLGIVWM-AHHDC-FFANECGTICI-PQMAEMRPPETTPYEI--DIIFMMF-WKE--HMSTTIL-DVVGMYRP-ATFSHWHDAHH-QCEPYLTPL-MCQSKLVFDAAFT--QVG-VKGVW-YHTEKLELMAGFNHM-K-FKKEEAQ---QSCFYWFQDCPDYDPPDAVRKTDEKHIRAHGEIWWLMRYYCMYHILHI-ASRHEWMHLRWDQACTNPGY--ELFE-F',s2)
 

        def test_ba5n(self):  # BA5N 	Find a Topological Ordering of a DAG 
            self.assertEqual([5, 4, 1, 2, 3],
                             topological_order({
                                 1 : [2],
                                 2 : [3],
                                 4 : [2],
                                 5 : [3]
                             }))
                    
    unittest.main()   