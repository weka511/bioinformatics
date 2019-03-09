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

import numpy as np

from reference_tables import createSimpleDNASubst


# create_distance_matrix
def create_distance_matrix(nrows,ncolumns):
    distances = []
    for i in range(nrows):
        distances.append([0]*ncolumns)
 
    return distances

def get_indel_cost(sigma,delta,i,j,di,dj,moves):
    return di,dj,sigma

def build_matrix(s,t,matrix,replace_score=createSimpleDNASubst(),indel_cost=1,get_indel_cost=get_indel_cost):
    def init_indel_costs():
        if isinstance(indel_cost,tuple):
            return indel_cost
        return indel_cost,None
    
    sigma,delta = init_indel_costs()
    moves = {}
    def score(pair):
        def reverse(pair):
            a,b=pair
            return (b,a)
        return replace_score[pair] if pair in replace_score else replace_score[reverse(pair)] 
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
                                               matrix[i-1][j-1] + score((s[i-1],t[j-1]))]
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
            index  = np.argmin(scores)
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

if __name__=='__main__':
    from Bio.SubsMat.MatrixInfo import blosum62
    score,s1,s2=align('PLEASANTLY','MEANLY',replace_score=blosum62,indel_cost=5)
    print (score,s1,s2)