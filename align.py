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

import numpy as np, sys

from reference_tables import createSimpleDNASubst


# create_distance_matrix
def create_distance_matrix(nrows,ncolumns,initial_value=-sys.float_info.max):
    distances = []
    for i in range(nrows):
        row  = []
        mrow = []
        for j in range(ncolumns):
            row.append(initial_value)  
            mrow.append([])
        distances.append(row)
 
    distances[0][0] = 0
    return distances

def build_matrix(s,t,matrix,replace_score=createSimpleDNASubst(),indel_cost=1):
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
                matrix[i][j] = matrix[i][j-1] - indel_cost
                moves[(i,j)]  = (0,j-1,0,-1)
            elif j==0:
                matrix[i][j] = matrix[i-1][j] - indel_cost
                moves[(i,j)]  = (i-1,0,-1,0)
            else:
                scores       = [matrix[i-1][j]   - indel_cost,
                                matrix[i][j-1]   - indel_cost,
                                matrix[i-1][j-1] + score((s[i-1],t[j-1]))]
                froms        = [(i-1,j,-1,0),
                                (i,j-1,0,-1),
                                (i-1,j-1,-1,-1)]
                index        = np.argmax(scores)
                matrix[i][j] = scores[index]
                moves[(i,j)]  = froms[index]

    return matrix,moves

def backtrack(s,t,matrix,moves):
    i     = len(s)
    j     = len(t)
    score = matrix[i][j]
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

    return score,s1[::-1],t1[::-1]

def align(s,t,replace_score=createSimpleDNASubst(),indel_cost=1,build_matrix=build_matrix,backtrack=backtrack):
    distances = create_distance_matrix(len(s)+1,len(t)+1)
    distances,moves = build_matrix(s,t,distances,replace_score=replace_score,indel_cost=indel_cost)
    return backtrack(s,t,distances,moves)

if __name__=='__main__':
    from Bio.SubsMat.MatrixInfo import blosum62
    score,s1,s2=align('PLEASANTLY','MEANLY',replace_score=blosum62,indel_cost=5)
    print (score,s1,s2)