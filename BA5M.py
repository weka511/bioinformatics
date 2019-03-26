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

# BA5M.py Find a Highest-Scoring Multiple Sequence Alignment 


from Bio.SubsMat.MatrixInfo import blosum62
from numpy import argmax

def FindHighestScoringMultipleSequenceAlignment (u,v,w,score=lambda x,y,z: 1 if x==y and y==z else 0):
    s = [[[0 for i in range(len(w)+1)] for j in range(len(v)+1)] for k in range(len(u)+1) ]
    moves = {}
    for i in range(1,len(u)+1):
        for j in range(1,len(v)+1):
            for k in range(1,len(w)+1):
                scores     = [
                    s[i-1][j][k]     + score(u[i-1], '-',    '-'),
                    s[i][j-1][k]     + score('-',    v[j-1], '-'),
                    s[i][j][k-1]     + score('-',    '-',    w[k-1]),
                    s[i][j-1][k-1]   + score('-',    v[j-1], w[k-1]),
                    s[i-1][j][k-1]   + score(u[i-1], '-',    w[k-1]),
                    s[i-1][j-1][k]   + score(u[i-1], v[j-1], '-'),
                    s[i-1][j-1][k-1] + score(u[i-1], v[j-1], w[k-1])
                ]
                possible_moves = [
                    (-1,  0,   0),
                    (0,  -1,   0),
                    (0,  0,   -1),
                    (0,  -1,  -1),
                    (-1, 0,  -1),
                    (0,  -1, -1),
                    (-1, -1, -1),
                ]
                index          = argmax(scores)
                s[i][j][k]     = scores[index]
                moves[(i,j,k)] = possible_moves[index]
    i  = len(u)
    j  = len(v)
    k  = len(w)
    u1 = []
    v1 = []
    w1 = []
    while i>0 and j>0 and k>0:
        di,dj,dk = moves[(i,j,k)]
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
    while j>0:
        j-=1
        v1.append(v[j])
    while k>0:
        k-=1
        w1.append(w[k])    
    return s[len(u)][len(v)][len(w)],''.join(u1[::-1]),''.join(v1[::-1]),''.join(w1[::-1])             
            

if __name__=='__main__':
    from helpers import create_strings    
    s,u,v,w = FindHighestScoringMultipleSequenceAlignment('ATATCCG','TCCGA','ATGTACTG')
    print (s)
    print (u)
    print (v)
    print (w)
    s1,u1,v1,w1 = FindHighestScoringMultipleSequenceAlignment('TGTTTAAAAATGTCCGCAACCATTTC',
                                                              'GATATAAAACAGGGATAACTGCAATGG',
                                                              'CCTGCTACTTTATGCCGTCTCCATATGCG')
    print (s1)
    print (u1)
    print (v1)
    print (w1)    
 