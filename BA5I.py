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

# BA5I Find a Highest-Scoring Overlap Alignment of Two Strings

from align import create_distance_matrix
from numpy import argmax

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
            #print (distances[i])
        
        i = len(v)
        j = argmax(distances[i])
        distance = distances[i][j]
        v1 = []
        w1 = []
        while i>0 and j>0:
            i1,j1 = path[(i,j)]
            v1.append(v[i1] if i1<i else '-')
            w1.append(w[j1] if j1<j else '-')
            i,j=i1,j1
    
        return distance,v1[::-1],w1[::-1]
    
    score,u1,v1=dynamic_programming([vv for vv in v],[ww for ww in w])
    return score,''.join(u1),''.join(v1)

if __name__=='__main__':

    strings = []
 
    with open(r'C:\Users\Simon\Downloads\rosalind_ba5i(1).txt','r') as f:
        for line in f:
            strings.append(line.strip())

    d,s1,t1 = overlap_assignment(strings[0],strings[1])      
    print ('{0}'.format(d))
    print (s1)
    print (t1)    
    