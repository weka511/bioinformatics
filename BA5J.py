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

from align import create_distance_matrix
from Bio.SubsMat.MatrixInfo import blosum62

def san_kai(s,t, replace_score=blosum62,sigma=11,epsilon=1):
    def score(pair):
        def reverse(pair):
            a,b=pair
            return (b,a)
        return replace_score[pair] if pair in replace_score else replace_score[reverse(pair)]     
    lower  =  create_distance_matrix(len(s)+1,len(t)+1)
    middle =  create_distance_matrix(len(s)+1,len(t)+1)
    upper  =  create_distance_matrix(len(s)+1,len(t)+1)
    for i in range(len(s)+1):
        for j in range(len(t)+1):
            #if i==0 and j==0:
                #pass
            #elif i==0:
            if i==0:
                lower[i][j] = - (sigma + epsilon *(j-1))
                upper[i][j] = - (sigma + epsilon *(j-1)) 
 
            elif j==0:
                lower[i][j] = - (sigma + epsilon *(i-1))
                upper[i][j] = - (sigma + epsilon *(i-1))

            else:
                lower[i][j]  = max(lower[i-1][j]    - epsilon,
                                   middle[i-1][j]   - sigma)
                middle[i][j] = max(lower[i][j],
                                   middle[i-1][j-1] + score((s[i-1],t[i-1])) ,
                                   upper[i][j])     
                upper[i][j]  = max(upper[i][j-1]    - epsilon,
                                   middle[i][j-1]   - sigma)
        
    return middle[len(s)][len(t)],[],[]

def ba5j(s,t):
    score,s1,t1 = san_kai([s0 for s0 in s],[t0 for t0 in t])
    return score,''.join(s1),''.join(t1)

if __name__=='__main__':
    print(ba5j('PRTEINS','PRTWPSEIN'))