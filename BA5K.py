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

# BA5K.py Find a Middle Edge in an Alignment Graph in Linear Space

from Bio.SubsMat.MatrixInfo import blosum62
from numpy import argmax

def FindMiddleEdge(s,t,replace_score=blosum62,indel_cost=5):
    def score(pair):
        def reverse(pair):
            a,b=pair
            return (b,a)
        return replace_score[pair] if pair in replace_score else replace_score[reverse(pair)]     
    def update(j,col1,col2):
        col2[0] = j*indel_cost
        for i in range(1,len(s)+1):
            scores = [col1[i-1] + score((s[i-1],t[j-1])),
                      col2[i-1] + indel_cost,
                      col1[i]   + indel_cost
            ]
 
            best        = argmax(scores)
            col2[i]     = scores[best]
   
    column_A = [i*indel_cost for i in range(len(s)+1)]
    column_B = [0 for i in range(len(s)+1)]
    for j in range(1,len(t)//2):
        if j%2 ==1:
            update(j,column_A,column_B)
        else:
            update(j,column_B,column_A)
        
    return ((4, 3), (5, 4))

if __name__=='__main__':
    from helpers import create_strings    
    print (FindMiddleEdge('PLEASANTLY','MEASNLY'))