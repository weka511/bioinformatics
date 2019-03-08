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
from numpy import argmax

# https://web.stanford.edu/class/cs262/presentations/lecture2.pdf

def score(pair,replace_score=blosum62):
    def reverse(pair):
        a,b=pair
        return (b,a)
    return replace_score[pair] if pair in replace_score else replace_score[reverse(pair)]

def san_kai(s,t, replace_score=blosum62,sigma=11,epsilon=1):
    
    lower    = create_distance_matrix(len(s)+1,len(t)+1)
    middle   = create_distance_matrix(len(s)+1,len(t)+1)
    upper    = create_distance_matrix(len(s)+1,len(t)+1)
    #path  = create_distance_matrix(len(s)+1,len(t)+1)
    moves = {}
    lower[0][0]  = -float('inf')
    middle[0][0] = 0
    upper[0][0]  = -float('inf')
    for i in range(1,len(s)+1):
        lower[i][0]  = - (sigma + epsilon *(i-1))
        middle[i][0] = -float('inf')
        upper[i][0]  = -float('inf')
    for j in range(1,len(t)+1):
        lower[0][j]  = -float('inf')
        middle[0][j] = -float('inf')
        upper[0][j]  = - (sigma + epsilon *(j-1))
        
    for i in range(1,len(s)+1):
        for j in range(1,len(t)+1):
            lower[i][j]  = max(lower[i-1][j] - epsilon,
                               middle[i-1][j] - sigma)
            
            upper[i][j]  =  max(upper[i][j-1] - epsilon,
                               middle[i][j-1] - sigma)
            
            middle[i][j] = max(lower[i][j],
                               middle[i-1][j-1] + score((s[i-1],t[j-1])),
                               upper[i][j])
            
            #history      = []
            #F[i][j]      = V[i-1][j-1] + score((s[i-1],t[j-1]))
            #history.append((i-1,j-1,s[i-1],t[j-1]))
            
            #G[i][j]      = max(V[i-1][j] - sigma, 
                               #G[i-1][j] - epsilon)
            #history.append((i-1,j,s[i-1],'-'))
            
            #H[i][j]      = max(V[i][j-1] - sigma,
                               #H[i][j-1] - epsilon)
            #history.append((i,j-1,'-',t[j-1]))
            
            choices      = [lower[i][j], 
                            middle[i-1][j-1] + score((s[i-1],t[j-1])),
                            upper[i][j]]
            index        = argmax(choices)
            history=[(i-1,j,s[i-1],'-'),(i-1,j-1,s[i-1],t[j-1]),(i,j-1,'-',t[j-1])]
            #V[i][j]      = choices[index]
            moves[(i,j)] = history[index]
            #path[i][j]   = index
        #print (V[i])
        #print (path[i])
    i  = len(s)
    j  = len(t)
    ss = []
    ts = []

    while i>0 and j > 0:
        #direction = path[i][j]
        #if direction==0: #->i,->j
            #ss.append(s[i-1])
            #ts.append(t[j-1])
            #i-=1
            #j-=1
        #elif direction==1: #->i
            #ss.append(s[i-1])
            #ts.append('-')
            #i-=1
        #else:               #->j
            #ss.append('-')
            #ts.append(t[j-1])
            #j-=1
        i,j,s0,t0=moves[(i,j)]
        ss.append(s0)
        ts.append(t0)
    return middle[len(s)][len(t)],ss[::-1],ts[::-1]


def ba5j(s,t):
    score,s1,t1 = san_kai([s0 for s0 in s],[t0 for t0 in t])
    return score,''.join(s1),''.join(t1)

if __name__=='__main__':
    from helpers import create_strings    
    import os, os.path
    score,s,t=ba5j('PRTEINS','PRTWPSEIN')
    #score,s,t=ba5j('AHRQPQ','AHED')  # supposedly should align AHRQPQ AHE--D
    #score,s,t=ba5j('YHFDVPDCWAHRYWVENPQAIAQMEQICFNWFPSMMMKQPHVFKVDHHMSCRWLPIRGKKCSSCCTRMRVRTVWE',
                #'YHEDVAHEDAIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPIISATCARMRVRTVWE')    
#    144
#    YHFDVPDCWAHRYWVENPQAIAQM------E-QICFNWFPSMMMK-------QPHVF---KVDHHMSCRWLPIRGKKCSSCCTRMRVRTVWE
#    YHEDV----AH-----E-DAIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPI----ISATCARMRVRTVWE
# File has:
#    YHFDVPDCWAHRYWVENPQAIAQME-------QICFNWFPSMMMK-------QPHVFKV---DHHMSCRWLPIRGKKCSSCCTRMRVRTVWE
#    YHEDV----AHE------DAIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPI----ISATCARMRVRTVWE
#
# I tried replacing blossum62 to no avail

    #strings  = create_strings('ba5j',ext=4)
    #score,s,t=ba5j(strings[0],strings[1])
    print (score)
    print (s)
    print (t)        
        #with open('ba5j.txt','w') as o:
            #o.write('{0}\n'.format(score))
            #o.write('{0}\n'.format(s))
            #o.write('{0}\n'.format(t))