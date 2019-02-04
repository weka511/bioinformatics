'''
 Copyright (C) 2019 Greenweaves Software Limited

 This is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This software is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>
 
 Programs written for 2nd week of Bioinformatics: Introduction and Methods

'''

import sys

def align(x,
          y,
          s=[[ 2, -7, -5, -7],
             [-7,  2, -7, -5],
             [-5, -7,  2, -7],
             [-7, -5, -7,  2]],
          d=-5,
          bound=-sys.maxsize):
        
    def index(x):
        return 'ACGT'.index(x)
    
    def dp():
        F = [[0 for j in range(len(y)+1)] for i in range(len(x)+1)]
        for i in range(1,len(F)):
            F[i][0] = F[i-1][0] + d
        for j in range(1,len(F[0])):
            F[0][j] = F[0][j-1] + d    
        for i in range(1,len(F)):
            for j in range(1,len(F[0])):
                F[i][j] = max(F[i-1][j-1] + s[index(x[i-1])][index(y[j-1])],
                              F[i-1][j] + d,
                              F[i][j-1] + d,
                              bound)
        return F

    def traceback(F):
        i    = len(F)-1
        j    = len(F[0])-1
        path = []
        while True:
            predecessors = []
            indices      = []
            xs           = []
            ys           = []
 
            if i>0 and j>0:             
                predecessors.append(F[i-i][j-1]+s[index(x[i-1])][index(y[j-1])])
                indices.append((i-1,j-1))
                xs.append(x[i-1])
                ys.append(y[j-1]) 
            if i>0:
                predecessors.append(F[i-1][j]+d)
                indices.append((i-1,j))
                xs.append(x[i-1])
                ys.append('-')
            if j>0:
                predecessors.append(F[i][j-1]+d)
                indices.append((i,j-1))
                xs.append('-')
                ys.append(y[j-1])  
                
            m   = predecessors.index(max(predecessors))
            i,j = indices[m]
            path.append((i,j,xs[m],ys[m]))
            if i==0 and j==0:
                return path
     
    F = dp()
    for f in F:
        print (f)
        
    t = traceback(F)
    for (_,_,x,y) in t:
        print (x,y)
    return (F,t,''.join([x for (_,_,x,_) in t]),''.join([y for (_,_,Y,y) in t]))

if __name__=='__main__':
    F,path,s1,s2 = align('AGCT','AAGT')
    for f in F:
        print (f)
    print (path)
    print (s1)
    print (s2)
