#!/usr/bin/env python

# Copyright (C) 2017 Greenweaves Software Pty Ltd

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

import numpy as np

def soft_kmeans(k,m,beta,points,N=1,centres=None):
    def distance(p1,p2):
        return np.sqrt(sum((p1[i]-p2[i])**2 for i in range (m)))
    @np.vectorize
    def hidden_matrix(i,j):
        return np.exp(-beta*distance(centres[i],points[j]))
    def step(centres):
        ks=np.array(range(k))
        ns=np.array(range((len(points))))
        kk,mm=np.meshgrid(ns,ks,indexing='xy')
        numerator=hidden_matrix(mm,kk)
        denominator=np.sum(numerator,axis=0)
        matrix=np.divide(numerator,denominator)
        new_centres=[[] for i in range(k)]
        for i in range(k):
            for j in range(m):
                x_i_j=sum(matrix[i,l]*points[l][j] for l in range(len(points)))/sum(matrix[i,l] for l in range(len(points)))
                new_centres[i].append(x_i_j)
        return new_centres

    if centres==None:
        centres=points[:k]
    for i in range(N):
        centres=step(centres)
    return centres

if __name__=='__main__':
    m = -1
    k = -1
    beta = -1
    points=[]

    with open (r'C:\Users\Weka\Downloads\rosalind_ba8d.txt') as f:
    #with open('BA8D.txt') as f:
        for line in f:
            if k==-1:
                values=line.strip().split()
                k=int(values[0])
                m=int(values[1])
            elif beta==-1:
                beta=float(line.strip())
            else:
                points.append([float(v) for v in line.strip().split()])
        for pt in soft_kmeans(k,m,beta,points,N=100):
            print (' '.join('{0:.3f}'.format(p) for p in pt))
