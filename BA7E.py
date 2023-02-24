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

# BA7E Implement the Neighbor Joining Algorithm

from rosalind import Tree,read_matrix

def NeighborJoining(D,n,node_list=None):

    def remove(i,D):
        D_new=[]
        for j in range(len(D)):
            if j!=i:
                D_new.append([D[j][k] for k in range(len(D[j])) if k!=i])
        return D_new

    def create_DPrime(Total_distance):
        DPrime=[[0]*n for dummy in range(n)]
        for i in range(n):
            for j in range(i+1,n):
                DPrime[i][j]=(n-2)*D[i][j]-Total_distance[i]-Total_distance[j]
                DPrime[j][i]= DPrime[i][j]
        return DPrime

    def get_ij_minDPrime(DPrime):
        i=-1
        j=-1
        minD=float('inf')
        for ii in range(n):
            for jj in range(ii,n):
                if  DPrime[ii][jj]<minD:
                    i=ii
                    j=jj
                    minD= DPrime[i][j]
        return (i,j,minD)

    def createDelta(Total_distance):
        Delta=[[0]*n for dummy in range(n)]
        for i in range(n):
            for j in range(i+1,n):
                Delta[i][j]=(Total_distance[i]-Total_distance[j])/(n-2)
                Delta[j][i]=Delta[i][j]
        return Delta
    if node_list==None:
        node_list=list(range(n))

    if n==2:
        T=Tree()
        T.link(node_list[0],node_list[1],D[0][1])
        return T
    else:
        Total_distance=[sum(row) for row in D]
        DPrime=create_DPrime(Total_distance)
        (i,j,minD)=get_ij_minDPrime(DPrime)
        Delta=createDelta(Total_distance)
        limbLength_i=(D[i][j]+Delta[i][j])/2
        limbLength_j=(D[i][j]-Delta[i][j])/2
        new_row=[0.5*(D[k][i]+D[k][j]-D[i][j]) for k in range(n)]+[0]

        D.append(new_row)
        for l in range(n):
            D[l].append(new_row[l])

        m=node_list[-1]+1
        node_list.append(m)
        D=remove(max(i,j),D)
        D=remove(min(i,j),D)
        node_i=node_list[i]
        node_j=node_list[j]
        node_list.remove(node_i)
        node_list.remove(node_j)
        T=NeighborJoining(D,n-1,node_list)
        T.link(node_i,m,limbLength_i)
        T.link(node_j,m,limbLength_j)
        return T

if __name__=='__main__':
    params,D=read_matrix('c:/Users/Weka/Downloads/rosalind_ba7e.txt',conv=float)
    NeighborJoining(D,params[0]).print()
