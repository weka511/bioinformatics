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

# BA7D Implement UPGMA

from rosalind import Tree,read_matrix,DPrint 

def UPGMA(D, n):
    def find_two_closest_clusters():
        ii=-1
        jj=-1
        best_distance=float('inf')
        for i in range(len(D)):
            for j in range(i):
                if i in Clusters and j in Clusters and D[i][j]<best_distance:
                    ii=i
                    jj=j
                    best_distance=D[i][j]
        return (ii,jj,best_distance)    
    Clusters={}
    for i in range(n):
        Clusters[i]=[i]

    T=Tree(n)
    Age={}
    for node in T.get_nodes():
        Age[node]=0
    while len(Clusters)>1:
        def d(i,j):
            return sum([D[cl_i][cl_j] for cl_i in Clusters[i] for cl_j in Clusters[j]])/(len( Clusters[i])* len(Clusters[j])) \
                   if i in Clusters and j in Clusters \
                   else float('nan')
  
        i,j,distance=find_two_closest_clusters()
        node=T.next_node()

        T.link(node,i)
        T.link(node,j)        
        Clusters[node]=Clusters[i]+Clusters[j]
        Age[node]=D[i][j]/2
        del Clusters[i]
        del Clusters[j]

        row=[d(i,node) for i in range(len(D))] + [0.0]
        for k in range(len(D)):
            D[k].append(row[k])
        D.append(row)

    for node in T.nodes:
        T.edges[node]=[(e,abs(Age[node]-Age[e])) for e,W in T.edges[node]]
        
    return T

if __name__=='__main__':
    params,D=read_matrix('c:/Users/Weka/Downloads/rosalind_ba7d(1).txt',conv=float) 
    UPGMA(D, params[0]).print()


