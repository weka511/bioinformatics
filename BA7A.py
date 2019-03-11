# Copyright (C) 2017-2019 Greenweaves Software Limited

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

# BA7A Compute Distances Between Leaves

# ComputeDistancesBetweenLeaves

# Inputs:  n an integer n 
#          T the adjacency list of a weighted tree with n leaves.
#
# Returns: An n by n symmetric matrix of distannces between leaves 

def ComputeDistancesBetweenLeaves(n,T):

    # D
    # Recursively compute the distances between two nodes
    def D(i,j,path=[]):
        if i==j:  return 0
        d = float('inf')
        for node,weight in T[i]:
            if node==j:
                return weight
            if node in path:
                continue
            test = weight + D(node,j,path+[node])
            if test<d:
                d=test
        return d

    return [[D(i,j)for j in range(n)] for i in range(n) ]
    
if __name__=='__main__':
    from helpers import create_weighted_adjacency_list
    n,T = create_weighted_adjacency_list()           
    for ds in ComputeDistancesBetweenLeaves(n,T):
        print(' '.join(str(d) for d in ds))
