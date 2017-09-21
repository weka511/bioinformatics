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

# BA7A Compute Distances Between Leaves

def ComputeDistancesBetweenLeaves(n,T):
    def get_leaves():
        Froms = list(T.keys())
        Tos   = [a for v in list(T.values()) for (a,_) in v]
        print (Froms)
        print (Tos)
        X=Froms + Tos
        X.sort()
        Counts={}
        for x in X:
            if x in Counts:
                Counts[x]+=1
            else:
                Counts[x]=1
        return [node for node in Counts.keys() if Counts[node]==2]
    def D(i,j,path=[]):
        if i==j:
            return 0
        d=float('inf')
        for node,weight in T[i]:
            if node==j:
                return weight
            if node in path:
                continue
            test = weight + D(node,j,path+[node])
            if test<d:
                d=test
        return d
    Leaves=get_leaves()
    return [[D(i,j)for j in range(n)] for i in range(n) ]
    
if __name__=='__main__':
    T={
    0:[(4,11)],
    1:[(4,2)],
    2:[(5,6)],
    3:[(5,7)],
    4:[(0,11),(1,2),(5,4)],
    5:[(4,4),(3,7),(2,6)]
    }    
    print (ComputeDistancesBetweenLeaves(4,T))
