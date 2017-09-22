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

# BA7C Implement Additive Phylogeny
from BA7B import ComputeLimbLength

def AdditivePhylogeny(D,n):
    
    def find_ik():
        '''
        Find three leaves such that Di,k = Di,n + Dn,k
        '''
        for i in range(n-1):
            for k in range(n-1):
                if D[i][k]==D[i][n-1]+D[n-1][k]:
                    return (i,k)
        print ('not found')
        
    def Node(x,i,k):
        '''
        the (potentially new) node in T at distance x from i on the path between i and k
        '''
        candidates=[l for l in range(n) if D[i][l]==x]
        
        if len(candidates)==1:
            return candidates[0]
    
    def AddNode(v,limbLength):
        '''
         add leaf n back to T by creating a limb (v, n) of length limbLength
         '''
        T[v]=(n,limbLength)
    
            
    if n==2:
        return {0:(1,D[0][1])}
    else:
        limbLength=ComputeLimbLength(n-1,n-1,D)
        for j in range(n-1):
            D[n-1][j]-=limbLength
            D[j][n-1]=D[n-1][j]
        i,k=find_ik()
        x=D[i][n-1]
        T=AdditivePhylogeny(D,n-1)
        print('i={0},n={1},k={2},x={3}'.format(i,n,k,x),D)
        v=Node(x,i,k)
        print ('v={0}'.format(v))
        AddNode(v,limbLength)
        return T
    
if __name__=='__main__':
    n=4
    D=[[0,   13,  21,  22],
        [13,  0 ,  12,  13],
        [21,  12 , 0 ,  13],
        [22 , 13 , 13,  0]]
    print (AdditivePhylogeny(D,n))