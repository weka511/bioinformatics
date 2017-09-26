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

from rosalind import Tree,read_matrix  
    
def AdditivePhylogeny(D,n,N=-1):
    def find_ikn(DD):
        '''
        Find three leaves such that Di,k = Di,n + Dn,k
        '''
        for i in range(n):
            for k in range(n):
                if DD[i][k]==DD[i][n-1]+DD[n-1][k] and i!=k:
                    return(i,k,n-1,DD[i][n-1])
    
    def get_Position_v(traversal):
        d=0
        for l,w in traversal:
            d0=d
            d+=w
            if d==x: return (True,l,l,d0,d)
            if d>x: return (False,l_previous,l,d0,d)
            l_previous=l
            
        return (False, l_previous, l, d0,d)
        
    if N==-1:
        N=n
        
    if n==2:
        T=Tree(N)
        T.link(0,1,D[0][1])
        return T
    else:
        limbLength=ComputeLimbLength(n,n-1,D)
        
        D_bald=[d_row[:] for d_row in D]
        for j in range(n-1):   
            D_bald[n-1][j]-=limbLength
            D_bald[j][n-1]=D_bald[n-1][j]  
         
        i,k,node,x=find_ikn(D_bald)
        #x=D_bald[i][n-1]
        
        D_Trimmed=[D_bald[l][:-1] for l in range(n-1)]
        
        T=AdditivePhylogeny(D_Trimmed,n-1,N)
        # v= the (potentially new) node in T at distance x from i on the path between i and 
        found_k,traversal=T.traverse(i,k)
        path,weights=zip(*traversal)

        found,l0,l1,d,d0=get_Position_v(traversal)
         
        if found:
            v=l0  #Untested!
            T.link(node,v,limbLength)
        else:
            v=T.next_node()
            weight_i=ComputeLimbLength(n,i,D)
            weight_k=ComputeLimbLength(n,k,D)
            T.unlink(l0,l1)
            T.link(v,l0,x-d)
            T.link(v,l1,d0-x)
        # add leaf n back to T by creating a limb (v, n) of length limbLength
        T.link(node,v,limbLength)
        
        return T
 


if __name__=='__main__':
    params,D=read_matrix('c:/Users/Weka/Downloads/rosalind_ba7c(8).txt')           
    AdditivePhylogeny(D,params[0]).print()
    
