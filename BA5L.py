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

# BA5L.py Align Two Strings Using Linear Space

from align import FindMiddleEdge
from Bio.SubsMat.MatrixInfo import blosum62

# alignUsingLinearSpace
#
# Align Two Strings Using Linear Space
#
# Inputs: v
#         w

def alignUsingLinearSpace(v,w,replace_score=blosum62,indel_cost=5):
    # LinearSpaceAlignment
    #
    # Find longest path between a substring of v[top] v[bottom-1]
    # and w[left] and w[right-1]
    #
    # Inputs: top
    #         bottom
    #         left
    #         right
    
    def isRightOrDownRight(midEdge):
        return midEdge==RIGHT or midEdge==DOWNRIGHT
    
    def isDownOrDownRight(midEdge):
        return midEdge==DOWN or midEdge==DOWNRIGHT
    
    # MiddleNodeAndEdge
    #
    # An adapter which replaces MiddleNode and MiddleEdge in the pseudocode, and calls FindMiddleEdge
    def MiddleNodeAndEdge(top, bottom, left, right):
        ((i1,j1),(i2,j2)) = FindMiddleEdge(v[top:bottom],w[left:right],replace_score=replace_score,indel_cost=indel_cost)
        direction         = RIGHT if i1==i2 else DOWN if j1==j2 else DOWNRIGHT
        return j1,direction
    
    
    def  LinearSpaceAlignment(top, bottom, left, right):
        if left==right:
            return indel_cost*(bottom - top)
        if top==bottom:
            return indel_cost*(right-left)
        middle           = (left + right)//2
        midNode,midEdge  = MiddleNodeAndEdge(top, bottom, left, right)
  
        LinearSpaceAlignment(top, midNode, left, middle)
        # output midEdge
        if isRightOrDownRight(midEdge):
            middle += 1
        if isDownOrDownRight(midEdge):
            midNode+= 1
        LinearSpaceAlignment(midNode, bottom, middle, right) 
        
    RIGHT     = 0
    DOWN      = 1
    DOWNRIGHT = 2
    LinearSpaceAlignment(0,len(v),0,len(w))
    
if __name__=='__main__':
    from helpers import create_strings    
    print (alignUsingLinearSpace('PLEASANTLY','MEANLY'))
 