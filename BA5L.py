#    Copyright (C) 2019-2020 Greenweaves Software Limited
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
#    along with this program.  If not, see <http://www.gnu.org/licenses/>

#    BA5L Align Two Strings Using Linear Space

import argparse
import os
import time
from   helpers import read_strings
from   align import FindMiddleEdge
from   Bio.Align import substitution_matrices

# alignUsingLinearSpace
#
# Align Two Strings Using Linear Space
#
# Inputs: v
#         w
#         replace_score
#         indel_cost

def alignUsingLinearSpace(v,w,
                          replace_score = substitution_matrices.load("BLOSUM62"),
                          indel_cost    = 5):
    
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
    
    # LinearSpaceAlignment
    #
    # Find longest path between a substring of v[top] v[bottom-1]
    # and w[left] and w[right-1]
    #
    # Inputs: top
    #         bottom
    #         left
    #         right
    
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
    LinearSpaceAlignment(0,len(v)+1,0,len(w)+1)
        
    
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA5L.py Align Two Strings Using Linear Space')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',     default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (alignUsingLinearSpace('PLEASANTLY','MEANLY'))
     
    if args.extra:
        Input,Expected             = read_strings(f'data/linear_space_alignment.txt',init=0)
        print (alignUsingLinearSpace(Input[0],Input[1]))
        
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = None
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes    
