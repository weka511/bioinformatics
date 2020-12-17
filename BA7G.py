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

# BA7G Adapt SmallParsimony to Unrooted Trees  http://rosalind.info/problems/ba7g/


from rosalind import LabelledTree
from phylogeny import SmallParsimony
from BA7F import print_assignments

def AdaptSmallParsimonyToUnrootedTrees(N,T):
    '''
    When the position of the root in an evolutionary tree is unknown,
    we can simply assign the root to any edge that we like, apply
    SmallParsimony from "Implement SmallParsimony" to the resulting rooted tree, and then remove the root.
    It can be shown that this method provides a solution to the following problem.
    Small Parsimony in an Unrooted Tree Problem
    '''
    def assign_root():
        '''
        Assign a root to the tree.
        
        Initially I followed Igor Segota's solution from
        https://stepik.org/lesson/10335/step/12?course=Stepic-Interactive-Text-for-Week-3&unit=8301,
        but found that a random root  generally led to problems with the Small Parsimony algorithm.
        Using the last node as one half of the broken lenk works well.
        '''
        a=T.nodes[len(T.nodes)-1]
        b,_=T.edges[a][0]
        print ("Breaking at {0} {1}".format(a,b))
        T.unlink(a,b)
        c=T.next_node()
        T.link(c,a)
        T.link(c,b)
        return (a,b,c)
    
    a,b,root=assign_root()
    
    T.remove_backward_links(root)

    return a,b,root,T


                        
if __name__=='__main__':
    def parse(input):
        N=-1
        lines=[]
        for line in input:
            if N==-1:
                N=int(line.strip())
            else:
                lines.append(line.strip())
        return N,LabelledTree.parse(N,lines,bidirectional=True)
    
 
    with open('c:/Users/Weka/Downloads/rosalind_ba7g(3).txt') as f:
        N,T=parse(f)
        a,b,root,T1= AdaptSmallParsimonyToUnrootedTrees(N,T)
        score,assignments=SmallParsimony(T1)
        # This is the fixup at the emd of processing
        assignments.unlink(root,b)     
        assignments.unlink(root,a)
        assignments.link(a,b)
        
        print (score)
        print_assignments(assignments)