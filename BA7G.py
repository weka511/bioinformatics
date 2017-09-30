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

# BA7G Adapt SmallParsimony to Unrooted Trees  

#If anyone is stuck, Igor Segota solved this without changing any code from the previous challenge, only the input to it.

#I did this, as pointed out, by inserting a new node between a random pair of nodes
# (store this old pair of nodes as they'll be needed later). Then I removed old two
# (same but opposite) edges, replaced them with two new edges going away from the root.

#Then I did the breadth-first scan through the graph and deleted all 'backwards' facing edges. 
#This nice tree graph is now all good and ready for small_parsimony().
#After the small_parsimony is done (for all characters from sequences), I remove any reference to the root node 
# (I forgot to do this initially and was stuck). Then I add back an edge (*) between the old two nodes 
# with the corresponding hamming distance from these two removed nodes.

#(*) actually two because the output requires every edge in duplicate, going the  other way

from rosalind import LabelledTree

def AdaptSmallParsimonyToUnrootedTrees(N,T):
    pass

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
    
    input=[
        '4',
        'TCGGCCAA->4',
        '4->TCGGCCAA',
        'CCTGGCTG->4',
        '4->CCTGGCTG',
        'CACAGGAT->5',
        '5->CACAGGAT',
        'TGAGTACC->5',
        '5->TGAGTACC',
        '4->5',
        '5->4'    
    ]
    N,T=parse(input)
    T.print()
    print (AdaptSmallParsimonyToUnrootedTrees(N,T))
