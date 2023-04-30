#!/usr/bin/env python

# Copyright (C) 2017-2023 Greenweaves Software Pty Ltd

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

''' BA7G Adapt SmallParsimony to Unrooted Trees  http://rosalind.info/problems/ba7g/'''


from rosalind import LabelledTree
from phylogeny import SmallParsimony, AdaptSmallParsimonyToUnrootedTrees
from BA7F import print_assignments


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
        # This is the fixup at the end of processing
        assignments.unlink(root,b)
        assignments.unlink(root,a)
        assignments.link(a,b)

        print (score)
        print_assignments(assignments)
