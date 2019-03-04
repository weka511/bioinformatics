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
#
# fragile.py # code for Chapter 6 - Are there fragile regions in the human genome?

# BA6B Compute the Number of Breakpoints in a Permutation 

def getBreakPoints(P):
    return len( [(a,b) for (a,b) in zip([0]+P,P+[len(P)+1]) if b!=a+1] )

# BA6C Compute the 2-Break Distance Between a Pair of Genomes

def get_synteny_blocks(a):
    return [j for i in a for j in i] 

def count_synteny_blocks(a):
    return max(abs(b) for b in get_synteny_blocks(a))

# BA6F Implement Chromosome to Cycle

def ChromosomeToCycle(Chromosome):
    Nodes = []
    for i in Chromosome:
        if i> 0:
            Nodes.append(2*i-1)
            Nodes.append(2*i)
        else:
            Nodes.append(-2*i)
            Nodes.append(-2*i-1)
        
    return Nodes

# BA6H Implement ColoredEdges

def ColouredEdges(P):
    Edges = []
    for  Chromosome in P:
        Nodes = ChromosomeToCycle(Chromosome)
        it = iter(Nodes[1:]+[Nodes[0]])
        for i in it:
            Edges.append((i,next(it)))
    return Edges

