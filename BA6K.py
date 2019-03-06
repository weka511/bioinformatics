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
#    BA6K Implement 2-BreakOnGenome

from BA6J import get2BreakOnGenomeGraph
from fragile import GraphToGenome,ChromosomeToCycle

def ColouredEdges(P):
     Edges = []
     Nodes = ChromosomeToCycle(P)
     it = iter(Nodes[1:]+[Nodes[0]])
     for i in it:
          Edges.append((i,next(it)))
     return Edges

def perform_2_BreakOnGenome(P,i0,i1,j0,j1):
     def BlackEdges(P):
          Nodes = ChromosomeToCycle(P)
          Edges = []
          it = iter(Nodes[1:]+[Nodes[0]])
          for i in it:
               Edges.append((i,next(it)))          
          return Edges
     GenomeGraph = ColouredEdges(P) # BlackEdges(P) +
     GenomeGraph = get2BreakOnGenomeGraph(GenomeGraph,i0,i1,j0,j1)
     return GraphToGenome(GenomeGraph)

if __name__=='__main__':
     print (perform_2_BreakOnGenome([+1, -2, -4, +3], 1, 6, 3, 8))
    