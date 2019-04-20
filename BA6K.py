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
from fragile import GraphToGenome,ChromosomeToCycle,CycleToChromosome

def GraphToGenomeAll(GenomeGraph):
     def build_index():
          def insert(a,b):
               if not a in Index:
                    Index[a]=[]
               Index[a].append(b)
          Index = {}
          for (a,b) in GenomeGraph:
               insert(a,b)
               insert(b,a)
          return Index
     
     def extract_cycle(Index,a,b):
          Cycle = [a,b]
          while True:    # build one cycle
               nexts = Index[b]
               next_link = list(set(nexts) - set(Cycle))
               if len(next_link)>0:
                    b = next_link[0]
                    Cycle.append(b)
               else:
                    x,y = nexts
                    assert x== Cycle[0] or y == Cycle[0]
                    for c in Cycle:
                         del Index[c]
                    return Cycle
                
     Index  = build_index()
     Graph  = GenomeGraph[:]
     Genome = []
     while len(Graph)>0:
          a,b   = Graph[0]
          Cycle = extract_cycle(Index,a,b)
          Genome.append(Cycle)
          Graph = [(a,b) for (a,b) in Graph if not a in Cycle]
     return [CycleToChromosome(g) for g in Genome]
     
def perform_2_BreakOnGenome(P,i0,i1,j0,j1):
     def Edges(it):
          return [(i,next(it)) for i in it]
     def BlackEdges():
          return Edges(iter(Nodes))     
     def ColouredEdges():
          return Edges(iter(Nodes[1:]+[Nodes[0]]))
     Nodes = ChromosomeToCycle(P)
     GenomeGraph =  BlackEdges() + ColouredEdges()
     GenomeGraph = get2BreakOnGenomeGraph(GenomeGraph,i0,i1,j0,j1)
     return GraphToGenomeAll(GenomeGraph)

if __name__=='__main__':
     print (perform_2_BreakOnGenome([+1, -2, -4, +3], 1, 6, 3, 8))
    