# Copyright (C) 2019-2021 Greenweaves Software Limited

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

# BA6D Find a Shortest Transformation of One Genome into Another by 2-Breaks 

import argparse
import os
import time
from   helpers import read_strings

from fragile import ChromosomeToCycle,ColouredEdges,BlackEdges,get2BreakOnGenomeGraph,CycleToChromosome

# Find a shortest transformation of one genome into another by 2-breaks.

# Given: Two genomes with circular chromosomes on the same set of synteny blocks.

# Return: The sequence of genomes resulting from applying a shortest sequence of
#         2-breaks transforming one genome into the other.
def FindShortestTransformation(s,t,N=25,M=10):
     def mismatches(s,t):
          return sum([0 if a==b else 1 for (a,b) in zip(sorted(s),sorted(t))])
     def get2_breaks(Configuration):
          result = []
          for k in range(len(Configuration)):
               for l in range(k):
                    i0,j0 = Configuration[k]
                    i1,j1 = Configuration[l]
                    result.append((i0,j0,i1,j1))
          return result
     
     def FindShortestTransformationCycles(s,t):
          assert sorted(s)== sorted(t)
          ColouredS     = ColouredEdges(s)
          Blacks        = BlackEdges(s)
          ColouredT     = ColouredEdges(t)
          leader_board = [(ColouredS,mismatches(ColouredS,ColouredT),[])] #Configuration, score
          for _ in range(M):
               new_leaders = []
               for Configuration,_,path in leader_board:
                    for i0,i1,j0,j1 in get2_breaks(Configuration):
                         transformed = get2BreakOnGenomeGraph(Blacks + Configuration,i0,i1,j0,j1)
                         Coloured    = [tt for tt in list(set(transformed)) if not tt in Blacks]
                         new_leaders.append( (Coloured, mismatches(Coloured,ColouredT),path+[Coloured]))
               leader_board = sorted(new_leaders,key=lambda tuple:tuple[1])
               if len(leader_board)>N:
                    leader_board = leader_board[:N]
               #Config,score = leader_board[0]
               #print (Config,score)
               Config,score,path = leader_board[0]
               if score==0:
                    return path,Blacks
     
     return FindShortestTransformationCycles(ChromosomeToCycle(s),ChromosomeToCycle(t))

def CycleToChromosome1(Blacks,Coloured):
     Chromosome = []
     flattened = [y for x in Coloured for y in x]
     for a,b in Blacks:
          x = flattened.index(a)
          if x%2==0:
               a,b=b,a
          if a<b:
               Chromosome.append(b//2)
          else:
               Chromosome.append(-a//2)

     return Chromosome

if __name__=='__main__':
     start = time.time()
     parser = argparse.ArgumentParser('Find a Shortest Transformation of One Genome into Another by 2-Breaks')
     parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
     parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
     args = parser.parse_args()
     if args.sample:
          paths,Blacks = FindShortestTransformation([+1, -2, -3, +4],[+1, +2, -4, -3])
          #print (len(paths))
          for Coloured in paths:
               #print (Coloured)
               print (CycleToChromosome1(Blacks,Coloured))          



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
     print (f'Elapsed Time {minutes} m {seconds:.2f} s')  
  