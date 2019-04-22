# Copyright (C) 2019 Greenweaves Software Limited

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

# BA6D Find a Shortest Transformation of One Genome into Another by 2-Breaks 

from fragile import ChromosomeToCycle,ColouredEdges,BlackEdges,get2BreakOnGenomeGraph

def FindShortestTransformation(s,t,N=5,M=10):
     def mismatches(s,t):
          return sum([0 if a==b else 1 for (a,b) in zip(s,t)])
     def FindShortestTransformationCycles(s,t):
          print (s)
          print (t)
          print (mismatches(s,t))
          Nodes       = ChromosomeToCycle(s)
          check_nodes = ChromosomeToCycle(t)
          assert sorted(Nodes)== sorted(check_nodes)
          Coloured     = ColouredEdges(Nodes)
          Blacks       = BlackEdges(Nodes)
          ColouredT     = ColouredEdges(check_nodes)
          leader_board = [([],Coloured,mismatches(Coloured,ColouredT))] #2-breaks, 2-break(s), score
          for _ in range(M):
               print ('---')
               new_leaders = []
               for path,Coloured,_ in leader_board:
                    for k in range(len(Coloured)):
                         for l in range(k):
                              i0,i1 = Coloured[k]
                              j0,j1 = Coloured[l]
                              if len(path)>0:
                                   if (i0,i1,j0,j1)==path[-1]: continue
                              transformed = get2BreakOnGenomeGraph(Blacks + Coloured,i0,i1,j0,j1)
                              tt = sorted([t for t in transformed if not t in Blacks])
                              new_leaders.append((path+[(i0,i1,j0,j1)],
                                                  tt,
                                                  mismatches(tt,ColouredT)))
               leader_board = sorted(new_leaders)
               if len(leader_board)>N:
                    leader_board = leader_board[:N]
               path,_,score = leader_board[0]
               print (score,path)
          x=0
               
     return FindShortestTransformationCycles(ChromosomeToCycle(s),ChromosomeToCycle(t))

if __name__=='__main__':
     FindShortestTransformation([+1, -2, -3, +4],[+1, +2, -4, -3])