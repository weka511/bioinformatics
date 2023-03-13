#!/usr/bin/env python

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
#    BA5H Find a Highest-Scoring Fitting Alignment of Two Strings

from reference_tables import createSimpleDNASubst
from align import align,get_highest_scoring_global_alignment
import numpy as np
from helpers import create_strings

def build_matrix(s,t,matrix,replace_score=createSimpleDNASubst(),indel_cost=1,get_indel_cost=None):
      moves = {}
      def score(pair):
            def reverse(pair):
                  a,b=pair
                  return (b,a)
            return replace_score[pair] if pair in replace_score else replace_score[reverse(pair)]
      for i in range(len(s)+1):
            for j in range(len(t)+1):
                  if i==0 and j==0: pass
                  elif i==0:
                        matrix[i][j] = 0
                        moves[(i,j)] = (0,0,0,-1)
                  elif j==0:
                        matrix[i][j] = 0
                        moves[(i,j)]  =(0,0,-1,0)
                  else:
                        scores       = [matrix[i-1][j]   - indel_cost,
                                        matrix[i][j-1]   - indel_cost,
                                        matrix[i-1][j-1] + score((s[i-1],t[j-1]))]
                        froms        = [(i-1, j,   -1,  0),
                                        (i,   j-1,  0, -1),
                                        (i-1, j-1, -1, -1)]
                        index        = np.argmax(scores)
                        matrix[i][j] = scores[index]
                        moves[(i,j)] = froms[index]

      return matrix,moves

def backtrack(s,t,matrix,moves,showPath=False):

      score = max([matrix[i][-1] for i in range(len(s)+1)])
      i     = -1
      j     = len(t)
      for k in range(len(s)-1,-1,-1):
            if matrix[k][-1]==score:
                  i = k
                  break
      s1    = []
      t1    = []
      while i>0 or j>0:
            i,j,di,dj = moves[(i,j)]
            if di==0:
                  s1.append('-')
                  t1.append(t[j])
            elif dj==0:
                  s1.append(s[i])
                  t1.append('-')
            else:
                  s1.append(s[i])
                  t1.append(t[j])

      return score,s1[:-1][::-1],t1[:-1][::-1]


def ba5h(s,t):
      d,s1,t1 =align([s0 for s0 in s], [t0 for t0 in t],build_matrix=build_matrix,backtrack=backtrack)
      return (d,''.join(s1),''.join(t1))


if __name__=='__main__':
      import timeit
      start_time = timeit.default_timer()
      strings   = create_strings(ext=1)
      d,s1,t1 = ba5h(strings[0],strings[1])
      print ('Score = {0}'.format(d))
      print (s1)
      print (t1)
      print ('Elapsed Time = {0}'.format(timeit.default_timer() - start_time))
