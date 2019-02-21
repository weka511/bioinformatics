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

def ba5h(s,t,indel_cost=1,replace_cost=lambda a,b: 1,show_scores=False):
      
      def dynamic_programming(s,t):
            score=[[0 for j in range(len(t)+1)] for i in range(len(s)+1)]

            for j in range(len(t)+1):
                  score[0][j]=j
            for i in range(len(s)+1):
                  score[i][0]=i

            for i in range(1,len(s)+1):
                  for j in range(1,len(t)+1):
                        score[i][j] = min(
                              score[i-1][j]   + indel_cost,
                              score[i][j-1]   + indel_cost,
                              score[i-1][j-1] + (0 if s[i-1]==t[j-1] else replace_cost(s[i-1],t[j-1])))
                        
            if show_scores:
                  print (' ',['#']+t)
                  for i in range(len(s)+1):
                        print ((['#']+s)[i] ,score[i])
                  print (' ',['#']+t) 
                  
            return 0,[],[]
      
      d,s1,t1=dynamic_programming([s0 for s0 in s], [t0 for t0 in t])
      return d,''.join(s1),''.join(t1)

if __name__=='__main__':
      d,s1,t1 = ba5h('GTAGGCTTAAGGTTA','TAGATA',show_scores=True)
      print (d)
      print (s1)
      print (t1)
