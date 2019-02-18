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
# EDTA Edit Distance Alignment http://rosalind.info/problems/edta/

from edit import edit
import numpy as np

def edta(s,t,indel_cost=1,replace_cost=lambda a,b: 1):
    def extract(s,t,matrix):
        m  = len(matrix)-1
        n  = len(matrix[0])-1
        s1 = []
        t1 = []
        while m>0 and n>0:
            moves  = [(m-1,n),(m,n-1),(m-1,n-1)]
            scores = [matrix[m-1][n]+indel_cost,
                      matrix[m][n-1]+indel_cost,
                      matrix[m-1][n-1] + (0 if s[m-1]==t[n-1] else replace_cost(s[m-1],t[n-1]))]
            ss     = [s[m-1],'-',s[m-1]]
            ts     = ['-',t[n-1],t[n-1]]
            index  = np.argmin(scores)
            m,n    = moves[index]
            s1.append(ss[index])
            t1.append(ts[index])
        s1.reverse()
        t1.reverse()
        return ''.join(s1),''.join(t1)
    d,matrix = edit(s,t,indel_cost,replace_cost)
    s1,t1    = extract([s0 for s0 in s], [t0 for t0 in t],matrix)
    return (d,s1,t1)

if __name__=='__main__':
    print (edta('PRETTY','PRTTEIN'))
