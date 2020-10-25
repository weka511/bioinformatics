#    Copyright (C) 2020 Simon Crase
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
# SIMS Finding Mutated Motifs

import numpy as np

def get_score(c1,c2, match=1,mismatch=-1):
    return match if c1==c2 else mismatch

def score_string(s,t):
    return sum(get_score(c1,c2) for (c1,c2) in zip(s,t))

def sims(s,t, match=1,mismatch=-1):
    scores    = np.zeros((len(t),len(s)),dtype=int)
    link_back = {}
    for j in range(len(s)):
        scores[0,j]  = get_score(s[j],t[0],match=match,mismatch=mismatch)
    for i in range(1,len(t)):
        for j in range(i,len(s)):
            score_diag       = scores[i-1][j-1] + get_score(s[j],t[i],match=match,mismatch=mismatch)
            score_horizontal = scores[i-1][j] + mismatch
            if score_horizontal>score_diag:
                scores[i,j] = score_horizontal
                link_back[(i,j)] = (i,j-1)
            else:
                scores[i,j] = score_diag
                link_back[(i,j)] = (i-1,j-1)
    for i in range(len(scores)):
        print (scores[i])
    print (max(scores[-1]))
    j0      = np.argmax(scores[-1])
    i0      = len(scores)-1
    s_match = [s[j0]]
    t_match = [t[i0]]
    while (i0,j0) in link_back:
        i1,j1=link_back[(i0,j0)]
        if i0==i1:
            t_match.append('-')
        else:
            t_match.append(t[i1])
        s_match.append(s[j1])
        i0,j0= i1,j1
    return (max(scores[-1]),''.join(s_match[::-1]),''.join(t_match[::-1]))

if __name__=='__main__':
    print(score_string('ACCATAAGCCCTACGTG-CCG',
                'GCCGTCAGGC-TG-GTGTCCG'))
    
    print (
        sims('GCAAACCATAAGCCCTACGTGCCGCCTGTTTAAACTCGCGAACTGAATCTTCTGCTTCACGGTGAAAGTACCACAATGGTATCACACCCCAAGGAAAC',
         'GCCGTCAGGCTGGTGTCCG'))
    