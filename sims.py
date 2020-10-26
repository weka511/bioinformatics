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

def get_score(a,b, match=1,mismatch=-1):
    return match if a==b else mismatch

def score_string(s,t):
    return sum(get_score(a,b) for (a,b) in zip(s,t))

def sims(s,t, match=1,mismatch=-1):
    
    def dynamic_programming():
        scores      = np.zeros((len(t),len(s)),dtype=int)
        predecessor = {}
        for j in range(len(s)):
            scores[0]  = get_score(s[j],t[0],match=match,mismatch=mismatch)
        for i in range(1,len(t)):
            for j in range(i-1):
                scores[i,j]  = get_score(s[j],t[i],match=match,mismatch=mismatch)
            for j in range(i,len(s)):
                score_diag       = scores[i-1][j-1] + get_score(s[j],t[i],match=match,mismatch=mismatch)
                score_horizontal = scores[i][j-1]   + mismatch
                score_vertical   = scores[i-1][j]   + mismatch
                scores[i,j]      = max(score_diag,
                                       score_horizontal,
                                       score_vertical)
    
                predecessor[(i,j)] = (i-1, j-1) if scores[i,j] == score_diag       else \
                                     (i,   j-1) if scores[i,j] == score_horizontal else \
                                     (i-1, j)   if scores[i,j] == score_vertical   else \
                                     None
        return scores,predecessor

    def trace_back(scores,predecessor):
        i      = len(scores)-1
        j      = np.argmax(scores[i])
 
        print (f'Number of maxima={sum(1 for s in scores[i] if s>=scores[i,j])}')
        s_match = [s[j]]
        t_match = [t[i]]
        while (i,j) in predecessor:
            i_pre,j_pre = predecessor[(i,j)]
            if i==i_pre:
                assert j != j_pre
                t_match.append('-')
                s_match.append(s[j_pre])
            elif j==j_pre:
                t_match.append(t[i_pre])
                s_match.append('-')            
            else:
                t_match.append(t[i_pre])
                s_match.append(s[j_pre])
            print (f'{(i_pre,j_pre)}=>{(i,j)} s={s_match[-1]} t={t_match[-1]}')    
            i,j = i_pre,j_pre
            
        return (s_match,t_match)
    
    scores,predecessor = dynamic_programming()
    
    print (max(scores[-1]))
    print (scores[-1])

    s_match,t_match = trace_back(scores,predecessor)

    return (max(scores[-1]), ''.join(s_match[::-1]), ''.join(t_match[::-1]))

if __name__=='__main__':
     
    score,s1,t1 = sims('GCAAACCATAAGCCCTACGTGCCGCCTGTTTAAACTCGCGAACTGAATCTTCTGCTTCACGGTGAAAGTACCACAATGGTATCACACCCCAAGGAAAC',
                       'GCCGTCAGGCTGGTGTCCG')
    print (score)
    print (len(s1),s1)
    print (len(t1),len('GCCGTCAGGCTGGTGTCCG'),t1)
    print (score_string(s1,t1))
       