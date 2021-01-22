#   Copyright (C) 2020-2021 Greenweaves Software Limited

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

#  LAFF Local Alignment with Affine Gap Penalty

import argparse
import os
import time
from   helpers   import read_strings
from   numpy     import argmax,zeros
from   fasta     import FastaContent
from   Bio.Align import substitution_matrices
from   align     import reverse
from   Bio       import SeqIO
# laff
#
# Local Alignment with Affine Gap Penalty

def laff(s,t,replace_score=substitution_matrices.load("BLOSUM62"),sigma=11,epsilon=1):
    def score_pair(a,b):
        return int(replace_score[(a,b)] if (a,b) in replace_score else replace_score[(b,a)])
    def unwind_moves(moves,score,i,j):
        ss = []
        ts = []
    
        while i>0 and j > 0:
            i,j,s0,t0=moves[(i,j)]
            ss.append(s0)
            ts.append(t0)
        return score,reverse(ss),reverse(ts)
    
    m      = len(s)
    n      = len(t)
    lower  = zeros((m+1,n+1),dtype=int)
    middle = zeros((m+1,n+1),dtype=int)
    upper  = zeros((m+1,n+1),dtype=int)
    moves  = {}
    
    imax      = -1
    max_score = -1
    jmax      = -1
    
    for i in range(1,m+1):
        if i%10==0: print (f'{i} {m}')
        for j in range(1,n+1):
            lower[i,j]      = max(lower[i-1,j]  - epsilon,
                                  middle[i-1,j] - sigma)
            upper[i,j]      = max(upper[i,j-1]  - epsilon,
                                  middle[i,j-1] - sigma)
            possible_scores = [lower[i,j], 
                               upper[i,j], 
                               middle[i-1,j-1] + score_pair(s[i-1],t[j-1])]
            possible_moves = [(i-1, j,   s[i-1], '-'   ),  # Comes from lower
                              (i-1, j-1, s[i-1], t[j-1]),  # Comes from middle
                              (i,   j-1, '-',    t[j-1])]  # Comes from upper
            
            index          = argmax(possible_scores)
            middle[i,j]    = possible_scores[index]                             
            moves[(i,j)]   = possible_moves[index]
            
        j     = argmax(middle[i])
        score = middle[i,j]
        if score>max_score:
            imax      = i
            jmax      = j
            max_score = score
        
    return unwind_moves(moves,max_score,imax,jmax)

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('LAFF Local Alignment with Affine Gap Penalty')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        score,s1,t1 = laff('PLEASANTLY','MEANLY')
        print (score)
        print (s1)
        print (t1)    
        
    if args.rosalind:
        data = [str(rec.seq) for rec in SeqIO.parse(open(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt','r'),
                                                         'fasta')]                
        score,s,t = laff(data[0],data[1])

        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:               
            f.write(f'{score}\n')
            f.write(f'{s}\n')
            f.write(f'{t}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
