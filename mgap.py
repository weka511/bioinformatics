#   Copyright (C) 2020 Greenweaves Software Limited

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

#  MGAP Maximizing the Gap Symbols of an Optimal Alignment

import argparse
import os
import time
from   laff import read_fasta

# mgap
#
# For the computation of an alignment score generalizing the edit alignment score, 
# let m denote the score assigned to matched symbols, d denote the score assigned to mismatched non-gap symbols,
# and g denote the score assigned a symbol matched to a gap symbol '-' (i.e., gis a linear gap penalty).

# Given: Two DNA strings s and t,  each of length at most 5000 bp.
#
# Return: The maximum number of gap symbols that can appear in any maximum score alignment of s
# and t with score parameters satisfying m>0, d<0, and g<0.

def mgap(s,t,m0=1,d0=-1,g0=-1,Nm=2,Nd=2,Ng=2):
    def get_N(m=1,d=-1,g=-1):
        def dynamic_programming():
            scores = [[0 for j in range(len(t)+1)] for i in range(len(s)+1)]
    
            for j in range(len(t)+1):
                scores[0][j] = g * j
            for i in range(len(s)+1):
                scores[i][0] = g * i
    
            for i in range(1,len(s)+1):
                for j in range(1,len(t)+1):
                    scores[i][j] = max(
                        scores[i-1][j]   + g,
                        scores[i][j-1]   + g,
                        scores[i-1][j-1] + (m if s[i-1]==t[j-1] else d))
    
            return scores
        
        def backtrack(scores):
            gaps = 0
            i    = len(scores) - 1
            j    = len(scores[0]) -1
            while i>0 or j>0:
                if scores[i][j]==scores[i-1][j]   + g:
                    #print (s[i-1], '-')
                    i    -= 1
                    gaps += 1
                elif  scores[i][j]==scores[i][j-1]   + g:
                    #print ('-', t[j-1])
                    j    -= 1
                    gaps += 1
                elif  scores[i][j]==scores[i-1][j-1] + (m if s[i-1]==t[j-1] else d):
                    #print (s[i-1], t[j-1])
                    i  -= 1
                    j  -= 1
                else:
                    raise Exception(f'{i} {j}')
            return gaps    
        
        return backtrack(dynamic_programming())
    # Nm=2,Nd=2,Ng=2
    max_gaps = -1
    for i in range(Nm):
        for j in range(Nd):
            for k in range(Ng):
                gaps = get_N(m=m0+i, d=d0-j, g = g0-k)
                if gaps>max_gaps:
                    max_gaps=gaps
    return max_gaps

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('MGAP Maximizing the Gap Symbols of an Optimal Alignment')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='proces[1]s Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (mgap('AACGTA','ACACCTA'))
        
    if args.rosalind:
        Data      = read_fasta(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = mgap(Data[0],Data[1])
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
