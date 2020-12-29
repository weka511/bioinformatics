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

#  MREP Identifying Maximal Repeats

import argparse
import os
import time
from   helpers import read_strings
from   snp     import SuffixArray
from   numpy   import argsort

def mrep(w,ml=20):
    n       = len(w)
    r,p,LCP = SuffixArray(w,auxiliary=True) # 0..n-1, 0..n-1, 0..n-2
    S       = [u for u in range(len(LCP)) if LCP[u]<ml]
    S.append(-1)
    S.append(n-1)
    I       = argsort(LCP)
    initial = min([t for t in range(len(I)) if LCP[I[t]]>=ml])
 
    for t in range(initial,n-1):
        i   = I[t]
        p_i = max([j for j in S if j<i])+1
        n_i = min([j for j in S if j>i])
        S.append(i)
        if (p_i==0 or LCP[p_i-1]!=LCP[i]) and (n_i==n-1 or LCP[n_i]!=LCP[i]):
            if r[p_i]==0 or r[n_i]==0 or w[r[p_i]-1]!=w[r[n_i]-1] or p[r[n_i]-1]-p[r[p_i]-1]!=n_i-p_i:
                yield w[r[i]:r[i]+LCP[i]]


if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('MREP Identifying Maximal Repeats')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    
    if args.extra:
        for r in mrep('TAGTTAGCGAGA',ml=2):
            print (r)    
            
    if args.sample:
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in mrep('TAGAGATAGAATGGGTCCAGAGTTTTGTAATTTCCATGGGTCCAGAGTTTTGTAATTTATTATATAGAGATAGAATGGGTCCAGAGTTTTGTAATTTCCATGGGTCCAGAGTTTTGTAATTTAT'):
                print (line)
                f.write(f'{line}\n')        
         
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in mrep(Input[0]):
                print (line)
                f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
