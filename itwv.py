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

#  ITWV Finding Disjoint Motifs in a Gene 

import argparse
import os
import time
from   helpers import read_strings,format_list
from   graphs   import trie

def itwv(s,patterns):
    def can_interweave(s,t):
        J = len(s)
        K = len(t)
        for i in range(n-J-K):
            j = 0
            k = 0
        return 0
    def interweave(j,k):
        if can_interweave(patterns[j],patterns[k]):
            return 1
        else:
            return can_interweave(patterns[k],patterns[j])
    n = len(patterns)
    M = [[float('nan') for k in range(n)] for j in range(n)]
    for j in range(n):
        for k in range(j,n):
            M[j][k] = interweave(j,k)
            if k>j:
                M[k][j] = M[j][k]
    return M

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        for line in itwv('GACCACGGTT',
                    ['ACAG',
                     'GT',
                     'CCG']):
            print (format_list(line))
        
    

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = itwv(Input[0],Input[1:])
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                out = format_list(line)
                f.write(f'{out}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
