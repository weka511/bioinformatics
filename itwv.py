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
import itertools

# itwv
#
# Given: A text DNA string s of length at most 10 kbp, followed by a collection of n DNA strings of
#        length at most 10 bp acting as patterns.
#
# Return: An nxn matrix M for which M[j,k]==1 if the jth and kth pattern strings can be interwoven into s and Mj,k=0 otherwise.

def itwv(s,patterns):

    def can_match(p,q,l):
        ss = s_list[l:l+len(p)+len(q)]
        indices = set()
        matched = False
        for index in itertools.permutations([0]*len(p) + [1]*len(q)):
            if index in indices: continue
            indices.add(index)
            matched = True # assumption
            i  = 0
            j  = 0
            for k in range(len(p)+len(q)):
                if index[k]==0: # sample p
                    if ss[k]==p[i]:
                        i += 1
                    else:
                        matched = False
                        break    # from k
                else:           #sample q
                    if ss[k]==q[j]:
                        j+= 1
                    else:
                        matched = False
                        break  # from k
            if matched: return 1
        return 0
    
    def can_interweave(p,q):
        L = m - len(p) -len(q)
        for l in range(L+1):
            if can_match(p,q,l): return 1
        return 0
    
    s_list = list(s)
    p_list = [list(p) for p in patterns]
    n      = len(p_list)
    m      = len(s_list)
    return [[can_interweave(p_list[i],p_list[j]) for j in range(n)] for i in range(n)]

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('ITWV Finding Disjoint Motifs in a Gene ')
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
