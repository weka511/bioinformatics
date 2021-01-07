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
from   numpy   import zeros, amax 
# itwv
#
# Given: A text DNA string s of length at most 10 kbp, followed by a collection of n DNA strings of
#        length at most 10 bp acting as patterns.
#
# Return: An nxn matrix M for which M[j,k]==1 if the jth and kth pattern strings can be interwoven into s and Mj,k=0 otherwise.

def itwv2(s,u,v):
    def score(a,b):
        return 1 if a==b else 0
    def match(k):
        def update_score(i,j):
            if (i,j) in Closed: return
            #print (i,j, i+j)
            s0 = s[k+i+j-1]
            u0 = u[i-1]
            v0 = v[j-1]
            S[i,j] = max(S[i-1,j] + score(s0,u0),S[i,j-1] + score(s0,v0))
            Closed.add((i,j))
        S      = zeros((len(u)+1,len(v)+1),dtype=int)
        N      = max(len(u),len(v))+1
        Closed = set()
        Closed.add((0,0))
        for n in range(1,N+1):
            #print (f'n={n}')
            for i in range(0,min(n,len(u)+1)):
                for j in range(0,min(n,len(v)+1)):
                    update_score(i,j)

        return amax(S)==len(u)+len(v)
    
    for i in range(len(s)-len(u)-len(v)+1):
        if s[i]==u[0] or s[i]==v[0]:
            if match(i)>0:
                return 1
    return 0

def itwv(s,patterns):        
    return [[itwv2(s,u,v)  for v in patterns] for u in patterns]



if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('ITWV Finding Disjoint Motifs in a Gene ')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        #print (itwv2('GACCACGGTT','GT','GT'))
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
