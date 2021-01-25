#    Copyright (C) 2019-2021 Greenweaves Software Limited
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
#    along with this program.  If not, see <http://www.gnu.org/licenses/>

#    OAP Overlap Alignment

import argparse
import os
import time
import sys
from   laff import read_fasta

# oap
#
# Given: Two DNA strings s and t, each having length at most 10 kbp.
#
# Return: The score of an optimal overlap alignment of s and t, followed by an alignment of a suffix  of s
#         and a prefix  of t achieving this optimal score. 

def oap(s,t,
        match_bonus   = 1,
        mismatch_cost = 2,
        indel_cost    = 2):
    def score(v,w):
        return match_bonus if v==w else -mismatch_cost
        
    def dynamic_programming(s,t):
        
        distances = [[0 for j in range(len(t)+1)] for i in range(len(s)+1)]
        for i in range(1,len(s)+1):
            for j in range(1,len(t)+1):
                distances[i][j]  = max(
                                        distances[i-1][j]   - indel_cost,
                                        distances[i][j-1]   - indel_cost,
                                        distances[i-1][j-1] + score(s[i-1],t[j-1]))
        
        # Begin backtracking. Since we want a suffix of s, start in the last row.
        # For some reason I need to use the last j that matches maximum score
        
        i = len(s)
        distance = max(distances[i])
        for j in range(len(t),0,-1):
            if distances[i][j] == distance: break
 
        s1       = []
        t1       = []
        while j>0:                  # we want a prefix of t
            if distances[i][j]==distances[i-1][j]   - indel_cost:
                i1,j1 = (i-1,j)
            elif distances[i][j]==distances[i][j-1]   - indel_cost:
                i1,j1 = (i,j-1)
            elif distances[i][j]==distances[i-1][j-1] + score(s[i-1],t[j-1]):
                i1,j1 = (i-1,j-1)
            else:
                raise Exception(f'This cannot possible happen {i} {j}!')
             
            s1.append(s[i1] if i1<i else '-')
            t1.append(t[j1] if j1<j else '-')
   
            i,j = i1,j1
            
        return distance,s1[::-1],t1[::-1]
    
    score,u1,s1=dynamic_programming(s,t)
    return score,''.join(u1),''.join(s1)

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('OAP Overlap Alignment')
    parser.add_argument('--sample',    default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind',  default=False, action='store_true', help='process Rosalind dataset')
    parser.add_argument('--version',   default=False, action='store_true', help='Get version of python')
    args = parser.parse_args()
    
    if args.version:
        print (f'{sys.version}')
        
    if args.sample:
        d,s1,t1 = oap('CTAAGGGATTCCGGTAATTAGACAG',
                      'ATAGACCATATGTCAGTGACTGTGTAA')
        
        print ('{0}'.format(d))
        print (s1)
        print (t1)
        
    if args.rosalind:
        Data      = read_fasta(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')          
        d,s1,t1   = oap(Data[0],Data[1])      
        print ('{0}'.format(d))
        print (s1)
        print (t1)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as o:
            o.write('{0}\n'.format(d))
            o.write('{0}\n'.format(s1))
            o.write('{0}\n'.format(t1))    

    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')        