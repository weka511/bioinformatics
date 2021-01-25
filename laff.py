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
#
# Rewritten to use pypy - cannot use numpy, etc.

import argparse
import os
import time
import sys

def create_blosum62() :
 
    amino_acids = 'A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y'.split()
           
    scores = [
              4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2,
              0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2,
              -2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3,
              -1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2,
              -2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3,
              0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3,
              -2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2,
              -1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1,
              -1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2,
              -1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1,
              -1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1,
              -2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2,
              -1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3,
              -1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1,
              -1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2,
              1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2,
              0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2,
              0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1,
              -3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2,
              -2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7    ]
    
    i = 0
    Product = {}
    for x in amino_acids:
        for y in amino_acids:
            Product[(x,y)] = scores[i]
            i += 1
    return Product

# reverse
#
# Input:    a string, 
# Return:   new string with characters in reverse order

def reverse(chars):
    return ''.join(chars[i] for i in range(len(chars)-1,-1,-1))

# laff
#
# Local Alignment with Affine Gap Penalty

def laff(s,t,
         replace_score = create_blosum62(),
         sigma         = 11,
         epsilon       = 1,
         frequency     = 0):

    # get_max_score
    
    def get_max_score():
        imax      = -1
        max_score = -1
        jmax      = -1        
        for i in range(1,m+1):
            for j in range(1,n+1):
                if middle[i][j] > max_score:
                    imax      = i
                    jmax      = j
                    max_score = middle[i][j]    
        return max_score,imax,jmax

    def unwind_moves(score,i,j):
        print (f'Unwinding from: {score}, {i}, {j}')
        ss = []
        tt = []
        assert score == middle[i][j],f'score={score}, middle[{i}][{j}]={middle[i][j]}'
        while i>0 and j > 0:
            if middle[i][j] == middle[i-1][j-1] + replace_score[(s[i-1],t[j-1])]:
                i -= 1
                j -= 1
                ss.append(s[i])
                tt.append(t[j])
            elif middle[i][j]==upper[i][j]: # See Yury Grushetsky's hint http://rosalind.info/problems/laff/questions/
                j -= 1
                tt.append(t[j])                
            elif middle[i][j]==lower[i][j]: # See Yury Grushetsky's hint http://rosalind.info/problems/laff/questions/
                i -= 1
                ss.append(s[i])
            else:
                raise Exception(f'This cannot possible happen {i} {j}!')
            
            if middle[i][j]==0:
                print (f'Exiting at {i} {j}')
                break
            
        return score,reverse(ss),reverse(tt)
    
    m      = len(s)
    n      = len(t)
 
    lower  = [[0 for j in range(n+1)] for i in range((m+1)) ]
    middle = [[0 for j in range(n+1)] for i in range((m+1)) ]
    upper  = [[0 for j in range(n+1)] for i in range((m+1)) ]

    start     = time.time()
    for i in range(1,m+1):
        
        if frequency>0 and i%frequency==0:
            T = time.time()-start
            print (f'Row={i}/{m}, elapsed={T:.1f} sec., time/step={T/i:.4f} sec., ETA={m*T/i-T:.0f} sec.')
            
        for j in range(1,n+1):
            lower[i][j]      = max(lower[i-1][j]  - epsilon,
                                   middle[i-1][j] - sigma,
                                   0)
            upper[i][j]      = max(upper[i][j-1]  - epsilon,
                                   middle[i][j-1] - sigma,
                                   0)
            middle[i][j]     = max(lower[i][j], 
                                   upper[i][j], 
                                   middle[i-1][j-1] + replace_score[(s[i-1],t[j-1])],
                                   0)                             
    
    max_score,imax,jmax = get_max_score()
   
    return unwind_moves(max_score,imax,jmax)

def read_fasta(name):
    Data   = []
    Record = []
    with open(name,'r') as input:
        for line in input:
            text = line.strip()
            if text[0]=='>':
                if len(Record)>0:
                    Data.append(''.join(Record))
                Record = []
            else:
                Record.append(text)
    Data.append(''.join(Record))
    return Data

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('LAFF Local Alignment with Affine Gap Penalty')
    parser.add_argument('--sample',    default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind',  default=False, action='store_true', help='process Rosalind dataset')
    parser.add_argument('--version',   default=False, action='store_true', help='Get version of python')
    parser.add_argument('--frequency', default=100,   type=int,            help='Number of iteration per progress tick' )
    args = parser.parse_args()
    
    if args.version:
        print (f'{sys.version}')
        
    if args.sample:
        score,s1,t1 = laff('PLEASANTLY','MEANLY')
        print (score)
        print (s1)
        print (t1)    
        
    if args.rosalind:
        Data      = read_fasta(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        score,s,t = laff(Data[0],Data[1],frequency=args.frequency)
        print (score)
        print (s)
        print (t)   
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:               
            f.write(f'{score}\n')
            f.write(f'{s}\n')
            f.write(f'{t}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
