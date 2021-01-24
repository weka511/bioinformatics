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

def create_blosum62() :
    Product = {}
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
        -2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7
    ]
    i = 0
    for x in amino_acids:
        for y in amino_acids:
            Product[(x,y)] = scores[i]
            i += 1
    return Product
    
#from   helpers   import read_strings
#from   numpy     import argmax,zeros
#from   fasta     import FastaContent
#from   Bio.Align import substitution_matrices
#from   align     import reverse
#from   Bio       import SeqIO
# laff
#
# Local Alignment with Affine Gap Penalty

def laff(s,t,replace_score=create_blosum62(),sigma=11,epsilon=1):
 
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
    print (m,n)
    lower  = [[0 for j in range(n+1)] for i in range(2) ]
    middle = [[0 for j in range(n+1)] for i in range(2) ]
    upper  = [[0 for j in range(n+1)] for i in range(2) ]
    #moves  = {}
    
    imax      = -1
    max_score = -1
    jmax      = -1
    start = time.time()
    for i in range(1,m+1):
        
        i_to      = i%2
        i_from    = 1 - i_to
        if i%100==0:
            print (f'{i} {m} {int(time.time()-start)}')
        for j in range(1,n+1):
            lower[i_to][j]      = max(lower[i_from][j]  - epsilon,
                                  middle[i_from][j] - sigma)
            upper[i_to][j]      = max(upper[i_to][j-1]  - epsilon,
                                  middle[i_to][j-1] - sigma)
            possible_scores = [lower[i_to][j], 
                               upper[i_to][j], 
                               middle[i_from][j-1] + replace_score[(s[i-1],t[j-1])]]
            #possible_moves = [(i-1, j,   s[i-1], '-'   ),  # Comes from lower
                              #(i-1, j-1, s[i-1], t[j-1]),  # Comes from middle
                              #(i,   j-1, '-',    t[j-1])]  # Comes from upper
            
            max_score      = max(possible_scores)
            for index in range(len(possible_scores)):
                if possible_scores[index]==max_score:
                    break
            middle[i_to][j] = possible_scores[index]                             
            #moves[(i,j)]   = possible_moves[index]
        
        score =  max(middle[i_to])
        for j in range(1,n+1):
            if middle[i_to][j]  == score: break
            
        #score = middle[i_to,j]
        if score>max_score:
            imax      = i
            jmax      = j
            max_score = score
        
    return unwind_moves(None,max_score,imax,jmax)

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
        Data = []
        Record = []
        with open ('data/rosalind_laff.txt') as input:
            i = -1
            for line in input:
                text = line.strip()
                #print(text)
                if text[0]=='>':
                    i += 1
                    if len(Record)>0:
                        Data.append(''.join(Record))
                    Record = []
                else:
                    Record.append(text)
        Data.append(''.join(Record))
        x=0
        #data = [str(rec.seq) for rec in SeqIO.parse(open(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt','r'),
                                                         #'fasta')]                
        score,s,t = laff(Data[0],Data[1])

        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:               
            f.write(f'{score}\n')
            f.write(f'{s}\n')
            f.write(f'{t}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
