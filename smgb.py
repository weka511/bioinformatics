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

#  SMGB Semiglobal Alignment

import argparse
import os
import time
from   helpers import read_strings
from   numpy import argmax
from   fasta import FastaContent

def smgb(s,t,match=+1,mismatch=-1,indel=-1):
    def backtrack(index,n):
        s1 = []
        t1 = []
        i  = index
        j  = n
        while i>0 and j>0:
            step = argmax([0,
                            scores[i-1][j]   + indel,
                            scores[i][j-1]   + indel,
                            scores[i-1][j-1] + (match if s[i-1]==t[j-1] else mismatch)])
            if step==0:
                print ('zero')
                break
            elif step==1:
                i-=1
                s1.append(s[i])
                t1.append('-')
            elif step==2:
                j-=1
                s1.append('-')
                t1.append(t[j])
            else:
                i-=1
                j-=1
                s1.append(s[i])
                t1.append(t[j])
        return s1,t1
    def reverse(chars):
        return ''.join(c for c in chars[::-1])
    m      = len(s)
    n      = len(t)
    if m<n: # assume t shorter than s
        return smgb(t,s,match=match,mismatch=mismatch,indel=indel)
    scores = [[0 for j in range(n+1)] for i in range(m+1)]
    for i in range(1,m+1):
        for j in range(1,n+1):
            scores[i][j] = max(0,
                               scores[i-1][j]   + indel,
                               scores[i][j-1]   + indel,
                               scores[i-1][j-1] + (match if s[i-1]==t[j-1] else mismatch))
        #print (''.join([f'{score:3d}' for score in scores[i]]) )
    last_column = [row[-1] for row in scores]
    index       = argmax(last_column)
    s1,t1       = backtrack(index,n)
    return last_column[index],reverse(s1),reverse(t1)
    

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('SMGB Semiglobal Alignment')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        score,s1,t1 = smgb('CAGCACTTGGATTCTCGG','CAGCGTGG')
        print (score)
        print (s1)
        print (t1)
        
        
    if args.rosalind:
        Input       = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        fasta       = FastaContent(Input)
        a,s         = fasta[0]
        b,t         = fasta[1]
        score,s1,t1 = smgb(s,t)

        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            print (score)
            print (s1)
            print (t1)                
            f.write(f'{score}\n')
            f.write(f'{s1}\n')
            f.write(f'{t1}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
