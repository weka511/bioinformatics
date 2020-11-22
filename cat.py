#    Copyright (C) 2019-2020 Greenweaves Software Limited
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
#    along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>
#
#    cat 	Catalan Numbers and RNA Secondary Structures (WIP)

import argparse
import os
import time
from helpers import read_strings


def count_matchings(seq):
    def partition(indices,i,j):
        I1 = []
        I2 = []
        for k in indices:
            if k==i: continue
            if k==j: continue
            if i<k and k <j:
                I1.append(k)
            else:
                I2.append(k)
        return (I1,I2)
    
    def count(indices):
        result = 0
        if 0 != sum(seq[i] for i in indices if abs(seq[i])==1): return 0
        if 0 != sum(seq[i] for i in indices if abs(seq[i])==2): return 0
        if len(indices)==0: return 1
        if len(indices)==2: return 1
        i = min(indices)
        for j in range(i+1,max(indices)+1,2):
            I1,I2  = partition(indices,i,j)
            count1 = count(I1)
            count2 = count(I2)
            result += (count1*count2)
        return result
    
    return count(list(range(len(seq))))
    
def cat(s):
    to_int = {'A':+1, 'U':-1, 'G':+2, 'C':-2}
    return count_matchings ([to_int[c] for c in s])
        
if __name__=='__main__':
    
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (cat('AUAU'))
        print (cat('UAGCGUGAUCAC'))
        print (cat('CGGCUGCUACGCGUAAGCCGGCUGCUACGCGUAAGC'))
  
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = cat(''.join(line for line in Input[1:]))
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')
                
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
