# Copyright (C) 2020 Greenweaves Software Limited

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

#  PDPL Creating a Restriction Map

import argparse
import math
import os
import random
import time
from   helpers import read_strings
        
def pdpl(L):
    def get_n(m):  # Solve simple quadratic Diophantine equation: n*(n-1) ==2*m
        return (1+math.isqrt(8*m+1))//2
    
    def compare(list1,list2):
        return sum([i1!=i2 for i1,i2 in zip(list1,list2)])
    
    def get_signature(L):
        signature = {}
        for i in L:
            if i in signature:
                signature[i]+=1
            else:
                signature[i]=1
        return signature  
    
    def get_diffs(X):
        return sorted([b-a for a in X for b in X if b>a])  
    
    def isCompatible(X):
        sx = get_signature(get_diffs(X))
        for x,count in sx.items():
            if x not in signature or count>signature[x]:
                return False
        return True
    
    m = len(L)
    n = get_n(len(L))
    assert n*(n-1)==2*m
    signature = get_signature(L)
    X = [0,L[-1]]
    for i in L:
        if isCompatible(X+[i]):             
            X.append(i)
    print (m,n,len(X))
    LL = get_diffs(X)
    print (compare(L,LL))
    return sorted(X)

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:  
        print (pdpl([2, 2, 3, 3, 4, 5, 6, 7, 8, 10]))
          
    if args.rosalind:
        Input     = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        converted = [int(i) for i in Input[0].split()]
        Result    = pdpl(sorted(converted))
        Formatted = ' '.join(str(r) for r in Result)
        print (Formatted)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
                f.write(f'{Formatted}\n')
                
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
