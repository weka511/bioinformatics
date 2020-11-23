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

#  Combinatorics The mathematics of counting objects.

import argparse
import os
import time
from   helpers import read_strings

# catalan
#
# Calculate Catalan numbers
#
# Verified against list from http://mathforum.org/advanced/robertd/catalan.html
#
# 1 1
# 2 2
# 3 5
# 4 14
# 5 42
# 6 132
# 7 429
# 8 1430
# 9 4862
# 10 16796
# 11 58786
# 12 208012
# 13 742900
# 14 2674440
# 15 9694845
#
# Verified against Online Encyclopedia of Integer Sequences (https://oeis.org/A000108) for case:
# 30 3814986502092304

class Catalan:
    def __init__(self):
        self.c = [1]
        
    def get(self,n):
        for m in range(len(self.c),n+1):
            self.c.append(sum([self.c[k]*self.c[m-1-k] for k in range(m)]))
        return self.c[n]

class Motzkin:
    def __init__(self):
        selfs=[1,1]
        
    def get(n):
        for n0 in range(1,n):
            self.s.append(self.s[-1]+sum([self.s[k]*self.s[n0-1-k] for k in range(n0)]))
        return self.s[n]

#    cat 	Catalan Numbers and RNA Secondary Structures (WIP)

# partition
#
# Split set into two partitions, one between i and and j, one outside

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

def count_perfect_matchings(seq):

    def count(indices):
        key = str(indices)
        if key in cache:
            return cache[key]
        
        result = 0
        if 0 != sum(seq[i] for i in indices if abs(seq[i])==1): return 0
        if 0 != sum(seq[i] for i in indices if abs(seq[i])==2): return 0
        if len(indices)==0: return 1
        if len(indices)==2: return 1
        i = min(indices)
        for j in range(i+1,max(indices)+1,2):
            if seq[i] + seq[j]!=0: continue # Is i-j a valid split of the data?
            I1,I2  = partition(indices,i,j)
            result += count(I1)*count(I2)
            
        cache[str(indices)]= result
        return result
    
    cache = {}
    return count(list(range(len(seq))))
    
def catmotz(s,counter=count_perfect_matchings):
    to_int = {'A':+1, 'U':-1, 'G':+2, 'C':-2}
    return counter([to_int[c] for c in s])
    
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        catalan = Catalan()
        for i in range(31):
            print (i,catalan.get(i))
 
        
    

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = None
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
