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

#  RNAS Wobble Bonding and RNA Secondary Structures

import argparse
import os
import time
from   helpers import read_strings
from   combinatorics import catmotz, partition

def count_wobble(seq):
    def wrapped_count(indices):
        def valid_bonds():
            j      = min(indices)
            for k in range(j+4,max(indices)+1):
                if seq[j] + seq[k]==0 or (seq[j]<0 and seq[k]<0): # Bond or wobble bond
                    yield (j,k)
                    
        def count():
            n      = len(indices)
            if n<2: return 1
            #j      = min(indices)          
            count1 = wrapped_count(indices[1:]) #  If first node is not involved in a matching            
            count2 = 0                            #  If first node is involved in a matching
            #for k in range(j+4,max(indices)+1):
                #if seq[j] + seq[k]==0 or (seq[j]<0 and seq[k]<0): # Bond or wobble bond
            for j,k in valid_bonds():
                I1,I2   = partition(indices,j,k)
                count21 = wrapped_count(I1)
                count22 = wrapped_count(I2)
                count2 += (count21*count22)                
            return count1 + count2
        key = str(indices)
        if not key in cache:
            cache[key] = count()        
        return cache[key]
    cache = {}
    return wrapped_count(list(range(len(seq))))         

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (catmotz('CGAUGCUAG',
                       counter=count_wobble))        
        print (catmotz('AUGCUAGUACGGAGCGAGUCUAGCGAGCGAUGUCGUGAGUACUAUAUAUGCGCAUAAGCCACGU',
                       counter=count_wobble))
 
        
    

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
