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

#  BA9G Construct the Suffix Array of a String

import argparse
import os
import time
from   helpers import read_strings

# Suffix Array
#
# In 1993, Udi Manber and Gene Myers introduced suffix arrays as a memory-efficient alternative to suffix trees. 
# To construct SuffixArray(Text), we first sort all suffixes of Text lexicographically, assuming that "$" 
# $ comes first in the alphabet. The suffix array is the list of starting positions of these sorted suffixes.
#
# I have used a naive algorithm, as it appears to be adequate for the test data

def SuffixArray(s):
    suffixes = [(s[i:],i) for i in range(len(s))]
    suffixes.sort()
    return [i for (_,i) in suffixes]

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9G Construct the Suffix Array of a String')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra', default=False, action='store_true', help='process Extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (SuffixArray('panamabananas$'))
        print (SuffixArray('AACGATAGCGGTAGA$'))
    
    if args.extra:
        Input,Expected  = read_strings('data/SuffixArray.txt',init=0)
 
        Actual=SuffixArray(Input[0])
 
        E = [int(s) for s in Expected[0].split(',')]
        count = 0
        if len(E)==len(Actual):
            for a,b in zip(E,Actual):
                if a!=b:
                    print (f'Mismatch {a} {b}')
                    count+=1
        else:
            print("Lengths don't match")
        if count==0:
            print ('Matched')
        
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result =SuffixArray(Input[0])
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            line = ', '.join(str(i) for i in Result)
            f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
