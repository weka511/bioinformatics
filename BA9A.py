# Copyright (C) 2019-2020 Greenweaves Software Limited

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

# BA9A Construct a Trie from a Collection of Patterns 

import argparse
import time
from helpers import read_strings
from snp import Trie
                  
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9D Find the Longest Repeat in a String')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        trie = Trie(['ATAGA','ATC','GAT'])
        for a,b,c in trie.Edges():
            print (f'{a}->{b}:{c}') 
    
    if args.extra:
        Input,Expected  = read_strings('data/TrieConstruction.txt',init=0)
        trie = Trie(Input)
        Actual = sorted([f'{a}->{b}:{c}' for a,b,c in trie.Edges()])
        Expected.sort()
        print (len(Expected),len(Actual))
        diffs = [(e,a) for e,a in zip(Expected,Actual) if e!=a]
        print (diffs)
  
    if args.rosalind:
        Input  = read_strings('data/rosalind_ba9a.txt')
        trie   = Trie(Input)
        Result = sorted([f'{a}->{b}:{c}' for a,b,c in trie.Edges()])
        print (Result)
        with open('ba9a.txt','w') as f:
            for line in Result:
                f.write(f'{line}\n')
                
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
        
