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


# BA9B Implement TrieMatching

import argparse
import os
import time
from helpers import read_strings
from snp import Trie,create_trie


def MatchPrefix(Text,Trie):
    i      = iter(Text)
    symbol = next(i)
    v      = min(Trie.keys())
    path   = [v]
    while True:
        if len(Trie[v])==0:
            return path
        elif symbol in Trie[v]:
            w      = Trie[v][symbol]
            try:
                symbol = next(i)
                v      = w
                path.append(v)
            except StopIteration:
                if len(Trie[w])==0:
                    return path
                else:
                    return []
        else:
            return []

def MatchAll(Text,Trie):
    return [i for i in range(len(Text)) if MatchPrefix(Text[i:],Trie)]


if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9B Implement TrieMatching')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:    
        Trie = create_trie(['ATCG','GGGT'])
        print (MatchAll('AATCGGGTTCAATCGGGGT',Trie))
    
    if args.extra:
        Input,Expected  = read_strings('data/TrieMatching.txt',init=0)
        Trie            = create_trie(Input[1:])
        Actual          = MatchAll(Input[0],Trie)
        Expected        = [int(e) for e in Expected[0].split()]
        print (len(Expected),len(Actual))
        diffs = [(e,a) for e,a in zip(Expected,Actual) if e!=a]
        print (diffs)
        
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        Trie   = create_trie(Input[1:])
        Result = MatchAll(Input[0],Trie)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            Output = ' '.join(str(n) for n in Result)
            print (Output)
            f.write(f'{Output}\n')
                    
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s') 
    