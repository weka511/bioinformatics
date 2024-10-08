#!/usr/bin/env python

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
import os
import time
from   helpers import read_strings
from   snp import create_trie

def convert_trie(Trie):
    return [(key,node,symbol)for key,Edge in Trie.items() for symbol,node in Edge.items()]

def format(key,node,symbol):
    return f'{key}->{node}:{symbol}'

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9A Construct a Trie from a Collection of Patterns')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        Trie = create_trie(['ATAGA','ATC','GAT'],root=0)
        for key,node,symbol in convert_trie(Trie):
            print (format(key,node,symbol))


    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        Trie = create_trie(Input,root=0)

        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for key,node,symbol in convert_trie(Trie):
                print (format(key,node,symbol))
                f.write(f'{format(key,node,symbol)}\n')

    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
