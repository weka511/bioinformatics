#!/usr/bin/env python
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

#  MMCH Maximum Matchings and RNA Secondary Structures

import argparse
import os
import time
from   helpers import read_strings
from math import factorial
import decimal
def mmch(s):
    def comb(n,r):
        return factorial(n) // factorial(n-r)
    def prod(a,b):
        return int(comb(max(counts[a],counts[b]), min(counts[a],counts[b])))

    counts = {'A':0,'C':0,'G':0,'U':0}
    for c in s:
        counts[c]+=1
    return prod('A','U') * prod('C','G')

if __name__=='__main__':
    decimal.getcontext().prec=100
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--sample2',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--sample3',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (mmch('AUGCUUC'))
    if args.sample2:
        print (mmch('CAGCGUGAUCAC'))
    if args.sample3:
        print (mmch('CAGCGUGAUCACCAGCGUGAUCAC'))


    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')

        Result = mmch(''.join(s for s in Input[1:]))
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')

    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
