#!/usr/bin/env python

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

# 2SAT 2-Satisfiability

import argparse
import os
import time
from   helpers import read_strings,create_list,parse_graphs
from   graphs import two_sat

# to_int
#
# Used to parse a line of input
#
# Convert a line of text representing a single number to an int
# and a line of text representing a list to a list of ints

def to_int(s):
    parts=s.split()
    if len(parts)==1:
        return int(s)
    else:
        return [int(p) for p in parts]

# create_sets
#
# Create sets of data for 2SAT problem. Each set is in the form given in problem statement, e.g.
# n,m              i.e. number of variables, number of clauses
# clause 1
# clause 2
# ...
# clause m
#
# Output: list of form [(n,m, (...), (...) ...],,,]  where each inner tuple is a clause, each out tuple a set

def create_sets(data):
    product = []
    i       = 1
    while i<len(data):
        n,m = data[i]
        product.append((n,m, [(xx[0],xx[1]) for xx in data[i+1:i+1+m]]))
        i = i+1+m
    assert len(product)==data[0]
    return product

# Format
#
# Format Solution for display as a single line of numbers

def Format(status,Solution):
    return f'{status} {" ".join(str(sol) for sol in Solution)}'

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('2SAT')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:

        data = [2,
                [2, 4],
                [1 ,2],
                [-1, 2],
                [1, -2],
                [-1, -2],

                [3, 4],
                [1, 2],
                [2, 3],
                [-1, -2],
                [-2, -3]]

        for problem in create_sets(data):
            status,Solution = two_sat(problem )
            print (Format(status,Solution))



    if args.rosalind:
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
            Sets   = create_sets([to_int(s) for s in Input if len(s)>0])
            for problem in Sets:
                status,Solution = two_sat(problem )
                print (Format(status,Solution))
                f.write (f'{Format(status,Solution)}\n')

    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
