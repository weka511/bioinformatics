#!/usr/bin/env python

#    Copyright (C) 2019-2021 Greenweaves Software Limited
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
#    along with this program.  If not, see <http://www.gnu.org/licenses/>
#
#    BA6C Compute the 2-Break Distance Between a Pair of Genomes

import argparse
import os
import time
from   helpers import read_strings
from   fragile import d2break

def parse_permutations(text):
    Permutations = []
    depth        = 0
    i0           = None
    for i in range(len(text)):
        if text[i] =='(':
            depth += 1
            i0     = i
        else:
            if depth==0: continue
            if text[i]==')':
                depth -= 1
                i1     = i
                Permutations.append([int(j) for j in text[i0+1:i1-1].split()])

    return Permutations

if __name__=='__main__':

    start = time.time()
    parser = argparse.ArgumentParser('Find a Shortest Transformation of One Genome into Another by 2-Breaks')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (d2break([[+1, +2, +3, +4, +5, +6]],
                       [[+1, -3, -6, -5],[+2, -4]]))

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        Perms  = [parse_permutations(i) for i in Input]
        print (d2break(Perms[0],Perms[1]))

    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')

