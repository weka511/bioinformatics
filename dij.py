#!/usr/bin/env python

#    Copyright (C) 2019-2024 Greenweaves Software Limited
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
#    along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

'''
DIJ  Dijkstra's Algorithm: compute single-source shortest distances
                            in a directed graph with positive edge weights.
'''

from argparse import ArgumentParser
from os.path import basename
from time import time
from helpers import read_strings
from graphs import dij
from helpers import create_list

if __name__=='__main__':
    start = time()
    parser = ArgumentParser(__doc__)
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (dij([[6 ,10],
                   [1, 2, 4],
                   [1, 3, 2],
                   [2, 3, 3],
                   [6, 3, 2],
                   [3, 5, 5],
                   [5, 4, 1],
                   [3, 2, 1],
                   [2, 4, 2],
                   [2, 5, 3]]))

    if args.rosalind:
        with open(f'{basename(__file__).split(".")[0]}.txt','w') as f:
            Solution = ' '.join([str(int(i)) for i in dij(create_list(path='./data'))])
            print (Solution)
            f.writelines(f'{Solution}\n')

    elapsed = time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
