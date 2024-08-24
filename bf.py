#!/usr/bin/env python

# Copyright (C) 2020-2024 Greenweaves Software Limited

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

'''BF Bellman-Ford Algorithm'''

from argparse import ArgumentParser
from os.path  import basename
from time     import time

from graphs import bf
from helpers import read_strings,create_list

def Format(dists ):
    return ' '.join(str(d) if d<float('inf') else 'x' for d in dists )

if __name__=='__main__':
    start = time()
    parser = ArgumentParser(__doc__)
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        edges = [(9, 13),
                 (1, 2, 10),
                 (3, 2, 1),
                 (3, 4, 1),
                 (4, 5, 3),
                 (5, 6, -1),
                 (7, 6, -1),
                 (8, 7, 1),
                 (1, 8, 8),
                 (7, 2, -4),
                 (2, 6, 2),
                 (6, 3, -2),
                 (9, 5 ,-10),
                 (9, 4, 7)]
        _,dists,_ = bf(edges)
        print (Format(dists))


    if args.rosalind:
        _,dists,_ = bf(create_list(path='./data'))
        print (Format(dists))

        with open(f'{basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Format(dists)}\n')

    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
