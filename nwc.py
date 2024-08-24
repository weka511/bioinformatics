#!/usr/bin/env python

# Copyright (C) 2019-2024 Greenweaves Software Limited

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

'''NWC Negative Weight Cycle'''

from argparse import ArgumentParser
from os.path import basename
from time import time
from helpers import read_strings,create_list
from graphs import bf

def extract_data(Input):
    '''
    Parse data file

    Parameters:
        Input   Strings of raw data

    Returns:
       A collection of graphs in edge weight format
    '''
    def generate_rows():
        '''
        Iterate through rows of input, converting each to a list of ints
        '''
        for row in Input:
            yield [int(i) for i in row.split()]

    return [data if len(data)>1 else data[0] for data in generate_rows() if len(data) > 0]


def extract_graphs(data):
    '''
    Extract a set of graphs from input

    Parameters:
       data A collection of graphs in edge weight format
    '''
    Edges  = []
    for row in data[1:]:
        if len(row)==2:
            if len(Edges)>0:
                yield Edges
                Edges=[]
            Edges.append(row)
        else:
            assert len(row)==3
            Edges.append(row)
    yield Edges

if __name__=='__main__':
    start = time()
    parser = ArgumentParser(__doc__)
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print ([n for n,_,_ in [bf(edges) for edges in
                                extract_graphs( [2,
                                                 [4, 5],
                                                 [1, 4, 4],
                                                 [4, 2, 3],
                                                 [2, 3, 1],
                                                 [3, 1, 6],
                                                 [2, 1, -7],

                                                 [3, 4],
                                                 [1, 2, -8],
                                                 [2, 3, 20],
                                                 [3, 1, -1],
                                                 [3, 2, -30]])]])

    if args.rosalind:
        Input = read_strings(f'data/rosalind_{basename(__file__).split(".")[0]}.txt')
        Result = [n for n,_,_ in [bf(edges,s=0) for edges in extract_graphs(extract_data(Input))]]
        Formatted = ' '.join(str(r) for r in Result)
        print (Formatted)
        with open(f'{basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Formatted}\n')

    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
