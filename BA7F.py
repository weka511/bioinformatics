#!/usr/bin/env python

# Copyright (C) 2017-2023 Greenweaves Software Limited
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

'''BA7F Implement SmallParsimony'''

from argparse  import ArgumentParser
from os.path   import basename
from time      import time
from helpers   import read_strings
from rosalind  import LabeledTree,hamm
from phylogeny import SmallParsimony


def print_assignments(assignments):
    '''
    Used to format output
    '''
    assignments.nodes.sort()

    for node in assignments.nodes:
        if node in assignments.edges:
            for edge in assignments.edges[node]:
                end,weight=edge
                if node in assignments.labels and end in assignments.labels:
                    print ('{0}->{1}:{2}'.format(assignments.labels[node],
                                                 assignments.labels[end],
                                                 hamm(assignments.labels[node],assignments.labels[end])))

if __name__=='__main__':
    start = time()
    parser = ArgumentParser(__doc__)
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        N = 4
        T = LabeledTree.parse(4,
                               ['4->CAAATCCC',
                                '4->ATTGCGAC',
                                '5->CTGCGCTG',
                                '5->ATGGACGA',
                                '6->4',
                                '6->5'],
                               bidirectional=False
                               )
        score,assignments=SmallParsimony(T)
        print (score)
        print_assignments(assignments)


    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{basename(__file__).split(".")[0]}.txt')
        T = LabeledTree.parse(int(Input[0]),Input[1:],bidirectional=False)
        score,assignments = SmallParsimony(T)
        print (score)
        print_assignments(assignments)

    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
