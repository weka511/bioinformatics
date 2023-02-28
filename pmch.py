#!/usr/bin/env python
#    Copyright (C) 2017-2020  Greenweaves Software Limited
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
#    pmch 	Perfect Matchings and RNA Secondary Structures

from rosalind import RosalindException,verify_counts_complete_graph
import math
import argparse
import os
import time
from helpers import create_strings

#    NumberPerfectMatchings
#    Verify that string contains the same number of As as Us, and the same number of Cs as Gs;
#    hence, that it is possible to match all bases.

#    Then determin numbere of matches

def NumberPerfectMatchings(string):

    counts=verify_counts_complete_graph(string)
    return math.factorial(counts['A']) * math.factorial(counts['G'])

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (NumberPerfectMatchings('AGCUAGUCAU'))


    if args.extra:
        Input,Expected  = read_strings('data/....txt',init=0)
        ...

    if args.rosalind:

        Input  = create_strings(path='./data',fasta=True)
        Result = NumberPerfectMatchings(Input[0])
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')

    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')



