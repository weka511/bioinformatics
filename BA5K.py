#!/usr/bin/env python

#    Copyright (C) 2019-2020 Greenweaves Software Limited
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

#    BA5K Find a Middle Edge in an Alignment Graph in Linear Space

import argparse
import os
import time
from   helpers import read_strings
from   align import FindMiddleEdge

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA5K Find a Middle Edge in an Alignment Graph in Linear Space')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',     default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()

    if args.sample:
        print (FindMiddleEdge('PLEASANTLY','MEASNLY'))
        # (4, 3) (5, 4)

    if args.extra:
        Input,Expected = read_strings(f'data/middle_edge.txt',init=0)
        ((i,j),(k,l))  = FindMiddleEdge(Input[0],Input[1])
        print (f'Calculated: (({i},{j}),({k},{l}))')
        print (f'Expected {Expected[0]}')
        # Expect (512,510)(513,511)

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        ((i,j),(k,l))  = FindMiddleEdge(Input[0],Input[1])
        print (f'({i},{j}) ({k},{l})')
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'({i},{j}) ({k},{l})\n')

    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
