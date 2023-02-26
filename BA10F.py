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

#  BA10F 	Construct a Profile HMM with Pseudocounts

import argparse
import os
import time
from   helpers import read_strings
from   hmm import ConstructProfileHMM,formatEmission,formatTransition

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA10F 	Construct a Profile HMM with Pseudocounts ')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    parser.add_argument('--precision', default=3,                          help='Controls display of probabilities')

    args = parser.parse_args()
    if args.sample:
        Transition, Emission,StateNames = ConstructProfileHMM(0.358,
                                                         'ABCD',
                                                         ['ADA', 'ADA', 'AAA', 'ADC', '-DA', 'D-A'],
                                                         sigma=0.01)

        for row in formatTransition(Transition,StateNames,precision=args.precision):
            print (row)
        print ('--------')
        for row in formatEmission(Emission,StateNames,'ABCD',precision=args.precision):
            print (row)



    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')

        Transition,Emission,States = ConstructProfileHMM(float(Input[0]),
                                                         Input[2].split(),
                                                         Input[4:])

        Result = None
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                f.write(f'{line}\n')

    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
