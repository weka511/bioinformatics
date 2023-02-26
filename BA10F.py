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

'''Construct a Profile HMM with Pseudocounts'''

from argparse import ArgumentParser
from os.path import basename
from time import time
from helpers import read_strings
from hmm import ConstructProfileHMM,formatEmission,formatTransition

if __name__=='__main__':
    start = time()
    parser = ArgumentParser(__doc__)
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
        Input  = read_strings(f'data/rosalind_{basename(__file__).split(".")[0]}.txt')

        Transition,Emission,States = ConstructProfileHMM(float(Input[0].split()[0]),
                                                         Input[2].split(),
                                                         Input[4:],
                                                         sigma=float(Input[0].split()[1]))
        with open(f'{basename(__file__).split(".")[0]}.txt','w') as f:
            for row in formatTransition(Transition,States,precision=args.precision):
                    print (row)
                    f.write(f'{row}\n')
            print ('--------')
            f.write('--------\n')
            for row in formatEmission(Emission,States,Input[2].split(),precision=args.precision):
                print (row)
                f.write(f'{row}\n')

    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
