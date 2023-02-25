#!/usr/bin/env python

#   Copyright (C) 2020-2023 Simon Crase

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

'''
BA10D 	Compute the Probability of a String Emitted by an HMM
'''

from argparse import ArgumentParser
from os.path  import basename
from time     import time

import numpy as np

from helpers  import read_strings,create_hmm_from_strings
from hmm      import Likelihood

if __name__=='__main__':
    start = time()
    parser = ArgumentParser(__doc__)
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (Likelihood('xzyyzzyzyy',
                          'xyz',
                          'AB',
                          np.array([[0.303, 0.697],
                                    [0.831, 0.169]]),
                          np.array([[0.533, 0.065, 0.402],
                                    [0.342, 0.334, 0.324]])))

                          # {('A','A') : 0.303, ('A','B') :   0.697,
                           # ('B','A') : 0.831, ('B','B') :   0.169, },
                          # {('A','x') : 0.533, ('A','y') :   0.065, ('A','z') : 0.402,
                           # ('B','x') : 0.342, ('B','y') :   0.334,  ('B','z') : 0.324, }))


    if args.extra:
        Input,Expected                         = read_strings(f'data/OutcomeLikelihood.txt',init=0)
        xs,alphabet,States,Transition,Emission = create_hmm_from_strings(Input)
        print (Likelihood(xs,alphabet,States,Transition,Emission))

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{basename(__file__).split(".")[0]}.txt')
        xs,alphabet,States,Transition,Emission = create_hmm_from_strings(Input)
        Result = Likelihood(xs,alphabet,States,Transition,Emission)
        print (Result)
        with open(f'{basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')

    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
