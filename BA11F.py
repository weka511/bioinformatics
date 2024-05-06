#!/usr/bin/env python

#   Copyright (C) 2020-2024 Greenweaves Software Limited

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

'''  BA11F 	Find a Highest-Scoring Peptide in a Proteome against a Spectrum'''

from argparse import ArgumentParser
from os.path import  basename,splitext,join
from time import time
import numpy as np
from   helpers import read_strings

def IdentifyPeptide(S):
    '''
    Find a peptide from a proteome with maximum score against a spectrum.

        Inputs:
           S         A space-delimited spectral vector
           Proteome  An amino acid string Proteome

        Returns:
            A peptide in Proteome with maximum score against S.
    '''
    pass

if __name__=='__main__':
    start = time()
    parser = ArgumentParser('BA11F 	Find a Highest-Scoring Peptide in a Proteome against a Spectrum ')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (IdentifyPeptide([0, 0, 0, 4, -2, -3, -1, -7, 6, 5, 3, 2, 1, 9, 3, -8, 0, 3, 1, 2, 1, 8],
                               'XZZXZXXXZXZZXZXXZ'))

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{basename(__file__).split(".")[0]}.txt')

        Result = None
        print (Result)
        with open(f'{basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                f.write(f'{line}\n')

    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
