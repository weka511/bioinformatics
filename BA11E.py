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

'''BA11E 	Sequence a Peptide'''

from argparse import ArgumentParser
from os.path import  basename,splitext,join
from time import time
import numpy as np
from matplotlib.pyplot import figure, show
from helpers import read_strings
from spectrum import invert
from reference_tables import integer_masses, test_masses

def SequencePeptide(spectral, protein_masses = integer_masses):
    '''
    BA11E Sequence a Peptide

    Input: A spectral vector S.

    Returns: A peptide with maximum score against S. For masses with more than one amino acid, any choice may be used.
    '''
    inverse_masses = invert(protein_masses)
    masses = sorted([mass for mass in inverse_masses])
    Spectrum = np.array([0] + spectral)
    Scores = np.full_like(Spectrum,-np.inf,dtype=np.float64)
    m = len(Scores)
    Scores[0] = 0
    Choices = {}
    for i in range(1,m):
        Predecessors = [i-k for k in masses if i>=k]
        Candidates = [Scores[j] for j in Predecessors]
        if len(Predecessors)>0:
            k = np.argmax(Candidates)
            Scores[i] = Spectrum[i] + Candidates[k]
            Choices[i] = k

    imax = np.argmax(Scores)
    max_score = Scores[imax]
    n = len(Scores) - 1
    while Scores[n]<max_score:
        n -= 1
    peptide = []
    while n > 0:
        i = Choices[n]
        peptide.append(masses[i])
        n -= masses[i]
    return ''.join(inverse_masses[i] for i in peptide[::-1])

if __name__=='__main__':
    start = time()
    parser = ArgumentParser('BA11E 	Sequence a Peptide ')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (SequencePeptide([0, 0, 0, 4, -2, -3, -1, -7, 6, 5, 3, 2, 1, 9, 3, -8, 0, 3, 1, 2, 1, 0],
               protein_masses=test_masses))

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')

        Result = None
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                f.write(f'{line}\n')

    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
