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

''' BA11H 	Compute the Size of a Spectral Dictionary'''

from  argparse import ArgumentParser
from os.path import basename
from time import time
from helpers import read_strings
from reference_tables import  TEST_MASSES
from spectrum import SizeSpectralDictionary

if __name__=='__main__':
    start = time()
    parser = ArgumentParser(__doc__)
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (SizeSpectralDictionary([4, -3, -2, 3, 3, -4, 5, -3, -1, -1, 3, 4, 1, 3],1,8,
                                                    protein_masses = TEST_MASSES))



    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{basename(__file__).split(".")[0]}.txt')

        Result = SizeSpectralDictionary( [int(s) for s in Input[0].split()],int(Input[1]),int(Input[2]))
        print (Result)
        with open(f'{basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')

    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
