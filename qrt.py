#!/usr/bin/env python
# Copyright (C) 2020 Greenweaves Software Limited

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

''' QRT Incomplete Characters'''

from argparse import ArgumentParser
from os.path import basename
from time import time

from helpers import read_strings
from   phylogeny import qrt, expand_character



def format(quartet):
    return f'{{{quartet[0]}, {quartet[1]}}} {{{quartet[2]}, {quartet[3]}}}'

if __name__=='__main__':
    start = time()
    parser = ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        for quartet in qrt(
                        ['cat', 'dog', 'elephant', ' ostrich', 'mouse', 'rabbit', 'robot'],
                        [expand_character(character) for character in ['01xxx00', 'x11xx00', '111x00x']]):
            print (format(quartet))

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{basename(__file__).split(".")[0]}.txt')

        with open(f'{basename(__file__).split(".")[0]}.txt','w') as f:
            for quartet in qrt(Input[0].split(),
                               [expand_character(character) for character in Input[1:]]):
                Result = format(quartet)
                print (Result)
                f.write(f'{Result}\n')

    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
