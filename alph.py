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

#  ALPH  Alignment-Based Phylogeny

import argparse
import os
import time
from helpers import read_strings


from fasta   import FastaContent, fasta_out

from phylogeny import alph

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('ALPH  Alignment-Based Phylogeny')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        fc = FastaContent(['>ostrich',
                           'AC',
                           '>cat',
                           'CA',
                           '>duck',
                           'T-',
                           '>fly',
                           'GC',
                           '>elephant',
                           '-T',
                           '>pikachu',
                           'AA'
                           ])

        d,Assignment = alph('(((ostrich,cat)rat,(duck,fly)mouse)dog,(elephant,pikachu)hamster)robot;',fc.to_list())
        print (d)
        for label,String in Assignment:
            for line in fasta_out(label,String):
                print (line)


    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        fc     = FastaContent(Input[1:])
        d,Assignment = alph(Input[0],fc.to_list())
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            print (d)
            f.write(f'{d}\n')
            for label,String in Assignment:
                for line in fasta_out(label,String):
                    print (line)
                    f.write(f'{line}\n')

    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
