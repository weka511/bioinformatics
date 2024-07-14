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

''' SMGB Semiglobal Alignment'''

from argparse import ArgumentParser
from os.path import basename
from time import time
from helpers import read_strings
from fasta import FastaContent
from align import  smgb

if __name__=='__main__':
    start = time()
    parser = ArgumentParser('SMGB Semiglobal Alignment')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        score,s1,t1 = smgb('CAGCACTTGGATTCTCGG','CAGCGTGG')
        print (score)
        print (s1)
        print (t1)

    if args.rosalind:
        Input = read_strings(f'data/rosalind_{basename(__file__).split(".")[0]}.txt')
        fasta = FastaContent(Input)
        a,s = fasta[0]
        b,t = fasta[1]
        score,s1,t1 = smgb(s,t)

        with open(f'{basename(__file__).split(".")[0]}.txt','w') as f:
            print (score)
            print (s1)
            print (t1)
            f.write(f'{score}\n')
            f.write(f'{s1}\n')
            f.write(f'{t1}\n')

    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
