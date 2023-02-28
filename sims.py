#!/usr/bin/env python
#    Copyright (C) 2020 Simon Crase
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
#    along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>
#
# SIMS Finding Mutated Motifs

import argparse
import os
import time
from   helpers import read_strings
from   align import sims
from  fasta import FastaContent

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('KSIM Finding All Similar Motifs')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        score,s1,t1 = sims('GCAAACCATAAGCCCTACGTGCCGCCTGTTTAAACTCGCGAACTGAATCTTCTGCTTCACGGTGAAAGTACCACAATGGTATCACACCCCAAGGAAAC',
                           'GCCGTCAGGCTGGTGTCCG')
        print (score)
        print (len(s1),s1)
        print (len(t1),len('GCCGTCAGGCTGGTGTCCG'),t1)

    if args.rosalind:
        Input       = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        fc          = FastaContent(Input)
        score,s1,t1 = sims(fc[0][1],fc[1][1])
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            print (f'{score}')
            f.write(f'{score}\n')
            print (f'{s1}')
            f.write(f'{s1}\n')
            print (f'{t1}')
            f.write(f'{t1}\n')

    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')

