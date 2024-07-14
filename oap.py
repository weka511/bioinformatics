#!/usr/bin/env python
#    Copyright (C) 2019-2024 Greenweaves Software Limited
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
#    along with this program.  If not, see <http://www.gnu.org/licenses/>

'''  OAP Overlap Alignment'''

from argparse import ArgumentParser
from os.path import basename
from time import time
from sys import version
from helpers import read_fasta
from align import oap

if __name__=='__main__':
    start = time()
    parser = ArgumentParser('OAP Overlap Alignment')
    parser.add_argument('--sample',    default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind',  default=False, action='store_true', help='process Rosalind dataset')
    parser.add_argument('--version',   default=False, action='store_true', help='Get version of python')
    args = parser.parse_args()

    if args.version:
        print (f'{version}')

    if args.sample:
        d,s1,t1 = oap('CTAAGGGATTCCGGTAATTAGACAG',
                      'ATAGACCATATGTCAGTGACTGTGTAA')

        print ('{0}'.format(d))
        print (s1)
        print (t1)

    if args.rosalind:
        Data  = read_fasta(f'data/rosalind_{basename(__file__).split(".")[0]}.txt')
        d,s1,t1 = oap(Data[0],Data[1])
        print ('{0}'.format(d))
        print (s1)
        print (t1)
        with open(f'{basename(__file__).split(".")[0]}.txt','w') as o:
            o.write('{0}\n'.format(d))
            o.write('{0}\n'.format(s1))
            o.write('{0}\n'.format(t1))

    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
