#!/usr/bin/env python

#   Copyright (C) 2020-2023 Greenweaves Software Limited

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

'''BA10J 	Solve the Soft Decoding Problem'''

from argparse import ArgumentParser
from os.path  import basename
from time     import time
from helpers  import create_hmm_from_strings, read_strings, format_list
from hmm      import float2str, SoftDecode

if __name__=='__main__':
    start = time()
    parser = ArgumentParser(__doc__)
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    parser.add_argument('--precision', default=4,                          help='Controls display of probabilities')
    args = parser.parse_args()
    if args.sample:
        print ('A\tB')
        for Ps in SoftDecode('zyxxxxyxzz',
                             ['xyz'],
                             ['A','B'],
                             {('A','A'):  0.911, ('A','B'):  0.089,
                              ('B','A'): 0.228, ('B','B'):  0.772 },
                             {('A','x'):0.356,('A','y'):   0.191, ('A','z'):   0.453,
                              ('B','x'):   0.04, ('B','y'):    0.467, ('B','z'):   0.493   }):
            print ('\t'.join(float2str(p,precision=args.precision) for p in Ps))



    if args.rosalind:
        Input                                  = read_strings(f'data/rosalind_{basename(__file__).split(".")[0]}.txt')
        xs,alphabet,States,Transition,Emission = create_hmm_from_strings(Input,
                                                                         sep = '\t')

        with open(f'{basename(__file__).split(".")[0]}.txt','w') as f:
            print (format_list(States))
            f.write(f'{format_list(States)}\n')
            for line in SoftDecode(xs,alphabet,States,Transition,Emission):
                print (format_list(line))
                f.write(f'{format_list(line)}\n')

    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
