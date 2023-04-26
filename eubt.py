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

'''eubt Enumerating Unrooted Binary Trees'''

from argparse  import ArgumentParser
from os.path   import basename
from time      import time
from helpers   import read_strings
from phylogeny import UnrootedBinaryTree

if __name__=='__main__':
    start = time()
    parser = ArgumentParser(__doc__)
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        for tree in UnrootedBinaryTree.Enumerate('dog cat mouse elephant'.split()):
            print (f'{tree};\n')

    if args.rosalind:
        Input = read_strings(f'data/rosalind_{basename(__file__).split(".")[0]}.txt')
        with open(f'{basename(__file__).split(".")[0]}.txt','w') as f:
            for tree in UnrootedBinaryTree.Enumerate(Input[0].split()):
                print (f'{tree};')
                f.write(f'{tree};\n')

    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
