#!/usr/bin/env python

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

# ctbl Creating a Character Table  http://rosalind.info/problems/ctbl/

import argparse
import os
import time
from   helpers import read_strings
from   Bio import Phylo
from   io import StringIO
from   phylogeny import CharacterTable


def create_tree(treedata):
    handle   = StringIO(treedata)
    return Phylo.read(handle, "newick")

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        for ch in CharacterTable(create_tree('(dog,((elephant,mouse),robot),cat);')):
            print (ch)

    if args.rosalind:
        tree = Phylo.read(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt', 'newick')

        Result = CharacterTable(tree)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                print (line)
                f.write(f'{line}\n')

    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')

