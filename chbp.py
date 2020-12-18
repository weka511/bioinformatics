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

#  chbp Character-Based Phylogeny

import argparse
import os
import time
from   helpers import read_strings
from phylogeny import chbp


def expand(s):
     return [int(c) for c in s]

if __name__=='__main__':
     start = time.time()
     parser = argparse.ArgumentParser('Character-Based Phylogeny')
     parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
     parser.add_argument('--abcde',   default=False, action='store_true', help="process Jonathan Dursi's dataset")
     parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
     args = parser.parse_args()
     if args.sample:
          print (chbp(
               ['cat', 'dog', 'elephant', 'mouse', 'rabbit', 'rat'],
               [expand('011101'),
               expand('001101'),
               expand('001100')] ))
        
     if args.abcde:
          clade = chbp(
               ['a', 'b', 'c', 'd', 'e'],
               [expand('00000'),
                expand('00110'),
                expand('01000'),
                expand('01110'),
                expand('01111')] )
          print (f'{clade}')        
  
     if args.rosalind:
          Input           = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
          species         = Input[0].split()
          character_table = [expand(row) for row in Input[1:]]
          clade           = chbp(species,character_table)
          print (clade)
          with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
               f.write(f'{clade}\n')
                
     elapsed = time.time()-start
     minutes = int(elapsed/60)
     seconds = elapsed-60*minutes
     print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
