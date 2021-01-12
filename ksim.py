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

#  KSIM Finding All Similar Motifs 

import argparse
import os
import time
from   helpers import read_strings
from   numpy   import zeros
from   snp     import create_trie

# ksim
#
# Finding All Similar Motifs
#
# Given: A positive integer k,
#        a DNA string s of length at most 5 kbp representing a motif,
#        and a DNA string t of length at most 50 kbp representing a genome.
#
#       Return: All substrings t' of t such that the edit distance dE(s,t') is less than or equal to k.
#       Each substring should be encoded by a pair containing its location in t followed by its length.
#
# The following code is an attempt to see wther we can execure the Manhattan algorithm for rosalind data
# within 5 minutes (300 seconds). It looks as if this isn't possible.
#
# 4900 iterations along t axis required 456 sec, so all of th will need 3934 sec.

def ksim(k,s,t):
   Trie = create_trie([t])
   x=0
if __name__=='__main__':
   start = time.time()
   parser = argparse.ArgumentParser('KSIM Finding All Similar Motifs')
   parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
   parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
   args = parser.parse_args()
   if args.sample:
      for a,b in ksim(2,'ACGTAG','ACGGATCGGCATCGT'):
         print (a,b)
        
   if args.rosalind:
      Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')

      with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
         for a,b in ksim(int(Input[0]),Input[1],Input[2]):
            print (f'{a} {b}')
            f.write(f'{a} {b}\n')
                
   elapsed = time.time() - start
   minutes = int(elapsed/60)
   seconds = elapsed - 60*minutes
   print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
