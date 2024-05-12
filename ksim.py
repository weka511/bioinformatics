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

'''KSIM Finding All Similar Motifs'''

from argparse import ArgumentParser
from os.path import basename
from time import time
from helpers import read_strings
import numpy as np

def ksim(k,s,t):
   '''
   Finding All Similar Motifs

   Given: A positive integer k,
            a DNA string s of length at most 5 kbp representing a motif,
            and a DNA string t of length at most 50 kbp representing a genome.

           Return: All substrings t' of t such that the edit distance dE(s,t') is less than or equal to k.
           Each substring should be encoded by a pair containing its location in t followed by its length.

   '''
   m = len(s)
   n = len(t)
   matrix = np.full((m+1,n+1),-m*n,dtype=np.int64)
   matrix[0,:] = 0
   matrix[:,0] = 0
   indel = 1
   mismatch = 1
   for i in range(1,m+1):
      for j in range(1,n+1):
            matrix[i,j] = max(matrix[i-1,j] - indel,
                              matrix[i,j-1] - indel,
                              matrix[i-1,j-1] + (1 if s[i-1] == t[j-1] else - mismatch),
                              0)
   print (matrix)
   return [(-1,-1)]

if __name__=='__main__':
   start = time()
   parser = ArgumentParser(__doc__)
   parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
   parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
   args = parser.parse_args()
   if args.sample:
      for a,b in ksim(2,'ACGTAG','ACGGATCGGCATCGT'):
         print (a,b)

   if args.rosalind:
      Input  = read_strings(f'data/rosalind_{basename(__file__).split(".")[0]}.txt')

      with open(f'{basename(__file__).split(".")[0]}.txt','w') as f:
         for a,b in ksim(int(Input[0]),Input[1],Input[2]):
            print (f'{a} {b}')
            f.write(f'{a} {b}\n')

   elapsed = time() - start
   minutes = int(elapsed/60)
   seconds = elapsed - 60*minutes
   print (f'Elapsed Time {minutes} m {seconds:.2f} s')
