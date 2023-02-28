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

#  OSYM Isolating Symbols in Alignments

import argparse
import os
import time
from   helpers import read_strings

# osym
#
# Say that we have two strings s and t of respective lengths m and n and an alignment score.
# Let's define a matrix M corresponding to s and t by setting Mj,k equal to the maximum score of any
# alignment that aligns s[j] with t[k]. So each entry in M can be equal to at most the maximum
# score of any alignment of s and t
#
# Given: Two DNA strings s and t in FASTA format, each having length at most 1000 bp.
#
# Return: The maximum alignment score of a global alignment of s and t, followed by
# the sum of all elements of the matrix M corresponding to s and t that was defined above.
# Apply the mismatch score introduced in Finding a Motif with Modifications.

def osym(s,t, match=1,mismatch=-1):
  def get_score(a,b, match=1,mismatch=-1):
    return match if a==b else mismatch

  def score_string(s,t):
    return sum(get_score(a,b) for (a,b) in zip(s,t))
  return 0,0

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('OSYM Isolating Symbols in Alignments')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (osym('ATAGATA','ACAGGTA'))



    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')

        Result = None
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                f.write(f'{line}\n')

    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
