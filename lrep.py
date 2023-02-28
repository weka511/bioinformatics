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

#  lrep Finding the Longest Multiple Repeat

import argparse
import os
import time
from   helpers import read_strings
from   numpy import argsort

def lrep(text,target_count,edges):

    # extract_entries_from_edges
    #
    # Convert edges from string to list of lists:
    #       -from
    #       -to
    #       - start of substring
    #       - length of substring

    def extract_entries_from_edges():
        return [[int(e.replace('node','')) for e in edge.split()] for edge in edges]

    # standardize_entry
    #
    # 1. If string includes '$', remove it (so it matches same substring as internal)
    # 2. Change start of substring to have zero offset as in Python

    def standardize_entry(entry):
        link1,link2,pos,length = entry
        if pos+length>len(text):
            length -=1
        return [link1,link2,pos-1,length]

    # remove_duplicates
    def remove_duplicates(entries):
        Result = []
        i = 0
        while i < len(entries):
            j = i + 1
            while j<len(entries) and entries[i][LENGTH]==entries[j][LENGTH]:
                j +=1
            same_length = sorted(entries[i:j],key=lambda x:x[POS])
            Result.append(entries[i])
            for k in range(i+1,j):
                if Result[-1][POS]!=entries[k][POS]:
                    Result.append(entries[k])
            i = j
        return Result

    POS    = 2 # Index in entry
    LENGTH = 3 # Index in entry
    n      = len(text)

    # Sort entries into descending order by length
    entries = remove_duplicates(
                   sorted([standardize_entry(entry) for entry in extract_entries_from_edges()],
                          key=lambda x:-x[LENGTH]))
    i = 0
    while True:
        j       = i
        strings = []
        while entries[j][LENGTH] == entries[i][LENGTH]:
            strings.append(text[entries[j][POS]: entries[j][POS] + entries[j][LENGTH]])
            j += 1
        if len(strings)>=target_count:
            indices = argsort(strings)
            for k in range(len(strings)):
                candidate = strings[indices[k]]
                m         = k + 1
                while m<  len(strings) and strings[indices[m]]==candidate:
                    if m-k >= target_count-1:  return candidate
                    m += k
                k = m
        i = j+1


if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (lrep('CATACATAC$',
                    2,
                    ['node1 node2 1 1',
                     'node1 node7 2 1',
                     'node1 node14 3 3',
                     'node1 node17 10 1',
                     'node2 node3 2 4',
                     'node2 node6 10 1',
                     'node3 node4 6 5',
                     'node3 node5 10 1',
                     'node7 node8 3 3',
                     'node7 node11 5 1',
                     'node8 node9 6 5',
                     'node8 node10 10 1',
                     'node11 node12 6 5',
                     'node11 node13 10 1',
                     'node14 node15 6 5',
                     'node14 node16 10 1']))

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')

        Result = lrep(Input[0], int(Input[1]), Input[2:])
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')

    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
