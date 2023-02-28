#!/usr/bin/env python
# Copyright (C) 2020 Greenweaves Software Limited

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

#  PDPL Creating a Restriction Map

import argparse
import os
import time
from   helpers import read_strings

# pdpl
#
# Given: A multiset L containing n*(n-1)/s positive integers for some positive integer n
# Return: A set X containing n nonnegative integers such that L is set of positive differences.
#
# Idea: set X isn't unique, as we can add a constant offset to all elements. So we can always take
#       X to contain zero, and also max(L). If 0 is X, X must also contain many of the same elements
#       as L. We therefore try adding them one by one, then check to see whether this allows us to
#       generate an invalid value.

def pdpl(L):

    # compare
    #
    # Used to check solution-returns zero if the two lists are identical

    def compare(list1,list2):
        return sum([i1!=i2 for i1,i2 in zip(list1,list2)])

    # get_signature
    #
    # Takes a "census" of a string: all elemnts plus their counts

    def get_signature(L):
        signature = {}
        for i in L:
            if i in signature:
                signature[i]+=1
            else:
                signature[i]=1
        return signature

    # get_diffs
    #
    # Calculate difference multiset

    def get_diffs(X):
        return sorted([b-a for a in X for b in X if b>a])

    # isCompatible
    #
    # Determine whether differences of X are a subset of L
    def isCompatible(X):
        for x,count in get_signature(get_diffs(X)).items():
            if x not in signature: return False      # generated element that isn't in L
            if count>signature[x]: return False      # generated too many of some value that is in L
        return True

    signature = get_signature(L)
    X = [0,L[-1]]
    for i in L:
        if isCompatible(X+[i]):
            X.append(i)
    assert 0==compare(L,get_diffs(X)),'Inconsistent solution'
    return sorted(X)

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (pdpl([2, 2, 3, 3, 4, 5, 6, 7, 8, 10]))

    if args.rosalind:
        Input     = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        converted = [int(i) for i in Input[0].split()]
        Result    = pdpl(sorted(converted))
        Formatted = ' '.join(str(r) for r in Result)
        print (Formatted)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
                f.write(f'{Formatted}\n')

    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
