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

#  Combinatorics The mathematics of counting objects.

import argparse
import os
import time
from   helpers import read_strings



#    cat 	Catalan Numbers and RNA Secondary Structures
#    motz 	Motzkin Numbers and RNA Secondary Structures

# partition
#
# Ensure that matching are non-crossing by splitting set into two partitions,
# one between i and and j, one outside
#
#   Parameters:
#       indices The set (of indices) that are to be partitioned
#       i       One end of bond
#       j       The other end of bond

def partition(indices,i,j):
    I1 = []
    I2 = []

    for k in indices:
        if k==i: continue
        if k==j: continue
        if i<k and k <j:
            I1.append(k)
        else:
            I2.append(k)

    return (I1,I2)

# valid_bonds
#
# Used to iterate through valid pairs of bonds, all originating
# from first index

def valid_bonds(seq,indices):
    j = min(indices)
    for k in range(j+1,max(indices)+1):
        if seq[j] + seq[k]==0:
            yield(j,k)

# count_partitioned
#
# Split indices onto two subsets, corresponding to a bond

def count_partitioned(indices,i,j,count=lambda x:0):
    I1,I2  = partition(indices,i,j)
    return count(I1)*count(I2)

#    count_perfect_matchings
#
#    Used by CAT Catalan Numbers and RNA Secondary Structures
#    to count perfect matchings, i.e. every base must be linked to another

def count_perfect_matchings(seq):
    def wrapped_count(indices):
        def count(indices):
            if 0 != sum(seq[i] for i in indices if abs(seq[i])==1): return 0
            if 0 != sum(seq[i] for i in indices if abs(seq[i])==2): return 0
            if len(indices)==0: return 1
            if len(indices)==2: return 1
            return sum([count_partitioned(indices,i,j,count=count) for i,j in valid_bonds(seq,indices)])

        key = str(indices)
        if not key in cache:
            cache[key] = count(indices)
        return cache[key]
    cache = {}
    return wrapped_count(list(range(len(seq))))

# count_matchings
#
# Used by MOTZ 	Motzkin Numbers and RNA Secondary Structures to
# count possible matches, which need not be perfect. Includes the
# case where n nodes match



def count_matchings(seq):
    def wrapped_count(indices):
        def count(indices):

            n      = len(indices)
            if n<2: return 1
            return wrapped_count(indices[1:]) + \
                   sum([count_partitioned(indices,i,j,count=count) for i,j in valid_bonds(seq,indices)])

        key = str(indices)
        if not key in cache:
            cache[key] = count(indices)
        return cache[key]
    cache = {}
    return wrapped_count(list(range(len(seq))))

# catmotz
#
# Count entries in bonding graph (probleme CAT, MOTZ, and RNAS).
# Function starts be translating sting to list of ints. For CAT and MOTZ,
# a bond is possible iff the two endx of the bond sum to zero. RNAS allows
# UG also - both negative.
#
# Parameters:  s        The string
#              counter   Function that performs counting

def catmotz(s,
            counter=count_perfect_matchings,
            to_int = {'A': +1,
                      'U': -1,
                      'G': -2,
                      'C': +2}):
    return counter([to_int[c] for c in s])
