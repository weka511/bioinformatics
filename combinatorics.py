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

#    count_perfect_matchings
#
#    Used by CAT Catalan Numbers and RNA Secondary Structures
#    to count perfect matchings, i.e. every base must be linked to another

def count_perfect_matchings(seq):

    def count(indices):
        key = str(indices)
        if key in cache:
            return cache[key]
        
        result = 0
        if 0 != sum(seq[i] for i in indices if abs(seq[i])==1): return 0
        if 0 != sum(seq[i] for i in indices if abs(seq[i])==2): return 0
        if len(indices)==0: return 1
        if len(indices)==2: return 1
        i = min(indices)
        for j in range(i+1,max(indices)+1,2):
            if seq[i] + seq[j]!=0: continue # Is i-j a valid split of the data?
            I1,I2  = partition(indices,i,j)
            result += count(I1)*count(I2)
            
        cache[str(indices)]= result
        return result
    
    cache = {}
    return count(list(range(len(seq))))

# count_matchings
#
# Used by MOTZ 	Motzkin Numbers and RNA Secondary Structures to
# count possible matches, which need not be perfect. Includes the
# case where n nodes match

def count_matchings(seq):
    def wrapped_count(indices):
        def count():
            n      = len(indices)
            if n<2: return 1
            i      = min(indices)          
            count1 = wrapped_count(indices[1:]) #  If first node is not involved in a matching            
            count2 = 0                            #  If first node is involved in a matching
            for j in range(i+1,max(indices)+1):
                if seq[i] + seq[j]!=0: continue
                I1,I2   = partition(indices,i,j)
                count21 = wrapped_count(I1)
                count22 = wrapped_count(I2)
                count2 += (count21*count22)                
            return count1 + count2
        
        key = str(indices)
        if not key in cache:
            cache[key] = count()        
        return cache[key]
    cache = {}
    return wrapped_count(list(range(len(seq)))) 

# catmotz
#
# Count entries in bonding graph (probleme CAT, MOTZ, and RNAS).
# Function starts be translating sting to list of ints. Fow CAT and MOTZ,
# a bond it possible iff to tho end of the bond sum to zero. RNAS allows
# CG also - both potitive.
#
# Parameters:  s        The string
#              counter   Function that performs counting

def catmotz(s,counter=count_perfect_matchings):
    to_int = {'A':+1, 'U':-1, 'G':-2, 'C':+2}
    return counter([to_int[c] for c in s])
    