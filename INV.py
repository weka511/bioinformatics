# Copyright (C) 2017-2019 Greenweaves Software Limited

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

# INV Counting Inversions

import timeit, numpy as np

# inv
#
# Input: A positive integer and an array A[1..n] of integers
#
# Return: The number of inversions in A
#
# Divide and conquer
#
# Split array in half, solve for each array, then merge, counting swaps

def inv(n,A):
    
    # Merge the two halves of the original array
    
    def merge(head,tail):
        count       = 0
        sorted_list = []
        i           = 0
        j           = 0
        
        while i<len(head) and j < len(tail):
            if head[i]<=tail[j]:
                sorted_list.append(head[i])
                i+=1
            else:
                sorted_list.append(tail[j])
                j+=1
                count+=(len(head)-i)  
                
        return sorted_list + head[i:] + tail[j:],count
    
    # Split array, calculate count for each half, then merge
    
    def merge_sort(A):
        if len(A)==1: return A,0
            
        head,count_left          = merge_sort(A[:len(A)//2])
        tail,count_right         = merge_sort(A[len(A)//2:])
        sorted_list,count_merged = merge(head,tail)   
        
        return sorted_list,count_left + count_right + count_merged
    
    return merge_sort(A)


if __name__=='__main__':
    
    #print (inv(5,[-6, 1, 15, 8, 10]))
    
    with open('c:/Users/Simon/Downloads/rosalind_inv(3).txt') as f:
        start_time = timeit.default_timer()
        i=0
        n=0
        A=[]
        for line in f:
            text=line.strip()
            if i==0:
                n=int(text)
                print (n)
            else:
                A=[int(t) for t in text.split(' ')]
            i+=1
        _,nn=inv(n,A)
        print (nn)

        elapsed = timeit.default_timer() - start_time
        print (elapsed)