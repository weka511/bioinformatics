#!/usr/bin/env python
#    Copyright (C) 2019 Greenweaves Software Limited
#
#    This is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This software is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>
#
#    PS Partial sort

from hea import hea,siftUp,heapsort,is_heap


def ps(n,A,k):
    result = A[0:k]
    result = hea(k,result)
    j = k
    while j<n:
        result.append(A[j])
        for l in range(1,k):
            siftUp(0,l,result)
        result.pop(0)
        j+=1
    return heapsort(k,result)

if __name__=='__main__':
    #print (ps(10,[4, -6, 7, 8, -9, 100, 12, 13, 56, 17],3))
    with open('/Users/Simon/Downloads/rosalind_ps(3).txt') as f:
        n=int(f.readline().strip())
        A = [int(a) for a in f.readline().strip().split(' ')]
        k=int(f.readline().strip())
        print (' '.join(str(a) for a in ps(n,A,k)))
