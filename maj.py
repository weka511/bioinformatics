#!/usr/bin/env python
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
#    maj 	Majority Element

def maj(A):
    A.sort()
    first = len(A)//2
    last  = first +1
    match = A[first]
    count = 0
    while first>=0 and A[first]==match:
        first-=1
        count+=1
    while last<len(A) and A[last]==match:
        last+=1
        count+=1
    return match if 2 * count > len(A) else -1

if __name__=='__main__':
    m = -1
    k = -1
    As=[]

    with open (r'data\rosalind_maj.txt') as f:
        for line in f:
            if k==-1:
                values=line.strip().split()
                k=int(values[0])
                m=int(values[1])
            else:
                As.append([int(v) for v in line.strip().split()])

    print (' '.join([str(maj(As[i])) for i in range(k)]))

