# Copyright (C) 2017 Greenweaves Software Pty Ltd

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

import timeit

def inv(n,A):
    count=0
    for i in range(n):
        for j in range(i+1,n):
            if A[i]>A[j]:
                count+=1
    return count

#print (inv(5,[-6, 1, 15, 8, 10]))

with open('c:/Users/Weka/Downloads/rosalind_inv.txt') as f:
    start_time = timeit.default_timer()
    i=0
    n=0
    A=[]
    for line in f:
        text=line.strip()
        if i==0:
            n=int(text)
        else:
            A=[int(t) for t in text.split(' ')]
        i+=1
    nn=inv(n,A)
    elapsed = timeit.default_timer() - start_time
    print (nn)
    print (elapsed)