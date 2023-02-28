#!/usr/bin/env python
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

import random

def qs(A):
    def split(v):
        prefix=[]
        suffix=[]
        mid=[]
        for a in A:
            if a<v:
                prefix.append(a)
            elif a==v:
                mid.append(a)
            else:
                suffix.append(a)
        return (prefix,mid,suffix)

    if len(A)<2:
        return A
    v=random.choice(A)
    (prefix,mid,suffix)=split(v)
    return qs(prefix)+mid+qs(suffix)

#print (qs([5, -2, 4, 7, 8, -10, 11]))

if __name__=='__main__':
    import timeit
    start_time = timeit.default_timer()
    with open('c:/Users/Weka/Downloads/rosalind_qs.txt') as f:
        n=0
        A=[]
        i=0
        for line in f:
            text=line.strip()
            if i==0:
                n=int(text)
            elif i==1:
                A=[int(t) for t in text.split(' ')]
            i+=1
        print (qs(A))
        print ('Elapsed: {0} seconds'.format(timeit.default_timer() - start_time))
