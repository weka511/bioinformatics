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

def mer(n,A,m,B):
    i=0
    j=0
    C=[]
    while j<n and i<m:
        if A[j]<B[i]:
            C.append(A[j])
            j+=1
        else:
            C.append(B[i])
            i+=1
    while j<n:
        C.append(A[j])
        j+=1
    while i<m:
        C.append(B[i])
        i+=1
    return C

if __name__=='__main__':

    with open('c:/Users/Weka/Downloads/rosalind_mer.txt') as f:
        n=0
        A=[]
        m=0
        B=[]
        i=0
        for line in f:
            text=line.strip()
            if i==0:
                n=int(text)
            elif i==1:
                A=[int(t) for t in text.split(' ')]
            elif i==2:
                m=int(text)
            else:
                B=[int(t) for t in text.split(' ')]
            i+=1

    print (' '.join([str(r) for r in mer(n,A,m,B)]))
