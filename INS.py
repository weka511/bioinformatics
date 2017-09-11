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

def ins(n,A):
    swap=0
    for i in range(1,n):
        k=i
        while k>0 and A[k]<A[k-1]:
            A[k-1],A[k]=A[k],A[k-1]
            k-=1
            swap+=1

    return (swap,A)

with open('c:/Users/Weka/Downloads/rosalind_ins.txt') as f:
    i=0
    n=0
    a=[]
    for line in f:
        text=line.strip()
        if i==0:
            n=int(text)
        else:
            a=[int(t) for t in text.split(' ')]
        i+=1
    print (ins(n,a))
