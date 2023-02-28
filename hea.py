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
#    HEA buiding a heap

def parent(i):
    return (i-1)//2

def leftChild(i):
    return 2*i+1

def rightChild(i):
    return 2*i + 2

def swap(i,j,A):
    x    = A[i]
    y    = A[j]
    A[i] = y
    A[j] = x

def siftDown(start,end,A):
    root = start
    while leftChild(root)<=end:
        child     = leftChild(root)
        swap_with = root
        if A[swap_with]<A[child]:
            swap_with = child
        if child+1<=end and A[swap_with]<A[child+1]:
            swap_with = child+1
        if swap_with==root:
            return
        else:
            swap(root,swap_with,A)
            root = swap_with

def siftUp(start,end,A):
    child = end
    while child>start:
        par = parent(child)
        if A[par]<A[child]:
            swap(par,child,A)
            child = par
        else:
            return

def hea(n,A):
    start = parent(n - 1)
    while start >= 0:
        siftDown(start,n-1,A)
        start -=1

    return A

def heapsort(n,A):
    hea(n,A)
    end = n-1
    while end>0:
        swap(end,0,A)
        end-=1
        siftDown(0,end,A)
    return A

def is_heap(n,A):
    for i in range(1,n):
        if A[parent(i)]<A[i]:
            return i
    return 0

if __name__=='__main__':
    #print (heapsort(9,[2, 6, 7, 1, 3, 5, 4, 8, 9]))
    with open('/Users/Simon/Downloads/rosalind_hs(3).txt') as f:
        n=int(f.readline().strip())
        A = [int(a) for a in f.readline().strip().split(' ')]
        print (' '.join(str(a) for a in heapsort(n,A)))
