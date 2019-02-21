#    Copyright (C) 2018 Greenweaves Software Pty Ltd
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
#   3SUM

import numpy as np

def sum3(A,n):
    for p in range(n):
        for q in range(p+1,n):
            for r in range(q+1,n):
                if A[p] + A[q] + A[r] ==0:
                    return (p+1,q+1,r+1)
    return (-1,-1,-1)        

if __name__=='__main__':
    k = -1
    n = -1
    with open(r'c:/Users/Simon/Downloads/rosalind_3sum.txt') as f:
        for line in f:
            text = line.strip().split()
            if k==-1:
                k = int(text[0])
                n = int(text[1])
                print (k,n)
            else:
                A = np.zeros(n)
                i = 0
                for s in text:
                    A[i] = int(s)
                p,q,r = sum3(A,n)
                if p==-1:
                    print (-1)
                else:
                    print ('{0} {1} {2}'.format(p,q,r))
