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

import MER

def ms(n,A):
    if n<2:
        return A
    else:
        n1=n//2
        A1=A[0:n1]
        B1=A[n1:]
        m1=len(B1)
        return MER.mer(n1,ms(n1,A1),m1,ms(m1,B1))
    
#print (ms(10,[20, 19, 35, -18, 17, -20, 20, 1 ,4, 4]))

with open('c:/Users/Weka/Downloads/rosalind_ms(1).txt') as f:
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
    #print (n,a)
    #print ()
    print (' '.join([str(r) for r in ms(n,a)]))