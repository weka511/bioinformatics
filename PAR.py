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

def par(n,A):
    prefix=[]
    suffix=[A[0]]
    for a in A[1:]:
        if a<=A[0]:
            prefix.append(a)
        else:
            suffix.append(a)
    return prefix + suffix


if __name__=='__main__':

    with open('c:/Users/Weka/Downloads/rosalind_par.txt') as f:
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
   
    print (' '.join([str(r) for r in par(n,A)]))