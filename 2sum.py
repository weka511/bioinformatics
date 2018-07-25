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
#    maj 	Majority Element

def f2sum(A):
    def order(a,b):
        return (a,b) if a<b else (b,a)
    indices = sorted(range(len(A)), key=lambda k: A[k])
    AA = [A[i] for i in indices]
    first = 0
    last = len(indices)-1
    while AA[first]<=0 or AA[last]>=0: 
        if AA[first]<0:
            first+=1
        if AA[last]>0:
            last-=1
    if AA[first]>0:
        first-=1
    if AA[last]<0:
        last+=1
    while 0<=first and last < len(AA): 
        if -AA[first]==AA[last]:
            return order(indices[first]+1,indices[last]+1)
        if abs(AA[first])<abs(AA[last]):
            first-=1
        else:
            last+=1
    return -1

if __name__=='__main__':
    #k=4 
    #m=5
    #As=[
        #[2 ,-3, 4, 10, 5],
        #[8 ,2, 4, -2, -8],
        #[-5, 2, 3, 2, -4],
        #[5, 4, -5, 6, 8]
       #]     
    m = -1
    k = -1
    As=[]
 
    with open (r'C:\Users\Simon\Downloads\rosalind_2sum.txt') as f:    
        for line in f:
            if k==-1:
                values=line.strip().split()
                k=int(values[0])
                m=int(values[1])
            else:
                As.append([int(v) for v in line.strip().split()])    

    for A in As:
        print(f2sum(A))

