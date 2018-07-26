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
#   2SUM

def f2sum(A):
    def order(a,b):
        return (a,b) if a<b else (b,a)
    def best(pairs):
        x0,y0=pairs[0]
        for x1,y1 in pairs[1:]:
            if y1<y0:
                x0,y0=x1,y1
        return (x0,y0)
    indices = sorted(range(len(A)), key=lambda k: A[k])
    A_Sorted= [A[i] for i in indices]
    p = 0
    q = len(A_Sorted)-1
    pairs=[]
    while p < q and A_Sorted[p]<= 0 and A_Sorted[q]>=0:
        if A_Sorted[p] + A_Sorted[q]==0:
            pairs.append(order(indices[p],indices[q]))
            q-=1
            p+=1
        elif abs(A_Sorted[p])<abs( A_Sorted[q]):
            q-=1
        elif abs(A_Sorted[p])>abs( A_Sorted[q]):
            p+=1
 
    return best(pairs) if len(pairs)>0 else (-1,-1)

if __name__=='__main__':
    m = -1
    k = -1
    As=[]
 
    with open (r'C:\Users\Simon\Downloads\rosalind_2sum(5).txt') as f:    
        for line in f:
            if k==-1:
                values=line.strip().split()
                k=int(values[0])
                m=int(values[1])
            else:
                As.append([int(v) for v in line.strip().split()])    
                
    for A in As:
        i,j=f2sum(A)
        if i>-1:
            print (i+1,j+1)  #Rosalind used 1-based indices!
        else:
            print (i)

