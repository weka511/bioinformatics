'''
 GASM Genome Assembly Using Reads

 Copyright (C) 2017 Greenweaves Software Pty Ltd

 This is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This software is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

'''
import rosalind as r

def split(S,k):
    def split1(s):
        return[s[i:i+k] for i in range(len(s)-k+1)]
    if k==len(S[0]):
        return S[:] 
    return[spl for s in S for spl in split1(s)]

def get_run(E):
    def search(s):
        for a,b in E:
            if a==s:
                return a,b
        else:
            return ([],[])
    prefix,suffix=E[0]
    R=[prefix]
    while len(suffix)>0:
        R.append(suffix)
        E.remove((prefix,suffix))
        prefix,suffix=search(suffix)      
    return R,E
def gasm(S):
    shortest_cycle=None
    for k in range(len(S[0]),2,-1):
        Sk= split(S,k)
        B,E=r.dbru(Sk)
        while(len(E)>0):
            R,E=get_run(E)
            if R[0]==R[-1]:
                cycle=''.join([r[0] for r in R[0:-1]])
                if shortest_cycle==None or len(cycle)<len(shortest_cycle):
                    shortest_cycle=cycle

        if shortest_cycle!=None:
            return shortest_cycle
    return shortest_cycle


if __name__=='__main__':
    #print(gasm(['AATCT','TGTAA','GATTA','ACAGA']))
    S=[]
    with open('c:/Users/Weka/Downloads/rosalind_gasm.txt') as f:
        for line in f:
            S.append(line.strip())     
    print (gasm(S))    