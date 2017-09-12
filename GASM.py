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

def gasm(S):
    k=len(S[0])
    B,E=r.dbru(S)
    print (B)
    print (E)
    fragments=[]
    (prefix,suffix)=E[0]
    while len(fragments)<len(E):
        fragments.append(prefix)
        for (a,b) in E:
            if a==suffix:
                prefix,suffix=a,b
                break

    superstring=[]
    for i in range(len(S)):
        if i==0:
            superstring.append(fragments[i])
        else:
            superstring.append(fragments[i][-1])

    return ''.join(superstring)

if __name__=='__main__':
    #print(gasm(['AATCT','TGTAA','GATTA','ACAGA']))
    S=[]
    with open('c:/Users/Weka/Downloads/rosalind_gasm.txt') as f:
        for line in f:
            S.append(line.strip())     
    print (gasm(S))    