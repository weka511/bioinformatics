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
    def partition(B,E):
        B1=set()
        (prefix,suffix)=E[0]
        found=True
        while found:
            found=False
            B1.add(prefix)
            for (a,b) in E:
                if not found and a==suffix:
                    prefix,suffix=a,b
                    found=True
        B1.add(suffix)           
        
        B2=set([b for b in B if not b in B1])
        
        E1=[]
        E2=[]
        for prefix,suffix in E:
            if prefix in B1 or suffix in B1:
                E1.append((prefix,suffix))
            else:
                E2.append((prefix,suffix))
                
        return (list(B1),E1,B2,E2)

    k=len(S[0])
    B,E=r.dbru(S,include_revc=True)
    print (B)
    print (E)
    (B1,E1,B2,E2)=partition(B,E)
    print (B1)
    print (E1)
    print (B2)
    print (E2)
    fragments=[]
    

    return []

if __name__=='__main__':
    gasm(['AATCT','TGTAA','GATTA','ACAGA'])
    #S=[]
    #with open('c:/Users/Weka/Downloads/rosalind_gasm.txt') as f:
        #for line in f:
            #S.append(line.strip())     
    #print (gasm(S))    