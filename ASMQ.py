'''
 ASMQ Assessing Assembly Quality with N50 and N75

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

def asmq_n(S,N):
    S.sort()
    total=sum(s for s in S)
    target=N*total/100
    #print(total, target)
    ss=0
    n=len(S)-1
    while ss<target:
        #print (S[n])
        ss+=S[n]
        n-=1
    return S[n+1]
    
def asmq(S,N=50):
    return asmq_n([len(s) for s in S],N)

if __name__=='__main__':
    S=['GATTACA','TACTACTAC','ATTGAT','GAAGA']
    print (asmq(S),asmq(S,N=75))
    #S=[]
    #with open('c:/Users/Weka/Downloads/rosalind_asmq.txt') as f:
        #for line in f:
            #S.append(line.strip())     
    #print (asmq(S))    