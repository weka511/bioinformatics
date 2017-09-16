'''
 GREP Genome Assembly with Perfect Coverage and Repeats

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

from rosalind import dbru 

def create_lookup(B,E):
    F={}
    for b in B:
        F[b]=[f for (e,f) in E if e==b]
    return F

def grow(run,F):
    end = run[-1]
    if end in F:
        return [run + [rr] for rr in F[end]]
    else:
        return run

def is_cycle(r,F):
    return r[-1] in F

def is_complete(r,freqs):
    #print (format(r))
    #print (r)
    kmers=[]
    for i in range(len(r)):
        if i<len(r)-1:
            kmers.append(r[i]+r[i+1][0])
        else:
            kmers.append(r[i]+r[0][0])
    #print ('k',kmers)
    myfreqs={}
    for rr in kmers:
        if rr in myfreqs:
            myfreqs[rr]+=1
        else:
            myfreqs[rr]=1
            
    for k in myfreqs.keys():
        if (not k in freqs):
            print (k,format(r),myfreqs[k],'-')
            return False
        if myfreqs[k] != freqs[k]:
            print (k,format(r),myfreqs[k],freqs[k])
            return False
    return True

def format(r):
    return ''.join([rr[0] for rr in r] + [r[-1][-1]])

def grep(S):
    B,E=dbru(S,include_revc=False)
    F=create_lookup(B,E)
    Runs=[[S[0][0:-1],S[0][1:]]]

    N=len(Runs[0])
    for i in range(len(S)-3): #while len(Runs[0])<5:
        Runs=[g for run in Runs for g in grow(run,F)]
        Runs=[r for r in Runs if len(r)>N]
        N+=1

    cycles =[r for r in Runs if is_cycle(r,F)]

    freqs={}
    for key in S:
        freqs[key]=1+freqs[key] if key in freqs else 1
    #freqs['GC']+=1
    return [format(c) for c in cycles if is_complete(c,freqs)]
    
if __name__=='__main__':
    S=[
        'CAG',
        'AGT',
        'GTT',
        'TTT',
        'TTG',
        'TGG',
        'GGC',
        'GCG',
        'CGT',
        'GTT',
        'TTC',
        'TCA',
        'CAA',
        'AAT',
        'ATT',
        'TTC',
        'TCA'    
    ]
    
    #S=[]
    #with open('c:/Users/Weka/Downloads/rosalind_grep.txt') as f:
        #for line in f:
            #S.append(line.strip())     
    for s in grep(S):
        print (s)