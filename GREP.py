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

def grow(cycles,E):
    def grow_one(c):
        return [c+[suffix] for (prefix,suffix) in E if prefix==c[-1]]
    return [cycle for c in cycles for cycle in grow_one(c)]
def is_cycle(S):
    return S[0]==S[-1]
def fmt(cycle):
    return ''.join([c[0] for c in cycle])
def grep(S):
    B,E=dbru(S,include_revc=False)
    print (B)
    print (E)
    e=S[0]
    cycles=[]
    runs=[[e[0:-1],e[1:]]]
    while len(cycles)==0:# or len(cycles[-1])<16:
        runs=grow(runs,E)
        new_cycles=[r for r in runs if is_cycle(r) and len(r)==len(S)-1]
        for c in new_cycles:
            runs.remove(c)
            cycles.append(c)
    cc= [fmt(c) for c in cycles]
    return list(set(cc))
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