'''
Copyright (C) 2017-2019 Greenweaves Software Limited

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

REAR 	Reversal Distance
'''
import time
from fragile import GreedySorting, kReverse,isSorted
from numpy import argmax


def leaderBoardSort(S,N=5):
    def get_all_reversals(S):
        def reverse(i,j):
            def reverse_segment(S):
                return S[::-1]       
            return S[:i] + reverse_segment(S[i:j+1]) + S[j+1:]    
        return [reverse(i,j) for j in range(len(S)) for i in range(j)]    
    #def get_score(s):
        #return sum([abs(s[i+1]-s[i]) for i in range(len(s)-1)])
    def get_breakpoints(S):
        return [1 if S[i+1]-S[i]>0 else -1 for i in range(len(S)-1) if abs(S[i+1]-S[i])>1] 
    #def reversals(s):
        #result = []
        #for i in range(1,len(s)):
            #for j in range(i+1,len(s)-1):
                #result.append(s[:i]+s[i:j+1][::-1]+s[j+1:])
        #return result
                
    def create_leaders(leaders):
        permutation = [get_all_reversals(s) for s,_ in leaders]
        new_list    = [(p,len(get_breakpoints(p))) for ps in permutation for p in ps]
        #print (len(new_list),new_list)
        min_breakpoint_count = min([c for (p,c) in new_list])
        #print (min_breakpoint_count)
        result               = []
        while len(result)<N:
            result = result + [(p,c) for (p,c) in new_list if c==min_breakpoint_count]
            min_breakpoint_count+=1
        return result     
        #sorted_list = sorted(new_list, key=lambda tup: tup[1])
        #return sorted_list if len(sorted_list)<=N else sorted_list[:N]
        
    reversalDistance = 0
    leaders    = [(S,len(get_breakpoints(S)))]
    
    for k in range(1,len(S)+1):
        leaders = create_leaders(leaders)
        reversalDistance+=1
        s,b = leaders[0]
        if b==0:
            if s[0]<s[1]:
                return reversalDistance
            else:
                return reversalDistance+1
    return reversalDistance    

def rear(s1,s2,sort=leaderBoardSort):
    if isSorted(s2):
        if isSorted(s1):
            return 0
        else:
            return sort(s1)
    if isSorted(s1): return sort(s2)
    return sort([s2.index(p)+1 for p in s1])    

if __name__=='__main__':
    def parse(line):
        return [int(c) for c in line.strip().split()]
    
    start_time = time.time()
    i=0
    original=''


    data = [
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        [3, 1, 5, 2, 7, 4, 9, 6, 10, 8],
        
        [3, 10, 8, 2, 5, 4, 7, 1, 6, 9],
        [5, 2, 3, 1, 7, 4, 10, 8, 6, 9],
        
        [8, 6, 7, 9, 4, 1, 3, 10, 2, 5],
        [8, 2, 7, 6, 9, 1, 5, 3, 10, 4],
        
        [3, 9, 10, 4, 1, 8, 6, 7, 5, 2],
        [2, 9, 8, 5, 1, 7, 3, 4, 6, 10],
        
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]    
    ]
    
    #it = iter(data)
    #result = []
    #for a in it:
        #print ('----')
        #result.append(rear(a,next(it)))
    #print (result)
    with open (r'C:\Users\Simon\Downloads\rosalind_rear(2).txt') as f:
        result = []
        for line in f:
            if i%3==0:
                original=parse(line)
            if i%3==1:
                result.append(rear(original,parse(line)))
                print("--- {0} seconds ---".format(time.time() - start_time))
            i+=1
    
        print(result)
        