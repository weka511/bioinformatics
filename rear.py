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
from fragile import GreedySorting

def sorted(S):
    for i in range(1,len(S)):
        if S[i-1]>S[i]: return False
    return True

# snarfed from https://stackoverflow.com/questions/3755136/pythonic-way-to-check-if-a-list-is-sorted-or-not

def isSorted(x, key = lambda x: x): 
    return all([key(x[i]) <= key(x[i + 1]) for i in range(len(x) - 1)])

# http://www.cs.utoronto.ca/~brudno/csc2417_09/ElementaryHP.pdf
def bergeron(s):
    def get_oriented(S):
        return [(i,j,S[i],S[j]) for j in range(len(S)) for i in range(0,j) if abs(S[i]+S[j])==1 and S[i]*S[j]<0]
    def reverse_segment(S):
        return [-s for s in S[::-1]]    
    def reverse(S,pair):
        i,j,pi_i,pi_j = pair
        if pi_i+pi_j>0:
            return S[:i+1] + reverse_segment(S[i+1:j+1]) + S[j+1:]
        else:
            return S[:i+2] + reverse_segment(S[i+2:j] + S[j:])    
    def get_score(S):
        return len(get_oriented(S))
    framed = [0] + s + [len(s)+1]
    d = 0
    for i in range(20):
        if isSorted(framed):
            return d,framed[1:-2]
        oriented_pairs = get_oriented(framed)
        
        d+=1

def rear(s1,s2,sort=GreedySorting):

    if isSorted(s2): return sort(s1)
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
    
    it = iter(data)
    result = []
    for a in it:
        print ('----')
        result.append(rear(a,next(it)))
    print (result)
    #with open (r'C:\Users\Simon\Downloads\rosalind_rear.txt') as f:
        #result = []
        #for line in f:
            #if i%3==0:
                #original=parse(line)
            #if i%3==1:
                #result.append(rear(original,parse(line)))
                #print("--- {0} seconds ---".format(time.time() - start_time))
            #i+=1
    
        #print(result)
        