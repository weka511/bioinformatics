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


#def rear(s1,s2):
    #def reverse(s,i0,i1):
        #if i0<i1:
            #return s[:i0] + s[i0:i1][::-1] + s[i1:]
    #def generate_single_reversals(s):
        #return [reverse(s,i0,i1) for i0 in range(len(s)) for i1 in range(i0+1,len(s)+1)]
    #if s1==s2:
        #return 0
    #history=set()
    #history.add(str(s2))
    #candidates=[s2]
    #for depth in range(25):
        #print (depth, len(candidates), len(history))
        #for s in candidates:
            #if s1==s:
                #return depth        
        #next_generation=[]
        #for s in candidates:
            #for rr in generate_single_reversals(s):
                #key= str(rr)
                #if not key in history:
                    #next_generation.append(rr)
                    #history.add(key)
        ##for s in next_generation:
            ##if s1==s:
                ##return depth            
        #candidates=next_generation
    #return None

def rear(s1,s2):
    return -1



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
        result.append(rear(a,next(it)))
    print (result)
 #   with open (r'C:\Users\Simon\Downloads\rosalind_rear.txt') as f:
#    with open (r'rear.txt') as f:
        #for line in f:
            #if i%3==0:
                #original=parse(line)
            #if i%3==1:
                #result.append(rear(original,parse(line)))
                #print("--- {0} seconds ---".format(time.time() - start_time))
            #i+=1
    
    #print(result)
        