#    Copyright (C) 2019-2020 Greenweaves Software Limited
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
#    along with this program.  If not, see <http://www.gnu.org/licenses/>
#
#    LING   Linguistic Complexity of a Genome

import argparse
import os
import time
from   helpers import read_strings

# ling_naive

def ling_naive(string,a=4):
    def possible(k):
        return a if k==1 else min(a**k,n-k+1)
    
    def actual(k):
        return len({string[i:i+k]:True for i in range(n-k+1)})
    
    n   = len(string)
    sub = [actual(k) for k in range(1,n+1)]
    m   = [possible(k) for k in range(1,n+1)]

    return sum(sub)/sum(m)

class Node:
    def __init__(self,k=None):
        self.Edges = {}

# explore
#
# k      4^k
# 8	65,536
# 9	262,144
#       
def explore(string,a=4):
    def possible(k):
        return a if k==1 else min(a**k,n-k+1) 
    
    n       = len(string)    
    Root    = Node(k=0)
    Current = [Root]
    
    Nodes=[]
    
    for k in range(10):
        for l in range(possible(k)):
            Nodes.append(Node())
        
    for i in range(n):
        if len(Root.Edges)==m: break
        if string[i] not in Root.Edges:
            Root.Edges[string[i]] = Node(k=1)
           
    for k in range(2,10):
        m = possible(k)
        Current = [edge for node in Current for edge in node.Edges.values()]
        x=0
    for i in range(n-k+1):
        if len(substrings)>=m: break
        substrings.add(string[i:i+k])
        print (len(substrings),substrings)
        
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (ling_naive('ATTTGGATT'))

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        explore(Input[0])
        #Result = ling(Input[0])
        #print (Result)
        #with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            #for line in Result:
                #f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s') 
