# Copyright (C) 2020 Greenweaves Software Limited

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

# SDAG Shortest Paths in DAG

import argparse
import os
import time
from   helpers import read_strings
from   align   import topological_order

def create_adjacency(edges):
    m,n       = edges[0]

    product = {}
    weights = {}
    
    for a in range(1,m+1):
        product[a]   = []
        
    for a,b,w in edges[1:]:
        product[a].append(b) 
        weights[(a,b)] = w
    
    for a in product.keys():
        product[a]=sorted(list(set(product[a])))  
        
    return m,n,product,weights

def sdag(m,adjacency,weights):
    t    = topological_order(adjacency.copy())
    D    = [None]*(m+1)
    D[1] = 0
    for i  in t:
        if D[i] == None: continue
        for j in adjacency[i]:
            if D[j]==None:
                D[j] = D[i] +weights[(i,j)]
            else:
                trial = D[i] +weights[(i,j)]
                if trial < D[j]:
                    D[j] = trial                 
    return D[1:]

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    
    if args.sample:
        Edges = [
            [5, 6],
            [2, 3, 4],
            [4, 3, -2],
            [1, 4, 1],
            [1, 5, -3],
            [2, 4, -2],
            [5, 4, 1]        
        ]
        
  
        
        m,n,adjacency,weights = create_adjacency(Edges)
        Lengths = sdag(m,adjacency,weights)
        print (' '.join(str(l) if l!= None else 'x' for l in Lengths))
        
    if args.extra:
        Input,Expected  = read_strings('data/....txt',init=0)
        trie = Trie(Input)
        Actual = None
        Expected.sort()
        print (len(Expected),len(Actual))
        diffs = [(e,a) for e,a in zip(Expected,Actual) if e!=a]
        print (diffs)
  
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = None
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                f.write(f'{line}\n')
                
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
