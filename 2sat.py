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

# 2SAT 2-Satisfiability

import argparse
import os
import time
from helpers import read_strings,create_list,parse_graphs
from graphs import scc,create_adj 
import tarjan

def to_int(s):
    parts=s.split()
    if len(parts)==1:
        return int(s)
    else:
        return [int(p) for p in parts]
    
def two_sat(problem):
    def represent(i):
        return 2*i if i>0 else 2*(-i)+1
    def tneserper(i):
        return i//2 if i%2==0 else -((i-1)//2)
    def convert(adj,op=tneserper):
        result={}
        for head,tail in adj.items():
            result[op(head)] = [op(i) for i in tail]
        return result
    def collect(index,adj,path=[]):
        if index in path: return path
        path.append(index)
        for i in adj[index]:
            path = collect(i,adj,path)
        return path
    
    n,m,clauses = problem
    edges       = [(represent(-a),represent(b)) for a,b in clauses] + [(represent(-b),represent(a)) for a,b in clauses]
    scc         = tarjan.tarjan(create_adj([[2*n+1,len(edges)]] + edges)) 
    scc_conv    = [[tneserper(node) for node in component] for component in scc]
    for component in scc_conv:
        for i in range(len(component)):
            for j in range(i+1,len(component)):
                if component[i]==-component[j]:
                    return 0,[]        
    scc_new = sorted(scc_conv,reverse=True,key=lambda x:len(x))
    return 1,scc_new[0]

def create_sets(data):
    product = []
    i    = 1
    while i<len(data):
        n,m = data[i]
        product.append((n,m, [(xx[0],xx[1]) for xx in data[i+1:i+1+m]]))
        i = i+1+m
    assert len(product)==data[0]
    return product
 
def Format(status,Solution):
    return f'{status} {" ".join(str(sol) for sol in Solution)}'
     
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:

        data = [2,
                [2, 4],
                [1 ,2],
                [-1, 2],
                [1, -2],
                [-1, -2],
                
                [3, 4],
                [1, 2],
                [2, 3],
                [-1, -2],
                [-2, -3]]
        
        for problem in create_sets(data):
            status,Solution = two_sat(problem )
            print (Format(status,Solution))

 
      
    if args.rosalind:
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
            Data = [to_int(s) for s in Input if len(s)>0]
            Sets = create_sets(Data)
            for problem in Sets:
                status,Solution = two_sat(problem )
                print (Format(status,Solution))
                f.write (f'{Format(status,Solution)}\n')                
            #print (Sets[0])
            #print (Data)
        #with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            #with open(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt') as g:
                #for problem in parse_graphs(g):
                    #print (problem)
                    #status,Solution = two_sat(problem )
                    #print (status,Solution)
                    #f.write (Format(status,Solution))
            #print (format_list([bip(g) for g in parse_graphs(f)]))
             
        #Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        
        #with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            #for line in Input:
                #print (line)
                #f.write(f'{line}\n')
                
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
