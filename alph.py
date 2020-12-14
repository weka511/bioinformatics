#   Copyright (C) 2020 Greenweaves Software Limited

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

#  ALPH  Alignment-Based Phylogeny 

import argparse
import os
import time
from helpers import read_strings
from newick  import newick_to_adjacency_list
from numpy   import argmin

def alph(T,Alignment,Alphabet=['A','T','C','G','-']):
    
    def create_fixed_alignments():
        Leaves = {}
        k      = None
        for i in range(0,len(Alignment),2):
            Leaves[Alignment[i][1:]] = Alignment[i+1]
            if k==None:
                k=len(Alignment[i+1])
            else:
                assert k==len(Alignment[i+1])
        return k,Leaves
    
    def SmallParsimony(l):
    
        def delta(i,j):
            return 0 if i==j else 1
    
        def is_ripe(v):
            for child in Adj[v]:
                if not Tagged[child]: return False
            return True 
        
        def find_ripe(Open):
            Ripe   = []
            Unripe = []
            for v in Open:
                if is_ripe(v):
                    Ripe.append(v)
                else:
                    Unripe.append(v)
            return Ripe,Unripe
        
        def get_distance(v,k):
            return sum([min([s[child][i] + delta(i,k)  for i in range(len(Alphabet))]) for child in Adj[v]])
 
        def backtrack(root,s):
            score = 0
            Open  = [root]
            ks    = {}
            while len(Open)>0:
                v = Open.pop(0)
                for child in Adj[v]:
                    Open.append(child)
                index = argmin([s[v][k] for k in range(len(Alphabet))])
                if score==0:
                    score += s[v][index]
                ks[v] = index
            return score,ks
        
        s      = {}
        Tagged = {}
        Open   = []
        for v in Adj.keys():
            if v in Leaves:
                char      = Leaves[v][l]
                s[v]      = [0 if Alphabet[k]==char else float('inf') for k in range(len(Alphabet))]
                Tagged[v] = True
            else:
                Tagged[v] = False
                Open.append(v)
        Ripe,Open = find_ripe(Open)
        while len(Ripe)>0:
            for v in Ripe:
                s[v] = [get_distance(v,k) for k in range(len(Alphabet))]
                #print (v, s[v])
                Tagged[v] = True
            Ripe,Open = find_ripe(Open)
        assert len(Open)==0
        return backtrack(v,s)
    
 
        
    Adj     = newick_to_adjacency_list(T)   
    L,Leaves = create_fixed_alignments()
    assert len([node for node,value in Adj.items() if len(value)==0 and node not in Leaves])==0
    Assignment = {a:[] for a in Adj.keys()}
    
    d = 0
    for l in range(L):
        score,ks = SmallParsimony(l)
        d       += score
        for v,index in ks.items():
            Assignment[v].append(Alphabet[index])
            
    return d,[(f'>{a}',''.join(b)) for a,b in Assignment.items() if len(Adj[a])==0]
    

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('ALPH  Alignment-Based Phylogeny')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        d,Assignment = alph('(((ostrich,cat)rat,(duck,fly)mouse)dog,(elephant,pikachu)hamster)robot;',
                    ['>ostrich',
                     'AC',
                     '>cat',
                     'CA',
                     '>duck',
                     'T-',
                     '>fly',
                     'GC',
                     '>elephant',
                     '-T',
                     '>pikachu',
                     'AA'
                     ])
        print (d)
        for label,String in Assignment:
            print (label)
            print (String)        
        
  
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        d,Assignment = alph(Input[0],Input[1])
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            print (d)
            f.write(f'{line}\n')
            for label,String in Assignment:
                print (label)
                f.write(f'{label}\n')
                print (String)
                f.write(f'{Sringt}\n')                
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
