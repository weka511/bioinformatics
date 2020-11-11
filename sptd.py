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

#  SPTD Phylogeny Comparison with Split Distance

import argparse
import os
import time
from   helpers import read_strings
from   phylogeny import parse,create_adj
    
def sptd(species,newick1,newick2):
    def replace_leaves(adj):
        return {parent:sorted([seiceps[child] if child in seiceps else child for child in children]) for parent,children in adj.items() }
    def edges(adj):
        for parent,children in adj.items():
            for child in children:
                if child >= n:
                    #print (parent,child)
                    yield parent,child
    def splits(adj,min_size=2,terminator=None):
        def find_leaves(node,path=[]):
            for child in adj[node]:
                if child<n:
                    path.append(child)
                else:
                    find_leaves(child,path=path)            
  
        for parent,child in edges(adj):
            s1 = []
            find_leaves(child,s1)#[leaf for leaf in find_leaves(child)]
            if len(s1)<min_size: continue
            s2 = [leaf for leaf in range(n) if not leaf in s1]
            yield sorted(s1),sorted(s2)
        if terminator!=None:
            yield terminator
           
    def ds(adj1,adj2):
        shared = 0
        splits1 = sorted([s for s,_ in splits(adj1)])
        splits2 = sorted([s for s,_ in splits(adj2)])
        k1 = 0
        k2 = 0
        i1 = splits1[k1]#next(splits1)
        i2 = splits2[k2]#next(splits2)      
        while k1<len(splits1) and k2<len(splits2): 
            #if len(i1)==0: break  
            #if len(i2)==0: break
   
            if i1==i2:
                shared+=1
                k1+=1
                k2+=1
                if k1<len(splits1) and k2<len(splits2):
                    i1 = splits1[k1]
                    i2 = splits2[k2]
                 #i1,_ = next(splits1)
                #i2,_ = next(splits2)
            elif i1<i2:
                k1+=1
                if k1<len(splits1):
                    i1 = splits1[k1]            
                #i1,_ = next(splits1)         
            else:
                k2+=1
                if k2<len(splits2):
                    i2 = splits2[k2]  

                #i2,_ = next(splits2) 
        return 2*(n-3)- 2* shared
    
    n       = len(species)
    seiceps = {species[i]:i for i in range(n)}
    tree1   = replace_leaves(create_adj(parse(newick1,start=n)))
    tree2   = replace_leaves(create_adj(parse(newick2,start=n)))

    return ds(tree1,tree2)

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print(sptd('dog rat elephant mouse cat rabbit'.split(),
             '(rat,(dog,cat),(rabbit,(elephant,mouse)));',
             '(rat,(cat,dog),(elephant,(mouse,rabbit)));)'))
        
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = sptd(Input[0].split(), Input[1], Input[2])
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')
                
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
