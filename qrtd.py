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

#  QRTD  Quartet Distance

import argparse
import os
import time
from   helpers   import read_strings
from   phylogeny import  parse

def create_adj(tree,indices=None):
    adj = {}
    def dfs(tree):
        id       = tree['id']
        name     = tree['name']
        children = tree['children']
        parentid = tree['parentid']
        if len(name)==0:
            adj[id]=[]
        if parentid>-1:
            adj[parentid].append(id if len(name)==0 else indices[name])     
        for child in children:
            dfs(child)
    dfs(tree)
    return adj

def create_edges(adj,n=None):
        return [(a,b) for a,children in adj.items() for b in children if b>=n]

# extract_quartets

def extract_quartets(edges,adj,n=None):
    
    # get_leaves
    #
    # Get all leaves in graph 'adj' that are below some specified node
    def get_leaves(node):
        # dfs
        #
        # Conduct depth first search, looking for leaves
        def dfs(u):
            if u>=n:
                for v in adj[u]:
                    dfs(v)                
            else:
                leaves.append(u)
                  
        leaves = []
        dfs(node)
        return leaves
    
    def split(a,b):
        s_b = leaves[b]
        s_a = [s for s in all_leaves if s not in s_b]
        return [(s_a[j],s_a[i],s_b[l],s_b[k]) for i in range(len(s_a))
                for j in range(i) 
                for k in range(len(s_b)) 
                for l in range(k)]
    
    all_leaves     = list(range(n))
    leaves         = {x:sorted(get_leaves(x)) for x in adj.keys()}
    internal_edges = [(a,b) for a,b in edges if b>=n]
    return [q  for a,b in internal_edges for q in split(a,b)]

def get_matches(quartets1,quartets2):
    i       = 0
    j       = 0
    matches = 0
    while i < len(quartets1) and j < len(quartets2):
        if quartets1[i]==quartets2[j]:
            matches += 1
            i       += 1
            j       += 1
        elif quartets1[i]<quartets2[j]:
            i += 1
        else:  # quartets1[i]>quartets2[j]
            j+=1
            
    return matches
            
# qrtd
#
# Given: A list containing n taxa  and two unrooted binary trees T1 and T2 on the given taxa. 
#        Both T1 and T2 are given in Newick format.
#
# Return: The quartet distance dq(T1,T2)

def qrtd(species,T1,T2):
    n         = len(species)
    indices   = {species[i]:i for i in range(n)}
    tree1     = parse(T1,start=n)
    adj1      = create_adj(tree1,indices=indices)
    edges1    = create_edges(adj1,n=n)
    quartets1 = sorted(set(extract_quartets(edges1,adj1,n=n)))
    print (len(quartets1), quartets1)
    
    tree2     = parse(T2,start=n)
    adj2      = create_adj(tree2,indices=indices)
    edges2    = create_edges(adj2,n=n)
    quartets2 = sorted(set(extract_quartets(edges2,adj2,n=n)))
    print (len(quartets2),quartets2)           
    return len(quartets1) + len(quartets2) - 2*get_matches(quartets1,quartets2)

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print(qrtd(
            'A B C D E'.split(),
            '(A,C,((B,D),E));',
            '(C,(B,D),(A,E));'            
        ))
        
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = qrtd(Input[0].split(), Input[1], Input[2])
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
