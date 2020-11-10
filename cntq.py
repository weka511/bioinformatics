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

#  cntq Counting Quartets

import argparse
import os
import time
from   helpers import read_strings#, Tree

import re

# snarfed from https://stackoverflow.com/questions/51373300/how-to-convert-newick-tree-format-to-a-tree-like-hierarchical-object
def parse(newick):
    tokens = re.findall(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")
    
    def recurse(nextid = 0, parentid = -1): # one node
        thisid = nextid;
        children = []

        name, length, delim, ch = tokens.pop(0)
        if ch == "(":
            while ch in "(,":
                node, ch, nextid = recurse(nextid+1, thisid)
                children.append(node)
            name, length, delim, ch = tokens.pop(0)
        return {"id": thisid, "name": name, "length": float(length) if length else None, 
                "parentid": parentid, "children": children}, delim, nextid

    return recurse()[0]

def create_adj(tree):
    adj = {}
    def bfs(tree):
        id       = tree['id']
        name     = tree['name']
        children = tree['children']
        parentid = tree['parentid']
        if len(name)==0:
            adj[id]=[]
        if parentid>-1:
            adj[parentid].append(id if len(name)==0 else name)     
        for child in children:
            bfs(child)
    bfs(tree)
    return adj

def cntq(n,newick):
    def bfs(subtree,leaves):
        for node in adj[subtree]:
            if type(node)==str:
                leaves.append(node)
            else:
                bfs(node,leaves)
    def pair(leaves):
        for i in range(len(leaves)):
            for j in range(i+1,len(leaves)):
                yield [leaves[i],leaves[j]] if leaves[i]<leaves[j] else [leaves[j],leaves[i]]
                
    adj = create_adj(parse(newick))
    taxa = [leaf for children in adj.values() for leaf in children if type(leaf)==str]
    splitting_edges = [(key,child) for key,value in adj.items() for child in value if type(child)==int]
    Quartets        = []
    for _,node in splitting_edges:
        leaves = []
        bfs(node,leaves)
        other_leaves = [leaf for leaf in taxa if leaf not in leaves]
        for pair1 in pair(leaves):
            for pair2 in pair(other_leaves):
                quartet = pair1 + pair2 if pair1[0]<pair2[0] else pair2 + pair1
                Quartets.append(quartet)
    Quartets.sort()
    Unique =[Quartets[0]]
    for i in range(1,len(Quartets)):
        if Quartets[i]!=Unique[-1]:
            Unique.append(Quartets[i])
    return len(Unique),Unique

# cntq is desparately slow for read problems (e.g. n=4525). This quick hack works - ignores actual tree!

def q(n):
    return (n*(n-1)*(n-2)*(n-3))/(4*3*2*1)


if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (cntq(6,
                    '(lobster,(cat,dog),(caterpillar,(elephant,mouse)));'))
        
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = cntq(int(Input[0]), Input[1])
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
                f.write(f'{Result}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
