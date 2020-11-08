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

#  eubt Enumerating Unrooted Binary Trees

import argparse
import os
import time
from   helpers import read_strings
from   phylogeny import NumberBinaryTrees

def create_adj(edges):
    product = {}
    for a,b in edges:
        if a in product:
            product[a].append(b)
        else:
            product[a] = [b]
    return product

def bfs(adj,node=0):
    children = adj[node]
    newick = []
    for child in children:
        if type(child)==int:
            newick.append(bfs(adj,node=child))
        else:
            newick.append(child)
    representation = ','.join(newick) 
    return f'({representation})'

def newick(edges):
    return f'{bfs(create_adj(edges))};'

    
# EnumerateUnrootedBinaryTrees
#
# Given: A collection of species names representing n taxa.
#
# Return: A list containing all unrooted binary trees whose leaves are these n
#         taxa. Trees should be given in Newick format, with one tree on each line; 
#         the order of the trees is unimportant.

def EnumerateUnrootedBinaryTrees(species):
    def permutations(xs):
        if len(xs)==1:
            return [xs]
        else:
            Result = []
            for x in xs:
                rest = permutations([xx for xx in xs if xx !=x])
                for perm in rest:
                    Result.append([x]+perm)
            return Result
                
    def selector():
        perms = permutations(list(range(n)))
        for perm in perms:
            yield perm
    n              = len(species)
    N              = NumberBinaryTrees(n,rooted=False)
    nInternalNodes = n - 2 
    InternalNodes  = list(range(nInternalNodes))
    sel            = selector()
    for i in range(N):
        edges     = [(i,i+1) for i in range(nInternalNodes-1)]
        selection = next(sel)
        for i in range(nInternalNodes):
            edges.append((i,species[selection[i]]))
        edges.append((0,species[selection[-2]]))
        edges.append((nInternalNodes-1,species[selection[-1]]))
        yield newick(edges)
    

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        for tree in EnumerateUnrootedBinaryTrees('dog cat mouse elephant'.split()):
            print (tree)
         
    if args.rosalind:
        Input = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for tree in EnumerateUnrootedBinaryTrees(Input[0].split()):
                print (tree)
                f.write(f'{tree}\n')
                
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
