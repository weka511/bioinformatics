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

class Graph:
    def __init__(self,adj):
        self.adj = adj
    
    def __str__(self):
        return self.bfs_newick()
    
    def bfs_newick(self,node=0):
        children = self.adj[node]
        newick = []
        for child in children:
            if type(child)==int:
                newick.append(self.bfs_newick(node=child))
            else:
                newick.append(child)
        representation = ','.join(newick) 
        return f'({representation})'
    
    def Edges(self):
        for a,b in self.adj.items():
            for c in b:
                yield a,c


def insert(species,edge,graph):
    nextNode = max(list(graph.adj.keys())) + 1
    n1,n2    = edge
    adj      = {}
    if type(n2)==int:
        adj[nextNode] = [species,n2]
        for node,links in graph.adj.items():
            if node==n1:
                adj[node] = [nextNode if ll==n2 else ll for ll in links]
            else:
                adj[node] = links        
    else:
        adj[nextNode] = [species,n2]
        for node,links in graph.adj.items():
            if node==n1:
                adj[node] = [nextNode if ll==n2 else ll for ll in links]
            else:
                adj[node] = links
    return Graph(adj)
    
# EnumerateUnrootedBinaryTrees
#
# Given: A collection of species names representing n taxa.
#
# Return: A list containing all unrooted binary trees whose leaves are these n
#         taxa. Trees should be given in Newick format, with one tree on each line; 
#         the order of the trees is unimportant.
#
# Idea: all rooted trees with a given number of leaves are isomorphic if we
#       ignore the labels of the leaves and nodes. Therfore it is enough to 
#       build a tree with 3 leaves, and keep adding one leaf at a time in all available positions.

def EnumerateUnrootedBinaryTrees(species):
    n              = len(species) 
    def enumerate(n):
        if n==3:
            return [Graph({0:[species[0], species[1], species[2]]})]
        else:
            Enumeration = []
            for graph in enumerate(n-1):
                for edge in graph.Edges():
                    Enumeration.append(insert(species[n-1],edge,graph))
            return Enumeration
        
    return enumerate(n)
    
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        for tree in EnumerateUnrootedBinaryTrees('dog cat mouse elephant unicorn'.split()):
            print (f'{tree};\n')
         
    if args.rosalind:
        Input = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for tree in EnumerateUnrootedBinaryTrees(Input[0].split()):
                print (f'{tree};')
                f.write(f'{tree};\n')
                
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
