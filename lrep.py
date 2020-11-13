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

#  lrep Finding the Longest Multiple Repeat 

import argparse
import os
import time
from   helpers import read_strings

class Node:
    seq = 0
    def __init__(self):
        self.Edges =  {}
        self.Label =  None
        self.seq   =  Node.seq
        Node.seq   += 1
        
    def isLeaf(self):
        return len(self.Edges)==0
    
    def bfs(self,
            prefix        = '',
            visitLeaf     = lambda prefix,node:print (f'{prefix}{node.Label}'),
            visitInternal = lambda prefix,symbol,node:print (f'{prefix}{symbol}')):
        if self.isLeaf():
            visitLeaf(prefix,self)
        else:
            for symbol,edge in self.Edges.items():
                visitInternal(prefix,symbol,self)
                edge.EndingNode.bfs(prefix + '-')
                           
class Edge:
    def __init__(self,EndingNode,Position):
        self.EndingNode = EndingNode
        self.Position   = Position
        
def create_suffix_trie(Text):
    Nodes = []
    Trie = Node()
    Nodes.append(Trie)
    
    for i in range(len(Text)):
        currentNode = Trie
        for j in range(i,len(Text)):
            currentSymbol = Text[j]
            if currentSymbol in currentNode.Edges:
                currentNode = currentNode.Edges[currentSymbol].EndingNode
            else:
                newNode                          = Node()
                Nodes.append(newNode)
                newEdge                          = Edge(newNode,j)
                currentNode.Edges[currentSymbol] = newEdge
                currentNode                      = newNode
        if currentNode.isLeaf():
            currentNode.Label = i
        
    return Trie,Nodes
 
def create_suffix_tree(Text):
    def create_branches(Nodes):
        OpenNodes  = {node.seq: node for node in Nodes if len(node.Edges)==1}
        
        for seq in sorted([seq for seq in OpenNodes.keys()]):
            if seq in OpenNodes:
                node   = OpenNodes.pop(seq)
                Branch = [Node]
                keys   = []
                while len(node.Edges)==1:
                    key    = list(node.Edges.keys())[0]
                    node   = node.Edges[key].EndingNode
                    Branch.append(node)
                    keys.append(key)
                    OpenNodes.pop(node.seq,'')
                
                yield Branch,keys
    
    Trie,Nodes = create_suffix_trie(Text)
    
    for Branch,keys in create_branches(Nodes):
        text      = ''.join(keys)
        firstNode = Branch[0]
        lastNode  = Branch[-1]
        print (text,[key for key in lastNode.Edges.keys()])
        if len(lastNode.Edges)==0:
            pass
    x=0
    
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        Trie,_ = create_suffix_trie('GTCCGAAGCTCCGG$')
        Trie.bfs()
        #create_suffix_tree('GTCCGAAGCTCCGG$')
        
    

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = None
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
