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
    
    def get_symbols(self):
        return list(self.Edges.keys())
    
    def bfs(self,
            prefix        = '',
            predecessor   = None,
            visitLeaf     = lambda prefix,node:None,
            visitInternal = lambda prefix,symbol,node,predecessor:None):
        if self.isLeaf():
            visitLeaf(prefix,self)
        else:
            for symbol,edge in self.Edges.items():
                visitInternal(prefix,symbol,self,predecessor)
                edge.EndingNode.bfs(prefix        = prefix + '-',
                                    predecessor   = self,
                                    visitLeaf     = visitLeaf,
                                    visitInternal = visitInternal)
                           
class Edge:
    def __init__(self,EndingNode,Position):
        self.EndingNode = EndingNode
        self.Position   = Position
        
def create_suffix_trie(Text):
    Trie = Node()
    
    for i in range(len(Text)):
        currentNode = Trie
        for j in range(i,len(Text)):
            currentSymbol = Text[j]
            if currentSymbol in currentNode.Edges:
                currentNode = currentNode.Edges[currentSymbol].EndingNode
            else:
                newNode                          = Node()
                newEdge                          = Edge(newNode,j)
                currentNode.Edges[currentSymbol] = newEdge
                currentNode                      = newNode
        if currentNode.isLeaf():
            currentNode.Label = i
        
    return Trie
 
def create_branches(Trie):      
    def visitInternal(prefix,symbol,node,predecessor,Branches):
        if len(node.Edges.items())==1:
            if len(predecessor.Edges)>1:
                #print (f'<<{predecessor.get_symbols()}')
                if len(Branches[-1])>0:
                    Branches.append([])
                Predecessors.append(predecessor)
                #Branches[-1].append(predecessor)
            Branches[-1].append(node)
            #b = ''.join(s for n in Branches[-1] for s in n.get_symbols())
            #print (f'{node.seq} {node.get_symbols()} {b}')
        else:
            if len(Branches[-1])>0:
                Branches.append([])
            
    def visitLeaf(prefix,node,Branches):
        if len(Branches[-1])>0:
            Branches.append([])
    
    Branches      = [[]]
    Predecessors  = []
    Trie.bfs(visitLeaf     = lambda prefix,node:visitLeaf (prefix,node,Branches),
             visitInternal = lambda prefix,symbol,node,predecessor: visitInternal(prefix,symbol,node,predecessor,Branches))
    #print( len(Predecessors),len(Branches))
    #for b in Branches:
        #print (len(b))
    if len(Branches[-1])==0:
        Branches.pop()
    assert len(Predecessors)==len(Branches) 
    return Branches

def create_suffix_tree(Trie):
    pass   
    
 

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        Trie = create_suffix_trie('GTCCGAAGCTCCGG$')
        Trie.bfs(visitLeaf     = lambda prefix,node:print (f'{prefix}{node.Label}'),
                 visitInternal = lambda prefix,symbol,node,predecessor:print (f'{prefix}{symbol}'))
        print ('==========')
        for Branch in create_branches(Trie):
            symbols = [s for node in Branch for s in node.Edges.keys()]
            if len(symbols)>0:
                print ('-'.join(symbols))
        x=0        
        #Tree = create_suffix_tree(Trie)
        
    

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
