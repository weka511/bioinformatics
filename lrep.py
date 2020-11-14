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

#class Node:
    #seq = 0
    #def __init__(self):
        #self.Edges =  {}
        #self.Label =  None
        #self.seq   =  Node.seq
        #Node.seq   += 1
        
    #def isLeaf(self):
        #return len(self.Edges)==0
    
    #def get_symbols(self):
        #return list(self.Edges.keys())
    
    #def bfs(self,
            #prefix        = '',
            #predecessor   = None,
            #visitLeaf     = lambda prefix,node:None,
            #visitInternal = lambda prefix,symbol,node,predecessor:None):
        #if self.isLeaf():
            #visitLeaf(prefix,self)
        #else:
            #for symbol,edge in self.Edges.items():
                #visitInternal(prefix,symbol,self,predecessor)
                #edge.EndingNode.bfs(prefix        = prefix + '-',
                                    #predecessor   = self,
                                    #visitLeaf     = visitLeaf,
                                    #visitInternal = visitInternal)
                           
#class Edge:
    #def __init__(self,EndingNode,Position):
        #self.EndingNode = EndingNode
        #self.Position   = Position
        
#def create_suffix_trie(Text):
    #Trie = Node()
    
    #for i in range(len(Text)):
        #currentNode = Trie
        #for j in range(i,len(Text)):
            #currentSymbol = Text[j]
            #if currentSymbol in currentNode.Edges:
                #currentNode = currentNode.Edges[currentSymbol].EndingNode
            #else:
                #newNode                          = Node()
                #newEdge                          = Edge(newNode,j)
                #currentNode.Edges[currentSymbol] = newEdge
                #currentNode                      = newNode
        #if currentNode.isLeaf():
            #currentNode.Label = i
        
    #return Trie
 
#def create_branches(Trie):
    
    #def visitInternal(prefix,symbol,node,predecessor,Branches):
        #if len(node.Edges.items())==1:
            #if len(predecessor.Edges)>1:
                #if len(Branches[-1])>0:
                    #Branches.append([])
                #Predecessors.append((predecessor,symbol))
            #Branches[-1].append(node)

        #else:
            #if len(Branches[-1])>0:
                #Branches.append([])
            
    #def visitLeaf(prefix,node,Branches):
        #if len(Branches[-1])>0:
            #Branches.append([])
    
    #Branches      = [[]]
    #Predecessors  = []
    #Trie.bfs(visitLeaf     = lambda prefix,node:visitLeaf (prefix,node,Branches),
             #visitInternal = lambda prefix,symbol,node,predecessor: visitInternal(prefix,symbol,node,predecessor,Branches))

    #if len(Branches[-1])==0:
        #Branches.pop()
    #assert len(Predecessors)==len(Branches) 
    #return Branches,Predecessors

#def create_suffix_tree(Trie):
    #pass   

class Node:
    def __init__(self,seq):
        self.seq        = seq
        self.edges      = []
    def bfs(self,
            Nodes,
            text,
            prefix        = '',
            visitLeaf     = lambda prefix:None,
            visitInternal = lambda prefix,symbol:None):
        if len(self.edges)==0:
            visitLeaf(prefix)
        else:
            for next_pos,s,l in self.edges:
                visitInternal(prefix,text[s:s+l])
                Nodes[next_pos].bfs(Nodes,
                                    text,
                                    prefix        = prefix + '-',
                                    visitLeaf     = visitLeaf,
                                    visitInternal = visitInternal)
        
def lrep(text,k,edges):
    entries    = sorted([[int(e.replace('node','')) for e in edge.split()] for edge in edges],key=lambda x:-x[-1])
    for i in range(len(entries)):
        _,_,pos_i,length_i = entries[i]
        if pos_i+length_i>len(text):
            length_i-=1
        for j in range(i,len(entries)):
            _,_,pos_j,length_j = entries[j]
            if pos_j==pos_i: continue
            if pos_j+length_j>len(text):
                length_j-=1        
            if text[pos_i-1:pos_i-1+length_i]==text[pos_j-1:pos_j-1+length_j]:
                while pos_i>1 and pos_j>1 and text[pos_i-2]==text[pos_j-2]:
                    pos_i    -= 1
                    pos_j    -= 1
                    length_i += 1
                return text[pos_i-1:pos_i-1+length_i]
    
    #def create_suffix_tree():
        #Nodes = []
        #for edge in edges:
            #entry = [int(e.replace('node','')) for e in edge.split()]
            #while len(Nodes)<max(entry[0],entry[1]):
                #Nodes.append(None)
            #if Nodes[entry[0]-1] ==None:
                #Nodes[entry[0]-1] = Node(entry[0]-1)
            #if Nodes[entry[1]-1] ==None:
                #Nodes[entry[1]-1] = Node(entry[1]-1)
            #Nodes[entry[0]-1].edges.append([entry[1]-1,entry[2]-1,entry[3]])
        #return Nodes
    #Nodes = create_suffix_tree()
    #Nodes[0].bfs(Nodes,s,
                 #visitInternal = lambda prefix,symbol:print (f'{prefix}{symbol}'))
    x=0

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (lrep('CATACATAC$',
                    2,
                    ['node1 node2 1 1',
                     'node1 node7 2 1',
                     'node1 node14 3 3',
                     'node1 node17 10 1',
                     'node2 node3 2 4',
                     'node2 node6 10 1',
                     'node3 node4 6 5',
                     'node3 node5 10 1',
                     'node7 node8 3 3',
                     'node7 node11 5 1',
                     'node8 node9 6 5',
                     'node8 node10 10 1',
                     'node11 node12 6 5',
                     'node11 node13 10 1',
                     'node14 node15 6 5',
                     'node14 node16 10 1']))
        #Trie = create_suffix_trie('GTCCGAAGCTCCGG$')
        #Trie.bfs(visitLeaf     = lambda prefix,node:print (f'{prefix}{node.Label}'),
                 #visitInternal = lambda prefix,symbol,node,predecessor:print (f'{prefix}{symbol}'))
        #print ('==========')
        #Branches,Predecessors = create_branches(Trie)
        #for Branch,pre in zip(Branches,Predecessors):
            #_,symbol = pre
            #symbols = [s for node in Branch for s in node.Edges.keys()]
            #if len(symbols)>0:
                #print (f'{symbol}-{"".join(symbols)}')
        #x=0        
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
