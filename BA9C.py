#  Copyright (C) 2019 Greenweaves Software Limited
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# BA9C Construct the Suffix Tree of a String 

import argparse
from helpers import read_strings
#from snp import ConstructSuffixTreeEdges   

           
class SuffixTree:
    
    class Node:
      
        def __init__(self):
            self.symbol = None
            self.edges  = {}
            self.label  = None
        def hasLabelledEdge(self,symbol):
            return symbol in self.edges
        def endEdge(self,symbol):
            return self.edges[symbol] 
        def isLeaf(self):
            return len(self.edges)==0
        def setLabel(self,j):
            self.label=j

        def traverse(self,path=[]):
            if len(self.edges)==0:
                print (f"{self.label}, {''.join(path)}" )
            else:
                for symbol,edge in self.edges.items():
                    edge.node.traverse(path=path+[symbol])
                
    class Edge:
        def __init__(self,node,position):
            self.node     = node
            self.position = position
        
    def __init__(self):
        self.root = self.Node()
        
    def build(self,text):
        for i in range(len(text)):
            currentNode = self.root
            for j in range(i,len(text)):
                currentSymbol = text[j]
                if currentNode.hasLabelledEdge(currentSymbol):
                    currentNode = currentNode.endEdge(currentSymbol).node
                else:
                    newNode                          = self.Node()
                    currentNode.edges[currentSymbol] = self.Edge(newNode,j)
                    currentNode                      = newNode
            if currentNode.isLeaf():
                currentNode.setLabel(i)
    
    def getEdges(self):
        self.root.traverse()
        x=0
    
    
# Each node has a symbol, position and labels of further edges


def check(Edges,Expected):
    print (f'Expected = {len(Expected)} edges, actual={len(Edges)}')
    mismatches = 0
    for a,b in zip(sorted(Edges),sorted(Expected)):
        if a!=b:
            mismatches+=1
            print (f'Expected {b}, was {a}')
    print(f'{0} mismatches')

def compare_edges(Edges,Expected):
    print (f'Expected = {len(Expected)} edges, actual={len(Edges)}')
    expected = iter(sorted(Expected))
    edges    = iter(sorted(Edges))
    exp      = next(expected)
    ed       = next(edges)
    while exp != '-' and ed !='-':
        if exp<ed:
            print('{0},{1}'.format(exp,'-'))
            exp = next(expected,'-')           
        elif ed<exp:
            print('{0},{1}'.format('-',ed))
            ed = next(edges,'-') 
        else:
            exp = next(expected,'-')
            ed = next(edges,'-')
            
if __name__=='__main__':
    parser = argparse.ArgumentParser('BA9C Construct the Suffix Tree of a String ')
    parser.add_argument('--sample', default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',  default=False, action='store_true', help='process extra dataset')
    args = parser.parse_args()
    if args.sample:
        tree = SuffixTree()
        tree.build('ATAAATG$')
        tree.getEdges()
        
        Edges = tree.Edges
        print('---')
        for edge in Edges:
            print(edge)
        print('---')
        Expected = [
            'AAATG$',
            'G$',
            'T',
            'ATG$',
            'TG$',
            'A',
            'A',
            'AAATG$',
            'G$',
            'T',
            'G$',
            '$'        
        ]
        check(Edges,Expected)
        
    if args.extra:
        Input,Expected  = read_strings('data/SuffixTreeConstruction.txt',init=0)
        Edges           = ConstructSuffixTreeEdges(Input[0])   
        
        compare_edges(Edges,Expected) 
            
            