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

#  SUFF Creating a Suffix Tree

import argparse
import os
import time
from   helpers import read_strings

class Node:
    def __init__(self):
        self.edges     = {}
        self.positions = {}
        self.label     = None
        
    def bfs(self,
            prefix        = '',
            visitLeaf     = lambda label,prefix:print (f'{prefix} {label}'),
            visitInternal = lambda node,symbol,position,prefix: print (f'{prefix} {symbol} {position}')):
        if len(self.edges)==0:
            visitLeaf(self.label,prefix+'-')
        else:
            for symbol,node in self.edges.items():
                visitInternal(node,symbol,self.positions[symbol],prefix)
                node.bfs(prefix+'-')
                
def suff(string):
    def createTrie():
        Trie = Node()
        for i in range(len(string)):
            currentNode = Trie
            for j in range(i,len(string)):
                currentSymbol = string[j]
                if currentSymbol in currentNode.edges:
                    currentNode = currentNode.edges[currentSymbol]
                else:
                    newNode                              = Node()
                    currentNode.edges[currentSymbol]     = newNode
                    currentNode.positions[currentSymbol] = j
                    currentNode                          = newNode
            if len(currentNode.edges)==0:
                currentNode.label = i
        return Trie
    
    Trie = createTrie()
    Trie.bfs()
    
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        for substring in suff('ATAAATG$'):
            print (substring)
        
    

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
