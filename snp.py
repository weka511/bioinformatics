# Copyright (C) 2019-2020 Greenweaves Software Limited

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

# BA9B Implement TrieMatching



def PrefixTreeMatching(Text,Tree):
    
    def isLeaf(v):
        found = False
        for a,b,c in Tree:
            if a==v:
                return False
            if b==v:
                found = True
        return found
    
    def edge(v,symbol):
        for a,b,c in Tree:
            if a==v and c ==symbol:
                return b
        return -1
    def pattern(v):
        result=[]
        return result
            
    i = 0
    v = 0
    while True:
        if isLeaf(v):
            #print (i,v)
            return True
            #pattern(v)
        else:
            w=edge(v,Text[i])
            if w>-1:
                i+=1
                v=w
            else:
                return False
            
def MatchTries(s,t):
    return [i for i in range(len(s)-1) if PrefixTreeMatching(s[i:],t) ]


# BA9C Construct the Suffix Tree of a String

# ConstructModifiedSuffixTrie
#
# Construct Modified Suffix Trie, as described in 
# Charging Station Bioinformatics Algorithms Vol II page 165
#
# See also https://sandipanweb.wordpress.com/2017/05/10/suffix-tree-construction-and-the-longest-repeated-substring-problem-in-python/


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

        def print(self,path=[]):
            if len(self.edges)==0:
                print (f"{self.label}, {''.join(path)}" )
            else:
                for symbol,edge in self.edges.items():
                    edge.node.print(path=path+[symbol])
 
        def collectEdges(self,path=[],accumulator=[]):
            #print (self.symbol,self.label)
            if len(self.edges)==0:
                accumulator.append(path+[self.symbol])
            elif len(self.edges)==1:
                for symbol,edge in self.edges.items():
                    edge.node.collectEdges(path+[symbol],accumulator=accumulator)                
                #key = next(iter(self.edges)) 
                #self.edges[key].node.collectEdges(path+[self.symbol],accumulator=accumulator)
            else:
                if len(path)>0:
                    accumulator.append(path+[self.symbol])
                for symbol,edge in self.edges.items():
                    edge.node.collectEdges([symbol],accumulator=accumulator)
            return accumulator
        
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
    
    def print(self):
        self.root.print()
        x=0
    
    def collectEdges(self):
        return [ ''.join(run[:-1]) for run in self.root.collectEdges()]

# BA9I Construct the Burrows-Wheeler Transform of a String

def BurrowsWheeler(s):
    def cyclicPermutation(s):
        return s[1:] + s[0:1]

    perms = []
    perm = s
    for i in range(len(s)):
        perm = cyclicPermutation(perm)
        perms.append(perm)
    perms.sort()
    bwt = [perm[-1] for perm in perms]
    return ''.join(bwt)

