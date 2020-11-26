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
# along with this program.  If not, see <http://www.gnu.org/licenses/>

# Code for: Chapter 9, how do we locate disease causing mutations?
# BA9A Construct a Trie from a Collection of Patterns 

# Trie - a trie, represented as a collection of Nodes and Edges
# 
#    The trie has a single root node with indegree 0, denoted root.
#    Each edge of Trie(Patterns) is labeled with a letter of the alphabet.
#    Edges leading out of a given node have distinct labels.
#    Every string in Patterns is spelled out by concatenating the letters along some path from the root downward.
#    Every path from the root to a leaf, or node with outdegree 0, spells a string from Patterns.

class Trie:
    class Node:
        def __init__(self,id):
            self.edges = {}
            self.id    = id
            
        def add(self,edge):
            self.edges[edge.symbol] = edge
            
        def hasEdge(self,symbol):
            return symbol in self.edges
        
        def getEdge(self,symbol):
            return self.edges[symbol]        
        
        def Edges(self):
            for symbol,edge in self.edges.items():
                yield self.id, edge.end.id, edge.symbol
                yield from edge.end.Edges()
                
        def isLeaf(self):        
            return len(self.edges)==0
        
    class Edge:
        def __init__(self,end=None,symbol=None):
            self.end    = end
            self.symbol = symbol
            
    def __init__(self,Patterns):
        self.nodeCounter = 0
        self.root        = self.Node(self.nodeCounter)
        self.nodeCounter +=1
        
        for Pattern in Patterns: # Build Trie - http://rosalind.info/problems/ba9a/
            currentNode = self.root
            for currentSymbol in Pattern:
                if currentNode.hasEdge(currentSymbol):
                    currentNode = currentNode.getEdge(currentSymbol).end
                else:
                    newNode = self.Node(self.nodeCounter)
                    self.nodeCounter+=1
                    currentNode.add(self.Edge(end=newNode,symbol=currentSymbol))
                    currentNode = newNode

    def Edges(self):
        return self.root.Edges()
    
    def MatchAll(self,text,minimum_length=1):
        return [i for i in range(len(text)) if self.Match(text[i:]) and len(text[i:])>minimum_length]
    
    def Match(self,text):
        def characters():
            i = 0
            while i<len(text):
                yield text[i]
                i += 1
                
        symbols = characters()
        symbol = next(symbols)
        v      = self.root
        path   = []
        while True:
            if v.isLeaf():
                return path
            elif v.hasEdge(symbol):
                edge   = v.getEdge(symbol)
                path.append(edge.symbol)
                v      = edge.end
                try:
                    symbol = next(symbols)
                except StopIteration:
                    return path
            else:
                return []
                
            
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
            if len(self.edges)==0:
                accumulator.append(path+[self.symbol])
            elif len(self.edges)==1:
                for symbol,edge in self.edges.items():
                    edge.node.collectEdges(path+[symbol],accumulator=accumulator)                
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
        return [ ''.join(run[:-1]) for run in self.root.collectEdges(accumulator=[])]

# BA9D Find the Longest Repeat in a String
# BA9E Find the Longest Substring Shared by Two Strings

def sort_by_descending_length(edges):
    return [e for _,e in sorted([(len(edge),edge) for edge in edges],reverse=True)]

def FindLongestRepeat(string1,string2=None):        
    tree = SuffixTree()
    string = string1 if string2 == None else string1 + string2
    tree.build(string)
    edges = sort_by_descending_length(tree.collectEdges())
    for edge in edges:
        index1 = string.find(edge,0,-1 if string2==None else len(string1))
        if index1==-1: continue
        index2 = string.find(edge,(index1+1) if string2==None else len(string1))
        if index2==-1: continue
        # try to extend edge -- first to the left
        i = 1
        while string[index1-i]==string[index2-i]:
            i+=1
        i -= 1
        # then to the right
        j  = len(edge)
        while index2+j<len(string) and string[index1+j]==string[index2+j]:
            j+=1

        return string[index1-i:index1+j]

# BA9G Construct the Suffix Array of a String

# Suffix Array
#
# In 1993, Udi Manber and Gene Myers introduced suffix arrays as a memory-efficient alternative to suffix trees. 
# To construct SuffixArray(Text), we first sort all suffixes of Text lexicographically, assuming that "$" 
# $ comes first in the alphabet. The suffix array is the list of starting positions of these sorted suffixes.
#
# I have used a naive algorithm, as it appears to be adequate for the test data

def SuffixArray(s):
    suffixes = [(s[i:],i) for i in range(len(s))]
    suffixes.sort()
    return [i for (_,i) in suffixes]
    
# BA9I Construct the Burrows-Wheeler Transform of a String
#
# Given a string Text, form all possible cyclic rotations of Text; a cyclic
# rotation is defined by chopping off a suffix from the end of Text and appending
# this suffix to the beginning of Text. Next - similarly to suffix arrays - order
# all the cyclic rotations of Text lexicographically to form a |Text| x |Text| matrix of symbols
# that we call the Burrows-Wheeler matrix and denote by M(Text).
#
# Note that the first column of M(Text) contains the symbols of Text ordered lexicographically.
# In turn, the second column of M(Text) contains the second symbols of all cyclic rotations of Text,
# and so it too represents a (different) rearrangement of symbols from Text. 
# The same reasoning applies to show that any column of M(Text) is some rearrangement 
# of the symbols of Text. We are interested in the last column of M(Text), 
# called the Burrows-Wheeler transform of Text, or BWT(Text).
def BurrowsWheeler(s):
    def cyclicPermutation(s):
        return s[1:] + s[0:1]

    perms = []
    perm  = s
    for i in range(len(s)):
        perm = cyclicPermutation(perm)
        perms.append(perm)
    perms.sort()
    bwt = [perm[-1] for perm in perms]
    return ''.join(bwt)

# BA9J Reconstruct a String from its Burrows-Wheeler Transform

def getN(ch,row,column):
    n = 0
    i = 0
    while i<=row:
        if ch==column[i]:
            n+=1
        i+=1
    return n

def get_char(ch,column,seq):
    pos   = 0
    count = 0
    for i in range(len(column)):
        if column[i]==ch:
            pos = i
            count+=1
            if count==seq:
                return pos,count
            
def InverseBWT(string):
    lastColumn  = [a for a in string]
    firstColumn = sorted(lastColumn)
    Result      = [firstColumn[0]]
    ch          = min(string)
    seq         = 1
    while len(Result)<len(string):
        row,count    = get_char(ch,lastColumn,seq)
        ch           = firstColumn[row]
        seq          = getN(ch,row,firstColumn)
        Result.append(ch)
 
    return ''.join(Result[1:]+Result[0:1])

# BA9K Generate the Last-to-First Mapping of a String

def LastToFirst(Transform,i):
    last      = Transform
    first     = sorted(Transform)
    ch        = last[i-1]
    n         = getN(ch,i,last)
    pos,count = get_char(ch,first,n)
    return pos