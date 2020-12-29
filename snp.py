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


from rosalind   import RosalindException
from deprecated import deprecated
from numpy      import empty_like,arange

# BA9A Construct a Trie from a Collection of Patterns 
# Trie - a trie, represented as a collection of Nodes and Edges
# 
# create_trie
#    The trie has a single root node with indegree 0, denoted root.
#    Each edge of Trie(Patterns) is labeled with a letter of the alphabet.
#    Edges leading out of a given node have distinct labels.
#    Every string in Patterns is spelled out by concatenating the letters along some path from the root downward.
#    Every path from the root to a leaf, or node with outdegree 0, spells a string from Patterns.
#
# Parameters:
#    Patterns   List of patterns
#    root       Controls whther root is 0, 1, or some other value
#
# Representation:
#      Trie is represented as a dict or dicts
#      {node: {symbol:node}}
def create_trie(Patterns,root=1):
    next_node = root + 1
    Trie      = {root:{}}
    for Pattern in Patterns:
        currentNode = root
        for currentSymbol in Pattern:
            if currentSymbol in Trie[currentNode]:
                currentNode = Trie[currentNode][currentSymbol]
            else:
                new_node                         = next_node
                Trie[new_node]                   = {}
                Trie[currentNode][currentSymbol] = new_node
                currentNode                      = new_node
                next_node                       += 1

    return Trie


                
            
# BA9B Implement TrieMatching


def MatchPrefix(Text,Trie):
    i      = iter(Text)
    symbol = next(i)
    v      = min(Trie.keys())
    path   = [v]
    while True:
        if len(Trie[v])==0:
            return path
        elif symbol in Trie[v]:
            w      = Trie[v][symbol]
            try:
                symbol = next(i)
                v      = w
                path.append(v)
            except StopIteration:
                if len(Trie[w])==0:
                    return path
                else:
                    return []
        else:
            return []

def MatchAll(Text,Trie):
    return [i for i in range(len(Text)) if MatchPrefix(Text[i:],Trie)]




# BA9C Construct the Suffix Tree of a String

# ConstructModifiedSuffixTrie
#
# Construct Modified Suffix Trie, as described in 
# Charging Station Bioinformatics Algorithms Vol II page 165
#
# See also https://sandipanweb.wordpress.com/2017/05/10/suffix-tree-construction-and-the-longest-repeated-substring-problem-in-python/


class SuffixTree:
    next_seq = 0
    class Node:
        def __init__(self):
            self.symbol         = None
            self.edges          = {}
            self.label          = None
            self.seq            = SuffixTree.next_seq
            SuffixTree.next_seq += 1
 
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
                    
        def create_adj(self,adj={}):
            if len(self.edges)==0:
                adj[self.seq] = {}
            else:
                adj[self.seq] = [edge.node.seq for edge in self.edges.values()]
                for edge in self.edges.values():
                    edge.node.create_adj(adj)
 
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
        Leaves   = []
        Internal = []
        Nodes    = {}
        for i in range(len(text)):
            currentNode = self.root
            for j in range(i,len(text)):
                currentSymbol = text[j]
                if currentNode.hasLabelledEdge(currentSymbol):
                    currentNode = currentNode.endEdge(currentSymbol).node
                else:
                    newNode                          = self.Node()
                    Nodes[newNode.seq]               = newNode
                    currentNode.edges[currentSymbol] = self.Edge(newNode,j)
                    currentNode                      = newNode
                    Internal.append(newNode)                    
            if currentNode.isLeaf():
                currentNode.setLabel(i)
                Leaves.append(currentNode)
  
                
        return Leaves,Internal,Nodes
    
    def print(self):
        self.root.print()
        
    def create_adj(self):
        adj = {}
        self.root.create_adj(adj)
        return adj
    
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

def SuffixArray(s,auxiliary=False):
    r = [i for (_,i) in sorted([(s[i:],i) for i in range(len(s))],
                               key=lambda x:x[0])]
    if auxiliary:
        n    = len(s)
        p    = empty_like(r)
        p[r] = arange(len(p), dtype=p.dtype)   # https://stackoverflow.com/questions/9185768/inverting-permutations-in-python
        LCP  = []
        for i in range(len(r)-1):
            i0     = r[i]
            i1     = r[i+1]
            LCP.append(0)
            for j in range(n-max(i0,i1)):
                if s[i0+j] == s[i1+j]:
                    LCP[-1] += 1
                else:
                    break
        return (r,p,LCP)
    else:
        return r
        
        

#  BA9H 	Pattern Matching with the Suffix Array

# MatchOnePatternUsingSuffixArray
#
# Used by MatchPatternsUsingSuffixArray to match one pattern
#
def MatchOnePatternUsingSuffixArray(Text,Pattern,SuffixArray):
    minIndex = 0
    maxIndex = len(Text)
    while minIndex < maxIndex:
        midIndex = (minIndex + maxIndex)//2
        if Pattern>Text[SuffixArray[midIndex]:]:
            minIndex = midIndex + 1
        else:
            maxIndex = midIndex
            
    first    = minIndex
    maxIndex = len(Text)
    while minIndex<maxIndex:
        midIndex = (minIndex + maxIndex)//2
        if Text[SuffixArray[midIndex]:].startswith(Pattern):
            minIndex = midIndex + 1
        else:
            maxIndex = midIndex
            
    last = maxIndex
    return (first,last)

# MatchPatternsUsingSuffixArray
#
# Given: A string Text and a collection of strings Patterns.
#
# Return: All starting positions in Text where a string from Patterns appears as a substring.

def MatchPatternsUsingSuffixArray(Text,Patterns):
    sfa = SuffixArray(Text)
    Matches = []
    for Pattern in Patterns:
        first,last = MatchOnePatternUsingSuffixArray(Text,Pattern,sfa)
        for i in range(first,last):
            Matches.append(sfa[i])
    return sorted(Matches)

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
    ch        = last[i]              # 0 based
    n         = getN(ch,i,last)
    pos,count = get_char(ch,first,n)
    return pos

# BA9L 	Implement BWMatching
#
# Given: A string BWT(Text), followed by a collection of strings Patterns.
#
# Return: A list of integers, where the i-th integer corresponds to the number
#          of substring matches of the i-th member of Patterns in Text.

def BW_Match(LastColumn,Patterns):
    def LastColumnContains(symbol,top,bottom):
        topIndex    = None
        bottomIndex = None
        N           = len(LastColumn)
        
        for i in range(top,N):
            if LastColumn[i]==symbol:
                bottomIndex = i
                break
            
        for i in range(bottom,-1,-1):
            if LastColumn[i]==symbol:
                topIndex = i
                break
            
        return topIndex,bottomIndex
    
    def Match(Pattern):
        top    = 0
        bottom = len(LastColumn)-1
        while top <= bottom:
            if len(Pattern) > 0:
                symbol  = Pattern[-1]
                Pattern = Pattern[:-1]
                topIndex,bottomIndex = LastColumnContains(symbol,top,bottom)
                if type(topIndex)==int and type(bottomIndex)==int:
                    top    = LastToFirst(LastColumn,bottomIndex)
                    bottom = LastToFirst(LastColumn,topIndex)
                else:
                    return 0
            else:
                return bottom - top + 1
        return 0
    
    return [Match(Pattern) for Pattern in Patterns]

#  BA9M 	Implement BetterBWMatching 


def BetterBWMatching(LastColumn,Patterns):
    
    def Count(symbol,i):
        return sum([1 for ch in LastColumn[:i] if ch==symbol]) 

    def FirstOccurence(ch):
        for i in range(len(FirstColumn)):
            if FirstColumn[i]==ch:
                return i
     
    def LastColumnContains(symbol,top,bottom):
        topIndex    = None
        bottomIndex = None
        N           = len(LastColumn)
        
        for i in range(top,N):
            if LastColumn[i]==symbol:
                bottomIndex = i
                break
            
        for i in range(bottom,-1,-1):
            if LastColumn[i]==symbol:
                topIndex = i
                break
            
        return topIndex,bottomIndex
    
    def Match(Pattern):
        top    = 0
        bottom = len(LastColumn) - 1
        while top <= bottom:
            if len(Pattern) > 0:
                symbol  = Pattern[-1]
                Pattern = Pattern[:-1]
                topIndex,bottomIndex = LastColumnContains(symbol,top,bottom)
                if type(topIndex)==int and type(bottomIndex)==int:
                    top    = FirstOccurences[symbol] + Count(symbol,top)
                    bottom = FirstOccurences[symbol] + Count(symbol,bottom+1) - 1
                else:
                    return 0
            else:
                return bottom - top + 1
        return 0    
    
 
    FirstColumn = sorted(LastColumn)
    FirstOccurences = {ch:FirstOccurence(ch) for ch in FirstColumn}
    return [Match(Pattern) for Pattern in Patterns]

#  BA9P 	Implement TreeColoring 
#
#  Colour the internal nodes of a suffix tree given colours of the leaves.
#
#  Given: An adjacency list, followed by colour labels for leaf nodes.
#
#  Return: Colour labels for all nodes.
#
# Implementation: I have coded pure red [True,False], blue [False,True], purple [True,True],
#                 as this makes it easy to generalize to n primary colours

def ColourTree(adj,colours):
    # Make sure that all leaves have been coloured
    for node,children in adj.items():
        if len(children)==0 and node not in colours:
            raise RosalindException(f'Leaf {node} has not been coloured')
        
    n          = len(next(iter(colours.values())))   # Number of primary colours
    Coloured   = {node:colour for node,colour in colours.items()}
    Discovered = [node for node in adj.keys() if node not in Coloured]
    iteration  = 0
    while True:
        Ripe = []
        Rest = []
        for node in Discovered:
            if all(child in Coloured for child in adj[node]):
                Ripe.append(node)
            else:
                Rest.append(node)
        Discovered = Rest
        #Ripe       = [node for node in Discovered if all(child in Coloured for child in adj[node])]
        if iteration%1==0:
            print (f'{iteration}: Discovered={len(Discovered)}, Ripe={len(Ripe)}')
        iteration += 1
        #Discovered = [node for node in Discovered if node not in Ripe]
        if len(Ripe)==0:
            assert len(Discovered)==0,'There are unprocessed elements: inconsistent data?'
            return Coloured

        for node in Ripe:
            Coloured [node] = []
            for i in range(n):
                Coloured [node].append(any([Coloured[child][i] for child in adj [node]]))

#  BA9Q 	Construct the Partial Suffix Array of a String 
#
# Given: A string Text and a positive integer K.
#
# Return: SuffixArrayK(Text), in the form of a list of ordered pairs (i, SuffixArray(i)) 
#                             for all nonempty entries in the partial suffix array.

def PartialSuffixArray(String,K) :
    Suffixes = SuffixArray(String)
    return [(i,Suffixes[i]) for i in range(len(Suffixes)) if Suffixes[i]%K ==0]

### Deprecated code ###
                
# BA9A Construct a Trie from a Collection of Patterns 
@deprecated(reason="Use snp.create_tree instead")
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
            
@deprecated(reason='Use snp.MatchPrefix')
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
@deprecated(reason='Use snp.MatchAll')            
def MatchTries(s,t):
    return [i for i in range(len(s)-1) if PrefixTreeMatching(s[i:],t) ]

if __name__=='__main__':
    print (SuffixArray('banana$',auxiliary=True))