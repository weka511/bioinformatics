#!/usr/bin/env python

#   Copyright (C) 2020-2024 Greenweaves Software Limited

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

'''Code for: Chapter 9, how do we locate disease causing mutations?'''

from unittest import TestCase, main, skip
from deprecated import deprecated
import numpy as np
from rosalind import RosalindException

def create_trie(Patterns,root=1):
    '''
    BA9A Construct a Trie from a Collection of Patterns
    Trie - a trie, represented as a collection of Nodes and Edges

    create_trie
       The trie has a single root node with indegree 0, denoted root.
       Each edge of Trie(Patterns) is labeled with a letter of the alphabet.
       Edges leading out of a given node have distinct labels.
       Every string in Patterns is spelled out by concatenating the letters along some path from the root downward.
       Every path from the root to a leaf, or node with outdegree 0, spells a string from Patterns.

    Parameters:
       Patterns   List of patterns
       root       Controls whther root is 0, 1, or some other value

    Representation:
         Trie is represented as a dict or dicts
         {node: {symbol:node}}
    '''
    next_node = root + 1
    Trie = {root:{}}
    for Pattern in Patterns:
        currentNode = root
        for currentSymbol in Pattern:
            if currentSymbol in Trie[currentNode]:
                currentNode = Trie[currentNode][currentSymbol]
            else:
                new_node = next_node
                Trie[new_node] = {}
                Trie[currentNode][currentSymbol] = new_node
                currentNode = new_node
                next_node += 1

    return Trie


def MatchPrefix(Text,Trie):
    '''
    BA9B Implement TrieMatching

    MatchPrefix

    Given: A string Text and a collection of strings Patterns.

    Return: First starting positions in Text where a string from Patterns appears as a substring.
    '''
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

# MatchAll
#
# Given: A string Text and a collection of strings Patterns.
#
# Return: All starting positions in Text where a string from Patterns appears as a substring.

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

def sort_by_descending_length(edges):
    return [e for _,e in sorted([(len(edge),edge) for edge in edges],reverse=True)]

# FindLongestRepeat
#
# Given: A string Text.
#
# Return: A longest substring of Text that appears in Text more than once.

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

# BA9E Find the Longest Substring Shared by Two Strings

# LongestSharedSubstring
#
# Find the Longest Substring Shared by Two Strings
#
# Given: Strings Text1 and Text2.
#
# Return: The longest substring that occurs in both Text1 and Text2.
#
# See https://en.wikipedia.org/wiki/Longest_common_substring_problem

def LongestSharedSubstring(s,t):
    # find_straddling
    #
    # Find indices of LCPs that straddle end of s,
    # i.e. one string is from s, the other from t

    def find_straddling():
        previous_from_s = False
        Pairs           = []
        for i in range(len(lcp)):
            from_s          = r[i]<len(s)+2
            if i>0 and previous_from_s != from_s:
                Pairs.append(i)
            previous_from_s = from_s
        return Pairs

    text           = s + '$' + t + '#'
    r,_,lcp        = SuffixArray(text,auxiliary=True,padLCP=True)
    candidate_LCPs = [lcp[i] if i in  find_straddling() else 0 for i in range(len(lcp))]
    index          = np.argmax(candidate_LCPs)
    return text[r[index]:r[index]+candidate_LCPs[index]]


def SuffixArray(s,auxiliary=False,padLCP=False):
    '''
    BA9G Construct the Suffix Array of a String

    Suffix Array

    In 1993, Udi Manber and Gene Myers introduced suffix arrays as a memory-efficient alternative to suffix trees.
    To construct SuffixArray(Text), we first sort all suffixes of Text lexicographically, assuming that "$"
    comes first in the alphabet. The suffix array is the list of starting positions of these sorted suffixes.

    I have used a naive algorithm, as it appears to be adequate for the test data

    This can also create the auxiliary arrays used by Efficient repeat finding via suffix arrays
    Veronica Becher, Alejandro Deymonnaz, and Pablo Ariel Heiber
    '''
    r = [i for (_,i) in sorted([(s[i:],i) for i in range(len(s))],
                               key=lambda x:x[0])]
    if auxiliary:
        n    = len(s)
        p    = np.empty_like(r)
        p[r] = np.arange(len(p), dtype=p.dtype)   # https://stackoverflow.com/questions/9185768/inverting-permutations-in-python
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
        if padLCP:
            LCP.insert(0,0)
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

    return ''.join([row[-1] for row in sorted([s[k:] + s[0:k] for k in range(len(s))])])

# BA9J Reconstruct a String from its Burrows-Wheeler Transform
# BA9K Generate the Last-to-First Mapping of a String

# getN
#
# Used by InverseBWT and LastToFirst: given a character that is known to appear in a
# specific position in a column, deterine whether is is the first, 2nd, ..., occurrence
#
# Parameters:
#     ch
#     row
#     column
def getN(ch,row,column):
    n = 0
    i = 0
    while i<=row:
        if ch==column[i]:
            n+=1
        i+=1
    return n

# get_char
#
# Used by InverseBWT and LastToFirst: find a specific occurrence of chracter in column
#
# Parameters:
#     ch      Character to search for
#     column  Column to be searched
#     seq     Number of occurrence - 1 for first occurence, 2 for 2nd, etc.

def get_char(ch,column,seq):
    count = 0
    for pos in range(len(column)):
        if column[pos]==ch:
            count +=1
            if count==seq:
                return pos

# Inverse BWT
#
# Reconstruct a String from its Burrows-Wheeler Transform
#
# See pages 139-149 of Compeau Pevzner, Volume II

def InverseBWT(BWT):
    lastColumn  = [a for a in BWT]
    firstColumn = sorted(lastColumn)
    Result      = [firstColumn[0]]
    ch          = min(BWT)        # $ in examples in textbook
    seq         = 1               # Keep track of whether this is 1st, 2nd,... occurence
                                  # Of course there is only one occurrence of $, but this
                                  # isn't necessarly true for other characters
    while len(Result) < len(BWT):
        row = get_char(ch,lastColumn,seq)  # Find row number for seq-th occurence of character in last coumn
        ch  = firstColumn[row]             # Find coresponding symbol in first column
                                           # this is the next character arising from a cyclic permutation
        seq = getN(ch,row,firstColumn)     # Indicates whether 1st, 2nd, ... occurrence
        Result.append(ch)

    return ''.join(Result[1:]+Result[0:1])

# BA9K Generate the Last-to-First Mapping of a String
#
#  LastToFirst
#
#  Parameters:
#      BWT
#      i
#
#  Return:
#    The position LastToFirst(i) in FirstColumn in the Burrows-Wheeler matrix if LastColumn = BWT.

def LastToFirst(BWT,i):
    first = sorted(BWT)         # First column = last column, sorted
    ch    = BWT[i]              # ith character, 0 based
    n     = getN(ch,i,BWT)      # character is nth occurrence
    return get_char(ch,first,n) # nth occurence in first column

# ColumnContains
#
# Used by Match to search the last column for last symbol in Pattern
#
# Parameters:
#     symbol   Current last symbol in Pattern
#     top      The first row in last column that matches current suffix (excluding symbol)
#     bottom   The last row in last column that matches current suffix (excluding symbol)
#
# Returns:
#     topIndex     Index in first column of first row that matches current suffix (including symbol)
#     bottomIndex  Index in first column of last row that matches current suffix (including symbol)

def ColumnContains(Column,symbol,top,bottom):
    topIndex    = None
    bottomIndex = None

    for i in range(top,len(Column)):
        if Column[i]==symbol:
            bottomIndex = i
            break

    for i in range(bottom,-1,-1):
        if Column[i]==symbol:
            topIndex = i
            break

    return topIndex,bottomIndex

# getCount
#
# Count number of occurrences of symbol in first i positions of last column.
def getCount(Column, symbol, i):
    return sum([1 for ch in Column[:i] if ch==symbol])

# getFirstOccurrence
#
# Find first occurrence of symbol in first column

def getFirstOccurrence(Column,symbol):
    for i in range(len(Column)):
        if Column[i]==symbol:
            return i

# BA9L 	Implement BWMatching
#
# Given: A string BWT(Text), followed by a collection of strings Patterns.
#
# Return: A list of integers, where the i-th integer corresponds to the number
#          of substring matches of the i-th member of Patterns in Text.

def BW_Match(LastColumn,Patterns):

    # Match
    #
    # Used by BW_Match to match one occurrence of Pattern in text
    def Match(Pattern):
        top    = 0                  # When we iterate, top will be the first row in last column that matches current suffix
        bottom = len(LastColumn)-1  # bottom will be the last row in last column that matches current suffix
        i      = len(Pattern)       # Current position within pattern
        while top <= bottom:
            if i >0 :
                i                   -= 1
                topIndex,bottomIndex = ColumnContains(LastColumn,Pattern[i],top,bottom)
                if type(topIndex)==int and type(bottomIndex)==int:
                    top    = LastToFirst(LastColumn, bottomIndex)
                    bottom = LastToFirst(LastColumn, topIndex)
                else:
                    return 0
            else:
                return bottom - top + 1
        return 0

    return [Match(Pattern) for Pattern in Patterns]

#  BA9M 	Implement BetterBWMatching
#
# Given: A string BWT(Text), followed by a collection of strings Patterns.
#
# Return: A list of integers, where the i-th integer corresponds to the number
#          of substring matches of the i-th member of Patterns in Text.

def BetterBWMatching(LastColumn,Patterns):

    # Match
    #
    # Used by BetterBWMatching to match one occurrence of Pattern in text
    #
    # Parameters:
    #      Pattern
    #      FirstOccurrences
    def Match(Pattern,FirstOccurrences):
        top    = 0
        bottom = len(LastColumn) - 1
        while top <= bottom:
            if len(Pattern) > 0:
                symbol  = Pattern[-1]
                Pattern = Pattern[:-1]
                topIndex,bottomIndex = ColumnContains(LastColumn,symbol,top,bottom)
                if type(topIndex)==int and type(bottomIndex)==int:
                    top    = FirstOccurrences[symbol] + getCount(LastColumn,symbol,top)
                    bottom = FirstOccurrences[symbol] + getCount(LastColumn,symbol,bottom+1) - 1
                else:
                    return 0
            else:
                return bottom - top + 1

        return 0


    FirstColumn      = sorted(LastColumn)
    FirstOccurrences = {ch:getFirstOccurrence(FirstColumn,ch) for ch in FirstColumn}
    return [Match(Pattern,FirstOccurrences) for Pattern in Patterns]



def FindApproximateMatches(Text,Patterns,d):
    '''
    BA9O Find all approximate occurrences of a collection of patterns in a text.
    Given: A string Text, a collection of strings Patterns, and an integer d.

    Return: All positions in Text where a string from Patterns appears as a substring with at most d mismatches.
    '''
    for p in Patterns:
        pass


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
        # if iteration%1==0:
            # print (f'{iteration}: Discovered={len(Discovered)}, Ripe={len(Ripe)}')
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

# mrep

#  MREP Identifying Maximal Repeats
#
# Based on Efficient repeat finding via suffix arrays
# Ver�nica Becher, Alejandro Deymonnaz, and Pablo Ariel Heiber

def mrep(w,ml=20):
    # create_LCP_already_seen
    #
    # Create list of LCP values that have already been seen

    def create_LCP_already_seen():
        Product = [u for u in range(len(LCP)) if LCP[u]<ml]
        Product.append(-1)
        Product.append(n-1)
        return Product

    n       = len(w)
    r,p,LCP = SuffixArray(w,auxiliary=True)
    S       = create_LCP_already_seen()
    I       = np.np.argsort(LCP)

    for t in range(min([t for t in range(len(I)) if LCP[I[t]]>=ml]),
                   n-1):
        i   = I[t]
        p_i = max([j for j in S if j<i]) + 1
        n_i = min([j for j in S if j>i])
        S.append(i)
        if (p_i==0 or LCP[p_i-1]!=LCP[i]) and (n_i==n-1 or LCP[n_i]!=LCP[i]):                          # Maximal on right?
            if r[p_i]==0 or r[n_i]==0 or w[r[p_i]-1]!=w[r[n_i]-1] or p[r[n_i]-1]-p[r[p_i]-1]!=n_i-p_i: # Maximal on left?
                yield w[r[i]:r[i]+LCP[i]]                                                              # banzai


# BA9N 	Find All Occurrences of a Collection of Patterns in a String

def EvenBetterBWMatching(Text,Patterns,K=10):
    # get_entry
    #
    # Search partial suffix array for a value whose position matches row

    def get_entry(row):
        for i,value in PSA:
            if i==row: return value
            if i>row: return None

    # find_position
    #
    # Locate match within Text
    #
    def find_position(row):
        steps = 0
        pos = get_entry(row)
        while pos==None:
            predecessor = LastColumn[row]                              # Predecessor of 1st character of match (cyclic)
            occurrence  = getN(predecessor,row,LastColumn)             # 1st, 2nd, ... occurence
            row         = get_char(predecessor,FirstColumn,occurrence) # Find occurrence in 1st column
            steps += 1
            pos = get_entry(row)
        return steps + pos

    # Match
    #
    # Used by BetterBWMatching to match one occurrence of Pattern in text
    #
    # Parameters:
    #      Pattern
    #      FirstOccurrences
    def Match(Pattern,FirstOccurrences):
        top    = 0
        bottom = len(LastColumn) - 1
        while top <= bottom:
            if len(Pattern) > 0:
                symbol  = Pattern[-1]
                Pattern = Pattern[:-1]
                topIndex,bottomIndex = ColumnContains(LastColumn,symbol,top,bottom)
                if type(topIndex)==int and type(bottomIndex)==int:
                    top    = FirstOccurrences[symbol] + getCount(LastColumn,symbol,top)
                    bottom = FirstOccurrences[symbol] + getCount(LastColumn,symbol,bottom+1) - 1
                else:
                    return []
            else:
                return [find_position(i) for i in range(top,bottom+1)]

        return []

    PSA              = PartialSuffixArray(Text+'$',K)
    LastColumn       = BurrowsWheeler(Text+'$')
    FirstColumn      = sorted(LastColumn)
    FirstOccurrences = {ch:getFirstOccurrence(FirstColumn,ch) for ch in FirstColumn}
    return [match for Pattern in Patterns for match in Match(Pattern,FirstOccurrences)]

if __name__=='__main__':

    class TestSNP(TestCase):
        '''Test cases for: Chapter 9, how do we locate disease causing mutations?'''

        def test_ba9a(self):
            ''' BA9A Construct a Trie from a Collection of Patterns'''
            Trie = create_trie(['ATAGA','ATC','GAT'],root=0)
            self.assertEqual(10,len(list(Trie.keys())))
            trie0 = Trie[0]
            self.assertEqual(1,trie0['A'])
            self.assertEqual(7,trie0['G'])
            trie1 = Trie[1]
            self.assertEqual(2,trie1['T'])
            trie2 = Trie[2]
            self.assertEqual(3,trie2['A'])
            self.assertEqual(6,trie2['C'])
            trie3 = Trie[3]
            self.assertEqual(4,trie3['G'])
            trie4 = Trie[4]
            self.assertEqual(5,trie4['A'])
            self.assertEqual(0,len(Trie[5]))
            self.assertEqual(0,len(Trie[6]))
            trie7 = Trie[7]
            self.assertEqual(8,trie7['A'])
            trie8 = Trie[8]
            self.assertEqual(9,trie8['T'])
            self.assertEqual(0,len(Trie[9]))

        def test_ba9b(self):
            '''BA9B Implement TrieMatching'''
            Trie = create_trie(['ATCG','GGGT'])
            self.assertEqual([1,4,11,15],MatchAll('AATCGGGTTCAATCGGGGT',Trie))


        def test_ba9c(self):
            '''BA9C Construct the Suffix Tree of a String'''
            tree = SuffixTree()
            tree.build('ATAAATG$')
            Edges = tree.collectEdges()
            self.assertEqual(12,len(Edges))
            self.assertIn('AAATG$',Edges)
            self.assertIn('G$',Edges)
            self.assertIn('T',Edges)
            self.assertIn('ATG$',Edges)
            self.assertIn('TG$',Edges)
            self.assertIn('A',Edges)
            self.assertIn('A',Edges)
            self.assertIn('AAATG$',Edges)
            self.assertIn('G$',Edges)
            self.assertIn('T',Edges)
            self.assertIn('G$',Edges)
            self.assertIn('$',Edges)

        def test_ba9d(self):
            '''BA9D Find the Longest Repeat in a String'''
            self.assertEqual('TATCGTT', FindLongestRepeat('ATATCGTTTTATCGTT'))

        def test_ba9e(self):
            '''BA9E Find the Longest Substring Shared by Two Strings'''
            self.assertEqual('AGA',LongestSharedSubstring('TCGGTAGATTGCGCCCACTC','AGGGGCTCGCAGTGTAAGAA'))

        @skip('TBP')
        def test_ba9f(self):
            ''''''
        @skip('TBP')
        def test_ba9g(self):
            ''''''
        @skip('TBP')
        def test_ba9h(self):
            ''''''
        @skip('TBP')
        def test_ba9i(self):
            ''''''
        @skip('TBP')
        def test_ba9j(self):
            ''''''
        @skip('TBP')
        def test_ba9k(self):
            ''''''
        @skip('TBP')
        def test_ba9l(self):
            ''''''

        def test_ba9m(self):
            '''BA9M Implement BetterBWMatching'''
            self.assertEqual([1,2,1],
                             BetterBWMatching('GGCGCCGC$TAGTCACACACGCCGTA',
                                              ['ACC', 'CCG', 'CAG']))

        def test_ba9n(self):
            '''BA9N Find All Occurrences of a Collection of Patterns in a String'''
            self.assertEqual([11, 1, 15, 4],
                             EvenBetterBWMatching('AATCGGGTTCAATCGGGGT', ['ATCG', 'GGGT'], K=5))

        @skip('TBP')
        def test_ba9o(self):
            '''BA9O Find All Approximate Occurrences of a Collection of Patterns in a String'''
            self.assertEqual([2, 4, 4, 6, 7, 8, 9],
                             FindApproximateMatches('ACATGCTACTTT', ['ATT', 'GCC', 'GCTA', 'TATT'], 1))


        @skip('TBP')
        def test_ba9q(self):
            ''''''

        def test_ba9p(self):
            ''' BA9P 	Implement TreeColoring'''
            def colour2list(colour):
                if colour=='red': return [True,False]
                if colour=='blue': return [False,True]
            Colours = {node: colour2list(c) for (node),c in  {
                                  0: 'red',
                                  1: 'red',
                                  3: 'blue',
                                  4: 'blue',
                                  6: 'red'}.items()}
            Coloured = ColourTree({
                                  0 : [],
                                  1 : [],
                                  2 : [0,1],
                                  3 : [],
                                  4 : [],
                                  5 : [3,2],
                                  6 : [],
                                  7 : [4,5,6]},Colours)


        @skip('TBP')
        def test_ba9r(self):
            '''BA9R Construct a Suffix Tree from a Suffix Array'''
            Tree = SuffixArray2Tree('GTAGT$',
                             [5, 2, 3, 0, 4, 1],
                             [0, 0, 0, 2, 0, 1])
            z=0
    main()
