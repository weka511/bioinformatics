# Copyright (C) 2019 Greenweaves Software Limited

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

# BA9C Construct the Suffix Tree of a String 

# ConstructSuffixTree
#
# Construct the suffix tree of a string.
#
# Inputs: A string Text.
#
# Return: The strings labeling the edges of SuffixTree(Text).

from rosalind import trie

def ConstructModifiedSuffixTrie(Text):
    def edgeLabelled(currentNode, currentSymbol):
        if currentNode in Trie:
            for symbol,pos,node in Trie[currentNode]:
                if symbol==currentSymbol:
                    return node
        return None
        
    Trie   = {0:[]}
    labels = {}
    n      = 1
    for i in range(len(Text)):
        currentNode = 0
        for j in range(i,len(Text)):
            currentSymbol = Text[j]
            next_node = edgeLabelled(currentNode, currentSymbol)
            if next_node != None:
                currentNode = next_node
            else:
                Trie[n] = []
                Trie[currentNode].append((currentSymbol,j,n))
                currentNode = n
                n+=1
        if len(Trie[currentNode])==0:
            labels[currentNode] = i
    return Trie,labels

def ConstructSuffixTree(Text):
    Trie,labels = ConstructModifiedSuffixTrie(Text)
    return Trie

if __name__=='__main__':
    print (ConstructSuffixTree('ATAAATG$'))
    
    trie,labels = ConstructModifiedSuffixTrie('panamabananas$')
    print (trie)
    print (labels)