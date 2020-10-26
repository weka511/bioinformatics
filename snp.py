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


def ConstructSuffixTreeEdges(s):
    def explore(suffixes):
        if len(suffixes)==0: return
        prefixes = list(set([s[0] for s in suffixes]))
        if len(prefixes)==len(suffixes):
            for s in suffixes:
                Edges.append(s)
        else:
            if len(prefixes)==1:
                for k in range(2,min([len(s) for s in suffixes])):
                    extended_prefixes = list(set([s[0:k] for s in suffixes]))
                    if len(extended_prefixes)==1:
                        prefixes = extended_prefixes
                    else:
                        break
                    
            for p in prefixes:
                subset = [s[len(prefixes[0]):] for s in suffixes if s[0:len(prefixes[0])]==p ]
                if len(subset)>1:
                    Edges.append(p)
                    explore(subset)
                elif len(subset)==1:
                    Edges.append(p+subset[0])
                #else:
                    #Edges.append(p)
                    
     
    next_node = 0                
    Edges     = []
    explore (sorted([s[i:] for i in range(len(s))]))
    return Edges

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

