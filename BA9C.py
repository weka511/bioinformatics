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

# ConstructModifiedSuffixTrie
#
# Construct Modified Suffix Trie, as described in 
# Charging Station Bioinformatics Algorithms Vol II page 165
#
# See also https://sandipanweb.wordpress.com/2017/05/10/suffix-tree-construction-and-the-longest-repeated-substring-problem-in-python/

def ConstructModifiedSuffixTrie(Text):
    
    # edgeLabelled
    #
    # Find edge labelled with symbol
    #
    # Inputs: currentNode
    #         currentSymbol
    #
    # Returns: node that edge links to, or None
    
    def edgeLabelled(currentNode, currentSymbol):
        if currentNode in Trie:
            for symbol,_,node in Trie[currentNode]:
                if symbol==currentSymbol:
                    return node
        return None
    
    # isLeaf
    #
    # Verify that node is a leaf of Trie
    
    def isLeaf(node):
        return len(Trie[node])==0
    
    ROOT            = 0    
    Trie            = {ROOT:[]}
    labels          = {}
    next_avail_node = 1
    
    for i in range(len(Text)):
        currentNode = ROOT
        for j in range(i,len(Text)):
            currentSymbol = Text[j]
            end_node_labelled_edge = edgeLabelled(currentNode, currentSymbol)
            if end_node_labelled_edge != None:
                currentNode = end_node_labelled_edge
            else:
                assert(not next_avail_node in Trie)
                Trie[next_avail_node] = []
                Trie[currentNode].append((currentSymbol,j,next_avail_node))
                currentNode = next_avail_node
                next_avail_node+=1
        if isLeaf(currentNode):
            labels[currentNode] = i
    return Trie,labels

# ConstructSuffixTree
#
# Construct the suffix tree of a string, as described in Charging Station
# Bioinfornmatics Algorithms Voll II page 167
#
# Inputs: A string Text.
#
# Return: The strings labeling the edges of SuffixTree(Text).

def ConstructSuffixTree(Text):
                
    def get_non_branching_paths(Trie):
        def explore(node=0,symbol='',path=[],symbols=[]):
            if len(Trie[node])==0:
                #if len(path)>0:
                Paths.append(path + [node])
                Strings.append(''.join(symbols +[symbol]))
            elif len(Trie[node])==1:
                next_symbol,_,next_node=Trie[node][0]
                return explore(node=next_node,symbol=next_symbol,path=path+[node],symbols=symbols+[symbol])
            else:
                for next_symbol,_,next_node in Trie[node]:
                    explore(node=next_node,symbol=next_symbol)
 
        
        Paths     = []
        Strings   = []
        #explored  = set()
        
        explore()
        
        #for node in sorted(Trie.keys()):
            #path = explore(node)
            #if len(path)>0:
                #Paths.append(path)
        print (len(Paths))
        return Paths,Strings
        
    ROOT        = 0         
    Trie,labels = ConstructModifiedSuffixTrie(Text)
    for path in get_non_branching_paths(Trie):
        symbols = []
        first   = True
        j0      = None
        for node in path:
            if len(Trie[node])>0:
                s,j,_ = Trie[node][0]
                if first:
                    j0     = j
                    first = False
                symbols.append(s)
            del Trie[node]
        Trie[path[0]]= [(''.join(symbols),j0,-1)]
        
    return Trie

# PrintTrie
#
# Print Trie

def PrintTrie(Trie):
    
    # dfs
    #
    # Depth First Search, printing nodes as we go
    
    def dfs(Trie,currentNode=0,depth=1):
        for symbol,pos,next_node in Trie[currentNode]:
            if next_node<0: break
            padding = ''.join(['-']*depth)
            print ('{0} {1}, {2}, {3}'.format(padding,symbol,pos,next_node))
            dfs(Trie,next_node,depth+1)
    print ('root')
    dfs(Trie)

if __name__=='__main__':
    #for k,v in ConstructSuffixTree('ATAAATG$').items():
        #print (k,v)    
    #Trie,labels=ConstructModifiedSuffixTrie('panamabananas$')
    #Trie,labels=ConstructModifiedSuffixTrie('ATAAATG$')
    #PrintTrie(Trie)
    #Trie1 = ConstructSuffixTree('panamabananas$')
    #PrintTrie(Trie1)
    #trie,labels = ConstructModifiedSuffixTrie('panamabananas$')
    #print (trie)
    #print (labels)
    
    Trie2 = ConstructSuffixTree('ATAAATG$')
    PrintTrie(Trie2)    