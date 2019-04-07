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
# Construct Modified Suffix Trie, as described in Charging Station
# Bioinfornmatics Algorithms Voll II page 165

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
                Trie[next_avail_node] = []
                Trie[currentNode].append((currentSymbol,j,next_avail_node))
                currentNode = next_avail_node
                next_avail_node+=1
        if len(Trie[currentNode])==0:
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
    # dfs
    #
    # Search Trie for non-branching paths
    #
    # Inputs: currentNode
    #         path
    def dfs(currentNode,path=[]):
        if len(Trie[currentNode])<2: # Non-branching or end
            path = path + [currentNode]
            if len(Trie[currentNode])==0: #end of branch
                symbols=[]
                for step in path:
                    if (len(Trie[step]))>0:
                        symbol,_,_=Trie[step][0]
                        symbols.append(symbol)
                    #else:
                        #symbols.append(str(labels[step]))
                    
                return ''.join(symbols)
        else:
            path=[]
        for symbol,pos,next_node in Trie[currentNode]:
            symbols=dfs(next_node,path)
            if symbols!=None:
                print (symbols)
                Trie[currentNode]=[(symbols,pos,next_node)]   
            
    Trie,labels = ConstructModifiedSuffixTrie(Text)
    Paths=[]
    dfs(0)
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
            padding = ''.join(['-']*depth)
            print ('{0} {1}, {2}, {3}'.format(padding,symbol,pos,next_node))
            dfs(Trie,next_node,depth+1)
    print ('root')
    dfs(Trie)

if __name__=='__main__':
    #for k,v in ConstructSuffixTree('ATAAATG$').items():
        #print (k,v)    
    #Trie,labels=ConstructModifiedSuffixTrie('panamabananas$')
    #PrintTrie(Trie)
    #Trie1 = ConstructSuffixTree('panamabananas$')
    #PrintTrie(Trie1)
    #trie,labels = ConstructModifiedSuffixTrie('panamabananas$')
    #print (trie)
    #print (labels)
    
    Trie2 = ConstructSuffixTree('ATAAATG$')
    PrintTrie(Trie2)    