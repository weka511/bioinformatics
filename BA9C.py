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
            if end_node_labelled_edge == None:
                assert(not next_avail_node in Trie)
                Trie[next_avail_node] = []
                Trie[currentNode].append((currentSymbol,j,next_avail_node))
                currentNode = next_avail_node
                next_avail_node+=1                
            else:
                currentNode = end_node_labelled_edge
    
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
        def explore(node=0,symbol='',path=[],symbols=[],origin=0,index=-1):
            if len(Trie[node])==0:
                Paths.append(path + [node])
                Strings.append(''.join(symbols +[symbol]))
                Origins.append((origin,index))
            elif len(Trie[node])==1:
                next_symbol,_,next_node=Trie[node][0]
                explore(node=next_node,symbol=next_symbol,path=path+[node],symbols=symbols+[symbol],origin=origin,index=index)
            else:
                trie = Trie[node]
                for i in range(len(trie)):
                    next_symbol,_,next_node = trie[i]
                    explore(node=next_node,symbol=next_symbol,origin=node,index=i)
 
        Origins   = []
        Paths     = []
        Strings   = [] 
        explore()
        return Paths,Strings,Origins
        
    ROOT                  = 0         
    Trie,labels           = ConstructModifiedSuffixTrie(Text)
    Paths,Strings,Origins = get_non_branching_paths(Trie)
    for path,string,origin in zip(Paths,Strings,Origins):
        if len(path)<2: break
        node,index        = origin
        s,j,n             = Trie[node][index]
        Trie[node][index] = (string,j,-1)
        for del_node in path[0:]:
            del Trie[del_node]
        
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
            if next_node>0: 
                dfs(Trie,next_node,depth+1)
    print ('root')
    dfs(Trie)
    
def PrintEdges(Trie):
    
    # dfs
    #
    # Depth First Search, printing nodes as we go
    
    def dfs(Trie,currentNode=0):
        for symbol,pos,next_node in Trie[currentNode]:
            print ('{0}'.format(symbol))
            if next_node>0: 
                dfs(Trie,next_node)
    dfs(Trie)    

def ConstructSuffixTreeEdges(s):
    def explore(suffixes):
        if len(suffixes)==0: return
        prefixes = list(set([s[0] for s in suffixes]))
        if len(prefixes)==len(suffixes):
            for s in suffixes:
                edges.append(s)
        else:
            for p in prefixes:
                subset = [s[1:] for s in suffixes if s[0]==p and len(s)>len(p)]
                if len(subset)>1:
                    edges.append(p)
                    explore(subset)
                elif len(subset)==1:
                    edges.append(subset[0])
    edges=[]
    explore (sorted([s[i:] for i in range(len(s))]))
    return edges

if __name__=='__main__':
    for edge in ConstructSuffixTreeEdges('ATAAATG$'):
        print(edge)
    #PrintEdges(Trie)
    #Trie1 = ConstructSuffixTree(
        #'ATTTCAACGGTCTAAACACGGGGGGTGCGATGAGGAGGGTAATCTGGGACTCACTCCCAACATCACTGCGGCTCGGAACTCCTACGACTGTTCACCAAGCCAATGAAGTTTTATTTCAGTCATCAACGCGTCCGTGTAGCATTAGATATCAAAGGTGGTTTTGGTGCTGTGGTCGGCCGACCAATGTAACAAGCTCATCCTCTTCCCAGTGCGATAGGTTCGAATCGATATTAATATGCAAACAAAACCGACGTCCAAGCGACAACTATTTCAGCGTGAATGCCAGGCAAGAGCCGCCCCAAACCTTCCTCAACTCTCGCACCCACCTGACTCTCCTTACATTTGGGCTTCGGCGTTTAAGGCGCCGTTAACGTCGTTGTCAGGCGAGAATAGCCTATACGTCCCTATCTATGACTTTTTCTGATACATCTTCTGAGGAGGTATTGGTATAGGCTTGCATACACGTACATAAAGAGGATAAAATTCCGAGTGGCATTACAGACATGGCTCTAAGCTGGATTGGAAGAAATCTGTATGGTTTAGGAGGGTGGGTAGTCATAAGATCGCATTAGTAGCAACGGCCGGGTTATCTTTTGCTGTTCTGGCCTTTTCCGGTTTGAACCTATGCTAGCGTCAACATACACAGTTCTCTTTTGAGCCGGCAAGATCTCAATATACCGTCCGCCCCTGGCCTAGAGCTCTATTATATGCACTACGTGAGGAATGACTTGGCAAAAATCTGTCAGTGGCGCAACGTTGAAATGGAGAAACCCATGGTTTTACCAACTTTCTGTCGTGTAAATTTGCAGATCTCACCATGTCCGCCGCACCTTCCTTTGGATTAACTTCCCAGTGAGGATCACTCGGCCATCACGTCTGGAATCATACGAACACACGATAGCTCGCGTTTACGCC$'
    #)
    #PrintEdges(Trie1)    