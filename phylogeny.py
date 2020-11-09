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

# Phylogeny -- http://rosalind.info/problems/topics/phylogeny/


#  tree -- Completing a Tree 
#
# Given: A positive integer n (n<=1000) and an adjacency list corresponding to a graph on n nodes that contains no cycles.
#
# Return: The minimum number of edges that can be added to the graph to produce a tree.
#         This is the number of independent components - 1
def CompleteTree(n,adj):
    # create_twigs
    #
    # Build dictionary tow show which node is linked to what
    
    def create_twigs():        
        twigs = {i:set() for i in range(1,n+1)}
        for a,b in adj:
            twigs[a].add(b)
            twigs[b].add(a)
        return twigs
    
    # find_component
    #
    # Find one component of graph
    
    def find_component(start):
        component = []      # The component being built 
        todo      = set()   # Nodes being considered for inclusion
        todo.add(start)
        while len(todo)>0:
            current = todo.pop()
            component.append(current)
            for node in twigs[current]:
                if node not in component:
                    todo.add(node)
        for c in component:
            del twigs[c]
        return component    

    twigs = create_twigs()
    components = []
    while len(twigs)>0:
        components.append(find_component(list(twigs.keys())[0]))
    return len(components)-1

def chbp(species,character_table):
    pass

# cstr 
#
#  Creating a Character Table from Genetic Strings  http://rosalind.info/problems/cstr/

def cstr(strings):
    def trivial(split):
        if len(split)<2: return True
        for k,v in split.items():
            if v<2: return True
        return False
    
    choices = [[] for s in strings[0]]
    counts  = [{} for s in strings[0]]
    for i in range(len(strings[0])):
        for s in strings:
            if not s[i] in choices[i]:
                choices[i].append(s[i])
            if s[i] in counts[i]:
                counts[i][s[i]] += 1
            else:
                counts[i][s[i]] = 1
    
    splits=[]
    for i in range(len(strings[0])):
        split = {}
        for c in choices[i]:
            split[c] = 0
        for s in strings:
            for c in choices[i]:
                if s[i]==c:
                    split[c]+=1
        splits.append(split)
    result=[]
    for i in range(len(strings[0])):
        character = []
        split = splits[i]
        if not trivial(split):
            chs = list(split.keys())
            for s in strings:
                character.append('0' if s[i]==chs[0] else '1')  
            result.append(''.join(character))
    return result

# ctbl Creating a Character Table  http://rosalind.info/problems/ctbl/

def CharacterTable(tree):
    def create_character(split_species):
        character=[]
        for s in species:
            character.append(1 if s in split_species else 0)
        return ''.join([str(c) for c in character])
    
    species=[spec.name for spec in tree.find_elements(terminal=True)]
    species.sort()

    clades=[clade for clade in tree.find_clades(terminal=False)]
    # we iterate over all Clades except the root
    return [create_character([spec.name for spec in split.find_elements(terminal=True)]) for split in clades[1:]]

# NumberBinaryTrees
#
# cunr  Counting Unrooted Binary Trees
# root  Counting Rooted Binary Trees
# See http://carrot.mcb.uconn.edu/~olgazh/bioinf2010/class16.html

def NumberBinaryTrees(n,rooted=True):
    N = 1
    m = 2*n-3 if rooted else 2*n-5
    while m>1:
        N *=m
        m -= 2
    return N

class UnrootedBinaryTree:
    @classmethod
    # EnumerateUnrootedBinaryTrees
    #
    # Given: A collection of species names representing n taxa.
    #
    # Return: A list containing all unrooted binary trees whose leaves are these n
    #         taxa. Trees should be given in Newick format, with one tree on each line; 
    #         the order of the trees is unimportant.
    #
    # Idea: all rooted trees with a given number of leaves are isomorphic if we
    #       ignore the labels of the leaves and nodes. Therfore it is enough to 
    #       build a tree with 3 leaves, and keep adding one leaf at a time in all available positions.
    
    def Enumerate(cls,species):
        def enumerate(n):
            if n==3:
                return [cls({0:[species[0], species[1], species[2]]})]
            else:
                return [cls.insert(species[n-1],edge,graph) for graph in enumerate(n-1) for edge in graph.Edges()] 
            
        return enumerate(len(species))
    
    # insert
    #
    # Create a rooted tree by adding one nre inernal node and a leaf to a specifed edge
    @classmethod
    def insert(cls,species,edge,graph):
        nextNode = max(list(graph.adj.keys())) + 1
        n1,n2    = edge
        adj      = {nextNode: [species,n2]}
 
        for node,links in graph.adj.items():
                adj[node] = [nextNode if ll==n2 else ll for ll in links] if node==n1 else links        
 
        return cls(adj)
    
    def __init__(self,adj):
        self.adj = adj
    
    def __str__(self):
        return self.bfs_newick()
    
    # bfs_newick
    #
    # Create Newick representation by best first search
    
    def bfs_newick(self,node=0):
        newick = []
        for child in  self.adj[node]:
            if type(child)==int:
                newick.append(self.bfs_newick(node=child))
            else:
                newick.append(child)
        representation = ','.join(newick) 
        return f'({representation})'
    
    def Edges(self):
        for a,b in self.adj.items():
            for c in b:
                yield a,c