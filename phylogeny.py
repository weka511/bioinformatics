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

import re
import numpy as np

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

#  qrt Incomplete Characters
#
# Given: A partial character table C
#
# Return: The collection of all quartets that can be inferred from the splits corresponding to the underlying characters of C

def qrt(taxa,characters):
    
    def tuples(n):
        for i in range(n):
            for j in range(n):
                if i==j: continue
                for k in range(n):
                    if k in [i,j]: continue
                    for l in range(n):
                        if l in [i,j,k]: continue
                        if i<j and k<l and i<k:
                            yield i,j,k,l
                        
    def isConsistent(selector):
        for char in characters:
            character = [char[i] for i in selector]
            if any(c is None for c in character): continue
            if character[0]==character[1] and character[2]==character[3] and character[0]!=character[2]: return True
        return False
    
    for (i,j,k,l) in tuples(len(taxa)):     
        selector =  [i,j,k,l]
        if isConsistent(selector):
            yield [taxa[m] for m in selector]
            
# snarfed from https://stackoverflow.com/questions/51373300/how-to-convert-newick-tree-format-to-a-tree-like-hierarchical-object
def parse(newick,start=0):
    tokens = re.findall(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")
    
    def recurse(nextid = start, parentid = -1): # one node
        thisid = nextid;
        children = []

        name, length, delim, ch = tokens.pop(0)
        if ch == "(":
            while ch in "(,":
                node, ch, nextid = recurse(nextid+1, thisid)
                children.append(node)
            name, length, delim, ch = tokens.pop(0)
        return {"id": thisid, "name": name, "length": float(length) if length else None, 
                "parentid": parentid, "children": children}, delim, nextid

    return recurse()[0]

def create_adj(tree):
    adj = {}
    def bfs(tree):
        id       = tree['id']
        name     = tree['name']
        children = tree['children']
        parentid = tree['parentid']
        if len(name)==0:
            adj[id]=[]
        if parentid>-1:
            adj[parentid].append(id if len(name)==0 else name)     
        for child in children:
            bfs(child)
    bfs(tree)
    return adj

#  SPTD Phylogeny Comparison with Split Distance

def sptd(species,newick1,newick2):
    def replace_leaves(adj):
        return {parent:sorted([seiceps[child] if child in seiceps else child for child in children]) for parent,children in adj.items() }
    def edges(adj):
        for parent,children in adj.items():
            for child in children:
                if child >= n:
                    yield parent,child
    def splits(adj,min_size=2):
        def find_leaves(node,path=[]):
            for child in adj[node]:
                if child<n:
                    path.append(child)
                else:
                    find_leaves(child,path=path)            
  
        for parent,child in edges(adj):
            s1 = []
            find_leaves(child,s1)#[leaf for leaf in find_leaves(child)]
            if len(s1)<min_size: continue
            s2 = [leaf for leaf in range(n) if not leaf in s1]
            yield sorted(s1),sorted(s2)
  
           
    def ds(adj1,adj2):
        shared  = 0
        splits1 = sorted([s for s,_ in splits(adj1)])
        splits2 = sorted([s for s,_ in splits(adj2)])
        k1      = 0
        k2      = 0
        i1      = splits1[k1]
        i2      = splits2[k2]      
        while k1<len(splits1) and k2<len(splits2): 
            if i1==i2:
                shared += 1
                k1     += 1
                k2     += 1
                if k1<len(splits1) and k2<len(splits2):
                    i1 = splits1[k1]
                    i2 = splits2[k2]
  
            elif i1<i2:
                k1+=1
                if k1<len(splits1):
                    i1 = splits1[k1]                  
            else:
                k2+=1
                if k2<len(splits2):
                    i2 = splits2[k2]  

        return 2*(n-3)- 2* shared
    
    n       = len(species)
    seiceps = {species[i]:i for i in range(n)}

    return ds(replace_leaves(create_adj(parse(newick1,start=n))),
              replace_leaves(create_adj(parse(newick2,start=n))))

# MEND Inferring Genotype from a Pedigree
#
# Given: A rooted binary tree T in Newick format encoding an individual's pedigree 
#       for a Mendelian factor whose alleles are A (dominant) and a (recessive).
#
#       Return: Three numbers between 0 and 1, corresponding to the respective probabilities
#       that the individual at the root of T will exhibit the "AA", "Aa" and "aa" genotypes.

def mend(node):

    # combine
    #
    # Combine two genomes with known probabilities - work out proabilites in next generation
    #
    # NB: the tree is a pedigree, not a phylogeny: the root is the descendent!
    
    def combine(f1,f2):
        return np.sum([[f*f1[i]*f2[j] for f in factors[i][j]] for i in range(n) for j in range(n)],
                      axis=0)
    
    # Probability of each combination in the initial generation, when we know the genome
    
    frequencies = {
        'aa': (0,0,1),
        'Aa': (0,1,0),
        'AA': (1,0,0)
    }
    
    # Probabilty of each combination when we combine two genomes
    
    factors=[#          AA                  Aa/aA                  aa
                [ [1.0, 0.0, 0.0],  [0.50, 0.50, 0.00],  [0.0, 1.0, 0.0] ], #AA
                [ [0.5, 0.5, 0.0],  [0.25, 0.50, 0.25],  [0.0, 0.5, 0.5] ], #Aa/aA
                [ [0.0, 1.0, 0.0],  [0.00, 0.50, 0.50],  [0.0, 0.0, 1.0] ]  #aa
    ]
    
    n = len(frequencies)  # Number of combinations
    
    # If we are at a leaf, we have a known ancestor
    if len(node.nodes)==0:
        try: 
            return frequencies['Aa' if node.name=='aA' else node.name]
        except KeyError:
            return (0,0)
        
    parent_freqs = [mend(parent) for parent in node.nodes]
    parent_freqs = [pp for pp in parent_freqs if len(pp)==n]
    return combine(parent_freqs[0],parent_freqs[1])