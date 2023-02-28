#!/usr/bin/env python
#   Copyright (C) 2020-2021 Greenweaves Software Limited

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
from   rosalind import LabelledTree
from   random   import randrange
from   newick   import newick_to_adjacency_list
from   numpy    import argmin,argmax
from   fasta    import FastaContent
from   helpers  import  flatten

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
    def dfs(tree):
        id       = tree['id']
        name     = tree['name']
        children = tree['children']
        parentid = tree['parentid']
        if len(name)==0:
            adj[id]=[]
        if parentid>-1:
            adj[parentid].append(id if len(name)==0 else name)
        for child in children:
            dfs(child)
    dfs(tree)
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

# SmallParsimony
#
# Find the most parsimonious labeling of the internal nodes of a rooted tree.
#
# Given: An integer n followed by an adjacency list for a rooted binary tree with n leaves labeled by DNA strings.
#
# Return: The minimum parsimony score of this tree, followed by the adjacency list of the tree
#         corresponding to labeling internal nodes by DNA strings in order to minimize the parsimony score of the tree.

def SmallParsimony(T,alphabet='ATGC'):

    # SmallParsimonyC Solve small parsimony for one character

    def SmallParsimonyC(Character):

        # get_ripe
        #
        # Returns: a node that is ready fpr preocssing

        def get_ripe():
            for v in T.get_nodes():
                if not processed[v] and v in T.edges:
                    for e,_ in T.edges[v]:
                        if e>v: continue
                        if not processed[e]: break

                    return v
            return None

        #  calculate_s
        #     Calculate score if node v is set to a specified symbol
        #     Parameters:
        #        symbol The symbol, e.g. 'A', not the index in alphabet
        #        v      The node

        def calculate_s(symbol,v):
            # delta
            #
            # Complement of Kronecker delta
            def delta(i):
                return 0 if symbol==alphabet[i] else 1

            def get_min(e):
                return min(s[e][i]+delta(i) for i in range(len(alphabet)))

            return sum([get_min(e) for e,_ in T.edges[v]])

        # update_assignments
        #
        # Parameters:
        #     v
        #     s
        def update_assignments(v,s):
            if not v in assignments.labels:
                assignments.labels[v]=''
            index = 0
            min_s = float('inf')
            for i in range(len(s)):
                if s[i] < min_s:
                    min_s = s[i]
                    index = i
            assignments.set_label(v,assignments.labels[v]+alphabet[index])
            return alphabet[index]

        # backtrack
        #
        # Process internal node of tree top down, starting from root

        def backtrack(v, current_assignment):
            for v_next,_ in T.edges[v]:
                if T.is_leaf(v_next): continue
                if not v_next in assignments.labels:
                    assignments.labels[v_next]=''
                min_score = min([s[v_next][i] for i in range(len(alphabet))])
                indices   = [i for i in range(len(alphabet)) if s[v_next][i]==min_score ]
                matched   = False
                for i in indices:
                    if alphabet[i]==current_assignment:

                        matched = True
                        assignments.set_label(v_next,assignments.labels[v_next]+current_assignment)
                        backtrack(v_next,current_assignment)
                if not matched:
                    # Black magic alert: I am not clear why the introduction of random numbers
                    # helps here. Maybe it stops the tree being biased towatds the first strings
                    # in the alphabet.
                    next_assignment = alphabet[indices[randrange(0,(len(indices)))]]
                    assignments.set_label(v_next,assignments.labels[v_next]+next_assignment)
                    backtrack(v_next,next_assignment)

        processed = {}
        s         = {}

        # Compute scores for a leaves, and mark internal notes unprocessed
        for v in T.get_nodes():
            if T.is_leaf(v):
                processed[v]=True
                s[v]        = [0 if symbol==Character[v] else float('inf') for symbol in alphabet]
            else:
                processed[v]=False

        # Process ripe (unprocessed, but whose children have been processed)
        # until there are none left
        # Keep track of last node as we will use it to start backtracking
        v = get_ripe()
        while not v == None:
            processed[v] = True
            s[v]         = [calculate_s(symbol,v) for symbol in alphabet ]
            v_last       = v
            v            = get_ripe()

        backtrack(v_last,update_assignments(v_last,s[v_last]))
        return min([s[v_last][c] for c in range(len(alphabet))])

    assignments = LabelledTree(T.N)
    assignments.initialize_from(T)

    return sum([SmallParsimonyC([v[i] for l,v in T.labels.items()]) for i in range(len(T.labels[0]))]),assignments

# alph
#
# Given: A rooted binary tree T on n  species, given in Newick format, followed by a multiple alignment of m
#        augmented DNA strings having the same length (at most 300 bp) corresponding to the species
#        and given in FASTA format.
#
# Return: The minimum possible value of dH(T), followed by a collection of DNA strings to be
#         assigned to the internal nodes of T that will minimize dH(T).

def alph(T,Alignment,Alphabet=['A','T','C','G','-']):

    # create_fixed_alignments
    #
    # Extract dictionary of leaves from Alignment
    #
    # Returns: length of any string in alignment, plus dictionary of leaves
    def create_fixed_alignments():
        Leaves = {}
        k      = None
        for i in range(0,len(Alignment),2):
            Leaves[Alignment[i]] = Alignment[i+1]
            if k==None:
                k = len(Alignment[i+1])
            else:
                assert k==len(Alignment[i+1]),f'Alignments should all have same length.'

        return k,Leaves

    # SmallParsimony
    #
    # This is the Small Parsimony algorihm from Pevzner and Compeau, which
    # processes a single character
    #
    # Parameters:
    #      l       Index of character in Alignment
    # Returns:     Score of best assignment, plus an assignment of character that provides this score

    def SmallParsimony(l):

        # is_ripe
        #
        # Determine whether now is ready for processing
        # A ripe node is one that hasn't been processed,
        # but its children have

        def is_ripe(v):
            for child in Adj[v]:
                if not Tag[child]: return False
            return True

        # find_ripe
        #
        # Find list of nodes that are ready to be processed
        #
        # Input:   A list of nodes
        # Returns: Two lists, those ready for processing, and those which are not

        def find_ripe(Nodes):
            Ripe   = []
            Unripe = []
            for v in Nodes:
                if is_ripe(v):
                    Ripe.append(v)
                else:
                    Unripe.append(v)
            return Ripe,Unripe

        # delta
        #
        # The delta function from Pevzner and Compeau: not the Kronecker delta

        def delta(i,j):
            return 0 if i==j else 1

        # get_distance
        #
        # Get total distance of node from its children assuming one trial assignmnet
        #
        # Parameters:
        #     v       Current node
        #     k       Index of character for trial

        def get_distance(v,k):

            # best_alignment
            #
            # Find best alignment with child (measured by varying child's index) given
            # the current choice of character in this node
            #
            # Parameters:
            #      k       Trial alignmnet for this node

            def best_alignment(child):
                return min([s[child][i] + delta(i,k)  for i in range(len(Alphabet))])

            return sum([best_alignment(child) for child in Adj[v]])

        # backtrack
        #
        # Perform a depth first search through all nodes to determive alignmant.
        # Parameters:
        #     root    Root node
        #     s       Scores for all possible best assignments to all nodes
        # Returns:
        #     score    Score of best assignment,
        #     ks       For each node the assignment of character that provides this score
        #              represented an an index into alphabet
        #
        #
        #       Comment by  Robert Goldberg-Alberts.
        #       The Backtrack portion of the code consists of a breath first tracking through the tree from
        #       the root in a left to right fashion through the nodes (sons and daughters)
        #       row after row until you finally reach the leaves. The leaves already have values assigned to them from the data
        #       At the root, determine whether the value of the node is A, C, T, G by taking the minimum value of the
        #       four numbers created for the root set. Break ties by selecting from the ties at random.
        #       After that, for subsequent nodes take the minimum of each value at a node and determine if
        #       there are ties at the minimum. Check to see if the ancestor parent of that node has a value
        #       that is contained in the eligible nucleotides from the node. If it IS contained there force the
        #       ancestor value for that node.
        #       Continue in that fashion until all the internal nodes have values assigned to them.

        def backtrack(root,s):
            def dfs(node,k,parent_score):

                def match(i,j,child_scores):
                    return parent_score == child_scores[0][i] + child_scores[1][j]

                if len(Adj[node])==0: return
                children           = Adj[node]
                child_scores_delta = [[s[child][i] + delta(i,k)  for i in range(len(Alphabet))] for child in children]
                child_scores_raw   = [[s[child][i]  for i in range(len(Alphabet))] for child in children]
                candidates         = [(i,j,child_scores_raw) for i in range(len(Alphabet)) for j in range(len(Alphabet)) \
                                      if match(i,j,child_scores_delta)]
                selection          = candidates[randrange(len(candidates))]
                scores_children    = [selection[2][i][selection[i]] for i in range(len(children))]
                for i in range(len(children)):
                    ks[children[i]] = selection[i]

                for i in range(len(children)):
                    dfs(children[i],ks[children[i]],scores_children[i])

            ks       = {}
            index    = argmin(s[root])
            score    = s[root][index]
            ks[root] = index
            dfs(root,index,score)
            return score, ks




        s      = {}           # Scores for nodes
        Tag    = {}           # Nodes that have been processed
        ToBeProcessed   = []  # Nodes that have yet to be processed

        # Partition nodes into two groups: leaves are easily processed,
        # the others are all marked as unprocessed
        for v in Adj.keys():
            if v in Leaves:
                char      = Leaves[v][l]
                s[v]      = [0 if Alphabet[k]==char else float('inf') for k in range(len(Alphabet))]
                Tag[v]    = True
            else:
                Tag[v] = False
                ToBeProcessed.append(v)

        Ripe,ToBeProcessed = find_ripe(ToBeProcessed)
        while len(Ripe)>0:
            for v in Ripe:
                s[v]      = [get_distance(v,k) for k in range(len(Alphabet))]
                Tag[v]    = True
            Ripe,ToBeProcessed = find_ripe(ToBeProcessed)

        assert len(ToBeProcessed)==0,'If there are no ripe nodes, ToBeProcessed should be exhausted'
        return backtrack(v,s)

    Adj        = newick_to_adjacency_list(T)
    L,Leaves   = create_fixed_alignments()
    Assignment = {a:[] for a in Adj.keys()}
    d          = 0
    assert len([node for node,value in Adj.items() if len(value)==0 and node not in Leaves])==0,\
           f'Some nodes are leaves, but have no strings in alignment'

    for l in range(L):
        score,ks = SmallParsimony(l)
        d       += score
        for v,index in ks.items():
            Assignment[v].append(Alphabet[index])

    return d,[(f'{a}',''.join(b)) for a,b in Assignment.items() if len(Adj[a])!=0]

#  chbp Character-Based Phylogeny
#
#  Strategy: sort character table on entropy, then use each character to divide clades into two.

def chbp(species,character_table):
    # Clade
    #
    # This class represents one clade or taxon

    class Clade:
        def __init__(self,taxa):
            self.taxa = [s for s in taxa]

        def is_trivial(self):
            return len(self.taxa)==0

        def is_singleton(self):
            return len(self.taxa)==1

        # newick
        #
        # Convert to string in Newick format

        def newick(self):
            def conv(taxon):
                if type(taxon)==int:
                    return species[taxon]
                else:
                    return taxon.newick()
            if self.is_singleton():
                return conv(self.taxa[0])
            else:
                return '(' + ','.join(conv(taxon) for taxon in self.taxa) +')'

        # split
        #
        # Split clade in two using character: list of taxa is replaced by two clades
        #
        # Returns True if clade has been split into two non-trivial clades
        #         False if at least one clade would be trivial--in which case clade is unchanged
        #
        def split(self,character):
            left  = []
            right = []
            for i in self.taxa:
                if character[i]==0:
                    left.append(i)
                else:
                    right.append(i)
            leftTaxon  = Clade(left)
            rightTaxon = Clade(right)
            if leftTaxon.is_trivial(): return False
            if rightTaxon.is_trivial(): return False
            self.taxa = [leftTaxon,rightTaxon]
            return True

        # splitAll
        #
        # Split clade using character table
        def splitAll(self,characters,depth=0):
            if depth<len(characters):
                if self.split(characters[depth]):
                    for taxon in self.taxa:
                        taxon.splitAll(characters,depth+1)
                else:
                    self.splitAll(characters,depth+1)


    # Calculate entropy of a single character

    def get_entropy(freq):
        if freq==0 or freq==n: return 0
        p1 = freq/n
        p2 = 1-p1
        return - p1 *np.log(p1) - p2 * np.log(p2)



    n               = len(species)
    entropies       = [get_entropy(sum(char)) for char in character_table]
    entropy_indices = np.argsort(entropies)
    characters      = [character_table[i] for i in entropy_indices[::-1]]
    indices         = list(range(len(species)))
    root            = Clade(indices)
    root.splitAll(characters)
    return f'{root.newick()};'

#  RSUB  	Identifying Reversing Substitutions
#
# Given: A rooted binary tree T with labeled nodes in Newick format, followed by a collection of at most
#        100 DNA strings in FASTA format whose labels correspond to the labels of T.
#
#        We will assume that the DNA strings have the same length, which does not exceed 400 bp).
#
# Return: A list of all reversing substitutions in T (in any order), with each substitution encoded by the following three items:
#
#    the name of the species in which the symbol is first changed, followed by the name of the species in which it changes back to its original state
#    the position in the string at which the reversing substitution occurs; and
#    the reversing substitution in the form original_symbol->substituted_symbol->reverted_symbol.

def rsub(T,Assignments):

    # find_path
    #
    # Find path from the root down to a specified leaf

    def find_path(leaf):
        Path   = [leaf]
        parent = Parents[leaf]
        while len(parent)>0:
            Path.append(parent)
            if parent in Parents:
                parent = Parents[parent]
            else:
                break
        return Path[::-1]

    # FindReversingSubstitutions
    #
    # Find reversion substitutions in one specified path trhough tree,
    # affecting a specified position in the strings
    #
    # Parameters: Path    Path to be searched
    #             pos     position in tree
    # Strategy:  build up history of changes, and search back whenever a change is detected.
    def FindReversingSubstitutions(Path,pos):
        History   = [Characters[Path[0]][pos]]
        Names     = Path[0:1]
        Reverses  = []
        for taxon in Path[1:]:
            current = Characters[taxon][pos]
            if current==History[-1]: continue
            History.append(current)
            Names.append(taxon)
            if len(History)>2 and History[-3]==History[-1]: # we have a reverse
                Reverses.append((Names[-2],Names[-1],pos+1,History[-3],History[-2],History[-1]))

        return Reverses

    # create_parents

    # Invert Ajacency list to we have the parent of each child

    def create_parents(Adj):
        Product = {node:[] for node in flatten(Adj.values())}
        for parent,children in Adj.items():
            for child in children:
                Product[child] = parent
        return Product

    # get_unique
    #
    # Convert list of lists into a single list and remove duplicate elements

    def get_unique(list_of_lists):
        return list(set(flatten(list_of_lists)))

    Adj,root = newick_to_adjacency_list(T,return_root=True)
    fc       = FastaContent(Assignments)
    Characters = fc.to_dict()            # So we can find character for each species
    _,string = fc[0]
    m        = len(string)
    Parents  = create_parents(Adj)
    Paths    = [find_path(node) for node in flatten(Adj.values()) if len(Adj[node])==0]

    # Build list of unique reversals.
    return get_unique([subst for subst in [FindReversingSubstitutions(path,pos) for path in Paths for pos in range(m)] if len(subst)>0])

# cset A submatrix of a matrix M is a matrix formed by selecting rows and columns from M and
# taking only those entries found at the intersections of the selected rows and columns.
# We may also think of a submatrix as formed by deleting the remaining rows and columns from M
#
# Given: An inconsistent character table C on at most 100 taxa.
#
# Return: A submatrix of C representing a consistent character table on the same taxa
#         and formed by deleting a single row of C.

def cset(table):
    # get_split
    #
    # Used to split indices of character (row) into two groups, one for each allele
    # First we yield all indices corresponding to 0, then those to 1

    def get_splits(character):
        for allele in [0,1]:
            yield set(i for i, c in enumerate(character) if c == allele)

    # conflicts_with
    #
    # Determine whether two characters are in conflict
    # We iterate through all the splits of each character.
    # If any pair of splits consists of two disjoint subsets,
    # the characters are compatible.
    def conflicts_with(c1, c2):
        for part1 in get_splits(c1):
            for part2 in get_splits(c2):
                if len(part1.intersection(part2)) == 0: return False

        return True

    n         = len(table)
    Conflicts = [0 for _ in range(n)]   # Count number of times each row conflicts with another
    for i in range(n):
        for j in range(i+1,n):
            if conflicts_with(table[i],table[j]):
                Conflicts[i] += 1
                Conflicts[j] += 1

    return [table[row] for row in range(n) if row!=argmax(Conflicts)]

#  cntq Counting Quartets

def cntq(n,newick):

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

    def bfs(subtree,leaves):
        for node in adj[subtree]:
            if type(node)==str:
                leaves.append(node)
            else:
                bfs(node,leaves)
    def pair(leaves):
        for i in range(len(leaves)):
            for j in range(i+1,len(leaves)):
                yield [leaves[i],leaves[j]] if leaves[i]<leaves[j] else [leaves[j],leaves[i]]

    adj             = create_adj(parse(newick))
    taxa            = [leaf for children in adj.values() for leaf in children if type(leaf)==str]
    splitting_edges = [(key,child) for key,value in adj.items() for child in value if type(child)==int]
    Quartets        = []
    for _,node in splitting_edges:
        leaves       = []
        bfs(node,leaves)
        other_leaves = [leaf for leaf in taxa if leaf not in leaves]
        for pair1 in pair(leaves):
            for pair2 in pair(other_leaves):
                quartet = pair1 + pair2 if pair1[0]<pair2[0] else pair2 + pair1
                Quartets.append(quartet)
    Quartets.sort()
    Unique =[Quartets[0]]
    for i in range(1,len(Quartets)):
        if Quartets[i]!=Unique[-1]:
            Unique.append(Quartets[i])
    return len(Unique),Unique
