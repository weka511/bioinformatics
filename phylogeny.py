#!/usr/bin/env python
#   Copyright (C) 2020-2023 Greenweaves Software Limited

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

'''Phylogeny -- http://rosalind.info/problems/topics/phylogeny/'''


from   unittest import TestCase, main, skip
from   io       import StringIO

import numpy    as np
from   Bio      import Phylo

from   rosalind import LabelledTree
from   random   import randrange
from   newick   import newick_to_adjacency_list, Parser, Tokenizer, Hierarchy
from   fasta    import FastaContent, fasta_out
from   helpers  import flatten, expand


def CompleteTree(n,adj):
    '''
    tree -- Completing a Tree

    Given: A positive integer n (n<=1000) and an adjacency list corresponding to a graph on n nodes that contains no cycles.

    Return: The minimum number of edges that can be added to the graph to produce a tree.
            This is the number of independent components - 1

    '''
    # Yury Zavarin pointed out:
    #    "The main idea is that the particular edges don't matter. The only thing
    #     that matters is the number of edges as a tree on n vertices always has n−1 edges....
    #     We simply count the number of edges, and subtract it from n−1." https://rosalind.info/problems/chbp/
    return  n-1 - len(adj)

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

def create_character(included,species):
    '''
    create_character

    Produce an array of 1s for each species that belongs to the split

    Parameters:
        included
        species

    '''
    def fmt(c):
        return '1' if c else '0'
    return ''.join([fmt(s in included) for s in species])

def CharacterTable(tree):
    '''
    ctbl Creating a Character Table  http://rosalind.info/problems/ctbl/
    '''

    species = sorted([spec.name for spec in tree.find_elements(terminal=True)])
    clades  = [clade for clade in tree.find_clades(terminal=False)]

    # we iterate over all Clades except the root
    return [create_character([spec.name for spec in split.find_elements(terminal=True)],species) for split in clades[1:]]

def NumberBinaryTrees(n,rooted=True):
    '''
    NumberBinaryTrees

    cunr  Counting Unrooted Binary Trees
    root  Counting Rooted Binary Trees
    See http://carrot.mcb.uconn.edu/~olgazh/bioinf2010/class16.html
    '''
    N = 1
    m = 2*n-3 if rooted else 2*n-5
    while m>1:
        N *=m
        m -= 2
    return N

class UnrootedBinaryTree:
    @classmethod
    def Enumerate(cls,species):
        '''
        EnumerateUnrootedBinaryTrees

        Given: A collection of species names representing n taxa.

        Return: A list containing all unrooted binary trees whose leaves are these n
                taxa. Trees should be given in Newick format, with one tree on each line;
                the order of the trees is unimportant.

        Idea: all rooted trees with a given number of leaves are isomorphic if we
              ignore the labels of the leaves and nodes. Therfore it is enough to
              build a tree with 3 leaves, and keep adding one leaf at a time in all available positions.
        '''
        def enumerate(n):
            if n==3:
                return [cls({0:[species[0], species[1], species[2]]})]
            else:
                return [cls.insert(species[n-1],edge,graph) for graph in enumerate(n-1) for edge in graph.Edges()]

        return enumerate(len(species))


    @classmethod
    def insert(cls,species,edge,graph):
        '''
        insert

        Create a rooted tree by adding one internal node and a leaf to a specifed edge
        '''
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


    def bfs_newick(self,node=0):
        '''
        bfs_newick

        Create Newick representation by best first search
        '''
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



def qrt(taxa,characters):
    '''
     qrt Incomplete Characters

     Given: A partial character table C

     Return: The collection of all quartets that can be inferred from the splits corresponding to the underlying characters of C
    '''
    def quartets(n):
        '''Generate all possible quarters made from elements less that n'''
        for i in range(n):
            for j in range(n):
                if i==j: continue
                for k in range(n):
                    if k in [i,j]: continue
                    for l in range(n):
                        if l in [i,j,k]: continue
                        if i<j and k<l and i<k:
                            yield ((i,j),(k,l))

    def isConsistent(selector):
        '''Verify that a selector is consistent: no element is null, and no two elements are equal'''
        for char in characters:
            character = [char[i] for pair in selector for i in pair]
            if (all(c is not None for c in character) and
                character[0]==character[1]            and
                character[2]==character[3]            and
                character[0]!=character[2]):
                return True
        return False

    for selector in quartets(len(taxa)):
        if isConsistent(selector):
            ((i,j),(k,l)) = selector
            yield ((taxa[i],taxa[j]),(taxa[k],taxa[l]))



def qrtd(species,T1,T2):
    '''
    qrtd  Quartet Distance

    Given: A list containing n taxa  and two unrooted binary trees T1 and T2 on the given taxa.
           Both T1 and T2 are given in Newick format.

    Return: The quartet distance dq(T1,T2)
    '''
    species_index = {i:species[i] for i in range(len(species))}
    tree1         = Phylo.read(StringIO(T1), "newick")
    tree2         = Phylo.read(StringIO(T2), "newick")

    z=0





def sptd(species,newick1,newick2):
    '''
    SPTD Phylogeny Comparison with Split Distance

    Given: A collection of at most 3,000 species taxa and two unrooted binary trees T1 and T2 on these taxa in Newick format.

    Return: The split distance dsplit(T1,T2), the number of nontrivial splits
            contained in one unrooted binary tree but not the other.
.
    '''
    def to_edges(tree):
        '''
        Iterate through edges
        '''
        for parent in tree.find_clades():
            for child in parent.clades:
                yield (parent,child)

    def create_non_trivial_splits(tree):
        '''
        create_non_trivial_splits

        We can specify a split by giving one componenet, since the other contains all the other leaves.
        '''
        product = []
        for (parent,child) in to_edges(tree):
            leaves = [clade.name for clade in child.find_clades() if clade.name]
            if len(leaves)>1 and len(leaves)<n-1:  # Non trivial iff both components have length 2 or more
                leaves.sort()
                product.append(leaves)
        product.sort()
        return product

    def ds(splits1,splits2):
        '''
        ds

        This is the workhorse that determines number of shared splits

        Parameters:
            splits1     Splits from one tree, each represented as the list of
                        leaves defining one component of the split. The inner and outer lists
                        are both sorted into ascending order.
            splits2     Splits from the other tree, ordered similaly to splits1

        Returns:
            Count of splits shared by the two trees
        '''
        n_shared = 0
        k1       = 0
        k2       = 0
        i1       = splits1[k1]
        i2       = splits2[k2]
        while k1<len(splits1) and k2<len(splits2):
            if i1==i2:
                n_shared += 1
                k1       += 1
                k2       += 1
                if k1<len(splits1) and k2<len(splits2):
                    i1 = splits1[k1]
                    i2 = splits2[k2]
            elif i1<i2:
                k1 += 1
                if k1<len(splits1):
                    i1 = splits1[k1]
            else:
                assert i1>i2
                k2 += 1
                if k2<len(splits2):
                    i2 = splits2[k2]

        return n_shared

    n             = len(species)
    species_index = {species[i]:i for i in range(n)}
    #  Any unrooted binary tree on n taxa must have n−3 nontrivial splits,
    #  so the two trees have 2(n-3) nontrivial splits
    #  whih comprises non-shared + shared, which must be counted twice
    return 2*(n-3)- 2* ds(create_non_trivial_splits(Phylo.read(StringIO(newick1), "newick")),
                          create_non_trivial_splits(Phylo.read(StringIO(newick2), "newick")))



def mend(node):
    '''
     MEND Inferring Genotype from a Pedigree

     Given: A rooted binary tree T in Newick format encoding an individual's pedigree
           for a Mendelian factor whose alleles are A (dominant) and a (recessive).

           Return: Three numbers between 0 and 1, corresponding to the respective probabilities
           that the individual at the root of T will exhibit the "AA", "Aa" and "aa" genotypes.
    '''

    def combine(f1,f2):
        '''
        combine

        Combine two genomes with known probabilities - work out proabilites in next generation

        NB: the tree is a pedigree, not a phylogeny: the root is the descendent!
        '''
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


def alph(T,Alignment,
         Alphabet=['A','T','C','G','-']):
    '''
    alph

    Given: A rooted binary tree T on n  species, given in Newick format, followed by a multiple alignment of m
           augmented DNA strings having the same length (at most 300 bp) corresponding to the species
           and given in FASTA format.

    Return: The minimum possible value of dH(T), followed by a collection of DNA strings to be
            assigned to the internal nodes of T that will minimize dH(T).

    '''

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
            index    = np.argmin(s[root])
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
    Assignment = {a:[] for a,_ in Adj.items()}
    d          = 0
    assert len([node for node,value in Adj.items() if len(value)==0 and node not in Leaves])==0,\
           f'Some nodes are leaves, but have no strings in alignment'

    for l in range(L):
        score,ks = SmallParsimony(l)
        d       += score
        for v,index in ks.items():
            Assignment[v].append(Alphabet[index])

    return d,[(f'{a}',''.join(b)) for a,b in Assignment.items() if len(Adj[a])!=0]



def chbp(species,character_table):
    '''
    chbp Character-Based Phylogeny

    Strategy: sort character table on entropy, then use each character to divide clades into two.
    '''
    class Clade:
        '''
        Clade

        This class represents one clade or taxon
        '''
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


def rsub(T,Assignments):
    '''
    RSUB  	Identifying Reversing Substitutions

    Given: A rooted binary tree T with labeled nodes in Newick format, followed by a collection of at most
           100 DNA strings in FASTA format whose labels correspond to the labels of T.

           We will assume that the DNA strings have the same length, which does not exceed 400 bp).

    Return: A list of all reversing substitutions in T (in any order), with each substitution encoded by the following three items:

       the name of the species in which the symbol is first changed, followed by the name of the species in which it changes back to its original state
       the position in the string at which the reversing substitution occurs; and
       the reversing substitution in the form original_symbol->substituted_symbol->reverted_symbol.
    '''

    def find_path(leaf):
        Path   = [leaf]
        '''
        find_path

        Find path from the root down to a specified leaf

        '''
        parent = Parents[leaf]
        while len(parent)>0:
            Path.append(parent)
            if parent in Parents:
                parent = Parents[parent]
            else:
                break
        return Path[::-1]


    def FindReversingSubstitutions(Path,pos):
        '''
        FindReversingSubstitutions

         Find reversion substitutions in one specified path trhough tree,
         affecting a specified position in the strings

         Parameters: Path    Path to be searched
                     pos     position in tree
         Strategy:  build up history of changes, and search back whenever a change is detected.
        '''
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



    def create_parents(Adj):
        '''
        create_parents

        Invert Adjacency list to we have the parent of each child
        '''
        Product = {node:[] for node in flatten(Adj.values())}
        for parent,children in Adj.items():
            for child in children:
                Product[child] = parent
        return Product

    def get_unique(list_of_lists):
        '''
        get_unique

        Convert list of lists into a single list and remove duplicate elements
        '''
        return list(set(flatten(list_of_lists)))

    Adj,root = newick_to_adjacency_list(T,return_root=True)
    fc       = FastaContent(Assignments)
    Characters = fc.to_dict()            # So we can find character for each species
    _,string = fc[0]
    m        = len(string)
    Parents  = create_parents(Adj)
    Paths    = [find_path(node) for node in flatten(Adj.values()) if len(Adj[node])==0]

    # Build list of unique reversals.
    return get_unique([subst for subst in [FindReversingSubstitutions(path,pos) for path in Paths for pos in range(m)]
                       if len(subst)>0])



def cset(table):
    '''
    CSET Fixing an Inconsistent Character Set

    A submatrix of a matrix M is a matrix formed by selecting rows and columns from M and
    taking only those entries found at the intersections of the selected rows and columns.
    We may also think of a submatrix as formed by deleting the remaining rows and columns from M

    Given: An inconsistent character table C on at most 100 taxa.

    Return: A submatrix of C representing a consistent character table on the same taxa
            and formed by deleting a single row of C.
    '''


    def get_splits(character):
        '''
        get_split

        Used to split indices of character (row) into two groups, one for each allele
        First we yield all indices corresponding to 0, then those to 1
        '''
        for allele in [0,1]:
            yield set(i for i, c in enumerate(character) if c == allele)


    def conflicts_with(c1, c2):
        '''
        conflicts_with

        Determine whether two characters are in conflict
        We iterate through all the splits of each character.
        If any pair of splits consists of two disjoint subsets,
        the characters are compatible.
        '''
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

    return [table[row] for row in range(n) if row!=np.argmax(Conflicts)]



def cntq(n,newick):
    '''CNTQ Counting Quartets'''
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

    adj             = Hierarchy(newick).create_adj()
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

def expand_character(s):
    return [int(c) if c.isdigit() else None for c in s]

if __name__=='__main__':
    class PhylogenyTestCase(TestCase):

        def test_alph(self):
            fc = FastaContent(['>ostrich',
                               'AC',
                               '>cat',
                               'CA',
                               '>duck',
                               'T-',
                               '>fly',
                               'GC',
                               '>elephant',
                               '-T',
                               '>pikachu',
                               'AA'
                               ])

            d,Assignment = alph('(((ostrich,cat)rat,(duck,fly)mouse)dog,(elephant,pikachu)hamster)robot;',fc.to_list())
            self.assertEqual(8,d)   # matches Rosalind; reset doesn't match, but gives same distance


        def test_chbp(self):
            '''CHBP Character-Based Phylogeny '''
            tree = chbp(['cat', 'dog', 'elephant', 'mouse', 'rabbit', 'rat'],
                        [expand('011101'),
                         expand('001101'),
                         expand('001100')] )
            self.assertEqual('(((cat,rabbit),dog),(rat,(elephant,mouse)));',tree)

        def test_cntq(self):
            '''CNTQ Counting Quartets'''
            n,_ = cntq(6,'(lobster,(cat,dog),(caterpillar,(elephant,mouse)));')
            self.assertEqual(15,n)


        def test_cset(self):
            '''CSET Fixing an Inconsistent Character Set'''
            submatrix = cset([expand('100001'),
                              expand('000110'),
                              expand('111000'),
                              expand('100111')])
            self.assertEqual(3,len(submatrix))
            self.assertIn([0,0,0,1,1,0],submatrix)
            self.assertIn( [1, 1, 1, 0, 0, 0],submatrix)
            self.assertIn([1,0,0,1,1,1],submatrix)

        @skip('#129')
        def test_cstr(self):
            character_table = cstr(['ATGCTACC',
                  'CGTTTACC',
                  'ATTCGACC',
                  'AGTCTCCC',
                  'CGTCTATC'
                  ])
            self.assertEqual(2,len(character_table))
            self.assertIn('10110',character_table)
            self.assertIn('10100',character_table)

        def test_ctbl(self):
            '''
            ctbl Creating a Character Table  http://rosalind.info/problems/ctbl/
            '''
            character_table = CharacterTable(Phylo.read(StringIO('(dog,((elephant,mouse),robot),cat);'), 'newick'))
            self.assertEqual(2,len(character_table))
            self.assertIn('00111',character_table)
            self.assertIn('00110',character_table)

        def test_mend(self):
            '''MEND Inferring Genotype from a Pedigree'''
            newick_parser = Parser(Tokenizer())
            tree,_        = newick_parser.parse('((((Aa,aa),(Aa,Aa)),((aa,aa),(aa,AA))),Aa);')
            P             = mend(tree)
            self.assertAlmostEqual(0.156, P[0], places=3)
            self.assertAlmostEqual(0.5,   P[1], places=3)
            self.assertAlmostEqual(0.344, P[2], places=3)

        def test_qrt1(self):
            '''qrt Incomplete Characters'''
            quartets = list(qrt(
                                ['cat', 'dog', 'elephant', 'ostrich', 'mouse', 'rabbit', 'robot'],
                                [expand_character(character) for character in [
                                    '01xxx00',
                                    'x11xx00',
                                    '111x00x']]))
            self.assertEqual(4, len(quartets))
            self.assertIn((('cat', 'dog'),  ('mouse', 'rabbit')),quartets)
            self.assertIn((('cat',  'elephant'),  ('mouse', 'rabbit')),quartets)
            self.assertIn((('dog', 'elephant'), ('mouse', 'rabbit')),quartets)
            self.assertIn((( 'dog', 'elephant'), ('rabbit', 'robot')),quartets)

        def test_qrt2(self):
            '''qrt Incomplete Characters'''
            quartets = list(qrt(['Acanthis_azureus',
                                 'Ahaetulla_solitaria',
                                 'Apodora_classicus',
                                 'Dafila_coturnix',
                                 'Eirenis_gallicus',
                                 'Eudrornias_Bernicla',
                                 'Leiocephalus_edulis',
                                 'Leptobrachium_bicoloratum',
                                 'Minipterus_carnifex',
                                 'Paraphysa_ruthveni',
                                 'Pareas_gallinago',
                                 'Pelodytes_nigropalmatus',
                                 'Phrynohyas_leiosoma',
                                 'Rhabdophis_noctua',
                                 'Spermophilus_bobac',
                                 'Terpsihone_completus',
                                 'Thecla_castaneus'],
                                [expand_character(character) for character in [
                                    'x11x0111100x11011',
                                    '01xxxxxxxxxxxxxxx',
                                    'xxxxxxxxxx1xxxx0x',
                                    'x101x0x1x0010101x',
                                    '1010xx0011x000100',
                                    '111x0x1xxx0111xx1',
                                    '1x111x0x11xx11x01',
                                    'xx01001x0xx1x001x',
                                    '111111111x0111011',
                                    'xxx100xxxxxxxxx1x',
                                    '01x1xxx100x1x1x11',
                                    '01x00x0000xxx1x0x',
                                    '0xxxx0x1xxxx0xxxx']]))
        @skip('#46')
        def test_qrtd(self):
            '''qrtd Quartet Distance'''
            self.assertEqual(4,qrtd('A B C D E'.split(),
                                    '(A,C,((B,D),E));',
                                    '(C,(B,D),(A,E));'
            ))

        def test_rsub(self):
            ''' RSUB  	Identifying Reversing Substitutions'''
            substitutions = list(rsub('(((ostrich,cat)rat,mouse)dog,elephant)robot;',
                                ['>robot',
                                 'AATTG',
                                 '>dog',
                                 'GGGCA',
                                 '>mouse',
                                 'AAGAC',
                                 '>rat',
                                 'GTTGT',
                                 '>cat',
                                 'GAGGC',
                                 '>ostrich',
                                 'GTGTC',
                                 '>elephant',
                                 'AATTC']))
            self.assertEqual(5,len(substitutions))
            self.assertIn(('dog', 'rat', 3, 'T', 'G', 'T'),substitutions)
            self.assertIn(('dog', 'mouse', 2, 'A', 'G', 'A'),substitutions)
            self.assertIn(('dog', 'mouse', 1, 'A', 'G', 'A'),substitutions)
            self.assertIn(('rat', 'ostrich', 3, 'G', 'T', 'G'),substitutions)
            self.assertIn(('rat', 'cat', 3, 'G', 'T', 'G'),substitutions)


        def test_sptd(self):
            '''
            SPTD Phylogeny Comparison with Split Distance
            '''
            self.assertEqual(2,
                             sptd('dog rat elephant mouse cat rabbit'.split(),
                                  '(rat,(dog,cat),(rabbit,(elephant,mouse)));',
                                  '(rat,(cat,dog),(elephant,(mouse,rabbit)));'))

        def test_tree1(self):
            '''
            tree -- Completing a Tree
            '''
            self.assertEqual(3,CompleteTree(10,
                                            [(1, 2),
                                             (2, 8),
                                             (4, 10),
                                             (5, 9),
                                             (6, 10),
                                             (7, 9)]))

        def test_tree2(self):
            '''
            tree -- Completing a Tree
            '''
            self.assertEqual(61,CompleteTree(977,[(10, 28),
                                                  (269, 428),
                                                  (109, 878),
                                                  (110, 708),
                                                  (799, 778),
                                                  (375, 223),
                                                  (729, 81),
                                                  (254, 106),
                                                  (317, 310),
                                                  (368, 788),
                                                  (821, 312),
                                                  (596, 56),
                                                  (929, 88),
                                                  (776, 74),
                                                  (162, 41),
                                                  (14, 90),
                                                  (130, 104),
                                                  (333, 735),
                                                  (668, 680),
                                                  (28, 328),
                                                  (18, 461),
                                                  (639, 519),
                                                  (10, 31),
                                                  (452, 253),
                                                  (493, 449),
                                                  (728, 55),
                                                  (360, 479),
                                                  (533, 201),
                                                  (310, 251),
                                                  (235, 98),
                                                  (418, 690),
                                                  (570, 410),
                                                  (541, 136),
                                                  (41, 57),
                                                  (57, 146),
                                                  (759, 42),
                                                  (15, 14),
                                                  (404, 23),
                                                  (831, 303),
                                                  (351, 122),
                                                  (702, 594),
                                                  (661, 739),
                                                  (38, 373),
                                                  (847, 920),
                                                  (549, 120),
                                                  (167, 272),
                                                  (281, 224),
                                                  (381, 41),
                                                  (122, 182),
                                                  (59, 48),
                                                  (397, 614),
                                                  (447, 725),
                                                  (815, 971),
                                                  (677, 468),
                                                  (939, 576),
                                                  (389, 477),
                                                  (95, 141),
                                                  (344, 319),
                                                  (33, 487),
                                                  (99, 176),
                                                  (59, 559),
                                                  (498, 334),
                                                  (881, 21),
                                                  (74, 172),
                                                  (802, 335),
                                                  (619, 421),
                                                  (887, 610),
                                                  (377, 552),
                                                  (270, 16),
                                                  (405, 223),
                                                  (204, 527),
                                                  (204, 184),
                                                  (426, 417),
                                                  (561, 363),
                                                  (764, 438),
                                                  (276, 312),
                                                  (112, 76),
                                                  (624, 206),
                                                  (724, 740),
                                                  (515, 673),
                                                  (323, 113),
                                                  (337, 19),
                                                  (417, 233),
                                                  (20, 12),
                                                  (924, 176),
                                                  (554, 711),
                                                  (561, 698),
                                                  (693, 902),
                                                  (11, 494),
                                                  (294, 301),
                                                  (31, 282),
                                                  (264, 421),
                                                  (443, 431),
                                                  (542, 477),
                                                  (890, 869),
                                                  (251, 525),
                                                  (33, 900),
                                                  (547, 98),
                                                  (857, 352),
                                                  (843, 670),
                                                  (23, 113),
                                                  (540, 248),
                                                  (377, 968),
                                                  (18, 774),
                                                  (128, 64),
                                                  (61, 133),
                                                  (293, 475),
                                                  (453, 219),
                                                  (250, 91),
                                                  (891, 226),
                                                  (810, 41),
                                                  (76, 100),
                                                  (315, 212),
                                                  (318, 915),
                                                  (630, 232),
                                                  (535, 732),
                                                  (328, 341),
                                                  (419, 719),
                                                  (818, 13),
                                                  (319, 38),
                                                  (203, 25),
                                                  (14, 69),
                                                  (156, 36),
                                                  (50, 189),
                                                  (808, 325),
                                                  (45, 114),
                                                  (397, 556),
                                                  (145, 271),
                                                  (43, 108),
                                                  (259, 30),
                                                  (171, 633),
                                                  (365, 347),
                                                  (227, 355),
                                                  (168, 141),
                                                  (748, 341),
                                                  (144, 240),
                                                  (363, 96),
                                                  (364, 610),
                                                  (540, 560),
                                                  (935, 460),
                                                  (547, 566),
                                                  (927, 128),
                                                  (863, 951),
                                                  (712, 160),
                                                  (308, 510),
                                                  (793, 778),
                                                  (34, 38),
                                                  (922, 118),
                                                  (581, 672),
                                                  (63, 61),
                                                  (262, 195),
                                                  (506, 943),
                                                  (811, 638),
                                                  (331, 822),
                                                  (731, 778),
                                                  (110, 260),
                                                  (435, 173),
                                                  (765, 191),
                                                  (380, 26),
                                                  (617, 30),
                                                  (635, 582),
                                                  (218, 176),
                                                  (936, 28),
                                                  (44, 31),
                                                  (482, 492),
                                                  (153, 5),
                                                  (597, 121),
                                                  (56, 238),
                                                  (643, 214),
                                                  (782, 493),
                                                  (440, 675),
                                                  (952, 12),
                                                  (81, 247),
                                                  (269, 133),
                                                  (844, 139),
                                                  (28, 511),
                                                  (109, 488),
                                                  (215, 263),
                                                  (112, 258),
                                                  (573, 83),
                                                  (167, 433),
                                                  (357, 627),
                                                  (645, 259),
                                                  (24, 97),
                                                  (512, 418),
                                                  (80, 11),
                                                  (503, 518),
                                                  (247, 694),
                                                  (359, 465),
                                                  (14, 2),
                                                  (71, 96),
                                                  (32, 302),
                                                  (323, 526),
                                                  (554, 794),
                                                  (91, 72),
                                                  (67, 838),
                                                  (26, 39),
                                                  (973, 297),
                                                  (37, 7),
                                                  (158, 514),
                                                  (471, 640),
                                                  (424, 733),
                                                  (470, 241),
                                                  (580, 730),
                                                  (29, 152),
                                                  (7, 5),
                                                  (451, 918),
                                                  (62, 83),
                                                  (253, 19),
                                                  (853, 398),
                                                  (43, 5),
                                                  (279, 258),
                                                  (534, 247),
                                                  (27, 382),
                                                  (159, 290),
                                                  (609, 547),
                                                  (537, 903),
                                                  (340, 124),
                                                  (291, 305),
                                                  (732, 946),
                                                  (655, 290),
                                                  (189, 190),
                                                  (99, 150),
                                                  (629, 892),
                                                  (613, 306),
                                                  (447, 302),
                                                  (384, 93),
                                                  (544, 278),
                                                  (529, 267),
                                                  (788, 894),
                                                  (291, 942),
                                                  (414, 281),
                                                  (193, 199),
                                                  (703, 71),
                                                  (372, 809),
                                                  (536, 31),
                                                  (615, 591),
                                                  (874, 792),
                                                  (130, 163),
                                                  (206, 46),
                                                  (467, 909),
                                                  (257, 571),
                                                  (977, 167),
                                                  (291, 123),
                                                  (603, 234),
                                                  (278, 87),
                                                  (341, 557),
                                                  (955, 119),
                                                  (226, 187),
                                                  (159, 160),
                                                  (260, 958),
                                                  (110, 813),
                                                  (616, 683),
                                                  (727, 600),
                                                  (626, 415),
                                                  (29, 71),
                                                  (7, 16),
                                                  (35, 391),
                                                  (205, 110),
                                                  (413, 186),
                                                  (941, 896),
                                                  (517, 726),
                                                  (200, 159),
                                                  (667, 627),
                                                  (62, 695),
                                                  (78, 580),
                                                  (287, 411),
                                                  (264, 786),
                                                  (706, 261),
                                                  (304, 297),
                                                  (49, 20),
                                                  (380, 440),
                                                  (144, 111),
                                                  (30, 135),
                                                  (31, 75),
                                                  (329, 819),
                                                  (61, 170),
                                                  (717, 646),
                                                  (238, 537),
                                                  (446, 127),
                                                  (767, 319),
                                                  (654, 701),
                                                  (34, 140),
                                                  (5, 115),
                                                  (962, 734),
                                                  (581, 781),
                                                  (523, 372),
                                                  (112, 286),
                                                  (631, 47),
                                                  (244, 32),
                                                  (112, 476),
                                                  (149, 48),
                                                  (1, 21),
                                                  (183, 187),
                                                  (56, 448),
                                                  (848, 1),
                                                  (305, 393),
                                                  (787, 324),
                                                  (407, 151),
                                                  (360, 861),
                                                  (93, 473),
                                                  (20, 792),
                                                  (491, 245),
                                                  (41, 295),
                                                  (408, 194),
                                                  (416, 125),
                                                  (752, 289),
                                                  (738, 309),
                                                  (576, 222),
                                                  (5, 8),
                                                  (159, 193),
                                                  (456, 509),
                                                  (121, 96),
                                                  (19, 371),
                                                  (870, 913),
                                                  (419, 346),
                                                  (207, 178),
                                                  (266, 548),
                                                  (292, 39),
                                                  (402, 21),
                                                  (230, 111),
                                                  (842, 55),
                                                  (286, 406),
                                                  (203, 251),
                                                  (961, 213),
                                                  (215, 543),
                                                  (56, 127),
                                                  (764, 836),
                                                  (948, 816),
                                                  (161, 503),
                                                  (380, 574),
                                                  (524, 486),
                                                  (716, 416),
                                                  (249, 332),
                                                  (197, 93),
                                                  (464, 241),
                                                  (608, 457),
                                                  (567, 179),
                                                  (17, 1),
                                                  (351, 590),
                                                  (349, 345),
                                                  (10, 23),
                                                  (87, 58),
                                                  (11, 107),
                                                  (975, 440),
                                                  (246, 125),
                                                  (464, 806),
                                                  (456, 259),
                                                  (846, 198),
                                                  (103, 9),
                                                  (19, 180),
                                                  (361, 139),
                                                  (108, 136),
                                                  (217, 357),
                                                  (250, 684),
                                                  (4, 3),
                                                  (756, 111),
                                                  (970, 906),
                                                  (12, 9),
                                                  (41, 15),
                                                  (52, 747),
                                                  (345, 118),
                                                  (944, 824),
                                                  (134, 115),
                                                  (860, 82),
                                                  (99, 47),
                                                  (417, 957),
                                                  (953, 200),
                                                  (339, 761),
                                                  (386, 75),
                                                  (451, 638),
                                                  (572, 744),
                                                  (150, 214),
                                                  (23, 167),
                                                  (227, 273),
                                                  (127, 817),
                                                  (102, 62),
                                                  (785, 375),
                                                  (383, 855),
                                                  (252, 89),
                                                  (270, 480),
                                                  (954, 754),
                                                  (668, 240),
                                                  (35, 21),
                                                  (149, 369),
                                                  (58, 665),
                                                  (522, 507),
                                                  (516, 314),
                                                  (581, 197),
                                                  (484, 800),
                                                  (745, 856),
                                                  (839, 360),
                                                  (383, 288),
                                                  (34, 9),
                                                  (84, 594),
                                                  (720, 240),
                                                  (421, 931),
                                                  (55, 29),
                                                  (101, 60),
                                                  (713, 647),
                                                  (302, 530),
                                                  (326, 459),
                                                  (78, 19),
                                                  (128, 158),
                                                  (297, 124),
                                                  (657, 12),
                                                  (155, 575),
                                                  (592, 305),
                                                  (875, 35),
                                                  (458, 318),
                                                  (33, 20),
                                                  (6, 66),
                                                  (938, 826),
                                                  (115, 823),
                                                  (1, 2),
                                                  (463, 193),
                                                  (762, 96),
                                                  (424, 70),
                                                  (772, 23),
                                                  (129, 415),
                                                  (52, 21),
                                                  (268, 398),
                                                  (497, 406),
                                                  (96, 621),
                                                  (225, 688),
                                                  (60, 58),
                                                  (949, 564),
                                                  (374, 219),
                                                  (515, 108),
                                                  (253, 358),
                                                  (75, 118),
                                                  (550, 106),
                                                  (395, 252),
                                                  (1, 805),
                                                  (910, 592),
                                                  (46, 2),
                                                  (750, 15),
                                                  (514, 945),
                                                  (496, 317),
                                                  (211, 359),
                                                  (737, 239),
                                                  (316, 163),
                                                  (114, 209),
                                                  (733, 889),
                                                  (112, 965),
                                                  (97, 215),
                                                  (55, 352),
                                                  (647, 607),
                                                  (539, 482),
                                                  (45, 148),
                                                  (300, 232),
                                                  (321, 628),
                                                  (140, 432),
                                                  (362, 28),
                                                  (908, 770),
                                                  (916, 280),
                                                  (636, 399),
                                                  (143, 692),
                                                  (796, 411),
                                                  (143, 25),
                                                  (484, 745),
                                                  (689, 73),
                                                  (641, 288),
                                                  (558, 28),
                                                  (194, 42),
                                                  (758, 31),
                                                  (19, 27),
                                                  (307, 49),
                                                  (105, 44),
                                                  (175, 9),
                                                  (233, 343),
                                                  (44, 661),
                                                  (284, 513),
                                                  (106, 306),
                                                  (868, 323),
                                                  (736, 454),
                                                  (44, 124),
                                                  (272, 827),
                                                  (428, 588),
                                                  (921, 901),
                                                  (647, 841),
                                                  (41, 390),
                                                  (427, 106),
                                                  (192, 930),
                                                  (50, 92),
                                                  (709, 202),
                                                  (877, 412),
                                                  (262, 392),
                                                  (658, 372),
                                                  (10, 53),
                                                  (335, 403),
                                                  (178, 116),
                                                  (60, 700),
                                                  (606, 755),
                                                  (847, 505),
                                                  (225, 4),
                                                  (231, 120),
                                                  (387, 676),
                                                  (789, 571),
                                                  (378, 56),
                                                  (16, 394),
                                                  (46, 117),
                                                  (417, 743),
                                                  (508, 474),
                                                  (8, 73),
                                                  (333, 67),
                                                  (331, 125),
                                                  (572, 651),
                                                  (399, 259),
                                                  (274, 283),
                                                  (876, 654),
                                                  (201, 59),
                                                  (330, 265),
                                                  (220, 768),
                                                  (506, 262),
                                                  (466, 276),
                                                  (34, 50),
                                                  (670, 212),
                                                  (593, 622),
                                                  (676, 829),
                                                  (654, 423),
                                                  (178, 401),
                                                  (445, 142),
                                                  (249, 154),
                                                  (489, 458),
                                                  (68, 137),
                                                  (84, 62),
                                                  (263, 618),
                                                  (208, 40),
                                                  (142, 90),
                                                  (187, 659),
                                                  (798, 50),
                                                  (966, 571),
                                                  (139, 780),
                                                  (452, 455),
                                                  (186, 170),
                                                  (542, 820),
                                                  (365, 431),
                                                  (686, 3),
                                                  (335, 438),
                                                  (27, 109),
                                                  (294, 4),
                                                  (1, 74),
                                                  (499, 476),
                                                  (280, 664),
                                                  (886, 778),
                                                  (14, 26),
                                                  (29, 171),
                                                  (210, 245),
                                                  (842, 865),
                                                  (650, 153),
                                                  (48, 30),
                                                  (173, 232),
                                                  (334, 577),
                                                  (791, 438),
                                                  (41, 532),
                                                  (144, 195),
                                                  (685, 153),
                                                  (797, 485),
                                                  (242, 235),
                                                  (545, 151),
                                                  (324, 172),
                                                  (4, 10),
                                                  (385, 751),
                                                  (553, 416),
                                                  (173, 385),
                                                  (372, 742),
                                                  (131, 217),
                                                  (919, 299),
                                                  (531, 262),
                                                  (790, 694),
                                                  (302, 652),
                                                  (181, 485),
                                                  (637, 317),
                                                  (625, 591),
                                                  (5, 2),
                                                  (65, 28),
                                                  (274, 368),
                                                  (511, 582),
                                                  (336, 327),
                                                  (377, 94),
                                                  (959, 495),
                                                  (241, 118),
                                                  (2, 229),
                                                  (830, 790),
                                                  (165, 142),
                                                  (36, 7),
                                                  (662, 576),
                                                  (224, 141),
                                                  (671, 312),
                                                  (220, 101),
                                                  (771, 135),
                                                  (1, 3),
                                                  (8, 13),
                                                  (112, 233),
                                                  (449, 413),
                                                  (257, 9),
                                                  (9, 221),
                                                  (481, 826),
                                                  (562, 51),
                                                  (803, 651),
                                                  (233, 264),
                                                  (357, 563),
                                                  (29, 4),
                                                  (289, 252),
                                                  (743, 828),
                                                  (134, 285),
                                                  (236, 284),
                                                  (882, 335),
                                                  (7, 11),
                                                  (858, 748),
                                                  (76, 9),
                                                  (265, 51),
                                                  (440, 601),
                                                  (620, 872),
                                                  (972, 887),
                                                  (268, 255),
                                                  (36, 68),
                                                  (4, 770),
                                                  (81, 60),
                                                  (88, 34),
                                                  (967, 391),
                                                  (339, 102),
                                                  (395, 591),
                                                  (95, 928),
                                                  (212, 148),
                                                  (486, 478),
                                                  (19, 106),
                                                  (750, 884),
                                                  (324, 517),
                                                  (815, 48),
                                                  (87, 222),
                                                  (633, 899),
                                                  (741, 604),
                                                  (472, 58),
                                                  (933, 636),
                                                  (337, 502),
                                                  (5, 216),
                                                  (425, 233),
                                                  (148, 870),
                                                  (165, 474),
                                                  (318, 292),
                                                  (29, 188),
                                                  (420, 368),
                                                  (607, 185),
                                                  (430, 62),
                                                  (681, 236),
                                                  (204, 423),
                                                  (641, 779),
                                                  (568, 492),
                                                  (141, 663),
                                                  (914, 602),
                                                  (26, 219),
                                                  (30, 32),
                                                  (103, 192),
                                                  (213, 118),
                                                  (682, 368),
                                                  (4, 19),
                                                  (6, 54),
                                                  (459, 833),
                                                  (334, 148),
                                                  (731, 96),
                                                  (478, 189),
                                                  (225, 602),
                                                  (175, 723),
                                                  (862, 129),
                                                  (108, 298),
                                                  (196, 35),
                                                  (72, 871),
                                                  (147, 82),
                                                  (484, 400),
                                                  (167, 507),
                                                  (118, 409),
                                                  (25, 179),
                                                  (926, 769),
                                                  (299, 40),
                                                  (947, 254),
                                                  (325, 816),
                                                  (329, 51),
                                                  (13, 18),
                                                  (905, 770),
                                                  (444, 154),
                                                  (276, 451),
                                                  (328, 612),
                                                  (573, 904),
                                                  (565, 773),
                                                  (312, 495),
                                                  (82, 166),
                                                  (367, 555),
                                                  (34, 629),
                                                  (78, 155),
                                                  (674, 53),
                                                  (932, 357),
                                                  (142, 441),
                                                  (898, 651),
                                                  (276, 22),
                                                  (206, 228),
                                                  (256, 1),
                                                  (52, 467),
                                                  (236, 16),
                                                  (535, 56),
                                                  (925, 422),
                                                  (23, 151),
                                                  (105, 138),
                                                  (182, 754),
                                                  (699, 964),
                                                  (445, 599),
                                                  (485, 584),
                                                  (418, 201),
                                                  (537, 888),
                                                  (13, 72),
                                                  (44, 47),
                                                  (13, 280),
                                                  (300, 976),
                                                  (410, 107),
                                                  (185, 107),
                                                  (159, 121),
                                                  (74, 77),
                                                  (354, 693),
                                                  (440, 734),
                                                  (255, 587),
                                                  (234, 151),
                                                  (828, 885),
                                                  (828, 873),
                                                  (669, 511),
                                                  (25, 94),
                                                  (186, 705),
                                                  (583, 329),
                                                  (606, 288),
                                                  (232, 604),
                                                  (322, 177),
                                                  (126, 60),
                                                  (3, 840),
                                                  (98, 923),
                                                  (806, 940),
                                                  (867, 183),
                                                  (25, 5),
                                                  (287, 103),
                                                  (488, 893),
                                                  (250, 366),
                                                  (41, 327),
                                                  (159, 579),
                                                  (400, 75),
                                                  (360, 283),
                                                  (338, 35),
                                                  (98, 49),
                                                  (835, 170),
                                                  (824, 634),
                                                  (329, 348),
                                                  (361, 501),
                                                  (851, 456),
                                                  (269, 308),
                                                  (623, 216),
                                                  (412, 272),
                                                  (176, 211),
                                                  (134, 145),
                                                  (153, 266),
                                                  (504, 460),
                                                  (783, 34),
                                                  (27, 223),
                                                  (177, 18),
                                                  (600, 795),
                                                  (228, 883),
                                                  (187, 303),
                                                  (832, 653),
                                                  (79, 123),
                                                  (52, 56),
                                                  (642, 46),
                                                  (934, 18),
                                                  (675, 969),
                                                  (662, 879),
                                                  (807, 359),
                                                  (911, 184),
                                                  (956, 528),
                                                  (9, 6),
                                                  (763, 690),
                                                  (678, 287),
                                                  (185, 198),
                                                  (870, 950),
                                                  (293, 26),
                                                  (255, 237),
                                                  (10, 896),
                                                  (226, 313),
                                                  (90, 267),
                                                  (660, 481),
                                                  (25, 120),
                                                  (132, 54),
                                                  (217, 462),
                                                  (40, 33),
                                                  (457, 237),
                                                  (30, 70),
                                                  (284, 595),
                                                  (603, 679),
                                                  (13, 67),
                                                  (469, 70),
                                                  (50, 184),
                                                  (866, 399),
                                                  (129, 119),
                                                  (170, 367),
                                                  (804, 240),
                                                  (852, 456),
                                                  (449, 691),
                                                  (869, 682),
                                                  (554, 444),
                                                  (766, 696),
                                                  (116, 59),
                                                  (82, 32),
                                                  (897, 853),
                                                  (644, 412),
                                                  (434, 201),
                                                  (307, 687),
                                                  (89, 44),
                                                  (715, 657),
                                                  (203, 863),
                                                  (907, 66),
                                                  (58, 64),
                                                  (5, 288),
                                                  (277, 15),
                                                  (291, 326),
                                                  (460, 586),
                                                  (152, 183),
                                                  (69, 161),
                                                  (104, 3),
                                                  (697, 314),
                                                  (845, 296),
                                                  (320, 104),
                                                  (530, 611),
                                                  (435, 634),
                                                  (148, 376),
                                                  (415, 699),
                                                  (209, 718),
                                                  (880, 665),
                                                  (93, 69),
                                                  (279, 859),
                                                  (592, 912),
                                                  (11, 51),
                                                  (572, 87),
                                                  (45, 11),
                                                  (174, 116),
                                                  (756, 895),
                                                  (130, 309),
                                                  (850, 708),
                                                  (468, 36),
                                                  (556, 656),
                                                  (598, 25),
                                                  (917, 557),
                                                  (40, 125),
                                                  (298, 593),
                                                  (114, 154),
                                                  (437, 490),
                                                  (53, 111),
                                                  (653, 571),
                                                  (112, 812),
                                                  (648, 546),
                                                  (546, 361),
                                                  (121, 429),
                                                  (379, 203),
                                                  (906, 29),
                                                  (18, 834),
                                                  (632, 317),
                                                  (85, 12),
                                                  (528, 864),
                                                  (241, 721),
                                                  (814, 418),
                                                  (251, 963),
                                                  (439, 221),
                                                  (27, 396),
                                                  (194, 237),
                                                  (460, 321),
                                                  (350, 81),
                                                  (432, 722),
                                                  (101, 122),
                                                  (441, 589),
                                                  (552, 854),
                                                  (347, 345),
                                                  (42, 26),
                                                  (64, 519),
                                                  (353, 138),
                                                  (360, 757),
                                                  (494, 538),
                                                  (346, 116),
                                                  (119, 181),
                                                  (124, 837),
                                                  (521, 232),
                                                  (616, 487),
                                                  (24, 11),
                                                  (139, 52),
                                                  (157, 72),
                                                  (614, 801),
                                                  (117, 202),
                                                  (296, 241),
                                                  (27, 79),
                                                  (974, 385),
                                                  (74, 86),
                                                  (649, 241),
                                                  (572, 646),
                                                  (324, 746),
                                                  (753, 606),
                                                  (564, 319),
                                                  (270, 335),
                                                  (34, 58),
                                                  (605, 130),
                                                  (109, 173),
                                                  (257, 354),
                                                  (620, 484),
                                                  (144, 749),
                                                  (6, 2),
                                                  (321, 138),
                                                  (696, 481),
                                                  (62, 53),
                                                  (244, 356),
                                                  (422, 53),
                                                  (37, 389),
                                                  (295, 825),
                                                  (415, 769) ]))

    main()
