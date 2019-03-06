#    Copyright (C) 2019 Greenweaves Software Limited
#
#    This is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This software is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>
#
# fragile.py # code for Chapter 6 - Are there fragile regions in the human genome?

from helpers import create_frequency_table
from rosalind import revc,subs

# BA6B Compute the Number of Breakpoints in a Permutation 
#
# getBreakPoints
#
# Parameters:    P list representing a permutation
#
# Returns:  list of pairs of elements that are not adjacent, in extended permutation
#           formed by bracketing P between 0 and len(P) +1.

def getBreakPoints(P):
    return len( [(a,b) for (a,b) in zip([0]+P,P+[len(P)+1]) if b!=a+1] )

# BA6C Compute the 2-Break Distance Between a Pair of Genomes

#  get_synteny_blocks
#
# Input: a genome with circular chromosomes, e.g. 
#        (+1 +2 +3 +4 +5 +6), or
#        (+1 -3 -6 -5)(+2 -4)
#
# Returns: list of synteny blocks, e.g. [1,2,3,4,5,6]

def get_synteny_blocks(a):
    return [j for i in a for j in i] 

#  count_synteny_blocks
#
# Input: a genome with circular chromosomes, e.g. 
#        (+1 +2 +3 +4 +5 +6), or
#        (+1 -3 -6 -5)(+2 -4)
#
# Returns: Number of synteny blocks

def count_synteny_blocks(a):
    return max(abs(b) for b in get_synteny_blocks(a))

# d2break
#
# Find the 2-break distance between two genomes.
#
# Input: Two genomes with circular chromosomes on the same set of synteny blocks.
#
# Return: The 2-break distance between these two genomes.

def d2break(a,b):
    # update   Used to add one edge to the coloured graph
    #
    # Inputs: adjacency    Used to lookup all edges for one node
    #         edges    List of edges
    #         x        Defines start of one edge
    #         y        Defines end of one edge x->y
    #         max_node Hightest node number encountered
    #
    # Returns: Hightest node number encountered
    
    def update(adjacency,edges,x,y,max_node):
        edges.append((x,y))
        
        if x in adjacency:
            adjacency[x].append(y)
        else:
            adjacency[x]=[y] 
            
        if x>max_node:
            max_node=x
        return max_node

    # build_cycle
    #
    # Build one cycle from adjacency list,
    # and delete the edges gthat have been used
    #
    # Inputs:
    #         start     First node of graph
    #         adjacency Adjacency list for noded and edges
    #
    # Returns: newly constructed graph
    #       
    
    def build_cycle(start,adjacency):
        cycle=[]
        ins = [start]
        while len(ins)>0:
            j = ins[0]
            if not j in cycle:
                cycle.append(j)
                for link in adjacency[j]:
                    if not link in cycle:
                        ins.append(link)
            ins.pop(0)
        for i in cycle:
            adjacency.pop(i)
        return cycle
    
    # Verify that the two genomes share the same synteny blocks
    
    blocks = get_synteny_blocks(a)
    n      = count_synteny_blocks(a)
    nb     = count_synteny_blocks(b)
    assert (n==nb),'Mismatched synteny blocks {0} != {1}'.format(n,nb)

    # Combines coloured edges from both chromosomes
    
    edges     = []
    adjacency = {}
    max_node  = -1
    for (x,y) in ColouredEdges(a) + ColouredEdges(b):       
        max_node = update(adjacency,edges,x,y,max_node)
        max_node = update(adjacency,edges,y,x,max_node)

    # Count cycles
    # We'll do this by building each cycle and pruning the ajacency list
    count = 0
    for i in range(1,max_node+1):
        if i in adjacency:
            build_cycle(i,adjacency)
            count +=1
            
    return n - count

# BA6E	Find All Shared k-mers of a Pair of Strings
#
# We say that a k-mer is shared by two genomes if either the k-mer or its
# reverse complement appears in each genome. In Figure 1 are four pairs of
# 3-mers that are shared by "AAACTCATC" and "TTTCAAATC".
#
# A shared k-mer can be represented by an ordered pair (x, y), where x is the
# starting position of the k-mer in the first genome and y is the starting
# position of the k-mer in the second genome. For the genomes "AAACTCATC" and 
# "TTTCAAATC", these shared k-mers are (0,4), (0,0), (4,2), and (6,6).
#
# Input: An integer k and two strings.
#
# Return: All k-mers shared by these strings, in the form of ordered pairs
#         (x, y) corresponding to starting positions of these k-mers in the
#         respective strings.

def find_shared_kmers(k,s1,s2):
    def create_matches():
        freq1   = create_frequency_table(s1,k)
        freq2   = create_frequency_table(s2,k)
        matches = []
        for kmer in freq1:
            if kmer in freq2:
                matches.append((kmer,kmer))
            remk=revc(kmer)
            if remk in freq2:
                matches.append((kmer,remk))
        return matches
    index_pairs=[]
    for k1,k2 in create_matches():
        for i1 in subs(s1,k1):
            for i2 in subs(s2,k2):
                index_pairs.append((i1-1,i2-1))
    return sorted(index_pairs)

# BA6F Implement Chromosome to Cycle
#
# ChromosomeToCycle
#
# Inputs: Chromosome  A chromsome represented as signed synteny blocks e.g. [+1, -2, -3, +4]
#
# Returns: A string representing block x as 2x(head), and 2x-1 (tail). E.g. [1, 2, 4, 3, 6, 5, 7, 8]

def ChromosomeToCycle(Chromosome):
    Nodes = []
    for i in Chromosome:
        if i> 0:
            Nodes.append(2*i-1)
            Nodes.append(2*i)
        else:
            Nodes.append(-2*i)
            Nodes.append(-2*i-1)
        
    return Nodes

# BA6G Implement Cycle to Chromosome
#
# CycleToChromosome
#
# Invert ChromosomeToCycle
#
# Input:  Nodes - fraph in form 2x,2x-1 e.g. [1, 2, 4, 3, 6, 5, 7, 8]
#
# Returns: Chromosome e.g. [+1, -2, -3, +4]

def CycleToChromosome(Nodes):
    Chromosome = []
    it =iter(Nodes)
    for i in it:
        a,b = (i, next(it))
        if a<b:
            Chromosome.append(b//2)
        else:
            Chromosome.append(-a//2)
        
    return Chromosome

# BA6H Implement ColoredEdges

def ColouredEdges(P):
    Edges = []
    for  Chromosome in P:
        Nodes = ChromosomeToCycle(Chromosome)
        it = iter(Nodes[1:]+[Nodes[0]])
        for i in it:
            Edges.append((i,next(it)))
    return Edges

# BA6I Implement Graph to Genome

def GraphToGenome(GenomeGraph):
    def diff(a,b):
        _,x=a
        y,_=b
        return abs(x-y)
    def build_cycle(pair,dcycles):
        result = [pair]
        while pair in dcycles:
            pair = dcycles[pair]
            result.append(pair)
        return result
    def black_edges(cycle):
        result=[]
        for i in range(len(cycle)):
            a,_ = cycle[i]
            _,b = cycle[i-1]
            result.append((b,a))
        return result
 
    extract = [(a,b,diff(a,b)) for (a,b) in zip([GenomeGraph[-1]]+GenomeGraph[0:],GenomeGraph)]
    gaps    = [(a,b) for (a,b,diff) in extract if diff>1]
    cycles  = [(a,b) for (a,b,diff) in extract if diff==1]
    dcycles = {}
    for (a,b) in cycles:
        dcycles[a]=b
    P       = [build_cycle(pair,dcycles) for _,pair in gaps]
    Q       = [black_edges(p) for p in P]
    next_node = 1
    R = []
    for q in Q:
        r = []
        for a,b in q:
            r.append(next_node if a<b else -next_node)
            next_node+=1
        R.append(r)
    return R
