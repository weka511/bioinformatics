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

def get_synteny_blocks(a):
    return [j for i in a for j in i] 

def count_synteny_blocks(a):
    return max(abs(b) for b in get_synteny_blocks(a))

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
