#!/usr/bin/env python

#   (C) 2017-2024 Greenweaves Software Limited

#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

''' Rosalind utilities and simple problems'''

from re import compile
from unittest import main, skip, TestCase
from deprecated import deprecated
import numpy as np
from numpy.testing import assert_array_equal
from scipy.special import comb
from fasta import FastaContent
from reference_tables import CODON_TABLE, BASES, INTEGER_MASSES, AMINO_ACIDS

def iterate_markov(e,p,g):
    '''
    Step through multiple iterations of  Markov chain

    Parameters:
        e         State
        p         Transion matrix
        g         Number of steps
    '''
    for i in range(g):
        e = np.dot(e,p)
    return e

def create_wf_transition_matrix(n):
    '''
    create transition matrix (Feller, page 380)

    Parameters:
        n        Specified dimension of matrix (2*n + 1)
    Returns:
        p[j][k] = probabilty of a transition from j to k,
        where j and k is number of recessives
    '''
    product = np.empty((2*n+1,2*n+1))
    for j in range(2*n+1):
        for k in range(2*n+1):
            product[j,k] =  comb(2*n,k) * (j/(2*n))**k * (1-j/(2*n))**(2*n-k)

    return product

def create_1hot_vector(N,m):
    '''
    Create vector of zeros, except for on specifix element set to zero
    Parameters:
        m
        N
    '''
    e = np.zeros((2*N+1))
    e[2*N-m] =1
    return e

def triplets(dna):
    '''Extract codons (triplets of bases) from a string of DNA'''
    return [dna[i:i+3] for i in range(0,len(dna),3)]

def k_mers(k, bases = list(BASES)):
    '''
    Calculate all possible strings of bases of specified length

    Parameters:
        k          Length of each string
        bases      Specifies characters
    '''
    if k<=0:
        return ['']

    return [ks + b for ks in k_mers(k-1) for b in bases]

def create_frequency_table(string,k):
    '''
    Used to build a table of all the k-mers in a string, with their frequencies
    '''
    frequencies = {}
    for kmer in [string[i:i+k] for i in range(len(string)-k+1)]:
        if kmer not in frequencies:
            frequencies[kmer]=0
        frequencies[kmer]+=1
    return frequencies

def count_subset(s,subset):
    '''
   count_subset

   Count occurences of specific characters

   Parameters:
          s       String in which occurences are to be counted
          subset  String specifying subset of characters to be counted

   Returns:  counts in the same sequence as in subset
    '''
    return [s.count(c) for c in subset]

def grph(fasta,k):
    '''
    GRPH	Overlap Graphs

    A graph whose nodes have all been labeled can be represented by an adjacency list,
    in which each row of the list contains the two node labels corresponding to a unique edge.

    A directed graph (or digraph) is a graph containing directed edges, each of
    which has an orientation. That is, a directed edge is represented by an arrow
    instead of a line segment; the starting and ending nodes of an edge form its
    tail and head, respectively. The directed edge with tail v and head w is
    represented by (v,w) (but not by (w,v)). A directed loop is a directed edge
    of the form (v,v).

    For a collection of strings and a positive integer k, the overlap graph for
    the strings is a directed graph Ok in which each string is represented by a node,
    and string s is connected to string t with a directed edge when there is a
    length k suffix of s that matches a length k prefix of t, as long as sÃ¢â€°Â t;
    we demand sÃ¢â€°Â t to prevent directed loops in the overlap graph
    (although directed cycles may be present).

    Input : A collection of DNA strings in FASTA format having total length at most 10 kbp.

    Return: The adjacency list corresponding to O3. You may return edges in any order.
    '''
    graph=[]
    for name_s,s in fasta:
        for name_t,t in fasta:
            if s!=t and s[-k:]==t[:k]:
                graph.append((name_s,name_t))

    return graph

def revc(dna):
    '''Complementing a Strand of DNA '''
    return dna.translate({
            ord('A'): 'T',
            ord('C'): 'G',
            ord('G'):'C',
            ord('T'): 'A'})[::-1]

def dbru(S,
         include_revc=True):
    '''DBRU Constructing a De Bruijn Graph'''

    def union(S):
        U = set(S)
        if include_revc:
            for s in S:
                U.add(revc(s))
        return U

    def create_nodes(E):
        '''
           Create the nodes of the De Bruijn Graph: these are all kmers present as a substring of a (k+1)mer
           present in S or revc(S)
        '''
        product = set()
        for (a,b) in E:
            product.add(a)
            product.add(b)
        return list(product)

    E = [(e[0:-1],e[1:]) for e in union(S)]
    return (create_nodes(E),E)

@deprecated("I don't think this is actually being used")
def distance(p1,p2):
    return sum([(p1[i]-p2[i])**2 for i in range(len(p1))])


def hamm(s,t):
    '''
    HAMM	Counting Point Mutations
    BA1G	Compute the Hamming Distance Between Two Strings

    Given two strings s and t of equal length, the Hamming distance between s and t,
    denoted dH(s,t), is the number of corresponding symbols that differ in s and t.

    Input: Two DNA strings s and t of equal length (not exceeding 1 kbp).
    Return: The Hamming distance dH(s,t).
    '''
    return len([a for (a,b) in zip(s,t) if a!=b])


class Tree(object):
    '''
    Tree

    This class represents an undirected, weighted tree
    '''

    def __init__(self,
                 N= -1,
                 bidirectional = True):
        '''
        Tree

      Inputs
          N              Number of nodes
          bidirectional
        '''
        self.nodes         = list(range(N))
        self.edges         = {}
        self.bidirectional = bidirectional
        self.N             = N

    def link(self,start,end,weight=1):
        '''
        link

        Link two nodes

            Inputs: start
                    end
                    weight
        '''
        self.half_link(start,end,weight)
        if self.bidirectional:
            self.half_link(end,start,weight)

    def unlink(self,i,k):
        '''
        unlink

        Break link between two nodes
        '''
        try:
            self.half_unlink(i,k)
            if self.bidirectional:
                self.half_unlink(k,i)
        except KeyError:
            print ('Could not unlink {0} from {1}'.format(i,k))
            self.print()

    def half_link(self,a,b,weight=1):
        '''
        half_link

        Create a unidirectional link between two nodes
        '''
        if not a in self.nodes:
            self.nodes.append(a)
        if a in self.edges:
            self.edges[a]=[(b0,w0) for (b0,w0) in self.edges[a] if b0 != b] + [(b,weight)]
        else:
            self.edges[a]=[(b,weight)]

    def half_unlink(self,a,b):
        '''
        half_unlink

        Break link between two nodes
        '''
        links=[(e,w) for (e,w) in self.edges[a] if e != b]
        if len(links)<len(self.edges[a]):
            self.edges[a]=links
        else:
            print ('Could not unlink {0} from {1}'.format(a,b))
            self.print()

    def are_linked(self,a,b):
        '''
        are_linked

        Verify that two nodes are linked
        '''
        return len([e for (e,w) in self.edges[a] if e == b])>0

    def generate_adjacency(self):
        self.nodes.sort()
        for node in self.nodes:
            if node in self.edges:
                for edge in self.edges[node]:
                    end,weight=edge
                    yield node,end,weight

    def print_adjacency(self,includeNodes=False):
        print('-----------------')
        self.nodes.sort()
        if includeNodes:
            print (self.nodes)
        for node in self.nodes:
            if node in self.edges:
                for edge in self.edges[node]:
                    end,weight=edge
                    print ('{0}->{1}:{2}'.format(node,end,weight))



    def traverse(self,i,k,path=[],weights=[]):
        '''
        Traverse a tree, looking for a path between two specified nodes

        Parameters:
            i        Start of path
            k        End of path
            path     Used during recursion to construct path
            weights  Used during recursion to construct weights

        Returns:
            A list of tuples, [(i,0), (n1,w1), ... (k,w)], where the path is [i,n1,...k]
            and the weights [0,w1,...w]
        '''
        if not i in self.edges:
            return (False,[])
        else:
            if len(path) == 0:
                path = [i]
                weights = [0]

            for j,w in self.edges[i]:
                if j in path: continue
                path1 = path + [j]
                weights1 = weights + [w]
                if j==k:
                    return (True,list(zip(path1,weights1)))
                else:
                    found_k,test=self.traverse(j,k,path1,weights1)
                    if found_k:
                        return (found_k,test)
            return (False,[])

    def __len__(self):
        return len(self.nodes)

    def get_nodes(self):
        for node in self.nodes:
            yield(node)

    def initialize_from(self,T):
        for node in T.nodes:
            if not node in self.nodes:
                self.nodes.append(node)
                if node in T.edges:
                    for a,w in T.edges[node]:
                        self.link(node,a,w)

    def get_links(self):
        return [(a,b,w) for a in self.nodes for (b,w) in self.edges[a] if a in self.edges]

    def remove_backward_links(self,root):
        self.bidirectional = False
        for (child,_) in self.edges[root]:
            self.half_unlink(child,root)
            self.remove_backward_links(child)

    def create_topological_order(self):
        '''
        Return a set of all nodes, ordered so descendents prede ancestors
        '''
        def partition_nodes():
            '''
            Partition nodes into leaves (already sorted) and the rest
            '''
            Sorted = []
            ToBeSorted = []
            for node in self.nodes:
                if node in self.edges:
                    ToBeSorted.append(node)
                else:
                    Sorted.append(node)
            return Sorted, ToBeSorted

        def sort_remaining_nodes(Sorted, ToBeSorted):
            '''
            sort_remaining_nodes

            Repeaedly move nodes whose children have been processed
            until there are none left
            '''
            while len(ToBeSorted)>0:
                Unripe = []
                for node in ToBeSorted:
                    if all ([child in Sorted for child,_ in self.edges[node]]):
                        Sorted.append(node)
                    else:
                        Unripe.append(node)
                assert len(Unripe) < len(ToBeSorted)
                ToBeSorted = Unripe
            return Sorted

        Sorted, ToBeSorted = partition_nodes()
        return sort_remaining_nodes(Sorted, ToBeSorted)

class LabeledTree(Tree):
    '''
    LabeledTree

    A labeled tree is a tree in which each vertex is given a unique label.
    '''
    @staticmethod
    def parse(N,lines,
              letters='ATGC',
              bidirectional=True):
        '''
        parse

        Used to create a labeled tree from a list having the following form
            ['4->CAAATCCC',
            '4->ATTGCGAC',
            '5->CTGCGCTG',
            '5->ATGGACGA',
            '6->4',
            '6->5']
        '''
        def get_leaf(string):
            for index,label in T.labels.items():
                if string==label:
                    return index
            return None

        T  = LabeledTree(bidirectional=bidirectional)
        pattern = compile('(([0-9]+)|([{0}]+))->(([0-9]+)|([{0}]+))'.format(letters))

        for line in lines:
            m = pattern.match(line)
            if m:
                i=None

                if m.group(2)==None:
                    i=get_leaf(m.group(3))
                    if i==None:
                        i=len(T.labels)
                    T.set_label(i,m.group(3))

                    if not i in T.leaves:
                        T.leaves.append(i)

                if m.group(3)==None:
                    i=int(m.group(2))

                if m.group(5)==None:
                    j=get_leaf(m.group(6))
                    if j==None:
                        j=len(T.labels)
                        if not j in T.leaves:
                            T.leaves.append(j)
                        T.set_label(j,m.group(6))
                    T.link(i,j)

                if m.group(6)==None:
                    j=int(m.group(5))
                    T.link(i,j)

        return T

    def __init__(self,N=-1,bidirectional=True):
        super().__init__(N,bidirectional=bidirectional)
        self.labels={}
        self.leaves=[]

    def is_leaf(self,v):
        return v in self.leaves

    def set_label(self,node,label):
        if not node in self.nodes:
            self.nodes.append(node)
        self.labels[node]=label

    def initialize_from(self,T):
        super().initialize_from(T)
        for node,label in T.labels.items():
            self.set_label(node,label)


def read_strings(file_name):
    '''
    Read a bunch of strings from file, e.g. reads
    '''
    with open(file_name) as f:
        return [line.strip() for line in f ]

def read_list(file_name):
    '''
    Read a single list of numbers from file
    '''
    with open(file_name) as f:
        return [int(n) for n in f.read().split()]



def write_list(numeric_list,out=None):
    text=' '.join(str(l) for l in numeric_list )
    if out==None:
        print (text)
    else:
        with open(out,'w') as outfile:
            outfile.write(text)

def DPrint(D):
    print ('=====================')
    for drow in D:
        print (', '.join([str(d) for d in drow]))

class RosalindException(Exception):
    '''
    Provide a common class for all exceptions thrown by the library
    '''
    pass


def verify_counts_complete_graph(string):
    '''
    Verify that string contains the same number of As as Us, and the same number of Cs as Gs;
    hence, that it is possible to match all bases.
    '''
    counts={
        'A' : 0,
        'U' : 0,
        'C' : 0,
        'G' : 0
    }
    for c in string:
        counts[c]+=1
    if counts['A']==counts['U'] and counts['C']==counts['G']:
        return counts
    else:
        raise RosalindException('Mismatched counts {0}/{1} or {2}/{3}'.format(counts['A'],counts['U'],counts['C'],counts['G']))



def revc(dna):
    '''
    REVC	Complementing a Strand of DNA
    BA1C	Find the Reverse Complement of a String
           In DNA strings, symbols 'A' and 'T' are complements of each other,
           as are 'C' and 'G'.
           The reverse complement of a DNA string s is the string sc formed by
           reversing the symbols of s, then taking the complement of each symbol
           (e.g., the reverse complement of "GTCA" is "TGAC").

           Input:  A DNA string s of length at most 1000 bp.
           Return: The reverse complement sc of s.
    '''
    def translate(seq,table):
        return ''.join([table[x] for x in seq])
    return translate(dna,{'A': 'T', 'C': 'G', 'G':'C', 'T': 'A'})[::-1]


def subs(s,t):
    '''
    SUBS	Finding a Motif in DNA

    Given two strings s and t, t is a substring of s if t is contained as a
    contiguous collection of symbols in s (as a result, t must be no longer than s).

    The position of a symbol in a string is the total number of symbols found to its left,
    including itself (e.g., the positions of all occurrences of 'U' in "AUGCUUCAGAAAGGUCUUACG"
    are 2, 5, 6, 15, 17, and 18). The symbol at position i of s is denoted by s[i].

    A substring of s can be represented as s[j:k], where j and k represent the starting
    and ending positions of the substring in s; for example,
    if s = "AUGCUUCAGAAAGGUCUUACG", then s[2:5] = "UGCU".

    The location of a substring s[j:k] is its beginning position j; note that t
    will have multiple locations in s if it occurs more than once as a substring of s.

    Input: Two DNA strings s and t (each of length at most 1 kbp).

    Return: All locations of t as a substring of s.
    '''
    matches = []
    start   = 0
    while start>-1:
        start = s.find(t,start)
        if start>-1:
            start += 1
            matches.append(start)
    return matches

def fib(n,k):
    '''FIB	Rabbits and Recurrence Relations'''
    cache=[]
    for i in range(n):
        if i<2:
            cache.append(1)
        else:
            cache.append(cache[i-1]+k*cache[i-2])
    return cache[n-1]



def fibd(n,m):
    '''
    FIBD	Mortal Fibonacci Rabbits

    Input: Positive integers nÃ¢â€°Â¤100 and mÃ¢â€°Â¤20.

    Return: The total number of pairs of rabbits that will remain after the n-th
    month if all rabbits live for m months.
    '''
    def helper(p):
        return [sum([p[i] for i in range(1,len(p))])] +\
               [p[i] for i in range(0,len(p)-1)]

    state=[1]
    while len(state)<m:
        state.append(0)
    for i in range(n-1):
        state=helper(state)
    return sum(state)

def count_nucleotides(s):
    '''
    DNA	  Counting DNA Nucleotides

    Input:  A DNA string s of length at most 1000 nt.
    Return: Four integers (separated by spaces) counting the respective number
            of times that the symbols 'A', 'C', 'G', and 'T' occur in s.
    '''
    return count_subset(s,'ACGT')



def dna_to_rna(dna):
    '''
    RNA	Transcribing DNA into RNA
    An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.
    Given a DNA string t corresponding to a coding strand, its transcribed RNA
    string u is formed by replacing all occurrences of 'T' in t with 'U' in u.

    Input: A DNA string t having length at most 1000 nt.
    Return: The transcribed RNA string of t.
    '''
    return dna.replace('T','U')


def gc(fasta):
    '''
     GC	Computing GC Content

       The GC-content of a DNA string is given by the percentage of symbols in the string
       that are 'C' or 'G'. For example, the GC-content of "AGCTATAG" is 37.5%.
       Note that the reverse complement of any DNA string has the same GC-content.

       Input: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

       Return: The ID of the string having the highest GC-content, followed by
       the GC-content of that string.
    '''
    gcs = [(k,100*float(sum(count_subset(s,'GC')))/len(s)) for (k,s) in fasta]
    return gcs[np.argmax([gc_content for _,gc_content in gcs])]




def sseq(fasta,bases=list(BASES)):
    '''
    SSEQ  Finding a spliced motif
    Input: Two DNA strings s and t (each of length at most 1 kbp) in FASTA format.

    Return: One collection of indices of s in which the symbols of t appear
    as a subsequence of s. If multiple solutions exist, you may return any one.
    '''
    def find_all(text,motif):
        def find_occurenences():
            occurrences = {}
            for base in bases:
                occurrences[base]=[]
            for i in range(len(text)):
                occurrences[text[i]].append(i)
            return occurrences

        occurrences = find_occurenences()
        result = [occurrences[motif[0]][0]]
        for m in motif[1:]:
            for i in range(len(occurrences[m])):
                j = occurrences[m][i]
                if j>result[-1]:
                    result.append(j)
                    break
        return [a + 1 for a in result]

    _,text=fasta[0]
    _,motif=fasta[1]
    return find_all([t for t in text],[m for m in motif])

def lcsm(fasta):
    '''
    LCSM	Finding a Shared Motif

    A common substring of a collection of strings is a substring of every member
    of the collection. We say that a common substring is a longest common substring
    if there does not exist a longer common substring. For example, "CG" is a common
    substring of "ACGTACGT" and "AACCGGTATA", but it is not as long as possible;
    in this case, "GTA" is a longest common substring of "ACGTACGT" and "AACCGTATA".

    Note that the longest common substring is not necessarily unique; for a simple
    example, "AA" and "CC" are both longest common substrings of "AACC" and "CCAA".

    Input: A collection of k DNA strings of length at most 1 kbp each in FASTA format.

    Return: A longest common substring of the collection.
    (If multiple solutions exist, you may return any single solution.)
    '''
    def matches(kmer,strings):
        for s in strings:
            if not kmer in s:
                return False
        return True
    strings=[value for (_,value) in fasta]
    best=['A','C','G','T']
    while len(best[0])<min(len(string) for  string in strings):
        nn=[]
        for s in best:
            for t in ['A','C','G','T']:
                st=s+t
                if matches(st,strings):
                    nn.append(st)
        if len(nn)>0:
            best=nn
        else:
            return best
    return best




def mrna(string,modulus=1000000,table=CODON_TABLE):
    '''
    MRNA	Inferring mRNA from Protein

    Given: A protein string of length at most 1000 aa.

    Return: The total number of different RNA strings from which the protein could
    have been translated, modulo 1,000,000.
    (Don't neglect the importance of the stop codon in protein translation.)
    '''
    nodoc={}
    for key in list(set(table.values())):
        nodoc[key]=0
    for k,v in table.items():
        nodoc[v]+=1
    product=1
    for c in string:
        product*=nodoc[c]
        product=product%modulus
    return (product*nodoc[';'])%modulus


def iprb(k,m,n):
    '''
    IPRB	Mendel's First Law

    Input: Three positive integers k, m, and n, representing a population
           containing k+m+n organisms: k individuals are homozygous dominant for a factor,
           m are heterozygous, and n are homozygous recessive.

    Return: The probability that two randomly selected mating organisms will
    produce an individual possessing a dominant allele (and thus displaying
    the dominant phenotype). Assume that any two organisms can mate.
    '''
    def P_both_offspring_recessive():
        '''
        Calculate probability that both offspring are recessive. There are 4 cases to consider.
        '''
        P1 = n*(n-1)/((k+m+n)*(k+m+n-1))     # Both parents homozygous for recessive
        W1  = 1.0                            # This guarantees offspring recessive
        P2 = n * m /((k+m+n)*(k+m+n-1))      # 1st parent homozygous for recessive, 2nd heterozygous
        W2 =0.5                              # So P(offspring recessive)=0.5
        P3 = m * n /((k+m+n)*(k+m+n-1))      # 1st parent heterozygous, 2nd homozygous for recessive
        W3 = 0.5                             # So P(offspring recessive)=0.5
        P4 = m * (m-1) /((k+m+n)*(k+m+n-1))  # Both heterozygous.
        W4 = 0.25                            # So P(offspring recessive)=0.25
        return  W1 * P1 +  W2 *P2 +  W3 * P3 + W4 * P4

    return 1.0 - P_both_offspring_recessive()



def lia(k,n):
    '''
    LIA 	Independent Alleles

    Input: Two positive integers k (k≤7) and N (N≤2**k). In this problem, we begin
    with Tom, who in the 0th generation has genotype Aa Bb. Tom has two children
    in the 1st generation, each of whom has two children, and so on. Each organism
    always mates with an organism having genotype Aa Bb.

    Return: The probability that at least N Aa Bb organisms will belong to the
    k-th generation of Tom's family tree (don't count the Aa Bb mates at each
    level). Assume that Mendel's second law holds for the factors.
    '''
    transition_probabilities = [[2/4,1/4,0],[2/4,2/4,2/4],[0,1/4,2/4]]
    k_probabilities = [0,1,0]
    for kk in range(k-1):
        new_probabilities=[0,0,0]
        for j in range(3):
            for i in range(3):
                new_probabilities[j] += transition_probabilities[j][i]*k_probabilities[i]
        k_probabilities = new_probabilities

    probability = 0
    prob_individual = k_probabilities[1]**2
    for nn in range(n,2**k+1):
        n1 = 2**k-nn
        probability += comb(2**k,nn) *(1-prob_individual)**n1*prob_individual**nn
    return probability

def rstr(n,x,string):
    '''
    RSTR 	Matching Random Motifs

    Given: A positive integer N<100000, a number x between 0 and 1, and
    a DNA string s of length at most 10 bp.

    Return: The probability that if N random DNA strings having the same length
    as s are constructed with GC-content x (see Introduction to Random Strings),
    then at least one of the strings equals s. We allow for the same random
    string to be created more than once.

    NB "GC-content x" is interpreted as P(G)=P(C)=0.5*x, and similarly for A &T
   '''

    def prob_char(c):
        return 0.5*x if c in ['C','G'] else 0.5*(1-x)
    probability = np.exp(sum([np.log(prob_char(c)) for c in string]))
    return 1- (1-probability)**n

# CONS	Consensus and Profile
#
# Say that we have a collection of DNA strings, all having the same length n.
# Their profile matrix is a  matrix P in which P[1,j] represents the number of
# times that 'A' occurs in the jth position of one of the strings, P[2,j] represents
# the number of times that C occurs in the jth position, and so on.
#
# A consensus string c is a string of length n formed from our collection by
# taking the most common symbol at each position; the jth symbol of c therefore
# corresponds to the symbol having the maximum value in the j-th column of the
# profile matrix. Of course, there may be more than one most common symbol,
# leading to multiple possible consensus strings.
# Input: A collection of at most 10 DNA strings of equal length (at most 1 kbp)
# in FASTA format.
#
# Return: A consensus string and profile matrix for the collection.
#        (If several possible consensus strings exist, then you may return any one of them.)

def cons(fasta):
    (_,string)=fasta[0]
    n=len(string)

    def create_profile():
        product={}
        for c in ['A','C','G','T']:
            product[c] = np.zeros((n))
        for _,string in fasta:
            for i in range(len(string)):
                row=product[string[i]]
                row[i]+=1
        return product

    def create_consensus(profile):
        def most_common_in_column(i):
            best_candidate = ''
            maximum_count = -1
            for candidated in ['A','C','G','T']:
                if profile[candidated][i]>maximum_count:
                    best_candidate = candidated
                    maximum_count = profile[candidated][i]
            return best_candidate

        return ''.join([most_common_in_column(i) for i in range(n)])

    profile=create_profile()
    return (create_consensus(profile),profile)

#IEV	Calculating Expected Offspring
#
# Input: six positive integers, each of which does not exceed 20,000. The
#        integers correspond to the number of couples in a population possessing each
#        genotype pairing for a given factor. In order, the six given integers
#        represent the number of couples having the following genotypes:
#    AA-Aa
#    AA-aa
#    Aa-Aa
#    Aa-aa
#
# Return: The expected number of offspring displaying the dominant phenotype
# in the next generation, under the assumption that every couple has exactly
# two offspring.

def iev(ncopies):
    return 2*sum(n*prob for (n,prob) in zip(ncopies,[1,1,1,0.75,0.5,0]))


def revp(fasta,len1=4,len2=12):
    '''
     REVP	Locating Restriction Sites

     A DNA string is a reverse palindrome if it is equal to its reverse complement.
     For instance, GCATGC is a reverse palindrome because its reverse complement is GCATGC.

     Input: A DNA string of length at most 1 kbp in FASTA format.

     Return: The position and length of every reverse palindrome in the string
             having length between 4 and 12.
             You may return these pairs in any order.
    '''

    def is_palindrome(dna,i,half_length):
        '''Test a substring to see whether it is a palindrome'''
        return (len(dna[i:i+half_length])==half_length and
               dna[i:i+half_length]==revc(dna[i+half_length:i+2*half_length]))

    def find_palindrome(dna,half_length):
        return [(i+1,2*half_length) for i in range(0,len(dna)+1)
                if is_palindrome(dna,i,half_length)]

    def extend(palindromes,half_length):
        return [(i-1,2*half_length) for (i,_) in palindromes
                if is_palindrome(dna,i-2,half_length)]

    (_,dna)     = fasta[0]
    palindromes = find_palindrome(dna,len1//2)
    extension   = palindromes
    for half_length in range(len1//2,len2//2):
        extension = extend(extension,half_length+1)
        if len(extension)==0:
            return palindromes
        else:
            palindromes = palindromes + extension
            latest      = extension

    return palindromes


def random_genome(s,a):
    '''
     PROB 	Introduction to Random Strings

     Input: A DNA string s of length at most 100 bp and an array A containing
            at most 20 numbers between 0 and 1.

     Return: An array B having the same length as A in which B[k] represents the
             common logarithm of the probability that a random string constructed
             with the GC-content found in A[k] will match s exactly.
    '''
    def log_probability(prob_gc):
        log_prob = {
            'G' : np.log10(0.5*prob_gc),
            'C' : np.log10(0.5*prob_gc),
            'A' : np.log10(0.5*(1-prob_gc)),
            'T' : np.log10(0.5*(1-prob_gc))
        }
        return sum([log_prob[ch] for ch in s])
    return [log_probability(prob_gc) for prob_gc in a]



def longestIncreasingSubsequence(N,X):
    '''
    LGIS 	Longest Increasing Subsequence
    A subsequence of a permutation is a collection of elements of the permutation
    in the order in which they appear.
    For example, (5, 3, 4) is a subsequence of (5, 1, 3, 4, 2).

    A subsequence is increasing if the elements of the subsequence increase,
    and decreasing if the elements decrease. For example, given the permutation
    (8, 2, 1, 6, 5, 7, 4, 3, 9), an increasing subsequence is (2, 6, 7, 9),
    and a decreasing subsequence is (8, 6, 5, 4, 3). You may verify that these
    two subsequences are as long as possible.

    Input:
        A positive integer n>=10000 followed by a permutation X of length n.

    Return:
        A longest increasing subsequence of X, followed by a longest
        decreasing subsequence of X.

    Uses algorithm at https://en.wikipedia.org/wiki/Longest_increasing_subsequence
    '''
    def longestMonotoneSubsequence(ascending):
        def ordered(x,y):
            if ascending:
                return x<y
            else:
                return y<x
        P = np.zeros((N),dtype=int)
        M = np.zeros((N+1),dtype=int)
        L = 0
        for i in range(N):
            lo = 1
            hi = L
            while lo <= hi:
                mid = int(np.ceil((lo+hi)/2))
                if ordered(X[M[mid]],X[i]):
                    lo = mid + 1
                else:
                    hi = mid - 1
            newL = lo
            P[i] = M[newL-1]
            M[newL] = i
            if newL > L:
                L = newL

        S = np.zeros((L))
        k = M[L]
        for i in range(L-1,-1,-1):
            S[i] = X[k]
            k = P[k]

        return S

    return (longestMonotoneSubsequence(True),longestMonotoneSubsequence(False))


def get_reading_frames(fasta):
    '''
    ORF 	Open Reading Frames

    Either strand of a DNA double helix can serve as the coding strand for RNA
    transcription. Hence, a given DNA string implies six total reading frames,
    or ways in which the same region of DNA can be translated into amino acids:
    three reading frames result from reading the string itself,
    whereas three more result from reading its reverse complement.

    An open reading frame (ORF) is one which starts from the start codon and ends
    by stop codon, without any other stop codons in between. Thus, a candidate
    protein string is derived by translating an open reading frame into amino
    acids until a stop codon is reached.

    Input: A DNA string s of length at most 1 kbp in FASTA format.

    Return: Every distinct candidate protein string that can be translated from
            ORFs of s. Strings can be returned in any order.

    '''
    def get_start_symbols(peptide,char='M'):
        result=[]
        index=peptide.find(char)
        while index>-1:
            result.append(index)
            index=peptide.find(char,index+1)
        return result

    def read_one_strand(rna):
        peptide=''.join([CODON_TABLE[codon]            \
                         for codon in triplets(rna)     \
                         if len(codon)==3])
        starts=get_start_symbols(peptide)
        ends=get_start_symbols(peptide,';')
        result=[]
        for start in starts:
            for i in range(len(ends)):
                if start<ends[i] and (i==0 or ends[i-1]<start):
                    result.append(peptide[start:ends[i]])
        return result

    def read_one_direction(dna):
        return [val\
                for sublist in [read_one_strand(dna_to_rna(dna)[i:])\
                                for i in range(3)]\
                for val in sublist]

    (_,dna)=fasta.pairs[0]
    return list(set(read_one_direction(dna) + read_one_direction(revc(dna))))

def get_transition_transversion_ratio(fasta):
    '''
    TRAN Transitions and Transversions

    For DNA strings s1 and s2 having the same length, their
    transition/transversion ratio R(s1,s2) is the ratio of the total
    number of transitions to the total number of transversions, where
    symbol substitutions are inferred from mismatched corresponding symbols
    as when calculating Hamming distance (see Counting Point Mutations).

    Input: Two DNA strings s1 and s2 of equal length (at most 1 kbp).

    Return: The transition/transversion ratio R(s1,s2).
    '''
    def is_purine(x):
        '''
        Determine whether a base is a purine

        Parameters:
            x       A DNA base

        Returns:
            True for a purine, False for a pyrimidine
        '''
        return x.upper() in ['A', 'G']

    def is_transition(x,y):
        '''
        Determine whether a mutation is a transition or a transversion

        Parameters:
            x        A DNA base
            y        A mutated DNA base (we assume y != x)

        Returns:
            True if x and y are both purines or both pyrimidines

        '''
        return is_purine(x) == is_purine(y)

    (_,a) = fasta[0]
    (_,b) = fasta[1]
    n_transitions = 0
    n_transversions = 0
    for (x,y) in zip(a,b):
        if x != y:
            if is_transition(x,y):
                n_transitions += 1
            else:
                n_transversions += 1
    return n_transitions/n_transversions


def aspc(n,m):
    '''ASPC 	Introduction to Alternative Splicing'''
    return sum (comb(n,k) for k in range(m,n+1))%1000000


def pper(n,k):
    '''PPER 	Partial Permutations'''
    return n if k==1 else n*pper(n-1,k-1)%1000000

def indc(n):
    '''INDC 	Independent Segregation of Chromosomes'''
    mult=1.0
    for i in range(2*n):
        mult*=0.5
    def p(k):
        return comb(2*n,k)
    def p_cumulative(k):
        return sum(p(kk) for kk in range(k,2*n+1))*mult
    return [np.log10(p_cumulative(k+1)) for k in range(2*n)]


def afrq(ps):
    '''AFRQ 	Counting Disease Carriers'''
    def p_recessive(p):
        return 2*np.sqrt(p)-p
    return [p_recessive(p) for p in ps]



def wfmd(N,m,g,k):
    '''
    WFMD 	The Wright-Fisher Model of Genetic Drift

    Parameters:
        N   Number of diploid individuals in a population
        m   Initial number of copies of a dominant allele
        g   Number of generations
        k   Number of copies of a recessive allele observed after ge generations

    Return:
        The probability that in a population of N diploid individuals
        initially possessing m copies of a dominant allele, we will observe after
        g generations at least k copies of a recessive allele.
    Assume the Wright-Fisher model
    '''
    return iterate_markov(create_1hot_vector(N,m),create_wf_transition_matrix(N),g)[k:].sum()

def ebin(n,P):
    '''EBIN 	Wright-Fisher's Expected Behavior'''
    return [ n*p for p in P]

def foun(N,m,A):
    '''
    FOUN 	The Founder Effect and Genetic Drift
    '''
    def prob(g,a):
        if a==0:
            return 1
        final = iterate_markov(
            create_1hot_vector(N,2*N-a),
            create_wf_transition_matrix(N),
            g)
        return final[0]

    result=[]
    for i in range(m):
        result.append([np.log10(prob(i+1,a))  for a in A])
    return result


def sexl(A):
    '''SEXL 	Sex-Linked Inheritance'''
    return [2*x*(1-x) for x in A]

def sign(n):
    '''SIGN 	Enumerating Oriented Gene Orderings'''
    def expand(p):
        def expanded(bits):
            return [pp if bb==0 else -pp  for pp,bb in zip(p,bits)]
        return [expanded(binary(i,n)) for i in range(2**n)]
    return flatten([expand(p) for p in perm(n)])

def eval (n,s,A):
    '''
    EVAL 	Expected Number of Restriction Sites
    '''
    mult = n - len(s) + 1
    def probability_of_match(a):
        def prob_match(c):
            return 0.5*(a if c=='C' or c=='G' else 1-a)
        return np.exp(sum([np.log(prob_match(c)) for c in s]))
    return [mult*probability_of_match(a) for a in A]

#PERM	Enumerating Gene Orders

def perm(n):
    def perms(symbols):
        if len(symbols)==1:
            return [symbols]
        else:
            result=[]
            for s in symbols:
                r=symbols.copy()
                r.remove(s)
                ps=perms(r)
                for pp in ps:
                    ppp=pp.copy()
                    ppp.insert(0,s)
                    result.append(ppp)
            return result
    return perms([n+1 for n in range(n)])

def lexf(alphabet,k):
    '''
     LEXF	Enumerating k-mers Lexicographically

     Input: A positive integer n

     Return: The total number of permutations of length n, followed by a list of
              all such permutations (in any order)
    '''
    if k<=0:
        return ['']
    else:
        result=[]
        for ks in lexf(alphabet,k-1):
            for l in alphabet.split(' '):
                result.append(ks+l)
    return result



def lexv(alphabet,k):
    '''
    LEXV 	Ordering Strings of Varying Length Lexicographically

    Input: A collection of at most 10 symbols defining an ordered alphabet,
           and a positive integer n (nÃ‚Â¡ÃƒÅ“10).

    Return: All strings of length n that can be formed from the alphabet,
            ordered lexicographically.
    '''

    if k<=0:
        return ['']
    elif k==1:
        return alphabet.split(' ')
    else:
        result=[]
        previous=lexv(alphabet,k-1)
        for letter in alphabet.split(' '):
            result.append(letter)
            for string in previous:
                result.append(letter+string)
    return result


if __name__=='__main__':

    class TestCaseRosalind(TestCase):
        '''Tests for problems that are not from the textbook track'''
        def test_fib(self):
            '''FIB	Rabbits and Recurrence Relations'''
            self.assertEqual(19,fib(5,3))
            self.assertEqual(875089148811941,fib(35, 5))

        def test_fibd(self):
            ''' FIBD	Mortal Fibonacci Rabbits'''
            self.assertEqual(4,fibd(6,3))

        def test_dna(self):
            '''DNA	  Counting DNA Nucleotides'''
            self.assertEqual([20, 12, 17, 21],
                             count_nucleotides(
                                 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGT'
                                 'GGATTAAAAAAAGAGTGTCTGATAGCAGC'))

        def test_rna(self):
            '''RNA	Transcribing DNA into RNA'''
            self.assertEqual('GAUGGAACUUGACUACGUAAAUU',
                             dna_to_rna('GATGGAACTTGACTACGTAAATT'))

        def test_revc(self):
            '''Complementing a Strand of DNA '''
            self.assertEqual('ACCGGGTTTT',revc('AAAACCCGGT'))

        def test_gc(self):
            ''' GC	Computing GC Content'''
            string='''>Rosalind_6404
            CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
            TCCCACTAATAATTCTGAGG
            >Rosalind_5959
            CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
            ATATCCATTTGTCAGCAGACACGC
            >Rosalind_0808
            CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
            TGGGAACCTGCGGGCAGTAGGTGGAAT'''
            fasta = FastaContent(string.split('\n'))
            r     = gc(fasta)
            self.assertEqual('Rosalind_0808',r[0])
            self.assertAlmostEqual(60.919540,r[1],places=5)

        def test_sseq(self):
            string='''>Rosalind_14
            ACGTACGTGACG
            >Rosalind_18
            GTA'''
            fasta = FastaContent(string.split('\n'))
            ss    = sseq(fasta)

            self.assertEqual(3,len(ss))
            self.assertEqual(3,ss[0])
            self.assertEqual(4,ss[1])
            self.assertEqual(5,ss[2])

        def test_hamm(self):
            ''' HAMM	Counting Point Mutations'''
            self.assertEqual(7,hamm('GAGCCTACTAACGGGAT','CATCGTAATGACGGCCT'))

        def test_subs(self):
            '''SUBS	Finding a Motif in DNA'''
            self.assertEqual([2,4,10],subs('GATATATGCATATACTT','ATAT'))

        def test_lcsm(self):
            '''LCSM	Finding a Shared Motif'''
            string='''>Rosalind_1
            GATTACA
            >Rosalind_2
            TAGACCA
            >Rosalind_3
            ATACA'''
            fasta=FastaContent(string.split('\n'))
            self.assertEqual(set(['TA', 'CA', 'AC']),set(lcsm(fasta)))

        def test_mrna(self):
            self.assertEqual(12,mrna("MA"))

        def test_iprb(self):
            '''
            IPRB	Mendel's First Law
            '''
            self.assertAlmostEqual(0.78333,iprb(2,2,2),places=5)

        def test_lia(self): # LIA 	Independent Alleles
            self.assertAlmostEqual(0.684,lia(2,1),places=3)

        def test_rstr(self):
            ''' RSTR 	Matching Random Motifs'''
            self.assertAlmostEqual(0.689,
                                   rstr(90000, 0.6,'ATAGCCGA'),
                                   places=3)

        def test_cons(self):
            string='''>Rosalind_1
            ATCCAGCT
            >Rosalind_2
            GGGCAACT
            >Rosalind_3
            ATGGATCT
            >Rosalind_4
            AAGCAACC
            >Rosalind_5
            TTGGAACT
            >Rosalind_6
            ATGCCATT
            >Rosalind_7
            ATGGCACT'''
            consensus,profile=cons(FastaContent(string.split('\n')))
            self.assertEqual('ATGCAACT',consensus)
            assert_array_equal(profile['T'],[1,5, 0, 0, 0, 1, 1, 6])
            assert_array_equal(profile['G'],[1, 1, 6, 3, 0, 1, 0, 0])
            assert_array_equal(profile['A'],[5, 1, 0, 0, 5, 5, 0, 0])
            assert_array_equal(profile['C'],[0, 0, 1, 4, 2, 0, 6, 1])

        def test_grph(self):
            string='''>Rosalind_0498
            AAATAAA
            >Rosalind_2391
            AAATTTT
            >Rosalind_2323
            TTTTCCC
            >Rosalind_0442
            AAATCCC
            >Rosalind_5013
            GGGTGGG'''
            graph=grph(FastaContent(string.split('\n')),3)
            self.assertEqual(3,len(graph))
            self.assertIn(('Rosalind_0498', 'Rosalind_2391'),graph)
            self.assertIn(('Rosalind_0498', 'Rosalind_0442'),graph)
            self.assertIn(('Rosalind_2391', 'Rosalind_2323'),graph)



        def test_iev(self):
            self.assertEqual(3.5,iev([1,0,0,1,0,1]))

        def test_revp(self):
            string='''>Rosalind_24
            TCAATGCATGCGGGTCTATATGCAT'''
            palindromes=revp(FastaContent(string.split('\n')))
            self.assertIn((4, 6),palindromes)
            self.assertIn((5, 4),palindromes)
            self.assertIn((6, 6),palindromes)
            self.assertIn((7, 4),palindromes)
            self.assertIn((17, 4),palindromes)
            self.assertIn((18, 4),palindromes)
            self.assertIn((20, 6),palindromes)
            self.assertIn((21, 4),palindromes)
            self.assertEqual(8,len(palindromes))



        def test_tran(self):
            '''TRAN Transitions and Transversions'''
            string='''>Rosalind_0209
            GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGA
            AGTACGGGCATCAACCCAGTT
            >Rosalind_2200
            TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGC
            GGTACGAGTGTTCCTTTGGGT'''
            fasta = FastaContent(string.split('\n'))
            self.assertAlmostEqual(1.21428571429,get_transition_transversion_ratio(fasta),places=5)



        def test_aspc(self):
            '''ASPC 	Introduction to Alternative Splicing'''
            self.assertEqual(42,aspc(6,3))


        def test_pper(self):
            ''' PPER 	Partial Permutations'''
            self.assertEqual(51200,pper(21,7))

        def test_indc(self):
            '''INDC 	Independent Segregation of Chromosomes'''
            ii=indc(5)
            self.assertAlmostEqual(0.000,ii[0],3)
            self.assertAlmostEqual(-0.004,ii[1],2)
            self.assertAlmostEqual(-0.024,ii[2],2)
            self.assertAlmostEqual(-0.082,ii[3],3)
            self.assertAlmostEqual(-0.206,ii[4],2)
            self.assertAlmostEqual(-0.424,ii[5],3)
            self.assertAlmostEqual(-0.765,ii[6],3)
            self.assertAlmostEqual(-1.262,ii[7],3)
            self.assertAlmostEqual(-1.969,ii[8],3)
            self.assertAlmostEqual(-3.010,ii[9],3)

        def test_afrq(self):
            '''AFRQ 	Counting Disease Carriers'''
            aa=afrq([0.1, 0.25, 0.5])
            self.assertAlmostEqual(0.532,aa[0],3)
            self.assertAlmostEqual(0.75,aa[1],3)
            self.assertAlmostEqual(0.914,aa[2],3)

        def test_wfmd(self):
            '''
             WFMD 	The Wright-Fisher Model of Genetic Drift
            '''
            self.assertAlmostEqual(0.772,wfmd(4, 6, 2, 1),3)

        def test_ebin(self):
            '''
             EBIN 	Wright-Fisher's Expected Behavior
            '''
            B = ebin(17,[0.1, 0.2, 0.3])
            self.assertAlmostEqual(1.7,B[0],3)
            self.assertAlmostEqual(3.4,B[1],3)
            self.assertAlmostEqual( 5.1,B[2],3)

        def test_foun(self):
            ''' FOUN 	The Founder Effect and Genetic Drift'''
            B = foun(4,3,[0,1,2])
            self.assertAlmostEqual(0.0,              B[0][0],5)
            self.assertAlmostEqual(-0.463935575821,  B[0][1],5)
            self.assertAlmostEqual(-0.999509892866,  B[0][2],5)
            self.assertAlmostEqual(0.0,              B[1][0],5)
            self.assertAlmostEqual(-0.301424998891,  B[1][1],5)
            self.assertAlmostEqual(-0.641668367342,  B[1][2],5)
            self.assertAlmostEqual(0.0,              B[2][0],5)
            self.assertAlmostEqual(-0.229066698008,  B[2][1],5)
            self.assertAlmostEqual(-0.485798552456,  B[2][2],5)

        def test_sexl(self):
            '''SEXL 	Sex-Linked Inheritance'''
            B = sexl([0.1, 0.5, 0.8])
            self.assertAlmostEqual(0.18,B[0],3)
            self.assertAlmostEqual(0.5,B[1],3)
            self.assertAlmostEqual(0.32,B[2],3)

        def test_eval(self):
            '''
             EVAL 	Expected Number of Restriction Sites
            '''
            B = eval(10,'AG',[0.25, 0.5, 0.75])
            self.assertAlmostEqual(0.422,B[0],3)
            self.assertAlmostEqual(0.563,B[1],3)
            self.assertAlmostEqual(0.422,B[2],3)

        def test_longestIncreasingSubsequence1(self):
            '''LGIS 	Longest Increasing Subsequence'''
            (a,d)=longestIncreasingSubsequence(5,[5, 1, 4, 2, 3])
            assert_array_equal([1,2,3],a)
            assert_array_equal([5,4,3],d)

        def test_longestIncreasingSubsequence2(self):
            '''LGIS 	Longest Increasing Subsequence'''
            (a,d)=longestIncreasingSubsequence(9,[8, 2, 1, 6, 5, 7, 4, 3, 9])
            assert_array_equal([1, 5, 7, 9],a)
            assert_array_equal([8, 6, 5, 4, 3],d)

        def test_orf(self):
            string='''>Rosalind_99
            AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'''
            fasta=FastaContent(string.split('\n'))
            peptides=get_reading_frames(fasta)
            self.assertIn('MLLGSFRLIPKETLIQVAGSSPCNLS',peptides)
            self.assertIn('M',peptides)
            self.assertIn('MGMTPRLGLESLLE',peptides)
            self.assertIn('MTPRLGLESLLE',peptides)
            self.assertEqual(4,len(peptides))

        def test_dbru(self):
            '''DBRU Constructing a De Bruijn Graph'''
            _,debruijn= dbru(['TGAT',
                      'CATG',
                      'TCAT',
                      'ATGC',
                      'CATC',
                      'CATC'])
            self.assertEqual(9,len(debruijn))
            self.assertIn(('ATC', 'TCA'),debruijn)
            self.assertIn(('ATG', 'TGA'),debruijn)
            self.assertIn(('ATG', 'TGC'),debruijn)
            self.assertIn(('CAT', 'ATC'),debruijn)
            self.assertIn(('CAT', 'ATG'),debruijn)
            self.assertIn(('GAT', 'ATG'),debruijn)
            self.assertIn(('TCA', 'CAT'),debruijn)
            self.assertIn(('GCA', 'CAT'),debruijn)
            self.assertIn(('TGA', 'GAT'),debruijn)

    main()
