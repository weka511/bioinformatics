'''
 Rosalind utilities

 Copyright (C) 2017 Greenweaves Software Pty Ltd, (c) 2019 Greenweaves Software Limited

 This is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This software is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

 As is the case with point mutations, the most common type of sequencing
 error occurs when a single nucleotide from a read is interpreted incorrectly.

 Given: A collection of up to 1000 reads of equal length (at most 50 bp) in
 FASTA format. Some of these reads were generated with a single-nucleotide error.
 For each read s in the dataset, one of the following applies:
    s was correctly sequenced and appears in the dataset at least twice
	   (possibly as a reverse complement);
    s is incorrect, it appears in the dataset exactly once, and its
	  Hamming distance is 1 with respect to exactly one correct read 
	  in the dataset (or its reverse complement).

 Return: A list of all corrections in the form "[old read]->[new read]". 
 (Each correction must be a single symbol substitution, and you may return the corrections in any order.)
'''

import re,functools
from helpers import translate,count_subset,create_frequency_table,best,triplets,binomial_coefficients,prod,zeroes,k_mers
from fasta import FastaContent
from reference_tables import codon_table,skew_step,bases

def revc(dna):
    return dna.translate({
            ord('A'): 'T',
            ord('C'): 'G',
            ord('G'):'C', 
            ord('T'): 'A'})[::-1]

def dbru(S,include_revc=True):
    '''DBRU Constructing a De Bruijn Graph'''
    def union(S):                    
        U=set(S)
        for s in S:
            s_revc=revc(s)
            if include_revc and not s_revc in U:
                U.add(s_revc)
        return U
    def nodes(E):
        B=[]
        for (a,b) in E:
            if not a in B:
                B.append(a)
            if not b in B:
                B.append(b)
        return B
    E= [(e[0:-1],e[1:]) for e in union(S)]
    return (nodes(E),E)

def distance(p1,p2):
    return sum([(p1[i]-p2[i])**2 for i in range(len(p1))])

# HAMM	Counting Point Mutations
# BA1G	Compute the Hamming Distance Between Two Strings
#
# Given two strings s and t of equal length, the Hamming distance between s and t,
# denoted dH(s,t), is the number of corresponding symbols that differ in s and t.
#
# Input: Two DNA strings s and t of equal length (not exceeding 1 kbp).
# Return: The Hamming distance dH(s,t).

def hamm(s,t):
    return len([a for (a,b) in zip(s,t) if a!=b])

class Tree(object):
    '''
    Undirected, weighted tree
    '''
    def __init__(self,N=-1,bidirectional=True):
        self.nodes=list(range(N))
        self.edges={}
        self.bidirectional=bidirectional
        self.N = N
    def link(self,start,end,weight=1): 
        self.half_link(start,end,weight)
        if self.bidirectional:
            self.half_link(end,start,weight)
        
    def unlink(self,i,k):
        try:
            self.half_unlink(i,k)
            if self.bidirectional:
                self.half_unlink(k,i)
        except KeyError:
            print ('Could not unlink {0} from {1}'.format(i,k))
            self.print()        
        
    def half_link(self,a,b,weight=1):
        if not a in self.nodes:
            self.nodes.append(a)        
        if a in self.edges:               
            self.edges[a]=[(b0,w0) for (b0,w0) in self.edges[a] if b0 != b] + [(b,weight)]
        else:
            self.edges[a]=[(b,weight)]
    
    def half_unlink(self,a,b):
        links=[(e,w) for (e,w) in self.edges[a] if e != b]
        if len(links)<len(self.edges[a]):
            self.edges[a]=links
        else:
            print ('Could not unlink {0} from {1}'.format(a,b))
            self.print()
    
    def are_linked(self,a,b):
        return len([e for (e,w) in self.edges[a] if e == b])>0
        
    def print(self,includeNodes=False):
        print('-----------------')
        self.nodes.sort()
        if includeNodes:
            print (self.nodes)
        for node in self.nodes:
            if node in self.edges:
                for edge in self.edges[node]:
                    end,weight=edge
                    print ('{0}->{1}:{2}'.format(node,end,weight))
                    
    def next_node(self):
        return len(self.nodes)
    
    def traverse(self,i,k,path=[],weights=[]):
        if not i in self.edges: return (False,[])
        if len(path)==0:
            path=[i]
            weights=[0]
        
        for j,w in self.edges[i]:
            if j in path: continue
            path1=path + [j]
            weights1=weights+[w]
            if j==k:
                return (True,list(zip(path1,weights1)))
            else:
                found_k,test=self.traverse(j,k,path1,weights1)
                if found_k: 
                    return (found_k,test)
        return (False,[])  
    
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
            
class LabelledTree(Tree):
        
    @staticmethod
    def parse(N,lines,letters='ATGC',bidirectional=True):
        def get_leaf(string):
            for index,label in T.labels.items():
                if string==label:
                    return index
            return None
        
        T=LabelledTree(bidirectional=bidirectional)
        pattern=re.compile('(([0-9]+)|([{0}]+))->(([0-9]+)|([{0}]+))'.format(letters))
        
        for line in lines:
            m=pattern.match(line)
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

def read_matrix(file_name,conv=int,len_params=1):
    params=[]
    D=[]
    with open(file_name) as f:
        for line in f:
            if len(params)<len_params:
                params.append(int(line.strip()))
            else:
                D.append([conv(s) for s in line.strip().split()])
    return (params,D)

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

#def factorial(n):
#    return reduce(operator.mul,range(1,n+1))

# Verify that string contains the same number of As as Us, and the same number of Cs as Gs;
# hence, that it is possible to match all bases.

def verify_counts_complete_graph(string):
    counts={
        'A':0,
        'U':0,
        'C':0,
        'G':0
    }
    for c in string:
        counts[c]+=1
    if counts['A']==counts['U'] and counts['C']==counts['G']:
        return counts
    else:
        raise RosalindException('Mismatched counts {0}/{1} or {2}/{3}'.format(counts['A'],counts['U'],counts['C'],counts['G']))
    
# REVC	Complementing a Strand of DNA
# BA1C	Find the Reverse Complement of a String
#       In DNA strings, symbols 'A' and 'T' are complements of each other,
#       as are 'C' and 'G'.
#       The reverse complement of a DNA string s is the string sc formed by
#       reversing the symbols of s, then taking the complement of each symbol
#       (e.g., the reverse complement of "GTCA" is "TGAC").
#
#       Input:  A DNA string s of length at most 1000 bp.
#       Return: The reverse complement sc of s.

def revc(dna):
    return translate(dna,{'A': 'T', 'C': 'G', 'G':'C', 'T': 'A'})[::-1]

#SUBS	Finding a Motif in DNA
#
# Given two strings s and t, t is a substring of s if t is contained as a 
# contiguous collection of symbols in s (as a result, t must be no longer than s).
#
# The position of a symbol in a string is the total number of symbols found to its left, 
# including itself (e.g., the positions of all occurrences of 'U' in "AUGCUUCAGAAAGGUCUUACG"
# are 2, 5, 6, 15, 17, and 18). The symbol at position i of s is denoted by s[i].
#
# A substring of s can be represented as s[j:k], where j and k represent the starting
# and ending positions of the substring in s; for example,
# if s = "AUGCUUCAGAAAGGUCUUACG", then s[2:5] = "UGCU".

# The location of a substring s[j:k] is its beginning position j; note that t
# will have multiple locations in s if it occurs more than once as a substring of s.
#
# Input: Two DNA strings s and t (each of length at most 1 kbp).
#
# Return: All locations of t as a substring of s.

def subs(s,t):
    matches=[]
    start=0
    while start>-1:
        start=s.find(t,start)
        if start>-1:
            start+=1
            matches.append(start)
    return matches

#TRIE  Pattern matching

# trie
#
# Given: A list of at most 100 DNA strings of length at most 100 bp, none of which is a prefix of another.
#         one_based Indicates whther numnering of nodes should start at 1 or zero
#
# Return: The adjacency list corresponding to the trie T
#         for these patterns, in the following format. 
#         If T has n nodes, first label the root with 1 and then label the remaining nodes
#         with the integers 2 through n in any order you like.
#         Each edge of the adjacency list of T will be encoded by a triple
#         containing the integer representing the edge's parent node, followed by the integer
#         representing the edge's child node, and finally the symbol labeling the edge.

def trie(strings,one_based=True):
    def find_string_in_adjacency_list(string, adjacency_list):
        index=0
        parent=0
        path=[]
        for cc in string:
            matched=False
            while index<len(adjacency_list) and not matched:
                a,b,c=adjacency_list[index]
                if a==parent and c==cc:
                    matched=True
                    path.append(index)
                    parent=b
                else:
                    index+=1
        return (parent,path)
    
    def create_suffix(parent,path,string,b):
        result=[]
        a = parent
        for i in range(len(path),len(string)):
            b+=1
            result.append((a,b,string[i]))
            a=b

        return result
    
    def merge_string_with_adjacency_list(adjacency_list,string):
        parent,path=find_string_in_adjacency_list(string, adjacency_list)
        return adjacency_list+create_suffix(parent,path,string,len(adjacency_list))
    
    incr = 1 if one_based else 0
    
    def increment_indices(adjacency_list):
        return [(a+incr,b+incr,c) for (a,b,c) in adjacency_list] 
    
    return increment_indices(functools.reduce(merge_string_with_adjacency_list,strings,[]))

# FIB	Rabbits and Recurrence Relations
#
# Input: Positive integers nÃ¢â€°Â¤40 and kÃ¢â€°Â¤5.
#
# Return: The total number of rabbit pairs that will be present after n months
# if we begin with 1 pair and in each generation, every pair of reproduction-age 
# rabbits produces a litter of k rabbit pairs (instead of only 1 pair).
#
# When finding the n-th term of a sequence defined by a recurrence relation, 
# we can simply use the recurrence relation to generate terms for progressively
# larger values of n. This problem introduces us to the computational technique
# of dynamic programming, which successively builds up solutions by using the
# answers to smaller cases.

def fib(n,k):
    cache=[]
    for i in range(n):
        if i<2:
            cache.append(1)
        else:
            cache.append(cache[i-1]+k*cache[i-2])
    return cache[n-1]

# FIBD	Mortal Fibonacci Rabbits
#
# Input: Positive integers nÃ¢â€°Â¤100 and mÃ¢â€°Â¤20.
#
# Return: The total number of pairs of rabbits that will remain after the n-th
# month if all rabbits live for m months.

def fibd(n,m):
    def helper(p):
        return [sum([p[i] for i in range(1,len(p))])] +\
               [p[i] for i in range(0,len(p)-1)]
        
    state=[1]
    while len(state)<m:
        state.append(0)
    for i in range(n-1):
        state=helper(state)
    return sum(state)

# GRPH	Overlap Graphs
#
# A graph whose nodes have all been labeled can be represented by an adjacency list,
# in which each row of the list contains the two node labels corresponding to a unique edge.
#
# A directed graph (or digraph) is a graph containing directed edges, each of
# which has an orientation. That is, a directed edge is represented by an arrow
# instead of a line segment; the starting and ending nodes of an edge form its
# tail and head, respectively. The directed edge with tail v and head w is
# represented by (v,w) (but not by (w,v)). A directed loop is a directed edge
# of the form (v,v).
#
# For a collection of strings and a positive integer k, the overlap graph for
# the strings is a directed graph Ok in which each string is represented by a node,
# and string s is connected to string t with a directed edge when there is a
# length k suffix of s that matches a length k prefix of t, as long as sÃ¢â€°Â t;
# we demand sÃ¢â€°Â t to prevent directed loops in the overlap graph 
# (although directed cycles may be present).
#
# Input : A collection of DNA strings in FASTA format having total length at most 10 kbp.
#
# Return: The adjacency list corresponding to O3. You may return edges in any order.

def grph(fasta,k):
    graph=[]
    for name_s,s in fasta:
        for name_t,t in fasta:
            if s!=t and s[-k:]==t[:k]:
                graph.append((name_s,name_t))
            
    return graph

# SPLC	RNA Splicing
#
# After identifying the exons and introns of an RNA string, we only need to
# delete the introns and concatenate the exons to form a new string
# ready for translation.
#
# Input: A DNA string s (of length at most 1 kbp) and a collection of substrings
#        of s acting as introns. All strings are given in FASTA format.
#
# Return: A protein string resulting from transcribing and translating the 
#         exons of s. (Note: Only one solution will exist for the dataset provided.)

def splc(fasta):
    (_,dna)=fasta[0]
    for i in range(1,len(fasta)):
        (l,intron)=fasta[i]
        fragments=dna.split(intron)
        if len(fragments)>1:
            dna=''.join(fragments)
    return prot(dna_to_rna(dna))
    
# DNA	  Counting DNA Nucleotides
#
# Input:  A DNA string s of length at most 1000 nt.
# Return: Four integers (separated by spaces) counting the respective number
#         of times that the symbols 'A', 'C', 'G', and 'T' occur in s.

def count_nucleotides(s):
    return count_subset(s,'ACGT')

# RNA	Transcribing DNA into RNA
# An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.
# Given a DNA string t corresponding to a coding strand, its transcribed RNA
# string u is formed by replacing all occurrences of 'T' in t with 'U' in u.
#
# Input: A DNA string t having length at most 1000 nt.
# Return: The transcribed RNA string of t.

def dna_to_rna(dna):
    return translate(dna,{'A': 'A', 'C': 'C', 'G':'G', 'T': 'U'})

### Where in the Genome does DNA raplication begin? ###

# BA1A	Compute the Number of Times a Pattern Appears in a Text
#
# We define Count(Text, Pattern) as the number of times that a k-mer Pattern 
# appears as a substring of Text. For example,
#
# Count(ACAACTATGCATACTATCGGGAACTATCCT,ACTAT)=3.
#
# We note that Count(CGATATATCCATAG, ATA) is equal to 3 (not 2) since we should 
# account for overlapping occurrences of Pattern in Text.
#
# To compute Count(Text, Pattern), our plan is to slide a window down Text,
# checking whether each k-mer substring of Text matches Pattern. We will
# therefore refer to the k-mer starting at position i of Text as Text(i, k).
# Throughout this book, we will often use 0-based indexing, meaning that we
# count starting at 0 instead of 1. In this case, Text begins at position 0
# and ends at position |Text| - 1 (|Text| denotes the number of symbols in Text).
# For example, if Text = GACCATACTG, then Text(4, 3) = ATA. Note that the last
# k-mer of Text begins at position |Text| - k, e.g., the last 3-mer of
# GACCATACTG starts at position 10 - 3 = 7.
#
# Input: {DNA strings}} Text and Pattern.
#
# Return: Count(Text, Pattern).

def countOccurrences(pattern,string):
    return len(findOccurences(pattern,string))

# BA1B	Find the Most Frequent Words in a String
#
# We say that Pattern is a most frequent k-mer in Text if it maximizes
# Count(Text, Pattern) among all k-mers. For example, "ACTAT" is a most 
# frequent 5-mer in "ACAACTATGCATCACTATCGGGAACTATCCT", and "ATA" is a most
# frequent 3-mer of "CGATATATCCATAG".
#
# Input: A DNA string Text and an integer k.
#
# Output: All most frequent k-mers in Text (in any order).

def find_most_frequent_words(string,k):
    most_frequent_words=[]
    max_k=-1
    for kmer,count in create_frequency_table(string,k).items():
        if count>max_k:
            max_k=count
            most_frequent_words=[]
        if count==max_k:
            most_frequent_words.append(kmer)
    return most_frequent_words

# BA1D	Find All Occurrences of a Pattern in a String
#
# Input: Strings Pattern and Genome.
#
# Return: All starting positions in Genome where Pattern appears as a substring.
#         Use 0-based indexing.
def findOccurences(pattern,string):
    return [pos-1 for pos in subs(string,pattern)]


# BA1E	Find Patterns Forming Clumps in a String
#
# Given integers L and t, a string Pattern forms an (L, t)-clump inside a
# (larger) string Genome if there is an interval of Genome of length L in which
# Pattern appears at least t times. For example, TGCA forms a (25,3)-clump in
# the following Genome: gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac.
#
# Input: A string Genome, and integers k, L, and t.
#
# Return: All distinct k-mers forming (L, t)-clumps in Genome.

def findClumps(genome,k,L,t):
    def update_patterns(frequencies):
        for (kmer,count) in frequencies.items():
            if count>=t:
                patterns.append(kmer) 
                
    patterns=[]
    frequencies=create_frequency_table(genome[0:L],k)
    update_patterns(frequencies)
    for i in range(1,len(genome)-L+1):
        head=genome[i-1:i+k-1]
        frequencies[head]-=1
        tail=genome[i+L-k:i+L]
        if tail in frequencies:
            frequencies[tail]+=1
        else:
            frequencies[tail]=1        
        update_patterns(frequencies)
        
    return list(set(patterns))

# BA1F	Find a Position in a Genome Minimizing the Skew
#
# Define the skew of a DNA string Genome, denoted Skew(Genome), as the
# difference between the total number of occurrences of 'G' and 'C' in Genome.
#
# Input: A DNA string Genome.
#
# Return: All integer(s) i minimizing Skew(Prefixi (Text)) over all values of
#         i (from 0 to |Genome|).

def find_minimum_skew(genome):
    positions=[]
    min_skew=2
    skew=0
    pos=0
    for nucleotide in genome:
        pos+=1
        skew+=skew_step[nucleotide]
        if min_skew>skew:
            min_skew=skew
            positions=[pos]
        elif min_skew==skew:
            positions.append(pos)

    return positions, min_skew

# BA1H	Find All Approximate Occurrences of a Pattern in a String
#
# Input: Strings Pattern and Text along with an integer d.
#
# Return: All starting positions where Pattern appears as a substring of Text
#         with at most d mismatches.

def findApproximateOccurrences(pattern,text,d):
    return [i\
            for i in range(len(text)-len(pattern)+1)\
            if hamm(pattern,text[i:i+len(pattern)])<=d]

# helper for BA1I
def find_mismatches(pattern,text,d):
    return findApproximateOccurrences(pattern,text,d)

# helper for BA1J
def find_mismatches_and_rc(pattern,text,d):
    return findApproximateOccurrences(pattern,text,d) + \
           findApproximateOccurrences(revc(pattern),text,d)

# BA1I	Find the Most Frequent Words with Mismatches in a String
# BA1J	Find Frequent Words with Mismatches and Reverse Complements
#
# A most frequent k-mer with up to d mismatches in Text is simply a string
# Pattern maximizing Countd(Text, Pattern) among all k-mers. Note that Pattern
# does not need to actually appear as a substring of Text; for example, AAAAA 
# is the most frequent 5-mer with 1 mismatch in AACAAGCTGATAAACATTTAAAGAG,
# even though AAAAA does not appear exactly in this string.
#
# Input: A string Text as well as integers k and d.
#
# Return: All most frequent k-mers with up to d mismatches in Text.

def findMostFrequentWordsWithMismatches(text,k,d,find=find_mismatches):
    matches=[]
    max_count=-1
    for pattern in k_mers(k):
        count=len(find(pattern,text,d))
        if count>max_count:
            max_count=count
            matches=[]
        if count==max_count:
            matches.append(pattern)
            
    return (max_count,matches)

# BA1K	Generate the Frequency Array of a Strings
#
# Given an integer k, we define the frequency array of a string Text as an 
# array of length 4**k, where the i-th element of the array holds the number of
# times that the i-th k-mer (in the lexicographic order) appears in Text.
#
# Input: A DNA string Text and an integer k.
#
# Return: The frequency array of k-mers in Text.

def generateFrequencyArray(text,k):
    frequencies=[]
    for i in range(4**k):
        frequencies.append(0)
    for i in range(len(text)-k+1):
        frequencies[patternToNumber(text[i:i+k])]+=1
    return frequencies

# BA1L	Implement PatternToNumber

def patternToNumber(kmer):
    n=0
    for letter in kmer:
        n*=4
        n+=bases.find(letter)
    return n

# BA1M	Implement NumberToPattern

def numberToPattern(n,k):
    pattern=''
    nn=n
    for i in range(k):
        pattern=pattern+bases[nn%4]
        nn//=4
    return pattern[::-1]

# BA1N	Generate the d-Neighborhood of a String
#
# The d-neighborhood Neighbors(Pattern, d) is the set of all k-mers whose
# Hamming distance from Pattern does not exceed d.
#
# Input: A DNA string Pattern and an integer d.
#
# Return: The collection of strings Neighbors(Pattern, d).

def generate_dNeighborhood(pattern,d):
    def neighbours(p):
        neighbours=[]
        for i in range(len(p)):
            for ch in bases:
                neighbours.append(p[0:i]+ch+p[i+1:])
        return neighbours
    if d==0:
        return [pattern]
    else:
        neighbourhood=[]
        start=generate_dNeighborhood(pattern,d-1)
        for p in start:
            for n in neighbours(p):
                neighbourhood.append(n)
        return list(set(start+neighbourhood))


# GC	Computing GC Content
#
#       The GC-content of a DNA string is given by the percentage of symbols in the string 
#       that are 'C' or 'G'. For example, the GC-content of "AGCTATAG" is 37.5%.
#       Note that the reverse complement of any DNA string has the same GC-content.
#
#       Input: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).
#
#       Return: The ID of the string having the highest GC-content, followed by
#       the GC-content of that string. 

def gc(fasta):
    return best(\
        [(k,100*float(sum(count_subset(s,'GC')))/len(s)) for (k,s) in fasta]\
    )

#SSEQ  Finding a spliced motif
# Input: Two DNA strings s and t (each of length at most 1 kbp) in FASTA format.
# 
# Return: One collection of indices of s in which the symbols of t appear
# as a subsequence of s. If multiple solutions exist, you may return any one.

def sseq(fasta,bases=['A','C','G','T']):
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

# LCSM	Finding a Shared Motif
#
# A common substring of a collection of strings is a substring of every member
# of the collection. We say that a common substring is a longest common substring
# if there does not exist a longer common substring. For example, "CG" is a common
# substring of "ACGTACGT" and "AACCGGTATA", but it is not as long as possible;
# in this case, "GTA" is a longest common substring of "ACGTACGT" and "AACCGTATA".
#
# Note that the longest common substring is not necessarily unique; for a simple
# example, "AA" and "CC" are both longest common substrings of "AACC" and "CCAA".
#
# Input: A collection of k DNA strings of length at most 1 kbp each in FASTA format.
#
# Return: A longest common substring of the collection.
# (If multiple solutions exist, you may return any single solution.)

def lcsm(fasta):
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

def hamm(s,t):
    return len([a for (a,b) in zip(s,t) if a!=b])


# BA4A	Translate an RNA String into an Amino Acid String
# PROT Translating RNA into Protein
#
# The 20 commonly occurring amino acids are abbreviated by using 20 letters from
# the Roman alphabet (all letters except for B, J, O, U, X, and Z). Protein strings
# are constructed from these 20 symbols. Henceforth, the term genetic string will
# incorporate protein strings along with DNA strings and RNA strings.
#
# The RNA codon table dictates the details regarding the encoding of specific
# codons into the amino acid alphabet.
#
# Input: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).
#
# Return: The protein string encoded by s.
#
# NB: I have allowed an extra parameter to deal with alternatives, such as the 
# Mitochundrial codes (http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

def prot(rna,table=codon_table):
    return ''.join([table[codon] for codon in triplets(rna) if table[codon]!=';'])

# MRNA	Inferring mRNA from Protein
#
# Given: A protein string of length at most 1000 aa.
#
# Return: The total number of different RNA strings from which the protein could
# have been translated, modulo 1,000,000.
# (Don't neglect the importance of the stop codon in protein translation.)

def mrna(string,modulus=1000000,table=codon_table):
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

# IPRB	Mendel's First Law
#
# Input: Three positive integers k, m, and n, representing a population
#        containing k+m+n organisms: k individuals are homozygous dominant for a factor,
#        m are heterozygous, and n are homozygous recessive.
#
# Return: The probability that two randomly selected mating organisms will 
# produce an individual possessing a dominant allele (and thus displaying 
# the dominant phenotype). Assume that any two organisms can mate.

def iprb(k,m,n):
    return 1.0 -                                \
           1.0 * n* (n-1)/((k+m+n)*(k+m+n-1)) - \
           0.5 * n * m /((k+m+n)*(k+m+n-1))   - \
           0.5 * m * n /((k+m+n)*(k+m+n-1))   - \
           0.25 * m * (m-1) /((k+m+n)*(k+m+n-1))

# LIA 	Independent Alleles 
#
# Input: Two positive integers k (k≤7) and N (N≤2**k). In this problem, we begin
# with Tom, who in the 0th generation has genotype Aa Bb. Tom has two children
# in the 1st generation, each of whom has two children, and so on. Each organism
# always mates with an organism having genotype Aa Bb.
#
# Return: The probability that at least N Aa Bb organisms will belong to the
# k-th generation of Tom's family tree (don't count the Aa Bb mates at each
# level). Assume that Mendel's second law holds for the factors.

def lia(k,n):
    transition_probabilities=[[2/4,1/4,0],[2/4,2/4,2/4],[0,1/4,2/4]]
    k_probabilities=[0,1,0]
    for kk in range(k-1):
        new_probabilities=[0,0,0]
        for j in range(3):
            for i in range(3):
                new_probabilities[j]+=transition_probabilities[j][i]*k_probabilities[i]
        k_probabilities=new_probabilities
    counts=binomial_coefficients(2**k)
    probability=0
    prob_individual=k_probabilities[1]**2
    for nn in range(n,2**k+1):
        n1=2**k-nn
        probability+=counts[nn]*(1-prob_individual)**n1*prob_individual**nn
    return probability

# RSTR 	Matching Random Motifs
#
# Given: A positive integer N<100000, a number x between 0 and 1, and 
# a DNA string s of length at most 10 bp.
#
# Return: The probability that if N random DNA strings having the same length
# as s are constructed with GC-content x (see Introduction to Random Strings),
# then at least one of the strings equals s. We allow for the same random
# string to be created more than once.
#
# NB "GC-content x" is interpreted as P(G)=P(C)=0.5*x, and similarly for A &T

def rstr(n,x,string):

    def prob_char(c):
        return 0.5*x if c in ['C','G'] else 0.5*(1-x)
    probability=prod([prob_char(c) for c in string])
    return 1-(1-probability)**n

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
            product[c]=zeroes(n)
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
    
if __name__=='__main__':
 
    import unittest
    
    class TestRosalind(unittest.TestCase):
        
        def test_fib(self):
            self.assertEqual(19,fib(5,3))
            self.assertEqual(875089148811941,fib(35, 5))
    
    
    
        def test_fibd(self):
            self.assertEqual(4,fibd(6,3))
            
        
            
        def test_dna(self):
            self.assertEqual([20, 12, 17, 21],\
                             count_nucleotides(\
                                 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGT'\
                                 'GGATTAAAAAAAGAGTGTCTGATAGCAGC'))
        
        def test_rna(self):
            self.assertEqual('GAUGGAACUUGACUACGUAAAUU',dna_to_rna('GATGGAACTTGACTACGTAAATT'))
    
        def test_revc(self):
            self.assertEqual('ACCGGGTTTT',revc('AAAACCCGGT')) 
            
        def test_gc(self):
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
            
            # HAMM	Counting Point Mutations
            # BA1G	Compute the Hamming Distance Between Two Strings
            #
            # Given two strings s and t of equal length, the Hamming distance between s and t,
            # denoted dH(s,t), is the number of corresponding symbols that differ in s and t.
            #
            # Input: Two DNA strings s and t of equal length (not exceeding 1 kbp).
            # Return: The Hamming distance dH(s,t).
            
         
        def test_hamm(self):
            self.assertEqual(7,hamm('GAGCCTACTAACGGGAT','CATCGTAATGACGGCCT'))
            
        def test_subs(self):
            self.assertEqual([2,4,10],subs('GATATATGCATATACTT','ATAT')) 
            
        def test_lcsm(self):
            string='''>Rosalind_1
            GATTACA
            >Rosalind_2
            TAGACCA
            >Rosalind_3
            ATACA'''
            fasta=FastaContent(string.split('\n'))
            self.assertEqual(set(['TA', 'CA', 'AC']),set(lcsm(fasta)))
            
        def test_prot(self):
            self.assertEqual('MAMAPRTEINSTRING',prot('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'))        
          
        def test_mrna(self):
            self.assertEqual(12,mrna("MA"))
            
        def test_iprb(self):
            self.assertAlmostEqual(0.78333,iprb(2,2,2),places=5)  
    
        def test_lia(self): # LIA 	Independent Alleles 
            self.assertAlmostEqual(0.684,lia(2,1),places=3)
            
        def test_rstr(self):
            self.assertAlmostEqual(0.689,rstr(90000, 0.6,'ATAGCCGA'),places=3) # RSTR 	Matching Random Motifs
            
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
            self.assertEqual(profile['T'],[1,5, 0, 0, 0, 1, 1, 6])
            self.assertEqual(profile['G'],[1, 1, 6, 3, 0, 1, 0, 0])
            self.assertEqual(profile['A'],[5, 1, 0, 0, 5, 5, 0, 0])
            self.assertEqual(profile['C'],[0, 0, 1, 4, 2, 0, 6, 1]) 
            
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
            
        def test_splc(self):
            string='''>Rosalind_10
            ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG
            >Rosalind_12
            ATCGGTCGAA
            >Rosalind_15
            ATCGGTCGAGCGTGT'''
            self.assertEqual('MVYIADKQHVASREAYGHMFKVCA',splc(FastaContent(string.split('\n'))))    
            
    ### Where in the Genome does DNA replication begin? ###
    
    class Test_1_Replication(unittest.TestCase):                    
        def test_ba1b(self):
            most_frequent_words=find_most_frequent_words('ACGTTGCATGTCGCATGATGCATGAGAGCT',4)
            self.assertIn('CATG',most_frequent_words)
            self.assertIn('GCAT',most_frequent_words)
            most_frequent_words= find_most_frequent_words('ACATCTGGCGCCCATCGCCCATCTGACGGTTCTGACGGTTCTGACGGTTAGCAGCTCTGACGGTTCTGACGGTTCTGACGGTTCTGACGGTTCTGACGGTTTCGCGACCGCCCATCTGACGGTTCGCCCATACATCTGGTCGCGACCTGACGGTTACATCTGGCGCCCATCTGACGGTTCGCCCATAGCAGCTCTGACGGTTTCGCGACCGCCCATACATCTGGCTGACGGTTTCGCGACTCGCGACACATCTGGAGCAGCTACATCTGGTCGCGACCTGACGGTTAGCAGCTCGCCCATAGCAGCTCTGACGGTTAGCAGCTCTGACGGTTACATCTGGACATCTGGCGCCCATCGCCCATTCGCGACACATCTGGAGCAGCTCGCCCATTCGCGACTCGCGACAGCAGCTACATCTGGTCGCGACACATCTGGCGCCCATCTGACGGTTACATCTGGTCGCGACACATCTGGCGCCCATTCGCGACACATCTGGACATCTGGAGCAGCTACATCTGGTCGCGACAGCAGCTACATCTGGTCGCGACCTGACGGTTACATCTGGCTGACGGTTCGCCCATCGCCCATCGCCCATTCGCGACAGCAGCTTCGCGACCTGACGGTTACATCTGGCGCCCATACATCTGGTCGCGACAGCAGCTTCGCGACTCGCGACCTGACGGTTTCGCGACACATCTGGCTGACGGTTCTGACGGTTCGCCCATAGCAGCTCTGACGGTTCTGACGGTTAGCAGCTACATCTGGAGCAGCTCGCCCATCGCCCATAGCAGCTAGCAGCTTCGCGACACATCTGGCGCCCATCTGACGGTTTCGCGACAGCAGCT',14)
            self.assertIn('CGGTTCTGACGGTT',most_frequent_words)
            self.assertIn('GACGGTTCTGACGG',most_frequent_words)
            self.assertIn('TGACGGTTCTGACG',most_frequent_words)
            self.assertIn('ACGGTTCTGACGGT',most_frequent_words)
            self.assertIn('CTGACGGTTCTGAC',most_frequent_words)
            most_frequent_words=find_most_frequent_words('CGGAAGCGAGATTCGCGTGGCGTGATTCCGGCGGGCGTGGAGAAGCGAGATTCATTCAAGCCGGGAGGCGTGGCGTGGCGTGGCGTGCGGATTCAAGCCGGCGGGCGTGATTCGAGCGGCGGATTCGAGATTCCGGGCGTGCGGGCGTGAAGCGCGTGGAGGAGGCGTGGCGTGCGGGAGGAGAAGCGAGAAGCCGGATTCAAGCAAGCATTCCGGCGGGAGATTCGCGTGGAGGCGTGGAGGCGTGGAGGCGTGCGGCGGGAGATTCAAGCCGGATTCGCGTGGAGAAGCGAGAAGCGCGTGCGGAAGCGAGGAGGAGAAGCATTCGCGTGATTCCGGGAGATTCAAGCATTCGCGTGCGGCGGGAGATTCAAGCGAGGAGGCGTGAAGCAAGCAAGCAAGCGCGTGGCGTGCGGCGGGAGAAGCAAGCGCGTGATTCGAGCGGGCGTGCGGAAGCGAGCGG',12)
            self.assertIn('CGGCGGGAGATT',most_frequent_words)
            self.assertIn('CGGGAGATTCAA',most_frequent_words) 
            self.assertIn('CGTGCGGCGGGA',most_frequent_words)
            self.assertIn('CGTGCGGCGGGA',most_frequent_words)  
            self.assertIn('CGTGGAGGCGTG',most_frequent_words)
            self.assertIn('CGTGGCGTGCGG',most_frequent_words)             
            self.assertIn('GCGTGCGGCGGG',most_frequent_words)
            self.assertIn('GCGTGGAGGCGT',most_frequent_words)  
            self.assertIn('GCGTGGCGTGCG',most_frequent_words)
            self.assertIn('GGAGAAGCGAGA',most_frequent_words)  
            self.assertIn('GGAGATTCAAGC',most_frequent_words)
            self.assertIn('GGCGGGAGATTC',most_frequent_words) 
            self.assertIn('GGGAGATTCAAG',most_frequent_words)
            self.assertIn('GTGCGGCGGGAG',most_frequent_words)
            self.assertIn('TGCGGCGGGAGA',most_frequent_words) 
            
        def test_ba1e(self):
            clumps=findClumps('CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAAC'\
                              'ACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC',5,75,4)
            self.assertIn('CGACA',clumps)
            self.assertIn('GAAGA',clumps)
            self.assertIn('AATGT',clumps)        
            
        def test_ba1f(self):
            positions,_=find_minimum_skew(\
                'CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCT'\
                'TGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG')
            self.assertIn(53,positions)
            self.assertIn(97,positions)        
         
        def test_ba1g(self):
            self.assertEqual(7,hamm('GAGCCTACTAACGGGAT','CATCGTAATGACGGCCT'))
            
        def test_ba1h(self):
            self.assertEqual([6, 7, 26, 27, 78],
                             findApproximateOccurrences('ATTCTGGA',\
                                                        'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC',\
                                                           3)) 
            
        def test_ba1i(self):
            _,matches=findMostFrequentWordsWithMismatches('ACGTTGCATGTCGCATGATGCATGAGAGCT',4,1)
            self.assertIn('GATG',matches)
            self.assertIn('ATGC',matches)
            self.assertIn('ATGT',matches)
               
        def test_ba1k(self):
            self.assertEqual([2,1,0,0,0,0,2,2,1,2,1,0,0,1,1,0],
                             generateFrequencyArray('ACGCGGCTCTGAAA',2))
            
        def test_ba1n(self):
            neighbours=generate_dNeighborhood('ACG',1)
            self.assertIn('CCG', neighbours)
            self.assertIn('TCG', neighbours)
            self.assertIn('GCG', neighbours)
            self.assertIn('AAG', neighbours)
            self.assertIn('ATG', neighbours)
            self.assertIn('AGG', neighbours)
            self.assertIn('ACA', neighbours)
            self.assertIn('ACC', neighbours)
            self.assertIn('ACT', neighbours)
            self.assertIn('ACG', neighbours)        
            
    unittest.main(exit=False) 