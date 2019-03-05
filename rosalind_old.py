# $HeadURL: https://server/svn/sandbox/trunk/rosalind/rosalind.py $
# $LastChangedDate: 2016-10-02 09:36:14 +1300 (Sun, 02 Oct 2016) $
# $LastChangedRevision: 1030 $

# Copyright (C) 2015-2016 Greenweaves Software Pty Ltd

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

# This file contains a collection of functions to solve the problems
# at rosalind.info.

import helpers as rh, reference_tables as rrt, fasta as f, math, sys, numpy,random, functools

class RosalindException(Exception):
    def __init__( self, message ):
        Exception.__init__(self, message)
        
## Rosalind functions

# DNA	  Counting DNA Nucleotides
#
# Input:  A DNA string s of length at most 1000 nt.
# Return: Four integers (separated by spaces) counting the respective number
#         of times that the symbols 'A', 'C', 'G', and 'T' occur in s.

def count_nucleotides(s):
    return rh.count_subset(s,'ACGT')

# RNA	Transcribing DNA into RNA
# An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.
# Given a DNA string t corresponding to a coding strand, its transcribed RNA
# string u is formed by replacing all occurrences of 'T' in t with 'U' in u.
#
# Input: A DNA string t having length at most 1000 nt.
# Return: The transcribed RNA string of t.

def dna_to_rna(dna):
    return rh.translate(dna,{'A': 'A', 'C': 'C', 'G':'G', 'T': 'U'})



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
    return rh.best(\
        [(k,100*float(sum(rh.count_subset(s,'GC')))/len(s)) for (k,s) in fasta]\
    )

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



# MRNA	Inferring mRNA from Protein
#
# Given: A protein string of length at most 1000 aa.
#
# Return: The total number of different RNA strings from which the protein could
# have been translated, modulo 1,000,000.
# (Don't neglect the importance of the stop codon in protein translation.)

def mrna(string,modulus=1000000,table=rrt.codon_table):
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

# LEXF	Enumerating k-mers Lexicographically
#
# Input: A positive integer nÃ‚Â¡ÃƒÅ“7.
# 
# Return: The total number of permutations of length n, followed by a list of
# all such permutations (in any order)

def lexf(alphabet,k):
    if k<=0:
        return ['']
    else:
        result=[]
        for ks in lexf(alphabet,k-1):
            for l in alphabet.split(' '):
                result.append(ks+l)
    return result

# LEXV 	Ordering Strings of Varying Length Lexicographically
#
# Input: A collection of at most 10 symbols defining an ordered alphabet, 
#        and a positive integer n (nÃ‚Â¡ÃƒÅ“10).
#
# Return: All strings of length n that can be formed from the alphabet,
#         ordered lexicographically.

def lexv(alphabet,k):
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
    counts=rh.binomial_coefficients(2**k)
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
    probability=rh.prod([prob_char(c) for c in string])
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
            product[c]=rh.zeroes(n)
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

# REVP	Locating Restriction Sites
#
# A DNA string is a reverse palindrome if it is equal to its reverse complement.
# For instance, GCATGC is a reverse palindrome because its reverse complement is GCATGC.
#
# Input: A DNA string of length at most 1 kbp in FASTA format.
#
# Return: The position and length of every reverse palindrome in the string 
#         having length between 4 and 12.
#         You may return these pairs in any order.

def revp(fasta,len1=4,len2=12):
    
    # Test a substring to see whether it is a palindrome
    def is_palindrome(dna,i,half_length):
        return len(dna[i:i+half_length])==half_length and \
               dna[i:i+half_length]==revc(dna[i+half_length:i+2*half_length])
    
    def find_palindrome(dna,half_length):
        return [(i+1,2*half_length) for i in range(0,len(dna)+1)\
                if is_palindrome(dna,i,half_length)]
    
    def extend(palindromes,half_length):
        return [(i-1,2*half_length) for (i,_) in palindromes\
                if is_palindrome(dna,i-2,half_length)]
     
    (_,dna)=fasta[0]
    palindromes=find_palindrome(dna,len1//2)
    extension=palindromes
    for half_length in range(len1//2,len2//2):
        extension=extend(extension,half_length+1)
        if len(extension)==0:
            return palindromes
        else:
            palindromes=palindromes + extension
            latest=extension
    return palindromes

    
# PROB 	Introduction to Random Strings
#
# Input: A DNA string s of length at most 100 bp and an array A containing
#        at most 20 numbers between 0 and 1.
#
# Return: An array B having the same length as A in which B[k] represents the
#         common logarithm of the probability that a random string constructed
#         with the GC-content found in A[k] will match s exactly.

def random_genome(s,a):
    def log_probability(prob_gc):
        log_prob={
            'G' : math.log10(0.5*prob_gc),
            'C' : math.log10(0.5*prob_gc),
            'A' : math.log10(0.5*(1-prob_gc)),
            'T' : math.log10(0.5*(1-prob_gc))
        }
        return sum([log_prob[ch] for ch in s])
    return [log_probability(prob_gc) for prob_gc in a]

# LGIS 	Longest Increasing Subsequence
# A subsequence of a permutation is a collection of elements of the permutation
# in the order in which they appear.
# For example, (5, 3, 4) is a subsequence of (5, 1, 3, 4, 2).
#
# A subsequence is increasing if the elements of the subsequence increase, 
# and decreasing if the elements decrease. For example, given the permutation
# (8, 2, 1, 6, 5, 7, 4, 3, 9), an increasing subsequence is (2, 6, 7, 9),
# and a decreasing subsequence is (8, 6, 5, 4, 3). You may verify that these 
# two subsequences are as long as possible.
#
# Input: A positive integer n>=10000 followed by a permutation X of length n.
#
# Return: A longest increasing subsequence of X, followed by a longest
#         decreasing subsequence of X.
#
# Uses algorithm at https://en.wikipedia.org/wiki/Longest_increasing_subsequence

def longestIncreasingSubsequence(N,X):
    def longestMonotoneSubsequence(ascending):
        def ordered(x,y):
            if ascending:
                return x<y
            else:
                return y<x
        P = rh.zeroes(N)
        M = rh.zeroes(N+1)
        L = 0
        for i in range(N):
            lo = 1
            hi = L
            while lo<=hi:
                mid = math.ceil((lo+hi)/2)
                if ordered(X[M[mid]],X[i]):
                    lo = mid+1
                else:
                    hi = mid-1                
            newL = lo    
            P[i] = M[newL-1]
            M[newL] = i
            if newL > L:
                L = newL
                
        S = rh.zeroes(L)
        k = M[L]
        for i in range(L-1,-1,-1):
            S[i] = X[k]
            k = P[k]
                
        return S        
    
    return (longestMonotoneSubsequence(True),longestMonotoneSubsequence(False))

# ORF 	Open Reading Frames
#
# Either strand of a DNA double helix can serve as the coding strand for RNA
# transcription. Hence, a given DNA string implies six total reading frames,
# or ways in which the same region of DNA can be translated into amino acids: 
# three reading frames result from reading the string itself,
# whereas three more result from reading its reverse complement.

# An open reading frame (ORF) is one which starts from the start codon and ends
# by stop codon, without any other stop codons in between. Thus, a candidate
# protein string is derived by translating an open reading frame into amino
# acids until a stop codon is reached.
#
# Input: A DNA string s of length at most 1 kbp in FASTA format.
#
# Return: Every distinct candidate protein string that can be translated from 
#         ORFs of s. Strings can be returned in any order.

def get_reading_frames(fasta):
    def get_start_symbols(peptide,char='M'):
        result=[]
        index=peptide.find(char)
        while index>-1:
            result.append(index)
            index=peptide.find(char,index+1)
        return result
    
    def read_one_strand(rna):  
        peptide=''.join([rrt.codon_table[codon]            \
                         for codon in rh.triplets(rna)     \
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

# KMER  Generalizing GC-Content
#
# Input: A DNA string s in FASTA format (having length at most 100 kbp).
#
# Return: The 4-mer composition of s.

def lexig(k,fasta):
    (a,b)=fasta[0]
    counts=rh.zeroes(4**k)
    for index in [patternToNumber(kmer) for kmer in kmer_composition(k,b)]:
        counts[index]+=1
    return counts

#def rear(pairs):
    #def reversal_distance(a,b):
        #return 0 if a==b else 1
    #return [reversal_distance(a,b) for (a,b) in pairs]

#TRAN Transitions and Transversions
#
# For DNA strings s1 and s2 having the same length, their
# transition/transversion ratio R(s1,s2) is the ratio of the total
# number of transitions to the total number of transversions, where
# symbol substitutions are inferred from mismatched corresponding symbols
# as when calculating Hamming distance (see Counting Point Mutations).
#
# Input: Two DNA strings s1 and s2 of equal length (at most 1 kbp).
#
# Return: The transition/transversion ratio R(s1,s2).
def tran(fasta):
    def is_transition(x,y):
        return                      \
               x=='A' and y=='G' or \
               x=='G' and y=='A' or \
               x=='C' and y=='T' or \
               x=='T' and y=='C'
    (_,a)=fasta[0]
    (_,b)=fasta[1]
    n_transitions=0
    n_transversions=0
    for (x,y) in zip(a,b):
        if x!=y:
            if is_transition(x,y):
                n_transitions+=1
            else:
                n_transversions+=1
    return n_transitions/n_transversions

# PDST 	Creating a Distance Matrix 
#
# Input: A collection of nn (n≤10n≤10) DNA strings s1,…,sns1,…,sn of equal 
# length (at most 1 kbp). Strings are given in FASTA format.
#
# Return: The matrix DD corresponding to the p-distance dpdp on the given 
# strings. As always, note that your answer is allowed an absolute error of 0.001.

def distance_matrix(fasta):
    def get_string(i):
        _,string=fasta[i]
        return string
    # For two strings s and t of equal length, the p-distance between them
    # is the proportion of corresponding symbols that differ between s and t.
    def p_distance(s,t):
        return hamm(s,t)/len(s)
    def row(i):
        return [p_distance(get_string(i),get_string(j)) for j in range(len(fasta))]
    return [row(i) for i in range(len(fasta))]

# ASPC 	Introduction to Alternative Splicing 
def aspc(n,m):
    c=rh.create_binomial(n+1)
    return sum (c[rh.binomial_index(n,k)] for k in range(m,n+1))%1000000

# PPER 	Partial Permutations 
def pper(n,k):
    return n if k==1 else n*pper(n-1,k-1)%1000000

#INDC 	Independent Segregation of Chromosomes 
def indc(n):
    c=rh.create_binomial(2*n+1)
    mult=1.0
    for i in range(2*n):
        mult*=0.5
    def p(k):
        return c[rh.binomial_index(2*n,k)]
    def p_cumulative(k):
        return sum(p(kk) for kk in range(k,2*n+1))*mult
    return [math.log10(p_cumulative(k+1)) for k in range(2*n)]

# AFRQ 	Counting Disease Carriers
def afrq(ps):
    def p_recessive(p):
        return 2*math.sqrt(p)-p
    return [p_recessive(p) for p in ps]

# WFMD 	The Wright-Fisher Model of Genetic Drift 
#
# Return:  The probability that in a population of N diploid individuals
# initially possessing m copies of a dominant allele, we will observe after
# g generations at least k copies of a recessive allele.
# Assume the Wright-Fisher model
def wfmd(n,m,g,k):
    
    def accumulate(e):
        return sum([e[j] for j in range(k,2*n+1)])

    return  accumulate(
       rh.iterate_markov(
           rh.create_wf_initial_probabilites(n,m),
           rh.create_wf_transition_matrix(n),g,n))

# EBIN 	Wright-Fisher's Expected Behavior 
def ebin(n,P):
    def expected(p):
        return n*p
    return [expected(p) for p in P]

# FOUN 	The Founder Effect and Genetic Drift 
def foun(N,m,A):
    def prob(g,a):
        if a==0:
            return 1
        final=rh.iterate_markov(
           rh.create_wf_initial_probabilites(N,2*N-a),
           rh.create_wf_transition_matrix(N),
           g,
           N)
        return final[0]
    result=[]
    for i in range(m):
        result.append([math.log10(prob(i+1,a))  for a in A]) 
    return result

# SEXL 	Sex-Linked Inheritance 

def sexl(A):
    return [2*x*(1-x) for x in A]

#SIGN 	Enumerating Oriented Gene Orderings 

def sign(n):
    def expand(p):
        def expanded(bits):
            return [pp if bb==0 else -pp  for pp,bb in zip(p,bits)]
        return [expanded(rh.binary(i,n)) for i in range(2**n)]
    return rh.flatten([expand(p) for p in perm(n)])

# EVAL 	Expected Number of Restriction Sites 

def eval (n,s,A):
    mult=n-len(s)+1
    def probability_of_match(a):
        def prob_match(c):
            return 0.5*(a if c=='C' or c=='G' else 1-a)
        return rh.prod([prob_match(c) for c in s])
    return [mult*probability_of_match(a) for a in A]

#SSEQ  Finding a spliced motif
# Input: Two DNA strings s and t (each of length at most 1 kbp) in FASTA format.
# 
# Return: One collection of indices of s in which the symbols of t appear
# as a subsequence of s. If multiple solutions exist, you may return any one.

def sseq(fasta):
    def find_indices_of_motif(text,motif,start=0,target=None):
        if target==None:
            target=len(motif)
        for i in range(start,len(text)):
            if motif[0]==text[i]:
                if len(motif)==1:
                    return [i]
                else:
                    indices=find_indices_of_motif(text,motif[1:],i+1,target-1)
                    if len(indices)==target-1:
                        return indices + [i]
        return[]
    _,text=fasta[0]
    _,motif=fasta[1]
    return [i+1 for i in find_indices_of_motif(text,motif)[-1::-1]]

#TRIE  Pattern matching

def trie(strings):
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
    
    def increment_indices(adjacency_list):
        return [(a+1,b+1,c) for (a,b,c) in adjacency_list] 
    
    return increment_indices(functools.reduce(merge_string_with_adjacency_list,strings,[]))


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
    for kmer,count in rh.create_frequency_table(string,k).items():
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
    frequencies=rh.create_frequency_table(genome[0:L],k)
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
        skew+=rrt.skew_step[nucleotide]
        if min_skew>skew:
            min_skew=skew
            positions=[pos]
        elif min_skew==skew:
            positions.append(pos)

    return positions, min_skew

def create_skews(genome):
    skews=[]
    skew=0
    for nucleotide in genome:
        skew+=rrt.skew_step[nucleotide]
        skews.append(skew)       
    return skews

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

# helper for BA1I
def find_mismatches(pattern,text,d):
    return findApproximateOccurrences(pattern,text,d)

# helper for BA1J
def find_mismatches_and_rc(pattern,text,d):
    return findApproximateOccurrences(pattern,text,d) + \
           findApproximateOccurrences(revc(pattern),text,d)

def findMostFrequentWordsWithMismatches(text,k,d,find=find_mismatches):
    matches=[]
    max_count=-1
    for pattern in rh.k_mers(k):
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
        n+=rrt.bases.find(letter)
    return n

# BA1M	Implement NumberToPattern

def numberToPattern(n,k):
    pattern=''
    nn=n
    for i in range(k):
        pattern=pattern+rh.bases[nn%4]
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
            for ch in rrt.bases:
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

### Which DNA Patterns play the role of molecular clocks

# BA2A 	Implement MotifEnumeration
#
# Input: Integers k and d, followed by a collection of strings Dna.
#
# Return: All (k, d)-motifs in Dna.

def enumerateMotifs(k,d,dna):
    approximate_matches={}
    kmers= rh.k_mers(k)
    
    def near_matches(kmer):
        if not kmer in approximate_matches:
            approximate_matches[kmer]=[kk for kk in kmers if hamm(kk,kmer)<=d]
        return approximate_matches[kmer]
    def good_enough(pattern):
        def match(string):
            for i in range(len(string)-k+1):
                kmer=string[i:i+k]
                if hamm(kmer,pattern)<=d:
                    return True
            return False
        for string in dna:
            if not match(string):
                return False
        return True
                    
    patterns=[]
    for string in dna:
        for i in range(len(string)-k+1):
            kmer=string[i:i+k]
            for pattern in near_matches(kmer):
                if good_enough(pattern):
                    patterns.append(pattern)
    return ' '.join(sorted(list(set(patterns))))

# BA2B 	Find a Median String 
#
# Input: An integer k and a collection of strings Dna.
#
# Return: A k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers
# Pattern. (If multiple answers exist, you may return any one.)

def medianString(k,dna):
    def findClosest(d): 
        distance=sys.float_info.max
        closest=None
        for k_mer in rh.k_mers(k):
            if distance>d(k_mer,dna):
                distance=d(k_mer,dna)
                closest=k_mer
        return closest
    return findClosest(distanceBetweenPatternAndStrings)


# BA2C 	Find a Profile-most Probable k-mer in a String	
#
# Input: A string Text, an integer k, and a 4 × k matrix Profile.
#
# Return: A Profile-most probable k-mer in Text. 

def mostProbable(text,n,profile):
    # probability of kmer given profile
    def prob(kmer):
        p=1
        for j in range(n):
            i=rrt.bases.find(kmer[j])
            p*=profile[i][j]
        return p
    
    def findMostProbable():
        probability=-1
        best=[]
        for i in range(len(text)-n+1):
            p=prob(text[i:i+n])
            if probability<p:
                probability=p
                best=text[i:i+n]
        return best
    
    return findMostProbable()

# BA2D 	Implement GreedyMotifSearch	
# BA2E 	Implement GreedyMotifSearch with Pseudocounts
#
# Input: Integers k and t, followed by a collection of strings Dna.
#        Optional parameter pseudo_counts specifies whether pseudo counts are to be used
#
# Return: A collection of strings BestMotifs resulting from running 
# GreedyMotifSearch(Dna, k, t). If at any step you find more than one
# Profile-most probable k-mer in a given string, use the one occurring first.
def greedyMotifSearch(k,t,dna,pseudo_counts=False):
    
    # Create an array containing the count of occurences of
    # each base at each position, summed over all motifs
    def count_occurrences_of_bases(
        motifs,\
        initialise_counts=numpy.ones if pseudo_counts else numpy.zeros
        ):
        matrix = initialise_counts((len(rrt.bases),k),dtype=int)
        for kmer in motifs:
            for j in range(k):
                i=rrt.bases.find(kmer[j])
                matrix[i,j]+=1
        return matrix
    
    def profile(motifs):
        return count_occurrences_of_bases(motifs)/float(len(motifs))
    
    def score(motifs):
        matrix=count_occurrences_of_bases(motifs)
        total=0
        for j in range(k):
            m=0
            for i in range(len(rrt.bases)):
                if m<matrix[i,j]:
                    m=matrix[i,j]
            total+=(len(rrt.bases)-m)
        return total
    
    bestMotifs=[genome[0:k] for genome in dna]
    for motif in [dna[0][i:i+k] for i in range(len(dna[0])-k+1)]:
        motifs=[motif]
        for i in range(1,t):
            motifs.append(mostProbable(dna[i],k,profile(motifs)))
        if score(motifs)<score(bestMotifs):
            bestMotifs=motifs
    return bestMotifs	  	  	 

	
	  	  	 
# BA2F 	Implement RandomizedMotifSearch	

def randomized_motif_search(k,t,dna,eps=1):
    def score(k,motifs):
        total=0
        for j in range(k):
            counts=numpy.zeros(len(rrt.bases),dtype=numpy.int32)
            for motif in motifs:
                i=rrt.bases.find(motif[j])
                counts[i]+=1
            max=-1
            ii=-1
            for i in range(len(rrt.bases)):
                if max<counts[i]:
                    ii=i
                    max=counts[ii]
            for i in range(len(rrt.bases)):
                if i!=ii:
                    total+=counts[i]
        return total    
    def counts(motifs):
        matrix=numpy.ones((len(rrt.bases),k),dtype=int)
        for i in range(len(rrt.bases)):
            for j in range(k):
                matrix[i,j]*=eps
        for kmer in motifs:
            for j in range(k):
                i=rrt.bases.find(kmer[j])
                matrix[i,j]+=1
        return matrix
    def Motifs(profile,dna):
        def get_k(Profile):
            return len(Profile[0])        
        def prob(kmer):
            p=1
            for j in range(k):
                i=rrt.bases.find(kmer[j])
                p*=profile[i][j]
            return p    
        k=get_k(profile)
        motifs=[]
        for s in dna:
            max_probability=-1
            most_probable_kmer=''
            for kmer in [s[i:i+k].upper() for i in range(len(s)-k+1)]:
                #if len(kmer)<k:
                    #print kmer,s,i
                if max_probability<prob(kmer):
                    max_probability=prob(kmer)
                    most_probable_kmer=kmer
            motifs.append(most_probable_kmer)
        return motifs    
    def Profile(motifs):
        matrix=counts(motifs)
        probabilities=numpy.zeros((len(rrt.bases),k),dtype=float)
        for i in range(len(rrt.bases)):
            for j in range(k):
                probabilities[i,j]=matrix[i,j]/float(len(motifs))
        return probabilities
  
    def random_kmer(string):
        i=random.randint(0,len(string)-k)
        return string[i:i+k]
 
    motifs=[]
    
    for i in range(t):
        motifs.append(random_kmer(dna[i]))
    bestMotifs=motifs
    while True:
        profile=Profile(motifs)
        motifs = Motifs(profile, dna)
        if score(k,motifs) < score(k,bestMotifs):
            bestMotifs = motifs
        else:
            return (score(k,bestMotifs),bestMotifs)

def randomized_motif_search_driver(k,t,dna,N=1000):
    best=sys.float_info.max
    mm=[]
    for i in range(N):
        (sc,motifs) =randomized_motif_search(k,t,dna)
        if sc<best:
            best=sc
            mm=motifs
        if i%100==0:
            print (i,best)
            for m in mm:
                print (m)
    return (best,mm)

	  	  	 
# BA2G 	Implement GibbsSampler	
#
# Input: Integers k, t, and N, followed by a collection of strings Dna.
#
# Return: The strings BestMotifs resulting from running 
# GibbsSampler(Dna, k, t, N) with 20 random starts.
# Remember to use pseudocounts!

def gibbs(k,t,n,dna,eps=1):
    def score(k,motifs):
        total=0
        for j in range(k):
            counts=numpy.zeros(len(rrt.bases),dtype=numpy.int32)
            for motif in motifs:
                i=rrt.bases.find(motif[j])
                counts[i]+=1
            max=-1
            ii=-1
            for i in range(len(rrt.bases)):
                if max<counts[i]:
                    ii=i
                    max=counts[ii]
            for i in range(len(rrt.bases)):
                if i!=ii:
                    total+=counts[i]
        return total    
    def random_kmer(string):
        i=random.randint(0,len(string)-k)
        return string[i:i+k]
    def dropOneMotif(motifs,i):
        return [motifs[j] for j in range(len(motifs)) if j!=i]
    def counts(motifs):
        matrix=numpy.ones((len(rrt.bases),k),dtype=int)
        for i in range(len(rrt.bases)):
            for j in range(k):
                matrix[i,j]*=eps
        for kmer in motifs:
            for j in range(k):
                i=rrt.bases.find(kmer[j])
                matrix[i,j]+=1
        return matrix    
    def Profile(motifs):
        matrix=counts(motifs)
        probabilities=numpy.zeros((len(rrt.bases),k),dtype=float)
        for i in range(len(rrt.bases)):
            for j in range(k):
                probabilities[i,j]=matrix[i,j]/float(len(motifs))
        return probabilities
    
    def probability(kmer,profile):
        p=1
        for j in range(len(kmer)):
            i=rrt.bases.find(kmer[j])
            p*=profile[i][j]
        return p
    
    def accumulate(probabilities):
        total=0
        cumulative=[]
        for p in probabilities:
            total+=p
            cumulative.append(total)
        return cumulative
    
    def generate(probabilities):
        accumulated=accumulate(probabilities)
        rr=accumulated[len(accumulated)-1]*random.random()
        i=0
        while accumulated[i]<=rr:
            i+=1
        return i
    
    motifs=[]
    
    for i in range(t):
        motifs.append(random_kmer(dna[i]))
    bestMotifs=motifs
    
    trace=[]
    best_score=sys.float_info.max
    for j in range(n):
        i=random.randint(0,t-1)
        profile=Profile(dropOneMotif(motifs,i))
        motif_index=generate([probability(dna[i][ll:ll+k],profile)\
                              for ll in range(len(dna[i])-k)])
        motifs[i]=dna[i][motif_index:motif_index+k]
        sc=score(k,motifs)
        if  sc< best_score:
            best_score=sc
            bestMotifs = motifs
        trace.append(best_score)
                   
    return (score(k,bestMotifs),bestMotifs,trace)

	  	  	 
# BA2H 	Implement DistanceBetweenPatternAndStrings 

def distanceBetweenPatternAndStrings (pattern,dna):
    # Extend Hamming distance to work with string of unequal length
    def hamming(pattern,genome):
        return min([hamm(pattern,genome[i:i+len(pattern)]) 
                    for i in range(len(genome)-len(pattern)+1)])        
    return sum([hamming(pattern,motif) for motif in dna])
    
### How do we assemble genomes?

# BA3A	Generate the k-mer Composition of a String
#
# Input: An integer k and a string Text.
#
# Return: Compositionk(Text) (the k-mers can be provided in any order).

def kmer_composition(k,dna):
    return [dna[i:i+k] for i in range(1+len(dna)-k)]   

# BA3B	Reconstruct a String from its Genome Path
#
# Input: A sequence of k-mers Pattern1, ... , Patternn such that the last k - 1
# symbols of Pattern[i] are equal to the first k - 1 symbols of Pattern[i+1]
# for i from 1 to n-1.
#
# Return: A string Text of length k+n-1 where the i-th k-mer in Text is equal
# to Pattern[i] for all i.

def reconstruct_as_list(k,n,fragments,extract=lambda fragments,i: fragments[i]):
    result=[]   
    for i in range(0,n,k):
        result.append(extract(fragments,i))
    target_len=n+k-1
    actual_len=len(result)*k
    overlap=target_len-actual_len
    if overlap>0:
        result.append(fragments[-1][k-overlap:k])
    return result

def reconstruct(fragments):
    return ''.join(reconstruct_as_list(len(fragments[0]),len(fragments),fragments))

# BA3C	Construct the Overlap Graph of a Collection of k-mers
#
# Construct the overlap graph of a collection of k-mers.
#
# Input: A collection Patterns of k-mers.
#
# Return: The overlap graph Overlap(Patterns), in the form of an adjacency list.

def grph_kmers(strings):
    kk=len(strings[0])-1
    graph=[]
    for s in strings:
        for t in strings:
            if s!=t and s[-kk:]==t[:kk]:
                graph.append((s,t))
            
    return graph

# BA3D 	Construct the De Bruijn Graph of a String 
#
# Given: An integer k and a string Text.
#
# Return:DeBruijnk(Text), in the form of an adjacency list.

def deBruijn(k,text):
    kmers=kmer_composition(k-1,text)
    def standardize(ll):
        lll=list(set(ll))
        lll.sort()
        return lll
    
    pathgraph=grph_kmers(kmers)

    deBruijn_dict={}
    for [a,b] in pathgraph:
        if not a in deBruijn_dict:
            deBruijn_dict[a]=[]
        deBruijn_dict[a].append(b)
    graph= [(a,standardize(deBruijn_dict[a])) for a in deBruijn_dict.keys()]
    graph.sort()
    return graph

#BA3E 	Construct the De Bruijn Graph of a Collection of k-mers
#
# Input: A collection of k-mers Patterns.
#
# Return: The de Bruijn graph DeBruijn(Patterns), in the form of an adjacency list.

def deBruijn_collection(pattern,\
                        head=lambda kmer: kmer[0:-1],
                        tail=lambda kmer: kmer[1:]):
    graph={}
    k=len(pattern[0])
    for kmer in pattern:
        if not head(kmer) in graph:
            graph[head(kmer)]=[]
        graph[head(kmer)].append(tail(kmer))
    for kmer in graph.keys():
        graph[kmer].sort()
    return graph
 
# BA3F 	Find an Eulerian Cycle in a Graph 
#
# Input: An Eulerian directed graph, in the form of an adjacency list.
#
# Return: An Eulerian cycle in this graph.

def find_eulerian_cycle(graph):
    def create_unexplored_edges():
        edges=[]
        for a in graph.keys():
            for b in graph[a]:
                edges.append((a,b))
        return edges
    
    def find_next(node):
        for succ in graph[node]:
            if (node,succ) in unexplored:
                return succ
        return None
    
    def find_branch(cycle):
        pos=0
        for node in cycle:
            succ=find_next(node)
            if succ!=None:
                return (pos,node)
            pos+=1
            
    def find_cycle(cycle0,node):
        cycle=cycle0
        succ=find_next(node)
        while succ!=None:
            cycle.append(succ)
            unexplored.remove((node,succ))
            node=succ
            succ=find_next(node)
        return cycle
    
    unexplored= create_unexplored_edges()
    target_edges=len(unexplored)
    node,_=unexplored[0]
    cycle=find_cycle([node],node)
    while len(unexplored)>0:
        (pos,branch)=find_branch(cycle) 
        cycle=rh.rotate(cycle,pos)
        cycle=find_cycle(cycle,cycle[len(cycle)-1])
    return cycle

# BA3G 	Find an Eulerian Path in a Graph
#
# Input: A directed graph that contains an Eulerian path, where the graph
# is given in the form of an adjacency list.
#
# Return: An Eulerian path in this graph.

def find_eulerian_path(graph):
    def nodes():
        return list(set(list(graph.keys())+ \
                        [out for outs in graph.values() for out in outs]))        
    def get_start():
        start=[]
        for candidate in graph.keys():
            in_count=0
            for _,ins in graph.items():
                in_count+=sum([1 for i in ins if i==candidate])
            if in_count<len(graph[candidate]):
                start.append(candidate)
        return start
    
    def get_finish():
        finish=[]
        closed=[]
        for outs in graph.values():
            for candidate in outs:
                if not candidate in closed:
                    in_count=0
                    for os in graph.values():
                        for o in os:
                            if o==candidate:
                                in_count+=1
                    if not candidate in graph or in_count>len(graph[candidate]):
                        finish.append(candidate)                
                    closed.append(candidate) 
        return finish
    
    def adjust(path,start,finish):
        i=len(path)
        while i>0 and (path[0]!=start[0] or path[-1]!=finish[0]):
            path=path[1:]+[path[0]]
            i-=1
        return path

    start=get_start()
    finish=get_finish()
    if len(start)==1 and len(finish)==1:
        graph[finish[0]]=start[0:1]
        return adjust(find_eulerian_cycle(graph)[0:-1],start,finish)
    raise RosalindException(\
        'Start %(start)s and finish %(finish)s should have one element each'%locals())

# BA3H 	Reconstruct a String from its k-mer Composition 
#
# Input: An integer k followed by a list of k-mers Patterns.
#
# Return: A string Text with k-mer composition equal to Patterns. (If multiple
# answers exist, you may return any one.)

def reconstruct_from_kmers(k,patterns):
    return reconstruct(find_eulerian_path(deBruijn_collection(patterns)))

# BA3I 	Find a k-Universal Circular String 
# something off - I have taken this out of tests, even though the
# website accepts may answers. A small test case appears to give differnt
# answwrs on successive runs - maybe iterating through dict is not deterministic?

def k_universal_circular_string(k):
    def bits(i):
        return ('{0:b}'.format(i)).zfill(k)
    patterns=[bits(i) for i in range(2**k)]
    result= reconstruct(find_eulerian_cycle(deBruijn_collection(patterns)))[k-1:]
    for pattern in patterns:
        if not pattern in result:
            raise RosalindException('%(pattern)s %(result)s'%locals())   
    return result

# BA3J 	Reconstruct a String from its Paired Composition
#
# Input: Integers k and d followed by a collection of paired k-mers PairedReads.
#
# Return: A string Text with (k, d)-mer composition equal to PairedReads.
# (If multiple answers exist, you may return any one.)
def reconstruct_from_paired_kmers(k,d,patterns):
    def create_pair(string):
        [a,b]=string.split('|')
        return (a,b)
    def prefix(pair):
        (a,b)=pair
        return (a[0:-1],b[0:-1])
    def suffix(pair):
        (a,b)=pair
        return (a[1:],b[1:])
    def reconstruct_from_graph(path):
        def extract(pairs,i):
            (a,_)=pairs[i]
            return a        
        head_reconstruction=                         \
            ''.join(reconstruct_as_list(k-1,\
                                        len(path),\
                                        path,
                                        extract)[:-1])
        expected_length=len(patterns)+2*k+d-1
        deficit=expected_length-len(head_reconstruction)
        _,tail_reconstruction=path[-1]
        i=-2
        while len(tail_reconstruction)<deficit:
            _,tail=path[i]
            tail_reconstruction=tail[0]+tail_reconstruction
            i-=1
        return head_reconstruction+tail_reconstruction
    return reconstruct_from_graph(                            \
        find_eulerian_path(                                   \
            deBruijn_collection(                              \
                [create_pair(string) for string in patterns], \
                prefix,                                       \
                suffix)))

# BA3K 	Generate Contigs from a Collection of Reads 

def create_contigs(patterns):  
    contigs=[]
    for path in non_branching_paths(deBruijn_collection(patterns)):
        contig=path[0]
        for p in path[1:]:
            contig=contig+p[-1]
        contigs.append(contig)   
    return contigs

#BA3L 	Construct a String Spelled by a Gapped Genome Path
#   how does this differ from BA4J?

def construct_from_gapped(k,d,patterns):
    return reconstruct_from_paired_kmers(k,d,patterns)

#BA3M 	Generate All Maximal Non-Branching Paths in a Graph 

def non_branching_paths(graph):
    def add_counts(graph):
        counts={}
        for node,outs in graph.items():
            counts[node]=(0,len(outs))
        for node,outs in graph.items():
            for out in outs:
                if not out in counts:
                    counts[out]=(0,0)
                x,y=counts[out]
                counts[out]=x+1,y
        return counts
    def make_cycle(start,graph):
        def standardize(cycle):
            mm=min(cycle)
            ii=cycle.index(mm)
            return cycle[ii:]+cycle[:ii]+[mm]            
        cycle=[start]
        node=start

        while node in graph and len(graph[node])==1:
            succ=graph[node][0]
            if succ==start:
                return standardize(cycle)
            else:
                cycle.append(succ)
                node=succ

        return []
    
    def isolated_cycles(nodes,graph):
        cycles=[]
        for node in graph.keys():
            cycle=make_cycle(node,graph)
            if len(cycle)>0 and not cycle in cycles:
                cycles.append(cycle)
        return cycles
    
    paths=[]
    nodes=add_counts(graph)

    for v in nodes:
        (ins,outs)=nodes[v]
        if ins!=1 or outs!=1:
            if outs>0:
                for w in graph[v]:
                    nbp=[v,w]
                    w_in,w_out=nodes[w]
                    while w_in==1 and w_out==1:
                        u=graph[w][0]
                        nbp.append(u)
                        w=u
                        w_in,w_out=nodes[w]
                    paths.append(nbp)      
    return isolated_cycles(nodes,graph)+paths


### How do we sequence antibodies?

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

def prot(rna,table=rrt.codon_table):
    return ''.join([table[codon] for codon in rh.triplets(rna) if table[codon]!=';'])

# BA4B	Find Substrings of a Genome Encoding a Given Amino Acid String
#
# There are three different ways to divide a DNA string into codons for
# translation, one starting at each of the first three starting positions of
# the string. These different ways of dividing a DNA string into codons are
# called reading frames. Since DNA is double-stranded, a genome has six reading
# frames (three on each strand).
#
# We say that a DNA string Pattern encodes an amino acid string Peptide if
# the RNA string transcribed from either Pattern or its reverse complement
# Pattern translates into Peptide.
#
# Input: A DNA string Text and an amino acid string Peptide.
#
# Return: All substrings of Text encoding Peptide (if any such substrings exist)

def findEncodings(text,peptide):
    def encodes(dna):
        try:
            return prot(dna_to_rna(dna))==peptide
        except KeyError:
            return False    
    candidates=[text[i:i+3*len(peptide)] for i in range(len(text)-3)]
    return [rna for rna in candidates if encodes(rna) or encodes(revc(rna))]
            
# BA4C	Generate the Theoretical Spectrum of a Cyclic Peptide
#
# The workhorse of peptide sequencing is the mass spectrometer, an expensive
# molecular scale that shatters molecules into pieces and then weighs the 
# resulting fragments. The mass spectrometer measures the mass of a molecule
# in daltons (Da); 1 Da is approximately equal to the mass of a single nuclear
# particle (i.e., a proton or neutron).
#
# We will approximate the mass of a molecule by simply adding the number of
# protons and neutrons found in the molecule’s constituent atoms, which yields
# the molecule’s integer mass. For example, the amino acid "Gly", which has 
# chemical formula C2H3ON, has an integer mass of 57,
# since 2·12 + 3·1 + 1·16 + 1·14 = 57. Yet 1 Da is not exactly equal to the mass
# of a proton/neutron, and we may need to account for different naturally
# occurring isotopes of each atom when weighing a molecule. As a result, 
# amino acids typically have non-integer masses (e.g., "Gly" has total mass
# equal to approximately 57.02 Da); for simplicity, however, we will work with
# the integer mass table given in Figure 1.
#
# The theoretical spectrum of a cyclic peptide Peptide, denoted
# Cyclospectrum(Peptide), is the collection of all of the masses of its
# subpeptides, in addition to the mass 0 and the mass of the entire peptide.
# We will assume that the theoretical spectrum can contain duplicate elements,
# as is the case for "NQEL", where "NQ" and "EL" have the same mass.
#
# Input: An amino acid string Peptide.
#
# Return: Cyclospectrum(Peptide).

def cycloSpectrum(peptide,mass=rrt.integer_masses):
    def get_pairs(index_range):
        n=len(index_range)
        return [(i,j) for i in index_range for j in range(i,i+n) if j!=i]
    augmented_peptide=peptide+peptide   # allows easy extraction of substrings
                                        # fromcyclic peptide
    spectrum=[rh.get_mass(augmented_peptide[a:b],mass)\
              for (a,b) in get_pairs(range(len(peptide)))]
    spectrum.append(rh.get_mass('',mass))
    spectrum.append(rh.get_mass(peptide,mass)) 
    spectrum.sort()
    return spectrum

# BA4D	Compute the Number of Peptides of Given Total Mass
#
# In Generate the Theoretical Spectrum of a Cyclic Peptide, we generated the
# theoretical spectrum of a known cyclic peptide. Although this task is
# relatively easy, our aim in mass spectrometry is to solve the reverse problem:
# we must reconstruct an unknown peptide from its experimental spectrum. 
# We will start by assuming that a biologist is lucky enough to generate an
# ideal experimental spectrum Spectrum, which is one coinciding with the
# peptide’s theoretical spectrum. Can we reconstruct a peptide whose
# theoretical spectrum is Spectrum?
#
# Denote the total mass of an amino acid string Peptide as Mass(Peptide).
# In mass spectrometry experiments, whereas the peptide that generated a 
# spectrum is unknown, the peptide’s mass is typically known and is denoted
# ParentMass(Spectrum). Of course, given an ideal experimental spectrum,
# Mass(Peptide) is given by the largest mass in the spectrum.
#
# A brute force approach to reconstructing a peptide from its theoretical
# spectrum would generate all possible peptides whose mass is equal to 
# ParentMass(Spectrum) and then check which of these peptides has theoretical
# spectra matching Spectrum. However, we should be concerned about the running
# time of such an approach: how many peptides are there having mass equal
# to ParentMass(Spectrum)?
#
# Input: An integer m.
#
# Return: The number of linear peptides having integer mass m.
#
# NB, treat peptide as a vector of masses, so amino acids with the same
# mass are the same

def count_peptides_linear(total_mass):
    cache=[]
    masses=list(set(rrt.integer_masses.values()))
    for target_mass in range(total_mass+1):
        total=0
        for amino_acid_mass in masses:
            residual_mass=target_mass-amino_acid_mass
            if residual_mass==0:
                total+=1
            elif residual_mass>0:
                total+=cache[residual_mass]
        cache.append(total)
        
    return total


def get_weight(peptide):
    return sum(rrt.amino_acids[amino_acid].mon_mass for amino_acid in peptide)

def expand(peptides,masses):
    return [peptide+[mass] for peptide in peptides for mass in masses]
 
def mass(peptide):
    return sum([weight for weight in peptide])

def parentMass(spectrum):
    return max(spectrum)

# BA4E 	Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum
#
# Input: A collection of (possibly repeated) integers Spectrum corresponding
#       to an ideal experimental spectrum.
#
# Return: An amino acid string Peptide such that Cyclospectrum(Peptide) = 
#        Spectrum (if such a string exists).

def find_cyclopeptide_sequence(spectrum):
        
    def consistent(peptide,spectrum):
        def count(element,spect):
            return len ([s for s in spect if s==element])
        peptide_spectrum=linearSpectrum(peptide)
        for element in peptide_spectrum:
            if count(element,peptide_spectrum)>count(element,spectrum):
                return False
        return True
    
    def linearSpectrum(peptide):
        def get_pairs():
            return [(i,j) for i in range(len(peptide)) for j in range(i+1,len(peptide))]
        
        result=[sum(peptide[a:b]) for (a,b) in get_pairs()]
        result.append(0)
        result.append(sum(peptide)) 
        result.sort()
        return result 
        
    def cycloSpectrum(peptide):
        def get_pairs(index_range):
            n=len(index_range)
            return [(i,j) for i in index_range for j in range(i,i+n) if j!=i]
        augmented_peptide=peptide+peptide
        result=[sum(augmented_peptide[a:b]) for (a,b) in get_pairs(range(len(peptide)))]
        result.append(0)
        result.append(sum(peptide)) 
        result.sort()
        return result  
      
    peptides=[[]]
    output=[]
    masses=list(set(rrt.integer_masses.values()))
    
    while len(peptides)>0:
        next_peptides=[]
        for peptide in expand(peptides,masses):
            if mass(peptide) == parentMass(spectrum):
                if cycloSpectrum(peptide) == spectrum:
                    output.append(peptide)
            else:
                if consistent(peptide,spectrum):
                    next_peptides.append(peptide)    
        peptides=next_peptides
    return output

# BA4F 	Compute the Score of a Cyclic Peptide Against a Spectrum 
def score(peptide,spectrum,spect_from_peptide=cycloSpectrum):
    return rh.countMatchesInSpectra(spect_from_peptide(peptide),spectrum)

# BA4G 	Implement LeaderboardCyclopeptideSequencing 
def leaderPeptide(n,                                            \
                  spectrum,                                     \
                  masses=list(set(rrt.integer_masses.values())),\
                  spect1=rh.linearSpectrum,                     \
                  spect2=rh.linearSpectrum):
   
    
    leaderPeptide=[]
    leaderBoard=[leaderPeptide]
    
    while len(leaderBoard)>0:
        newBoard=[]
        for peptide in expand(leaderBoard,masses):
            if mass(peptide)==parentMass(spectrum):
                if score(peptide,spectrum,spect2)>\
                   score(leaderPeptide,spectrum,spect2):
                    leaderPeptide=peptide
                newBoard.append(peptide)
            elif mass(peptide)>parentMass(spectrum):
                pass #peptide will be dropped from leader board
            else:
                newBoard.append(peptide)
        leaderBoard=trim(newBoard, spectrum, n,spect1)
    return leaderPeptide

# BA4H 	Generate the Convolution of a Spectrum 
def convolution (spectrum):
    def create_counts(diffs):
        counts={}
        for diff in diffs:
            if not diff in counts:
                counts[diff]=0
            counts[diff]+=1
        return counts        
    diffs=[abs(spectrum[i]-spectrum[j])  \
           for i in range(len(spectrum))  \
           for j in range(i+1,len(spectrum)) if spectrum[i]!=spectrum[j]]
    return sorted([(diff,count)                                      \
                   for diff,count in create_counts(diffs).items()],  \
                  key=lambda x: (-x[1], x[0]))

# BA4I 	Implement ConvolutionCyclopeptideSequencing
#
# Given: An integer M, an integer N, and a collection of
# (possibly repeated) integers Spectrum.
#
# Return: A cyclic peptide LeaderPeptide with amino acids taken only from the
# top M elements (and ties) of the convolution of Spectrum that fall between
# 57 and 200, and where the size of Leaderboard is restricted to the top N (and ties).
#
# NB: I had to sort spectrum to pass the testcase in the textbook.

def convolutionCyclopeptideSequencing(m,n,spectrum,low_mass=57,high_mass=200):
    def get_masses_from_spectrum():
        masses=[]
        last_count=0
        for mass,count in convolution(spectrum):
            if low_mass<=mass and mass<=high_mass:
                if len(masses)<m:
                    masses.append(mass)
                    last_count=count
                else:
                    if count==last_count:
                        masses.append(mass)
                    else:
                        return masses
        return masses
    
    return leaderPeptide(n,spectrum,get_masses_from_spectrum(),spect2=rh.cycloSpectrum1)

# BA4L 	Trim a Peptide Leaderboard	
#
# Input: A leaderboard of linear peptides Leaderboard, a linear spectrum
# Spectrum, and an integer N.
#
# Return: The top N peptides from Leaderboard scored against Spectrum.
# Remember to use LinearScore.

def trim(leaderBoard, spectrum,n,spectrum_generator=rh.linearSpectrum):
    if len(leaderBoard)<n:
        return leaderBoard
    peptides_with_scores=[\
        (score(peptide,spectrum,spectrum_generator),peptide)\
        for peptide in leaderBoard]
    peptides_with_scores.sort(reverse=True)   
    (cutoff,_)= peptides_with_scores[n-1]
    return [peptide                                     \
            for (score,peptide) in peptides_with_scores \
            if score>=cutoff]

# Adapter for 'trim', so it will work with peptides as strings

def trim_for_strings(leaderBoard, spectrum,n,\
                     spectrum_generator=rh.linearSpectrum,\
                     masses=rrt.integer_masses):
    numeric_peptides=[]
    trans={}
    for peptide in leaderBoard:
        num=[masses[amino_acid] for amino_acid in peptide]
        numeric_peptides.append(num)
        trans[tuple(num)]=peptide
    trimmed=trim(numeric_peptides, spectrum,n,spectrum_generator)
    return [trans [tuple(peptide)] for peptide in trimmed]

# BA4M 	Solve the Turnpike Problem 
def turnpike(differences):
    def diffs(seq):
        return sorted([x-y for x in seq for y in seq if x>y]) 
    def match(rs,ss):
        for r,s in zip(rs,ss):
            if r!=s:
                return False
        return True    
    def get_length():
        n_zeroes=0
        i=0
        j=len(differences)-1
        index_first_pos=-1
        while i<j:
            if differences[i]+differences[j]!=0:
                raise RosalindException('Mismatch %(i)d and %(j)d'%locals())
            if differences[i]==0:
                n_zeroes+=1
            else:
                index_first_pos=j
            i+=1
            j-=1
        n_zeroes = 2*n_zeroes if len(differences)%2==0 else 2*n_zeroes+1
        if n_zeroes*n_zeroes==len(differences):
            return (index_first_pos,n_zeroes)
        else:
            raise RosalindException('Mismatched lengths')
        

    index_first_pos,length_result =get_length()
    largest=differences[-1]
    indices=[]
    while len(indices)<length_result-2:
        indices.append(index_first_pos+1)
    
    while True:
        result=[0]
        for index in indices:
            result.append(differences[index])
        result.append(largest)
       
        if match(diffs(result),differences[index_first_pos:]):
            return result
        j=len(indices)-1
        while j>0:
            if indices[j]<len(differences)-1:
                indices[j]+=1
                break
            else:
                indices[j]=index_first_pos+1
                j-=1
            
        
# BA4J 	Generate the Theoretical Spectrum of a Linear Peptide 
def linearSpectrumFromString(peptide):
    return rh.linearSpectrum([rrt.integer_masses[a] for a in peptide])

# BA4K 	Compute the Score of a Linear Peptide 
def linearScore(peptide,spectrum):
    return rh.countMatchesInSpectra(linearSpectrumFromString(peptide),spectrum)

def convolution_expanded(spectrum):
    return [diff for (diff,count) in convolution (spectrum) for i in range(count) ]
    
# BA5A 	Find the Minimum Number of Coins Needed to Make Change 	
#
# Input: An integer money and an array Coins of positive integers.
#
# Return: The minimum number of coins with denominations Coins that changes money.
#
# http://rosalind.info/problems/ba5a/

def number_of_coins(money,coins):
    number = [0]                           # We will use Dynamic Programming, and solve 
                                           # the problem for each amount up to and including money
    for m in range(1,money+1):             # solve for m
        nn = sys.float_info.max            # Number of coins: assume that we haven't solved
        for coin in coins:                 # Find a coin such that we can make change using it
                                           # plus a previoudly comuted value
            if m>=coin:
                if number[m-coin]+1<nn:
                    nn = number[m-coin]+1
        number.append(nn)
    return number[money]

# BA5B 	Find the Length of a Longest Path in a Manhattan-like Grid 
#
# Input: Integers n and m, followed by an n*(m+1) matrix Down and an
#        (n+1)*m matrix Right. The two matrices are separated by the "-" symbol.
#
# Return: The length of a longest path from source (0, 0) to sink (n, m)
#        in the n*m rectangular grid whose edges are defined by the matrices
#        Down and Right.
#
# http://rosalind.info/problems/ba5a/

def longest_manhattan_path(n,m,down,right):
    s=[]
    for i in range(n+1):
        s.append(rh.zeroes(m+1))

    for i in range(1,n+1):
        s[i][0]=s[i-1][0]+down[i-1][0]
        
    for j in range(1,m+1):
        s[0][j]=s[0][j-1]+right[0][j-1]
        
    for i in range(1,n+1):    
        for j in range(1,m+1):
            s[i][j]=max(s[i-1][j]+down[i-1][j],s[i][j-1]+right[i][j-1])

    return s[n][m]

# BA5C 	Find a Longest Common Subsequence of Two Strings
#
# Input: Two strings.
#
# Return: A longest common subsequence of these strings.
#
# http://rosalind.info/problems/ba5a/

def longest_common_subsequence(string1,string2):

    # Calculate longest path through "map" defined by the two strings
    #
    # At each point we have a state, s, which is defined as follows
    # 1. Count of matches (==path lebth from (0,0) to here
    # 2. horizonal position of predecessor of this point
    # 3. vertical position of predecessor of this point
    # 4. If this point corresponds to a matching character (in both strings)
    #    this should be that chracter, otherwise an empty string
    
    def longest_path():
        
        # Used to update distamce at each node     
        def new_s(i,j,s):
            count_horizonal,_,_,_=s[i-1][j]
            count_vertical,_,_,_=s[i][j-1]
            count_diagonal,_,_,_=s[i-1][j-1]
            if string1[i-1]==string2[j-1]:
                count_diagonal+=1
            count=max(count_horizonal,count_vertical,count_diagonal)
            if count_diagonal==count:
                return (count,\
                        i-1,\
                        j-1,
                        string1[i-1] if string1[i-1]==string2[j-1] else '')
            elif count_vertical==count:
                return (count,i,j-1,'')
            else: #count_horizonal==count
                return (count,i-1,j,'')
            
        m1=len(string1)+1
        m2=len(string2)+1
        s=[]
        for i in range(m1):
            ss=[]
            for j in range(m2):
                ss.append((0,-1,-1,''))
            s.append(ss)
        for i in range(1,m1):    
            for j in range(1,m2):
                s[i][j]=new_s(i,j,s)
        return s
    
    # Build string from status
    
    def construct_string(s):
        i=len(string1)
        j=len(string2)
        result=[]
        while i>-1 and j>-1:
            _,i,j,chars=s[i][j]
            result.append(chars)
        return ''.join(result[::-1])
    
    return construct_string(longest_path())

# BA5D 	Find the Longest Path in a DAG  	
#
# Input: An integer representing the source node of a graph, followed by an integer
#        representing the sink node of the graph, followed by an edge-weighted graph. 
#        The graph is represented by a modified adjacency list in which the notation "0->1:7"
#        indicates that an edge connects node 0 to node 1 with weight 7.
#
# Return: The length of a longest path in the graph, followed by a longest path. 
#         (If multiple longest paths exist, you may return any one.)
#
# http://rosalind.info/problems/ba5d/

def longest_path(source,sink,graph):
    def initialize_s():
        s={}
        for a,b,_ in graph:
            s[a]=-sys.float_info.max
            s[b]=-sys.float_info.max
        s[source]=0
        return s
    
    def create_adjacency_list():
        adjacency_list={}
        for a,b,w in graph:
            if not a in adjacency_list:
                adjacency_list[a]=[]
            adjacency_list[a].append(b)
        return adjacency_list
   
    def create_weights():
        weights={}
        for a,b,w in graph:
            weights[(a,b)]=w
        return weights
        
    def calculate_distances(ordering):
        s=initialize_s()
        weights=create_weights()        
        predecessors={}
        for b in ordering:
            for a in ordering:
                if a==b:
                    break
                
                new_s=max(s[b],s[a]+(weights[(a,b)] if (a,b) in weights else 0))
                if new_s>s[b]:
                    s[b]=new_s
                    predecessors[b]=a
        return (s,predecessors)
   
    def create_path(predecessors):
        path=[sink]
        node=sink
        while node in predecessors:
            node=predecessors[node]
            path.append(node)
        return path
    
    s,predecessors=calculate_distances(topological_order(create_adjacency_list()))
    
    return (s[sink],create_path(predecessors)[::-1])
    
# BA5E 	Find a Highest-Scoring Alignment of Two Strings
# BA5F 	Find a Highest-Scoring Local Alignment of Two Strings  
#
# common code

# create_distance_matrix
def create_distance_matrix(nrows,ncolumns,initial_value=-sys.float_info.max):
    s=[]
    for i in range(nrows):
        row=[]
        for j in range(ncolumns):
            row.append(initial_value)        
        s.append(row)
    s[0][0]=0
    return s

def calculate_scores_for_alignment(s,string1, string2, weights,sigma,init_predecessors=None):
    predecessors={}
    for i in range(len(string1)+1):
        for j in range(len(string2)+1):
            predecessors[(i,j)]=init_predecessors
            if i>0:
                s_new=s[i-1][j]-sigma
                if s_new>s[i][j]:
                    s[i][j]=s_new
                    predecessors[(i,j)]=(-1,0,i-1,j)
            if j>0:
                s_new=s[i][j-1]-sigma
                if s_new>s[i][j]:
                    s[i][j]=s_new
                    predecessors[(i,j)]=(0,-1,i,j-1)            
                if i>0:
                    s_new=s[i-1][j-1]+weights[(string1[i-1],string2[j-1])]
                    if s_new>s[i][j]:
                        s[i][j]=s_new
                        predecessors[(i,j)]=(-1,-1,i-1,j-1)
    return (s,predecessors)

def create_alignment(string1, string2,s_predecessors,i_start=-1,j_start=-1):
    s,predecessors = s_predecessors
    result1        = []
    result2        = []
    i              = len(string1) if i_start==-1 else i_start
    j              = len(string2) if j_start==-1 else j_start
    while i>0 or j>0:
        x,y,i,j=predecessors[(i,j)]
        if x==-1 and y==0:
            result1.append(string1[i])
            result2.append('-')
        elif x==0 and y==-1:
            result1.append('-')
            result2.append(string2[j])
        elif x==-1 and y==-1:
            result1.append(string1[i])
            result2.append(string2[j])
        
    return (s[len(string1)][len(string2)],\
            ''.join(result1[::-1]),       \
            ''.join(result2[::-1]))

# BA5E 	Find a Highest-Scoring Alignment of Two Strings
# Find the highest-scoring alignment between two strings using a scoring matrix.
#
# Input: Two amino acid strings.
#
# Return: The maximum alignment score of these strings followed by an
#         alignment achieving this maximum score. Use the BLOSUM62 scoring matrix
#         and indel penalty σ = 5. (If multiple alignments achieving the maximum 
#         score exist, you may return any one.)

def highest_scoring_global_alignment(string1,string2,weights=rrt.createBLOSUM62(),sigma=5):
    return create_alignment(string1,
                            string2,
                            calculate_scores_for_alignment(
                                create_distance_matrix(len(string1)+1,\
                                                       len(string2)+1),\
                                string1,\
                                string2,\
                                weights,\
                                sigma))


# BA5F 	Find a Highest-Scoring Local Alignment of Two Strings 
#
# Input: Two amino acid strings.
#
# Return: The maximum score of a local alignment of the strings, followed by
# a local alignment of these strings achieving the maximum score. Use the
# PAM250 scoring matrix and indel penalty  5. (If multiple local alignments
# achieving the maximum score exist, you may return any one.)

def highest_scoring_local_alignment(string1,string2,weights=rrt.createPAM250(),sigma=5):
    def find_best_substring(s):
        s_right=0    
        i_best=-1
        j_best=-1
        predecessor=(0,0,0,0)
        for i in range(len(string1)+1):
            for j in range(len(string2)+1):
                if s[i][j]>s_right:
                    i_best=i
                    j_best=j
                    s_right=s[i][j]
                    predecessor=(rh.sign(len(string1)-i_best),\
                                 rh.sign(len(string2)-j_best),\
                                 i_best,\
                                 j_best)
        if s_right>s[len(string1)][len(string2)]:
            s[len(string1)][len(string2)]=s_right
        else:
            i_best=len(string1)
            j_best=len(string2)
            predecessor=(0,0,0,0)
        return (i_best,j_best,s,predecessor)
   
    s,predecessors=calculate_scores_for_alignment(
        create_distance_matrix(len(string1)+1,len(string2)+1,0),
        string1,
        string2,
        weights,
        sigma,
        (0,0,0,0))
   
    i_best,j_best,s,predecessor=find_best_substring(s)
    if predecessor!=(0,0,0,0):
        predecessors[(len(string1),len(string2))]=predecessor
    
    return create_alignment(string1,string2,(s,predecessors),i_best,j_best)
	  	  	 
# BA5G 	Compute the Edit Distance Between Two Strings 	 	
	  	  	 
# BA5H 	Find a Highest-Scoring Fitting Alignment of Two Strings 	 	
	  	  	 
# BA5I 	Find a Highest-Scoring Overlap Alignment of Two Strings 	 	
	  	  	 
# BA5J 	Align Two Strings Using Affine Gap Penalties 	 	
	  	  	 
# BA5K 	Find a Middle Edge in an Alignment Graph in Linear Space 	 	
	  	  	 
# BA5L 	Align Two Strings Using Linear Space 	 	
	  	  	 
# BA5M 	Find a Highest-Scoring Multiple Sequence Alignment 	 	
	  	  	 
# BA5N 	Find a Topological Ordering of a DAG 
#
# Input: The adjacency list of a graph (with nodes represented by integers).
#
# Return: A topological ordering of this graph.

def topological_order(graph):
    def number_incoming(node):
        n=0
        for out in graph.values():
            if node in out:
                n+=1
        return n
    
    ordering=[]
    candidates=[node for node in graph.keys() if number_incoming(node)==0]
    while len(candidates)>0:
        a=candidates.pop()
        ordering.append(a)
        if a in graph:
            bs=[b for b in graph[a]]
            del graph[a]
            for b in bs:
                if number_incoming(b)==0:
                    candidates.append(b)
                    
    if len(graph)>0:
        raise RosalindException('Input graph is not a DAG')
    
    return ordering





### Test cases ###

if __name__=='__main__':
 
    import unittest,fragile
    
    class TestRosalind(unittest.TestCase):
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
            fasta=f.FastaContent(string.split('\n'))
            r=gc(fasta)
            self.assertEqual('Rosalind_0808',r[0])
            self.assertAlmostEqual(60.919540,r[1],places=5)
  
        def test_sseq(self):
            string='''>Rosalind_14
            ACGTACGTGACG
            >Rosalind_18
            GTA'''
            fasta=f.FastaContent(string.split('\n'))
            ss=sseq(fasta)
            #print (ss)
            self.assertEqual(3,len(ss))
            self.assertEqual(3,ss[0])
            self.assertEqual(8,ss[1])
            self.assertEqual(10,ss[2])
                
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
            fasta=f.FastaContent(string.split('\n'))
            self.assertEqual(set(['TA', 'CA', 'AC']),set(lcsm(fasta)))
            
        def test_prot(self):
            self.assertEqual('MAMAPRTEINSTRING',prot('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'))
            
        def test_mrna(self):
            self.assertEqual(12,mrna("MA"))
            
        def test_iprb(self):
            self.assertAlmostEqual(0.78333,iprb(2,2,2),places=5)

# LIA 	Independent Alleles 

        def test_lia(self):
            self.assertAlmostEqual(0.684,lia(2,1),places=3)

# RSTR 	Matching Random Motifs        
        def test_rstr(self):
            self.assertAlmostEqual(0.689,rstr(90000, 0.6,'ATAGCCGA'),places=3)
            
        def test_fib(self):
            self.assertEqual(19,fib(5,3))
            self.assertEqual(875089148811941,fib(35, 5))
            
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
            consensus,profile=cons(f.FastaContent(string.split('\n')))
            self.assertEqual('ATGCAACT',consensus)
            self.assertEqual(profile['T'],[1,5, 0, 0, 0, 1, 1, 6])
            self.assertEqual(profile['G'],[1, 1, 6, 3, 0, 1, 0, 0])
            self.assertEqual(profile['A'],[5, 1, 0, 0, 5, 5, 0, 0])
            self.assertEqual(profile['C'],[0, 0, 1, 4, 2, 0, 6, 1])
         
        def test_fibd(self):
            self.assertEqual(4,fibd(6,3))
        
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
            graph=grph(f.FastaContent(string.split('\n')),3)
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
            self.assertEqual('MVYIADKQHVASREAYGHMFKVCA',splc(f.FastaContent(string.split('\n'))))
            
        def test_iev(self):
            self.assertEqual(3.5,iev([1,0,0,1,0,1]))
            
        def test_revp(self):
            string='''>Rosalind_24
            TCAATGCATGCGGGTCTATATGCAT'''
            palindromes=revp(f.FastaContent(string.split('\n')))
            self.assertIn((4, 6),palindromes)
            self.assertIn((5, 4),palindromes)
            self.assertIn((6, 6),palindromes)
            self.assertIn((7, 4),palindromes)
            self.assertIn((17, 4),palindromes)
            self.assertIn((18, 4),palindromes)
            self.assertIn((20, 6),palindromes)
            self.assertIn((21, 4),palindromes)
            self.assertEqual(8,len(palindromes))

# KMER  Generalizing GC-Content
        def test_kmer(self):
            string='''>>Rosalind_6431
            CTTCGAAAGTTTGGGCCGAGTCTTACAGTCGGTCTTGAAGCAAAGTAACGAACTCCACGG
            CCCTGACTACCGAACCAGTTGTGAGTACTCAACTGGGTGAGAGTGCAGTCCCTATTGAGT
            TTCCGAGACTCACCGGGATTTTCGATCCAGCCTCAGTCCAGTCTTGTGGCCAACTCACCA
            AATGACGTTGGAATATCCCTGTCTAGCTCACGCAGTACTTAGTAAGAGGTCGCTGCAGCG
            GGGCAAGGAGATCGGAAAATGTGCTCTATATGCGACTAAAGCTCCTAACTTACACGTAGA
            CTTGCCCGTGTTAAAAACTCGGCTCACATGCTGTCTGCGGCTGGCTGTATACAGTATCTA
            CCTAATACCCTTCAGTTCGCCGCACAAAAGCTGGGAGTTACCGCGGAAATCACAG'''
            fasta=f.FastaContent(string.split('\n'))
            r=lexig(4,fasta)
            self.assertEqual( \
                [4, 1, 4, 3, 0, 1, 1, 5, 1, 3, 1, 2, 2, 1, 2, 0, 1, 1, 3, 1, 2,\
                 1, 3, 1, 1, 1, 1, 2, 2, 5, 1, 3, 0, 2, 2, 1, 1, 1, 1, 3, 1, 0,\
                 0, 1, 5, 5, 1, 5, 0, 2, 0, 2, 1, 2, 1, 1, 1, 2, 0, 1, 0, 0, 1,\
                 1, 3, 2, 1, 0, 3, 2, 3, 0, 0, 2, 0, 8, 0, 0, 1, 0, 2, 1, 3, 0,\
                 0, 0, 1, 4, 3, 2, 1, 1, 3, 1, 2, 1, 3, 1, 2, 1, 2, 1, 1, 1, 2,\
                 3, 2, 1, 1, 0, 1, 1, 3, 2, 1, 2, 6, 2, 1, 1, 1, 2, 3, 3, 3, 2,\
                 3, 0, 3, 2, 1, 1, 0, 0, 1, 4, 3, 0, 1, 5, 0, 2, 0, 1, 2, 1, 3,\
                 0, 1, 2, 2, 1, 1, 0, 3, 0, 0, 4, 5, 0, 3, 0, 2, 1, 1, 3, 0, 3,\
                 2, 2, 1, 1, 0, 2, 1, 0, 2, 2, 1, 2, 0, 2, 2, 5, 2, 2, 1, 1, 2,\
                 1, 2, 2, 2, 2, 1, 1, 3, 4, 0, 2, 1, 1, 0, 1, 2, 2, 1, 1, 1, 5,\
                 2, 0, 3, 2, 1, 1, 2, 2, 3, 0, 3, 0, 1, 3, 1, 2, 3, 0, 2, 1, 2,\
                 2, 1, 2, 3, 0, 1, 2, 3, 1, 1, 3, 1, 0, 1, 1, 3, 0, 2, 1, 2, 2,\
                 0, 2, 1, 1],
                r)  

        #def test_rear(self):
            #self.assertEqual([9, 4, 5, 7, 0],
                        #rear([
                            #([1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                            #[3, 1, 5, 2, 7, 4, 9, 6, 10, 8]),
                            
                            #([3, 10, 8, 2, 5, 4, 7, 1, 6, 9],
                            #[5, 2, 3, 1, 7, 4, 10, 8, 6, 9]),
                            
                            #([8, 6, 7, 9, 4, 1, 3, 10, 2, 5],
                            #[8, 2, 7, 6, 9, 1, 5, 3, 10, 4]),
                            
                            #([3, 9, 10, 4, 1, 8, 6, 7, 5, 2],
                            #[2, 9, 8, 5, 1, 7, 3, 4, 6, 10]),
                            
                            #([1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                            #[1, 2, 3, 4, 5, 6, 7, 8, 9, 10])                            
                        #]))
        
        def test_tran(self):
            string='''>Rosalind_0209
            GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGA
            AGTACGGGCATCAACCCAGTT
            >Rosalind_2200
            TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGC
            GGTACGAGTGTTCCTTTGGGT'''
            fasta=f.FastaContent(string.split('\n'))            
            self.assertAlmostEqual(1.21428571429,tran(fasta),places=5)

# PDST 	Creating a Distance Matrix 
        def test_pdst(self):
            string='''>Rosalind_9499
            TTTCCATTTA
            >Rosalind_0942
            GATTCATTTC
            >Rosalind_6568
            TTTCCATTTT
            >Rosalind_1833
            GTTCCATTTA'''
            fasta=f.FastaContent(string.split('\n'))            
            self.assertEqual([[0.00000, 0.40000, 0.10000, 0.10000],
                              [0.40000, 0.00000, 0.40000, 0.30000],
                              [0.10000, 0.40000, 0.00000, 0.20000],
                              [0.10000, 0.30000, 0.20000, 0.00000]],
                             distance_matrix(fasta))

# ASPC 	Introduction to Alternative Splicing 

        def test_aspc(self):
            self.assertEqual(42,aspc(6,3))
 
# PPER 	Partial Permutations        
        def test_pper(self):
            self.assertEqual(51200,pper(21,7))            

#INDC 	Independent Segregation of Chromosomes 
        def test_indc(self):
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

# AFRQ 	Counting Disease Carriers
        def test_afrq(self):
            aa=afrq([0.1, 0.25, 0.5])
            self.assertAlmostEqual(0.532,aa[0],3)
            self.assertAlmostEqual(0.75,aa[1],3)
            self.assertAlmostEqual(0.914,aa[2],3)

# WFMD 	The Wright-Fisher Model of Genetic Drift
        def test_wfmd(self):
            self.assertAlmostEqual(0.772,wfmd(4, 6, 2, 1),3)

# EBIN 	Wright-Fisher's Expected Behavior
        def test_ebin(self):
            B=ebin(17,[0.1, 0.2, 0.3])
            self.assertAlmostEqual(1.7,B[0],3)
            self.assertAlmostEqual(3.4,B[1],3)
            self.assertAlmostEqual( 5.1,B[2],3)

# FOUN 	The Founder Effect and Genetic Drift
        def test_foun(self):
            B=foun(4,3,[0,1,2])
            self.assertAlmostEqual(0.0,              B[0][0],5)
            self.assertAlmostEqual(-0.463935575821,  B[0][1],5)
            self.assertAlmostEqual(-0.999509892866,  B[0][2],5)
            self.assertAlmostEqual(0.0,              B[1][0],5)
            self.assertAlmostEqual(-0.301424998891,  B[1][1],5)
            self.assertAlmostEqual(-0.641668367342,  B[1][2],5) 
            self.assertAlmostEqual(0.0,              B[2][0],5)
            self.assertAlmostEqual(-0.229066698008,  B[2][1],5)
            self.assertAlmostEqual(-0.485798552456,  B[2][2],5)            
 
# SEXL 	Sex-Linked Inheritance 
        def test_sexl(self):
            B=sexl([0.1, 0.5, 0.8])
            self.assertAlmostEqual(0.18,B[0],3)
            self.assertAlmostEqual(0.5,B[1],3)
            self.assertAlmostEqual(0.32,B[2],3)

# EVAL 	Expected Number of Restriction Sites   
        def test_eval(self):
            B=eval(10,'AG',[0.25, 0.5, 0.75])
            self.assertAlmostEqual(0.422,B[0],3)
            self.assertAlmostEqual(0.563,B[1],3)
            self.assertAlmostEqual(0.422,B[2],3)
            
### Where in the Genome does DNA raplication begin? ###
            
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
 
        def test_ba1f(self):
            positions,_=find_minimum_skew(\
                'CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCT'\
                'TGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG')
            self.assertIn(53,positions)
            self.assertIn(97,positions)
        
        def test_ba6e(self):
            pairs=fragile.find_shared_kmers(3,'AAACTCATC','TTTCAAATC')
            self.assertIn((0, 4),pairs)
            self.assertIn((0, 0),pairs)
            self.assertIn((4, 2),pairs)
            self.assertIn((6, 6),pairs)
    
        def test_ba1e(self):
            clumps=findClumps('CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAAC'\
                              'ACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC',5,75,4)
            self.assertIn('CGACA',clumps)
            self.assertIn('GAAGA',clumps)
            self.assertIn('AATGT',clumps)
        
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

# BA2A 	Implement MotifEnumeration

        def test_ba2a(self):
            self.assertEqual('ATA ATT GTT TTT',
                             enumerateMotifs(3,
                                   1,
                                   ['ATTTGGC',
                                    'TGCCTTA',
                                    'CGGTATC',
                                    'GAAAATT']))

# BA2B 	Find a Median String 

        def test_ba2b(self):
            self.assertEqual('GAC',
                             medianString(3,[
                                 'AAATTGACGCAT',
                                 'GACGACCACGTT',
                                 'CGTCAGCGCCTG',
                                 'GCTGAGCACCGG',
                                 'AGTACGGGACAG'                               
                             ]))

          
            
# BA2C 	Find a Profile-most Probable k-mer in a String
        def test_ba2c(self):
            self.assertEqual(
                'CCGAG',
                mostProbable(
                    'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT',\
                    5,
                    [[0.2, 0.2, 0.3, 0.2, 0.3],
                     [0.4, 0.3, 0.1, 0.5, 0.1],
                     [0.3, 0.3, 0.5, 0.2, 0.4],
                     [0.1, 0.2, 0.1, 0.1, 0.2]]))

# BA2D 	Implement GreedyMotifSearch
        def test_ba2d(self):
            motifs=greedyMotifSearch(3,
                                     5,
                                     [
                                         'GGCGTTCAGGCA',
                                         'AAGAATCAGTCA',
                                         'CAAGGAGTTCGC',
                                         'CACGTCAATCAC',
                                         'CAATAATATTCG'                
                                     ])
            self.assertEqual(['CAA','CAA','CAA','CAG','CAG'],sorted(motifs))

       
# BA2E 	Implement GreedyMotifSearch with Pseudocounts  
        def test_ba2e(self):
            motifs=greedyMotifSearch(3,
                                     5,
                                     [
                                         'GGCGTTCAGGCA',
                                         'AAGAATCAGTCA',
                                         'CAAGGAGTTCGC',
                                         'CACGTCAATCAC',
                                         'CAATAATATTCG'                
                                     ],
                                     pseudo_counts=True)
            self.assertEqual(['ATC','ATC','TTC','TTC','TTC'],sorted(motifs))

# BA2F 	Implement RandomizedMotifSearch	
        #def test_ba2f(self):
            #(c,x)=randomized_motif_search_driver(8, 5,[
                #'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
                #'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
                #'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
                #'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
                #'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'],100000)
            #print (c)
            #self.assertIn('TCTCGGGG',x)
            #self.assertIn('CCAAGGTG',x)
            #self.assertIn('TACAGGCG',x)
            #self.assertIn('TTCAGGTG',x)
            #self.assertIn('TCCACGTG',x)

# BA2G 	Implement GibbsSampler	
        def test_ba2g(self):
            x=gibbs(8, 5, 100,[
                'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
                'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
                'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
                'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
                'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'])

# BA2H 	Implement DistanceBetweenPatternAndStrings 
        def test_ba2h(self):
            self.assertEqual(
                5,
                distanceBetweenPatternAndStrings('AAA',
                                                 ['TTACCTTAAC',
                                                  'GATATCTGTC',
                                                  'ACGGCGTTCG',
                                                  'CCCTAAAGAG',
                                                  'CGTCAGAGGT']))
            
# BA3B	Reconstruct a String from its Genome Path
        def test_ba3b(self):
            self.assertEqual("ACCGAAGCT",
                             reconstruct(["ACCGA",
                                          "CCGAA",
                                          "CGAAG",
                                          "GAAGC",
                                          "AAGCT"] ))
            
# BA3C	Construct the Overlap Graph of a Collection of k-mers
        def test_ba3c(self):
            graph=grph_kmers([
                'ATGCG',
                'GCATG',
                'CATGC',
                'AGGCA',
                'GGCAT'        
            ])
            self.assertEqual(4,len(graph))
            self.assertIn(('AGGCA','GGCAT'),graph)
            self.assertIn(('CATGC','ATGCG'),graph)
            self.assertIn(('GCATG','CATGC'),graph)
            self.assertIn(('GGCAT','GCATG'),graph)

# BA3D 	Construct the De Bruijn Graph of a String 
        def test_ba3d(self):
            graph=deBruijn(4,'AAGATTCTCTAC')
            self.assertIn(('AAG',['AGA']),graph)
            self.assertIn(('AGA',['GAT']),graph)
            self.assertIn(('ATT',['TTC']),graph)
            self.assertIn(('CTA',['TAC']),graph)
            self.assertIn(('CTC',['TCT']),graph)
            self.assertIn(('GAT',['ATT']),graph)
            self.assertIn(('TCT',['CTA','CTC']),graph)
            self.assertIn(('TTC',['TCT']),graph)            
            self.assertEqual(8,len(graph))

#BA3E 	Construct the De Bruijn Graph of a Collection of k-mers
        def test_ba3e(self):
            graph=deBruijn_collection(
                ['GAGG',
                 'CAGG',
                 'GGGG',
                 'GGGA',
                 'CAGG',
                 'AGGG',
                 'GGAG'])
            self.assertEqual(graph['AGG'],['GGG'])
            self.assertEqual(graph['CAG'],['AGG','AGG'])
            self.assertEqual(graph['GAG'],['AGG'])
            self.assertEqual(graph['GGA'],['GAG'])
            self.assertEqual(graph['GGG'],['GGA','GGG'])   
            self.assertEqual(5,len(graph.keys()))

        def assertCyclicEqual(self,a,b):
            self.assertEqual(len(a),len(b))
            for i in range(len(a)):
                if a==rh.rotate(b,i):
                    return
            self.fail('Cycles %(a)s and %(b)s do not match'%locals())
            
#BA3F 	Find an Eulerian Cycle in a Graph
        def test_ba3f(self):
            self.assertCyclicEqual([6,8,7,9,6,5,4,2,1,0,3,2,6],
                             find_eulerian_cycle(
                                 {0 : [3],
                                  1 : [0],
                                  2 : [1, 6],
                                  3 : [2],
                                  4 : [2],
                                  5 : [4],
                                  6 : [5, 8],
                                  7 : [9],
                                  8 : [7],
                                  9 : [6]}
                             ))

#BA3G 	Find an Eulerian Path in a Graph 
        def test_ba3g(self):
            self.assertEqual([6,7,8,9,6,3,0,2,1,3,4],find_eulerian_path(
                {
                    0 : [2],
                    1 : [3],
                    2 : [1],
                    3 : [0,4],
                    6 : [3,7],
                    7 : [8],
                    8 : [9],
                    9 : [6]                
                }
            ))

# BA3H 	Reconstruct a String from its k-mer Composition 
        def test_ba3h(self):
            self.assertEqual('GGCTTACCA',
                             reconstruct_from_kmers(
                                 4,
                                 ['CTTA',
                                  'ACCA',
                                  'TACC',
                                  'GGCT',
                                  'GCTT',
                                  'TTAC']))

# BA3I 	Find a k-Universal Circular String 
        #def test_ba3i(self):
            #self.assertEqual('0110010100001111',k_universal_circular_string(4))

# BA3J 	Reconstruct a String from its Paired Composition 
        def test_ba3j(self):
            self.assertEqual("GTGGTCGTGAGATGTTGA",
                             reconstruct_from_paired_kmers(
                                 4,
                                 2,
                                ['GAGA|TTGA',
                                'TCGT|GATG',
                                'CGTG|ATGT',
                                'TGGT|TGAG',
                                'GTGA|TGTT',
                                'GTGG|GTGA',
                                'TGAG|GTTG',
                                'GGTC|GAGA',
                                'GTCG|AGAT']))
            self.assertEqual('CACCGATACTGATTCTGAAGCTT',
                             reconstruct_from_paired_kmers(3,
                                           1,
                                           [
                                              'ACC|ATA', #(ACC|ATA) 
                                              'ACT|ATT', #(ACT|ATT)
                                              'ATA|TGA', #(ATA|TGA)
                                              'ATT|TGA', #(ATT|TGA)
                                              'CAC|GAT', #(CAC|GAT)
                                              'CCG|TAC', #(CCG|TAC)
                                              'CGA|ACT', #(CGA|ACT)
                                              'CTG|AGC', #(CTG|AGC)
                                              'CTG|TTC', #(CTG|TTC)
                                              'GAA|CTT', #(GAA|CTT)
                                              'GAT|CTG', #(GAT|CTG)
                                              'GAT|CTG', #(GAT|CTG)
                                              'TAC|GAT', #(TAC|GAT)
                                              'TCT|AAG', #(TCT|AAG)
                                              'TGA|GCT', #(TGA|GCT)
                                              'TGA|TCT', #(TGA|TCT)
                                              'TTC|GAA'  #(TTC|GAA)
                                           ]))

# BA3K 	Generate Contigs from a Collection of Reads
        def test_ba3k(self):
            contigs=create_contigs([
                'ATG',
                'ATG',
                'TGT',
                'TGG',
                'CAT',
                'GGA',
                'GAT',
                'AGA'            
            ])
            self.assertIn('AGA',contigs)
            self.assertIn('ATG',contigs)
            self.assertIn('ATG',contigs)
            self.assertIn('CAT',contigs)
            self.assertIn('GAT',contigs)
            self.assertIn('TGGA',contigs)
            self.assertIn('TGT',contigs)            

#BA3L 	Construct a String Spelled by a Gapped Genome Path 
        def test_ba3l(self):
            self.assertEqual('GACCGAGCGCCGGA',
                             construct_from_gapped(4,
                                                   2,
                                                   [
                                                       'GACC|GCGC',
                                                       'ACCG|CGCC',
                                                       'CCGA|GCCG',
                                                       'CGAG|CCGG',
                                                       'GAGC|CGGA'                             
                             ]))
    
#BA3M 	Generate All Maximal Non-Branching Paths in a Graph 
        def test_ba3m(self):   
            paths=non_branching_paths({
                1 :[2],
                2 : [3],
                3 : [4,5],
                6 : [7],
                7 : [6]})
            self.assertIn([1 ,2 ,3],paths)
            self.assertIn([3, 4],paths)
            self.assertIn([3 ,5],paths)
            self.assertIn([6,7,6 ],paths) 
            self.assertEqual(4,len(paths))

# BA4B	Find Substrings of a Genome Encoding a Given Amino Acid String
        def test_ba4b(self):
            encodings=findEncodings('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA','MA')
            self.assertEqual(3,len(encodings))
            self.assertIn('ATGGCC',encodings)
            self.assertIn('GGCCAT',encodings)
            self.assertIn('ATGGCC',encodings)

# BA4C	Generate the Theoretical Spectrum of a Cyclic Peptide
        def test_ba4c(self):
            self.assertEqual([0,113,114,128,129,227,242,242,257,355,356,370,371,484],
                             cycloSpectrum('LEQN'))
        
        def test_ba4d(self):
            self.assertEqual(14712706211,count_peptides_linear(1024))
        
        def test_chain(self):
            self.assertAlmostEqual(821.392,get_weight('SKADYEK'),places=3)
        

                
        def test_ba4e(self):
            seq=find_cyclopeptide_sequence([0,113,128,186,241,299,314,427])
            self.assertEqual(6,len(seq))
            self.assertIn([186,128,113],seq)
            self.assertIn([186,113,128],seq)
            self.assertIn([128,186,113],seq)
            self.assertIn([128,113,186],seq)
            self.assertIn([113,186,128],seq)
            self.assertIn([113,128,186],seq)
 
# BA4F 	Compute the Score of a Cyclic Peptide Against a Spectrum 
        def test_ba4f(self):
            self.assertEqual(
                11,
                score(
                    'NQEL',
                    [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]))

# BA4G 	Implement LeaderboardCyclopeptideSequencing 
        def test_ba4g(self):
            self.assertEqual([129, 71, 147, 113],
                             leaderPeptide(
                                 10,
                                 [0, 71, 113, 129, 147, 200, 218,\
                                  260, 313, 331, 347, 389, 460]))            

# BA4H 	Generate the Convolution of a Spectrum 
        def test_ba4h(self):
            self.assertEqual(
                [137, 137, 186, 186, 49, 323],
                convolution_expanded([0, 137, 186, 323]))
            
        #def test_random(self):
            #self.assertEqual([-5.737, -5.217, -5.263, -5.360, -5.958, -6.628, -7.009],
                             #random_genome('CGATACAA',[0.129, 0.287, 0.423, 0.476, 0.641, 0.742, 0.783]))


# BA4I 	Implement ConvolutionCyclopeptideSequencing 
        def test_ba4i(self):
            self.assertEqual([99,71,137,57,72,57],
                             convolutionCyclopeptideSequencing(
                                 20,
                                 60,
                                 [57, 57, 71, 99, 129, 137,\
                                  170, 186, 194, 208, 228,\
                                  265, 285, 299, 307, 323,\
                                  356, 364, 394, 422, 493])) 
            
# BA4J 	Generate the Theoretical Spectrum of a Linear Peptide 
        def test_ba4j(self):
            self.assertEqual(
                [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484],
                linearSpectrumFromString('NQEL'))     

# BA4L 	Trim a Peptide Leaderboard
        def test_ba4l(self):
            self.assertEqual(['LAST', 'ALST'],
                             trim_for_strings(['LAST', 'ALST', 'TLLT', 'TQAS'],
                                  [0,71,87,101,113,158,184,188,259,271,372],
                                  2))

# BA4M 	Solve the Turnpike Problem

        def test_ba4m(self):        
            self.assertEqual([0, 2, 4,7, 10],
                             turnpike([
                                 -10,-8, -7, -6, -5,
                                 -4, -3, -3, -2, -2,
                                 0, 0, 0, 0, 0, 
                                 2, 2, 3, 3, 4,
                                 5, 6, 7, 8, 10
                             ]))

# BA5A 	Find the Minimum Number of Coins Needed to Make Change 	
        def test_ba5a(self):
            self.assertEqual(2,number_of_coins(40,[1,5,10,20,25,50]))
            self.assertEqual(338,number_of_coins(8074,[24,13,12,7,5,3,1]))
            
# BA5B 	Find the Length of a Longest Path in a Manhattan-like Grid  	
        def test_ba5b(self):
            self.assertEqual(34,\
                             longest_manhattan_path(4,\
                                                    4,\
                                                    [[1, 0, 2, 4, 3],\
                                                     [4, 6, 5, 2, 1],\
                                                     [4, 4, 5, 2, 1],\
                                                     [5, 6, 8, 5, 3]],\
                                                    [[3, 2, 4, 0],\
                                                     [3, 2, 4, 2],\
                                                     [0, 7, 3, 3],\
                                                     [3, 3, 0, 2],\
                                                     [1, 3, 2, 2]]))
            self.assertEqual(84,\
                             longest_manhattan_path(17,\
                                                    9,\
                                                    [[2,3,4,0,3,1,1,1,1,1],
                                                     [4,2,3,4,3,3,0,4,1,1],\
                                                     [4,4,0,1,4,3,2,0,2,2],\
                                                     [4,3,0,3,4,4,3,2,4,4],\
                                                     [0,1,0,1,3,0,3,0,3,4],\
                                                     [3,2,4,4,4,3,1,0,0,0],\
                                                     [3,4,3,1,2,3,0,0,4,0],\
                                                     [2,4,3,4,1,2,0,3,2,0],\
                                                     [1,4,4,1,4,4,3,1,1,4],\
                                                     [3,1,2,2,3,3,0,4,0,0],\
                                                     [0,2,1,4,1,3,1,3,1,0],\
                                                     [1,0,4,0,4,3,3,2,3,1],\
                                                     [2,0,0,4,3,4,0,3,0,0],\
                                                     [4,1,0,4,3,2,1,1,3,1],\
                                                     [2,4,4,3,3,4,0,0,4,3],\
                                                     [1,0,2,3,3,0,4,0,2,0],\
                                                     [3,1,0,3,2,3,2,2,1,4]],\
                                                    [[1,0,4,4,3,3,1,0,4],\
                                                     [0,2,0,3,3,0,1,2,1],\
                                                     [3,2,3,1,1,4,2,4,4],\
                                                     [1,3,4,4,2,1,1,1,4],\
                                                     [1,4,2,2,3,1,3,2,3],\
                                                     [0,3,1,0,1,0,4,1,4],\
                                                     [1,3,4,4,1,0,3,2,1],\
                                                     [2,3,1,2,3,2,2,2,3],\
                                                     [3,2,1,4,0,2,4,2,4],\
                                                     [4,0,2,0,1,3,1,4,4],\
                                                     [1,3,0,2,2,1,0,3,2],\
                                                     [1,4,0,4,4,1,2,4,2],\
                                                     [0,2,4,3,4,0,3,2,2],\
                                                     [2,3,4,4,0,4,3,0,4],\
                                                     [1,0,4,1,3,3,1,4,2],\
                                                     [4,3,4,3,2,3,2,2,0],\
                                                     [0,1,2,2,4,4,2,4,2],\
                                                     [2,3,1,4,4,3,4,0,3]]))            
                         
# BA5C 	Find a Longest Common Subsequence of Two Strings  	
        def test_ba5c(self):
            self.assertEqual('ACCTTG',
                             longest_common_subsequence('AACCTTGG',
                                                        'ACACTGTGA'))            
# BA5D 	Find the Longest Path in a DAG  	
        def test_ba5d(self):
            n,path=longest_path(0,4,
                                [
                                    (0,1,7),
                                    (0,2,4),
                                    (2,3,2),
                                    (1,4,1),
                                    (3,4,3)
                                ])
            self.assertEqual(9,n)
            self.assertEqual([0,2,3,4],path)
            n,path=longest_path(11,18,
                                [(21,33,29),
                                 (5,19,8),
                                 (11,17,25),
                                 (10,35,33),
                                 (14,27,19),
                                 (6,27,20),
                                 (6,20,32),
                                 (14,23,4),
                                 (7,8,29),
                                 (12,18,9),
                                 (20,27,16),
                                 (15,35,39),
                                 (1,27,2),
                                 (9,19,37),
                                 (12,14,20),
                                 (30,34,9),
                                 (2,9,5),
                                 (24,31,37),
                                 (8,18,33),
                                 (8,13,20),
                                 (17,23,20),
                                 (6,34,13),
                                 (9,29,20),
                                 (16,29,12),
                                 (1,22,8),
                                 (2,32,27),
                                 (13,26,25),
                                 (3,21,2),
                                 (13,20,22),
                                 (3,23,20),
                                 (8,9,13),
                                 (12,16,30),
                                 (3,29,22),
                                 (1,3,2),
                                 (19,35,7),
                                 (14,19,36),
                                 (10,21,38),
                                 (19,30,38),
                                 (23,32,3),
                                 (5,6,19),
                                 (9,33,2),
                                 (6,10,17),
                                 (21,24,8),
                                 (9,15,36),
                                 (10,27,17),
                                 (25,33,20),
                                 (9,11,1),
                                 (0,10,39),
                                 (6,31,39),
                                 (6,16,24),
                                 (0,6,34),
                                 (14,15,30),
                                 (1,7,0),
                                 (0,22,3),
                                 (4,30,3),
                                 (4,15,20),
                                 (0,21,18),
                                 (0,26,28),
                                 (22,27,20),
                                 (3,32,32),
                                 (3,33,10),
                                 (5,28,13),
                                 (17,32,33),
                                 (7,13,25),
                                 (17,18,37),
                                 (8,31,33),
                                 (7,18,20),
                                 (20,35,20)
                                 ])
                                            
            self.assertEqual(62,n)
            self.assertEqual([11,17,18],path)  
            
# BA5E 	Find a Highest-Scoring Alignment of Two Strings  	
        def test_ba5e(self):
            score,s1,s2=highest_scoring_global_alignment('PLEASANTLY','MEANLY')
            self.assertEqual(8,score)
            self.assertEqual('PLEASANTLY',s1)
            self.assertEqual('-MEA--N-LY',s2)
            
# BA5F 	Find a Highest-Scoring Local Alignment of Two Strings 
        def test_ba5f(self):
            score,s1,s2=highest_scoring_local_alignment('MEANLY','PENALTY')
            self.assertEqual(15,score)
            self.assertEqual('EANL-Y',s1)
            self.assertEqual('ENALTY',s2)
            score,s1,s2=highest_scoring_local_alignment(\
                'AMTAFRYRQGNPRYVKHFAYEIRLSHIWLLTQMPWEFVMGIKMPEDVFQHWRVYSVCTAEPMRSDETYEQKPKPMAKWSGMTIMYQAGIIRQPPRGDRGVSDRNYSQCGKQNQAQLDNNPTWTKYEIEWRVQILPPGAGVFEGDNGQNQCLCPNWAWEQPCQWGALHSNEQYPNRIHLWAPMSKLHIKIEKSSYNRNAQFPNRCMYECEFPSYREQVDSCHYENVQIAFTIFSGAEQKRKFCSCHFWSNFIDQAVFSTGLIPWCYRRDDHSAFFMPNWNKQYKHPQLQFRVAGEGTQCRPFYTREMFTKVSAWRIAGRFAGPYERHHDAHLELWYQHHKVRTGQQLGIIWNNRDKTRNPCPFSAYYNKLPWWKINQNAFYNCLQNIAHSTHDETHEFNPVKCIDWLQGTMVPTECKKGFVHEKCECYRNPGPPLHDMYHQMEDIFGVRFDCLTGWKHLSDYNPCQERRNINDFYIFAYEIAPAVKNLVLSPQPLADATKKCAFNYTPLDQSPVVIACKWYIHQPICMLLIVLICAMDKYNAHMIVIRTTEGQQPMHACRMTEGPGMCMKEPLVTFTLPAQWQWPNHEFKYVYMYVLNYHLSQYTYTDEGHAGGQHYSFNVAVDVGMAWGHNRCYCQPACYSQQETQTRTIDYEKWQYMKHQAFKWGLWFCEQERHAWFKGQNRCEMFTAKMTRMGADSNLDQYKLMLAQNYEEQWEQPIMECGMSEIIEIDPPYRSELIFTFWPFCTYSPWQNLIKCRCNNVIEEMDQCVPLTFIGFGVKQAGGIQAWAFYKEEWTSTYYLMCQCMKSDKAQYPYEIILFWMQPMDTGEQEPPQQNMWIFLPHSWFFDWCCNAPWSEICSSRHDHGQCQDAFYPCELFTVFDDIFTAEPVVCSCFYDDPM',\
                'WQEKAVDGTVPSRHQYREKEDRQGNEIGKEFRRGPQVCEYSCNSHSCGWMPIFCIVCMSYVAFYCGLEYPMSRKTAKSQFIEWCDWFCFNHWTNWAPLSIVRTSVAFAVWGHCWYPCGGVCKTNRCKDDFCGRWRKALFAEGPRDWKCCKNDLQNWNPQYSQGTRNTKRMVATTNQTMIEWKQSHIFETWLFCHVIIEYNWSAFWMWMNRNEAFNSIIKSGYPKLLLTQYPLSQGSTPIVKPLIRRDQGKFWAWAQMWWFREPTNIPTADYCHSWWQSRADLQNDRDMGPEADASFYVEFWYWVRCAARTYGQQLGIIWNNRLKTRNPCPYSADGIQNKENYVFWWKNMCTKSHIAFYYCLQNVAHYTHDVTAEFNPVKCIDWLQGHMVLSSWFKYNTECKKLFVHEKCECYRMFCGVVEDIFGVRFHTGWKHLSTAKPVPHVCVYNPSVQERRNINDFYIFYEIAPAVKNLVLSAQPLHDYTKKCAFNYTPITITRIISTRNQIIWAHVVIACQFYSPHQMLLIELAMDKYCADMNVRRSTEGHQPMHACRSTFGPGMAAKEPLVTFTLVAFWQWPNHEFQYVYMYTEDKIIQIGPHLSNGCEMVEYCVDCYAKRPCYRAYSAEAQYWRMITEAEDYSYKTRNAIAATATVRGQYCHPFRWLGIVWMAHHDCFFANECGTICIPQMAEMRPPETTPYEIDIIFMMFWKEHMSTTILDVVGMYRPATFSHWHDAHHQCEPYLTPLMCQSKLVFDAAFTQVGVKGVWYHTEKLELMAGFNHMKFKKEEAQQSCFYWFQDCPDYDPPDAVRKTDEKHIRAHGEIWWLMRYYCMYHILHIASRHEWMHLRWDQACTNPGYELFEFIPWVLRRYVVYDKIRYNYSYRNSASMEFV')
            self.assertEqual(1062,score)
            self.maxDiff=None
            self.assertEqual('YQAGIIRQPPRGD-RGVSDRNYSQCGKQ-NQ-AQLDNNPTWTKYEIEWRVQI-LPPGAGVFEGDNGQNQCLCPNW--A-W-EQPCQW----GALHS-NEQYPNRIHLWAPMSKLHIKIEKSSYN-RNAQ-FPNRCMYECE-FPSY-REQVDSCHYENVQIAF-TIFSGAEQKRKFCSCHFWSNFIDQAVFSTGLI-PWCYRRDDHSAFFMPNWNKQ--YKHPQLQFRVAGEGTQCRPFYTREMFTKVSAWRIAGRFAGPYERHHDAHLELWY-QHHKVRT-GQQLGIIWNNRDKTRNPCPFSAY-Y-NK--LP-WWK-I-NQ-N-AFYNCLQNIAHSTHDETHEFNPVKCIDWLQGTMV-P------TECKKGFVHEKCECYRNPGPPLHDMYHQMEDIFGVRFDCLTGWKHLS------D---YNPC-QERRNINDFYIFAYEIAPAVKNLVLSPQPLADATKKCAFNYTPLDQSPVVIACK---WYIHQPI-CMLL----IVLIC-AMDKYNAHMIVIRTTEGQQPMHACRMTEGPGMCMKEPLVTFTLPAQWQWPNHEFKYVYMYVLNYHLSQYTYTDEGHAGGQHYSFNVAVDVGMAWGHNRCYCQPACYSQQETQTRTIDYEKWQYMKHQAFKWGLWFCEQER-HA--WFKGQNRCEMFTAKMTRMGADSNLDQYKLMLAQNYEEQWEQPIMECGMSEIIEIDPPYRSELIFTFWPFCTYSPWQNLIKCRCNNVIEEMDQCVP-LTF-IGFGVKQAGGIQA-WAFYKE--EWTSTYYLMCQCMKSDKAQYPYEIILFWMQ--P-MDTGE--QEPPQQNMWIFLPHSWFFDWCCNAPWSEICSSRHD--H---GQ-CQDAFYPCELFTVF',s1)
            self.assertEqual('Y-P-MSRKTAKSQFIEWCDW-F--CFNHWTNWAPLSIVRTSVAFAV-W-GHCWYPCG-GVCKTNRCKDD-FCGRWRKALFAEGPRDWKCCKNDLQNWNPQYSQGTR--NTK-RMVATTNQTMIEWKQSHIFETW-LF-CHVIIEYNWSAF-W-MWMNRNEAFNSIIKSGYPKLLL-T-QY-P-L-SQG--STPIVKPL-IRRD-QGKFW-A-WAQMWWFREPT-NIPTA-D-Y-CHSW--WQ--SR-ADLQ-NDRDMGP-EADASFYVEFWYWVRCAARTYGQQLGIIWNNRLKTRNPCPYSADGIQNKENYVFWWKNMCTKSHIAFYYCLQNVAHYTHDVTAEFNPVKCIDWLQGHMVLSSWFKYNTECKKLFVHEKCECYRM----FCGV---VEDIFGVRFH--TGWKHLSTAKPVPHVCVYNPSVQERRNINDFYIF-YEIAPAVKNLVLSAQPLHDYTKKCAFNYTPITITRIISTRNQIIW-AHVVIACQFYSPHQMLLIELAMDKYCADMNVRRSTEGHQPMHACRSTFGPGMAAKEPLVTFTLVAFWQWPNHEFQYVYMYTED-KIIQIG-PHLSN-GCEMVEYCVDC-YAK-RPCYRAYSAEAQYWRMITEAEDYSYKTRNAIAATATVRGQ-YCHPFRWLGIVWM-AHHDC-FFANECGTICI-PQMAEMRPPETTPYEI--DIIFMMF-WKE--HMSTTIL-DVVGMYRP-ATFSHWHDAHH-QCEPYLTPL-MCQSKLVFDAAFT--QVG-VKGVW-YHTEKLELMAGFNHM-K-FKKEEAQ---QSCFYWFQDCPDYDPPDAVRKTDEKHIRAHGEIWWLMRYYCMYHILHI-ASRHEWMHLRWDQACTNPGY--ELFE-F',s2)
                         
# BA5G 	Compute the Edit Distance Between Two Strings 	 	
                         
# BA5H 	Find a Highest-Scoring Fitting Alignment of Two Strings 	 	
                         
# BA5I 	Find a Highest-Scoring Overlap Alignment of Two Strings 	 	
                         
# BA5J 	Align Two Strings Using Affine Gap Penalties 	 	
                         
# BA5K 	Find a Middle Edge in an Alignment Graph in Linear Space 	 	
                         
# BA5L 	Align Two Strings Using Linear Space 	 	
                         
# BA5M 	Find a Highest-Scoring Multiple Sequence Alignment 	 	
                         
# BA5N 	Find a Topological Ordering of a DAG 
        def test_ba5n(self):
            self.assertEqual([5, 4, 1, 2, 3],
                             topological_order({
                                 1 : [2],
                                 2 : [3],
                                 4 : [2],
                                 5 : [3]
                             }))
            
#LGIS 	Longest Increasing Subsequence
        def test_longestIncreasingSubsequence(self):
            (a,d)=longestIncreasingSubsequence(5,[5, 1, 4, 2, 3])
            self.assertEqual([1,2,3],a)
            self.assertEqual([5,4,3],d)
            (a,d)=longestIncreasingSubsequence(9,[8, 2, 1, 6, 5, 7, 4, 3, 9])
            self.assertEqual([1, 5, 7, 9],a)
            self.assertEqual([8, 6, 5, 4, 3],d)
        
        def test_orf(self):
            string='''>Rosalind_99
            AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'''
            fasta=f.FastaContent(string.split('\n'))            
            peptides=get_reading_frames(fasta)
            self.assertIn('MLLGSFRLIPKETLIQVAGSSPCNLS',peptides)
            self.assertIn('M',peptides)
            self.assertIn('MGMTPRLGLESLLE',peptides)
            self.assertIn('MTPRLGLESLLE',peptides)            
            self.assertEqual(4,len(peptides))
            #for peptide in get_reading_frames(f.FastaFile('./data/rosalind_orf(3).txt') ):
                #print (peptide)

        #def test_wfmd(self):
            #self.assertAlmostEqual(0.772, wfmd(4,6,2,1),3)
            
        #def test_q5(self):
            #for peptide in ['MTAI',
                            #'MLAT',
                            #'MAIT',
                            #'MTAL',
                            #'TALM',
                            #'TLAM']:
                #sp=cycloSpectrum(peptide)
                #if sp==[0, 71, 101, 113, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416]:
                    #print (peptide, sp)
         
        #def test_q6(self):
            ## 0 71 99 101 103 128 129 199 200 204 227 230 231 298 303 328 330 332 333
            #spectrum=[0 ,71, 99, 101, 103, 128, 129, 199, 200, 204, 227, 230, 231,\
                #298, 303, 328, 330, 332, 333]
            #print (spectrum)
            #for peptide in ['TVQ',
                            #'CTV',
                            #'AVQ',
                            #'AQV',
                            #'TCE',
                            #'VAQ']:
                #masses=[rrt.integer_masses[a] for a in peptide]
                #masses.sort()
                #ls=rh.linearSpectrum(masses)
                #if rh.consistent(masses,spectrum):
                    #print (peptide)
                #else:
                    #print (peptide,masses,ls)
                    
    unittest.main(exit=False)
