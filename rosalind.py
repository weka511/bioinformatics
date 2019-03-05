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

import re
from helpers import translate

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