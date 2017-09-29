'''
 Rosalind utilities

 Copyright (C) 2017 Greenweaves Software Pty Ltd

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

class Tree(object):
    '''
    Undirected, weighted tree
    '''
    def __init__(self,N=-1):
        self.nodes=list(range(N))
        self.edges={}
        
    def link(self,start,end,weight=1): 
        self.half_link(start,end,weight)
        self.half_link(end,start,weight)
        
    def unlink(self,i,k):
        try:
            self.half_unlink(i,k)
            self.half_unlink(k,i)
        except KeyError:
            print ('Could not unlink {0} from {1}'.format(i,k))
            self.print()        
        
    def half_link(self,a,b,weight=1):
        if not a in self.nodes:
            self.nodes.append(a)        
        if a in self.edges:
            self.edges[a].append((b,weight))
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
            
class LabelledTree(Tree):
        
    @staticmethod
    def parse(N,lines,letters='ATGC'):
        T=LabelledTree(4)
        pattern=re.compile('([0-9]+)->(([0-9]+)|([{0}]+))'.format(letters))
        
        for line in lines:
            m=pattern.match(line)
            if m:
                i=int(m.group(1))
                if m.group(3)==None:
                    j=len(T.leaves)
                    T.half_link(i,j)
                    T.leaves.append(j)
                    T.labels[j]=m.group(4)
                if m.group(4)==None:
                    j=int(m.group(3))
                    T.half_link(i,j)
        return T
    
    def __init__(self,N=-1):
        super().__init__(N)
        self.labels={}
        self.leaves=[]
        
    def is_leaf(self,v):
        return v in self.leaves
    
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