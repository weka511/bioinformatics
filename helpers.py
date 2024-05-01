#!/usr/bin/env python
# Copyright (C) 2015-2023 Simon Crase

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a cnp.argmaxy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

# This file contains a collection of functions to solve the problems
# at rosalind.info.

'''
Helper functions for Rosalind to facilitate testing
'''

from os.path import basename, join, expanduser
from re      import compile, match
from sys     import argv

from Bio     import SeqIO
import numpy as np


def read_strings(test_data,
                 init=1      # Set to zero to handle format used in extra dataset
                 ):
    '''
    read_strings

    Read an input file and return set of strings
    '''
    inputs   = []
    expected = []
    with open(test_data) as f:
        state = init
        for line in f:
            content=line.strip()
            if state==0:
                if content.startswith('Input'):
                    state=1
            elif state==1:
                if content.startswith('Output'):
                    state=2
                else:
                    inputs.append(content)
            else:
                expected.append(content)
    return inputs if init==1 else (inputs,expected)

def print_profile(profile):
    for key,value in profile.items():
        line=key
        sep=': '
        for v in value:
            line=line+sep
            line = line + str(v)
            sep = ' '
        print (line)


def print_peptide(seq):
    for p in seq:
        l=''
        for x in p:
            if len(l)>0:
                l=l+'-'
            l=l+str(x)
        print (l)




def read_list(file_name):
    with open(file_name) as f:
        return [line.strip() for line in f]

def print_adjacency_graph(graph):
    for (a,b) in graph:
        print ('%(a)s -> %(b)s'%locals())


def format_list(list):
    '''
    format_list

    Format a list of numbers for use as a Rosalind solution

    '''
    return ' '.join([str(l) for l in list])



def create_adjacency(edges,back=True,self=True):
    '''
    Create adjacency list, with both forward and backward links

     Inputs:   edges    Edges in graph
               back     Indicates whether to create back links (undirected graph)
               self     Indicates whether self links, a->a, are to be included

     Outputs:  Adjacency list, including loop a->a
    '''
    m,n       = edges[0]
    product   = {}

    for a in range(1,m+1):
        product[a]   = [a] if self else []

    for a,b in edges[1:]:
        product[a].append(b)
        if back:
            product[b].append(a)

    for a in range(1,m+1):
        product[a]=sorted(list(set(product[a])))

    return m,n,product


def print_adjacency_graph2(graph):
    for (a,b) in graph:
        succ=''
        for bb in b:
            if len(succ)>0:
                succ=succ+','
            succ=succ+bb
        print ('%(a)s -> %(succ)s'%locals())

def print_adjacency_graph3(graph):
    for (a,b) in graph.items():
        succ=''
        for bb in b:
            if len(succ)>0:
                succ=succ+','
            succ=succ+bb
        print ('%(a)s -> %(succ)s'%locals())


# Count occurences of a character in a string
#
# Input: a character and a string
#
# Return: Nmber of times character occurs in string

def count_occurences(c,string):
    count=0
    start=0
    while True:
        start=string.find(c,start)
        if start<0:
            return count
        count+=1
        start+=1



def print_strings(strings,sorted=False):
    if sorted:
        strings.sort()
    for s in strings:
        print (s)

# convert a number i to a list of binary bits of length n

def binary(i,n):
    dividend=i
    bits=[]
    for j in range(n):
        next_dividend=dividend//2
        bits.append(dividend-2*next_dividend)
        dividend=next_dividend
    return bits[::-1]



def flatten(x):
    '''flatten a list of lists into a single list'''
    return [inner for outer in x for inner in outer]

def print_matrix(B,out='./print_matrix.txt'):
    with open(out,'w') as f:
        for row in B:
            line=""
            for x in row:
                if len(line)>0:
                    line=line+' '
                line=line+"{:15.12f}".format(x)
            f.write(line+'\n')


class Newick(object):
    '''
    Newick -- Represent a tree in Newlick format

    Tree --> Subtree ";" | Branch ";"
    Subtree --> Leaf | Internal
    Leaf --> Name
    Internal --> "(" BranchSet ")" Name
    BranchSet --> Branch | Branch "," BranchSet
    Branch --> Subtree Length
    Name --> empty | string
    Length --> empty | ":" number
    '''
    @staticmethod
    def find_commas(string):
        commas=[]
        bracket_level=0
        for i in range(len(string)):
            if string[i]==',' and bracket_level==0:
                commas.append(i)
            if string[i]=='(':
                bracket_level+=1
            if string[i]==')':
                bracket_level-=1
        return commas

    @staticmethod
    def find_balanced_brackets(string):
        bracket_level=0
        for i in range(len(string)):
            if string[i]=='(':
                bracket_level+=1
            if string[i]==')':
                bracket_level-=1
                if bracket_level==0:
                    return i

class Tree(Newick):
    def __init__(self,subtree):
        self.subtree=subtree
    def __str__(self):
        return f'{self.subtree};'

    @staticmethod
    def parse(string):
        s=string.replace(' ','')
        branch=Branch.parse(s[0:-1])
        if branch and s[-1]==';':
            return Tree(branch)
        subtree=SubTree.parse(s[0:-1])
        if subtree and s[-1]==';':
            return Tree(subtree)

class SubTree(Newick):
    def __init__(self,subtreeOrInternal):
        self.subtreeOrInternal=subtreeOrInternal
    def __str__(self):
        return str(self.subtreeOrInternal)
    @staticmethod
    def parse(string):
        internal=Internal.parse(string)
        if internal:
            return SubTree(internal)
        leaf=Leaf.parse(string)
        if leaf:
            return SubTree(leaf)

class Internal(Newick):
    def __init__(self,branchSet,name):
        self.branchSet=branchSet
        self.name=name
    def __str__(self):
        return '{}{}'.format(str(self.branchSet), str(self.name))
    @staticmethod
    def parse(string):
        if len(string)>1 and string[0]=='(':
            rh=Newick.find_balanced_brackets(string)
            internal_string=string[1:rh]
            remainder=string[rh+1:].strip()
            branch_set=BranchSet.parse(internal_string)
            if branch_set:
                name=Name.parse(remainder)
                return Internal(branch_set,name)

class BranchSet(Newick):
    def __init__(self,branches):
        self.branches=[b for b in branches]
    def __str__(self):
        return '({})'.format(','.join([str(b) for b in self.branches]))
    @staticmethod
    def parse(string):
        commas=Newick.find_commas(string)
        if len(commas)==0:
            return BranchSet([Branch.parse(string)])
        start=0
        branch_set=[]
        for pos in commas:
            branch=Branch.parse(string[start:pos])
            if branch:
                branch_set.append(branch)
            else:
                return None
            start=pos+1
        branch=Branch.parse(string[start:])
        if branch:
            branch_set.append(branch)
            return BranchSet(branch_set)

class Branch(Newick):
    def __init__(self,subtree,length=None):
        self.subtree=subtree
        self.length=length
    def __str__(self):
        return str(self.subtree)                         \
               if self.length==None or self.length.number==None else  \
               '{}:{}'.format(self.subtree, self.length)
    @staticmethod
    def parse(string):
        parts=string.split(':')
        if len(parts)>1:
            return Branch(SubTree.parse(':'.join(parts[0:-1])),Length(parts[-1]))
        return Branch(SubTree.parse(string))

class Leaf(Newick):
    def __init__(self,name):
        self.name=name
    def __str__(self):
        return str(self.name)
    @staticmethod
    def parse(string):
        return Leaf(Name.parse(string))

class Name(object):
    def __init__(self,name=None):
        self.name=name
    def __str__(self):
        return '' if self.name==None else '{}'.format(str(self.name))
    @staticmethod
    def parse(string):
        if len(string)==0:
            return Name()
        matches=match('[A-Za-z_]*',string)
        if matches:
            return Name(matches.group(0))

class Length(object):
    def __init__(self,number=None):
        self.number=number
    def parse(string):
        return None if number==None else float(number)
    def __str__(self):
        return '' if self.number==None else str(self.number)


def parse_graph(f):
    product = []
    for line in f:
        parts =line.strip().split()
        product.append((int(parts[0]),int(parts[1])))
    return product

def parse_graphs(f):
    product = []
    graph   = []
    state   = 0
    n       = -1
    for line in f:
        ll =line.strip()
        if state==0:
            if len(ll)>0:
                n = int(ll)
                state = 1
                continue
        if state==1:
            if len(ll)==0:
                state = 2
                graph = []
                continue
        if state == 2:
            if len(ll)==0:
                product.append(graph)
                graph=[]
            else:
                lll =ll.split()
                graph.append((int(lll[0]),int(lll[1])))
    product.append(graph)
    assert len(product)==n,'{0} {1}'.format(len(product),n)
    for graph in product:
        a,b=graph[0]
        assert len(graph)==b+1,'{0} {1}'.format(len(graph),b)
    return product



def create_strings(problem = basename(argv[0]).split('.')[0],
                   path    = join(expanduser('~'),'Downloads'),
                   ext     = None,
                   fasta   = False,
                   name    = None):
    '''
    create_strings
    Read test data from file

    Inputs:   problem  Name of Rosalind problem
                       - defaults to name of script callling this method
              path     Location of test data
              ext      Ext used to identify additional datasets - 1, 2, 3...
              fasta    File is in FASTA format, so skip FASTA ids

    Returns: A list of strings, one for each row
    '''
    def get_file_name():
        label        = problem if ext==None else '{0}({1})'.format(problem,ext)
        base         = 'rosalind_{0}.txt'.format(label) if name==None else name+'.txt'
        return join(path,base)

    product      = []
    with open(get_file_name(),'r') as f:
        if fasta:
            for record in SeqIO.parse(f,'fasta'):
                product.append(str(record.seq))
        else:
            for line in f:
                product.append(line.strip())
    return product

def create_list(problem  = basename(argv[0]).split('.')[0],
                path     = join(expanduser('~'),'Downloads'),
                ext      = None,
                fasta    = False,
                name     = None):
    '''
    create_list

    Read test data from file

    Inputs:   problem  Name of Rosalind problem
                       - defaults to name of script callling this method
              path     Location of test data
              ext      Ext used to identify additional datasets - 1, 2, 3...
              fasta    File is in FASTA format, so skip FASTA ids

    Returns: A list of lists, one for each row
    '''
    product = []
    for row in create_strings(problem = problem,
                              path    = path,
                              ext     = ext,
                              fasta   = fasta,
                              name    = name):
        if len(row)>0:
            product.append([int(s) for s in row.split(" ")])

    return product



def create_weighted_adjacency_list(problem=basename(argv[0]).split('.')[0],
                                   path=join(expanduser('~'),'Downloads'),
                                   ext=None,
                                   name=None):
    '''
    create_weighted_adjacency_list

    Read test data from file

    Inputs:   problem  Name of Rosalind problem
                       - defaults to name of script callling this method
              path     Location of test data
              ext      Ext used to identify additional datasets - 1, 2, 3...
              fasta    File is in FASTA format, so skip FASTA ids

    Returns: A weighted adjacency list - e.g. http://rosalind.info/problems/ba7a/
    '''
    n = -1
    T = {}
    p = compile('([0-9]+)->([0-9]+):([.0-9]+)')

    for line in create_strings(problem=problem, path=path, ext=ext, name=name):
        if n == -1:
            n = int(line)
        else:
            m = p.match(line)
            if  not int(m.group(1)) in T:
                T[int(m.group(1))] = [(int(m.group(2)),int(m.group(3)))]
            else:
                T[int(m.group(1))].append((int(m.group(2)),int(m.group(3))))
    return n,T

def read_matrix(problem=basename(argv[0]).split('.')[0],
                path=join(expanduser('~'),'Downloads'),
                ext=None,
                conv=int,
                len_params=1,
                name=None):
    '''
     read_matrix

     Read test data from file

     Inputs:   problem    Name of Rosalind problem
                          - defaults to name of script callling this method
               path       Location of test data
               ext        Ext used to identify additional datasets - 1, 2, 3...
               fasta      File is in FASTA format, so skip FASTA ids
               len_params Number of parameters that precede matrix

     Returns:  Paramters, followed by matrix - e.g. http://rosalind.info/problems/ba7b/
    '''
    params=[]
    D=[]

    for line in create_strings(problem=problem,path=path,ext=ext,name=name):
        if len(params)<len_params:
            params.append(int(line))
        else:
            D.append([conv(s) for s in line.split()])

    return (params,D)

def create_hmm(problem=basename(argv[0]).split('.')[0],
               path = join(expanduser('~'),'Downloads'),
               ext  = None,
               name = None):
    '''
    Create arguments for a Hidden Markov Model from strings read from a file
    '''
    return create_hmm_from_strings(create_strings(problem=problem,path=path,ext=ext,name=name))

def create_hmm_from_strings(strings,sep=' '):
    xs         = strings[0]
    alphabet   = strings[2].replace(sep,'')
    States     = strings[4].replace(sep,'')
    n          = len(States)
    m          = len(alphabet)
    Transition = np.zeros((n,n))
    i          = 0
    while i<n:
        items = strings[7+i].split()
        for j in range(1,len(items)):
            Transition[i,j-1] = float(items[j])
        i+=1
    i+=9

    Emission   = np.zeros((n,m))
    k          = 0
    while i<len(strings):
        items = strings[i].split()
        for j in range(1,len(items)):
            Emission[k,j-1] = float(items[j])
        i+=1
        k+=1
    return (xs,alphabet,States,Transition,Emission)



def expand(s):
    '''
    expand

    Convert a string of digits to an integer list
    '''
    return [int(c) for c in s]
