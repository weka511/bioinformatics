#!/usr/bin/env python
# Copyright (C) 2015-2024 Simon Crase

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

'''
Helper functions for Rosalind to facilitate testing with data from rosalind.info
'''

from os.path import basename, join, expanduser
from re      import compile, match
from sys     import argv
import numpy as np
from Bio     import SeqIO

def read_strings(test_data,
                 init=1      # Set to zero to handle format used in extra dataset
                 ):
    '''
    read_strings

    Read an input file and return set of strings
    '''
    inputs = []
    expected = []
    with open(test_data) as f:
        state = init
        for line in f:
            content = line.strip()
            if state == 0:
                if content.startswith('Input'):
                    state = 1
            elif state == 1:
                if content.startswith('Output'):
                    state = 2
                else:
                    inputs.append(content)
            else:
                expected.append(content)
    return inputs if init == 1 else (inputs,expected)



def read_list(file_name):
    '''
    Used to read a list of items from an input file
    '''
    with open(file_name) as f:
        return [line.strip() for line in f]


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

def flatten(x):
    '''flatten a list of lists into a single list'''
    return [inner for outer in x for inner in outer]

def parse_graph(f):
    '''Used to read a graph from a file'''
    product = []
    for line in f:
        parts = line.strip().split()
        product.append((int(parts[0]),int(parts[1])))
    return product

def parse_graphs(f):
    '''Used to read multiple graphs from a file'''
    product = []
    graph   = []
    state   = 0
    n       = -1
    for line in f:
        ll = line.strip()
        if state == 0:
            if len(ll) > 0:
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
