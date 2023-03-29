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

''' QRTD  Quartet Distance'''

from argparse           import ArgumentParser
from deprecated         import deprecated
from io                 import StringIO
from os.path            import basename
from time               import time
from Bio.Phylo          import read, draw
from Bio.Phylo.BaseTree import Clade, Tree
from matplotlib.pyplot  import figure, show
from helpers            import read_strings

class Component(Clade):
    def __init__(self,name):
        super().__init__(name=name)

class HTree(Tree):
    '''Tree associated with Hierarchical Decomposition'''
    def __init__(self,T):
        super().__init__()
        self.components  = {}
        self.leaves      = []
        for idx, clade in enumerate(T.find_clades()):
            if not clade.name:
                clade.name = str(idx)
            self.components[clade.name] = Component(clade.name)
            if clade.is_terminal():
                self.leaves.append(clade.name)
        edges_leaves = []
        edges_internal = []
        for a in T.find_clades():
            for b in a.clades:
                if b.is_terminal():
                    edges_leaves.append((a.name,b.name))
                else:
                    edges_internal.append((a.name,b.name))

        open_edges = edges_leaves + edges_internal
        self.agenda = []
        while len(open_edges)>0:
            skipped   = []
            closed     = []
            for a,b in open_edges:
                if not a in closed and not b in closed:
                    closed.append(a)
                    closed.append(b)
                    self.agenda.append((a,b))
                else:
                    skipped.append((a,b))

            open_edges = skipped.copy()

        z=0



def draw_tree(T,title,ax=None):
    draw(T,
        axes    = ax,
        do_show = False)
    ax.set_title(title)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])


@deprecated
def create_adj(tree,indices=None):
    adj = {}
    def dfs(tree):
        id       = tree['id']
        name     = tree['name']
        children = tree['children']
        parentid = tree['parentid']
        if len(name)==0:
            adj[id]=[]
        if parentid>-1:
            adj[parentid].append(id if len(name)==0 else indices[name])
        for child in children:
            dfs(child)
    dfs(tree)
    return adj
@deprecated
def create_edges(adj,n=None):
        return [(a,b) for a,children in adj.items() for b in children if b>=n]


@deprecated
def extract_quartets(edges,adj,n=None):

    # get_leaves
    #
    # Get all leaves in graph 'adj' that are below some specified node
    def get_leaves(node):
        # dfs
        #
        # Conduct depth first search, looking for leaves
        def dfs(u):
            if u>=n:
                for v in adj[u]:
                    dfs(v)
            else:
                leaves.append(u)

        leaves = []
        dfs(node)
        return leaves

    def split(a,b):
        s_b = leaves[b]
        s_a = [s for s in all_leaves if s not in s_b]
        return [(s_a[j],s_a[i],s_b[l],s_b[k]) for i in range(len(s_a))
                for j in range(i)
                for k in range(len(s_b))
                for l in range(k)]

    all_leaves     = list(range(n))
    leaves         = {x:sorted(get_leaves(x)) for x in adj.keys()}
    internal_edges = [(a,b) for a,b in edges if b>=n]
    return [q  for a,b in internal_edges for q in split(a,b)]
@deprecated
def get_matches(quartets1,quartets2):
    i       = 0
    j       = 0
    matches = 0
    while i < len(quartets1) and j < len(quartets2):
        if quartets1[i]==quartets2[j]:
            matches += 1
            i       += 1
            j       += 1
        elif quartets1[i]<quartets2[j]:
            i += 1
        else:  # quartets1[i]>quartets2[j]
            j+=1

    return matches

# qrtd
#
# Given: A list containing n taxa  and two unrooted binary trees T1 and T2 on the given taxa.
#        Both T1 and T2 are given in Newick format.
#
# Return: The quartet distance dq(T1,T2)
@deprecated
def qrtd(species,T1,T2):
    n         = len(species)
    indices   = {species[i]:i for i in range(n)}
    tree1     = parse(T1,start=n)
    adj1      = create_adj(tree1,indices=indices)
    edges1    = create_edges(adj1,n=n)
    quartets1 = extract_quartets(edges1,adj1,n=n)
    print (len(quartets1), quartets1)

    tree2     = parse(T2,start=n)
    adj2      = create_adj(tree2,indices=indices)
    edges2    = create_edges(adj2,n=n)
    quartets2 = extract_quartets(edges2,adj2,n=n)
    print (len(quartets2),quartets2)
    mismatches = 0
    for q in quartets1:
        if not q in quartets2:
            mismatches+=1
    for q in quartets2:
        if not q in quartets1:
            mismatches+=1
    return mismatches
    #return len(quartets1) + len(quartets2) - 2*get_matches(quartets1,quartets2)

if __name__=='__main__':
    start = time()
    parser = ArgumentParser(__doc__)
    parser.add_argument('--paper',    default=False, action='store_true', help='process dataset from paper')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    parser.add_argument('--show',     default=False, action='store_true', help='display plots')
    args = parser.parse_args()

    if args.paper:
        T1 = read(StringIO('a,((e,f),d),(b,c);'), 'newick')
        T2 = read(StringIO('a,((b,d),c),(e,f);'), 'newick')
        # HT1 = HTree(T1)
        HT2 = HTree(T2)

        fig = figure(figsize=(10,10))
        ax11 = fig.add_subplot(221)
        draw_tree(T1,'T1',ax11)

        ax12 = fig.add_subplot(222)
        draw_tree(T2,'T2',ax12)

        # ax21 = fig.add_subplot(223)
        # draw_tree(HT1,'HT1',ax21)

        ax22 = fig.add_subplot(224)
        draw_tree(HT2,HT2.agenda,ax22)


    if args.sample:
        print(qrtd(
            'A B C D E'.split(),
            '(A,C,((B,D),E));',
            '(C,(B,D),(A,E));'
        ))

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{basename(__file__).split(".")[0]}.txt')

        Result = qrtd(Input[0].split(), Input[1], Input[2])
        print (Result)
        with open(f'{basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')

    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()
