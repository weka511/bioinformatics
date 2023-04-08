#!/usr/bin/env python

#   Copyright (C) 2020-2023 Simon Crase, simon@greenweaves.nz

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

from   argparse           import ArgumentParser
from   deprecated         import deprecated
from   io                 import StringIO
from   os.path            import basename
from   time               import time
from   matplotlib.pyplot  import figure, show
import numpy              as     np
from   scipy.special      import binom
from   helpers            import read_strings
from   phylogeny          import qrtd

# class Colour:
    # '''The 3 colours mentioned in the paper'''
    # A = 0
    # B = 1
    # C = 2
    # n = 3
    # @staticmethod
    # def as_text(c):
        # '''
        # Used for determing the colour of leaves in plots
        # '''
        # if   c==Colour.A: return 'xkcd:red'
        # elif c==Colour.B: return 'xkcd:green'
        # elif c==Colour.C: return 'xkcd:blue'
        # else: return None

# class Species:
    # '''
    # The represents one Species; thr two trees are over the same set of species.
    # '''
    # def __init__(self,name=None,colour=None):
        # self.name   = name
        # self.colour = colour



# def F3(_tuple):
    # '''
    # F3

    # This is the F function described in Lemma 5, for a component consisting of a single degree 3 node.

    # '''
    # def generate_cycles(n=3):
        # basic_cycle = np.array(list(range(n)))
        # for i in range(n):
            # yield np.roll(basic_cycle,i)

    # def get_tuple_index(index,letter):
        # return 3*index+letter-1

    # def get_tuple_value(index,letter):
        # return _tuple[get_tuple_index(index,letter)]

    # def get_term(index,letter):
        # return binom(get_tuple_value(index[0],letter[0]),2)                                          *\
                        # (get_tuple_value(index[2],letter[1])*get_tuple_value(index[1],letter[2])     + \
                         # get_tuple_value(index[1],letter[1])*get_tuple_value(index[2],letter[2]))

    # return sum([get_term(index,letter) for index in generate_cycles() for letter in generate_cycles()])

# class ComponentClade(Clade):
    # '''
    # ComponentClade

    # This class represents a clade in the Hierarchical Decomposition Tree, and also
    # represents a component in the component decompostion
    # '''
    # def __init__(self,
                # name             = '',
                # clades           = [],
                # branch_length    = 1,
                # composition_type = None):
        # super().__init__(name        = name,
                       # clades        = clades,
                       # branch_length = branch_length)
        # self.composition_type = composition_type

    # @staticmethod
    # def get_label(clade):
        # '''
        # get_label

        # Label nodes on HT tree
        # '''
        # if clade.is_terminal():
            # return clade.name
        # else:
            # return ''.join(['i' for i in range(clade.composition_type)]) if clade.composition_type<4 else 'iv'



    # def decorate(self,S):

        # def get_element(species,target):
            # 1 if species.colour == target else 0

        # if self.composition_type==None:
            # if self.name in S:
                # species    = S[self.name]
                # self._tuple = np.array([get_element(species,colour) for colour in range(Colour.n)])
                # self.F     = lambda _: 0
            # else:
                # self._tuple = np.zeros((3))
                # self.F     = F3
        # else:
            # self._tuple = sum([child._tuple for child in self.clades])
            # if self.composition_type==1:
                # self.F = lambda a,b,c: 0      #TODO
            # elif self.composition_type==2:
                # self.F = lambda a,b,c: 0      #TODO
            # elif self.composition_type==3:
                # self.F = lambda a,b,c: 0      #TODO
            # else:
                # assert self.composition_type==4
                # self.F = lambda a,b,c: 0      #TODO

# class HTree(Tree):
    # def __init__(self, root=None, rooted=True, id=None, name=None):
        # super().__init__(root=root, rooted=rooted, id=id, name=name)

    # def decorate(self,S):
        # for clade in self.find_elements(order='postorder'):
            # if hasattr(clade,'composition_type'):
                # clade.decorate(S)
            # else:
                # print ('skip', clade)

    # def get_node_count(self,v):
        # pass

# class HTreeBuilder(Tree):
    # '''This class constructs the Tree associated with Hierarchical Decomposition'''

    # def build(self,T):
        # '''Build the tree'''
        # self.components  = {}
        # self.T           = T
        # self.top         = {}

        # self.create_initial_components()
        # open_edges = self.collect_edges()
        # while len(open_edges)>0:
            # open_edges = self.extract_open_edges(open_edges)

        # return HTree(self.get_root())


    # def create_initial_components(self):
        # '''
        # create_initial_components

        # Ensure all nodes have names, and wrap each node in a component
        # '''
        # for idx, clade in enumerate(self.T.find_clades()):
            # if not clade.name:
                # clade.name = str(idx)
            # self.components[clade.name] = ComponentClade(name=clade.name)


    # def collect_edges(self):
        # '''
        # collect_edges

        # Collect edges for first pass; prioritize edges that terminate in leaves.
        # '''
        # edges_leaves   = []
        # edges_internal = []
        # for a in self.T.find_clades():
            # for b in a.clades:
                # if b.is_terminal():
                    # edges_leaves.append((a.name,b.name))
                # else:
                    # edges_internal.append((a.name,b.name))

        # return edges_leaves + edges_internal

    # def get_type(self,key,clades):
        # '''
        # get_type

        # Assign a type to a node in tree
        # '''
        # if len(clades[0].clades)==0 and len(clades[1].clades)==0:              return 1
        # if all([hasattr(clade,'type') and clade.type==1 for clade in clades]): return 2
        # return 3

    # def extract_open_edges(self,open_edges):
        # '''
        # extract_open_edges

        # This is the heart build(...). It processes as many edges as possible, adding them
        # to the clade structure. The remaining edges are modified to point to the top of the
        # constructed clade that contains the original nodes.
        # '''
        # remaining_edges = []
        # for a,b in open_edges:
            # if a not in self.top and b not in self.top:
                # key = f'{a}-{b}'
                # self.top[a]          = key
                # self.top[b]          = key
                # clades               = [self.components[a],self.components[b]]
                # new_component        = ComponentClade(name           = key,
                                                    # clades           = clades,
                                                    # branch_length    = 1,
                                                    # composition_type = self.get_type(key,clades))

                # self.components[key] = new_component
            # else:
                # if a in self.top:
                    # a = self.top[a]
                # if b in self.top:
                    # b = self.top[b]
                # remaining_edges.append((a,b))
        # return remaining_edges

    # def get_root(self):
        # '''
        # get_root

        # Choose clade that will define the root node for tree
        # '''
        # names                  = [k for k,v in self.components.items()]
        # index_root             = np.argmax([len(name) for name in names])
        # root                   = self.components[names[index_root]]
        # root.composition_type  = 4
        # return root

# def ensure_all_nodes_have_names(tree):
    # '''
    # ensure_all_nodes_have_names

    # Attach a name to each node that doesn't alread have one
    # '''
    # for idx, clade in enumerate(tree.find_clades()):
        # if not clade.name:
            # clade.name = str(idx)


# def root_with_specified_leaf(T,S,index=0):
    # '''
    # Allow user to change root to some arbitrary leaf - see Figure 9
    # '''
    # leaves        = T.get_terminals()
    # T.root_with_outgroup([leaves[index]])
    # ensure_all_nodes_have_names(T)
    # root_name           = leaves[index].name
    # T.root.clades       = [clade for clade in T.root.clades if clade.name!=root_name]
    # T.root.name         = root_name
    # return root_name

# def link(T1,HT2,S):
    # '''
    # Create the links shown in Figure 9

    # Currently the link is via dictionary lookup on clade names.
    # '''
    # def link1(T):
        # return {clade.name:clade for clade in T.find_clades() if clade.name in S}
    # return link1(T1),link1(HT2)

# def cache_sizes(Tr1):
    # '''
    # cache_sizes

    # Used by Count(...) to determine which of two subtrees is larger or amaller

    # We process clades bottom up so we can use values that were calculated previously
    # '''
    # size   = {}
    # clades = list(Tr1.find_elements(order='postorder'))
    # for clade in clades[:-1]:
        # size[clade.name] = 0 if clade.is_terminal() else sum([(1 + size[child.name]) for child in clade.clades])

    # return size




# def colour_leaves(v,S,colour):
    # '''
    # colour_leaves

    # Set all leaves below a specified clade to specified colour
    # '''
    # for clade in v.find_elements(order='postorder'):
        # if clade.is_terminal():
            # S[clade.name].colour = colour


# def Count(species,v,HT2,sizes=[]):
    # '''Code from Figure 8'''
    # def get_subtree(v,small=True):
        # a = v.clades[0]
        # b = v.clades[1]
        # m = sizes[a.name]
        # n = sizes[b.name]
        # return m if small == (m<n) else n

    # if v.is_terminal():
        # v.colour = Colour.C
        # return 0
    # else:
        # colour_leaves(get_subtree(v), Colour.B)
        # x = HT2.get_node_count(v)
        # colour_leaves(get_subtree(v), Colour.C)
        # y = Count(species,get_subtree(v, small=False), HT2,sizes=sizes)
        # colour_leaves(get_subtree(v), Colour.A)
        # z = Count(species,get_subtree(v), HT2, sizes=sizes)
        # return x + y + z





# def draw_tree(T,title,
              # ax         = None,
              # label_func = str,
              # label_colors = lambda c:None):

    # draw(T,
        # axes       = ax,
        # do_show    = False,
        # label_func = label_func,
        # label_colors= label_colors)
    # ax.set_title(title)
    # ax.set_xlabel('')
    # ax.set_ylabel('')
    # ax.get_xaxis().set_ticks([])
    # ax.get_yaxis().set_ticks([])





# def qetd0(species,T1,T2):
    # S                   = {s:Species(s,colour=Colour.A) for s in species}
    # Tr1                 = read(StringIO(T1), 'newick')
    # Tr2                 = read(StringIO(T2), 'newick')
    # root_name           = root_with_specified_leaf(Tr1,S)
    # S[root_name].colour = Colour.C
    # Factory             = HTreeBuilder()
    # HT2                 = Factory.build(Tr2)
    # links1,links2       = link(Tr1,HT2,S)
    # HT2.decorate(S)
    # return Count(species,
                 # Tr1.root,
                 # HT2,
                 # sizes = cache_sizes(Tr1))

# def label_colors(name,S=[]):
    # '''
    # label_colors

    # Used to ensure that colour of dislayed label matches assigned colour
    # '''
    # if name in S:
        # return Colour.as_text(S[name].colour)
    # else:
        # return None


if __name__=='__main__':
    start = time()
    parser = ArgumentParser(__doc__)
    # parser.add_argument('--paper',    default=False, action='store_true', help='process dataset from paper')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    # parser.add_argument('--show',     default=False, action='store_true', help='display plots')
    # parser.add_argument('--tests',    default=False, action='store_true', help='exe')
    args = parser.parse_args()

    # if args.paper:
        # F3((1,2,3,4,5,6,7,8,9))
        # species = ['a', 'b', 'c', 'd','e', 'f']
        # T0      = read(StringIO('a,((e,f),d),(b,c);'), 'newick')
        # T1      = read(StringIO('a,((e,f),d),(b,c);'), 'newick')
        # T2      = read(StringIO('a,((b,d),c),(e,f);'), 'newick')
        # root_with_specified_leaf(T1, species, index = 3)

        # Factory = HTreeBuilder()
        # HT2     = Factory.build(T2)
        # S       = {s:Species(s,colour=Colour.A) for s in species}
        # if args.show:
            # fig = figure(figsize=(12,12))

            # ax11 = fig.add_subplot(221)
            # draw_tree(T0,'T1',ax=ax11)

            # ax12 = fig.add_subplot(222)
            # draw_tree(T2,'T2',
                      # ax           = ax12,
                      # label_colors = lambda name: label_colors (name,S=S))

            # colour_leaves(T1.root.clades[0].clades[1],S,Colour.C)
            # ax21 = fig.add_subplot(223)
            # draw_tree(T1,'T1 re rooted',
                      # ax           = ax21,
                      # label_colors = lambda name: label_colors (name,S=S))

            # ax22 = fig.add_subplot(224)
            # draw_tree(HT2,'HT2',
                      # ax           = ax22,
                      # label_func   = ComponentClade.get_label,
                      # label_colors = lambda name: label_colors (name,S=S))


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

    # if args.show:
        # show()
