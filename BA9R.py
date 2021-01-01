#   Copyright (C) 2020 Greenweaves Software Limited

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

#  BA9R 	Construct a Suffix Tree from a Suffix Array 

import argparse
import os
import time
from   helpers  import read_strings
from   rosalind import RosalindException
from   snp      import SuffixArray

# SuffixArray2Tree
#
# Construct a suffix tree from the suffix array and LCP array of a string.
#
# Given: A string Text, SuffixArray(Text), and LCP(Text).
#
# Return: The strings labeling the edges of SuffixTree(Text). (You may return these strings in any order.)

def SuffixArray2Tree(Text, SuffixArray, LCP, trace=False):
    
    # Node
    #
    # Represents a node in Suffix Tree
    #
    # Data Members:
    #    Edges
    #    entry       Will be set to the Edge that points into the Node
    #                once linked into tree, except Root. which points to self 
    #    count       Unique id for each Node FIXME do we really need this?
    #    descent
    #    label       Child nodes are labelled by position in suffix array   
    #
    class Node:
        count = 0   # Used to assign a unique id to each Node FIXME do we really need this?
        # __init__
        #
        # Parameters:
        #     root      Used to mark node as root
        
        def __init__(self,root=False,label=None):
            self.Edges   = []
            self.entry   = self if root else None
            self.count   = Node.count
            self.descent = 0 if root else None
            self.label   = label
            Node.count += 1
            
        # Attach an Edge to this Node
        
        def link(self,edge):
            self.Edges.append(edge)
            edge.link(self)
         
        # is_root
        #
        # Verify that node is a root
        def is_root(self):
            return self.entry == self
        
        # unlink
        #
        # Unlink rightmost edge from this node
        #
        # Returns: newly deleted node
        
        def unlink(self):
            return self.Edges.pop(-1)
        
        # gather_edge_labels
        #
        # Used to collect all edges of tree
        # This is just a simple depth-first search
        def gather_edge_labels(self):
            for edge in self.Edges:
                yield edge.label,edge.destination.label
                yield from edge.destination.gather_edge_labels()
                
        # find_attachment_point
        #
        # Starting with last node added, search up until descent <= limit. This will always 
        # return a value provided limit is from a legal LCP array, as DESCENT(root)==0
        
        def find_attachment_point(self,limit):
            v = self
            while v.descent>limit:
                v = v.entry.entry
            return v,v.descent
        
        def get_concatenated_edge_labels(self):
            Concat = []
            v = self
            while True:
                if v.is_root():
                    return ''.join(Concat)
                edge = v.entry
                Concat.insert(0,edge.label)
                v    = edge.entry
    # Edge
    #
    # Represents an Edge in Suffix Tree
    #
    # Data Members:
    #      label         Text used to label edge
    #      destination   Node that Edge points to
    #      entry         Node that Edge comes from
    
    class Edge():
        
        # __init__
        #
        # Create Edge and set destination Node to point back to edge
        
        def __init__(self,label,destination):
            self.label        = label
            self.destination  = destination
            self.entry        = None
            destination.entry = self
        
        # link
        #
        # Used by Node link to attach this Edge, and update the descent of the Node
        # that Edge points to
        
        def link(self,entry):
            self.entry               = entry
            self.destination.descent = entry.descent + len(self.label)
    
    # get_common_prefix
    #
    # Determine common prefix of two suffixes
    # Currently used only in assertions
    
    def get_common_prefix(suffix0,suffix1):
        for i in range(min(len(suffix0),len(suffix1))):
            if suffix0[i]!=suffix1[i]:
                return suffix0[:i]
           
    n    = len(Text)
    root = Node(root=True)
    
    # x, y, v, w are nodes
    # x will always be the last node that has been added
    # two letter combinations of nodes, e.g. vx, represent Edges
    
    x    = root
    
    for i in range(n): # cf Biometrics Algorithms 9R; my 'i' represents textbook's 'i+1'
        if trace:
            print ('----------------')
            print (f'i={i}, SuffixArray[{i}]={SuffixArray[i]}--{Text[SuffixArray[i]:]}/{Text[SuffixArray[i-1]:]} LCP[{i}]={LCP[i]}') 
            
        v,descent = x.find_attachment_point(LCP[i])
        if trace:
            print (f'Concatenated Edge labels: {v.get_concatenated_edge_labels()}')

        if descent == LCP[i]:                             # Concatenation of labels on path is equal to
                                                          # longest common prefix of suffixes corresponding to
                                                          # SuffixArray[i] and SuffixArray[i-1]
            x  = Node(label=SuffixArray[i])               # Insert suffix as a new leaf
            vx = Edge(Text[SuffixArray[i] + LCP[i]:], x)  # Edge label is suffix
            v.link(vx)                                    # Link back to V
            if trace:
                print (f'Added vx: {vx.label}')
                
        elif descent < LCP[i]:     # Concatenation of labels on path has fewer symbols than we need
                                   # to match the longest common prefix of suffixes corresponding to
                                   # SuffixArray[i] and SuffixArray[i-1]            
                                   # In this case we need to split a Node
            # Delete (v,w)
            vw        = v.unlink()
            w         = vw.destination
 
            # Add new node y and edge (v,y)
            y         = Node()
            new_label = Text[SuffixArray[i]+descent:]
            vy_label  = new_label[:LCP[i]-descent]
            vy        = Edge(vy_label, y)
            v.link(vy)
            if trace:
                print (f'descent={descent}. Deleting {vw.label}. Split: {new_label} vy: {vy.label}')  
                print (f'Common prefix of {Text[SuffixArray[i]:]} and {Text[SuffixArray[i-1]:]} is {get_common_prefix(Text[SuffixArray[i]:],Text[SuffixArray[i-1]:])}')
            #assert vy.label== get_common_prefix(Text[SuffixArray[i]:],Text[SuffixArray[i-1]:]), \
                              #f'i={i}: {vy.label}!= {get_common_prefix(Text[SuffixArray[i]:],Text[SuffixArray[i-1]:])} {Text[SuffixArray[i]:]} {Text[SuffixArray[i-1]:]}'
            
            # Verify that descent of y has been calculated correctly
            assert y.descent==LCP[i],f'{y.descent} != {LCP[i]}'
            
            # Connect w to y. This will set w to point back to edge
            
            yw        = Edge(Text[SuffixArray[i-1]+LCP[i]:], w)
            y.link(yw)
            #assert yw.label==vw.label[len(vy.label):],f'i={i}, {yw.label}!={vw.label[len(vy.label):]}'
            
            # Add x
            x         = Node(label=SuffixArray[i])
            yx        = Edge( Text[SuffixArray[i]+LCP[i]:], x)
            y.link(yx)
            
            if trace:
                print (f'Deleted {vw.label}. Split: {new_label}')
                print (f'vy: {vy.label}, yw: {yw.label}, yx: {yx.label}')
        else:
            raise RosalindException(f'Descent ({descent}) > LCP[{i}]({LCP[i]})')
           
        if trace:
            print ('Edges:')
            for edge_label,node_label in root.gather_edge_labels():
                print (f'\t{edge_label} {node_label if node_label != None else ""}')

    yield from root.gather_edge_labels()

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9R Construct a Suffix Tree from a Suffix Array ')
    parser.add_argument('--banana',   default=False, action='store_true', help='process bana example from textbook')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    parser.add_argument('--trace',    default=False, action='store_true', help='Print debugging information')
    args = parser.parse_args()
    
    if args.banana:
        r,p,LCP = SuffixArray('panamabananas$',auxiliary=True,padLCP=True)
        for edge,_ in SuffixArray2Tree('panamabananas$', r, LCP,trace=args.trace):
            print (edge)
            
    if args.sample:
        for edge,_ in SuffixArray2Tree('GTAGT$',
                                     [5, 2, 3, 0, 4, 1],
                                     [0, 0, 0, 2, 0, 1],
                                     trace=args.trace):
            print (edge)
        
    
    if args.extra:
        Input,Expected  = read_strings('data/SuffixTreeFromSuffixArray.txt',init=0)

        Result = [edge for edge,_ in SuffixArray2Tree(Input[0],
                                                    [int(s) for s in Input[1].split(',')],
                                                    [int(s) for s in Input[2].split(',')],
                                                    trace=args.trace)]
        print (f'Expected {len(Expected)} Edges, actual = {len(Result)}')
        Result.sort(key=lambda x: f'{len(x):04}{x}')
        Expected.sort(key=lambda x: f'{len(x):04}{x}')
        i = 0
        j = 0
        while i<len(Expected) and j<len(Result):
            if Expected[i]==Result[j]:
                i+=1
                j+=1
            elif Expected[i]<Result[j]:
                print (f'Expected {Expected[i]} < Actual {Result[j]}')
                i+=1
            else: # Expected[i]>Result[j]
                print (f'Expected {Expected[i]} > Actual {Result[j]}')
                j+=1
                
  
        while i<len(Expected):
            print (f'Expected {Expected[i]}, Actual ---')
            i+=1            
        
        while j<len(Result):
            print (f'Expected ---, Actual {Result[j]}')
            j+=1        
        
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for edge,_ in SuffixArray2Tree(Input[0],
                                         [int(s) for s in Input[1].split(',')],
                                         [int(s) for s in Input[2].split(',')],
                                         trace=args.trace):
                print (edge)
                f.write(f'{edge}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
