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
from   snp      import SuffixArray
from   rosalind import RosalindException

# SuffixArray2Tree
#
# Construct a suffix tree from the suffix array and LCP array of a string.
#
# Given: A string Text, SuffixArray(Text), and LCP(Text).
#
# Return: The strings labeling the edges of SuffixTree(Text). (You may return these strings in any order.)

def SuffixArray2Tree(Text, SuffixArray, LCP, Debug=False):
    
    # Node
    class Node:
        
        # __init__
        #
        # Paramters:
        #     root      Used to mark node as root
        
        def __init__(self,root=False):
            self.Edges  = []
            self.entry  = self if root else None
         
        # link
        
        def link(self,edge):
            self.Edges.append(edge)
            edge.link(self)
         
        # is_root
        #
        # Verify that node is a root
        def is_root(self):
            return self.entry == self
        
        def get_descent(self):
            Path    = []
            current = self
            while not current.is_root():
                edge = current.entry
                Path.insert(0,(current,len(edge.label)))
                current = edge.entry
            Path.insert(0,(root,0))
            Descents = []
            descent  = 0
            for node,contribution in Path:
                descent += contribution
                Descents.append((node,descent))
            return Descents
        
        def pop_rightmost_edge(self):
            return self.Edges.pop(-1)
        
        def gather_edges(self,Collection=[]):
            for edge in self.Edges:
                yield edge.label
                yield from edge.destination.gather_edges()
    # Edge
    
    class Edge():
        def __init__(self,label,destination):
            self.label        = label
            self.destination  = destination
            self.entry        = None
            destination.entry = self
            
        def link(self,entry):
            self.entry = entry
    
    def find_descent_le_LCP(last,lcp):
        for v,descent in last.get_descent():
            if descent <= lcp:
                return v,descent
                    
    n    = len(Text)
    root = Node(root=True)
    last = root
    
    for i in range(n):
        if i==-1:
            last = Node()
            root.link(Edge(Text[SuffixArray[i]:],last))
        else:
            v,descent = find_descent_le_LCP(last,LCP[i])   

            if descent == LCP[i]:
                x    = Node()
                v.link(Edge(Text[SuffixArray[i] + LCP[i]:], x))
                last = x         
            elif descent < LCP[i]:
                vw_edge_for_deletion = v.pop_rightmost_edge()
                y                    = Node()
                new_label            = Text[SuffixArray[i]+descent:]
                vy                   = Edge(new_label[:LCP[i]-descent],y)
                v.link(vy)
                yw                   = Edge(Text[SuffixArray[i-1]+LCP[i]:],vw_edge_for_deletion.destination)
                y.link(yw)
                x                    = Node()
                yx                   = Edge( Text[SuffixArray[i]+LCP[i]:],x)
                y.link(yx)
                last                 = x
                if Debug:
                    print (f'Deleted {vw_edge_for_deletion.label}. Split: {new_label}')
                    print (f'vy:{vy.label}, yw:{yw.label}, yx:{yx.label}')
            else:
                raise RosalindException(f'Descent ({descent}) > LCP[{i}]({LCP[i]})')
                
        if Debug:
            print ('----------------')
            print (f'{i}')
            for edge in root.gather_edges():
                print (edge)
    
    yield from root.gather_edges()

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9R Construct a Suffix Tree from a Suffix Array ')
    parser.add_argument('--banana',   default=False, action='store_true', help='process bana example from textbook')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    parser.add_argument('--debug',    default=False, action='store_true', help='Print debugging information')
    args = parser.parse_args()
    
    if args.banana:
        r,p,LCP = SuffixArray('panamabananas$',auxiliary=True)
        for edge in SuffixArray2Tree('panamabananas$', r, [0]+LCP,Debug=args.debug):
            print (edge)
            
    if args.sample:
        for edge in SuffixArray2Tree('GTAGT$',
                                     [5, 2, 3, 0, 4, 1],
                                     [0, 0, 0, 2, 0, 1],
                                     Debug=args.debug):
            print (edge)
        
    
    if args.extra:
        Input,Expected  = read_strings('data/SuffixTreeFromSuffixArray.txt',init=0)
        Expected.sort()
        Result = [edge for edge in SuffixArray2Tree(Input[0],
                                     [int(s) for s in Input[1].split(',')],
                                     [int(s) for s in Input[2].split(',')],
                                     Debug=args.debug)]
        Result.sort()
        print (f'Expected {len(Expected)} Edges, actual = {len(Result)}')
       
        
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for edge in SuffixArray2Tree(Input[0],
                                         [int(s) for s in Input[1].split(',')],
                                         [int(s) for s in Input[2].split(',')],
                                         Debug=args.debug):
                print (edge)
                f.write(f'{edge}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
