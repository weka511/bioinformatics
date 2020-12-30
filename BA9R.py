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

def SuffixArray2Tree(Text, SuffixArray, LCP, Debug=False):
    

    # Node
    class Node:
        count = 0
        # __init__
        #
        # Paramters:
        #     root      Used to mark node as root
        
        def __init__(self,root=False):
            self.Edges   = []
            self.entry   = self if root else None
            self.count   = Node.count
            self.descent = 0 if root else None
            Node.count += 1
        # link
        
        def link(self,edge):
            self.Edges.append(edge)
            edge.link(self)
         
        # is_root
        #
        # Verify that node is a root
        def is_root(self):
            return self.entry == self
        
   
        
        # pop_rightmost_edge
        #
        # Delete rightmost edge from this node
        #
        # Returns: newly deleted node
        
        def pop_rightmost_edge(self):
            return self.Edges.pop(-1)
        
        # gather_edges
        #
        # Used to collect all edges
        def gather_edges(self,Collection=[]):
            for edge in self.Edges:
                yield edge.label
                yield from edge.destination.gather_edges()
                
        # find_descent_LE_LCP
        #
        # Starting with last node added, search up until descent <= limit. This will always 
        # return a value provided limit is from a legal LCP arrau, as DESCENT(root)==0
        
        def find_descent_LE_LCP(self,limit):
            v = self
            while v.descent>limit:
                v = v.entry.entry
            return v,v.descent    
        
    # Edge
    
    class Edge():
        def __init__(self,label,destination):
            self.label        = label
            self.destination  = destination
            self.entry        = None
            destination.entry = self
            
        def link(self,entry):
            self.entry               = entry
            self.destination.descent = entry.descent + len(self.label)
    

    
    def get_common_prefix(sa0,sa1):
        for i in range(min(len(sa0),len(sa1))):
            if sa0[i]!=sa1[i]:
                return sa0[:i]
            
    n    = len(Text)
    root = Node(root=True)
    x    = root   # x is always the last node added
    
    for i in range(n): # cf Biometrics Algorithms 9R my 'i' represents their 'i+1'
        if Debug:
            print ('----------------')
            print (f'i={i}') 
            
        v,descent = x.find_descent_LE_LCP(LCP[i])   

        if descent == LCP[i]:
            x  = Node()
            vx = Edge(Text[SuffixArray[i] + LCP[i]:], x)
            v.link(vx)
            if Debug:
                print (f'Added vx: {vx.label}')
                
        elif descent < LCP[i]: 
            # Delete (v,w)
            vw_for_deletion = v.pop_rightmost_edge()
            w               = vw_for_deletion.destination
            ll              = vw_for_deletion.label
            
            # Add new node y and edge (v,y)
            y               = Node()
            new_label       = Text[SuffixArray[i]+descent:]
            vy_label        = new_label[:LCP[i]-descent]
            vy              = Edge(vy_label, y)
            v.link(vy)
            assert vy.label== get_common_prefix(Text[SuffixArray[i]:],Text[SuffixArray[i-1]:])
            
            # Verify that descent of y has been calculated correctly
            assert y.descent==LCP[i],f'{y.descent} != {LCP[i]}'
            
            # Connect w to y
            
            yw              = Edge(Text[SuffixArray[i-1]+LCP[i]:], w)
            y.link(yw)
            assert yw.label==vw_for_deletion.label[len(vy.label):],f'i={i}, {yw.label}!={vw_for_deletion.label[len(vy.label):]}'
            
            # Add x
            x               = Node()
            yx              = Edge( Text[SuffixArray[i]+LCP[i]:], x)
            y.link(yx)
            
            if Debug:
                print (f'Deleted {vw_for_deletion.label}. Split: {new_label}')
                print (f'vy: {vy.label}, yw: {yw.label}, yx: {yx.label}')
        else:
            raise RosalindException(f'Descent ({descent}) > LCP[{i}]({LCP[i]})')
           
        if Debug:
            print ('Edges:')
            for label in root.gather_edges():
                print (f'\t{label}')

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
        r,p,LCP = SuffixArray('panamabananas$',auxiliary=True,padLCP=True)
        for edge in SuffixArray2Tree('panamabananas$', r, LCP,Debug=args.debug):
            print (edge)
            
    if args.sample:
        for edge in SuffixArray2Tree('GTAGT$',
                                     [5, 2, 3, 0, 4, 1],
                                     [0, 0, 0, 2, 0, 1],
                                     Debug=args.debug):
            print (edge)
        
    
    if args.extra:
        Input,Expected  = read_strings('data/SuffixTreeFromSuffixArray.txt',init=0)

        Result = [edge for edge in SuffixArray2Tree(Input[0],
                                                    [int(s) for s in Input[1].split(',')],
                                                    [int(s) for s in Input[2].split(',')],
                                                    Debug=args.debug)]
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
