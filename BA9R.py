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
from   helpers import read_strings

def SuffixArray2Tree(Text, SuffixArray, LCP):
    class Node:
        def __init__(self):
            self.Edges  = []
            self.entry = None
        def link(self,edge):
            self.Edges.append(edge)
            edge.destination.entry = edge
            edge.entry             = self
        def is_root(self):
            return self.entry == None
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
            
    class Edge():
        def __init__(self,label,destination):
            self.label       = label
            self.destination = destination
            self.entry       = None
            
    n    = len(Text)
    root = Node()
    for i in range(n):
        if i==0:
            last = Node()
            root.link(Edge(Text[SuffixArray[i]:],last))
        else:
            for v,descent in last.get_descent():
                if descent <= LCP[i]:
                    break
            if descent == LCP[i]:
                x     = Node()
                v.link(Edge(Text[SuffixArray[i] + LCP[i]:],
                            x))
                last = x
            if descent < LCP[i]:
                vw = v.pop_rightmost_edge()
                label = vw.label
                y     = Node()
                new_label = []
                for j in range(n):
                    if SuffixArray[i]+j<n and SuffixArray[i-1]+j <n and Text[SuffixArray[i]+j]== Text[SuffixArray[i-1]+j]:
                        new_label.append(Text[SuffixArray[i]+j])
                x=0
                break            
    
    yield 'foo'

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9R Construct a Suffix Tree from a Suffix Array ')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        for edge in SuffixArray2Tree('GTAGT$',
                                     [5, 2, 3, 0, 4, 1],
                                     [0, 0, 0, 2, 0, 1]):
            print (edge)
        
    

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for edge in SuffixArray2Tree(Input[0],Input[1],Input[2]):
                print (edge)
                f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
