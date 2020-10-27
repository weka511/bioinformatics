#  Copyright (C) 2020 Greenweaves Software Limited
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# BA9D Find the Longest Repeat in a String

import argparse
from helpers import read_strings
from snp import SuffixTree   
import time

def FindLongestRepeat(string):
    def sort_by_length(edges,reverse=False):
        return [e for _,e in sorted([(len(edge),edge) for edge in edges],reverse=reverse)]
        
    tree = SuffixTree()
    tree.build(string)
    edges = sort_by_length(tree.collectEdges(),reverse=True)
    for edge in edges:
        index = string.find(edge)
        if index==-1: break
        if string.find(edge,index+1)>-1:
            return edge


if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9D Find the Longest Repeat in a String')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (FindLongestRepeat('ATATCGTTTTATCGTT'))
        
    if args.extra:
        pass
    if args.rosalind:
        pass
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')