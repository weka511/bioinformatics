#  Copyright (C) 2019 Greenweaves Software Limited
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
# BA9C Construct the Suffix Tree of a String 

import argparse
from helpers import read_strings
from snp import SuffixTree   

           

    
# Each node has a symbol, position and labels of further edges


def check(Edges,Expected):
    print (f'Expected = {len(Expected)} edges, actual={len(Edges)}')
    mismatches = 0
    for a,b in zip(sorted(Edges),sorted(Expected)):
        if a!=b:
            mismatches+=1
            print (f'Expected {b}, was {a}')
    print(f'{0} mismatches')

def compare_edges(Edges,Expected):
    print (f'Expected = {len(Expected)} edges, actual={len(Edges)}')
    expected = iter(sorted(Expected))
    edges    = iter(sorted(Edges))
    exp      = next(expected)
    ed       = next(edges)
    while exp != '-' and ed !='-':
        if exp<ed:
            print('{0},{1}'.format(exp,'-'))
            exp = next(expected,'-')           
        elif ed<exp:
            print('{0},{1}'.format('-',ed))
            ed = next(edges,'-') 
        else:
            exp = next(expected,'-')
            ed = next(edges,'-')
            
if __name__=='__main__':
    parser = argparse.ArgumentParser('BA9C Construct the Suffix Tree of a String ')
    parser.add_argument('--sample', default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',  default=False, action='store_true', help='process extra dataset')
    args = parser.parse_args()
    if args.sample:
        tree = SuffixTree()
        tree.build('ATAAATG$')
        for edge in tree.collectEdges():
            print (edge)
         
    if args.extra:
        Input,Expected  = read_strings('data/SuffixTreeConstruction.txt',init=0)
        tree = SuffixTree()
        tree.build(Input[0])
        #tree.print()
        Edges  = tree.collectEdges()   
        
        compare_edges(Edges,Expected) 
            
            