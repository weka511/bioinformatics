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

#  ALPH  Alignment-Based Phylogeny 

import argparse
import os
import time
from helpers import read_strings
from newick  import newick_to_adjacency_list

def alph(T,Alignment):
    def create_fixed_alignments():
        Leaves = {}
        for i in range(0,len(Alignment),2):
            Leaves[Alignment[i][1:]] = Alignment[i+1]
        return Leaves

    Adj   = newick_to_adjacency_list(T)   
    Fixed = create_fixed_alignments()
    assert len([node for node,value in Adj.items() if len(value)==0 and node not in Fixed])==0
    Assignment = {a:[] for a in Adj.keys()}
    
    d = 0
    return d,[(a,b) for a,b in Assignment.items() if len(Adj[a])==0]
    

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('ALPH  Alignment-Based Phylogeny')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (alph('(((ostrich,cat)rat,(duck,fly)mouse)dog,(elephant,pikachu)hamster)robot;',
                    ['>ostrich',
                     'AC',
                     '>cat',
                     'CA',
                     '>duck',
                     'T-',
                     '>fly',
                     'GC',
                     '>elephant',
                     '-T',
                     '>pikachu',
                     'AA'
                     ]))
        
  
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        d,Assignment = alph(Input[0],Input[1])
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            print (d)
            f.write(f'{line}\n')
            for label,String in Assignment:
                print (label)
                f.write(f'{label}\n')
                print (String)
                f.write(f'{Sringt}\n')                
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
