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

#  QRTD  Quartet Distance

import argparse
import os
import time
from   helpers import read_strings
from   phylogeny import create_adj, parse

def qrtd(species,newick1,newick2):
    
    def q(n):
        return (n*(n-1)*(n-2)*(n-3))/(4*3*2*1)  
    
    def replace_leaves(adj):
        return {parent:sorted([seiceps[child] if child in seiceps else child for child in children]) for parent,children in adj.items() }
 
    def dq(T1,T2):
        return q(count_leaves(T1)) +  q(count_leaves(T2)) -2 * shared_quartets(T1,T2)
    
    def count_leaves(adj):
        return sum([1 for children in adj.values() for child in children if child<n])
    
    def shared_quartets(T1,T2):
        return 0
    
    n       = len(species)
    seiceps = {species[i]:i for i in range(n)}
    tree1   = replace_leaves(create_adj(parse(newick1,start=n)))
    tree2   = replace_leaves(create_adj(parse(newick2,start=n)))
    return dq(tree1,tree2)

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print(qrtd(
            'A B C D E'.split(),
            '(A,C,((B,D),E));',
            '(C,(B,D),(A,E));'            
        ))
        
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = qrtd(Input[0].split(), Input[1], Input[2])
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
