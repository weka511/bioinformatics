#    Copyright (C) 2019-2020 Greenweaves Software Limited
#
#    This is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This software is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>
#
#    LING   Linguistic Complexity of a Genome

import argparse
import os
import time
from   helpers import read_strings
from   ukkonen import build,Node

#  ling
#
#  Given: A DNA string s of length at most 100 kbp.
#
#  Return: The linguistic complexity lc(s)
#
# This follows  Troyanskaya OG, Arbell O, Koren Y, Landau GM, Bolshoy A.
# Sequence complexity profiles of prokaryotic genomic sequences: a fast algorithm 
# for calculating linguistic complexity. Bioinformatics. 2002 May;18(5):679-88.
# doi: 10.1093/bioinformatics/18.5.679. PMID: 12050064.
#
# "Therefore the number of non-terminal characters represented on the edges of the suffix
#   tree T for string S is equal to the number of different subwords contained in S."

def ling(string,a=4):
    def possible(k):
        return a if k==1 else min(a**k,n-k+1)
    tree, pst = build(string, regularize=True)
       
    n   = len(string)
    m   = [possible(k) for k in range(1,n+1)]
    
    return Node.count(tree, pst, ed='#')  /sum(m)
        
def convert_to_indices(string,lookup = {'A':0,'C':1,'G':2,'T':3}):
    return [lookup[c] for c in string] 


        
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
        
    if args.sample:
        print (ling('ATTTGGATT'))

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')

        Result = ling(Input[0])
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result:.3f}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s') 
