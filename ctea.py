#   Copyright (C) 2020-2021 Greenweaves Software Limited

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

#  CTEA Counting Optimal Alignments

import argparse
import os
import time
from   helpers import read_strings
from   align import edit,topological_order
from   numpy import argmin
from   Bio import SeqIO

# ctea
#
# Given: Two protein strings s and t in FASTA format, each of length at most 1000 aa.
#
# Return: The total number of optimal alignments of s and t with respect to edit alignment score, modulo 134,217,727.

def ctea(s,t,indel_cost=1,replace_cost=lambda a,b: 1, mod=134217727):
    count        = 0
    Closed       = set()
    Open         = []
    Leaves       = []
    Adj          = {}
    def explore():
        mn = Open.pop(0)
        if mn in Closed: return
        Closed.add(mn)
        m,n = mn
        Adj[(m,n)] = []
        if m>0 and n>0:
            moves  = [(m-1,n),(m,n-1),(m-1,n-1)]
            scores = [matrix[m-1][n]+indel_cost,
                      matrix[m][n-1]+indel_cost,
                      matrix[m-1][n-1] + (0 if s[m-1]==t[n-1] else replace_cost(s[m-1],t[n-1]))]
            index  = argmin(scores)
            lowest = scores[index]
            candidates = [i for i in range(len(scores)) if scores[i]==lowest]
            Adj[(m,n)] = [moves[i] for i in candidates]
            
            for i in candidates:
                if i not in Closed:
                    Open.append((moves[i]))
        else:
            Leaves.append((m,n))
    
    def invert(Adj):
        Inverse = {a:[] for links in Adj.values() for a in links}
        for source,destinations in Adj.items():
            for destination in destinations:
                Inverse[destination].append(source)
        return Inverse
    
    def get_count_from(node,Inverse):
        N = 1
        parents = Inverse[node]
        return N
    
    _,matrix = edit(s,t,indel_cost,replace_cost)
    m,n = len(matrix)-1, len(matrix[0])-1
    Open.append((m,n))
    while len(Open)>0:
        explore() 
        
    ground          = (-1,-1)
    for leaf in Leaves:
        Adj[leaf] = [ground]    
    Inverse         = invert(Adj)

    C               = {ground:1}
    for node in topological_order(dict(Inverse)):
        if node == ground: continue
        C[node] = sum([C[x] for x in Adj[node]])       
 
    return C[(m,n)] %mod
    
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('CTEA Counting Optimal Alignments')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (ctea('PLEASANTLY','MEANLY'))
        
    if args.rosalind:
        inFile = open(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt','r')
        strings = []
        for record in SeqIO.parse(inFile,'fasta'):
            strings.append(str(record.seq))        
 
        Result = ctea(strings[0],strings[1])
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
