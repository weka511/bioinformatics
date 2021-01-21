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
from   align import edit
from   numpy import argmin
from   Bio import SeqIO

# ctea
#
# Given: Two protein strings s and t in FASTA format, each of length at most 1000 aa.
#
# Return: The total number of optimal alignments of s and t with respect to edit alignment score, modulo 134,217,727.

def ctea(s,t,indel_cost=1,replace_cost=lambda a,b: 1, mod=134217727):
 
    class Move:
        
        #count     = 0
        Processed = {}
        Leaves    = {}
        
        def __init__(self,m,n):
            self.m = m
            self.n = n
 
        def explore(self,path=[]):
            while self.m>0 and self.n>0:
                if (self.m,self.n) in Move.Processed:
                    Move.Processed[(self.m,self.n)] += 1
                    return
                Move.Processed[(self.m,self.n)] = 1
                moves        = [(self.m-1,self.n),(self.m,self.n-1),(self.m-1,self.n-1)]
                scores       = [matrix[self.m-1][self.n]+indel_cost,
                                matrix[self.m][self.n-1]+indel_cost,
                                matrix[self.m-1][self.n-1] + (0 if s[self.m-1]==t[self.n-1] else replace_cost(s[self.m-1],
                                                                                                              t[self.n-1]))]
                index        = argmin(scores)
                alternatives = [moves[i] for i in range(len(scores)) if scores[i]==scores[index]]
                if len(alternatives)>1:
                    for m,n in alternatives:
                        predecessor = Move(m,n)
                        predecessor.explore(path=path+[(m,n)])
                    return
                else:
                    self.m,self.n = alternatives[0]
                    continue
            Move.Leaves[(self.m,self.n)] = path
            #Move.count +=1
        
    _,matrix = edit(s,t,indel_cost,replace_cost)
    move = Move(len(matrix) - 1, len(matrix[0]) - 1)
    move.explore()
    count = 0
    for path in Move.Leaves.values():
        product = 1
        for z in path:
            if z in Move.Processed:
                product *= Move.Processed[z]
        count += product
    return count %mod
    
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
