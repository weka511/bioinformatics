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
#    along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>
#
#    cat 	Catalan Numbers and RNA Secondary Structures (WIP)

import argparse
import os
import time
from helpers import read_strings

from rosalind import RosalindException,verify_counts_complete_graph
from Bio import SeqIO
import os.path

# catalan
#
# Calculate Catalan numbers
#
# Verified against list from http://mathforum.org/advanced/robertd/catalan.html
#
# 1 1
# 2 2
# 3 5
# 4 14
# 5 42
# 6 132
# 7 429
# 8 1430
# 9 4862
# 10 16796
# 11 58786
# 12 208012
# 13 742900
# 14 2674440
# 15 9694845
#
# Verified against Online Encyclopedia of Integer Sequences (https://oeis.org/A000108) for case:
# 30 3814986502092304

def catalan(n,trace=False):
    c=[1]
    for n0 in range(1,n+1):
        c.append(sum([c[k]*c[n0-1-k] for k in range(n0)]))
        if trace:
            print (n0,c[-1])
    return c[n]

def motzkin(n):
    s=[1,1]
    for n0 in range(1,n):
        s.append(s[-1]+sum([s[k]*s[n0-1-k] for k in range(n0)]))
        #print (s[-1])
    return s[n]

# >seq01
#GCGCCGGCGGCCGCGUACAGUUAACAUUACGAUUAUGAUAAUAGCUUAUCGCUCAUAUGA
#GCAUCGCCGGGCGCAGUACAUUGCAUGCGCGGUAUUAUAAGCCCAUGAUCUAUAGUGCCG
#ACCAUGCAUAUAUAUAUUAGCGCGGCGAUCAUCGAUCGCGUACUCUGCAGAGAUUAGAUG
#CACGCGCUCGAGUGCGCGCACGUGCGCCGACGGCUGCCGGAUAUAUCAUAUAUUAUAGCC
#CGCGCGCGCGCCAUGGGC
#
# expect 412480

def cat(s):
    return motzkin(int(len(s)/2))

def perfect_matchings(s):
    matches = {
        'A':'U',
        'U':'A',
        'C':'G',
        'G':'C'
    }
    s = ''.join([c for c in s if c in 'AUCG']) # purge white spaces and other junk
    for i in range(len(s)):
        for k in range(i+1,len(s)):
            for j in range(k+1,len(s)):
                for l in range(j+1,len(s)):
                    print (i,k,j,l)

def read_fasta(file,path=r'C:\Users\Simon\Downloads'):
    with open(os.path.join(path,file),'rU') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return record.seq
        
if __name__=='__main__':
    
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (catalan(30))
        print (cat('AUAU')%1000000)
    
  
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = None
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                f.write(f'{line}\n')
                
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    
    
   
    #print (cat('CGGCUGCUACGCGUAAGCCGGCUGCUACGCGUAAGC')) # should be 736
    
    #print(cat(read_fasta('data.fna')))

