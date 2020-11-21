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


# >seq01
#GCGCCGGCGGCCGCGUACAGUUAACAUUACGAUUAUGAUAAUAGCUUAUCGCUCAUAUGA
#GCAUCGCCGGGCGCAGUACAUUGCAUGCGCGGUAUUAUAAGCCCAUGAUCUAUAGUGCCG
#ACCAUGCAUAUAUAUAUUAGCGCGGCGAUCAUCGAUCGCGUACUCUGCAGAGAUUAGAUG
#CACGCGCUCGAGUGCGCGCACGUGCGCCGACGGCUGCCGGAUAUAUCAUAUAUUAUAGCC
#CGCGCGCGCGCCAUGGGC
#
# expect 412480

def cat(s):
    pass

def perfect_matchings(s):
    def isCrossing(i,j,k,l):
        return i<k and k<j and j<l
    
    def doesMatch(a,b):
        return b==matches[a]
    
    matches = {
        'A':'U',
        'U':'A',
        'C':'G',
        'G':'C'
    }
    s = ''.join([c for c in s if c in 'AUCG']) # purge white spaces and other junk
    for i in range(len(s)):
        for k in range(i+1,len(s)):
            print ('------------')
            for j in range(i+1,len(s)):
                for l in range(k+1,len(s)):
                    if i in [j,k,l] or j in [k,l] or k==l: continue
                    if isCrossing(i,j,k,l): continue
                    if doesMatch(s[i],s[j]) and doesMatch(s[k],s[l]):
                        print (f'{{{i},{j}}}, {{{k},{l}}}')

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
        perfect_matchings('UAGCGUGAUCAC')
        #print (catalan(30))
        #print (cat('AUAU')%1000000)
    
  
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

