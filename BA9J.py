# Copyright (C) 2020 Greenweaves Software Limited

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

# BA9J Reconstruct a String from its Burrows-Wheeler Transform

import argparse
import os
import time
from helpers import read_strings
from snp import BurrowsWheeler

def InverseBWT(string):
    def getN(ch,row,column):
        n = 0
        i = 0
        while i<=row:
            if ch==column[i]:
                n+=1
            i+=1
        return n
    
    def get_char(ch,column,seq):
        pos   = 0
        count = 0
        for i in range(len(column)):
            if column[i]==ch:
                pos = i
                count+=1
                if count==seq:
                    return pos,count
            
    lastColumn  = [a for a in string]
    firstColumn = sorted(lastColumn)
    Result      = [firstColumn[0]]
    ch          = min(string)
    seq         = 1
    while len(Result)<len(string):
        row,count    = get_char(ch,lastColumn,seq)
        ch           = firstColumn[row]
        seq          = getN(ch,row,firstColumn)
        Result.append(ch)
        x=0
    return ''.join(Result[1:-1]+Result[0:1])

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9J Reconstruct a String from its Burrows-Wheeler Transform')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        #s = BurrowsWheeler('panamabananas$')
        #print (InverseBWT(s))
        print (InverseBWT('TTCCTAACG$A'))
    
    if args.extra:
        Input,Expected  = read_strings('data/InverseBWT.txt',init=0)
        trie = Trie(Input)
        Actual = None
        Expected.sort()
        print (len(Expected),len(Actual))
        diffs = [(e,a) for e,a in zip(Expected,Actual) if e!=a]
        print (diffs)
  
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
