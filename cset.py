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

#  CSET  	Fixing an Inconsistent Character Set 

import argparse
import os
import time
from helpers import read_strings,expand

def cset(matrix):
    def find_inconsistencies(i,j,skip=-1):
        values = [[],[],[],[]]
        for k in range(n):
            if k==skip: continue
            values[2*matrix[k][i] + matrix[k][j]].append(k)
        if all([len(values[l])>0 for l in range(4)]):
            return values
        else:
            return []
        
    m = len(matrix[0])
    n = len(matrix)
    inconsistencies = [find_inconsistencies(i,j) for i in range(m) for j in range(i)]
    candidates = set()
    for values in inconsistencies:
        if len(values)>0:
            for value in values:
                if len(value)==1:
                    candidates.add(value[0])
    
    for skip in candidates:
        inconsistencies = [find_inconsistencies(i,j,skip=skip) for i in range(m) for j in range(i)]
        for values in inconsistencies:
            if len(values)==0:
                return [matrix[row] for row in range(n) if row!=skip] 
            for value in values:
                if len(value)==0:
                    return [matrix[row] for row in range(n) if row!=skip] 
        #x=0
    #return [matrix[row] for row in range(n) if row!=candidate]        

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (cset([expand('100001'),
                     expand('000110'),
                     expand('111000'),
                     expand('100111')]))
  
    if args.rosalind:
        Input           = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        character_table = [expand(row) for row in Input]
        Result          = cset(character_table)
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for row in Result:
                f.write(f'{"".join(str(i) for i in row)}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
