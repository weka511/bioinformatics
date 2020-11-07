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

#  chbp Character-Based Phylogeny

import argparse
import os
import time
from   helpers import read_strings
import numpy as np

class Clade:
    seq = 0
    def __init__(self,symbols,character,weight=1):
        self.symbols   = symbols
        self.character = character
        self.seq       = Clade.seq
        self.weight    = weight
        Clade.seq      += 1
    def distance(self,other):
        return sum([abs(a/self.weight-b/other.weight) for a,b in zip(self.character,other.character)] )
    def __str__(self):
        return f'{self.seq}: {self.symbols} {self.character} {self.weight}]'
    
def chbp(species,character_table):
    clades = {}
    for clade in [Clade(species[i],
                        [character_table[j][i] for j in range(len(character_table))]) for i in range(len(species))]:
        clades[clade.seq] = clade
    while len(clades)>1:
        keys      = [(cl1,cl2) for cl1 in clades.keys() for cl2 in clades.keys() if cl1<cl2]
        distances = [ clades[i].distance(clades[j]) for i,j in keys]       
        index     = np.argmin(distances)
        i,j       = keys[index]
        clade1    = clades[i]
        clade2    = clades[j]
        merged    = Clade(f'({clade1.symbols},{clade2.symbols})',
                          [a+b for a,b in zip(clade1.character,clade2.character)],
                          weight=clade1.weight+clade2.weight)
        print (f'Merged {clade1.symbols} and {clade2.symbols}=>{merged}')
        clades[merged.seq] = merged
        del clades[clade1.seq]
        del clades[clade2.seq] 
    #v = list(clades.values())
    return list(clades.values())[0]

def expand_as_ints(s):
    return [int(c) for c in s]

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--abcde',   default=False, action='store_true', help="process Jonathan Dursi's dataset")
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        clade = chbp(
            ['cat', 'dog', 'elephant', 'mouse', 'rabbit', 'rat'],
            [expand_as_ints('011101'),
             expand_as_ints('001101'),
             expand_as_ints('001100')] )
        print (f'{clade.symbols};')
        
    if args.abcde:
        clade = chbp(
            ['a', 'b', 'c', 'd', 'e'],
            [expand_as_ints('00000'),
             expand_as_ints('00110'),
             expand_as_ints('01000'),
             expand_as_ints('01110'),
             expand_as_ints('01111')] )
        print (f'{clade.symbols};')        
  
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        species = Input[0].split()
        character_table = [expand_as_ints(row) for row in Input[1:]]
        clade = chbp(species,character_table)
        print (f'{clade.symbols};')
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{clade.symbols};\n')
                
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
