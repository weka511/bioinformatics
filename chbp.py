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

#  chbp Character-Based Phylogeny

import argparse
import os
import time
from   helpers import read_strings
import numpy as np

class Clade:
    @staticmethod
    def toString(Clades,species):
        def bfs(clade,suffix='',index=len(species)+1):
            if len(clade.children)>0:
                return f'({",".join([bfs(Clades[child],index=child) for child in clade.children])}){suffix}'
            if index<len(species):
                return species[index] 
        return bfs(Clades[-1],suffix=';')
        
    def __init__(self,index=None,character=None,children=[]):
        self.index     = index
        self.character = character
        self.fresh     = True
        self.children  = children
        self.weight    = max(1,len(children))
        
    def distance(self,other):
        return sum([abs(a-b)/(self.weight+other.weight) for a,b in zip(self.character,other.character)] )
        #a*self.weight-b*other.weight
        #self.weight+other.weight
        
    
def chbp(species,character_table):
    characters = sorted([([c[i] for c in character_table],i) for i in range(len(species))])
 
    x=0      
    Clades = [Clade(index=i, 
                    character=[c[i] for c in character_table]) for i in range(len(species))]
    while True:
        D     = [(i,j,Clades[i].distance(Clades[j]))     \
                        for i in range(len(Clades))      \
                        for j in range(i+1,len(Clades))  \
                        if Clades[i].fresh and Clades[j].fresh][::-1]
        index = np.argmin([d for _,_,d in D])
        i,j,_ = D[index]
        Clades.append(Clade(index=len(Clades),
                            children=[i,j],
                            character = np.mean([Clades[child].character for child in [i,j]],axis=0)))
        #if i<len(species):
            #print (species[i])
        #if j<len(species):
            #print (species[j])            
        Clades[i].fresh = False
        Clades[j].fresh = False
        if len(D)<2: break
        #print ()
    return Clade.toString(Clades,species)
    
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
        print (chbp(
                ['cat', 'dog', 'elephant', 'mouse', 'rabbit', 'rat'],
                [expand_as_ints('011101'),
                expand_as_ints('001101'),
                expand_as_ints('001100')] ))
        
    if args.abcde:
        clade = chbp(
            ['a', 'b', 'c', 'd', 'e'],
            [expand_as_ints('00000'),
             expand_as_ints('00110'),
             expand_as_ints('01000'),
             expand_as_ints('01110'),
             expand_as_ints('01111')] )
        print (f'{clade}')        
  
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        species = Input[0].split()
        character_table = [expand_as_ints(row) for row in Input[1:]]
        clade = chbp(species,character_table)
        print (clade)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{clade}\n')
                
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
