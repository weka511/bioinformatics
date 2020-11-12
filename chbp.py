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


class Clade():
    @staticmethod
    def toString(Root):
        def bfs(clade,suffix=''):
            if len(clade.children)>0:
                pp = [bfs(child) for child in clade.children]
                ppp=[pppp for pppp in pp if len(pppp)>0]
                ll = ",".join(ppp) 
                return f'({ll}){suffix}'
            else:
                return ','.join(clade.species)
            
        return bfs(Root,suffix=';')    
    def __init__(self,species):
        self.species = species
        self.children = []
    
    def bfs(self,prefix=''):
        print (f'{prefix}{self.species}')
        for child in self.children:
            child.bfs(prefix+'-')
  
    def newick(self,suffix=''):
        def wrap(value):
            return f'({value})' if len(value)>1 else value
        if len(self.children)==0:
            if len(self.species)==0:
                return ''
            elif len(self.species)==1:
                return self.species[0]
            else:
                return '(' + ','.join(self.species) + ')'
        values = [child.newick()  for child in self.children]
        non_trivial = [v for v in values if len(v)>0 and len(v[0])>0]
        if len(non_trivial)==0:
            return ''
        if len(non_trivial)==1:
            return non_trivial[0]
        return '(' + ','.join(non_trivial) + ')' + suffix
        #return ','.join([wrap(child.newick())  for child in self.children])
      
            
def chbp(species,character_table):
    def entropy(freq):
        p1 = freq/n
        p2 = 1-p1
        return - p1 *np.log(p1) - p2 * np.log(p2)
    
    n          = len(species)
    char_table_info = sorted(character_table,reverse=True,key=lambda char:entropy(sum(char)) )
    species_table = {species[i]:i for i in range(n)}
    Root = Clade(species)
    Current = [Root]
    for char in char_table_info:
        print (char)
        for clade in Current:
            Left = []
            Right = []
            for s in clade.species:
                i = species_table[s]
                if char[i]==0:
                    Left.append(s)
                else:
                    Right.append(s)
            clade.children = [Clade(Left),Clade(Right)]
        Current = [child for clade in Current for child in clade.children] 
    #Root.bfs()
    #s =  Root.newick(suffix=';')
    return Root.newick(suffix=';')#Clade.toString(Root)

 
    
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



    #class Clade:
        #@staticmethod
        #def toString(Clades,species):
            #def bfs(clade,suffix='',index=len(species)+1):
                #if len(clade.children)>0:
                    #return f'({",".join([bfs(Clades[child],index=child) for child in clade.children])}){suffix}'
                #if index<len(species):
                    #return species[index] 
            #return bfs(Clades[-1],suffix=';')
            
        #def __init__(self,index=None,character=None,children=[]):
            #self.index     = index
            #self.character = character
            #self.fresh     = True
            #self.children  = children
            #self.weight    = max(1,len(children))
            
        #def distance(self,other):
            #return sum([abs(a-b)/(self.weight+other.weight) for a,b in zip(self.character,other.character)] )
            ##a*self.weight-b*other.weight
            ##self.weight+other.weight
            
            
            #characters = [[c[i] for c in character_table] for i in range(len(species))]
            #freqs      = [entropy(sum(character))  for character in character_table]
            #x=0      
            #Clades = [Clade(index=i, 
            #character=[c[i] for c in character_table]) for i in range(len(species))]
            #while True:
            #D     = [(i,j,Clades[i].distance(Clades[j]))     \
            #for i in range(len(Clades))      \
            #for j in range(i+1,len(Clades))  \
            #if Clades[i].fresh and Clades[j].fresh][::-1]
            #index = np.argmin([d for _,_,d in D])
            #i,j,_ = D[index]
            #Clades.append(Clade(index=len(Clades),
            #children=[i,j],
            #character = np.mean([Clades[child].character for child in [i,j]],axis=0)))
            ##if i<len(species):
            ##print (species[i])
            ##if j<len(species):
            ##print (species[j])            
            #Clades[i].fresh = False
            #Clades[j].fresh = False
            #if len(D)<2: break
            ##print ()
            #return Clade.toString(Clades,species)            
