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


def chbp(species,character_table):
     class Clade:
          def __init__(self,taxa):
               self.taxa = [s for s in taxa]
          
          #def __str__(self):
               #return '(' + ','.join(str(taxon) for taxon in self.taxa) +')'
          
          def newick(self):
               def conv(taxon):
                    if type(taxon)==int:
                         return species[taxon]
                    else:
                         return taxon.newick()
               return '(' + ','.join(conv(taxon) for taxon in self.taxa) +')'
          
          def split(self,character):
               left  = []
               right = []
               for i in self.taxa:
                    if character[i]==0:
                         left.append(i)
                    else:
                         right.append(i)
                         
               self.taxa = [Clade(left),Clade(right)]
         
          def splitAll(self,characters,depth=0):
               if depth<len(characters):
                    self.split(characters[depth])
                    for taxon in self.taxa:
                         taxon.splitAll(characters,depth+1)
                         
          #__repr__ = __str__     # https://stackoverflow.com/questions/12448175/confused-about-str-on-list-in-python
          
     def entropy(freq):
          if freq==0 or freq==n: return 0
          p1 = freq/n
          p2 = 1-p1
          return - p1 *np.log(p1) - p2 * np.log(p2)
     
     
     
     n               = len(species)
     entropies       = [entropy(sum(char)) for char in character_table]
     entropy_indices = np.argsort(entropies)
     characters      = [character_table[i] for i in entropy_indices[::-1]]
     indices         = list(range(len(species)))
     root            = Clade(indices)
     root.splitAll(characters)
     print (root.newick())
     #left,right      = split(species,characters[0])
     x=0


 
    
def expand(s):
     return [int(c) for c in s]

if __name__=='__main__':
     start = time.time()
     parser = argparse.ArgumentParser('Character-Based Phylogeny')
     parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
     parser.add_argument('--abcde',   default=False, action='store_true', help="process Jonathan Dursi's dataset")
     parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
     args = parser.parse_args()
     if args.sample:
          print (chbp(
               ['cat', 'dog', 'elephant', 'mouse', 'rabbit', 'rat'],
               [expand('011101'),
               expand('001101'),
               expand('001100')] ))
        
     if args.abcde:
          clade = chbp(
               ['a', 'b', 'c', 'd', 'e'],
               [expand('00000'),
                expand('00110'),
                expand('01000'),
                expand('01110'),
                expand('01111')] )
          print (f'{clade}')        
  
     if args.rosalind:
          Input           = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
          species         = Input[0].split()
          character_table = [expand(row) for row in Input[1:]]
          clade           = chbp(species,character_table)
          print (clade)
          with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
               f.write(f'{clade}\n')
                
     elapsed = time.time()-start
     minutes = int(elapsed/60)
     seconds = elapsed-60*minutes
     print (f'Elapsed Time {minutes} m {seconds:.2f} s')    


#class Clade():
    #@staticmethod
    #def toString(Root):
        #def bfs(clade,suffix=''):
            #if len(clade.children)>0:
                #pp = [bfs(child) for child in clade.children]
                #ppp=[pppp for pppp in pp if len(pppp)>0]
                #ll = ",".join(ppp) 
                #return f'({ll}){suffix}'
            #else:
                #return ','.join(clade.species)
            
        #return bfs(Root,suffix=';')    
    #def __init__(self,species):
        #self.species = species
        #self.children = []
    
    #def bfs(self,prefix=''):
        #print (f'{prefix}{self.species}')
        #for child in self.children:
            #child.bfs(prefix+'-')
  
    #def newick(self,suffix=''):
        #def wrap(value):
            #return f'({value})' if len(value)>1 else value
        #if len(self.children)==0:
            #if len(self.species)==0:
                #return ''
            #elif len(self.species)==1:
                #return self.species[0]
            #else:
                #return '(' + ','.join(self.species) + ')'
        #values      = [child.newick()  for child in self.children]
        #non_trivial = [v for v in values if len(v)>0 and len(v[0])>0]
        #if len(non_trivial)==0:
            #return ''
        #if len(non_trivial)==1:
            #return non_trivial[0]
        #return '(' + ','.join(non_trivial) + ')' + suffix

#def chbp(species,character_table):
    #def entropy(freq):
        #p1 = freq/n
        #p2 = 1-p1
        #return - p1 *np.log(p1) - p2 * np.log(p2)
    
    #n               = len(species)
    #entropies       = [entropy(sum(char)) for char in character_table]
    #indices         = np.argsort(entropies)

    #Root            = Clade(list(range(n)))
    #Current         = [Root]
    #for k in range(len(character_table)):
        #character = character_table[indices[len(character_table)-k-1]]
        #l = 0
        #r = 0
        #for clade in Current:
            #Left = []
            #Right = []
            #for i in clade.species:
                #if character[i]==0:
                    #Left.append(i)
                #else:
                    #Right.append(i)
            #l=max(l,len(Left))
            #r=max(r,len(Right))
            #clade.children = [Clade(Left),Clade(Right)]
            #clade.species = None
        #Current = [child for clade in Current for child in clade.children] 
        #print (k,entropies[indices[len(character_table)-k-1]],l,r)
    #return Root.newick(suffix=';')
