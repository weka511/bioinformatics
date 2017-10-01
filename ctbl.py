# Copyright (C) 2017 Greenweaves Software Pty Ltd

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>
 

from Bio import Phylo
import random

def CharacterTable(tree,split_species_len=0.25):
    def create_character():
        character=[]
        for s in species:
            character.append(1 if s in split_species else 0)
        return ''.join([str(c) for c in character])
    
    #print(tree)
    species=[spec.name for spec in tree.find_elements(terminal=True)]
    species.sort()
    #print (species)
    clades=[clade for clade in tree.find_clades(terminal=False)]
    #print (clades)
    result=[]
    for split in clades[1:]:
        split_species=[spec.name for spec in split.find_elements(terminal=True)]
        #print (split_species)
        result.append(create_character())
    return result

if __name__=='__main__': 
    #tree = Phylo.read('ctbl.txt', 'newick')
    tree = Phylo.read('c:/Users/Weka/Downloads/rosalind_ctbl(3).txt', 'newick')
    for ch in CharacterTable(tree):
        print (ch)
