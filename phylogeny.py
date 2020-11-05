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

# Phylogeny -- http://rosalind.info/problems/topics/phylogeny/

def chbp(species,character_table):
    pass

# cstr 
#
#  Creating a Character Table from Genetic Strings  http://rosalind.info/problems/cstr/

def cstr(strings):
    def trivial(split):
        if len(split)<2: return True
        for k,v in split.items():
            if v<2: return True
        return False
    
    choices = [[] for s in strings[0]]
    counts  = [{} for s in strings[0]]
    for i in range(len(strings[0])):
        for s in strings:
            if not s[i] in choices[i]:
                choices[i].append(s[i])
            if s[i] in counts[i]:
                counts[i][s[i]] += 1
            else:
                counts[i][s[i]] = 1
    
    splits=[]
    for i in range(len(strings[0])):
        split = {}
        for c in choices[i]:
            split[c] = 0
        for s in strings:
            for c in choices[i]:
                if s[i]==c:
                    split[c]+=1
        splits.append(split)
    result=[]
    for i in range(len(strings[0])):
        character = []
        split = splits[i]
        if not trivial(split):
            chs = list(split.keys())
            for s in strings:
                character.append('0' if s[i]==chs[0] else '1')  
            result.append(''.join(character))
    return result

# ctbl Creating a Character Table  http://rosalind.info/problems/ctbl/

def CharacterTable(tree):
    def create_character(split_species):
        character=[]
        for s in species:
            character.append(1 if s in split_species else 0)
        return ''.join([str(c) for c in character])
    
    species=[spec.name for spec in tree.find_elements(terminal=True)]
    species.sort()

    clades=[clade for clade in tree.find_clades(terminal=False)]
    # we iterate over all Clades except the root
    return [create_character([spec.name for spec in split.find_elements(terminal=True)]) for split in clades[1:]]