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

# cstr Creating a Character Table from Genetic Strings  http://rosalind.info/problems/cstr/

from rosalind import hamm

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
 
if __name__=='__main__': 
    #for row in cstr([
        #'ATGCTACC',
        #'CGTTTACC',
        #'ATTCGACC',
        #'AGTCTCCC',
        #'CGTCTATC'    
        #]):
        #print (row)
    with open('c:/Users/Simon/Downloads/rosalind_cstr.txt') as f:
        strings=[]
        for line in f:
            strings.append(line.strip())
        for row in cstr(strings):
            print (row)
