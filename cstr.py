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
    def c(i):
        choices= list(set([s[i] for s in strings]))
        return [len([s[i] for s in strings if s[i]==choice]) for choice in choices]
    counts= ([c(i) for i in range(len(strings[0]))])
    print (counts)
    for i in range(len(strings)):
        for j in range(i):
            print (i,j,hamm(strings[i], strings[j]))
                   
if __name__=='__main__': 
    cstr([
        'ATGCTACC',
        'CGTTTACC',
        'ATTCGACC',
        'AGTCTCCC',
        'CGTCTATC'    
    ])
