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

import newick

frequencies={
    'aa':(0,0,1),
    'Aa': (0,1,0),
    'AA': (1,0,0)
}

factors=[
   [[1,0,0], [0.5,0.5,0], [0,1,0]],
   [[0.5,0.5,0],[0.25, 0.5, 0.25],[0, 0.5, 0.5]],
   [[0, 1,0],[0,0.5,0.5],[0,0,1]]
]

def combine(f1,f2):
    def total(i,contribs):
        return sum(c[i] for c in contribs)
    contribs=[]
    for i in range(3):
        for j in range(3):
            contribs.append([f*f1[i]*f2[j] for f in factors[i][j]])
    return [total(k,contribs) for k in range(3)]

def mend(node):
    if len(node.nodes)==0:
        try:
            if node.name=='aA':
                node.name=node.name[::-1]
            freqs=frequencies[node.name]
            return freqs
        except KeyError:
            return (0,0)
    parent_freqs = [mend(parent) for parent in node.nodes]
    parent_freqs=[pp for pp in parent_freqs if len(pp)==3]
    return combine(parent_freqs[0],parent_freqs[1])

tokenizer = newick.Tokenizer()
parser = newick.Parser(tokenizer)

with open (r'C:\Users\Weka\Downloads\rosalind_mend.txt') as f:
    for line in f:
        tree,_=parser.parse(line.strip())
        print (mend(tree))



