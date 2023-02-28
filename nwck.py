#!/usr/bin/env python
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

tokenizer = newick.Tokenizer()
parser = newick.Parser(tokenizer)

def get_path_to_root(clade,path=[]):
    if (len(path)==0):
        path=[(clade.id,clade.length)]
    ancestor = clade.parent
    if ancestor == None:
        return path
    pp=path[:]
    pp.append((ancestor.id,ancestor.length))
    return get_path_to_root(ancestor,pp)

def get_path(clade):
    return [clade.id]+get_path_to_root(clade)

def diff(path1,path2):
    def length(path):
        return sum(l for (_,l) in path)
    if len(path1)>len(path2):
        return diff(path2,path1)

    i=0
    while i<len(path1)  and path1[i][0]==path2[i][0]:
        i+=1

    return length(path1[i:]) +length(path2[i:])

with open (r'C:\Users\Weka\Downloads\rosalind_nkew.txt') as f:
    diffs=[]
    i = 0
    tree = None
    lookup = None
    for line in f:
        if i%3==0:
            print (line.strip())
            tree,lookup=parser.parse(line.strip())
        elif i%3==1:
            clades = line.strip().split()
            print (clades)
            ancestors=[get_path_to_root(lookup[clade])[::-1] for clade in clades]
            diffs.append(diff(ancestors[0],ancestors[1]))
        i+=1
    print (diffs)
