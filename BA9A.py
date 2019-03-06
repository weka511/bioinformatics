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

# BA9A Construct a Trie from a Collection of Patterns 

from rosalind import trie

if __name__=='__main__':
    #for a,b,c in trie(['ATAGA','ATC','GAT'],one_based=False):
        #print ('{0}->{1}:{2}'.format(a,b,c))
    with open('/Users/Simon/Downloads/rosalind_ba9a.txt') as f:
        strings=[]
        #f.readline()
        for line in f:
            #if line.startswith('Output'): break
            strings.append(line.strip())
        for a,b,c in trie(strings,one_based=False):
            print ('{0}->{1}:{2}'.format(a,b,c))        