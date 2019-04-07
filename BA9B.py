# Copyright (C) 2019 Greenweaves Software Limited

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

# BA9B Implement TrieMatching

from snp import MatchTries
from rosalind import trie

if __name__=='__main__':
    s = 'AATCGGGTTCAATCGGGGT'
    t = trie(['ATCG','GGGT'],one_based=False)
    print (MatchTries(s,t))
    #with open('/Users/Simon/Downloads/rosalind_ba9b.txt') as f:
        #strings=[]

        #text = f.readline().strip()
        #for line in f:
            ##if line.startswith('Output'): break
            #strings.append(line.strip())
        #matches = MatchTries(text,trie(strings,one_based=False))
        #print (' '.join([str(m) for m in matches]))