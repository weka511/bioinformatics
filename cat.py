#    Copyright (C) 2019 Greenweaves Software Limited
#
#    This is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This software is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>
#
#    cat 	Catalan Numbers and RNA Secondary Structures 

from rosalind import RosalindException,verify_counts_complete_graph

# 1 2 5 14

# Verified against list from http://mathforum.org/advanced/robertd/catalan.html
#
# 1 1
# 2 2
# 3 5
# 4 14
# 5 42
# 6 132
# 7 429
# 8 1430
# 9 4862
# 10 16796
# 11 58786
# 12 208012
# 13 742900
# 14 2674440
# 15 9694845

def catalan(n):
    c=[1]
    for n0 in range(1,n+1):
        c.append(sum([c[k]*c[n0-1-k] for k in range(n0)]))
    return c[n]

def cat(s):
    return catalan(int(len(s)/2))
    
if __name__=='__main__':
#    print (cat('AUAU')%1000000)
    for n in range(1,16):
        print (n,catalan(n))
