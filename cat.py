#    Copyright (C) 2017 Greenweaves Software Pty Ltd
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

def CatalanNumbers(N):
    cn=[1]
    for n in range(0,N+1):
        while len(cn)<n+1:
            cn.append(sum([cn[k-1]*cn[n-k] for k in range(1,n+1)]))
    return cn[-1]

def CountNonCrossingMatches(string):
    counts=verify_counts_complete_graph(string)
    
if __name__=='__main__':
    #print (CountNonCrossingMatches('AUAU'))
    print (CatalanNumbers(12))
