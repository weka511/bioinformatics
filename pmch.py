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
#    pmch 	Perfect Matchings and RNA Secondary Structures 

from rosalind import RosalindException
def NumberPerfectMatchings(string):
    def factorial(n):
        return n*factorial(n-1) if n>1 else 1
    counts={
        'A':0,
        'U':0,
        'C':0,
        'G':0
    }
    for c in string:
        counts[c]+=1
    if counts['A']==counts['U'] and counts['C']==counts['G']:
        return factorial(counts['A'])*factorial(counts['G'])
    else:
        raise RosalindException('Mismatched counts {0}/{1} or {2}/{3}'.format(counts['A'],counts['U'],counts['C'],counts['G']))
    
if __name__=='__main__':
    print (NumberPerfectMatchings('CGCUUGCAGGACAAGGUGCAUAUCUUCACCAUCAGCUAAACGAUAGCCGCUCACGCGGGAUGGUGCGUUCCGUG'))
