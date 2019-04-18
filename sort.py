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
#    SORT 	Sorting by Reversals

import time
from fragile import sort
from helpers import create_strings

if __name__=='__main__':
    def parse(line):
        return [int(c) for c in line.strip().split()]
    
    S = create_strings(ext=3)
    d,path = sort(S[0].split(' '),S[1].split(' '))
    print (d)
    for (i,j) in path:
        print ('{0} {1}'.format(i,j))