#!/usr/bin/env python

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

# BA5I Find a Highest-Scoring Overlap Alignment of Two Strings

from align import create_distance_matrix,overlap_assignment

if __name__=='__main__':
    from helpers import create_strings

    strings  = create_strings('ba5i',ext=1)
    d,s1,t1  = overlap_assignment(strings[0],strings[1])
    print ('{0}'.format(d))
    print (s1)
    print (t1)
