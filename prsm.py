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
#
#    prsm Matching a Spectrum to a Protein

from spectrum import prsm


if __name__=='__main__':
    from helpers import create_strings



    i = 0
    s = []
    R = []
    for ll in create_strings(ext=3):
        if i ==0:
            n = int(ll)
        elif i<n+1:
            s.append(ll)
        else:
            R.append(float(ll))
        i+=1

    m,s_max = prsm(s,R)

    print (m)
    print (s_max)
