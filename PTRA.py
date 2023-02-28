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

from Bio.Seq import translate
import sys

with open ('c:/Users/Weka/Downloads/rosalind_ptra(2).txt') as f:
    i=0
    for line in f:
        print (line.strip())
        if i==0:
            coding_dna=line.strip()
        if i==1:
            out=line.strip()
            for table in [1,2,3,4,5,6,9,10,11,12,13,14,15]:
                print (table)
                translated=translate(coding_dna,table=table,to_stop=True)
#                if len(translated)==len(out):
                if translated==out:
                    print (table,translated)
                    sys.exit()
                print ('Not matched: '+ translated)
        i+=1
