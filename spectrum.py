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
#    Utilities for mass sprctroscopy

from reference_tables import integer_masses

def SpectrumGraph(spectrum):
    def add(index=-1):
        value = spectrum[index] if index>-1 else 0
        for j in range(index+1,len(spectrum)):
            diff = spectrum[j]-value
            if diff in reversed:
                if not value in product:
                    product[value]=[]
                for protein in reversed[diff]:
                    product[value].append((spectrum[j],protein))
    reversed={}
    for k,v in integer_masses.items():
        if not v in reversed:
            reversed[v]=[k]
    product = {}
    add()
    for i in range(len(spectrum)):
        add(i)
    return product