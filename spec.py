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

from reference_tables import amino_acids
from bisect import bisect
def create_lookup(amino_acids=amino_acids):
    pairs = sorted([(abbrev,value.mon_mass) for abbrev,value in amino_acids.items()],
                   key =lambda x:x[1])
    masses = [mass for (_,mass) in pairs]
    return masses,pairs
    
def spectrum2protein(ms):
    masses,pairs = create_lookup()
    diffs = [m1-m0 for m0,m1 in zip(ms[:-1],ms[1:])]
    def get_abbrev(diff):
        index = bisect(masses,diff)
        m1 = masses[index]
        m0 = masses[(index-1) if index>0 else 0]
        if diff-m0 < m1-diff:
            index-=1
        abbrev,_ = pairs[index]
        return abbrev
    return ''.join([get_abbrev(diff) for diff in diffs])
  

if __name__=='__main__':
    for key,value in amino_acids.items():
        print (key,value.mon_mass)
    print (spectrum2protein([3524.8542,
                     3710.9335,
                     3841.974,
                     3970.0326,
                     4057.0646]))