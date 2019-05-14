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
#    full Inferring Peptide from Full Spectrum

from spectrum import create_lookup,get_abbrev
from reference_tables import amino_acids

def full(s,epsilon=0.00001):
    def create_candidates():
        masses,pairs = create_lookup()
        candidates={}
        for diff in [m1-m0 for m0 in s[1:] for m1 in s[1:] if m1>m0]:
            a = get_abbrev(diff,masses,pairs)
            m = amino_acids[a].mon_mass
            if abs(diff-amino_acids[a].mon_mass)<epsilon:
                if a in candidates:
                    candidates[a]+=1
                else:
                    candidates[a]=1
        return candidates
        
    def create_pairs():
        pairs = []
        for i in range(1,len(s)//2+1):
            assert(abs(s[i]+s[-i]-s[0])<epsilon)
            pairs.append((s[i],s[-i]))
        return pairs
    total = s[0]
    pairs = create_pairs()
    candidates= create_candidates()
 
    protein = []
    return ''.join(protein)

if __name__=='__main__':
    print (full([1988.21104821,
                 610.391039105,
                 738.485999105,
                 766.492149105,
                 863.544909105,
                 867.528589105,
                 992.587499105,
                 995.623549105,
                 1120.6824591,
                 1124.6661391,
                 1221.7188991,
                 1249.7250491,
                 1377.8200091
    ]))
