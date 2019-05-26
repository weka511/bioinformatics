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
#    sgra Using the Spectrum Graph to Infer Peptides

from reference_tables import amino_acids
from spectrum import create_lookup,get_abbrev

def sgra(L,Alphabet=amino_acids,epsilon=0.00001):
    def create_spectrum_graph():
        masses,pairs   = create_lookup()
        G = {}
        for u in L:
            G[u]=[]
        for u in L:
            for v in L:
                if u<v:
                    abbrev = get_abbrev(v-u,masses,pairs)
                    if abs(v-u-amino_acids[abbrev].mon_mass)<epsilon:
                        G[u].append(v)
                    
        return G
    G = create_spectrum_graph()


if __name__=='__main__':
    print (sgra([3524.8542,
                 3623.5245,
                 3710.9335,
                 3841.974,
                 3929.00603,
                 3970.0326,
                 4026.05879,
                 4057.0646,
                 4083.08025
                 ]))