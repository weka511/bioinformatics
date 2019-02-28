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
# BA6I Implement Graph to Genome

from BA6G import CycleToChromosome

def GraphToGenome(GenomeGraph):
    P = []
    for  Nodes in GenomeGraph:
        P.append( CycleToChromosome(Nodes))
    return P

if __name__=='__main__':
    print (GraphToGenome([(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)]))