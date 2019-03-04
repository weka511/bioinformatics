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
#    BA6C Compute the 2-Break Distance Between a Pair of Genomes

from fragile import count_synteny_blocks,get_synteny_blocks

def extract_all_edges(chromosome):
    def extract_edges(graph):
        return list(zip(graph,graph[1:]+graph[0:]))
    return [extract_edges(graph) for graph in chromosome]
    
def d2break(a,b):
    def cycles():
        return 0
    blocks = get_synteny_blocks(a)
    n      = count_synteny_blocks(a)
    nb     = count_synteny_blocks(b)
    assert (n==nb),'Mismatched synteny blocks {0} != {1}'.format(n,nb)
    e1     = extract_all_edges(a)
    e2     = extract_all_edges(b)
    print (e1)
    print (e2)
    return n-cycles()

if __name__=='__main__':
    print (
        d2break(
            [[+1, +2, +3, +4, +5, +6]],
            [[+1, -3, -6, -5],[+2, -4]]
    ))
