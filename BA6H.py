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
# BA6H Implement ColoredEdges

from BA6F import ChromosomeToCycle
def ColouredEdges(P):
    Edges = []
    for  Chromosome in P:
        #print (Chromosome)
        Nodes = ChromosomeToCycle(Chromosome)
        #print (Nodes)
        it = iter(Nodes[1:]+[Nodes[0]])
        for i in it:
            Edges.append((i,next(it)))
    return Edges

if __name__=='__main__':
    print (ColouredEdges([[+1 ,-2 ,-3],[+4 ,+5, -6]]))