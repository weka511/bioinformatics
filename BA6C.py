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

from fragile import count_synteny_blocks,get_synteny_blocks,ColouredEdges

#def extract_all_edges(chromosome):
    #def extract_edges(graph):
        #return list(zip(graph,graph[1:]+graph[0:]))
    #return [extract_edges(graph) for graph in chromosome]
    
def d2break(a,b):
    def update(index,edges,x,y,max_node):
        edges.append((x,y))
        if x in index:
            index[x].append(y)
        else:
            index[x]=[y] 
        if x>max_node:
            max_node=x
        return max_node
    
    def build_cycle(start,i,index,cycle):
        maybe_start = False
        cycle.append(i)
        for link in index[i]:
            if link==start:
                maybe_start = True
            else:
                if not link in cycle:
                    build_cycle(start,link,index,cycle)
                    maybe_start = False
        index.pop(i)
        
    blocks = get_synteny_blocks(a)
    n      = count_synteny_blocks(a)
    nb     = count_synteny_blocks(b)
    assert (n==nb),'Mismatched synteny blocks {0} != {1}'.format(n,nb)

    edges = []
    index = {}
    max_node = -1
    for (x,y) in ColouredEdges(a) + ColouredEdges(b):       
        max_node = update(index,edges,x,y,max_node)
        max_node = update(index,edges,y,x,max_node)
 
    print (edges)
    print (index)
    print (max_node)
    
    cycles = []
    for i in range(1,max_node+1):
        if i in index:
            cycles.append(build_cycle(i,i,index,[]))
    return n - len(cycles)

if __name__=='__main__':
    print (
        d2break(
            [[+1, +2, +3, +4, +5, +6]],
            [[+1, -3, -6, -5],[+2, -4]]
    ))
