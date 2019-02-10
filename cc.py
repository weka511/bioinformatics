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
#    cc 	Connected Components

# Create adjacency list, with both forward and backward links

def create_adjacency(graph):
    m,n       = graph[0]
    product   = {}

    for a in range(1,m+1):
        product[a]   = [a]
        
    for a,b in graph[1:]:
        product[a].append(b)
        product[b].append(a)
    for a in range(1,m+1):
        product[a]=sorted(list(set(product[a])))

    return m,n,product

def explore(a,adjacency,explored,component):
    explored.add(a)
    component.append(a)
    for b in adjacency[a]:
        if not b in explored:
            explore(b,adjacency,explored,component)
    return sorted(list(set(component)))
        
def cc(graph):          
    m,_,adjacency = create_adjacency(graph)
  
    explored   = set()

    components = {}  # The connected components, with one element as key
    for a,_ in adjacency.items():
        if not a in explored:
            component = explore(a,adjacency,explored,[])
            for c in component:
                components[c]=component
                    
    count = 0
    uniques = []
    duplicates = []
    for k,v in components.items():
        if k in uniques:
            duplicates.append(k)
        else:
            for vv in v:
                uniques.append(vv)
            count+=1
    for d in duplicates:
        del components[d]

#   a few sanity checks

    nodes = sorted(list(set([v for k,vs in components.items() for v in vs])))
    assert len(nodes) ==m, '{0} not {1}'.format(len(nodes), m) 
    v0=0
    for v in nodes:
        assert v == v0+1
        v0=v
        
    return count,components

if __name__=='__main__':
    graph=[]
    with open(r'C:\Users\Simon\Downloads\rosalind_cc(12).txt') as f:
        for line in f:
            ll =line.strip().split()
            graph.append((int(ll[0]),int(ll[1])))
        count,_ = cc(graph)
        print (count)
 