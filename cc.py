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

def create_adjacency(graph):
    m,n       = graph[0]
    product   = {}
    backlinks = {}
    for a in range(1,m+1):
        product[a]   = [a]
        backlinks[a] = [a]
    for a,b in graph[1:]:
        product[min(a,b)].append(max(a,b))
        backlinks[max(a,b)].append(min(a,b))
    for a in range(1,m+1):
        product[a].sort()
        if len(backlinks[a])>1:
            backlinks[a].sort()
        else:
            backlinks.pop(a)
    return m,n,product,backlinks

def explore(a,adjacency,explored,component):
#    print ("Exploring",a)
    explored.add(a)
    component.append(a)
    for b in adjacency[a]:
        if not b in explored:
            explore(b,adjacency,explored,component)
    return sorted(list(set(component)))
        
def cc(graph):          
    m,_,adjacency,backlinks = create_adjacency(graph)
  
    explored   = set()

    components = {}  # The connected components, with one element as key
    for a,_ in adjacency.items():
        if not a in explored:
            component = explore(a,adjacency,explored,[])
            for c in component:
                components[c]=component
# 826 [64, 88, 90, 142, 259, 306, 310, 362, 364, 826]
# 64 [26, 40, 42, 46, 58, 59, 61, 64, 72, 90, 92, 98, 100, 106, 107, 111, 128, 131, 147, 153, 154, 155, 168, 169, 177, 180, 193, 203, 210, 214, 227, 232, 236, 240, 251, 257, 259, 279, 295, 304, 305, 306, 310, 321, 330, 337, 342, 350, 352, 368, 376, 382, 385, 389, 392, 400, 402, 414, 427, 434, 461, 472, 484, 502, 538, 544, 550, 580, 582, 586, 602, 626, 677, 697, 715, 774, 788, 797, 798, 802, 814, 816]

    while True:
        merges = []
        for k,vs in components.items():
            for v in vs:
                if v in backlinks:
                    for b in backlinks[v]:
                        if not b in vs:
                            merges.append((k, backlinks[v]))
                            break
 
        print (len(merges))
        if len(merges)==0: break 
 
        merged=[]
        
        for k, backs in merges:
            if k in merged: continue
            merged.append(k)
            components[k] = sorted(list(set( components[k]+ backs)))
            for v in  components[k]:
                components[v]= components[k]
                merged.append(v)
                
        if len(merges)==0: break

                    
    count = 0
    uniques = []
    duplicates = []
    for k,v in components.items():
        if k in uniques:
            duplicates.append(k)
        else:
            for vv in v:
                uniques.append(vv)
            print(k,v)
            count+=1
    for d in duplicates:
        del components[d]
   
    nodes = sorted(list(set([v for k,vs in components.items() for v in vs])))
    assert len(nodes) ==m, '{0} not {1}'.format(len(nodes), m) 
    v0=0
    for v in nodes:
        assert v == v0+1
        v0=v
    return count

if __name__=='__main__':
    graph=[]
    with open(r'C:\Users\Simon\Downloads\rosalind_cc(10).txt') as f:
        for line in f:
            ll =line.strip().split()
            graph.append((int(ll[0]),int(ll[1])))
        c=cc(graph)
        print (c)
 