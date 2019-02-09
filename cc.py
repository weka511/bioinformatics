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
    m,n      = graph[0]
    product  = {}
    for a in range(1,m+1):
        product[a] = [a]
    for a0,b0 in graph[1:]:
        a = min(a0,b0)
        b = max(a0,b0)
        product[a].append(b)
    return product

def cc(graph):          
    adjacency = create_adjacency(graph)
  
    for a,bs in adjacency.items():
        print (a,bs)
 
    components = {}  # The connected components, with one element as key
    index      = {}  # Which component does node belong to?
        
    for a,bs in adjacency.items():
        # see whether any elements have been processed before,
        # is so, use the earlier component
        a0 = a
        for b in bs:
            if b in index:
                a0 = index[b]
                break
        # Add items one by one
        for b in bs:
            if not b in index:
                if not a0 in components:
                    components[a0] = [a0]
                if not b in components[a0]:
                    components[a0].append(b)
                    index[b] = a0

    # Verify every node in connected or is a singleton
    singletons=[c for c,d in components.items() if len(d)==1]
    for s in singletons:
        bs=adjacency[s]
        assert len(bs)==1
        for a,bs in adjacency.items():
            assert not s in bs or a == s
            
    print ('S',len(singletons),singletons)
    connected = list(set([item for _,d in components.items() if len(d)>1 for item in d]))
    print ('C',len(connected), connected)
    sum = set(singletons + connected)
    assert len(sum) == len(singletons) + len(connected),\
        "Singletons + connected should span set"
    
    for s in singletons:
        if s in connected:
            del components[s]
            assert False, "I think this is redundant"
            
    count = 0
    for c,d in components.items():
        print (c,d)
        count +=1
        
    return count,components

if __name__=='__main__':
    graph=[]
    with open(r'C:\Users\Simon\Downloads\rosalind_cc(5).txt') as f:
        for line in f:
            ll =line.strip().split()
            graph.append((int(ll[0]),int(ll[1])))
        c,_=cc(graph)
        print (c)
        
    #print (len(cc([(10 ,0),
#(1, 2),
#(2, 8),
#(4, 10),
#(5, 9),
#(6, 10),
#(7, 9)
        #]).keys()))