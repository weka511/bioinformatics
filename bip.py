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

# BIP Testing Bipartiteness

# A graph is bipartite if and only if it does not contain an odd cycle
# https://en.wikipedia.org/wiki/Bipartite_graph#Characterization

from helpers import format_list, create_adjacency

def bip(graph):

    def explore(node,adjacency,path,explored,edges):
         
        #print ('X',node)
        #explored.add(node)
        for linked in adjacency[node]:
            if linked==node: continue
            
            if  (node,linked) in edges or (linked,node) in edges: continue
           
            if len(path)>0 and linked==path[0]:
                cycle     = path + [node]
                found_odd = len(cycle)%2==1
                if found_odd:
                    print (cycle)
                explored.add(node)
                return found_odd
            else:
                if linked in explored:
                    explored.add(node)
                    return False
                if explore(linked,adjacency,path+[node],explored,edges+[(node,linked)]):
                    explored.add(node)
                    return True
        explored.add(node)        
        return False
   
    m,_,adjacency = create_adjacency(graph,back=True)
 
    cycles   = {}
    explored = set()
   
    for node in range(1,m+1):
        if not node in explored:
            print ("Explore ",node)
            try:
                if explore(node,adjacency,[],explored,[]): return -1
            except:
                return 1
    return 1


def parse_graphs(f):
    product = []        
    graph   = []
    state   = 0
    n       = -1
    for line in f:
        ll =line.strip()
        if state==0:
            if len(ll)>0:
                n = int(ll)
                state = 1
                continue
        if state==1:
            if len(ll)==0:
                state = 2
                graph = []
                continue
        if state == 2:
            if len(ll)==0:
                product.append(graph)
                graph=[]
            else:
                lll =ll.split()
                graph.append((int(lll[0]),int(lll[1])))
    product.append(graph) 
    assert len(product)==n,'{0} {1}'.format(len(product),n)
    for graph in product:
        a,b=graph[0]
        assert len(graph)==b+1,'{0} {1}'.format(len(graph),b)
    return product

if __name__=='__main__':
    simple = [[[3, 3],
               [1, 2],
               [3, 2],
               [3, 1]],
        
              [[4, 3],
               [1, 4],
               [3, 1],
               [1, 2]]]
    print (format_list([bip(g) for g in simple]))
    with open(r'C:\Users\Simon\Downloads\rosalind_bip(4).txt') as f:
        #parse_graphs(f)
        print (format_list([bip(g) for g in parse_graphs(f)]))
        # -1 1 1 -1 -1 1 -1 1 1 -1 -1 -1 1 1 1 1 1 1 1