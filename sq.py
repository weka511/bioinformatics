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

# SQ Square in a Graph
#

def sq(g):
    def explore(node,path=[]):
        if len(path)>4: return -1
 
        for next_node in adj[node]:
            if next_node==path[0] and len(path)==4:
                #print (path+[next_node])
                return 1
            if next_node in path: continue
            if explore(next_node,path+[next_node])==1: return 1
        return -1
    
        
    adj = {}
    for a,b in g[1:]:
        if not a in adj: adj[a]=[]
        adj[a].append(b)
        if not b in adj: adj[b]=[]
        adj[b].append(a)
    
    for node in adj.keys():
        if explore(node,path=[node])==1:
            return 1
        
    return -1

if __name__=='__main__':
    from helpers import create_strings
    
    k=2
    
    graphs = [[(4, 5),
               (3, 4),
               (4, 2),
               (3, 2),
               (3, 1),
               (1, 2)],
    
              [(4, 4),
                (1, 2),
                (3, 4),
                (2, 4),
                (4, 1)]]    

    print (' '.join([str(sq(g)) for g in graphs]))
    
    gs = []
    g  = []
    for line in create_strings(ext=1):
        if len(line)==0:
            if len(g)>0:
                gs.append(g)
            g=[]
            continue
        numbers = [int(s) for s in line.split(' ')]
        if len(numbers)==1:
            continue
        elif len(numbers)==2:
            
            g.append((numbers[0],numbers[1]))

    if len(g)>0:
        gs.append(g)           
 
        
    print (' '.join([str(sq(g)) for g in gs]))    
  