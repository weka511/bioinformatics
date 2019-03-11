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

def dij(g):
    n,_             = g[0]
    open_edges      = [edge for edge in g[1:]]
    D               = [-1 for i in range(n+1)]
    D[1]            = 0
    for i in range(n):
        for a,b,w in open_edges:
            if D[a]>-1:
                proposed_distance = D[a]+w
                if D[b]==-1 or D[b]>proposed_distance:
                    D[b]=proposed_distance
                
    return D[1:]

if __name__=='__main__':
    print (dij([
        (6, 10),
        (3, 4, 4),
        (1, 2, 4),
        (1, 3, 2),
        (2, 3, 3),
        (6, 3, 2),
        (3, 5, 5),
        (5, 4, 1),
        (3, 2, 1),
        (2, 4, 2),
        (2, 5, 3)  
    ]))
    
#def explore(start,adj,D):
    #for node,w in adj[explore]:
        #length_this_path = D[start-1] + w 
        #if D[node-1]==-1 or D[node-1]>length_this_path:
            #D[node-1] = length_this_path
        
#adj = {}
#n,_ = g[0]
#for i in range(1,n+1):
    #adj[i]=[]
#for (a,b,w) in g[1:]:
    #adj[a].append((b,w))    