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

def cc(graph):
    m,n        = graph[0]
    components = []
    outs       = []
    unexplored = graph[1:]
    count      = 0
    def explore(i,component):
        nonlocal unexplored, outs
        if not i in component:
            component.append(i)
        #print (i)
        outs.append(i)
        reachable = []
        ins       = []
        residue   = []
        for a,b in unexplored:
            if a == i:
                if not b in reachable:
                    reachable.append(b)
                    if not b in component:
                        component.append(b)                    
            else:
                residue.append((a,b))
        unexplored = residue
        for b in reachable:
            component= explore(b,component)
        return component
 
    for i in range(1,m):
        if not i in outs:
 #           print ("top:", i)
            print (explore(i,[]))
            count+=1
    return count

if __name__=='__main__':
    graph=[]
    with open(r'C:\Users\Simon\Downloads\rosalind_cc.txt') as f:
        for line in f:
            ll =line.strip().split()
            graph.append((int(ll[0]),int(ll[1])))
        print (cc(graph))
        
    #print(
        #cc([
            #(12, 13),
            #(1, 2),
            #(1, 5),
            #(5, 9),
            #(5, 10),
            #(9, 10),
            #(3, 4),
            #(3, 7),
            #(3, 8),
            #(4, 8),
            #(7, 11),
            #(8, 11),
            #(11, 12),
            #(8, 12)    
    #]))
    