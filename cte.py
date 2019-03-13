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

# CTE Shortest Cycle Through a Given Edge
#
# Returns: For each graph, output the length of a shortest cycle going
#           through the first specified edge if there is a cycle and "-1" otherwise.


def cte(edges,n=25):
    def create_adj():
        product = {}
        for a,b,w in edges[1:]:
            if not a in product:
                product[a]=[]
            product[a].append((b,w))
        return product
    adj = create_adj()        

    def explore(start,node,path=[],length=0):
        total_length = -1
        if not node in adj: return -1
        linked = adj[node]
        for succ,weight in linked:
            if succ==start:
                return length+weight
            if succ == node: continue
            if succ in path: continue
            total_length=explore(start,succ,path+[node],length+weight)
            if total_length>0:
                return total_length
        return -1

    explored   = set()
   
    a,b,w  = edges[1]
    max_length = explore(start=a,node=b,path=[a],length=w)
    print (max_length)
            
    #start,b,w          = edges[1]

    #partial_cycles = [([start,b],w)]
    #cycles         = []
    
    #while len(partial_cycles)>0:
        #partial_cycles.sort(key=lambda tup: len(tup[0]),reverse=True)
        #ccc,_ = partial_cycles[0]
        #print (len(ccc),len(partial_cycles))        
        #new_cycles=[]
        #for (cc,length) in partial_cycles[:n]:
            #if cc[-1] in adj:
                #for next_node,weight in adj[cc[-1]]:
                    #if not next_node in cc[:-1]:
                        #new_cycles.append((cc+[next_node],length+weight))
        
        ##ccc,_ = new_cycles[0]
        ##print (len(ccc),len(new_cycles))
        #leaders = []
        #for (cc,length) in new_cycles:
            ##print (cc[-1],cc[0])
            #if cc[-1]==cc[0]:
                #cycles.append((cc,length))
            ##elif not cc[-1] in cc[:-1]:
            #else:
                #leaders.append((cc,length))
        #for ccc in partial_cycles[n:]:
            #leaders.append(ccc)
        #partial_cycles = leaders
            ##else:
                ##print ("dropped",cc[-1],cc[:-1])
    #return min([w for (_,w) in cycles]) if len(cycles)>0 else -1


if __name__=='__main__':
    from helpers import create_strings

    #k = 2
    #gs = [[(4, 5),
           #(2, 4, 2),
           #(3, 2, 1),
           #(1, 4, 3),
           #(2, 1, 10),
           #(1, 3, 4)],

          #[(4, 5),
           #(3, 2, 1),
           #(2, 4, 2),
           #(4, 1, 3),
           #(2, 1, 10),
           #(1, 3, 4)]]
    #print (' '.join([str(cte(g)) for g in gs]))   
    gs = []
    g  = []
    for line in create_strings():
        numbers = [int(s) for s in line.split(' ')]
        if len(numbers)==1:
            pass
        elif len(numbers)==2:
            if len(g)>0:
                gs.append(g)
            g = [(numbers[0],numbers[1])]
        else:
            g.append((numbers[0],numbers[1],numbers[2]))
    if len(g)>0:
        gs.append(g)           
 
        
    print (' '.join([str(cte(g)) for g in gs])) # OK for 1,2
  