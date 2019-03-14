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

# HDAG Hamiltonian Path in DAG 
#

from helpers import create_adjacency

def hdag(graph):
    def explore(start,path=[]):
        if len(path)==m: return 1,path
        if start in adjacency:
            for next_node in adjacency[start]:
                if next_node in path: continue
                result,path = explore(next_node,path=path+[next_node])
                if result==1 and len(path)==m: return result,path
                if result==-1: return -1,[]
            pass
        else:
            return -1,[]
        return -1,[]
    
    m,_,adjacency = create_adjacency(graph,back=False,self=False)
    
    for start in range(1,m+1):
        result,path = explore(start,path=[start])
        if result>0:
            return result,path
        
    return  result,path

if __name__=='__main__':
    from helpers import create_strings
    
    k=2
    
    graphs = [[(3,3),
               (1,2),
               (2,3),
               (1, 3)],
    
              [(4, 3),
                (4,3),
                (3, 2),
                (4, 1)]]    

    for g in graphs:
        result,path=hdag(g)
        if result==-1:
            print (result)
        else:
            print (' '.join(str(x) for x in [result]+path) )
    
    #gs = []
    #g  = []
    #for line in create_strings(ext=1):
        #if len(line)==0:
            #if len(g)>0:
                #gs.append(g)
            #g=[]
            #continue
        #numbers = [int(s) for s in line.split(' ')]
        #if len(numbers)==1:
            #continue
        #elif len(numbers)==2:
            
            #g.append((numbers[0],numbers[1]))

    #if len(g)>0:
        #gs.append(g)           
 
        
    #print (' '.join([str(sq(g)) for g in gs]))    
  